/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd
-------------------------------------------------------------------------------
    License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "Myers.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "fvcSnGrad.H"
#include "fvcReconstruct.H"
#include "fvcFlux.H"
#include "fvm.H"
#include "zeroGradientFvPatchFields.H"
#include "mixedFvPatchFields.H"
#include "mappedFieldFvPatchField.H"
#include "mapDistribute.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"

// Sub-models
#include "filmThermoModel.H"
#include "filmViscosityModel.H"
#include "heatTransferModel.H"
#include "phaseChangeModel.H"
#include "filmRadiationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Myers, 0);

addToRunTimeSelectionTable(surfaceFilmRegionModel, Myers, mesh);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

wordList Myers::hsBoundaryTypes()
{
    wordList bTypes(T_.boundaryField().types());
    forAll(bTypes, patchi)
    {
        if
        (
            T_.boundaryField()[patchi].fixesValue()
         || isA<mixedFvPatchScalarField>(T_.boundaryField()[patchi])
         || isA<mappedFieldFvPatchField<scalar>>(T_.boundaryField()[patchi])
        )
        {
            bTypes[patchi] = fixedValueFvPatchField<scalar>::typeName;
        }
    }

    return bTypes;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Myers::read()
{
    // No additional properties to read
    return kinematicSingleLayer::read();
}


void Myers::resetPrimaryRegionSourceTerms()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    kinematicSingleLayer::resetPrimaryRegionSourceTerms();

    hsSpPrimary_ == dimensionedScalar(hsSp_.dimensions(), Zero);
}


void Myers::correctThermoFields()
{
    rho_ == filmThermo_->rho();
    sigma_ == filmThermo_->sigma();
    Cp_ == filmThermo_->Cp();
    kappa_ == filmThermo_->kappa();
}


void Myers::correctHsForMappedT()
{
    T_.correctBoundaryConditions();

    volScalarField::Boundary& hsBf = hs_.boundaryFieldRef();

    forAll(hsBf, patchi)
    {
        const fvPatchField<scalar>& Tp = T_.boundaryField()[patchi];
        if (isA<mappedFieldFvPatchField<scalar>>(Tp))
        {
            hsBf[patchi] == hs(Tp, patchi);
        }
    }
}


void Myers::updateSurfaceTemperatures()
{
    correctHsForMappedT();

    // Push boundary film temperature into wall temperature internal field
    for (label i=0; i<intCoupledPatchIDs_.size(); i++)
    {
        label patchi = intCoupledPatchIDs_[i];
        const polyPatch& pp = regionMesh().boundaryMesh()[patchi];
        UIndirectList<scalar>(Tw_, pp.faceCells()) =
            T_.boundaryField()[patchi];
    }
    Tw_.correctBoundaryConditions();

    // Update film surface temperature
    Ts_ = T_;
    Ts_.correctBoundaryConditions();
}


void Myers::transferPrimaryRegionThermoFields()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    kinematicSingleLayer::transferPrimaryRegionThermoFields();
    impSp_=rhoSp_;

    // Update primary region fields on local region via direct mapped (coupled)
    // boundary conditions
    TPrimary_.correctBoundaryConditions();
    forAll(YPrimary_, i)
    {
        YPrimary_[i].correctBoundaryConditions();
    }
    TParcels_.correctBoundaryConditions();
//    Info<<" min/max TParcels "<<gMin(TParcels_)<<", "<<gMax(TParcels_)<<endl;
//    Info<<" min/max Ts "<<gMin(Ts_)<<", "<<gMax(Ts_)<<endl;
//    Info<<" min/max Tw "<<gMin(Tw_)<<", "<<gMax(Tw_)<<endl;
}


void Myers::transferPrimaryRegionSourceFields()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    kinematicSingleLayer::transferPrimaryRegionSourceFields();

    volScalarField::Boundary& hsSpPrimaryBf =
        hsSpPrimary_.boundaryFieldRef();

    // Convert accumulated source terms into per unit area per unit time
    const scalar deltaT = time_.deltaTValue();
    forAll(hsSpPrimaryBf, patchi)
    {
        scalarField rpriMagSfdeltaT
        (
            (1.0/deltaT)/primaryMesh().magSf().boundaryField()[patchi]
        );

        hsSpPrimaryBf[patchi] *= rpriMagSfdeltaT;
    }

    // Retrieve the source fields from the primary region via direct mapped
    // (coupled) boundary conditions
    // - fields require transfer of values for both patch AND to push the
    //   values into the first layer of internal cells
//    Info<<"before min, max hsSp_ "<<gMin(hsSp_)<<", "<<gMax(hsSp_)<<nl;
    hsSp_.correctBoundaryConditions();
//    Info<<"after min, max hsSp_ "<<gMin(hsSp_)<<", "<<gMax(hsSp_)<<nl;
}


void Myers::correctAlpha()
{
    if (hydrophilic_)
    {
        const scalar hydrophilicDry = hydrophilicDryScale_*deltaWet_;
        const scalar hydrophilicWet = hydrophilicWetScale_*deltaWet_;

        forAll(alpha_, i)
        {
            if ((alpha_[i] < 0.5) && (delta_[i] > hydrophilicWet))
            {
                alpha_[i] = 1.0;
            }
            else if ((alpha_[i] > 0.5) && (delta_[i] < hydrophilicDry))
            {
                alpha_[i] = 0.0;
            }
        }

        alpha_.correctBoundaryConditions();
    }
    else
    {
        alpha_ ==
            pos0(delta_ - dimensionedScalar("deltaWet", dimLength, deltaWet_));
    }
}


void Myers::updateSubmodels()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    // Update heat transfer coefficient sub-models
    htcs_->correct();
    htcw_->correct();

    // Update radiation
    radiation_->correct();

    // Update injection model - mass returned is mass available for injection
    injection_.correct(availableMass_, cloudMassTrans_, cloudDiameterTrans_);

    phaseChange_->correct
    (
        time_.deltaTValue(),
        availableMass_,
        primaryMassTrans_,
        primaryEnergyTrans_
    );

    const volScalarField rMagSfDt((1/time().deltaT())/magSf());

    // Vapour recoil pressure
    pSp_ -= sqr(rMagSfDt*primaryMassTrans_)/(2*rhoPrimary_);

    // Update transfer model - mass returned is mass available for transfer
    transfer_.correct(availableMass_, primaryMassTrans_, primaryEnergyTrans_);

    // Update source fields
    rhoSp_ += rMagSfDt*(cloudMassTrans_ + primaryMassTrans_);
    hsSp_ += rMagSfDt*(cloudMassTrans_*hs_ + primaryEnergyTrans_);
  
    turbulence_->correct();
}


 tmp<fvScalarMatrix> Myers::q(volScalarField& hs) const
{
    return
    (
        // Heat-transfer to the primary region
      - fvm::Sp(htcs_->h()/Cp_, hs)
      + htcs_->h()*(hs/Cp_ + alpha_*(TPrimary_ - T_))

        // Heat-transfer to the wall
      - fvm::Sp(htcw_->h()/Cp_, hs)
      + htcw_->h()*(hs/Cp_ + alpha_*(Tw_- T_))
//      + htcw_->h()*(hs/Cp_ + (TPrimary_ - T_))
    );
}

/*tmp<volScalarField> Myers::qq(volScalarField& hs) const
{
    volScalarField r = htcw_->h()*(hs/Cp_ + (TPrimary_ - T_));
    forAll(r, i)
    {
	if (deltaIce_[i]<0.0001)
	{
	  r[i] = 0.0;
	}
    }
    return (r+htcs_->h()*(hs/Cp_ + alpha_*(TPrimary_ - T_)));
}*/


void Myers::solveContinuity()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

//    Info<<"SWIM cont rhoSp_ "<<gMin(rhoSp_)<<", "<<gMax(rhoSp_)<<nl;
	solve
    (
        fvm::ddt(deltaRho_)
      + fvc::div(phi_)
     ==
      - rhoSp_
      + iceSp_
    );	
    delta_ = deltaRho_ / rho_;
	// Bound film thickness by a minimum of zero
    delta_.max(1e-8);
}

/*
void Myers::solveThickness
(
    const volScalarField& pu,
    const volScalarField& pp//,
    //const fvVectorMatrix& UEqn
)
{
    if (debug)
    {
        InfoInFunction << endl;
    }   
	
	// employ simple coeff-based model	
	scalar Cf_ = turbulence_->coeffDict().get<scalar>("Cf");    
    volScalarField Cs("Cs", Cf_*rhoPrimary_ *mag(UPrimary_ - U_));
	volVectorField shearStress = (
									Cs * U_ + Cs * UPrimary_// surface contribution
									);
                                    
    U_ = 
       - delta_ * (
				fvc::reconstruct
				 (
				   fvc::interpolate(delta_)
				   * ( regionMesh().magSf()
						* (
							fvc::snGrad(pu, "snGrad(p)")
							+ fvc::snGrad(pp, "snGrad(p)")*fvc::interpolate(delta_)
							+ fvc::snGrad(delta_)*fvc::interpolate(pp)
						   )
					 )
				 ) 
			- delta_ * rho_*gTan()
			)/ (3 * mu_) 
       + delta_ *
			(
				shearStress
				//turbulence_->Su(U_) & U_ 
		    ) / (2 * mu_)
    ; 
    
    // Bound film thickness by a minimum of zero
    //delta_.max(0.0);

    // Remove any patch-normal components of velocity
    U_ -= nHat()*(nHat() & U_);

    U_.correctBoundaryConditions();

    // Update film wall and surface velocities
    updateSurfaceVelocities();
    
    phi_ = fvc::flux(deltaRho_ * U_);

    // Continuity check
    continuityCheck();
}
*/



void Myers::solveEnergy()
{
   /* Info << "OUTPUT ->" << nl;
    Info << "Tref = " << Tref << nl;
    Info << "Tw_ = " << Tw_ << nl;
    Info << "Tf_ = " << Tf_ << nl; 
    Info << "htcs_->h() = " << htcs_->h() << nl;
    Info << "htcw_->h() = " << htcw_->h() << nl;
    Info << "U8LWC = " << U8LWC << nl;*/
    
    //Info << "Cp_ = " << Cp_ << nl;
    
    //Info << "hsSp_ = " << hsSp_ << nl;
    //Info << "rhoSp_ = " << rhoSp_ << nl;
    //Info << "hsSp_/Cp_ + rhoSp_* Tref  = " << hsSp_/Cp_ + rhoSp_* Tref<< nl;
    
    /*if (min(rhoSp_) < 0*min(rhoSp_))
    {
        //Info << "hsSp_/min(rhoSp_) = " << hsSp_/min(rhoSp_) << nl;
        
        //Info << "T(hsSp_/min(rhoSp_)) = " << T(hsSp_/min(rhoSp_)) << nl;
    }*/
    
    //Info << - (hsSp_  + rhoSp_ * Cp_ * Tref) << nl;
    //Info << htcs_->h() * TPrimary_ << nl;
    //Info << rhoSp_*mag(UPrimary_)*mag(UPrimary_)/2 << nl;
    //Info << rhoSp_*mag(UPrimary_)*mag(UPrimary_)/2 + htcs_->h() * 0.8 * mag(UPrimary_)*mag(UPrimary_)/CpIce_ << nl;
    
    if (debug)
    {
        InfoInFunction << endl;
    }
	    
    dimensionedScalar Ca("Ca", dimEnergy/dimMass/dimTemperature, 1010);
//    dimensionedScalar Hscale("Hscale", dimLength, 0.0001);
    volScalarField q_0 = (
			    - (hsSp_  + rhoSp_ * Cp_ * (TParcels_ - Tf_))
			    + htcs_->h() * TPrimary_
                + htcs_->h()*mag(UPrimary_)*mag(UPrimary_)/2/Ca
			) / kappa_;
    
    volScalarField q_1 = (
			    -htcs_->h()
//				+ rhoSp_ * Cp_
			)  / kappa_;
			
	
    volScalarField q_0r = (
			    - (hsSp_  + rhoSp_ * Cp_ * (TParcels_ - Tf_))
			    - rhoSp_ * Lf_
			    + htcs_->h() * TPrimary_
                + htcs_->h()*mag(UPrimary_)*mag(UPrimary_)/2/Ca
			) / kappaIce_;
    
    volScalarField q_1r = (
			    - htcs_->h()
//				+ rhoSp_ * Cp_
			) / kappaIce_;
//    Info << "q_0 = " << q_0 << nl;
//    Info << "q_1 = " << q_1 << nl;
//    Info << "q_0r = " << q_0r << nl;
//    Info << "q_1r = " << q_1r << nl;
    
	Ts_ = Tw_ + (q_0r + q_1r * Tw_) / (1.0 - q_1r * deltaIce_) * deltaIce_;
    
    TIce = Ts_;
  
	
	//volScalarField 
    Twater = Tf_ + (q_0 + q_1 * Tf_ ) / (1.0 - q_1 * delta_) * delta_;    
    
    volScalarField Bg = (Tf_ - TPrimary_) / (q_0r + q_1r * Tf_);
//    Info << "Bg = " << Bg << nl;
    volScalarField Wh = pos(pos0(TPrimary_ - Tf_) + alpha_);
    forAll(Wh, i)
    {
	if (deltaIce_[i] < 0.0001)
	{
	    Wh[i] = 0.0;
	}
	else
	{
	    if (deltaIce_[i] < mag(Bg[i]))
		Wh[i] = 0.0;
	    else
		Wh[i] = 1.0;
	}
    }

    
    /*
    Info << "Ts_ = " << Ts_ << nl;
    Info << "Twater = " << Twater << nl;

    Info << "TPrimary_ = " << TPrimary_ << nl;
    Info << "Tf-Tw = " << Tf_ - TPrimary_ << nl;
    */

    iceSp_=
	(
		(1.0 - Wh) *
		    (
				rhoSp_
		    )
		- Wh * 
			(
				kappaIce_ * (Tf_ - TPrimary_) / ( (deltaIce_ + dimensionedScalar("smallDeltaIce", dimLength, SMALL)) * Lf_)
				- kappa_ / Lf_ * (q_0 + q_1 * Tf_ ) / (1.0 - q_1 * delta_)
			)
    );
//    Info << "deltaIce = " << deltaIce_ << nl;
//    Info << "iceSp = " << iceSp_ << nl;
    iceSp_.min(0.0);
    solve
    (
        fvm::ddt(rhoIce_, deltaIce_)
     ==
      - iceSp_
    );
//    deltaIce_.max(0.0);

    correctThermoFields();

    // Evaluate viscosity from user-model
    viscosity_->correct(pPrimary_, T_);

    // Update film wall and surface temperatures
    updateSurfaceTemperatures();
}
/*
 tmp<volScalarField> Myers::pu()
 {
     return tmp<volScalarField>
     (
         new volScalarField
         (
             IOobject
             (
                 typeName + ":pu",
                 time_.timeName(),
                 regionMesh(),
                 IOobject::NO_READ,
                 IOobject::NO_WRITE
             ),
             pPrimary_                  // pressure (mapped from primary region)
           - pSp_                           // accumulated particle impingement
           - fvc::laplacian(sigma_, delta_ + deltaIce_) // surface tension
         )
     );
 }
  
  
 tmp<volScalarField> Myers::pp()
 {
     return tmp<volScalarField>
     (
         new volScalarField
         (
             IOobject
             (
                 typeName + ":pp",
                 time_.timeName(),
                 regionMesh(),
                 IOobject::NO_READ,
                 IOobject::NO_WRITE
             ),
            -rho_*gNormClipped() // hydrostatic effect only
         )
     );
 }
*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Myers::Myers
(
    const word& modelType,
    const fvMesh& mesh,
    const dimensionedVector& g,
    const word& regionType,
    const bool readFields
)
:
    iceSingleLayer(modelType, mesh, g, regionType, false),
    kappaIce_("kappaIce_", dimensionSet(1, 1,-3,-1,0,0,0), readScalar(coeffs_.lookup("kappaIce"))),
    CpIce_("CpIce_", dimensionSet(0, 2,-2,-1,0,0,0), readScalar(coeffs_.lookup("CpIce"))),
    Twater
    (
        IOobject
        (
            "Twater", // Water film surface temperature 
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("Twater", dimTemperature, 273.16),
        zeroGradientFvPatchScalarField::typeName
    ),
    TIce
    (
        IOobject
        (
            "TIce", // Water film surface temperature 
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("TIce", dimTemperature, 273.16),
        zeroGradientFvPatchScalarField::typeName
    )
/*    ,
    bettaLWCW_
    (
        IOobject
        (
            "bettaLWCW", // Same name as T on primary region to enable mapping
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimDensity*dimVelocity, Zero),
        this->mappedFieldAndInternalPatchTypes<scalar>()
    )*/
{
    if (coeffs().readIfPresent("Tmin", Tmin_))
    {
        Info<< "    limiting minimum temperature to " << Tmin_ << endl;
    }

    if (coeffs().readIfPresent("Tmax", Tmax_))
    {
        Info<< "    limiting maximum temperature to " << Tmax_ << endl;
    }

    if (thermo_.hasMultiComponentCarrier())
    {
        YPrimary_.setSize(thermo_.carrier().species().size());

        forAll(thermo_.carrier().species(), i)
        {
            YPrimary_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        thermo_.carrier().species()[i],
                        time().timeName(),
                        regionMesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    regionMesh(),
                    dimensionedScalar(dimless, Zero),
                    pSp_.boundaryField().types()
                )
            );
        }
    }

    if (hydrophilic_)
    {
        coeffs_.lookup("hydrophilicDryScale") >> hydrophilicDryScale_;
        coeffs_.lookup("hydrophilicWetScale") >> hydrophilicWetScale_;
    }

    if (readFields)
    {
        transferPrimaryRegionThermoFields();

        correctAlpha();

        correctThermoFields();

        // Update derived fields
        hs_ == hs(T_); //-hs(Tf_);

        deltaRho_ == delta_*rho_;

        surfaceScalarField phi0
        (
            IOobject
            (
                "phi",
                time().timeName(),
                regionMesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE,
                false
            ),
            fvc::flux(deltaRho_*U_)
        );

        phi_ == phi0;

        // Evaluate viscosity from user-model
        viscosity_->correct(pPrimary_, T_);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Myers::~Myers()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Myers::addSources
(
    const label patchi,
    const label facei,
    const scalar massSource,
    const vector& momentumSource,
    const scalar pressureSource,
    const scalar energySource
)
{
    kinematicSingleLayer::addSources
    (
        patchi,
        facei,
        massSource,
        momentumSource,
        pressureSource,
        energySource
    );

    if (debug)
    {
        Info<< "    energy   = " << energySource << nl << endl;
    }

    hsSpPrimary_.boundaryFieldRef()[patchi][facei] -= energySource;
}


/*void Myers::preEvolveRegion()
{

//    bettaLWCW_.correctBoundaryConditions();
    const scalar deltaT = time_.deltaTValue();

    

	const scalar tmp(1.0);
    const vector tmpv(0.0,0.0,0.0);
    for (label i=0; i<primaryPatchIDs_.size(); i++)
    {
        label patchi = primaryPatchIDs_[i];
        const polyPatch& pp = regionMesh().boundaryMesh()[patchi];

        scalarField rpriMagSfdeltaT
        (
            deltaT*primaryMesh().magSf().boundaryField()[patchi]
        );
        
        
        scalarField massFlow = rhoSp_*rpriMagSfdeltaT;
        
        //Info << "rpriMagSfdeltaT" << rpriMagSfdeltaT<< endl;

        //Info << "bettaLWCW_" << bettaLWCW_ << endl;
  
         forAll(pp, k)
         {
             scalar energySource = 0 * massFlow[k] * Cp_[k]* (TParcels_[k] - 298.15);
             //Info << hsSpPrimary_.boundaryFieldRef()[pp.index()][k] << endl;
             //Info << pp[k] << endl;
             //Info << pp.index() << k << massFlow[k] << tmp * UPrimary_[k] << tmp * 0 << energySource << endl;
             addSources(pp.index(), k, massFlow[k], massFlow[k] * UPrimary_[k], tmp * 0, energySource);
         }
         //const label patchi = pp.index();
     }

    if (debug)
    {
        InfoInFunction << endl;
    }

    kinematicSingleLayer::preEvolveRegion();
    primaryEnergyTrans_ == dimensionedScalar(dimEnergy, Zero);
}
*/
void Myers::evolveRegion()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    // Solve continuity for deltaRho_
    //solveContinuity();

    // Update sub-models to provide updated source contributions
    updateSubmodels();    
    

    //Info <<"!!! MIN/MAX hsSp_" << min(hsSp_) << max(hsSp_) << endl;
    //Info <<"!!! MIN/MAX rhoSp_" << min(rhoSp_) << max(rhoSp_) << endl;
    //Info <<"!!! MIN/MAX hsSp_" << min(hsSp_) << max(hsSp_) << endl;

    // Explicit pressure source contribution
    tmp<volScalarField> tpu(this->pu());

    // Implicit pressure source coefficient
    tmp<volScalarField> tpp(this->pp());

    // Solve for momentum for U_
    tmp<fvVectorMatrix> UEqn = solveMomentum(tpu(), tpp());

    // Solve energy for hs_ - also updates thermo
    solveEnergy();
	
	// Explicit pressure source contribution
    //tmp<volScalarField> tpu(this->pu());

    // Implicit pressure source coefficient
    //tmp<volScalarField> tpp(this->pp());
	
	// Solve continuity for deltaRho_
    solveContinuity();

    // Solve thickness for delta_
    solveThickness(tpu(), tpp(), UEqn());

    // Update deltaRho_ with new delta_
    //deltaRho_ == delta_*rho_;

    // Update temperature using latest hs_
//    T_ == T(hs_);

	const scalar dt = time_.deltaTValue();
    const scalar curt = time_.time().value() - time_.time().startTime().value();
    if (time_.timeIndex() > time_.startTimeIndex() + 1)
    {
	betta_= (1.0 - dt / curt) * betta_ 
		- rhoSp_ / U8LWC * dt / (curt - dt);
    }
    else
    {
        betta_ += -rhoSp_ / U8LWC;
    }
}

const volScalarField& Myers::Cp() const
{
    return Cp_;
}


const volScalarField& Myers::kappa() const
{
    return kappa_;
}


const volScalarField& Myers::T() const
{
    return T_;
}


const volScalarField& Myers::Ts() const
{
    return Ts_;
}


const volScalarField& Myers::Tw() const
{
    return Tw_;
}


const volScalarField& Myers::hs() const
{
    return hs_;
}

const volScalarField& Myers::deltaIce() const
{
  return deltaIce_;
}


void Myers::info()
{
    thermoSingleLayer::info();

    Info<< indent << "min/max(deltaIce)  = "
        << gMin(deltaIce_) << ", "
        << gMax(deltaIce_) << nl
        << indent << "current ice mass   = "
        << gSum((rhoIce_*deltaIce_*magSf())()) << nl;

}


tmp<volScalarField::Internal> Myers::Srho() const
{
    tmp<volScalarField::Internal> tSrho
    (
        new volScalarField::Internal
        (
            IOobject
            (
                typeName + ":Srho",
                time().timeName(),
                primaryMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            primaryMesh(),
            dimensionedScalar(dimMass/dimVolume/dimTime, Zero)
        )
    );

    scalarField& Srho = tSrho.ref();
    const scalarField& V = primaryMesh().V();
    const scalar dt = time_.deltaTValue();

    forAll(intCoupledPatchIDs(), i)
    {
        const label filmPatchi = intCoupledPatchIDs()[i];

        scalarField patchMass =
            primaryMassTrans_.boundaryField()[filmPatchi];

        toPrimary(filmPatchi, patchMass);

        const label primaryPatchi = primaryPatchIDs()[i];
        const labelUList& cells =
            primaryMesh().boundaryMesh()[primaryPatchi].faceCells();

        forAll(patchMass, j)
        {
            Srho[cells[j]] += patchMass[j]/(V[cells[j]]*dt);
        }
    }

    return tSrho;
}


tmp<volScalarField::Internal> Myers::Srho
(
    const label i
) const
{
    const label vapId = thermo_.carrierId(filmThermo_->name());

    tmp<volScalarField::Internal> tSrho
    (
        new volScalarField::Internal
        (
            IOobject
            (
                typeName + ":Srho(" + Foam::name(i) + ")",
                time_.timeName(),
                primaryMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            primaryMesh(),
            dimensionedScalar(dimMass/dimVolume/dimTime, Zero)
        )
    );

    if (vapId == i)
    {
        scalarField& Srho = tSrho.ref();
        const scalarField& V = primaryMesh().V();
        const scalar dt = time().deltaTValue();

        forAll(intCoupledPatchIDs_, i)
        {
            const label filmPatchi = intCoupledPatchIDs_[i];

            scalarField patchMass =
                primaryMassTrans_.boundaryField()[filmPatchi];

            toPrimary(filmPatchi, patchMass);

            const label primaryPatchi = primaryPatchIDs()[i];
            const labelUList& cells =
                primaryMesh().boundaryMesh()[primaryPatchi].faceCells();

            forAll(patchMass, j)
            {
                Srho[cells[j]] += patchMass[j]/(V[cells[j]]*dt);
            }
        }
    }

    return tSrho;
}


tmp<volScalarField::Internal> Myers::Sh() const
{
    tmp<volScalarField::Internal> tSh
    (
        new volScalarField::Internal
        (
            IOobject
            (
                typeName + ":Sh",
                time().timeName(),
                primaryMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            primaryMesh(),
            dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero)
        )
    );

    scalarField& Sh = tSh.ref();
    const scalarField& V = primaryMesh().V();
    const scalar dt = time_.deltaTValue();

    forAll(intCoupledPatchIDs_, i)
    {
        const label filmPatchi = intCoupledPatchIDs_[i];

        scalarField patchEnergy =
            primaryEnergyTrans_.boundaryField()[filmPatchi];

        toPrimary(filmPatchi, patchEnergy);

        const label primaryPatchi = primaryPatchIDs()[i];
        const labelUList& cells =
            primaryMesh().boundaryMesh()[primaryPatchi].faceCells();

        forAll(patchEnergy, j)
        {
            Sh[cells[j]] += patchEnergy[j]/(V[cells[j]]*dt);
        }
    }

    return tSh;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
} // End namespace regionModels
} // End namespace surfaceFilmModels

// ************************************************************************* //
