/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
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

Application
    reactingParcelSWIMFoam

Group
    grpLagrangianSolvers

Description
    Transient solver for compressible, turbulent flow with a reacting,
    multiphase particle cloud, and surface film modelling.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "basicReactingMultiphaseCloud.H"
#include "surfaceFilmModel.H"
#include "rhoReactionThermo.H"
#include "CombustionModel.H"
#include "radiationModel.H"
#include "SLGThermo.H"
#include "fvOptions.H"
#include "pimpleControl.H"
#include "pressureControl.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
//    #include "createTimeControls.H"
    #include "createFields.H"
    
//    #include "createFieldRefs.H"
//    #include "createRegionControls.H"
//    #include "initContinuityErrs.H"

//    Info<<"0 "<<endl;
    const polyPatch& film=mesh.boundaryMesh()[filmPatchID];
//    Info<<"1 "<<intFilm<<endl;
    const labelUList& filmFaceCells=film.faceCells();
    const polyPatch& intFilm=(*rmesh).boundaryMesh()[intfilmPatchID];
//    Info<<"1 "<<intFilm<<endl;
    const labelUList& intFilmFaceCells=intFilm.faceCells();
//    Info<<"2 "<<intFilmFaceCells<<endl;
    forAll(intFilm, faceI)
    {
	const label faceID=intFilm.start()+faceI;
//        Info<<faceI<<" "<<intFilmFaceCells[faceI]<<endl;
//        Info<<faceI<<" "<<(*rmesh).faces()[faceID]<<endl;
	const face& intFace=(*rmesh).faces()[faceID];
	label findI=-1;
	forAll(film, fI)
	{
	    const label fID=film.start()+fI;
	    const face& priFace=mesh.faces()[fID];
	    bool flagall=true;
	    forAll(intFace, pointID) 
	    {
    	        bool flag=false;
    	        const point& intPoint=(*rmesh).points()[intFace[pointID]];
//    	        Info<<pointID<<" int "<<intPoint<<endl;;
		if (intFace.size()==priFace.size())
		{
		    forAll(priFace, pID)
		    {
			const point& priPoint=mesh.points()[priFace[pID]];
//			Info<<pID<<" pri "<<priPoint<<endl;
			if (fabs(intPoint.x()-priPoint.x())<1e-8) 
			{
//			    Info<<"x "<<intPoint.x()<<" "<<priPoint.x()<<endl;
			    if (fabs(intPoint.y()-priPoint.y())<1e-8)
			    {
//				Info<<"y "<<intPoint.y()<<" "<<priPoint.y()<<endl;
				if (fabs(intPoint.z()-priPoint.z())<1e-8) 
				{
//				    Info<<"z "<<intPoint.z()<<" "<<priPoint.z()<<endl;
				    flag=true;
				}
			    }
			}
//			Info<<flag<<endl;
			if (flag==true) break;
		    }
		}
		if (flag==false) {
		    flagall=false;
		    break;
		}
	    }
	    if (flagall==true) {
		findI=fI;
		break;
	    }
	}
//	Info<<faceI<<" "<<findI<<endl;
	if (findI>=0) {
	    filmCellDist[intFilmFaceCells[faceI]]=cellDist[filmFaceCells[findI]];
	    filmCellDecomposition[intFilmFaceCells[faceI]]=filmCellDist[intFilmFaceCells[faceI]];
	}
    }
    filmCellDist.correctBoundaryConditions();
    filmCellDist.write();
    filmCellDecomposition.write();
    
//    Info<<"2 "<<intFilmFaceCells<<endl;
/*    forAll(film, faceI)
    {
        Info<<faceI<<" "<<filmFaceCells[faceI]<<endl;
    }*/
/*    turbulence->validate();

    if (!LTS)
    {
        #include "compressibleCourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "compressibleCourantNo.H"
            #include "setMultiRegionDeltaT.H"
        }

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        parcels.evolve();
//	volScalarField& he = thermo.he();
	TParcels=dimensionedScalar("", dimTemperature, parcels.Tmin());
//(1.0/runTime.deltaT())*(parcels.Sh(he) & TParcels)*(1.0/parcels.constProps().Cp0());
        surfaceFilm.evolve();

        if (solvePrimaryRegion)
        {
            if (pimple.nCorrPIMPLE() <= 1)
            {
                #include "rhoEqn.H"
            }

            // --- PIMPLE loop
            while (pimple.loop())
            {
                #include "UEqn.H"
                #include "YEqn.H"
                #include "EEqn.H"

                // --- Pressure corrector loop
                while (pimple.correct())
                {
                    #include "pEqn.H"
                }

                if (pimple.turbCorr())
                {
                    turbulence->correct();
                }
            }

            rho = thermo.rho();
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }
*/
    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
