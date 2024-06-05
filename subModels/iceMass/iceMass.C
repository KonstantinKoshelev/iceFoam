/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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
    wallHeatFlux

Description
    Calculates and writes the heat flux for all patches as the boundary field
    of a volScalarField and also prints the integrated flux for all wall
    patches.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvcVolumeIntegrate.H"

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

// #include "harmonic.H"
// #include "incompressiblePhase.H"
// #include "capillarityModel.H"
// #include "relativePermeabilityModel.H"
#include "fixedValueFvPatchField.H"
#include "iceSingleLayer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
//    #include "addRegionOption.H"
//    #include "setRootCase.H"

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

//    #include "createTime.H
    regionModels::surfaceFilmModel& surfaceFilm = tsurfaceFilm();
    regionModels::surfaceFilmModels::iceSingleLayer* filmIce=dynamic_cast<regionModels::surfaceFilmModels::iceSingleLayer*> (&surfaceFilm);
    dimensionedScalar rhoIce("rhoIce_", dimensionSet(1,-3,0,0,0,0,0), readScalar(filmIce->coeffs().lookup("rhoIce")));
    Info << "rhoIce = " << rhoIce.value() << endl;

    instantList timeDirs = timeSelector::select0(runTime, args);
//    #include "createNamedMesh.H"

    Info<< "Time      Ice Mass" << endl;
    scalar volume0 = 0.0;
    scalar icemass = 0.0;
    bool first = true;
    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        mesh.readUpdate();
	if (first)
	{
	    volume0 = gSum(mesh.V());
	    first = false;
	}
	else
	{
    	    Info << runTime.timeName();
	    icemass = (volume0-gSum(mesh.V()))*rhoIce.value();
    	    Info<< "        " << icemass << endl;
	}
/*
        #include "createFields.H"

        Sb=epsScalar*(Sb-Sbmin);

        Info << "Qoil = " << fvc::domainIntegrate(Sb).value() << endl;
    
        Info<< endl;
*/
    }

    Info<< "End" << endl;

    return 0;
}

// ************************************************************************* //
