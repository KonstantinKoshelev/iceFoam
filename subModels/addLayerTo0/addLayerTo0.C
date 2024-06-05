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
    
    labelList cellDecCopy = cellDecomposition;
    
    forAll(mesh.cells(), cellI)
    {
	if (cellDecCopy[cellI]!=0) {
//	    Info<<cellI<<"  "<<cellDecomposition[cellI]<<endl;
	    bool flag = false;
	    const labelList& cc = mesh.cellCells()[cellI];
	    forAll(cc, ccI)
	    {
		if (cellDecCopy[cc[ccI]]==0) 
		{
		    flag = true;
		    break;
		}
	    }
	    if (flag)
		cellDecomposition[cellI] = 0;
	}
    }
    cellDecomposition.write();
    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
