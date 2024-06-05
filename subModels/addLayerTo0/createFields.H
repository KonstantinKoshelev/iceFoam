//*
#include "createRDeltaT.H"

Info<< "Reading thermophysical properties\n" << endl;
autoPtr<rhoReactionThermo> pThermo(rhoReactionThermo::New(mesh));
rhoReactionThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

SLGThermo slgThermo(mesh, thermo);

basicSpecieMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

const word inertSpecie(thermo.lookup("inertSpecie"));
if (!composition.species().found(inertSpecie))
{
    FatalIOErrorIn(args.executable().c_str(), thermo)
        << "Inert specie " << inertSpecie << " not found in available species "
        << composition.species()
        << exit(FatalIOError);
}

Info<< "Creating field rho\n" << endl;
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);
volScalarField TParcels
(
    IOobject
    (
        "TParcels",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    thermo.T()
);


volScalarField& p = thermo.p();

Info<< "\nReading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "compressibleCreatePhi.H"

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

Info<< "Creating combustion model\n" << endl;
autoPtr<CombustionModel<rhoReactionThermo>> combustion
(
    CombustionModel<rhoReactionThermo>::New(thermo, turbulence())
);

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Force p_rgh to be consistent with p
p_rgh = p - rho*gh;

pressureControl pressureControl(p, rho, pimple.dict(), false);

mesh.setFluxRequired(p_rgh.name());

Info<< "Creating multi-variate interpolation scheme\n" << endl;
multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(thermo.he());

volScalarField Qdot
(
    IOobject
    (
        "Qdot",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero)
);

#include "createDpdt.H"

#include "createK.H"

#include "createMRF.H"
#include "createRadiationModel.H"
#include "createClouds.H"
#include "createSurfaceFilmModel.H"
#include "createFvOptions.H"
//*/
    Info<< "\nReading cellDecomposition\n" << endl;
    labelIOList cellDecomposition
    (
        IOobject
        (
            "cellDecomposition",
	    runTime.constant(),     // dict is found in "constant"
            mesh,                   // registry for the dict            
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );
    
/*    volScalarField cellDist
    (
	IOobject
        (
	    "cellDist",
    	    runTime.timeName(),
    	    mesh,
    	    IOobject::MUST_READ,
    	    IOobject::NO_WRITE
        ),
	mesh
    );*/
//    Info<<cellDist<<endl;

    IOdictionary moveMeshesDict
    (
        IOobject
        (
            "moveMeshesDict",    // dictionary name
            runTime.constant(),     // dict is found in "constant"
            mesh,                   // registry for the dict
            IOobject::MUST_READ,    // must exist, otherwise failure
            IOobject::NO_WRITE      // dict is only read by the solver
        )
    );
    word mappedWallPatchName=word(moveMeshesDict.lookup("mappedWallPatchName")); //"region0_to_wallFilmRegion_wallFilmFaces";
    Info<<"Mapped patch: "<<mappedWallPatchName<<endl;
//    label filmPatchID=mesh.boundaryMesh().findPatchID(mappedWallPatchName);
//    vectorField& pDisp=refCast<vectorField>(pD.boundaryFieldRef()[filmPatchID]);
/*    Info<< "\nConstructing surface film model" << endl;

//    #include "readGravitationalAcceleration.H"
    autoPtr<regionModels::surfaceFilmModel> tsurfaceFilm
    (
	regionModels::surfaceFilmModel::New(mesh, g)
    );
*/
//    fvMesh* rmesh = runTime.getObjectPtr<fvMesh>("wallFilmRegion");
//    label intfilmPatchID=rmesh->boundaryMesh().findPatchID(mappedWallPatchName);
//    Info<<"b"<<endl;

    word filmTopPatchName=word(moveMeshesDict.lookup("filmTopPatchName")); //"wallFilmFaces_top";
    Info<<"Film top patch: "<<filmTopPatchName<<endl;
//    label inttopfilmPatchID=rmesh->boundaryMesh().findPatchID(filmTopPatchName);
//    Info<<"c"<<endl;
/*
    volScalarField filmCellDist
    (
	IOobject
        (
	    "filmCellDist",
    	    runTime.timeName(),
    	    *rmesh,
    	    IOobject::NO_READ,
    	    IOobject::NO_WRITE
        ),
	*rmesh,
    dimensionedScalar("filmCellDist", dimless, -1),
    zeroGradientFvPatchScalarField::typeName
    );

// Write decomposition as labelList for use with 'manual'
// decomposition method.
    labelIOList filmCellDecomposition
    (
        IOobject
        (
            "filmCellDecomposition",
            rmesh->facesInstance(),
            *rmesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        labelList(rmesh->nCells(), label(-1))
    );
*/
//                cellDecomposition.write();
