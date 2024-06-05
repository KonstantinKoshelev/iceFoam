{
    volScalarField& he = thermo.he();

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + mvConvection->fvmDiv(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), U),
                p,
                "div(phiv,p)"
            )
          : -dpdt
        )
      - fvm::laplacian(turbulence->alphaEff(), he)
     ==
        rho*(U&g)
      + parcels.Sh(he)
      + surfaceFilm.Sh()
      + radiation->Sh(thermo, he)
      + Qdot
      + fvOptions(rho, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();
    radiation->correct();

/*    scalarField& wallHeatFluxBf = wallHeatFlux.boundaryFieldRef()[filmPatchID];

    volScalarField T2=thermo.T();
    scalarField& T2Bf=T2.boundaryFieldRef()[filmPatchID];
    forAll(T2Bf, patchI)
    {
	T2Bf[patchI] = 273.15;
    }
    volScalarField he2=thermo.he(p, T2);

    const volScalarField::Boundary& heBf = he2.boundaryField();
    wallHeatFluxBf = heBf[filmPatchID].snGrad();
    const labelList& own=mesh.owner();
    const polyPatch& wallFilmPatch = mesh.boundaryMesh()[filmPatchID];
    forAll(wallFilmPatch, patchFaceI)
    {
	wallHeatFluxBf[patchFaceI] *= turbulence->alphaEff()()[own[wallFilmPatch.start()+patchFaceI]];
    }
*/
    Info<< "T gas min/max   = " << min(T).value() << ", "
        << max(T).value() << endl;
}
