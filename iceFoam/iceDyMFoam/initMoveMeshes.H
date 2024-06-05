    regionModels::surfaceFilmModels::iceSingleLayer* filmIce=dynamic_cast<regionModels::surfaceFilmModels::iceSingleLayer*> (&surfaceFilm);

    pointVectorField& pD=const_cast<pointVectorField&>(mesh.objectRegistry::lookupObject<pointVectorField>("pointDisplacement"));

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
    label filmPatchID=mesh.boundaryMesh().findPatchID(mappedWallPatchName);
    vectorField& pDisp=refCast<vectorField>(pD.boundaryFieldRef()[filmPatchID]);

    fvMesh* rmesh = runTime.getObjectPtr<fvMesh>("wallFilmRegion");
    label intfilmPatchID=rmesh->boundaryMesh().findPatchID(mappedWallPatchName);

    word filmTopPatchName=word(moveMeshesDict.lookup("filmTopPatchName")); //"wallFilmFaces_top";
    Info<<"Film top patch: "<<filmTopPatchName<<endl;
    label inttopfilmPatchID=rmesh->boundaryMesh().findPatchID(filmTopPatchName);
    pointField newPs(rmesh->points().size());
    pointField newPsNF(mesh.points().size());

/*
// cylinder radius
    scalar r=readScalar(moveMeshesDict.lookup("cylinderRadius")); //0.076; 
    Info<<"Cylinder radius: "<<r<<endl;
// film thickness (must be equal thickness in extrudeToRegionMeshDict)
    scalar filmThickness=readScalar(moveMeshesDict.lookup("filmThickness")); //0.005;
    Info<<"Film thickness: "<<filmThickness<<nl<<endl;
    scalar rtop=r-filmThickness;
*/    
    bool moveMeshes=readBool(moveMeshesDict.lookup("moveMeshes"));
    Info<<"Move meshes: "<<moveMeshes<<endl;
    
    const labelListList& meshPointPoints=mesh.pointPoints();
    const labelList& filmPoints=mesh.boundaryMesh()[filmPatchID].meshPoints();
    labelList pointsNearFilm(filmPoints.size());
    labelList pointsNearFilm2(filmPoints.size());
    labelList pointsNearFilm3(filmPoints.size());

    forAll(filmPoints, index) 
    {
        label i=filmPoints[index];
        pointsNearFilm[index]=-1;
        forAll (meshPointPoints[i], iPoint) 
        {
    	    bool flag = false;
    	    forAll (filmPoints, j)
    	    {
    		if (meshPointPoints[i][iPoint] == filmPoints[j]) 
    		{
    		    flag = true;
    		    break;
    		}
    	    }
    	    if (!flag)
    	    {
    		pointsNearFilm[index] = meshPointPoints[i][iPoint];
    		break;
    	    }
        }
    }
    forAll(pointsNearFilm, index)
    {
	label i = pointsNearFilm[index];
	label count = 0;
        pointsNearFilm2[index]=-1;
//	pointsNearFilm3[index]=-1;
        forAll (meshPointPoints[i], iPoint) 
        {
    	    bool flag = false;
    	    forAll (filmPoints, j)
    	    {
    		if (meshPointPoints[i][iPoint] == filmPoints[j]) 
    		{
    		    flag = true;
    		    break;
    		}
    	    }
    	    if (!flag)
    	    {
    		forAll(pointsNearFilm, k)
    		{
    		    if (meshPointPoints[i][iPoint] == pointsNearFilm[k]) 
    		    {
    			flag = true;
    			break;
    		    }
    		}
    	    }
    	    if (!flag)
    	    {
    		pointsNearFilm2[index] = meshPointPoints[i][iPoint];
    		count++;
    	    }
    	}
	if (count>1) {
	    pointsNearFilm2[index]=-1;
//	    pointsNearFilm3[index]=-1;
	}
    }
    forAll(pointsNearFilm2, index)
    {
	label i = pointsNearFilm2[index];
	label count = 0;
	pointsNearFilm3[index]=-1;
        if (pointsNearFilm2[index]>=0)
	{
        forAll (meshPointPoints[i], iPoint) 
        {
    	    bool flag = false;
    	    forAll (pointsNearFilm, j)
    	    {
    		if (meshPointPoints[i][iPoint] == pointsNearFilm[j]) 
    		{
    		    flag = true;
    		    break;
    		}
    	    }
    	    if (!flag)
    	    {
    		forAll(pointsNearFilm2, k)
    		{
    		    if (meshPointPoints[i][iPoint] == pointsNearFilm2[k]) 
    		    {
    			flag = true;
    			break;
    		    }
    		}
    	    }
    	    if (!flag)
    	    {
    		pointsNearFilm3[index] = meshPointPoints[i][iPoint];
    		count++;
    	    }
    	}
	}
	if (count>1) {
//	    pointsNearFilm2[index]=-1;
	    pointsNearFilm3[index]=-1;
	}
    }

    scalar angleCrit = moveMeshesDict.lookupOrDefault("angleCrit", 0.0);
    Info<<"angleCrit: "<<angleCrit<<" grad"<<endl;
    angleCrit = constant::mathematical::pi * angleCrit / 180.0;
    scalar angleWeight =  moveMeshesDict.lookupOrDefault("angleWeight", 0.3);
    Info<<"angleWeight: "<<angleWeight<<endl;
    scalar maxAngleWeight =  moveMeshesDict.lookupOrDefault("maxAngleWeight", 0.0);
    if (maxAngleWeight<angleWeight)
	maxAngleWeight = angleWeight;
    Info<<"maxAngleWeight: "<<maxAngleWeight<<endl;
    scalar maxPartDeltaIce = moveMeshesDict.lookupOrDefault("maxPartDeltaIce", 0.75);
    Info<<"maxPartDeltaIce: "<<maxPartDeltaIce<<endl;
    scalar magWeight = moveMeshesDict.lookupOrDefault("magWeight", 0.25);
    Info<<"magWeight: "<<magWeight<<endl;
    label reportOrtho = moveMeshesDict.lookupOrDefault("reportOrtho", 0);
    Info<<"reportOrtho: "<<reportOrtho<<endl;
    label iReportOrtho = 0;
//    Info<<pointsNearFilm<<endl;
    const labelListList& intMeshPointPoints=rmesh->pointPoints();
    const labelList& intfilmPoints=rmesh->boundaryMesh()[intfilmPatchID].meshPoints();
    labelList topFilmPoints(intfilmPoints.size());
    pointField filmNormal(intfilmPoints.size());
    labelList nearTopFilmPoints(intfilmPoints.size());

    forAll(intfilmPoints, index) 
    {
        label i=intfilmPoints[index];
        topFilmPoints[index]=-1;
        forAll (intMeshPointPoints[i], iPoint) 
        {
    	    bool flag = false;
    	    forAll (intfilmPoints, j)
    	    {
    		if (intMeshPointPoints[i][iPoint] == intfilmPoints[j]) 
    		{
    		    flag = true;
    		    break;
    		}
    	    }
    	    if (!flag)
    	    {
    		topFilmPoints[index] = intMeshPointPoints[i][iPoint];
    		break;
    	    }
        }
    }

    forAll(nearTopFilmPoints, index) 
    {
//        label i=nearTopFilmPoints[index];
        nearTopFilmPoints[index]=-1;
	forAll(intfilmPoints, jndex)
	{
	    if (intfilmPoints[jndex]==index)
	    {
		nearTopFilmPoints[index]=jndex;
		break;
	    }
	}
    }

//    Info<<topFilmPoints<<endl;
