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
    volScalarField& delta=const_cast<volScalarField&>((*rmesh).objectRegistry::lookupObject<volScalarField>("deltaf"));
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
    label moveMeshes=moveMeshesDict.lookupOrDefault("moveMeshes",1);
    label moveMeshesStep = 1;
    Info<<"Move meshes: "<<moveMeshes<<endl;
    
/*    const labelListList& meshPointPoints=mesh.pointPoints();
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
*/
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
// найти общие грани на границе внешней сетки и пленки
    const polyPatch& intFilm=(*rmesh).boundaryMesh()[intfilmPatchID];
    const polyPatch& film=mesh.boundaryMesh()[filmPatchID];
    labelList extFilmToInt(film.size());
    extFilmToInt=-1;
//    const label startFilm = film.start();
//    Info<<"1 "<<intFilm<<endl;
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
			if (fabs(intPoint.x()-priPoint.x())<1e-10) 
			{
//			    Info<<"x "<<intPoint.x()<<" "<<priPoint.x()<<endl;
			    if (fabs(intPoint.y()-priPoint.y())<1e-10)
			    {
//				Info<<"y "<<intPoint.y()<<" "<<priPoint.y()<<endl;
				if (fabs(intPoint.z()-priPoint.z())<1e-10) 
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
	if (findI>=0)
        {
	    extFilmToInt[findI]=faceI;
	}
    }
//    Info<<extFilmToInt<<endl;
/*    List<List<label> > spLabels((*rmesh).globalData().sharedPoints().size());
    List<List<label> > spProcessors((*rmesh).globalData().sharedPoints().size());
    forAll((*rmesh).globalData().sharedPoints(), isp)
    {
	const point& sp = (*rmesh).globalData().sharedPoints()[isp];
	spLabels[isp] = List<label>();
	spProcessors[isp] = List<label>();
	forAll((*rmesh).globalData().sharedPointLabels(), ispl)
	{
	    const point& localsp = (*rmesh).points()[(*rmesh).globalData().sharedPointLabels()[ispl]];
	if ((fabs(sp.x()-localsp.x())<1e-10) &&
	    (fabs(sp.y()-localsp.y())<1e-10) &&
	    (fabs(sp.z()-localsp.z())<1e-10))
	{
	    spLabels[isp].append((*rmesh).globalData().sharedPointLabels()[ispl]);
	    spProcessors[isp].append(UPstream::myProcNo());
	}
	}
    }
*/
    List<List<vector> > mySharedPoints(Pstream::nProcs());
    mySharedPoints[Pstream::myProcNo()] = List<vector>(0);
    labelList mySharedPointLabels(0);
    List<List<vector> > allPs(Pstream::nProcs());
    allPs[Pstream::myProcNo()] = List<vector>(rmesh->points().size(), vector::zero);
    forAll(rmesh->points(), ipoint)
    {
	allPs[Pstream::myProcNo()][ipoint] = rmesh->points()[ipoint];
    }
    Pstream::gatherList(allPs);
    Pstream::scatterList(allPs);
    if (Pstream::master())
    {
        for(label iProc=0; iProc<Pstream::nProcs(); iProc++)
	{
	    forAll(allPs[iProc], i)
	    {
		for(label j=i+1; j<allPs[iProc].size(); j++)
		{
		    if (allPs[iProc][i]==allPs[iProc][j])
		    {
			bool flag = true;
			forAll(mySharedPoints[Pstream::myProcNo()], k)
			{
			    if (allPs[iProc][i]==mySharedPoints[Pstream::myProcNo()][k])
			    {
				flag = false;
				break;
			    }
			}
			if (flag)
			{
			    label sz = mySharedPoints[Pstream::myProcNo()].size();
			    mySharedPoints[Pstream::myProcNo()].setSize(sz+1);
			    mySharedPoints[Pstream::myProcNo()][sz] = allPs[iProc][i];
			}
		    }
		}
		for(label jProc=iProc+1; jProc<Pstream::nProcs(); jProc++)
		forAll(allPs[jProc], j)
		{
		    if (allPs[iProc][i]==allPs[jProc][j])
		    {
			bool flag = true;
			forAll(mySharedPoints[Pstream::myProcNo()], k)
			{
			    if (allPs[iProc][i]==mySharedPoints[Pstream::myProcNo()][k])
			    {
				flag = false;
				break;
			    }
			}
			if (flag)
			{
			    label sz = mySharedPoints[Pstream::myProcNo()].size();
			    mySharedPoints[Pstream::myProcNo()].setSize(sz+1);
			    mySharedPoints[Pstream::myProcNo()][sz] = allPs[iProc][i];
			}
		    }
		}
	    }
	}
    }
    Pstream::gatherList(mySharedPoints);
    Pstream::scatterList(mySharedPoints);
    List<List<label> > spLabels(Pstream::nProcs());
    spLabels[Pstream::myProcNo()] = List<label>(mySharedPoints[0].size(), -1);
    forAll(mySharedPoints[0], isp)
    {
	point sp = mySharedPoints[0][isp];
	forAll((*rmesh).points(), ispl)
	{
	    const point& lp = (*rmesh).points()[ispl];
	    if ((fabs(sp.x()-lp.x())<1e-9) &&
	        (fabs(sp.y()-lp.y())<1e-9) &&
		(fabs(sp.z()-lp.z())<1e-9))
	    {
		spLabels[Pstream::myProcNo()][isp]=ispl;
	    }
	}
    }
    Pstream::gatherList(spLabels);
    Pstream::scatterList(spLabels);
//    Info<<"spLabels "<<spLabels<<endl;
    List<List<vector> > exSharedPoints(Pstream::nProcs());
    exSharedPoints[Pstream::myProcNo()] = List<vector>(0);
    labelList exSharedPointLabels(0);
    List<List<vector> > exallPs(Pstream::nProcs());
    exallPs[Pstream::myProcNo()] = List<vector>(mesh.boundaryMesh()[filmPatchID].meshPoints().size(), vector::zero);
    forAll(mesh.boundaryMesh()[filmPatchID].meshPoints(), ipoint)
    {
	exallPs[Pstream::myProcNo()][ipoint] = mesh.points()[mesh.boundaryMesh()[filmPatchID].meshPoints()[ipoint]];
    }
    Pstream::gatherList(exallPs);
    Pstream::scatterList(exallPs);
    if (Pstream::master())
    {
        for(label iProc=0; iProc<Pstream::nProcs(); iProc++)
	{
	    forAll(exallPs[iProc], i)
	    {
		for(label j=i+1; j<exallPs[iProc].size(); j++)
		{
		    if (exallPs[iProc][i]==exallPs[iProc][j])
		    {
			bool flag = true;
			forAll(exSharedPoints[Pstream::myProcNo()], k)
			{
			    if (exallPs[iProc][i]==exSharedPoints[Pstream::myProcNo()][k])
			    {
				flag = false;
				break;
			    }
			}
			if (flag)
			{
			    label sz = exSharedPoints[Pstream::myProcNo()].size();
			    exSharedPoints[Pstream::myProcNo()].setSize(sz+1);
			    exSharedPoints[Pstream::myProcNo()][sz] = exallPs[iProc][i];
			}
		    }
		}
		for(label jProc=iProc+1; jProc<Pstream::nProcs(); jProc++)
		forAll(exallPs[jProc], j)
		{
		    if (exallPs[iProc][i]==exallPs[jProc][j])
		    {
			bool flag = true;
			forAll(exSharedPoints[Pstream::myProcNo()], k)
			{
			    if (exallPs[iProc][i]==exSharedPoints[Pstream::myProcNo()][k])
			    {
				flag = false;
				break;
			    }
			}
			if (flag)
			{
			    label sz = exSharedPoints[Pstream::myProcNo()].size();
			    exSharedPoints[Pstream::myProcNo()].setSize(sz+1);
			    exSharedPoints[Pstream::myProcNo()][sz] = exallPs[iProc][i];
			}
		    }
		}
	    }
	}
    }
    Pstream::gatherList(exSharedPoints);
    Pstream::scatterList(exSharedPoints);
    List<List<label> > exspLabels(Pstream::nProcs());
    exspLabels[Pstream::myProcNo()] = List<label>(exSharedPoints[0].size(), -1);
    forAll(exSharedPoints[0], isp)
    {
	point sp = exSharedPoints[0][isp];
	forAll(mesh.boundaryMesh()[filmPatchID].meshPoints(), ispl)
	{
	    const point& lp = mesh.points()[mesh.boundaryMesh()[filmPatchID].meshPoints()[ispl]];
	    if ((fabs(sp.x()-lp.x())<1e-9) &&
	        (fabs(sp.y()-lp.y())<1e-9) &&
		(fabs(sp.z()-lp.z())<1e-9))
	    {
		exspLabels[Pstream::myProcNo()][isp]=ispl;
	    }
	}
    }
    Pstream::gatherList(exspLabels);
    Pstream::scatterList(exspLabels);
    Info<<"exspLabels "<<exspLabels<<endl;