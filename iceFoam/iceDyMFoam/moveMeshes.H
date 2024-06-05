    newPs = (*rmesh).points();
//- set-up interpolator
//	volPointInterpolation vpi(*rmesh);

//- perform interpolation
//	pointScalarField deltaIcep = vpi.interpolate(filmIce->deltaIce()-olddeltaIce);
    volScalarField dDeltaIce=filmIce->deltaIce()-olddeltaIce;
//	Info<<"new value deltaIce:"<<endl;
//	Info<<filmIce->deltaIce()<<endl;
//	Info<<"old value deltaIce:"<<endl;
//	Info<<deltaIcep<<endl;
	const labelList& intfilmPoints=(*rmesh).boundaryMesh()[intfilmPatchID].meshPoints();
	const labelListList& intfilmpointFaces=(*rmesh).boundaryMesh()[intfilmPatchID].pointFaces();
	const vectorField& intfilmnf=(*rmesh).boundary()[intfilmPatchID].Sf();

	forAll(intfilmPoints, index) 
	{
	    label i=intfilmPoints[index];
//	    if (deltaIcep[i]>0.0) {
//	    Info<<"index="<<index<<" i="<<i<<" deltaIce="<<deltaIcep[i]<<endl;
	    vector d(0.0, 0.0, 0.0);
	    label dicount = 0;
	    forAll(intfilmpointFaces[index], faceI) {
//		  Info<<"faceI="<<faceI<<" nf="<<intfilmnf[intfilmpointFaces[index][faceI]]<<" "<<intfilmpointFaces[index][faceI]<<" <="<<intfilmnf.size()<<endl;
		if (dDeltaIce.boundaryField()[intfilmPatchID][intfilmpointFaces[index][faceI]]>0.0) 
		{
//    		    Info<<"i "<<faceI<<" "<<dDeltaIce.boundaryField()[intfilmPatchID][intfilmpointFaces[index][faceI]]<<endl;
		    d+=dDeltaIce.boundaryField()[intfilmPatchID][intfilmpointFaces[index][faceI]]
		    *
		    intfilmnf[intfilmpointFaces[index][faceI]]
		    /
		    mag(intfilmnf[intfilmpointFaces[index][faceI]]);
		}
		dicount++;
	    }
//		Info<<"d="<<d<<endl;
//		Info<<"Point="<<newPs[i]<<endl;
		if (dicount>0)
		    newPs[i]+=d/dicount;
//		newPs[i].x()+=deltaIcep[i]*d.x()/mag(d);
//		newPs[i].z()+=deltaIcep[i]*d.z()/mag(d);
//    	    }
	}
	scalar filmThickness = 0.0;
//	(*rmesh).movePoints(newPs);
//        newPs = (*rmesh).points();

	const labelList& intfilmPointsNew=(*rmesh).boundaryMesh()[intfilmPatchID].meshPoints();
	const labelListList& intfilmpointFacesNew=(*rmesh).boundaryMesh()[intfilmPatchID].pointFaces();
	const vectorField& intfilmnfNew=(*rmesh).boundary()[intfilmPatchID].Sf();

	const labelList& inttopfilmPoints=(*rmesh).boundaryMesh()[inttopfilmPatchID].meshPoints();
//	const labelListList& inttopfilmpointFaces=(*rmesh).boundaryMesh()[inttopfilmPatchID].pointFaces();
//	const vectorField& inttopfilmnf=(*rmesh).boundary()[inttopfilmPatchID].Sf();
//	Pout<<"film "<<filmThickness<<endl;
	forAll(inttopfilmPoints, index)
	{
//	    Pout<<index<<endl;
	    if (index==0)
		filmThickness = mag((*rmesh).points()[intfilmPoints[nearTopFilmPoints[0]]] - (*rmesh).points()[(*rmesh).boundaryMesh()[inttopfilmPatchID].meshPoints()[0]]);
	    label i=inttopfilmPoints[index];
	    label jndex = nearTopFilmPoints[index];
	    label ifilm = intfilmPointsNew[jndex];
//	    Pout<<i<<" "<<jndex<<" "<<ifilm<<endl;
    	    vector d(0.0, 0.0, 0.0);
    	    label dicount = 0;
	    forAll(intfilmpointFacesNew[jndex], faceI) 
	    {
		{
		    d += intfilmnfNew[intfilmpointFacesNew[jndex][faceI]]
		    /
		    mag(intfilmnfNew[intfilmpointFacesNew[jndex][faceI]]);
		}
		dicount++;
	    }
	    if (dicount>0)
		newPs[i] = newPs[ifilm] - d/mag(d)*filmThickness;
//	    Info<<index<<" "<<mag(newPs[i]-(*rmesh).points()[i])<<" "<<newPs[ifilm]<<endl;
	}

	(*rmesh).movePoints(newPs);
	forAll(intfilmPoints, index) 
	{
	    label i=intfilmPoints[index];
	    filmNormal[index]=rmesh->points()[i]-rmesh->points()[topFilmPoints[index]];
	}

//- set-up interpolator
//	primitivePatchInterpolation priInterpolator(rmesh->boundaryMesh()[intfilmPatchID]);
//	scalarField dIcePatch = filmIce->deltaIce().boundaryField()[intfilmPatchID]-olddeltaIce.boundaryField()[intfilmPatchID];
//	scalarField dIcePoint = priInterpolator.faceToPointInterpolate(dIcePatch);
//	const scalarField& dIcePoint = deltaIcep.boundaryField()[intfilmPatchID];

	patchToPatchInterpolation patchInterpolator 
	(
	    rmesh->boundaryMesh()[intfilmPatchID], 
	    mesh.boundaryMesh()[filmPatchID],
	    intersection::FULL_RAY,
	    intersection::CONTACT_SPHERE
	);

//- perform interpolation
	scalarField iceValues = patchInterpolator.faceInterpolate(dDeltaIce.boundaryField()[intfilmPatchID].patchInternalField());
	pointField fromFilmNormal = patchInterpolator.pointInterpolate(filmNormal);
	scalarField iceFaceThickness;
	scalarField icePointThickness;
	scalar maxIceThickness = 0.0;
	if (maxAngleWeight>angleWeight)
	{
	    iceFaceThickness = patchInterpolator.faceInterpolate(filmIce->deltaIce().boundaryField()[intfilmPatchID].patchInternalField());
	    primitivePatchInterpolation priInterpolator(mesh.boundaryMesh()[filmPatchID]);
	    icePointThickness = priInterpolator.faceToPointInterpolate(iceFaceThickness);
	    maxIceThickness = gMax(icePointThickness);
	}
//	Info<<fromFilmNormal<<endl;
//	Info<<dDeltaIce.boundaryField()[intfilmPatchID].patchInternalField()<<endl;
//	Info<<"iceValues"<<endl;
//	Info<<iceValues<<endl;
//	bool moveman = false;

	const labelList& filmPoints=mesh.boundaryMesh()[filmPatchID].meshPoints();
	const labelListList& filmpointFaces=mesh.boundaryMesh()[filmPatchID].pointFaces();
	const vectorField& filmnf=mesh.boundary()[filmPatchID].Sf();

	forAll(filmPoints, index) 
	{
//	    if (iceValues[index]>0.0) {
	    vector d(0.0, 0.0, 0.0);
	    label dicount = 0;
	    forAll(filmpointFaces[index], faceI) 
	    {
		if (iceValues[filmpointFaces[index][faceI]]>0.0)
    		{
//    		    Info<<faceI<<" "<<iceValues[filmpointFaces[index][faceI]]<<endl;
		    d+=iceValues[filmpointFaces[index][faceI]]
		    *
		    filmnf[filmpointFaces[index][faceI]]
		    /
		    mag(filmnf[filmpointFaces[index][faceI]]);
		}
		dicount++;
	    }
	    if (dicount>0)
	    {
	        pDisp[index]-=d/dicount;
	    }
	}

    	// Do any mesh changes
        parcels.storeGlobalPositions();
	mesh.moving(true);
    	mesh.update();
	
	if (angleCrit>SMALL)
	{
	    newPsNF=mesh.points();
	    const vectorField& filmnf2=mesh.boundary()[filmPatchID].Sf();
	    scalar angle = 0.0;
	    scalar angle2 = 0.0;
	    forAll(filmPoints, index) 
	    {
		scalar localAngleWeight = angleWeight;
		if (maxAngleWeight>angleWeight)
		{
		    if (icePointThickness[index]<0.25*maxIceThickness) {
			localAngleWeight *= 0.5;
		    }
		    else if (icePointThickness[index]>maxPartDeltaIce*maxIceThickness)
		    {
			localAngleWeight += (icePointThickness[index]/maxIceThickness-maxPartDeltaIce)
			/
			(1.0-maxPartDeltaIce)
			*
			(maxAngleWeight-angleWeight);
		    }
		}
		label i=filmPoints[index];
		label ip=pointsNearFilm[index];
		vector p(mesh.points()[ip] - mesh.points()[i]);
    		{
		    vector d(0.0, 0.0, 0.0);
		    forAll(filmpointFaces[index], faceI) 
		    {
			d-=filmnf2[filmpointFaces[index][faceI]]/mag(filmnf2[filmpointFaces[index][faceI]]);
		    }
		    bool flagmove = false;
		    angle = Foam::vectorTools::radAngleBetween(p, d, SMALL*SMALL);

		    if (mag(angle)>=angleCrit)
		    {
//			Pout<<index<<" "<<radToDeg(angle)<<endl;
//			Pout<<mesh.points()[ip]<<endl;
//			Pout<<mesh.points()[i]<<endl;
//			if (!moveman) 
//			{
//			moveman = true;
//			}
			vector tempw=d/mag(d)*mag(p);
			tempw = (1.0-localAngleWeight)*p+localAngleWeight*tempw;
//			label ip2=pointsNearFilm2[index];
//			if (mag(tempw)<magWeight*mag(mesh.points()[ip2]-mesh.points()[i]))
//			    tempw *= magWeight*mag(mesh.points()[ip2]-mesh.points()[i])/mag(tempw);
			newPsNF[ip] = mesh.points()[i] + tempw;
			flagmove = true;
		    }
		    label ip2=pointsNearFilm2[index];
		    if (ip2>=0)
		    {
			if (flagmove)
			    p=newPsNF[ip] - mesh.points()[i];
			if (mag(p)<magWeight*mag(mesh.points()[ip2]-mesh.points()[i]))
			{
		    	    p *= magWeight*mag(mesh.points()[ip2]-mesh.points()[i])/mag(p);
		    	    newPsNF[ip] = mesh.points()[i] + p;
//		    	    moveman = true;
			}
		    }
		    label ip3=pointsNearFilm3[index];
		    if (ip3>=0)
		    {
		        vector p2(mesh.points()[ip2] - newPsNF[ip]);
			angle2 = Foam::vectorTools::radAngleBetween(p2, d, SMALL*SMALL);
			if (mag(angle2)>=angleCrit)
			{
//    			    moveman = true;
			    vector tempw2=d/mag(d)*mag(p2);
			    tempw2 = (1.0-0.5*localAngleWeight)*p2+0.5*localAngleWeight*tempw2;
			    newPsNF[ip2] = newPsNF[ip] + tempw2;
			}
		    }
    		}
	    }
//    	    if (moveman)
//    	    {
//    	    }
	    mesh.movePoints(newPsNF);
    	}
	mesh.moving(false);
    	parcels.updateMesh();

	if (reportOrtho>0)
	{
	    iReportOrtho++;
	    if (iReportOrtho>=reportOrtho)
	    {
	         iReportOrtho = 0;
    		faceSet faces(mesh, "nonOrthoFaces", mesh.nFaces()/100+1);
    		mesh.checkFaceOrthogonality(true, &faces);
    		label nFaces = returnReduce(faces.size(), sumOp<label>());
    		if (nFaces > 0)
    		{
        	    Info<<"non Ortho Faces: " <<faces<<endl;
    		}
    	    }
        }
        
