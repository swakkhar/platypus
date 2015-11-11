#include "cbls/idx.hh"
#include "psp/psp.hpp"

#include<fstream>

openKangarooSpace


int connectedPSP()
{
    char Protein [512] =  "HPHHPPHHHHPHHHPPHHPPHPHHHPHPHHPPHHPPPHPPPPPPPPHH" ;
    Idx Length = strlen(Protein);

    cout << "Length: " << Length << endl;
    Itn I = 10;
    Int tabuTenure = 4;
    Int MAX_CONTACT = 10000;

    kblock2<Int> latticePoints(strlen(Protein),3);
    Rnd theRnd(0);
    new_init(Protein, latticePoints, theRnd);

    Hdl handleSystem = Sys::def();
    EvalTi::def(handleSystem);

    HintTi::def(handleSystem);


    kb<Prm> VarX ( Length );   // kb is kangaroo block or a vector
    kb<Prm> VarY ( Length );    // it stores the handles for the vairables declared
    kb<Prm> VarZ ( Length );


    for(Idx tIdx = 0; tIdx < Length ; ++tIdx)
    {
        VarX[tIdx] = RngVi::def(handleSystem, 0, 2*Length);
        VarY[tIdx] = RngVi::def(handleSystem, 0, 2*Length);
        VarZ[tIdx] = RngVi::def(handleSystem, 0, 2*Length);
    }

    PtrPrms(PointsRefX,RngVi,handleSystem, VarX.items(), Length);
    PtrPrms(PointsRefY,RngVi,handleSystem, VarY.items(), Length);
    PtrPrms(PointsRefZ,RngVi,handleSystem, VarZ.items(), Length);



    kb<Prm> XYZ ( Length );

    for(Idx tIdx = 0; tIdx < Length ; ++tIdx)
    {
        XYZ [tIdx] = TpackXiFeVi::def(Xv, handleSystem, VarX[tIdx], VarY[tIdx], VarZ[tIdx]);
    }

    Prm AllDiff = Prm(AllDiffXiFcMiD::def(Xv, handleSystem, XYZ.items(), Length), Dd);


    kb<Prm> tHdl(1);
    tHdl[0] = AllDiff;

    Prm hardSum = Prm(SumXiHFcMiD::def(Xm|Hd,handleSystem,tHdl.items(),tHdl.size()),Hd);


    xb<Prm> DistCons;

    xb<Prm> Dists;
	for(Idx tIdx1 = 0; tIdx1 < Length; ++tIdx1)
	{
		if (Protein[tIdx1] == 'P') continue;
		for(Idx tIdx2 = tIdx1 + 2; tIdx2 < Length; ++tIdx2)
		{
			if (Protein[tIdx2] == 'P') continue;

            Prm pXYZ[3];

			pXYZ[0] = BsubXiFeVi::def(Xv, handleSystem, VarX[tIdx1], VarX[tIdx2]);
            pXYZ[1] = BsubXiFeVi::def(Xv, handleSystem, VarY[tIdx1], VarY[tIdx2]);
            pXYZ[2] = BsubXiFeVi::def(Xv, handleSystem, VarZ[tIdx1], VarZ[tIdx2]);

            Prm tXYZ[3];
            tXYZ[0] = BmultXiFeVi::def(Xv, handleSystem, pXYZ[0], pXYZ[0]);
            tXYZ[1] = BmultXiFeVi::def(Xv, handleSystem, pXYZ[1], pXYZ[1]);
            tXYZ[2] = BmultXiFeVi::def(Xv, handleSystem, pXYZ[2], pXYZ[2]);



			//tXYZ[0] = BdiffXiFeVi::def(Xv, handleSystem, VarX[tIdx1], VarX[tIdx2]);
			//tXYZ[1] = BdiffXiFeVi::def(Xv, handleSystem, VarY[tIdx1], VarY[tIdx2]);
			//tXYZ[2] = BdiffXiFeVi::def(Xv, handleSystem, VarZ[tIdx1], VarZ[tIdx2]);


            DistCons.annex(Prm(SumEquXiBiFeVi::def(Xv, handleSystem, tXYZ, 3, SumEquXiBiFeVi::bind(2)),Lv));

			Dists.annex(SumEqsXiBiFeVi::def(Xv, handleSystem, tXYZ, 3, SumEqsXiBiFeVi::bind(2)));
			//cout << tIdx1 << " " << tIdx2 << " " << Dists.size() << endl;
		}
	}
    // sum all the distances
	Prm DstSum = SumXiFeVi::def(Xv, handleSystem, Dists.items(), Dists.size());

    Prm DistConsSum = SumXiEFcMiD::def(Xv|El, handleSystem, DistCons.items(),DistCons.size());

    // Now we have to model the optimization part

    // we declare a range variable for setting the limit to the objective function
    Prm ObjLimit = RngVi::def(handleSystem,0, MAX_CONTACT);

    // now create another constraint that will limit the DistSum to this limit value

    Prm ObjCons = Prm(BltuXiFcMi::def(Xv, handleSystem, ObjLimit, DstSum), Lm);

    // now define a overall constraint that takes all constraints

    Prm PrmArr[2] = {hardSum, DistConsSum};
	Prm TopSum = SumXiOFcMiD::def(Xm|On, handleSystem, PrmArr, 2);


    // now we define a Tabu heap for the timebeing on the TopSum
    Prm Heap = TabuMaxHeapHiFrHi::def(handleSystem, TopSum);
    Hdl const VarSelcHdl = RankVarSp::def( handleSystem, Heap);
	Hdl const ValSelcHdl = MinValSdXi::def( Xm, handleSystem, TopSum);

    // we need another set of heap and selectors for the hard constraints

    Prm hardSumHeap = TabuMaxHeapHiFrHi::def(handleSystem, hardSum);
    Hdl const hardSumVarSelcHdl = RankVarSp::def( handleSystem, hardSumHeap);
	Hdl const hardSumValSelcHdl = MinValSdXi::def( Xm, handleSystem, hardSum);


    RefHdl(RankVarSpHard,RankVarSp,handleSystem,hardSumVarSelcHdl);

    Sys const & ss = Sys::ref(handleSystem);
    RefHdl(TopSumRef, SumXiOFcMiD, handleSystem, TopSum.TermHdl);
	RefPrm(DstSumRef, SumXiFeVi, handleSystem, DstSum);
	RefPrm(HardSumRef, SumXiHFcMiD, handleSystem, hardSum);
	//RefPrm(AllDiffRef, AllDiffXiFcMiD, handleSystem, AllDiff);


    Sys::setVarTabuLimit(handleSystem, tabuTenure);//theRnd(4,(Int)Length/4));

    kb<Itr> ItrVals(Length*3+1);
    for(Idx tIdx=0;tIdx < ItrVals.size()-1;tIdx+=3)
    {
        ItrVals[tIdx] =   latticePoints[tIdx/3][0];
        ItrVals[tIdx+1] = latticePoints[tIdx/3][1];
        ItrVals[tIdx+2] = latticePoints[tIdx/3][2];
        cout << ItrVals[tIdx] << " " << ItrVals[tIdx+1] << " " << ItrVals[tIdx+2] << endl;
    }
    ItrVals[Length*3] = 0;
//    Sys::execAnewRndVars(handleSystem,theRnd);
    Sys::execAnewItrVals(handleSystem,ItrVals.items());

    N sum=0;
    for(Idx tIdx1 = 0; tIdx1 < Length; ++tIdx1)
	{
		if (Protein[tIdx1] == 'P') continue;
		for(Idx tIdx2 = tIdx1 + 2; tIdx2 < Length; ++tIdx2)
		{
			if (Protein[tIdx2] == 'P') continue;



            N tempX= (PointsRefX[tIdx1]->CurrItr()-PointsRefX[tIdx2]->CurrItr())*(PointsRefX[tIdx1]->CurrItr()-PointsRefX[tIdx2]->CurrItr());
            N tempY = (PointsRefY[tIdx1]->CurrItr()-PointsRefY[tIdx2]->CurrItr())*(PointsRefY[tIdx1]->CurrItr()-PointsRefY[tIdx2]->CurrItr());
            N tempZ = (PointsRefZ[tIdx1]->CurrItr()-PointsRefZ[tIdx2]->CurrItr())*(PointsRefZ[tIdx1]->CurrItr()-PointsRefZ[tIdx2]->CurrItr());

            if(tempX+tempY+tempZ == 2)
            {
                ++sum;
            }
		}
	}

	cout << "Initial Sum :" << sum <<endl;

    Int globalbest = sum;
    Int hardThreshold = 0;
    Int countFeasible = 0;


    cout << "Initial Value of Obj:" << DstSumRef.ValueRec().CurrData() <<endl;


    // take values from the input file and initialize with them


    Sys::execDiffItrVal(handleSystem, ObjLimit.TermHdl,0);
    Sys::setVarLocked(handleSystem, ObjLimit.TermHdl, true);
    Sys::setVarTabu(handleSystem, ObjLimit.TermHdl,true);

    Sys::setVarLocked(handleSystem, VarX[0].TermHdl, true);
    Sys::setVarTabu(handleSystem, VarX[0].TermHdl,true);

    Sys::setVarLocked(handleSystem, VarY[0].TermHdl, true);
    Sys::setVarTabu(handleSystem, VarY[0].TermHdl,true);

    Sys::setVarLocked(handleSystem, VarZ[0].TermHdl, true);
    Sys::setVarTabu(handleSystem, VarZ[0].TermHdl,true);

    // the objlimit value should always be locked !!!!!

    cout << "Initial Top Sum:" << TopSumRef.MetricRec().CurrData() << endl;
    cout << "Initial DistSum Sum:" << DstSumRef.ValueRec().CurrData() << endl;
    cout << "Initial Hard Sum:" << HardSumRef.MetricRec().CurrData() << endl;

    while(ss.ExecClk() <= I) // only iteration
    {

        if( HardSumRef.MetricRec().CurrData() <= hardThreshold ) // top
        {
            cout << " Clk: " << ss.ExecClk() << " Upper: " << HardSumRef.MetricRec().CurrData()<< endl;
            Sys::doSelcPairExecDiff(handleSystem, VarSelcHdl, ValSelcHdl, theRnd);

            RefPrm(TopHeapRef,TabuMaxHeapHiFrHi,handleSystem,Heap);

            Sln tSln =  TopHeapRef.begin();
            cout << "Variables: ";
            while(tSln != TopHeapRef.end())
            {
                cout << tSln << " ";
                tSln=TopHeapRef.move(tSln);
            }
            cout << endl;



        }
        else
        {
            Sys::doSelcPairExecDiff(handleSystem, hardSumVarSelcHdl, hardSumValSelcHdl, theRnd);

            cout << " Clk: " << ss.ExecClk() << " Lower: " << HardSumRef.MetricRec().CurrData()<< endl;
            if(HardSumRef.MetricRec().CurrData() == 0)
            {
                ++countFeasible;
                cout << "obj: " << DstSumRef.ValueRec().CurrData() << endl;
            }
            if(HardSumRef.MetricRec().CurrData() == 0 && DstSumRef.ValueRec().CurrData() > globalbest)
            {
                globalbest =  DstSumRef.ValueRec().CurrData();

                Sys::setVarLocked(handleSystem, ObjLimit.TermHdl, false);
                Sys::execDiffItrVal(handleSystem, ObjLimit.TermHdl,globalbest);
                Sys::setVarLocked(handleSystem, ObjLimit.TermHdl, true);

                cout << "obj locked to: " << globalbest << endl;
                writeHPToFile("data.cml", Protein, PointsRefX, PointsRefY, PointsRefZ);
            }
        }

    }
    cerr << ss.ExecClk() << " Global Best: "<< globalbest << " No Of feasible Solutions: " << countFeasible << endl;


    return 0;
}

//int proteinSP()
//{
//
//    #define MAX_CONTACT 1000   // This is the upperbound set to the objective function, later we can change this
//
//    // start of code
//    char Protein [ 512 ] =  "HPHHPPHHHHPHHHPPHHPPHPHHHPHPHHPPHHPPPHPPPPPPPPHH" ;
////    char Protein [ 512 ] =  "HPHHPPHHHH" ;
//
//    N Length = strlen ( Protein );
//
//    Itr** latticePoints= new Itr*[Length];
//	for(Idx _i=0;_i<Length;_i++)
//		latticePoints[_i]=new Itr[3]; // 3 is the current dimension
//    inputFromFile("psp//input.txt", latticePoints,Length);
//
//    cout << "Protein: " << Protein << " Length: " << strlen ( Protein ) << endl;
//
//    // now we should define a system
//    Hdl handleSystem = Sys::def(); // returns a handle to the system
//    // we are considering integer variables only
//    EvalTi::def(handleSystem);
//
//    HintTi::def(handleSystem); // need to declare this since we need the hints
//
//    // now we have to declare the variables, we need protein length
//    // for each amino acids we need three variables X Y ans Z
//    kb<Prm> VarX ( Length );   // kb is kangaroo block or a vector
//    kb<Prm> VarY ( Length );    // it stores the handles for the vairables declared
//    kb<Prm> VarZ ( Length );
//
//
//
//    // take reference to all these variables
//
//
//    // these three arrays will contain the lattice point co-ordinates for all amino acid monomers
//
//    // at first we fix the first point to be a large positive constant equal to the lenth of the protein
//
//    // we have to set this value to VarX[0], VarY[0] and VarZ[0]
//    // we get handles for these three variables for a constant value
//    //VarX[0] = VarY[0] = VarZ[0] = Ci::def( handleSystem, Length );
//
//
//    // we have to specifiy domains for all these variables
//
//
//
//    for(Idx tIdx = 0; tIdx < Length ; ++tIdx)
//    {
//            VarX[tIdx] = RngVi::def(handleSystem, 0, 2*Length);
//            VarY[tIdx] = RngVi::def(handleSystem, 0, 2*Length);
//            VarZ[tIdx] = RngVi::def(handleSystem, 0, 2*Length);
//    }
//
//    PtrPrms(PointsRefX,RngVi,handleSystem, VarX.items(), Length);
//    PtrPrms(PointsRefY,RngVi,handleSystem, VarY.items(), Length);
//    PtrPrms(PointsRefZ,RngVi,handleSystem, VarZ.items(), Length);
//
//
//
//    // now we have to pack all the point triplets
//    //we define a new set of variables XYZ to indicate triples
//
//    kb<Prm> XYZ ( Length );
//
//    // we have to pack all the variables in to these XYZ triples
//
//    for(Idx tIdx = 0; tIdx < Length ; ++tIdx)
//    {
//        XYZ [tIdx] = TpackXiFeVi::def(Xv, handleSystem, VarX[tIdx], VarY[tIdx], VarZ[tIdx]);
//    }
//
//    // we must define all different for all these triples
//
//    Prm AllDiff = Prm(AllDiffXiFcMiD::def(Xv, handleSystem, XYZ.items(), Length), Dd);
//
//    // now we have to declare the FCC constraint, X+Y+Z mod 2 = 0 for all the points
//    kb<Prm> FccCons (Length);
//    for(Idx tIdx = 0; tIdx < Length ; ++tIdx)
//    {
//        // for each triplet create a new triplet
//        kb<Prm> Triplet(3); // size fixed
//        Triplet[0] = VarX[tIdx];
//        Triplet[1] = VarY[tIdx];
//        Triplet[2] = VarZ[tIdx];
//
//        //now define the fcc constraint
//        FccCons[tIdx] = Prm(SumModXiBiFeVi::def(Xv,handleSystem,Triplet.items(),3,SumModXiBiFeVi::bind(2)), Lv);
//        // probably we need a ground variable version
//    }
//
//    PtrPrms(FccConsRef,SumModXiBiFeVi,handleSystem,FccCons.items(),Length);
//
//    //we need a sum of all these fcc constraints
//    Prm SumFccCons  = SumXiEFcMiD::def(Xv|El,handleSystem,FccCons.items(),Length);
//
//
//
//    // now we need connectivity constraints
//    // for each of the consecutive amino acids the sum of the abs diff is exactly 2
//
//    Prm conConstant = Ci::def(handleSystem,1);
//    kb <Prm> extraConstraint(3*(Length - 1));
//
//    kb<Prm> ConnCons(Length - 1);
////    cout << "here" << endl;
//    for(Idx tIdx = 1; tIdx < Length ; tIdx++ )
//    {
//
//        kb<Prm> pXYZ(3);
//
//        pXYZ[0] = BsubXiFeVi::def(Xv, handleSystem, VarX[tIdx], VarX[tIdx - 1]);
//        pXYZ[1] = BsubXiFeVi::def(Xv, handleSystem, VarY[tIdx], VarY[tIdx - 1]);
//        pXYZ[2] = BsubXiFeVi::def(Xv, handleSystem, VarZ[tIdx], VarZ[tIdx - 1]);
//
//
//        kb<Prm> tXYZ(3);
//       // cout << " param " << pXYZ[0].TermHdl << endl;
//        tXYZ[0] = BmultXiFeVi::def(Xv, handleSystem, pXYZ[0], pXYZ[0]);
//       // cout << " param " << pXYZ[0].TermHdl << endl;
//        tXYZ[1] = BmultXiFeVi::def(Xv, handleSystem, pXYZ[1], pXYZ[1]);
//        tXYZ[2] = BmultXiFeVi::def(Xv, handleSystem, pXYZ[2], pXYZ[2]);
//
//        //tXYZ[0] = BdiffXiFeVi::def(Xv, handleSystem, VarX[tIdx], VarX[tIdx - 1]);
//		//tXYZ[1] = BdiffXiFeVi::def(Xv, handleSystem, VarY[tIdx], VarY[tIdx - 1]);
//		//tXYZ[2] = BdiffXiFeVi::def(Xv, handleSystem, VarZ[tIdx], VarZ[tIdx - 1]);
//
//        //kb<Prm> sXYZ(3);
//
//        extraConstraint[3*(tIdx-1)+0] = Prm (BleuXiFcMi::def(Xv, handleSystem, tXYZ[0],conConstant),Lm);
//        extraConstraint[3*(tIdx-1)+1] = Prm (BleuXiFcMi::def(Xv, handleSystem, tXYZ[1],conConstant),Lm);
//        extraConstraint[3*(tIdx-1)+2] = Prm (BleuXiFcMi::def(Xv, handleSystem, tXYZ[2],conConstant),Lm);
//
//        //extraConstraint[tIdx -1 ] = SumXiEFcMiD::def(Xv|El,handleSystem, tXYZ.items(),3);
//
//        ConnCons[tIdx - 1] = Prm(SumNesXiBiFeVi::def(Xv, handleSystem, tXYZ.items(), 3, SumNesXiBiFeVi::bind(2)),Lv);
////        ConnCons[tIdx - 1] = Prm(SumEquXiBiFeVi::def(Xv, handleSystem, tXYZ.items(), 3, SumEquXiBiFeVi::bind(2)),Lv);
//    }
////    cout << "here" << endl;
//    Prm SumConnCons  = SumXiEFcMiD::def(Xv|El,handleSystem,ConnCons.items(),Length-1);
//
//    Prm SumExtra = SumXiEFcMiD::def(Xm|El,handleSystem,extraConstraint.items(),extraConstraint.size());
//
//
//    kb<Prm> tHdl(3);
//    tHdl[0] = AllDiff;
//    //tHdl[1] = SumFccCons;
//    tHdl[1] = SumConnCons;
//    tHdl[2] = SumExtra;
//
//    Prm hardSum = Prm(SumXiHFcMiD::def(Xm|Hd,handleSystem,tHdl.items(),tHdl.size()),Hd);
//
//
//
//    //for debugging purpose we declare some Ref to the internal nodes
//
//    // now we have to add the objective function here
//    xb<Prm> DistCons;
//
//    xb<Prm> Dists;
//	for(Idx tIdx1 = 0; tIdx1 < Length; ++tIdx1)
//	{
//		if (Protein[tIdx1] == 'P') continue;
//		for(Idx tIdx2 = tIdx1 + 2; tIdx2 < Length; ++tIdx2)
//		{
//			if (Protein[tIdx2] == 'P') continue;
//
//            Prm pXYZ[3];
//
//			pXYZ[0] = BsubXiFeVi::def(Xv, handleSystem, VarX[tIdx1], VarX[tIdx2]);
//            pXYZ[1] = BsubXiFeVi::def(Xv, handleSystem, VarY[tIdx1], VarY[tIdx2]);
//            pXYZ[2] = BsubXiFeVi::def(Xv, handleSystem, VarZ[tIdx1], VarZ[tIdx2]);
//
//            Prm tXYZ[3];
//            tXYZ[0] = BmultXiFeVi::def(Xv, handleSystem, pXYZ[0], pXYZ[0]);
//            tXYZ[1] = BmultXiFeVi::def(Xv, handleSystem, pXYZ[1], pXYZ[1]);
//            tXYZ[2] = BmultXiFeVi::def(Xv, handleSystem, pXYZ[2], pXYZ[2]);
//
//
//
//			//tXYZ[0] = BdiffXiFeVi::def(Xv, handleSystem, VarX[tIdx1], VarX[tIdx2]);
//			//tXYZ[1] = BdiffXiFeVi::def(Xv, handleSystem, VarY[tIdx1], VarY[tIdx2]);
//			//tXYZ[2] = BdiffXiFeVi::def(Xv, handleSystem, VarZ[tIdx1], VarZ[tIdx2]);
//
//
//            DistCons.annex(Prm(SumEquXiBiFeVi::def(Xv, handleSystem, tXYZ, 3, SumEquXiBiFeVi::bind(2)),Lv));
//
//			Dists.annex(SumEqsXiBiFeVi::def(Xv, handleSystem, tXYZ, 3, SumEqsXiBiFeVi::bind(2)));
//			//cout << tIdx1 << " " << tIdx2 << " " << Dists.size() << endl;
//		}
//	}
//    // sum all the distances
//	Prm DstSum = SumXiFeVi::def(Xv, handleSystem, Dists.items(), Dists.size());
//
//    Prm DistConsSum = SumXiEFcMiD::def(Xv|El, handleSystem, DistCons.items(),DistCons.size());
//
//    // Now we have to model the optimization part
//
//    // we declare a range variable for setting the limit to the objective function
//    Prm ObjLimit = RngVi::def(handleSystem,0, MAX_CONTACT);
//
//    // now create another constraint that will limit the DistSum to this limit value
//
//    Prm ObjCons = Prm(BltuXiFcMi::def(Xv, handleSystem, ObjLimit, DstSum), Lm);
//
//    // now define a overall constraint that takes all constraints
//
//    Prm PrmArr[2] = {hardSum, DistConsSum};
//	Prm TopSum = SumXiOFcMiD::def(Xm|On, handleSystem, PrmArr, 2);
//
//
//    // now we define a Tabu heap for the timebeing on the TopSum
//    Prm Heap = TabuMaxHeapHiFrHi::def(handleSystem, TopSum);
//    Hdl const VarSelcHdl = RankVarSp::def( handleSystem, Heap);
//	Hdl const ValSelcHdl = MinValSdXi::def( Xm, handleSystem, TopSum);
//
//    // we need another set of heap and selectors for the hard constraints
//
//    Prm hardSumHeap = TabuMaxHeapHiFrHi::def(handleSystem, hardSum);
//    Hdl const hardSumVarSelcHdl = RankVarSp::def( handleSystem, hardSumHeap);
////	Hdl const hardSumValSelcHdl = MinValSdXi::def( Xm, handleSystem, hardSum);
//
//    Sys const & ss = Sys::ref(handleSystem);
//    RefHdl(TopSumRef, SumXiOFcMiD, handleSystem, TopSum.TermHdl);
//	RefPrm(DstSumRef, SumXiFeVi, handleSystem, DstSum);
//	RefPrm(HardSumRef, SumXiHFcMiD, handleSystem, hardSum);
//	//RefPrm(AllDiffRef, AllDiffXiFcMiD, handleSystem, AllDiff);
//
//    // now we can randomly initialize the values
//    Rnd theRnd;
//    cout << "The Random Seed:" << theRnd.Seed() << endl;
//    Sys::setVarTabuLimit(handleSystem, 3);//theRnd(4,(Int)Length/4));
//
//    kb<Itr> ItrVals(Length*3);
//    for(Idx tIdx=0;tIdx < ItrVals.size();tIdx+=3)
//    {
//        ItrVals[tIdx] =   latticePoints[tIdx/3][0];
//        ItrVals[tIdx+1] = latticePoints[tIdx/3][1];
//        ItrVals[tIdx+2] = latticePoints[tIdx/3][2];
//        //cout << ItrVals[tIdx] << " " << ItrVals[tIdx+1] << " " << ItrVals[tIdx+2] << endl;
//    }
//
//
//    //Sys::execAnewItrVals(handleSystem,ItrVals.items());
//
//    N sum=0;
//    for(Idx tIdx1 = 0; tIdx1 < Length; ++tIdx1)
//	{
//		if (Protein[tIdx1] == 'P') continue;
//		for(Idx tIdx2 = tIdx1 + 2; tIdx2 < Length; ++tIdx2)
//		{
//			if (Protein[tIdx2] == 'P') continue;
//
//
//
//            N tempX= (PointsRefX[tIdx1]->CurrItr()-PointsRefX[tIdx2]->CurrItr())*(PointsRefX[tIdx1]->CurrItr()-PointsRefX[tIdx2]->CurrItr());
//            N tempY = (PointsRefY[tIdx1]->CurrItr()-PointsRefY[tIdx2]->CurrItr())*(PointsRefY[tIdx1]->CurrItr()-PointsRefY[tIdx2]->CurrItr());
//            N tempZ = (PointsRefZ[tIdx1]->CurrItr()-PointsRefZ[tIdx2]->CurrItr())*(PointsRefZ[tIdx1]->CurrItr()-PointsRefZ[tIdx2]->CurrItr());
//
//            if(tempX+tempY+tempZ == 2)
//            {
//                ++sum;
//            }
//		}
//	}
//
////	cout << "Initial Sum :" << sum <<endl;
//
//
//
//
////    cout << "Initial Value of Obj:" << DstSumRef.ValueRec().CurrData() <<endl;
//    Sys::execAnewRndVars(handleSystem, theRnd);
//
//	Sys::execDiffItrVal(handleSystem,0,Length);
//	Sys::execDiffItrVal(handleSystem,1,Length);
//	Sys::execDiffItrVal(handleSystem,2,Length);
//    for(Idx tIdx= 0;tIdx < Length*3 - 3;tIdx+=3)
//    {
//        Int tempX = PointsRefX[tIdx/3]->ValueRec().CurrData();
//        Int tempY = PointsRefY[tIdx/3]->ValueRec().CurrData();
//        Int tempZ = PointsRefZ[tIdx/3]->ValueRec().CurrData();
//
////        Int tempX = Vi::ref(handleSystem, tIdx/3).ValueRec().CurrData();
////        Int tempY = Vi::ref(handleSystem, tIdx/3 + 1).ValueRec().CurrData();
////        Int tempZ = Vi::ref(handleSystem, tIdx/3 + 2).ValueRec().CurrData();
//		cerr << tempX << " " << tempY << " " << tempZ << endl;
//
//        Sys::execDiffItrVal(handleSystem,tIdx +3,tempX + uniform(theRnd,castZ(-1),castZ(1)));
//        Sys::execDiffItrVal(handleSystem,tIdx + 4,tempY + uniform(theRnd,castZ(-1),castZ(1)));
//        Sys::execDiffItrVal(handleSystem,tIdx + 5,tempZ + uniform(theRnd,castZ(-1),castZ(1)));
//    }
// //   Sys::execDiffItrVal(handleSystem,0,Length);
////    Sys::execDiffItrVal(handleSystem,1,Length);
////    Sys::execDiffItrVal(handleSystem,2,Length);
//
//    //writeHPToFile("data.cml", Protein, PointsRefX, PointsRefY, PointsRefZ);
//    //exit(0);
//
//
//
//    // take values from the input file and initialize with them
//
//
//    Sys::execDiffItrVal(handleSystem, ObjLimit.TermHdl,0);
//    Sys::setVarLocked(handleSystem, ObjLimit.TermHdl, true);
//    Sys::setVarTabu(handleSystem, ObjLimit.TermHdl,true);
//
//    //Sys::setVarLocked(handleSystem, VarX[0].TermHdl, true);
//    //Sys::setVarTabu(handleSystem, VarX[0].TermHdl,true);
//
//    //Sys::setVarLocked(handleSystem, VarY[0].TermHdl, true);
//    //Sys::setVarTabu(handleSystem, VarY[0].TermHdl,true);
//
//    //Sys::setVarLocked(handleSystem, VarZ[0].TermHdl, true);
//    //Sys::setVarTabu(handleSystem, VarZ[0].TermHdl,true);
//
//    // the objlimit value should always be locked !!!!!
//
//
//    Int hardThreshold = 0;
//
//    Int globalbest=0; // if the solution is feasible at initialization
//
//    //since all the propagators are set we can start the algorithm
//    Itn I = 1000000;
//
//    Int countFeasible =0;
//    //Bln MoveTypeHard = true;
//    while(ss.ExecClk() <= I) // only iteration
//    {
//        // check the packs
//        Int tFccSum = 0;
////        Int tConnSum = 0;
//            int tDisIdx=0;
//
//        Int tObjSum = 0;
//        Int tAllDiffCount = 0;
//        xset<Int>  diffDict(Length);
//        for(Idx tIdx = 0; tIdx < Length ; ++tIdx)
//        {
//            for(Idx tIdx2 = tIdx + 2; tIdx2 < Length; ++tIdx2)
//            {
//                if(Protein[tIdx]=='P'||Protein[tIdx2]=='P') continue;
//
//                //
//                Int tempX = PointsRefX[tIdx]->ValueRec().CurrData()-PointsRefX[tIdx2]->ValueRec().CurrData();
//                Int tempY = PointsRefY[tIdx]->ValueRec().CurrData()-PointsRefY[tIdx2]->ValueRec().CurrData();
//                Int tempZ = PointsRefZ[tIdx]->ValueRec().CurrData()-PointsRefZ[tIdx2]->ValueRec().CurrData();
//                Int sum=tempX*tempX+tempY*tempY+tempZ*tempZ;
//                if(SumEqsXiBiFeVi::ref(handleSystem,Dists[tDisIdx].TermHdl).ValueRec().CurrData() != (sum==2))
//                    cerr << "Error in Dist! "  <<
//                    SumEqsXiBiFeVi::ref(handleSystem,Dists[tDisIdx].TermHdl).ValueRec().CurrData()
//                    << " "  << sum <<  " " << (sum == 2) << endl;
//
//                if(sum==2)
//                {
//                    //cerr << " "  <<
//                    //SumEqsXiBiFeVi::ref(handleSystem,Dists[tDisIdx].TermHdl).ValueRec().CurrData()
//                    //<< " "  << sum <<  " " << (sum == 2) << endl;
//                    tObjSum++;
//                }
//                tDisIdx ++;
//                //cout << tIdx << " " << tIdx2 << " " << tDisIdx<< endl;
//
//            }
//
////            if(tIdx > 0 && SumEquXiBiFeVi::ref(handleSystem,ConnCons[tIdx-1].TermHdl).ValueRec().CurrData() !=
////
////               (bequ<Int,Int>::ifo(
////
////                usqr<Int,Int>::ifo(PointsRefX[tIdx]->ValueRec().CurrData()-PointsRefX[tIdx-1]->ValueRec().CurrData())+
////                usqr<Int,Int>::ifo(PointsRefY[tIdx]->ValueRec().CurrData()-PointsRefY[tIdx-1]->ValueRec().CurrData())+
////                usqr<Int,Int>::ifo(PointsRefZ[tIdx]->ValueRec().CurrData()-PointsRefZ[tIdx-1]->ValueRec().CurrData())
////
////               ,2)
////               ))
////
////                cerr << "Error in Conn !" << endl;
////                if(tIdx > 0)
////                tConnSum+=(bequ<Int,Int>::ifo(
////
////                usqr<Int,Int>::ifo(PointsRefX[tIdx]->ValueRec().CurrData()-PointsRefX[tIdx-1]->ValueRec().CurrData())+
////                usqr<Int,Int>::ifo(PointsRefY[tIdx]->ValueRec().CurrData()-PointsRefY[tIdx-1]->ValueRec().CurrData())+
////                usqr<Int,Int>::ifo(PointsRefZ[tIdx]->ValueRec().CurrData()-PointsRefZ[tIdx-1]->ValueRec().CurrData())
////
////               ,2));
//
////
////            if(tIdx > 0 && SumEquXiBiFeVi::ref(handleSystem,ConnCons[tIdx-1].TermHdl).ValueRec().CurrData() !=
////
////                (bequ<Int,Int>::ifo(
////                    ((usqr<Int,Int>::ifo(PointsRefX[tIdx]->ValueRec().CurrData()-PointsRefX[tIdx-1]->ValueRec().CurrData())+
////               usqr<Int,Int>::ifo(PointsRefY[tIdx]->ValueRec().CurrData()-PointsRefY[tIdx-1]->ValueRec().CurrData())+
////                usqr<Int,Int>::ifo(PointsRefZ[tIdx]->ValueRec().CurrData()-PointsRefZ[tIdx-1]->ValueRec().CurrData())),2))))
////                cerr << "Error in Conn !" << endl;
////            if(tIdx > 0)
////            tConnSum+=(bequ<Int,Int>::ifo(((usqr<Int,Int>::ifo(PointsRefX[tIdx]->ValueRec().CurrData()-PointsRefX[tIdx-1]->ValueRec().CurrData())+
////               usqr<Int,Int>::ifo(PointsRefY[tIdx]->ValueRec().CurrData()-PointsRefY[tIdx-1]->ValueRec().CurrData())+
////                usqr<Int,Int>::ifo(PointsRefZ[tIdx]->ValueRec().CurrData()-PointsRefZ[tIdx-1]->ValueRec().CurrData())),2)));
//
//            if((FeVi::ref(handleSystem, XYZ[tIdx].TermHdl).ValueRec().CurrData() !=
//            ipack3<Int,Int> ::ifo(PointsRefX[tIdx]->ValueRec().CurrData() ,
//                                  PointsRefY[tIdx]->ValueRec().CurrData(),
//                                  PointsRefZ[tIdx]->ValueRec().CurrData())))
//                cerr << "Error in pack !" << endl;
//
//
//            if(diffDict.append(ipack3<Int,Int> ::ifo(PointsRefX[tIdx]->ValueRec().CurrData() ,
//                                  PointsRefY[tIdx]->ValueRec().CurrData(),
//                                  PointsRefZ[tIdx]->ValueRec().CurrData()))==false)
//                                  tAllDiffCount++;
//
//            if(FccConsRef[tIdx]->ValueRec().CurrData()!=
//               ((PointsRefX[tIdx]->ValueRec().CurrData()+
//                PointsRefY[tIdx]->ValueRec().CurrData()+
//                PointsRefZ[tIdx]->ValueRec().CurrData()
//                )%2))
//                cerr << "Error in Fcc Constraint! " << endl;
//            tFccSum+=((PointsRefX[tIdx]->ValueRec().CurrData()+
//                PointsRefY[tIdx]->ValueRec().CurrData()+
//                PointsRefZ[tIdx]->ValueRec().CurrData()
//                )%2);
//
//        }
//
//        if(tAllDiffCount != AllDiffXiFcMiD::ref(handleSystem,AllDiff.TermHdl).MetricRec().CurrData())
//            cerr << "Error in AllDiff Count!" << endl;
//        if(DstSumRef.ValueRec().CurrData() != tObjSum)
//                cerr << "Error in Dist Sum! " << tObjSum << " " << DstSumRef.ValueRec().CurrData() <<endl;
//
//        if(SumXiEFcMiD::ref(handleSystem,SumFccCons.TermHdl).MetricRec().CurrData() !=
//           tFccSum)
//           cerr << "Error in Fcc Sum!" << endl;
////        if(SumXiEFcMiD::ref(handleSystem,SumConnCons.TermHdl).MetricRec().CurrData() !=
////           tConnSum )
////           cerr << "Error in Conc Sum" << endl;
////        //cout << tFccSum << endl ;
//
//
//        //if(HardSumRef.MetricRec().CurrData()!=tConnSum+tAllDiffCount+tFccSum)
//          //  cerr << "Error in hardsum!" << endl;
//
//
//
//
//        //if(  BltuXiFcMi::ref(handleSystem,ObjCons.TermHdl).MetricRec().CurrData()!=
//
//            //1+(RngVi::ref(handleSystem,ObjLimit.TermHdl).ValueRec().CurrData())-tObjSum)
//           //cerr << "Error in Obj Constraint!" << " " <<
//           //BltuXiFcMi::ref(handleSystem,ObjCons.TermHdl).MetricRec().CurrData() << " "<<
//           //abs(tObjSum - RngVi::ref(handleSystem,ObjLimit.TermHdl).ValueRec().CurrData()) << endl;
//
//        //Sys::execDiffRndVal(handleSystem,theRnd(3u,3*Length-1),theRnd);
//        cout << HardSumRef.MetricRec().CurrData() << " " <<
//        TopSumRef.MetricRec().CurrData() <<
//        " " << AllDiffXiFcMiD::ref(handleSystem,AllDiff.TermHdl).MetricRec().CurrData() <<
//        " " << SumXiEFcMiD::ref(handleSystem,SumConnCons.TermHdl).MetricRec().CurrData() <<
//        " " << SumXiEFcMiD::ref(handleSystem,SumFccCons.TermHdl).MetricRec().CurrData() <<
//        " " << SumXiEFcMiD::ref(handleSystem, DistConsSum.TermHdl).MetricRec().CurrData() << " >> ";
//        if( HardSumRef.MetricRec().CurrData() <= hardThreshold ) // top
//        {
//            Sys::doSelcPairExecDiff(handleSystem, VarSelcHdl, ValSelcHdl, theRnd);
//        }
//        else // hard
//        {
//<<<<<<< Updated upstream
///*			Sys::runSelc(handleSystem, hardSumVarSelcHdl, theRnd);
//            RefHdl(tVarSelc, RankVarSp, handleSystem, hardSumVarSelcHdl);
//            Hdl tVarHdl = tVarSelc.getSelcVars()[0];
////			cout << "selected " << tVarHdl << " ";
//
//			xb<Int> sVals;
////            xb<Itr> tVals;
//			tVarHdl = tVarHdl / 3 * 3;
////			if (tVarHdl < 3)
////			{
////				tVals.annex(PointsRefX[1]->ValueRec().CurrData() - 1);
////				tVals.annex(PointsRefX[1]->ValueRec().CurrData());
////				tVals.annex(PointsRefX[1]->ValueRec().CurrData() + 1);
////			}
////			else if (tVarHdl >= 3 * Length - 3)
////			{
////				tVals.annex(PointsRefX[Length - 2]->ValueRec().CurrData() - 1);
////				tVals.annex(PointsRefX[Length - 2]->ValueRec().CurrData());
////				tVals.annex(PointsRefX[Length - 2]->ValueRec().CurrData() + 1);
////			}
////			else
////			{
////				Int tTmp = PointsRefX[tVarHdl/3 - 1]->ValueRec().CurrData() + PointsRefX[tVarHdl/3 + 1]->ValueRec().CurrData();
////				tTmp /= 2;
////				tVals.annex(tTmp - 1);
////				tVals.annex(tTmp);
////				tVals.annex(tTmp + 1);
////			}
//			Idx minIdx = InvIdx;
//            Int minVal= MaxInt;
//            for(Idx tIdx = 0; tIdx < 2 * Length; ++tIdx )
////            for(Idx tIdx = 0; tIdx < tVals.size(); ++tIdx )
//            {
////                if(tVals.item(tIdx)>=0 && tVals[tIdx] <= 2*Length) // within the range
////                {
//                    Sys::simulDiffItrVal(handleSystem,tVarHdl,tIdx);
//                    Sys::simulIncr(handleSystem, hardSum.TermHdl);
//                    Int thisVal = HardSumRef.MetricRec().NextData(ss.SimulClk());
//                    if(thisVal < minVal)
//                    {
//                        minVal=thisVal;
//						minIdx= tIdx;
//						sVals.reset();
//						sVals.annex(tIdx);
//                    }
//					else if (thisVal == minVal)
//						sVals.annex(tIdx);
////                }
//=======
//			Sys::runSelc(handleSystem, hardSumVarSelcHdl, theRnd);
//            RefHdl(vaRSelc, RankVarSp, handleSystem, hardSumVarSelcHdl);
////            Selc * selector =(Selc *)(&vaRSelc);
////            selector->select(theRnd);
//          cout<< "Size " << vaRSelc.getSelcVars().size()<< " var:" << vaRSelc.getSelcVars()[0];
//            Idx varIdx=vaRSelc.getSelcVars()[0];
//			cout << " val: " << RngVi::ref(handleSystem, varIdx).ValueRec().CurrData() << endl;
//            xb<Itr> vals;
//            if(varIdx-3>=0)
//            {
//                NewPtrHdl(tmp1,RngVi,handleSystem,varIdx-3);
//                vals.annex(tmp1->ValueRec().CurrData()-1);
//                vals.annex(tmp1->ValueRec().CurrData());
//                vals.annex(tmp1->ValueRec().CurrData()+1);
//>>>>>>> Stashed changes
//            }
//			if (minIdx != InvIdx)
//				Sys::execDiffItrVal(handleSystem,tVarHdl,sVals[theRnd(sVals.size())]);
//
//			++tVarHdl;
////////			tVals.reset();
////////			if (tVarHdl < 3)
////////			{
////////				tVals.annex(PointsRefY[1]->ValueRec().CurrData() - 1);
////////				tVals.annex(PointsRefY[1]->ValueRec().CurrData());
////////				tVals.annex(PointsRefY[1]->ValueRec().CurrData() + 1);
////////			}
////////			else if (tVarHdl >= 3 * Length - 3)
////////			{
////////				tVals.annex(PointsRefY[Length - 2]->ValueRec().CurrData() - 1);
////////				tVals.annex(PointsRefY[Length - 2]->ValueRec().CurrData());
////////				tVals.annex(PointsRefY[Length - 2]->ValueRec().CurrData() + 1);
////////			}
////////			else
////////			{
////////				Int tTmp = PointsRefY[tVarHdl/3 - 1]->ValueRec().CurrData() + PointsRefY[tVarHdl/3 + 1]->ValueRec().CurrData();
////////				tTmp /= 2;
////////				tVals.annex(tTmp - 1);
////////				tVals.annex(tTmp);
////////				tVals.annex(tTmp + 1);
////////			}
////////
//			minIdx = InvIdx;
//            minVal= MaxInt;
//            for(Idx tIdx = 0; tIdx < 2 * Length; ++tIdx )
////            for(Idx tIdx = 0; tIdx < tVals.size(); ++tIdx )
//            {
//<<<<<<< Updated upstream
////                if(tVals.item(tIdx)>=0 && tVals[tIdx] <= 2*Length) // within the range
////                {
//                    Sys::simulDiffItrVal(handleSystem,tVarHdl,tIdx);
//                    Sys::simulIncr(handleSystem, hardSum.TermHdl);
//                    Int thisVal = HardSumRef.MetricRec().NextData(ss.SimulClk());
//                    if(thisVal < minVal)
//                    {
//                        minVal=thisVal;
//						minIdx = tIdx;
//						sVals.reset();
//						sVals.annex(tIdx);
//                    }
//					else if (thisVal == minVal)
//						sVals.annex(tIdx);
////                }
//            }
//			if (minIdx != InvIdx)
//				Sys::execDiffItrVal(handleSystem,tVarHdl,sVals[theRnd(sVals.size())]);
//
//			++tVarHdl;
////////			tVals.reset();
////////			if (tVarHdl < 3)
////////			{
////////				tVals.annex(PointsRefZ[1]->ValueRec().CurrData() - 1);
////////				tVals.annex(PointsRefZ[1]->ValueRec().CurrData());
////////				tVals.annex(PointsRefZ[1]->ValueRec().CurrData() + 1);
////////			}
////////			else if (tVarHdl >= 3 * Length - 3)
////////			{
////////				tVals.annex(PointsRefZ[Length - 2]->ValueRec().CurrData() - 1);
////////				tVals.annex(PointsRefZ[Length - 2]->ValueRec().CurrData());
////////				tVals.annex(PointsRefZ[Length - 2]->ValueRec().CurrData() + 1);
////////			}
////////			else
////////			{
////////				Int tTmp = PointsRefZ[tVarHdl/3 - 1]->ValueRec().CurrData() + PointsRefZ[tVarHdl/3 + 1]->ValueRec().CurrData();
////////				tTmp /= 2;
////////				tVals.annex(tTmp - 1);
////////				tVals.annex(tTmp);
////////				tVals.annex(tTmp + 1);
////////			}
//////
//			minIdx = InvIdx;
//            minVal= MaxInt;
//            for(Idx tIdx = 0; tIdx < 2 * Length; ++tIdx )
////            for(Idx tIdx = 0; tIdx < tVals.size(); ++tIdx )
//            {
////                if(tVals.item(tIdx)>=0 && tVals[tIdx] <= 2*Length) // within the range
////                {
//                    Sys::simulDiffItrVal(handleSystem,tVarHdl,tIdx);
//=======
//                NewPtrHdl(tmp1,RngVi,handleSystem,varIdx+3);
//                vals.annex(tmp1->ValueRec().CurrData()-1);
//                vals.annex(tmp1->ValueRec().CurrData());
//                vals.annex(tmp1->ValueRec().CurrData()+1);
//            }
//            vals.annex((vals[1]+vals[4])/2);
//
//            //cout << "Iterated Value # " << vals.size() << endl;
//
//
//            // now for each of these values find the best value change
//            Idx minVal=-1;
//            Idx minIdx=-1;
//            //cout<<minVal<<endl;
//			cerr << "value tried ";
//            for(Idx tIdx = 0; tIdx < vals.size(); ++tIdx )
//            {
//                if(vals.item(tIdx)>=0 && vals[tIdx] <= 2*Length) // within the range
//                {
//					cerr << vals.item(tIdx) << " ";
//                    Sys::simulDiffItrVal(handleSystem,varIdx,vals.items()[tIdx]);
//>>>>>>> Stashed changes
//                    Sys::simulIncr(handleSystem, hardSum.TermHdl);
//                    Int thisVal = HardSumRef.MetricRec().NextData(ss.SimulClk());
//                    if(thisVal < minVal)
//                    {
//                        minVal=thisVal;
//						minIdx = tIdx;
//						sVals.reset();
//						sVals.annex(tIdx);
//                    }
//					else if (thisVal == minVal)
//						sVals.annex(tIdx);
////                }
//            }
//			if (minIdx != InvIdx)
//				Sys::execDiffItrVal(handleSystem,tVarHdl,sVals[theRnd(sVals.size())]);
//*/
//
//<<<<<<< Updated upstream
//             Sys::doSelcPairExecDiff(handleSystem, hardSumVarSelcHdl, hardSumValSelcHdl, theRnd);
//        }
//=======
//			cerr << "minIdx " << vals.item(minIdx) << " " << minVal << endl;
//            // min Idx is the selected value
//            if(minIdx != InvIdx)
//            Sys::execDiffItrVal(handleSystem,varIdx,vals.item(6));
//
//
//            // now there should be only onne variable selected
//
//
//             //Sys::doSelcPairExecDiff(handleSystem, hardSumVarSelcHdl, hardSumValSelcHdl, theRnd);
//            // instead of choosing values by by the selector use a customized value selector
//            // to handle the connectivity constraint
//>>>>>>> Stashed changes
//
//        for(Idx tIdx = 1; tIdx <Length; ++tIdx)
//        {
//            Int tx1 = PointsRefX[tIdx - 1]->ValueRec().CurrData();
//            Int ty1 = PointsRefY[tIdx - 1]->ValueRec().CurrData();
//            Int tz1 = PointsRefZ[tIdx - 1]->ValueRec().CurrData();
//            Int tx2 = PointsRefX[tIdx]->ValueRec().CurrData();
//            Int ty2 = PointsRefY[tIdx]->ValueRec().CurrData();
//            Int tz2 = PointsRefZ[tIdx]->ValueRec().CurrData();
//			cout << (tx1 - tx2) * (tx1 - tx2) + (ty1 - ty2) * (ty1 - ty2) + (tz1 - tz2) * (tz1 - tz2) << " ";
//        }
//
////        for(Idx tIdx = 0; tIdx <Length; ++tIdx)
////        {
////            cout << PointsRefX[tIdx]->ValueRec().CurrData() << " "
////            << PointsRefY[tIdx]->ValueRec().CurrData() << " "
////            << PointsRefZ[tIdx]->ValueRec().CurrData() << " | ";
////        }
//        cout << " Clk: " << ss.ExecClk() << endl;
//        if(HardSumRef.MetricRec().CurrData() == 0)
//        {
//            ++countFeasible;
//            cout << "obj: " << DstSumRef.ValueRec().CurrData() << endl;
//        }
//        if(HardSumRef.MetricRec().CurrData() == 0 && DstSumRef.ValueRec().CurrData() > globalbest)
//        {
//            globalbest =  DstSumRef.ValueRec().CurrData();
//
//            Sys::setVarLocked(handleSystem, ObjLimit.TermHdl, false);
//            Sys::execDiffItrVal(handleSystem, ObjLimit.TermHdl,globalbest);
//            Sys::setVarLocked(handleSystem, ObjLimit.TermHdl, true);
//
//            cout << "obj locked to: " << globalbest << endl;
//            writeHPToFile("data.cml", Protein, PointsRefX, PointsRefY, PointsRefZ);
//        }
//    }
//
//
//    cerr << ss.ExecClk() << " Global Best: "<< globalbest << " No Of feasible Solutions: " << countFeasible << endl;
//
//
//    return 1;
//
//}
//


//////int proteinSPWithAbsDirModel()
//////{
//////
//////        #define ORIGIN 500
//////        char theProtein [ 512 ] =  "HPHHPPHHHHPHHHPPHHPPHPHHHPHPHHPPHHPPPHPPPPPPPPHH" ;
//////        N Length = strlen ( theProtein );
//////
//////        Cnt T = 0;//parseN(ArgV[1]);
//////
//////		Rnd theRnd;//1336536274
//////		cerr << "The Random Seed:" << theRnd.Seed() << endl;
//////
//////		Hdl tSys = Sys::def();
//////		EvalTi::def(tSys);
//////		HintTi::def(tSys);
//////
//////		kb<Prm> Vars(Length - 1);
//////		kb<Prm> Xvecs(Length - 1), Yvecs(Length - 1), Zvecs(Length - 1);
//////		Int Xs[12] = {0, 0, 0, 0, 	1, 1, -1, -1, 	1, 1, -1, -1};
//////		Int Ys[12] = {1, 1, -1, -1, 0, 0, 0, 0, 	1, -1, 1, -1};
//////		Int Zs[12] = {1, -1, 1, -1, 1, -1, 1, -1, 	0, 0, 0, 0};
//////		kb<Prm> Xpoints(Length), Ypoints(Length), Zpoints(Length), XYZ(Length);
//////		Xpoints[0] = Ypoints[0] = Zpoints[0] = Ci::def(tSys, ORIGIN);
//////		XYZ[0] = TpackXiFeVi::def(Xv, tSys, Xpoints[0], Ypoints[0], Zpoints[0]);
//////
//////		for(Idx tIdx = 0; tIdx < Length - 1; ++tIdx)
//////		{
//////			Vars[tIdx] = RngVi::def(tSys, 0, 11);
//////			Xvecs[tIdx] = UilXiKiFeVi::def(Xv, tSys, Vars[tIdx], UilXiKiFeVi::bind(Xs, 12));
//////			Yvecs[tIdx] = UilXiKiFeVi::def(Xv, tSys, Vars[tIdx], UilXiKiFeVi::bind(Ys, 12));
//////			Zvecs[tIdx] = UilXiKiFeVi::def(Xv, tSys, Vars[tIdx], UilXiKiFeVi::bind(Zs, 12));
//////
//////			Xpoints[tIdx + 1] = BaddXiFeVi::def(Xv, tSys, Xpoints[tIdx], Xvecs[tIdx]);
//////			Ypoints[tIdx + 1] = BaddXiFeVi::def(Xv, tSys, Ypoints[tIdx], Yvecs[tIdx]);
//////			Zpoints[tIdx + 1] = BaddXiFeVi::def(Xv, tSys, Zpoints[tIdx], Zvecs[tIdx]);
//////			XYZ[tIdx + 1] = TpackXiFeVi::def(Xv, tSys, Xpoints[tIdx + 1], Ypoints[tIdx + 1], Zpoints[tIdx + 1]);
//////		}
//////		Prm AllDiff = Prm(AllDiffXiFcMiD::def(Xv, tSys, XYZ.items(), Length), Dd);
//////
//////
//////
//////
//////
//////		xb<Prm> Dists;
//////		for(Idx tIdx1 = 0; tIdx1 < Length; ++tIdx1)
//////		{
//////			if (theProtein[tIdx1] == 'P') continue;
//////			for(Idx tIdx2 = tIdx1 + 2; tIdx2 < Length; ++tIdx2)
//////			{
//////				if (theProtein[tIdx2] == 'P') continue;
//////				Prm tXYZ[3];
//////				tXYZ[0] = BdiffXiFeVi::def(Xv, tSys, Xpoints[tIdx1], Xpoints[tIdx2]);
//////				tXYZ[1] = BdiffXiFeVi::def(Xv, tSys, Ypoints[tIdx1], Ypoints[tIdx2]);
//////				tXYZ[2] = BdiffXiFeVi::def(Xv, tSys, Zpoints[tIdx1], Zpoints[tIdx2]);
//////				Dists.annex(SumEqsXiBiFeVi::def(Xv, tSys, tXYZ, 3, SumEqsXiBiFeVi::bind(2)));
//////			}
//////		}
//////		Prm ObjVar = RngVi::def(tSys,0,1000);
//////		Prm DstSum = SumXiFeVi::def(Xv, tSys, Dists.items(), Dists.size());
//////		Prm ObjLim = Prm(BltuXiFcMi::def(Xv, tSys, ObjVar, DstSum), Lm);
//////
//////		Prm PrmArr[2] = {AllDiff, ObjLim};
//////		Prm TopSum = SumXiOFcMiD::def(Xm|On, tSys, PrmArr, 2);
//////
//////		Prm Heap = TabuMaxHeapHiFrHb::def(tSys, TopSum);
//////
//////		Hdl const VarSelcHdl = RankVarSp::def( tSys, Heap);
//////		Hdl const ValSelcHdl = MinValSdXi::def( Xm, tSys, TopSum);
//////
//////
//////		Prm hardSumHeap = TabuMaxHeapHiFrHi::def(tSys, AllDiff);
//////        Hdl const hardSumVarSelcHdl = RankVarSp::def( tSys, hardSumHeap);
//////        Hdl const hardSumValSelcHdl = MinValSdXi::def( Xm, tSys, AllDiff);
//////
//////
//////
//////        T = Length/4; //---length of tabu
//////
//////		Sys::setVarTabuLimit(tSys, T);
//////		Sys::execAnewRndVars(tSys, theRnd);
//////
//////
//////
//////
//////		//#if CompLazy
//////		//Sys::setSelcActive(tSys, VarSelcHdl, true);
//////		//Sys::setSelcActive(tSys, ValSelcHdl, true);
//////		//#endif
//////
//////		Sys const & ss = Sys::ref(tSys);
////////		RefHdl(TopSumRef, SumXiOFcMiD, tSys, TopSum.TermHdl);
//////		RefHdl(DstSumRef, SumXiFeVi, tSys, DstSum.TermHdl);
//////		RefHdl(AllDiffRef, AllDiffXiFcMiD, tSys, AllDiff.TermHdl);
//////
//////        Itn I = 1000000;//46;//81;//parseN(ArgV[2]);
//////
//////        Int hardThreshold = 15;
//////
//////        Int globalbest=0;
//////
//////        Sys::execDiffItrVal(tSys, ObjVar.TermHdl,0);
//////        Sys::setVarLocked(tSys, ObjVar.TermHdl, true);
//////        Sys::setVarTabu(tSys, ObjVar.TermHdl,true);
//////
//////
//////        while(ss.ExecClk() <= I) // only iteration
//////        {
//////            //cerr << AllDiffRef.MetricRec().CurrData() << " " << TopSumRef.MetricRec().CurrData() << endl;
//////
//////            if( AllDiffRef.MetricRec().CurrData() <= hardThreshold ) // top
//////            {
//////                Sys::doSelcPairExecDiff(tSys, VarSelcHdl, ValSelcHdl, theRnd);
//////            }
//////            else // hard
//////            {
//////                Sys::doSelcPairExecDiff(tSys, hardSumVarSelcHdl, hardSumValSelcHdl, theRnd);
//////            }
//////
//////            if(AllDiffRef.MetricRec().CurrData() == 0 && DstSumRef.ValueRec().CurrData() > globalbest)
//////            {
//////                globalbest =  DstSumRef.ValueRec().CurrData();
//////
//////                Sys::setVarLocked(tSys, ObjVar.TermHdl, false);
//////                Sys::execDiffItrVal(tSys, ObjVar.TermHdl,globalbest);
//////                Sys::setVarLocked(tSys, ObjVar.TermHdl, true);
//////
//////                cout << "obj locked to: " << globalbest << " " << ss.ExecClk() << endl;
//////
//////            }
//////
//////
//////
//////        }
//////
//////        return 1;
//////}

void inputFromFile(string filename, Itr **latticePoints, int length)
{
    std::ifstream file(filename.c_str());
	if(!file.is_open()){cout<<"File not found!"; exit(0);}
	for(int i=0;i<length;i++)
	{
		file>>latticePoints[i][0];
		file>>latticePoints[i][1];
		file>>latticePoints[i][2];
		latticePoints[i][0]+=length;
		latticePoints[i][1]+=length;
		latticePoints[i][2]+=length;

	}
	for(int i=0;i<length;i++)
		{
			cout<<latticePoints[i][0]<<" "<<latticePoints[i][1]<<" "<<latticePoints[i][2]<<endl;
		}
	file.close();

}

void writeHPToFile(const string filename, char* Protein, const RngVi** PointsX, const RngVi** PointsY, const RngVi** PointsZ)
{

	std::ofstream out(filename.c_str());
	 /*
	for(int i=0;i<protein->sequenceLength;i++)
	{
		out<<latticePoints[i][0]<<" "<<latticePoints[i][1]<<" "<<latticePoints[i][2]<<" "<<(char)('A'+protein->getProteinType(i))<<endl;;
	}
	out.close();*/
	int distance=3;
		int shift = strlen(Protein)*distance;	// shift output to display with positive dimensions only

				out	<<"<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>"
					<<"\n\n<!-- FOR BEST VIEWING e.g. IN JMOL USE AFTER LOADING THE SCRIPT -->"
					<<"\n<!-- select all; color bonds grey; -->"
					<<"\n\n<!-- created with CPSP-package : http://www.bioinf.uni-freiburg.de/sw/cpsp/ -->"
					<<"\n"
					<<"\n<list dictRef=\"cdk:model\" xmlns=\"http://www.xml-cml.org/schema\">"
					<<"\n  <list>"
					<<"\n    <molecule id=\"lattice-protein\">"
					<<"\n    <!-- model    : HP-model -->"
					<<"\n    <!-- sequence : " <<Protein<<" -->"
					<<"\n    <!-- encoding : H = Cl-atom, P = C-atom  -->"
					<<"\n    <!-- lattice  : " <<"fcc"<<" -->"
					<<"\n      <atomArray>";

//					out <<"\n        <atom id=\""<<1<<"_"<<Protein[0]<<"\""
//						<<" elementType=\"C"<<(Protein[0]=='H'?"l":"")<<"\""
//						<<" x3=\""<<(shift + (distance * strlen(Protein)))<<"\""
//						<<" y3=\""<<(shift + (distance * strlen(Protein)))<<"\""
//						<<" z3=\""<<(shift + (distance * strlen(Protein)))<<"\""
//						<<" />";
				for (Idx i = 0; i < strlen(Protein); ++i) {
					out <<"\n        <atom id=\""<<i+1<<"_"<<Protein[i]<<"\""
						<<" elementType=\"C"<<(Protein[i]=='H'?"l":"")<<"\""
						<<" x3=\""<<(shift + (distance * (*PointsX[i]).CurrItr()))<<"\""
						<<" y3=\""<<(shift + (distance * (*PointsY[i]).CurrItr()))<<"\""
						<<" z3=\""<<(shift + (distance * (*PointsZ[i]).CurrItr()))<<"\""
						<<" />";
				}
				out	<<"\n      </atomArray>"
					<<"\n      <bondArray>";
				for (Idx i = 1; i < strlen(Protein); i++) {
					out	<<"\n        <bond id=\"b"<<i<<"\" atomRefs2=\""<<i<<"_"<<Protein[i-1]<<" "<<i+1<<"_"<<Protein[i]<<"\" order=\"S\"/>";
				}
				out	<<"\n      </bondArray>"
					<<"\n    </molecule>"
					<<"\n  </list>"
					<<"\n</list>"
					<<"\n\n<!-- FOR BEST VIEWING e.g. IN JMOL USE THE SCRIPT -->"
					<<"\n<!-- select all; color bonds grey; -->"
					<<std::endl;
				out.flush();
				out.close();
}

void new_init(char* protein, kblock2<Int> & latticePoints, Rnd & theRnd)
{
	int ran;
	int ranf=0;
	int pTurn=0;
	int pUp=0;
	int up;
	int turn;
	int mv;
	int mv1;
	int mvUp;
	int mvP;
	int mvP_;
	int mv1_;
	int mvS;
	int *forbiden=(int*)calloc(13,sizeof(int));


//	for(int i=0;i<13;i++)
//	{
//		forbiden[i]=0;
//	}

	int pUpRate = 20; // initialization parameters
	int pTuRate = 15;

	// WE FIRSTLY INITIALIZE THE COORDINATES
	for (int i=1;i<=12;i++)
		if(((i-3)%4==0 && i>=3))
		{
			forbiden[i]=1;
			forbiden[i+1]=1;
		}
    int sequenceLength = strlen(protein);
	latticePoints[0][0]=sequenceLength;
	latticePoints[0][1]=sequenceLength;
	latticePoints[0][2]=sequenceLength;
	if ((latticePoints[0][0]+latticePoints[0][1]+latticePoints[0][2])%2 != 0)
        latticePoints[0][2]=(latticePoints[0][2]+1);//%bound; // it must be even!


	do{
		//ran=(rand()%12)+1;
		ran=uniform(theRnd,castZ(1),castZ(12));
	}while(forbiden[ran]==1);

	mv=ran;
	mv1=ran;

	if(mv1%2!=0)
		mv1_=mv1+1;
	else
		mv1_=mv1-1;

	updateForbiden(forbiden,mv1);

	do{
		//		ran=(rand()%12)+1;
		ran=uniform(theRnd,castZ(1),castZ(12));
	}while(forbiden[ran]==1);
	mvP=ran;

	if(mvP%2!=0) mvP_=mvP+1;
	else mvP_=mvP-1;

	forbiden[mvP]=1;

	forbiden[mvP_]=1;

	do{
//		ran=(rand()%12)+1;
		ran=uniform(theRnd,castZ(1),castZ(12));
	}while(forbiden[ran]==1);
	mvUp=ran;



	for(int i=1;i<sequenceLength;i++)
	{
//		ranf=rand()%101;
		ranf=uniform(theRnd,castZ(0),castZ(100));
		turn=0;
		up=0;
		//printf("ran %f, pT %f pU %f\n",ranf,pTurn,pUp);
		if(pTurn>ranf){
//			ranf=rand()%101;
			ranf=uniform(theRnd,castZ(0),castZ(100));
			if(pUp>ranf){
				mv=mvUp;
				pUp=0;
				up=1;
				pTurn=0;
			}else{
				pUp=pUp+pUpRate;
				mv=mvP;
				pTurn=0;
				turn=1;
			}
		}else{
			pTurn+=pTuRate;
		}

		for(int _i=0;_i<3;_i++)
		{
			latticePoints[i][_i]=latticePoints[i-1][_i]+fccNeighborData[transformMoves[mv]].vec[_i];
		}

  		if(up==1)
  		{
  			mv=mv1_;
  			//swapp moves
  			mvS=mv1;
  			mv1=mv1_;
  			mv1_=mvS;
  			mvS=mvP;
  			mvP=mvP_;
  			mvP_=mvS;
  		}
  		else
  		{
  			if(turn==1)
  			{
  				mv=mv1_;
  				//swapp moves
  				mvS=mv1;
  				mv1=mv1_;
  				mv1_=mvS;
  			}

  		}

  }



}

void updateForbiden(int *f, int m)
{

	if((m-1)%4<2)
	{
		for(int i=1;i<=12;i++)
			if(((i-3)%4==0 && i>=3))
			{

				f[i]=1;
				f[i+1]=1;
			}
	}
	else
	{
		for(int i=1;i<=12;i++)
			if(((i)%4==1))
			{
				f[i]=1;
				f[i+1]=1;
			}
	}

	if(m<5)
	{
		for(int i=1;i<=4;i++)
			f[i]=1;
	}
	else if(m<9)
	{
		for(int i=4;i<=8;i++)
    		  f[i]=1;
	}
	else if(m<13)
	{
		for(int i=8;i<=12;i++)
			f[i]=1;
	}

}



closeKangarooSpace



