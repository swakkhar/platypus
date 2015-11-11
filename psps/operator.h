#include "pspl/idx.hh"
#include<queue>

using namespace platypus;
using namespace kangaroo;
using namespace std;

void operatorSearch(int argC, char* argv[])
{

    if(argC<=3)
    {
        cout << "Not enough argument provided!" << endl;
        cout << "Usage:" << argv[0] << " <sequence> <timeout>" << endl;
        return;
    }

    block1<RPoint,xmm> allPoints;
        int start=0;
        int end=0;

    Int timeOut = parseN(argv[2]);
    Rnd theRnd;
    cout << "seed: " <<theRnd.Seed()<<endl;
    Space::Dimen(3);
    FccLattice fcc;
    Lattice::l(fcc);

    BarerraEnergy br;
    Energy::e(br);

    Protein p(20,argv[1]);
    Protein::p(p);

    PDBReader::readPDBFile(argv[3],allPoints,start,end);



	DoublePullMove pull;
    CrankShaftMove crank;
    SinglePullMove single;
    DiagonalMove dia;
    TriplePullMove triple;
    PushMove push;
    TiltMove tilt;
    TwistMove twist;
    LNSMove lns;

	Conf c(theRnd);

	EnergyObj bre(c,true);
    ACoreDistObj acdo(c,true);
    FreeFtCns free(c,true);

    RMSDObj rmsd(c,true,allPoints);
    //RepCntCns rep(c,true);

    //RandomValid rs(theRnd);
    ChainGrowth rs(theRnd);

    rs.compute(c);

    // hint system


	Sys const & tSys = Sys::refc(c.SysHdl());
    //
    block1<Prm,nmm> prmArr(3);
    prmArr.insert(Prm(Tf,free.FuncHdl(),HintMinAsHintMin));
    prmArr.insert(Prm(Tf,acdo.FuncHdl(),HintMinAsHintMin));
    prmArr.insert(Prm(Tf,bre.FuncHdl(),HintMinAsHintMin));
    //prmArr.insert(Prm(Tf,rep.FuncHdl(),HintMinAsHintMin));
    Prm hintSum=Prm(TermFunc,SumXiHFcMiHn::def(Xm|HintMin, tSys.SysHdl, prmArr.items(), prmArr.itemCount()));

    Prm HintHeapEn = Prm(Tf, Sv2TabuMaxHeapHiFrHi::def(tSys.SysHdl, hintSum));
    //Prm HintHeapEn = Prm(Tf, Sv2TabuMaxHeapHiFrHi::def(tSys.SysHdl, hintSum));
    //Prm HintHeapEn = Prm(Tf, Sv2TabuMaxHeapHiFrHi::def(tSys.SysHdl, Prm(Tf,free.FuncHdl())));
	Hdl const VarSelcHdl = RankedHintVar1Sp::def( tSys.SysHdl, HintHeapEn);

	c.initialise(rs.FullChange());

    R globalRMSDMin=0;

    Int globalMin= bre.ExecMetric();
    cout << "global best: "<<globalMin << endl;
    c.writeHPToFile("data.cml");
    exit(0);

    Int initWalkPV=15;
    Int walkPV=initWalkPV; // 5% walk probability


    TimeUtl timeUtl;
    R timeInit= timeUtl.getTime();

    block1<Change,xmm> partialChanges;

    typedef tuple2a<Int,Idx> mytuple;
    priority_queue<mytuple,vector<mytuple>,greater_equal<mytuple> > pQueue; // using default comparator

    Idx proteinLength = Protein::p().Length();
    //Idx varibaleSelectorIdx = 0;

    block1<Idx,kmm> tabuList(Protein::p().Length());
    Idx tabuLength = uniform(theRnd,4,(Z)Protein::p().Span()/8);
    // initialize tabu
    for(Idx tIdx=0;tIdx < Protein::p().Length(); ++tIdx)
        tabuList[tIdx]=0;

    Idx nonImp =0;
    Idx initMaxStable = 1000;
    Idx maxStable=initMaxStable;
    R mFactor=1.2;

    Change bestPoints;
    block1<Change,xmm> elitSet;

    while(timeUtl.getTime()-timeInit <= timeOut)
    {
        Idx stagLoop=0;
        partialChanges.clear();

        // at first random toss
        Int wp=uniform(theRnd,0,100);
        platypus::Pos tSelcPos = 0;
        Hdl selectedVar=0;
        if(wp<walkPV) //random walk
        {
            // select randomly a position that is not in the meta level tabu
            while(1)
            {
                tSelcPos=uniform(theRnd,(Z)0,(Z)proteinLength-1);
                if(tabuList[tSelcPos]< tSys.ExecClk())
                    break;
                stagLoop++;
                if(stagLoop>=proteinLength)
                {
                    break;
                }
            }
        }
        else
        {
            while(1)
            {
                    // select a variable from the variable selector
                Selc::ptrm(tSys.SysHdl,VarSelcHdl)->performSelection(theRnd);
                selectedVar= Selc::ptrm(c.SysHdl(), VarSelcHdl)->SelcVars()[0];
                tSelcPos= selectedVar/Space::Dimen();
                // if tPos is is meta-level tabu then consider next variable
                if(tabuList[tSelcPos]< tSys.ExecClk())
                    break;
                else
                    QcSv2Tabu::refm(c.SysHdl()).tabuVar(selectedVar);
                stagLoop++;
                if(stagLoop>=proteinLength)
                {
                    break;
                }
            }

        }
        if(stagLoop>=proteinLength)
        {
            // stagnation
            walkPV=walkPV*mFactor;
            nonImp=0;
            //clear tabu
            for(Idx tIdx =0; tIdx < proteinLength;++tIdx)
               {
                    tabuList[tIdx]=0;
               }

            // initialize with best saved points
            c.execute(bestPoints);

        }
        // now there is a danger in the inner loop it may live forever and find no moves
        // or it may find move but that may lead to stagnation, we need something to handle stagnation
        // restart - what type of restart? best position ? possible



        // now enumerate changes for all the operators on position tSelcPos


        //cout << "tPos: " << tSelcPos <<endl;



        /*dia.compute(c,tSelcPos);
        // copy all the changes into partialChanges
        for(Idx tIdx = 0; tIdx < dia.PartChanges().itemCount(); ++tIdx)
        {
            partialChanges.insertMem(dia.PartChanges()[tIdx]);
        }



        crank.compute(c,tSelcPos);
        // copy all the changes into partialChanges
        for(Idx tIdx = 0; tIdx < crank.PartChanges().itemCount(); ++tIdx)
        {
            partialChanges.insertMem(crank.PartChanges()[tIdx]);
        }

        triple.compute(c,tSelcPos);
        // copy all the changes into partialChanges
        for(Idx tIdx = 0; tIdx < triple.PartChanges().itemCount(); ++tIdx)
        {
            partialChanges.insertMem(triple.PartChanges()[tIdx]);
        }*/

        if(tSelcPos>0)
        {


        lns.reset();
        lns.addPosition(tSelcPos);
        lns.compute(c,tSelcPos);
        for(Idx tIdx = 0; tIdx < lns.PartChanges().itemCount(); ++tIdx)
        {
            partialChanges.insertMem(lns.PartChanges()[tIdx]);
        }

        lns.reset();
        lns.addPosition(tSelcPos);
        if(tSelcPos<Protein::p().Span())
        {
            lns.addPosition(tSelcPos+1);
            lns.compute(c,tSelcPos);
            for(Idx tIdx = 0; tIdx < lns.PartChanges().itemCount(); ++tIdx)
            {
                partialChanges.insertMem(lns.PartChanges()[tIdx]);
            }
        }
        lns.reset();
        lns.addPosition(tSelcPos);
        if(tSelcPos<Protein::p().Span()-1)
        {
            lns.addPosition(tSelcPos+1);
            lns.addPosition(tSelcPos+2);
            lns.compute(c,tSelcPos);
            for(Idx tIdx = 0; tIdx < lns.PartChanges().itemCount(); ++tIdx)
            {
                partialChanges.insertMem(lns.PartChanges()[tIdx]);
            }
        }
        }

        push.compute(c,tSelcPos);
        // copy all the changes into partialChanges
        for(Idx tIdx = 0; tIdx < push.PartChanges().itemCount(); ++tIdx)
        {
            partialChanges.insertMem(push.PartChanges()[tIdx]);
        }

        single.compute(c,tSelcPos);
        // copy all the changes into partialChanges
        for(Idx tIdx = 0; tIdx < single.PartChanges().itemCount(); ++tIdx)
        {
            partialChanges.insertMem(single.PartChanges()[tIdx]);
        }


        // also include all other moves
        pull.compute(c,tSelcPos);
        // copy all the changes into partialChanges
        for(Idx tIdx = 0; tIdx < pull.PartChanges().itemCount(); ++tIdx)
        {
            partialChanges.insertMem(pull.PartChanges()[tIdx]);
        }

        //cout << "total partial changes:" << partialChanges.itemCount() << endl;
        // now simulate the partial changes
        for(Idx tIdx = 0; tIdx < partialChanges.itemCount(); ++tIdx)
        {
            //cout << bre.SimulMetric() <<endl;
            c.simulate(partialChanges[tIdx]);

            int heu=uniform(theRnd,0,2);
            if(heu<=1)
            {
                Term::performSimulIncr(Func::ptrm(c.SysHdl(), bre.FuncHdl()));
                mytuple t(bre.SimulMetric(),tIdx);
                pQueue.push(t);
            }
            else
            {
                Term::performSimulIncr(Func::ptrm(c.SysHdl(), acdo.FuncHdl()));
                mytuple t(acdo.SimulMetric(),tIdx);
                pQueue.push(t);
            }

        }
        // now all the changes are in the que
        // randomly tie break from the top
        Idx selectedId=1;
        //Int selectedEn=0;
        if(!pQueue.empty())
        {
            mytuple top=pQueue.top();
            block1<mytuple,xmm> vector;
            while(!pQueue.empty()&&top==pQueue.top())
            {
                vector.insertMem(pQueue.top());
                pQueue.pop();
            }
            //cout << "vector count " << vector.itemCount() << endl;
            // now randomly select from the elements - random tie break
            Idx selRank=uniform(theRnd,(N)0,vector.itemCount()-1);
            selectedId=vector[selRank].Second;
        }
        else
        {
            //stagnation
            //cout << "Stagnation! " << endl<<endl;
            // probably select the next variable or randomly
            // tabu this variable
            if(wp>=walkPV)
                QcSv2Tabu::refm(c.SysHdl()).tabuVar(selectedVar);
            tabuList[tSelcPos]=tSys.ExecClk()+tabuLength;
            nonImp++;
            walkPV=walkPV*mFactor;
            continue;
        }
        //cout << "selected id:"<<selectedId  << "partial changes" << partialChanges.itemCount()<< endl;
        //cout << partialChanges[selectedId].size()<< endl;
        c.execute(partialChanges[selectedId]);
        tabuList[tSelcPos]=tSys.ExecClk()+tabuLength;
        if(bre.ExecMetric() < globalMin)
        {
            globalMin = bre.ExecMetric();
            globalRMSDMin=0.1*sqrt(2*rmsd.ExecMetric()/(proteinLength*(proteinLength-1)));
            cout << "global best: " << globalMin << " "<<tSys.ExecClk()<<endl;
            c.savePoints(bestPoints);
            elitSet.insertMem(bestPoints);
            nonImp=0;
            walkPV=initWalkPV;
            maxStable=initMaxStable;
        }
        else
        {
            nonImp++;
        }
        if(nonImp>maxStable)
        {
            // stagnation
            cout << "restart: " << bre.ExecMetric()<<" "<<tSys.ExecClk()<<endl;
            walkPV=walkPV*mFactor;
            maxStable=maxStable*mFactor;
            nonImp=0;
            //clear tabu
            for(Idx tIdx =0; tIdx < proteinLength;++tIdx)
               {
                    tabuList[tIdx]=0;
               }

            // initialize with best saved points from elit set // so each time a new
            // local minima is used for restart // not the global minima always
            bestPoints=elitSet[uniform(theRnd,(Dim)0,elitSet.itemCount()-1)];
            c.execute(bestPoints);
        }
        while(!pQueue.empty())pQueue.pop();

    }
    cout << "global min:" << globalMin << " "<<globalRMSDMin<<endl;
    for(Idx tIdx=0;tIdx<Protein::p().Length();++tIdx)
    {
           cout << c.Coord(tIdx).Comp(0) << " " << c.Coord(tIdx).Comp(1) << " " << c.Coord(tIdx).Comp(2) << endl;
    }

    return;
}

/*
void operatorSearch(int argC, char* argv[])
{

    if(argC<2)
    {
        cout << "Not enough argument provided!" << endl;
        cout << "Usage:" << argv[0] << " <sequence> <timeout>" << endl;
        return;
    }


    Int timeOut = parseN(argv[2]);
    Rnd theRnd;
    cout << "seed: " <<theRnd.Seed()<<endl;
    Space::Dimen(3);
    FccLattice fcc;
    Lattice::l(fcc);

    BarerraEnergy br;
    Energy::e(br);

    Protein p(20,argv[1]);
    Protein::p(p);


	DoublePullMove pull;
    CrankShaftMove crank;
    SinglePullMove single;
    DiagonalMove dia;
    TriplePullMove triple;
    PushMove push;
    TiltMove tilt;
    TwistMove twist;

	Conf c(theRnd);

	EnergyObj bre(c,true);
    ACoreDistObj acdo(c,true);
    //RandomSphere rs(theRnd);
    ChainGrowth rs(theRnd);

    rs.compute(c);

	Sys const & tSys = Sys::refc(c.SysHdl());
    Prm HintHeapEn = Prm(Tf, Sv2TabuMaxHeapHiFrHi::def(tSys.SysHdl, Prm(Tf,acdo.FuncHdl())));
	Hdl const VarSelcHdl = RankedHintVar1Sp::def( tSys.SysHdl, HintHeapEn);

	c.initialise(rs.FullChange());



    Int globalMin= bre.ExecMetric();
    cout << "global min:"<<globalMin << endl;

    Int initWalkPV=15;
    Int walkPV=initWalkPV; // 5% walk probability


    TimeUtl timeUtl;
    R timeInit= timeUtl.getTime();

    block1<Change,xmm> partialChanges;

    typedef tuple2a<Int,Idx> mytuple;
    priority_queue<mytuple,vector<mytuple>,greater_equal<mytuple> > pQueue; // using default comparator

    Idx protenLength = Protein::p().Length();
    //Idx varibaleSelectorIdx = 0;

    block1<Idx,kmm> tabuList(Protein::p().Length());
    Idx tabuLength = uniform(theRnd,4,(Z)Protein::p().Span()/8);
    // initialize tabu
    for(Idx tIdx=0;tIdx < Protein::p().Length(); ++tIdx)
        tabuList[tIdx]=0;

    Idx nonImp =0;
    Idx initMaxStable = 5000;
    Idx maxStable=initMaxStable;
    R mFactor=1.2;

    Change bestPoints;

    while(timeUtl.getTime()-timeInit <= timeOut)
    {
        Idx stagLoop=0;
        partialChanges.clear();

        // at first random toss
        Int wp=uniform(theRnd,0,100);
        platypus::Pos tSelcPos = 0;
        Hdl selectedVar=0;
        if(wp<walkPV) //random walk
        {
            // select randomly a position that is not in the meta level tabu
            while(1)
            {
                tSelcPos=uniform(theRnd,(Z)0,(Z)protenLength-1);
                if(tabuList[tSelcPos]< tSys.ExecClk())
                    break;
                stagLoop++;
                if(stagLoop>=protenLength)
                {
                    break;
                }
            }
        }
        else
        {
            while(1)
            {
                    // select a variable from the variable selector
                Selc::ptrm(tSys.SysHdl,VarSelcHdl)->performSelection(theRnd);
                selectedVar= Selc::ptrm(c.SysHdl(), VarSelcHdl)->SelcVars()[0];
                tSelcPos= selectedVar/Space::Dimen();
                // if tPos is is meta-level tabu then consider next variable
                if(tabuList[tSelcPos]< tSys.ExecClk())
                    break;
                else
                    QcSv2Tabu::refm(c.SysHdl()).tabuVar(selectedVar);
                stagLoop++;
                if(stagLoop>=protenLength)
                {
                    break;
                }
            }

        }
        if(stagLoop>=protenLength)
        {
            // stagnation
            walkPV=walkPV*mFactor;
            nonImp=0;
            //clear tabu
            for(Idx tIdx =0; tIdx < protenLength;++tIdx)
               {
                    tabuList[tIdx]=0;
               }

            // initialize with best saved points
            c.execute(bestPoints);

        }
        // now there is a danger in the inner loop it may live forever and find no moves
        // or it may find move but that may lead to stagnation, we need something to handle stagnation
        // restart - what type of restart? best position ? possible



        // now enumerate changes for all the operators on position tSelcPos


        //cout << "tPos: " << tSelcPos <<endl;



        dia.compute(c,tSelcPos);
        // copy all the changes into partialChanges
        for(Idx tIdx = 0; tIdx < dia.PartChanges().itemCount(); ++tIdx)
        {
            partialChanges.insertMem(dia.PartChanges()[tIdx]);
        }



        crank.compute(c,tSelcPos);
        // copy all the changes into partialChanges
        for(Idx tIdx = 0; tIdx < crank.PartChanges().itemCount(); ++tIdx)
        {
            partialChanges.insertMem(crank.PartChanges()[tIdx]);
        }

        triple.compute(c,tSelcPos);
        // copy all the changes into partialChanges
        for(Idx tIdx = 0; tIdx < triple.PartChanges().itemCount(); ++tIdx)
        {
            partialChanges.insertMem(triple.PartChanges()[tIdx]);
        }


        push.compute(c,tSelcPos);
        // copy all the changes into partialChanges
        for(Idx tIdx = 0; tIdx < push.PartChanges().itemCount(); ++tIdx)
        {
            partialChanges.insertMem(push.PartChanges()[tIdx]);
        }

        single.compute(c,tSelcPos);
        // copy all the changes into partialChanges
        for(Idx tIdx = 0; tIdx < single.PartChanges().itemCount(); ++tIdx)
        {
            partialChanges.insertMem(single.PartChanges()[tIdx]);
        }


        // also include all other moves
        pull.compute(c,tSelcPos);
        // copy all the changes into partialChanges
        for(Idx tIdx = 0; tIdx < pull.PartChanges().itemCount(); ++tIdx)
        {
            partialChanges.insertMem(pull.PartChanges()[tIdx]);
        }

        //cout << "total partial changes:" << partialChanges.itemCount() << endl;
        // now simulate the partial changes
        for(Idx tIdx = 0; tIdx < partialChanges.itemCount(); ++tIdx)
        {
            //cout << bre.SimulMetric() <<endl;
            c.simulate(partialChanges[tIdx]);
            Term::performSimulIncr(Func::ptrm(c.SysHdl(), acdo.FuncHdl()));
            mytuple t(acdo.SimulMetric(),tIdx);
            pQueue.push(t);
        }
        // now all the changes are in the que
        // randomly tie break from the top
        Idx selectedId=1;
        //Int selectedEn=0;
        if(!pQueue.empty())
        {
            mytuple top=pQueue.top();
            block1<mytuple,xmm> vector;
            while(!pQueue.empty()&&top==pQueue.top())
            {
                vector.insertMem(pQueue.top());
                pQueue.pop();
            }
            //cout << "vector count " << vector.itemCount() << endl;
            // now randomly select from the elements - random tie break
            Idx selRank=uniform(theRnd,(N)0,vector.itemCount()-1);
            selectedId=vector[selRank].Second;
        }
        else
        {
            //stagnation
            //cout << "Stagnation! " << endl<<endl;
            // probably select the next variable or randomly
            // tabu this variable
            if(wp>=walkPV)
                QcSv2Tabu::refm(c.SysHdl()).tabuVar(selectedVar);
            tabuList[tSelcPos]=tSys.ExecClk()+tabuLength;
            nonImp++;
            walkPV=walkPV*mFactor;
            continue;
        }
        //cout << "selected id:"<<selectedId  << "partial changes" << partialChanges.itemCount()<< endl;
        //cout << partialChanges[selectedId].size()<< endl;
        c.execute(partialChanges[selectedId]);
        tabuList[tSelcPos]=tSys.ExecClk()+tabuLength;
        if(bre.ExecMetric() < globalMin)
        {
            globalMin = bre.ExecMetric();
            cout << "global min:" << globalMin << " "<<tSys.ExecClk()<<endl;
            c.savePoints(bestPoints);
            nonImp=0;
            walkPV=initWalkPV;
        }
        else
        {
            nonImp++;
        }
        if(nonImp>maxStable)
        {
            // stagnation
            cout << "Stagnation! " << " "<<tSys.ExecClk()<<endl;
            walkPV=walkPV*mFactor;
            maxStable=maxStable*mFactor;
            nonImp=0;
            //clear tabu
            for(Idx tIdx =0; tIdx < protenLength;++tIdx)
               {
                    tabuList[tIdx]=0;
               }

            // initialize with best saved points
            c.execute(bestPoints);
        }
        while(!pQueue.empty())pQueue.pop();

    }
    cout << "global min:" << globalMin << endl;


    return;
}
*/
