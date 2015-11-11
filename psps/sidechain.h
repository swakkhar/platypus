#include "pspl/idx.hh"

using namespace kangaroo;
using namespace platypus;
using namespace std;

typedef block1<Dir,kmm> individual;

void encode(individual & dirs,Conf const &c)
{
    for(Idx tIdx=0;tIdx<Protein::p().Length()/2-1;++tIdx)
    {
        dirs[tIdx]=Lattice::l().Direction(c.Coord(tIdx+1)-c.Coord(tIdx));
    }
    for(Idx tIdx=Protein::p().Length()/2-1;tIdx<Protein::p().Span();++tIdx)
    {
        dirs[tIdx]=Lattice::l().Direction(c.Coord(tIdx+1)-c.Coord(tIdx+1-Protein::p().Length()/2));
    }
}

void decode(individual const dirs, Change &chg)
{
    chg.reset();
    Point tPoint=0;
    chg.add(0,tPoint);

    for(Idx tIdx=0;tIdx<Protein::p().Length()/2-1;++tIdx)
    {

        tPoint=tPoint+Lattice::l().DirVec(dirs[tIdx]);
        chg.add(tIdx+1,tPoint);
    }
    for(Idx tIdx=Protein::p().Length()/2-1;tIdx<Protein::p().Span();++tIdx)
    {
        tPoint=chg.Destination(tIdx+1-Protein::p().Length()/2)+Lattice::l().DirVec(dirs[tIdx]);
        chg.add(tIdx+1,tPoint);
    }


}
bool crossover(block1<Dir,kmm> const dirs1, block1<Dir,kmm> const dirs2,Idx pos,Change &chg)
{
    // if valid returns a Change

    //the crossover position pos is the half of the actual to take from both
    //backbone and sidechains depending on our representation
    hset<Point,xmmh> tPoints(Protein::p().Length());

    chg.reset();
    Point tPoint=0;
    chg.add(0,tPoint);
    tPoints.insertBll(tPoint);


    for(Idx tIdx=0;tIdx<Protein::p().Length()/2-1;++tIdx)
    {
        if(tIdx<pos)
            tPoint=tPoint+Lattice::l().DirVec(dirs1[tIdx]);
        else
            tPoint=tPoint+Lattice::l().DirVec(dirs2[tIdx]);
        if(tPoints.findBll(tPoint))
            return false;
        chg.add(tIdx+1,tPoint);
        tPoints.insertBll(tPoint);
    }

    for(Idx tIdx=Protein::p().Length()/2-1;tIdx<Protein::p().Span();++tIdx)
    {
        if(tIdx/2<pos)
            tPoint=chg.Destination(tIdx+1-Protein::p().Length()/2)+Lattice::l().DirVec(dirs1[tIdx]);
        else
            tPoint=chg.Destination(tIdx+1-Protein::p().Length()/2)+Lattice::l().DirVec(dirs2[tIdx]);
        if(tPoints.findBll(tPoint))
            return false;
        chg.add(tIdx+1,tPoint);
        tPoints.insertBll(tPoint);
    }
    return true;

}
void memeticSearch(int argC, char* argv[])
{
   if(argC<3)
    {
        cout << "Not enough argument provided!" << endl;
        cout << "Usage:" << argv[0] << " <sequence> <timeout>" << endl;
        return;
    }
    Int timeOut = parseN(argv[2]);
    Rnd theRnd;//(1367176523);
    cout << "seed: " <<theRnd.Seed()<<endl;
    Space::Dimen(3);
    CCLattice fcc;
    Lattice::l(fcc);

    HpEnergy hp;
    Energy::e(hp);

    // since this is a hp string, we have to copy all p in the beginning of the string
    // or take a double string as input

    // first one is convenient
    // how ever here the protein length is misleading, its actually conformation length

    // argv[1] == contains the sequence in HP

    int lengthHP=strlen(argv[1]);

    char fullSide[1000];
    // now copy the string
    int i=0;
    for(;i<lengthHP*2;i++)
    {
        if(i<lengthHP)
            fullSide[i]='P';
        else
            fullSide[i]=argv[1][i-lengthHP];
    }
    fullSide[i]=0;
    cout << fullSide << endl;

    Protein p(2,fullSide);
    Protein::p(p);


    LNSMoveSC lnsSC;
    DiagonalMoveSC diaSC;
    BackBonePullMove bbPull;

    Conf c(theRnd);
    c.setSideChain(true);



	EnergyObj hpsc(c,true);

    HCoreDistObj hcd(c,true);
    FreeFtCns free(c,true);

    SideChainCns scs(c,false);

    Sys const & tSys = Sys::refc(c.SysHdl());


    block1<Prm,nmm> prmArr(3);
    prmArr.insert(Prm(Tf,free.FuncHdl(),HintMinAsHintMin));
    prmArr.insert(Prm(Tf,hcd.FuncHdl(),HintMinAsHintMin));
    prmArr.insert(Prm(Tf,hpsc.FuncHdl(),HintMinAsHintMin));

    Prm hintSum=Prm(TermFunc,SumXiHFcMiHn::def(Xm|HintMin, tSys.SysHdl, prmArr.items(), prmArr.itemCount()));

    Prm HintHeapEn = Prm(Tf, Sv2TabuMaxHeapHiFrHi::def(tSys.SysHdl, hintSum));
    //Prm HintHeapEn = Prm(Tf, Sv2TabuMaxHeapHiFrHi::def(tSys.SysHdl, hintSum));
    //Prm HintHeapEn = Prm(Tf, Sv2TabuMaxHeapHiFrHi::def(tSys.SysHdl, Prm(Tf,free.FuncHdl())));
	Hdl const VarSelcHdl = RankedHintVar1Sp::def( tSys.SysHdl, HintHeapEn);



    /// START OF INITIALIZATION /////////////

    Idx populationSize=50;



    block1<individual,xmm> population(populationSize);
    block1<Int,xmm> fitness(populationSize);

    RandomValidSC rvsc(theRnd);
    rvsc.compute(c);
    c.initialise(rvsc.FullChange()); // dummy initialize


    individual globalBest;
    Int globalBestFitness=0;

    for(Idx tIdx=0;tIdx<populationSize;++tIdx)
    {
        individual temp(Protein::p().Span());
        // populate temp with initializations
        rvsc.compute(c);
        c.execute(rvsc.FullChange()); // dummy initialize
       c.writeHPSCToFile("data.cml");
    exit(0);
        encode(temp,c);
        population.insertMem(temp);
        fitness.insertMem(hpsc.ExecMetric());
        if(globalBestFitness>hpsc.ExecMetric())
        {
            globalBest=temp;
            globalBestFitness=hpsc.ExecMetric();
        }
    }
    cout << "global best: "<<globalBestFitness<< endl;


    // tabulist
    block1<Idx,kmm> tabuList(Protein::p().Length());
    Idx tabuLength = uniform(theRnd,4,(Z)Protein::p().Span()/8);

    typedef tuple2a<Int,Idx> mytuple;
    priority_queue<mytuple,vector<mytuple>,greater_equal<mytuple> > pQueue; // using default comparator



    TimeUtl timeUtl;
    R timeInit= timeUtl.getTime();
    while(timeUtl.getTime()-timeInit <= timeOut)
    {

        for(Idx tIdxOuter=0;tIdxOuter<populationSize;++tIdxOuter)
        {
            // for each of the ividual in population, apply operator and populate
            // alter population

            //cout << "Individual # " << tIdxOuter << endl;

            individual temp=population[tIdxOuter];
            Change tChange;
            decode(temp,tChange);
            c.execute(tChange);

            // for each of these conformations subsequent local search
            // a tabulist
            // initialize tabu
            for(Idx tIdx1=0;tIdx1 < Protein::p().Length(); ++tIdx1)
            tabuList[tIdx1]=0;

            platypus::Pos tSelcPos = 0;
            Hdl selectedVar=0;
            Idx proteinLength = Protein::p().Length();
            block1<platypus::Pos,xmm> curPositions;
            block1<Change,xmm> partialChanges;
            Idx nonImp=0;

            /// local search loop //////////////
            while(nonImp<1000)
            {
                Idx stagLoop=0;
                curPositions.clear();
                partialChanges.clear();
                lnsSC.reset();
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
                // out of this loop we either get one tSelcPos
                // or stagLoop is stuck
                if(stagLoop>=proteinLength)
                {
                    break; // stagnation // or random walk??
                }


                // operator?
                Idx windowType = uniform(theRnd,1,2);
                Idx moveType = uniform(theRnd,0,1);


                if(moveType==0 && (tSelcPos>0 && tSelcPos <Protein::p().Length()/2-1))
                {

                    // use Pull move

                    curPositions.insertMem(tSelcPos);
                    bbPull.compute(c,tSelcPos);
                    for(Idx tIdx2 = 0; tIdx2 < bbPull.PartChanges().itemCount(); ++tIdx2)
                    {
                        partialChanges.insertMem(bbPull.PartChanges()[tIdx2]);
                    }

                }
                else {
                    Idx windowSize=uniform(theRnd,1,3);

                    if(windowSize==1 && tSelcPos < Protein::p().Length()/2)
                    {
                        // diagonal
                        curPositions.insertMem(tSelcPos);
                        if(!c.isPosSideChain(tSelcPos))
                        {
                            curPositions.insertMem(tSelcPos+Protein::p().Length()/2);
                        }
                        diaSC.compute(c,tSelcPos);
                        for(Idx tIdx = 0; tIdx < diaSC.PartChanges().itemCount(); ++tIdx)
                        {
                            partialChanges.insertMem(diaSC.PartChanges()[tIdx]);
                        }
                    }
                    else
                    if(windowType==1 && tSelcPos>1+windowSize/2 && Protein::p().Span()-tSelcPos>1+windowSize/2)
                    {
                        // consecutive window
                        for(platypus::Pos tPos = 0; tPos < windowSize ; ++tPos)
                        {
                            lnsSC.addPosition(tSelcPos-windowSize/2+tPos);
                            curPositions.insertMem(tSelcPos-windowSize/2+tPos);
                            if(!c.isPosSideChain(tSelcPos-windowSize/2+tPos))
                            {
                                lnsSC.addPosition(tSelcPos-windowSize/2+tPos+Protein::p().Length()/2);
                                curPositions.insertMem(tSelcPos-windowSize/2+tPos+Protein::p().Length()/2);
                            }
                        }
            }
            else
            {
                // non-consecutive random window
                for(platypus::Pos tPos = 0; tPos < windowSize ; ++tPos)
                {
                    tSelcPos = uniform(theRnd, 1 , (Z)Protein::p().Span()-1);
                    if(!lnsSC.alreadyAdded(tSelcPos) && (tabuList[tSelcPos] < tSys.ExecClk()))
                    {
                        lnsSC.addPosition(tSelcPos);
                        curPositions.insertMem(tSelcPos);
                        //cout << tSelcPos << " " ;
                        if(!c.isPosSideChain(tSelcPos)&&!lnsSC.alreadyAdded(tSelcPos+Protein::p().Length()/2))
                        {
                            lnsSC.addPosition(tSelcPos+Protein::p().Length()/2);
                            curPositions.insertMem(tSelcPos+Protein::p().Length()/2);
                        }
                    }
                    else --tPos;

                }
            }
            lnsSC.compute(c, tSelcPos);
            for(Idx tIdx = 0; tIdx < lnsSC.PartChanges().itemCount(); ++tIdx)
            {
                    partialChanges.insertMem(lnsSC.PartChanges()[tIdx]);
            }


            }

            if(partialChanges.itemCount()==0)continue;



            //tMinEng = MaxInt;//, tSolCount = 0;
			//tMinIdx = MaxIdx;
            Idx heuristics = uniform(theRnd,0,3);
			for(Idx tIdx = 0; tIdx < partialChanges.itemCount(); ++tIdx)
			{
                c.simulate(partialChanges[tIdx]);

                if(heuristics<=1)
                {
                    Term::performSimulIncr(Func::ptrm(c.SysHdl(), hcd.FuncHdl()));
                    mytuple t(hcd.SimulMetric(),tIdx);
                    pQueue.push(t);
                    /*if (tMinIdx == MaxIdx || tMinEng >= hcd.SimulMetric())
                    {
                        tMinEng =hcd.SimulMetric();
                        tMinIdx = tIdx;
                    }*/
                }
                else
                {
                    Term::performSimulIncr(Func::ptrm(c.SysHdl(), hpsc.FuncHdl()));
                    mytuple t(hpsc.SimulMetric(),tIdx);
                    pQueue.push(t);
                    /*if (tMinIdx == MaxIdx || tMinEng >= hpsc.SimulMetric())
                    {
                        tMinEng =hpsc.SimulMetric();
                        tMinIdx = tIdx;
                    }*/
                }



            }
            //if (tMinIdx == MaxIdx) continue;
            Idx selectedId=1;
            if(!pQueue.empty())
            {
                mytuple top=pQueue.top();
                block1<mytuple,xmm> vector;
                while(!pQueue.empty()&&top==pQueue.top())
                {
                    vector.insertMem(pQueue.top());
                    pQueue.pop();
                }
                Idx selRank=uniform(theRnd,(N)0,vector.itemCount()-1);
                //cout << "selRank " << selRank << " selectedId: " << selectedId <<" ";
                selectedId=vector[selRank].Second;
            }
            else
            {
                for(Idx tIdx =0; tIdx < curPositions.itemCount();++tIdx)
                    {
                        tabuList[tIdx]=tSys.ExecClk()+tabuLength;
                    }

                continue;
            }
            //cout << "selected id:"<<selectedId  << "partial changes" << partChanges.itemCount()<< endl;
            //cout << partialChanges[selectedId].size()<< endl;
            c.execute(partialChanges[selectedId]);
            while(!pQueue.empty())pQueue.pop();
            //now after each of this execution update the global best
            if(globalBestFitness>hpsc.ExecMetric())
            {
                individual temp1(proteinLength-1);
                encode(temp1,c);
                globalBest=temp1;
                globalBestFitness=hpsc.ExecMetric();
                cout << "global best: "<<globalBestFitness<<endl;
            }
            else
            {
                nonImp++;
            }

            }
            /// local search loop //////////////
            /// after local search is done we need to save the latest conformation into the population
            ///Dont need the previous population
            /// discard them
            encode(temp,c);
            population[tIdxOuter]=temp;
            fitness[tIdxOuter]=hpsc.ExecMetric();
            /// then crossover

            block1<individual,xmm> alterPopulation;
            block1<Int,xmm> alterFitness;

            Idx tIdxAlt=0;
            while(tIdxAlt<50)
            {
                Idx t1=uniform(theRnd,(Z)0,(Z)populationSize-1);
                Idx t2=uniform(theRnd,(Z)0,(Z)populationSize-1);
                if(t1==t2) continue;

                Idx pos=uniform(theRnd,(Z)1,(Z)proteinLength/2-2);
                Change tempChange;
                //cout<<t1<< " "<< t2 << " "<<pos <<endl;
                bool cross = crossover(population[t1],population[t2],pos,tempChange);
                if(cross==false) continue;
                c.execute(tempChange);
                individual tempI(proteinLength-1);
                encode(tempI,c);
                alterPopulation.insertMem(temp);
                alterFitness.insertMem(hpsc.ExecMetric());
                tIdxAlt++;
                if(globalBestFitness>hpsc.ExecMetric())
                {
                    individual temp1(proteinLength-1);
                    encode(temp1,c);
                    globalBest=temp1;
                    globalBestFitness=hpsc.ExecMetric();
                    cout << "global best: "<<globalBestFitness<<endl;
                }
            }
            /// now we have got two populations; need to keep one
            while(tIdxAlt<populationSize*2)
            {
                alterPopulation.insertMem(population[tIdxAlt-populationSize]);
                alterFitness.insertMem(fitness[tIdxAlt-populationSize]);
                tIdxAlt++;
            }
            /// now sort

            for(Idx i=0;i<populationSize*2-1;++i)
                for(Idx j=i+1;j<populationSize*2;++j)
                {
                    if(alterFitness[i]>alterFitness[j])
                    {
                        Int tF=alterFitness[i];
                        alterFitness[i]=alterFitness[j];
                        alterFitness[j]=tF;

                        individual tI=alterPopulation[i];
                        alterPopulation[i]=alterPopulation[j];
                        alterPopulation[j]=tI;
                    }
                }
            /// time to select - just take best 50
            population[0]=globalBest;
            fitness[0]=globalBestFitness;

            for(Idx i=1;i<populationSize;++i)
            {
                population[i]=alterPopulation[i];
                fitness[i]=alterFitness[i];
            }


        }


    }
    cout <<"globalBest: "<< globalBestFitness << endl;


}

void scOperator(int argC, char* argv[])
{
    if(argC<3)
    {
        cout << "Not enough argument provided!" << endl;
        cout << "Usage:" << argv[0] << " <sequence> <timeout>" << endl;
        return;
    }
    Int timeOut = parseN(argv[2]);
    Rnd theRnd;//(1367176523);
    cout << "seed: " <<theRnd.Seed()<<endl;
    Space::Dimen(3);
    CCLattice fcc;
    Lattice::l(fcc);

    HpEnergy hp;
    Energy::e(hp);

    // since this is a hp string, we have to copy all p in the beginning of the string
    // or take a double string as input

    // first one is convenient
    // how ever here the protein length is misleading, its actually conformation length

    // argv[1] == contains the sequence in HP

    int lengthHP=strlen(argv[1]);

    char fullSide[1000];
    // now copy the string
    int i=0;
    for(;i<lengthHP*2;i++)
    {
        if(i<lengthHP)
            fullSide[i]='P';
        else
            fullSide[i]=argv[1][i-lengthHP];
    }
    fullSide[i]=0;
    cout << fullSide << endl;

    Protein p(2,fullSide);
    Protein::p(p);


    LNSMoveSC lnsSC;
    DiagonalMoveSC diaSC;
    BackBonePullMove bbPull;

    Conf c(theRnd);
    c.setSideChain(true);



	EnergyObj hpsc(c,true);

    HCoreDistObj hcd(c,true);
    FreeFtCns free(c,true);

	RandomValidSC rvsc(theRnd);
	rvsc.compute(c);


	Sys const & tSys = Sys::refc(c.SysHdl());




    block1<Prm,nmm> prmArr(3);
    prmArr.insert(Prm(Tf,free.FuncHdl(),HintMinAsHintMin));
    prmArr.insert(Prm(Tf,hcd.FuncHdl(),HintMinAsHintMin));
    prmArr.insert(Prm(Tf,hpsc.FuncHdl(),HintMinAsHintMin));

    Prm hintSum=Prm(TermFunc,SumXiHFcMiHn::def(Xm|HintMin, tSys.SysHdl, prmArr.items(), prmArr.itemCount()));

    Prm HintHeapEn = Prm(Tf, Sv2TabuMaxHeapHiFrHi::def(tSys.SysHdl, hintSum));
    //Prm HintHeapEn = Prm(Tf, Sv2TabuMaxHeapHiFrHi::def(tSys.SysHdl, hintSum));
    //Prm HintHeapEn = Prm(Tf, Sv2TabuMaxHeapHiFrHi::def(tSys.SysHdl, Prm(Tf,free.FuncHdl())));
	Hdl const VarSelcHdl = RankedHintVar1Sp::def( tSys.SysHdl, HintHeapEn);


    c.initialise(rvsc.FullChange());

    Int globalMin= hpsc.ExecMetric();
    cout << "global best: "<<globalMin << endl;


// tabulist
    block1<Idx,kmm> tabuList(Protein::p().Length());
    Idx tabuLength = uniform(theRnd,4,(Z)Protein::p().Span()/8);

    typedef tuple2a<Int,Idx> mytuple;
    priority_queue<mytuple,vector<mytuple>,greater_equal<mytuple> > pQueue; // using default comparator

    Change bestPoints;
     block1<Change,xmm> elitSet;

    Idx nonImp =0;
    TimeUtl timeUtl;
    R timeInit= timeUtl.getTime();
 for(Idx tIdx1=0;tIdx1 < Protein::p().Length(); ++tIdx1)
            tabuList[tIdx1]=0;

            platypus::Pos tSelcPos = 0;
            Hdl selectedVar=0;
            Idx proteinLength = Protein::p().Length();
            block1<platypus::Pos,xmm> curPositions;
            block1<Change,xmm> partialChanges;

            Int initWalkPV=15;
            Int walkPV=initWalkPV; // 5% walk probability
            R mFactor=1.2;
            Idx initMaxStable = 1000;
            Idx maxStable=initMaxStable;

            while(timeUtl.getTime()-timeInit <= timeOut)
            {
                Idx stagLoop=0;
                curPositions.clear();
                partialChanges.clear();
                lnsSC.reset();
                Int wp=uniform(theRnd,0,100);
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
                // out of this loop we either get one tSelcPos
                // or stagLoop is stuck
                if(stagLoop>=proteinLength)
                {
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


                // operator?
                Idx windowType = uniform(theRnd,1,2);
                Idx moveType = uniform(theRnd,0,1);


                if(moveType==0 && (tSelcPos>0 && tSelcPos <Protein::p().Length()/2-1))
                {

                    // use Pull move

                    curPositions.insertMem(tSelcPos);
                    bbPull.compute(c,tSelcPos);
                    for(Idx tIdx2 = 0; tIdx2 < bbPull.PartChanges().itemCount(); ++tIdx2)
                    {
                        partialChanges.insertMem(bbPull.PartChanges()[tIdx2]);
                    }

                }
                else {
                    Idx windowSize=uniform(theRnd,1,3);

                    if(windowSize==1 && tSelcPos < Protein::p().Length()/2)
                    {
                        // diagonal
                        curPositions.insertMem(tSelcPos);
                        if(!c.isPosSideChain(tSelcPos))
                        {
                            curPositions.insertMem(tSelcPos+Protein::p().Length()/2);
                        }
                        diaSC.compute(c,tSelcPos);
                        for(Idx tIdx = 0; tIdx < diaSC.PartChanges().itemCount(); ++tIdx)
                        {
                            partialChanges.insertMem(diaSC.PartChanges()[tIdx]);
                        }
                    }
                    else
                    if(windowType==1 && tSelcPos>1+windowSize/2 && Protein::p().Span()-tSelcPos>1+windowSize/2)
                    {
                        // consecutive window
                        for(platypus::Pos tPos = 0; tPos < windowSize ; ++tPos)
                        {
                            lnsSC.addPosition(tSelcPos-windowSize/2+tPos);
                            curPositions.insertMem(tSelcPos-windowSize/2+tPos);
                            if(!c.isPosSideChain(tSelcPos-windowSize/2+tPos))
                            {
                                lnsSC.addPosition(tSelcPos-windowSize/2+tPos+Protein::p().Length()/2);
                                curPositions.insertMem(tSelcPos-windowSize/2+tPos+Protein::p().Length()/2);
                            }
                        }
            }
            else
            {
                // non-consecutive random window
                for(platypus::Pos tPos = 0; tPos < windowSize ; ++tPos)
                {
                    tSelcPos = uniform(theRnd, 1 , (Z)Protein::p().Span()-1);
                    if(!lnsSC.alreadyAdded(tSelcPos) && (tabuList[tSelcPos] < tSys.ExecClk()))
                    {
                        lnsSC.addPosition(tSelcPos);
                        curPositions.insertMem(tSelcPos);
                        //cout << tSelcPos << " " ;
                        if(!c.isPosSideChain(tSelcPos)&&!lnsSC.alreadyAdded(tSelcPos+Protein::p().Length()/2))
                        {
                            lnsSC.addPosition(tSelcPos+Protein::p().Length()/2);
                            curPositions.insertMem(tSelcPos+Protein::p().Length()/2);
                        }
                    }
                    else --tPos;

                }
            }
            lnsSC.compute(c, tSelcPos);
            for(Idx tIdx = 0; tIdx < lnsSC.PartChanges().itemCount(); ++tIdx)
            {
                    partialChanges.insertMem(lnsSC.PartChanges()[tIdx]);
            }


            }

            if(partialChanges.itemCount()==0)continue;



            //tMinEng = MaxInt;//, tSolCount = 0;
			//tMinIdx = MaxIdx;
            Idx heuristics = uniform(theRnd,0,3);
			for(Idx tIdx = 0; tIdx < partialChanges.itemCount(); ++tIdx)
			{
                c.simulate(partialChanges[tIdx]);

                if(heuristics<=1)
                {
                    Term::performSimulIncr(Func::ptrm(c.SysHdl(), hcd.FuncHdl()));
                    mytuple t(hcd.SimulMetric(),tIdx);
                    pQueue.push(t);
                    /*if (tMinIdx == MaxIdx || tMinEng >= hcd.SimulMetric())
                    {
                        tMinEng =hcd.SimulMetric();
                        tMinIdx = tIdx;
                    }*/
                }

                else
                {
                    Term::performSimulIncr(Func::ptrm(c.SysHdl(), hpsc.FuncHdl()));
                    mytuple t(hpsc.SimulMetric(),tIdx);
                    pQueue.push(t);
                    /*if (tMinIdx == MaxIdx || tMinEng >= hpsc.SimulMetric())
                    {
                        tMinEng =hpsc.SimulMetric();
                        tMinIdx = tIdx;
                    }*/
                }



            }
            //if (tMinIdx == MaxIdx) continue;
            Idx selectedId=1;
            if(!pQueue.empty())
            {
                mytuple top=pQueue.top();
                block1<mytuple,xmm> vector;
                while(!pQueue.empty()&&top==pQueue.top())
                {
                    vector.insertMem(pQueue.top());
                    pQueue.pop();
                }
                Idx selRank=uniform(theRnd,(N)0,vector.itemCount()-1);
                //cout << "selRank " << selRank << " selectedId: " << selectedId <<" ";
                selectedId=vector[selRank].Second;
            }
            else
            {
                for(Idx tIdx =0; tIdx < curPositions.itemCount();++tIdx)
                    {
                        tabuList[tIdx]=tSys.ExecClk()+tabuLength;
                    }

                continue;
            }
            //cout << "selected id:"<<selectedId  << "partial changes" << partChanges.itemCount()<< endl;
            //cout << partialChanges[selectedId].size()<< endl;
            c.execute(partialChanges[selectedId]);
            while(!pQueue.empty())pQueue.pop();
            //now after each of this execution update the global best
            tabuList[tSelcPos]=tSys.ExecClk()+tabuLength;
            if(hpsc.ExecMetric() < globalMin)
            {
                globalMin = hpsc.ExecMetric();
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
                cout << "restart: " << hpsc.ExecMetric()<<" "<<tSys.ExecClk()<<endl;
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

        }
            /// local search loop //////////////


}
void sidechainSearch(int argC, char* argv[])
{

    if(argC<3)
    {
        cout << "Not enough argument provided!" << endl;
        cout << "Usage:" << argv[0] << " <sequence> <timeout>" << endl;
        return;
    }
    Int timeOut = parseN(argv[2]);
    Rnd theRnd;//(1367176523);
    cout << "seed: " <<theRnd.Seed()<<endl;
    Space::Dimen(3);
    CCLattice fcc;
    Lattice::l(fcc);

    HpEnergy hp;
    Energy::e(hp);

    // since this is a hp string, we have to copy all p in the beginning of the string
    // or take a double string as input

    // first one is convenient
    // how ever here the protein length is misleading, its actually conformation length

    // argv[1] == contains the sequence in HP

    int lengthHP=strlen(argv[1]);

    char fullSide[1000];
    // now copy the string
    int i=0;
    for(;i<lengthHP*2;i++)
    {
        if(i<lengthHP)
            fullSide[i]='P';
        else
            fullSide[i]=argv[1][i-lengthHP];
    }
    fullSide[i]=0;
    cout << fullSide << endl;

    Protein p(2,fullSide);
    Protein::p(p);


    LNSMoveSC lnsSC;
    DiagonalMoveSC diaSC;
    BackBonePullMove bbPull;

    Conf c(theRnd);
    c.setSideChain(true);



	EnergyObj hpsc(c,true);

    HCoreDistObj hcd(c,false);
    HhDistObj hhd(c,false);


    SideChainCns scs(c,false);

	RandomValidSC rvsc(theRnd);
	rvsc.compute(c);


	Sys const & tSys = Sys::refc(c.SysHdl());
	c.initialise(rvsc.FullChange());



    Int globalMin= hpsc.ExecMetric();
    cout << "global best: "<<globalMin << " "<<scs.ExecMetric()<<endl;
    //c.writeHPSCToFile("data.cml");

    TimeUtl timeUtl;
    R timeInit= timeUtl.getTime();

    // variables for tabu
        block1<Idx,kmm> tabuList(Protein::p().Length());
        Idx tabuLength = uniform(theRnd,4,(Z)Protein::p().Span()/16);
        // initialize tabu
        for(Idx tIdx=0;tIdx < Protein::p().Length(); ++tIdx)
            tabuList[tIdx]=0;

        block1<platypus::Pos,xmm> curPositions;

        Idx nonImp =0;
        Idx initMaxStable = 1000;
        Idx maxStable=initMaxStable;
        Idx initWindowSize = 1;
        Idx windowSize=initWindowSize;
        R mFactor=1.2;
        //Int     iteration=0;

        block1<Change,xmm> partialChanges;
        typedef tuple2a<Int,Idx> mytuple;
        priority_queue<mytuple,vector<mytuple>,greater_equal<mytuple> > pQueue;

        while(timeUtl.getTime()-timeInit <= timeOut/*++iteration<timeOut*/)
        {

             partialChanges.clear();

            for(Idx tIdx=0;tIdx<Protein::p().Length();++tIdx)
            {
                if(tIdx==Protein::p().Length()/2-1) continue;
                Point p1=c.Coord(tIdx);
                Idx pos2=tIdx<Protein::p().Length()/2?tIdx+1:tIdx-Protein::p().Length()/2;
                Point p2=c.Coord(pos2);
                if(!Lattice::l().areNeighbors(p1,p2))
                {
                    cout <<"Exiting!"<< tIdx << " " << pos2 << endl;
                    exit(0);
                }
            }

             //cout << "global best: "<<globalMin << " "<<scs.ExecMetric()<<endl;
            //windowSize=uniform(theRnd,4,(Idx)10);

            //Int tMinEng; Idx tMinIdx;

            curPositions.clear();

            lnsSC.reset();

            // add positions to the move
            platypus::Pos tSelcPos=0;

            Idx windowType = uniform(theRnd,1,2);
            Idx moveType = uniform(theRnd,0,1);


            if(moveType==0)
            {

            // use Diagonal move

                tSelcPos=uniform(theRnd,0,(Z)Protein::p().Length()/2-1);
                if(tabuList[tSelcPos]>=tSys.ExecClk())
                    continue;
                curPositions.insertMem(tSelcPos);
                bbPull.compute(c,tSelcPos);
                for(Idx tIdx = 0; tIdx < bbPull.PartChanges().itemCount(); ++tIdx)
                {
                    partialChanges.insertMem(bbPull.PartChanges()[tIdx]);
                }

            }
            else {

            if(windowSize==1)
            {
                // diagonal
                tSelcPos=uniform(theRnd,0,(Z)Protein::p().Length()/2-1);
                if(tabuList[tSelcPos]>=tSys.ExecClk())
                    continue;
                curPositions.insertMem(tSelcPos);
                if(!c.isPosSideChain(tSelcPos))
                {
                    curPositions.insertMem(tSelcPos+Protein::p().Length()/2);
                }
                diaSC.compute(c,tSelcPos);
                for(Idx tIdx = 0; tIdx < diaSC.PartChanges().itemCount(); ++tIdx)
                {
                    partialChanges.insertMem(diaSC.PartChanges()[tIdx]);
                }
            }
            else
            if(windowType==1)
            {


            // consecutive window

                tSelcPos = uniform(theRnd, 1+windowSize/2 , Protein::p().Span()-1-windowSize/2);
                B exitFor = false;
                for(platypus::Pos tPos = 0; tPos < windowSize ; ++tPos)
                {
                    if(tabuList[tSelcPos-windowSize/2+tPos] < tSys.ExecClk())
                    {
                        lnsSC.addPosition(tSelcPos-windowSize/2+tPos);
                        curPositions.insertMem(tSelcPos-windowSize/2+tPos);
                        if(!c.isPosSideChain(tSelcPos-windowSize/2+tPos))
                        {
                            lnsSC.addPosition(tSelcPos-windowSize/2+tPos+Protein::p().Length()/2);
                            curPositions.insertMem(tSelcPos-windowSize/2+tPos+Protein::p().Length()/2);
                        }
                    }
                    else
                    {
                        exitFor=true; break;
                    }

                }
                if(exitFor) continue;
            }
            else
            {


                // non-consecutive random window


                for(platypus::Pos tPos = 0; tPos < windowSize ; ++tPos)
                {
                    tSelcPos = uniform(theRnd, 1 , (Z)Protein::p().Span()-1);
                    if(!lnsSC.alreadyAdded(tSelcPos) && (tabuList[tSelcPos] < tSys.ExecClk()))
                    {
                        lnsSC.addPosition(tSelcPos);
                        curPositions.insertMem(tSelcPos);
                        //cout << tSelcPos << " " ;
                        if(!c.isPosSideChain(tSelcPos)&&!lnsSC.alreadyAdded(tSelcPos+Protein::p().Length()/2))
                        {
                            lnsSC.addPosition(tSelcPos+Protein::p().Length()/2);
                            curPositions.insertMem(tSelcPos+Protein::p().Length()/2);
                        }
                    }
                    else --tPos;

                }
            }
            /*cout <<iteration<< " type:" << windowType << " <<";
            for(Idx tIdx=0;tIdx<curPositions.itemCount();++tIdx)
            {
                cout << curPositions[tIdx] << " ";
            }
            cout<<endl;*/

            lnsSC.compute(c, tSelcPos);
            //cout << " window size:" << windowSize << " part Size : " << dm1.PartChanges().size()<<endl;
            for(Idx tIdx = 0; tIdx < lnsSC.PartChanges().itemCount(); ++tIdx)
            {
                    partialChanges.insertMem(lnsSC.PartChanges()[tIdx]);
            }


            }

            if(partialChanges.itemCount()==0)continue;



            //tMinEng = MaxInt;//, tSolCount = 0;
			//tMinIdx = MaxIdx;
            Idx heuristics = uniform(theRnd,1,3);
			for(Idx tIdx = 0; tIdx < partialChanges.itemCount(); ++tIdx)
			{
                c.simulate(partialChanges[tIdx]);

                if(heuristics==1)
                {
                    Term::performSimulIncr(Func::ptrm(c.SysHdl(), hcd.FuncHdl()));
                    mytuple t(hcd.SimulMetric(),tIdx);
                    pQueue.push(t);
                    /*if (tMinIdx == MaxIdx || tMinEng >= hcd.SimulMetric())
                    {
                        tMinEng =hcd.SimulMetric();
                        tMinIdx = tIdx;
                    }*/
                }
                else if(heuristics==2)
                {
                    Term::performSimulIncr(Func::ptrm(c.SysHdl(), hhd.FuncHdl()));
                    mytuple t(hhd.SimulMetric(),tIdx);
                    pQueue.push(t);
                    /*if (tMinIdx == MaxIdx || tMinEng >= hhd.SimulMetric())
                    {
                        tMinEng =hhd.SimulMetric();
                        tMinIdx = tIdx;
                    }*/
                }
                else
                {
                    Term::performSimulIncr(Func::ptrm(c.SysHdl(), hpsc.FuncHdl()));
                    mytuple t(hpsc.SimulMetric(),tIdx);
                    pQueue.push(t);
                    /*if (tMinIdx == MaxIdx || tMinEng >= hpsc.SimulMetric())
                    {
                        tMinEng =hpsc.SimulMetric();
                        tMinIdx = tIdx;
                    }*/
                }



            }
            //if (tMinIdx == MaxIdx) continue;
            Idx selectedId=1;
            if(!pQueue.empty())
            {
                mytuple top=pQueue.top();
                block1<mytuple,xmm> vector;
                while(!pQueue.empty()&&top==pQueue.top())
                {
                    vector.insertMem(pQueue.top());
                    pQueue.pop();
                }
                Idx selRank=uniform(theRnd,(N)0,vector.itemCount()-1);
                //cout << "selRank " << selRank << " selectedId: " << selectedId <<" ";
                selectedId=vector[selRank].Second;
            }
            else
            {
                for(Idx tIdx =0; tIdx < curPositions.itemCount();++tIdx)
                    {
                        tabuList[tIdx]=tSys.ExecClk()+tabuLength;
                    }

                nonImp++;
                continue;
            }
            //cout << "selected id:"<<selectedId  << "partial changes" << partChanges.itemCount()<< endl;
            //cout << partialChanges[selectedId].size()<< endl;
            c.execute(partialChanges[selectedId]);
            while(!pQueue.empty())pQueue.pop();


            //c.writeHPSCToFile("data.cml");

            if(hpsc.ExecMetric() < globalMin)
            {
                //c.writeHPSCToFile("dataB.cml");

                globalMin = hpsc.ExecMetric();
                windowSize = initWindowSize;
                nonImp = 0;
                cout << "new global best:" << globalMin << endl;
                maxStable = initMaxStable;
                for(Idx tIdx=0;tIdx < Protein::p().Length(); ++tIdx)
                    tabuList[tIdx]=0;
            }
            else
            {
                nonImp++;
            }

            if(nonImp>=maxStable)
            {
                cout << "Stagnation!" << endl;
                nonImp=0;
                windowSize++;
                maxStable=maxStable*mFactor;

			//if(tSys.ExecClk()%1000==0)
			{
                cout << tSys.ExecClk() <<" "<< timeUtl.getTime()-timeInit<< " " <<hhd.ExecMetric() << " " << hpsc.ExecMetric() << " "<< globalMin <<endl;
			}
			}
            // Update tabu list
            for(Idx tIdx =0; tIdx < curPositions.itemCount();++tIdx)
            {
                tabuList[tIdx]=tSys.ExecClk()+tabuLength;
            }


        }
        cout<< globalMin << endl;


    exit(0);



}
