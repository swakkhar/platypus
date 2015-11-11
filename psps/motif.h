#include "pspl/idx.hh"
#include<limits>
using namespace kangaroo;
using namespace platypus;
using namespace std;
typedef block1<Dir,kmm> individual;


void encodeM(individual & dirs,Conf const &theConf)
{
    Vector absVec;
    for(Cmp tCmp=0;tCmp<Protein::p().Span();++tCmp)
    {
        absVec= theConf.Coord(tCmp+1)-theConf.Coord(tCmp);
        // now we have got the absVec
        dirs[tCmp]=Lattice::l().Direction(absVec);
        //cout << absVec.Comp(0) << " "<<absVec.Comp(1) << " "<<absVec.Comp(2) << " "<<(int)mAbsDirs[tCmp]<<endl;
    }
}
void decodeM(individual const dirs, Change &chg)
{
    chg.reset();
    Point tPoint=0;
    chg.add(0,tPoint);
    for(Idx tIdx=0;tIdx<Protein::p().Span();++tIdx)
    {
        tPoint=tPoint+Lattice::l().DirVec(dirs[tIdx]);
        chg.add(tIdx+1,tPoint);
    }
}

bool crossoverM(block1<Dir,kmm> const dirs1, block1<Dir,kmm> const dirs2,Idx pos,Change &chg)
{
    chg.reset();
    Point tPoint=0;
    chg.add(0,tPoint);

    hset<Point,xmmh> tPoints(Protein::p().Length());

    tPoints.insertBll(tPoint);

    for(Idx tIdx=0;tIdx<Protein::p().Span();++tIdx)
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
    return true;

}


void motifSearch(int argC, char* argv[])
{

    if(argC<=4)
    {
        cout << "REMARK Not enough argument provided!" << endl;
        cout << "REMARK Usage:" << argv[0] << " <sequence> <timeout>" << endl;
        return;
    }
    Int timeOut = parseN(argv[2]);
    Rnd theRnd;
    cout << "REMARK seed: " <<theRnd.Seed()<<endl;
    Space::Dimen(3);
    FccLattice fcc;
    Lattice::l(fcc);

    //BarerraEnergy br;
    MjEnergy br;
    Energy::e(br);


    Protein p(20,argv[1]);
    Protein::p(p);


     block1<RPoint,xmm> allPoints;
        int start=0;
        int end=0;

    PDBReader::readPDBFile(argv[3],allPoints,start,end);

    Conf c(theRnd);
    c.setMotif(true);

    if(strlen(argv[4])!=strlen(argv[1]))
    {
        cout << "REMARK Encoded string and Protein sequence length mismatch! " << strlen(argv[4])<<" "<<strlen(argv[1]) <<endl;
        return;
    }

    c.insertMotifs(argv[4]);

    EnergyObj bre(c,false);



    RotationMove rotationM;
    MotifPullMove pullM;
    LNSMove lnsM;


/// START OF INITIALIZATION /////////////

    Idx populationSize=50;



    block1<individual,xmm> population(populationSize);
    block1<Int,xmm> fitness(populationSize);

    MotifInit rvsc(theRnd);
	rvsc.compute(c);

    ACoreDistObj aco(c,false);
    FreeFtCns free(c,false);
    RMSDObj rmsd(c,false,allPoints);

    Sys const & tSys = Sys::refc(c.SysHdl());
    c.initialise(rvsc.FullChange()); // dummy initialize






    individual globalBest;
    individual globalBestFine;
    Int globalBestFitness=0;
    R globalBestRMSD=0;
    R globalBestRMSDFine=std::numeric_limits<double>::max();

    Idx proteinLength = Protein::p().Length();

    for(Idx tIdx=0;tIdx<populationSize;++tIdx)
    {
        individual temp(Protein::p().Span());
        // populate temp with initializations
        rvsc.compute(c);
        c.execute(rvsc.FullChange()); // dummy initialize
        encodeM(temp,c);
        population.insertMem(temp);
        fitness.insertMem(bre.ExecMetric());
        if(globalBestFitness>bre.ExecMetric())
        {
            globalBest=temp;
            globalBestFine=temp;
            globalBestFitness=bre.ExecMetric();
            globalBestRMSD=0.1*sqrt(2*rmsd.ExecMetric()/(proteinLength*(proteinLength-1)));
        }
    }
    cout << "REMARK global best:"<<globalBestFitness<<" "<<globalBestRMSD<< endl;
    //c.writeMotifToFile("data.cml");

    // tabulist
    block1<Idx,kmm> tabuList(Protein::p().Length());
    Idx tabuLength = uniform(theRnd,4,(Z)Protein::p().Span()/8);

    typedef tuple2a<Int,Idx> mytuple;
    priority_queue<mytuple,vector<mytuple>,greater_equal<mytuple> > pQueue; // using default comparator


    TimeUtl timeUtl;
    R timeInit= timeUtl.getTime();

    while(timeUtl.getTime()-timeInit <= timeOut)
    {
        individual temp(Protein::p().Span());
        for(Idx tIdxOuter=0;tIdxOuter<populationSize;++tIdxOuter)
        {
            temp=population[tIdxOuter];
            Change tChange;
            decodeM(temp,tChange);
            c.execute(tChange);

            // now the conformation c contains the individual

            // for each of these conformations subsequent local search
            // a tabulist
            // initialize tabu
            for(Idx tIdx1=0;tIdx1 < Protein::p().Length(); ++tIdx1)
            tabuList[tIdx1]=0;

            platypus::Pos tSelcPos = 0;
            //Hdl selectedVar=0;

            block1<platypus::Pos,xmm> curPositions;
            block1<Change,xmm> partialChanges;
            Idx nonImp=0;
            //cout << endl;
            /// local search loop //////////////
            Idx noMoves=0;
            while(nonImp<1)
            {

                Idx moveType=uniform(theRnd,0,1);
                /*!
                    0 = pull move
                    1 = rotation
                    2 = lns move
                */

                Idx windowType=uniform(theRnd,0,1);
                // point selection is randomized

                curPositions.clear();
                partialChanges.clear();

                lnsM.reset();

                // all clear // now select the points according to the types
                Idx stagLoop=0;
                if(moveType==0)
                {
                    // select a random point that is not in the tabu list and motif
                    while(1)
                    {
                        tSelcPos=uniform(theRnd,(Idx)1,proteinLength-2); // leaving the end points
                        if(c.isPosMotif(tSelcPos)==false && tabuList[tSelcPos]< tSys.ExecClk())
                            break;
                        stagLoop++;
                        if(stagLoop>=proteinLength)
                            break;

                    }
                    if(stagLoop>=proteinLength)
                        break; // stagnation // or random walk??

                    // else apply the move

                    curPositions.insertMem(tSelcPos);
                    pullM.compute(c,tSelcPos);
                    for(Idx tIdx2 = 0; tIdx2 < pullM.PartChanges().itemCount(); ++tIdx2)
                        partialChanges.insertMem(pullM.PartChanges()[tIdx2]);

                }
                else if(moveType==1) // rotation move
                {

                    ///Idx accSize=0;
                    while(1)
                    {
                            tSelcPos=uniform(theRnd,(Idx)1,proteinLength-2); // no need to spare the end point, only beginning
                            if(c.isPosMotif(tSelcPos+1)/*&&c.isPosMotif(tSelcPos-1)*/)
                            {
                                stagLoop++;
                                if(stagLoop>=proteinLength)
                                    break;
                            }
                            else if(tabuList[tSelcPos]>= tSys.ExecClk())
                            {
                                stagLoop++;
                                if(stagLoop>=proteinLength)
                                    break;
                            }
                            else
                            {
                                break;
                            }

                    }
                    if(stagLoop>=proteinLength)
                        break; // stagnation // or random walk??

                    rotationM.compute(c,tSelcPos);
                    for(Idx tIdx2 = 0; tIdx2 < rotationM.PartChanges().itemCount(); ++tIdx2)
                        partialChanges.insertMem(rotationM.PartChanges()[tIdx2]);


                }
                else // lns move
                {
                    Idx moveSize=uniform(theRnd,1,3);
                    Idx accSize=0;
                    if(windowType==0) //consecutive
                    {
                        while(1)
                        {
                            lnsM.reset();

                            tSelcPos=uniform(theRnd,(Idx)1+moveSize/2,proteinLength-2-moveSize/2);
                            if(tabuList[tSelcPos]>= tSys.ExecClk())
                            {
                                stagLoop++;
                                if(stagLoop>=proteinLength)
                                    break;
                                continue;
                            }
                            ///now we have to see that the window does not belong to tabu list
                            /// or motif

                            B exitFor = false;
                            for(platypus::Pos tPos = 0; tPos < moveSize ; ++tPos)
                            {
                                if(!c.isPosMotif(tPos))
                                {
                                     cout << "moveSize "<<moveSize << " "<<tSelcPos-moveSize/2+tPos << " "<<lnsM.size()<<endl;
                                    lnsM.addPosition(tSelcPos-moveSize/2+tPos);
                                    curPositions.insertMem(tSelcPos-moveSize/2+tPos);
                                }
                                else
                                {
                                    exitFor=true; break;
                                }

                            }
                            if(exitFor)
                            {
                                stagLoop=proteinLength+1;
                            }
                            break;

                        }
                    }
                    else // non-consecutive
                    {
                        while(accSize<moveSize)
                    {
                            tSelcPos=uniform(theRnd,(Idx)1,proteinLength-2); // no need to spare the end point, only beginning
                            if(c.isPosMotif(tSelcPos))
                            {
                                stagLoop++;
                                if(stagLoop>=proteinLength)
                                    break;
                            }
                            else if(tabuList[tSelcPos]>= tSys.ExecClk())
                            {
                                stagLoop++;
                                if(stagLoop>=proteinLength)
                                    break;
                            }
                            else
                            {
                                // valid position
                                lnsM.addPosition(tSelcPos);
                                curPositions.insertMem(tSelcPos);
                                accSize++;
                            }

                    }
                    }
                    if(stagLoop>=proteinLength)
                        break; // stagnation // or random walk??

                    lnsM.compute(c,tSelcPos);
                    for(Idx tIdx2 = 0; tIdx2 < lnsM.PartChanges().itemCount(); ++tIdx2)
                        partialChanges.insertMem(lnsM.PartChanges()[tIdx2]);

                }

                if(partialChanges.itemCount()==0){

                    noMoves++;
                    if(noMoves>=proteinLength)
                    nonImp++;
                    continue;
                }
                else
                {
                    noMoves=0;
                }

                /// all moves genereated
                /// now make selection according to the heuristics
                Idx heuristics = uniform(theRnd,0,3);
                for(Idx tIdx = 0; tIdx < partialChanges.itemCount(); ++tIdx)
                {
                    c.simulate(partialChanges[tIdx]);

                    if(heuristics<=1)
                    {
                        Term::performSimulIncr(Func::ptrm(c.SysHdl(), bre.FuncHdl()));
                        mytuple t(bre.SimulMetric(),tIdx);
                        pQueue.push(t);

                    }
                    else if(heuristics<=2)
                    {
                        Term::performSimulIncr(Func::ptrm(c.SysHdl(), aco.FuncHdl()));
                        mytuple t(aco.SimulMetric(),tIdx);
                        pQueue.push(t);
                    }
                    else
                    {
                        Term::performSimulIncr(Func::ptrm(c.SysHdl(), free.FuncHdl()));
                        mytuple t(free.SimulMetric(),tIdx);
                        pQueue.push(t);
                    }
                } /// end of heuristic simulation

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
                //cout<<"here:" <<tSelcPos<<" "<<moveType<<" "<<" "<<selectedId<< " "<<partialChanges.itemCount()                                                                                                                                                                                                                                                                                                                         <<endl;
                /// prin the changes
                //cout <<"start print";
                /*for(Idx ttt=0;ttt<partialChanges[selectedId].size();++ttt)
                {
                    cout<<"("<<partialChanges[selectedId].Position(ttt) << " -- ";
                    cout << partialChanges[selectedId].Destination(ttt).Comp(0) << " "
                        << partialChanges[selectedId].Destination(ttt).Comp(1) << " "
                        <<partialChanges[selectedId].Destination(ttt).Comp(2)<<" ) ";
                }
                cout << endl;*/
                c.execute(partialChanges[selectedId]);
                 //c.writeMotifToFile("data.cml");
                /// DEBUG here ///////////////////////////////////////////


                /// DEBUG here ///////////////////////////////////////////


                while(!pQueue.empty())pQueue.pop();
                R rmsdVal=0.1*sqrt(2*rmsd.ExecMetric()/(proteinLength*(proteinLength-1)));
                if(globalBestRMSDFine>rmsdVal)
                {
                    individual temp1(proteinLength-1);
                    encodeM(temp1,c);
                    globalBestFine=temp1;
                    globalBestRMSDFine=rmsdVal;
                    cout << "REMARK "<<globalBestFitness<<" "<<globalBestRMSD<<" "<<globalBestRMSDFine<<" "<<timeUtl.getTime()-timeInit<<endl;
                }
                if(globalBestFitness>bre.ExecMetric())
                {
                    individual temp1(proteinLength-1);
                    encodeM(temp1,c);
                    globalBest=temp1;
                    globalBestFitness=bre.ExecMetric();
                    globalBestRMSD=0.1*sqrt(2*rmsd.ExecMetric()/(proteinLength*(proteinLength-1)));
                    cout << "REMARK "<<globalBestFitness<<" "<<globalBestRMSD<<" "<<globalBestRMSDFine<<" "<<timeUtl.getTime()-timeInit<<endl;
                    c.writeMotifToFile("data.cml");

                }
                else
                {
                    nonImp++;
                }


            } /// end of local search
            encodeM(temp,c);
            population[tIdxOuter]=temp;
            fitness[tIdxOuter]=bre.ExecMetric();


        }/// end of a generation
        ///after each generation crossover()
        ///cout << "REMARK end of generation"<<endl;
        block1<individual,xmm> alterPopulation;
        block1<Int,xmm> alterFitness;

        Idx tIdxAlt=0;
        while(tIdxAlt<50)
        {
            Idx t1=uniform(theRnd,(Z)0,(Z)populationSize-1);
            Idx t2=uniform(theRnd,(Z)0,(Z)populationSize-1);
            if(t1==t2) continue;

                Idx pos=uniform(theRnd,(Z)1,(Z)proteinLength-2);
                Change tempChange;
                /// if pos is a motif then crossover is not allowed
                if(c.isPosMotif(pos)&&c.isPosMotif(pos-1))
                    continue;
                //cout<<t1<< " "<< t2 << " "<<pos <<endl;

                bool cross = crossoverM(population[t1],population[t2],pos,tempChange);
                if(cross==false) continue;
                //cout << "REMARK success in crossover!"<<endl;
                c.execute(tempChange);
                individual tempI(proteinLength-1);
                encodeM(tempI,c);
                alterPopulation.insertMem(temp);
                alterFitness.insertMem(bre.ExecMetric());
                tIdxAlt++;
                R rmsdVal=0.1*sqrt(2*rmsd.ExecMetric()/(proteinLength*(proteinLength-1)));
                if(globalBestRMSDFine>rmsdVal)
                {
                    individual temp1(proteinLength-1);
                    encodeM(temp1,c);
                    globalBestFine=temp1;
                    globalBestRMSDFine=rmsdVal;
                    cout << "REMARK "<<globalBestFitness<<" "<<globalBestRMSD<<" "<<globalBestRMSDFine<<" "<<timeUtl.getTime()-timeInit<<endl;
                }
                if(globalBestFitness>bre.ExecMetric())
                {
                    individual temp1(proteinLength-1);
                    encodeM(temp1,c);
                    globalBest=temp1;
                    globalBestFitness=bre.ExecMetric();
                    globalBestRMSD=0.1*sqrt(2*rmsd.ExecMetric()/(proteinLength*(proteinLength-1)));
                    cout << "REMARK "<<globalBestFitness<<" "<<globalBestRMSD<<" "<<globalBestRMSDFine<<" "<<timeUtl.getTime()-timeInit<<endl;
                    c.writeMotifToFile("data.cml");
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

    }/// end of search

     cout << "REMARK "<<globalBestFitness<<" "<<globalBestRMSD<<" "<<globalBestRMSDFine<<endl;
    individual temp=globalBest;
    Change tChange;
    decodeM(temp,tChange);
    c.execute(tChange);
    c.writePDB("no file",-1);

}

