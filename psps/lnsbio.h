using namespace platypus;
#include "pspl/idx.hh"

void lnsmj(int argC, char* argv[])
{

        if(argC<3)
        {
            cout << "Not Enough parameters provided!" << endl;
            cout << "Usage:" << argv[0] << " <protein_seq> <timeout>" << endl;
            return;
        }
        Int timeOut = parseN(argv[2]);

		Rnd theRnd;//(timestamp);
        cout << "Time Stamp: " << theRnd.Seed() << endl;
		Space::Dimen(3);

		FccLattice fcc;
		Lattice::l(fcc);

		MjEnergy br;
		Energy::e(br);
        // PPPHPHHPHHPPPHPHPPPPHPHHPPHPHHHHHPPHHPPHHHHHHPPHPPHHPPHPHPHHHHHPHHPHHHPPPHHHPHHPPHPHPPHPPPHPPHPPHPPHHHPHHHPHPPHPHHPHHHHPHPHHHPHHHPPPPPPHHHHHHPPPPPPPPHHHPPHPHPPPHPHPHPHHPPHHPPPPHHHHHHPPPHHPPPPPHPPPHHPP

        Protein p(20,argv[1]);
		Protein::p(p);
		block1<Point,kmm> allPoints(Protein::p().Length());
        //ReadFile::readTxtFile(argv[3],allPoints);
        //BlockInit bi(allPoints);
        LNSMove lns;
        DoublePullMove pull;

		Conf c(theRnd);

		EnergyObj bre(c,true);
        //bi.compute(c);



        ContactOrderObj coo(c,true);
        ContactCountObj cco(c,true);
        MinusDistObj mdo(c,true);
        PlusDistObj pdo(c,true);
        CentroidDistObj cdo(c,true);

        ACoreDistObj aco(c,true);
        AffineDistObj ado(c,true);

        RandomSphere rs(theRnd);
        //cout<<"here"<<endl;
		rs.compute(c);
		c.initialise(rs.FullChange());
		//cout<<"also here"<<endl;
        //c.initialise(bi.FullChange());
   		//Sys::refc(c.SysHdl())setVarTabuLimit(c.SysHdl(), p.Length());

		Sys const & tSys = Sys::refc(c.SysHdl());

        Int globalMin= bre.ExecMetric();
        cout << globalMin << endl;
        //exit(0);
        // variables for tabu
        block1<Idx,kmm> tabuList(Protein::p().Length());
        Idx tabuLength = uniform(theRnd,4,(Z)Protein::p().Span()/8);
        // initialize tabu
        for(Idx tIdx=0;tIdx < Protein::p().Length(); ++tIdx)
            tabuList[tIdx]=0;

        block1<platypus::Pos,xmm> curPositions;
        TimeUtl timeUtl;
        R timeInit= timeUtl.getTime();
        Idx nonImp =0;
        Idx initMaxStable = 1000;
        Idx maxStable=initMaxStable;
        Idx initWindowSize = 1;
        Idx windowSize=initWindowSize;
        R mFactor=1.2;


        N lnsSelected=0;
        block1<Change,xmm> partChanges;
        int fileCounter=0;

        typedef tuple2a<Int,Idx> mytuple;
        priority_queue<mytuple,vector<mytuple>,greater_equal<mytuple> > pQueue;

        while(timeUtl.getTime()-timeInit <= timeOut)
        {
            //windowSize=uniform(theRnd,4,(Idx)10);
            //c.writeHPToFile("data.cml");
            //Int tMinEng; Idx tMinIdx;
            partChanges.clear();
            platypus::Pos tSelcPos=0;

            lnsSelected=uniform(theRnd,0,3);

            if(lnsSelected<=2)
            {

                curPositions.clear();

                lns.reset();

                // add positions to the move


                Idx windowType = uniform(theRnd,1,2);

                if(windowType==1)
                {


                // consecutive window

                    tSelcPos = uniform(theRnd, 1+windowSize/2 , Protein::p().Span()-1-windowSize/2);
                    B exitFor = false;
                    for(platypus::Pos tPos = 0; tPos < windowSize ; ++tPos)
                    {
                        if(tabuList[tSelcPos-windowSize/2+tPos] < tSys.ExecClk())
                        {
                            lns.addPosition(tSelcPos-windowSize/2+tPos);
                            curPositions.insertMem(tSelcPos-windowSize/2+tPos);
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
                        if(!lns.alreadyAdded(tSelcPos) && (tabuList[tSelcPos] < tSys.ExecClk()))
                        {
                            lns.addPosition(tSelcPos);
                            curPositions.insertMem(tSelcPos);
                            //cout << tSelcPos << " " ;
                        }
                        else --tPos;

                    }
                }


                lns.compute(c, tSelcPos);
                //cout << " window size:" << windowSize << " part Size : " << dm1.PartChanges().size()<<endl;
                if(lns.PartChanges().itemCount()==0) continue;
                partChanges=lns.PartChanges();
            }
            else
            {
                // random pull move
                while(1)
                {
                    tSelcPos = uniform(theRnd, 1 , (Z)Protein::p().Span()-1);
                    if(tabuList[tSelcPos] < tSys.ExecClk())
                    {
                        break;
                    }
                }
                pull.compute(c,tSelcPos);
                if(pull.PartChanges().itemCount()==0) continue;
                partChanges = pull.PartChanges();
            }

            /********** heuristics************/
            //tMinEng = MaxInt;//, tSolCount = 0;
			//tMinIdx = MaxIdx;
            Idx heuristics = uniform(theRnd,0,9);
			for(Idx tIdx = 0; tIdx < partChanges.itemCount(); ++tIdx)
			{
			//cout << tIdx << " ";
                c.simulate(partChanges[tIdx]);

                if(heuristics<4)
                {
                    Term::performSimulIncr(Func::ptrm(c.SysHdl(), bre.FuncHdl()));
                    mytuple t(bre.SimulMetric(),tIdx);
                    pQueue.push(t);
                    /*if (tMinIdx == MaxIdx || tMinEng >= bre.SimulMetric())
                    {
                        tMinEng =bre.SimulMetric();
                        tMinIdx = tIdx;
                    }*/
                }
                else if(heuristics<7)
                {
                    Term::performSimulIncr(Func::ptrm(c.SysHdl(), mdo.FuncHdl()));
                    mytuple t(mdo.SimulMetric(),tIdx);
                    pQueue.push(t);
                    /*if (tMinIdx == MaxIdx || tMinEng <= coo.SimulMetric())
                    {
                        tMinEng =coo.SimulMetric();
                        tMinIdx = tIdx;
                    }*/
                }
                else if(heuristics<8)
                {
                    Term::performSimulIncr(Func::ptrm(c.SysHdl(), coo.FuncHdl()));
                    mytuple t(-1*coo.SimulMetric(),tIdx);
                    pQueue.push(t);
                    /*if (tMinIdx == MaxIdx || tMinEng <= cco.SimulMetric())
                    {
                        tMinEng =mdo.SimulMetric();
                        tMinIdx = tIdx;
                    }*/
                }
                else if(heuristics<9)
                {
                    Term::performSimulIncr(Func::ptrm(c.SysHdl(), cco.FuncHdl()));
                    mytuple t(-1*cco.SimulMetric(),tIdx);
                    pQueue.push(t);
                    /*if (tMinIdx == MaxIdx || tMinEng >= cco.SimulMetric())
                    {
                        tMinEng =cco.SimulMetric();
                        tMinIdx = tIdx;
                    }*/
                }
                else if(heuristics == 9)
                {
                    Term::performSimulIncr(Func::ptrm(c.SysHdl(), pdo.FuncHdl()));
                    mytuple t(-1*pdo.SimulMetric(),tIdx);
                    pQueue.push(t);
                    /*if (tMinIdx == MaxIdx || tMinEng <= coo.SimulMetric())
                    {
                        tMinEng =coo.SimulMetric();
                        tMinIdx = tIdx;
                    }*/
                }



            }
            //if (tMinIdx == MaxIdx) continue;
            //cout << tMinIdx <<endl;

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
                if(lnsSelected<=2)
                {
                    for(Idx tIdx =0; tIdx < curPositions.itemCount();++tIdx)
                    {
                        tabuList[tIdx]=tSys.ExecClk()+tabuLength;
                    }
                }
                else
                {
                    tabuList[tSelcPos]=tSys.ExecClk()+tabuLength;
                }
                nonImp++;
                continue;
            }
            //cout << "selected id:"<<selectedId  << "partial changes" << partChanges.itemCount()<< endl;
            //cout << partialChanges[selectedId].size()<< endl;
            c.execute(partChanges[selectedId]);
            while(!pQueue.empty())pQueue.pop();

            if(bre.ExecMetric() < globalMin)
            {
                globalMin = bre.ExecMetric();
                windowSize = initWindowSize;
                nonImp = 0;
                maxStable = initMaxStable;
                for(Idx tIdx=0;tIdx < Protein::p().Length(); ++tIdx)
                    tabuList[tIdx]=0;


                std::stringstream strS;
                strS<<fileCounter;
                string tt="file"+strS.str()+".pdb";
                //c.writePDB(tt,globalMin/1000.0);
                fileCounter++;
            }
            else
            {
                nonImp++;
            }

            if(nonImp>=maxStable)
            {
                nonImp=0;
                windowSize++;
                maxStable=maxStable*mFactor;
                for(Idx tIdx =0; tIdx < curPositions.itemCount();++tIdx)
               {
                    tabuList[tIdx]=0;
               }

			//if(tSys.ExecClk()%1000==0)
			{
			    //cout << " window size:" << windowSize << " part Size : " << dm1.PartChanges().size()<<endl;

                cout << tSys.ExecClk() <<" "<<
                timeUtl.getTime()-timeInit<< " "
                << cco.ExecMetric() << " "
                << coo.ExecMetric() << " " <<
                bre.ExecMetric() << " "<< globalMin <<endl;
			}
			}
            // Update tabu list
            if(lnsSelected<=2)
            {

               for(Idx tIdx =0; tIdx < curPositions.itemCount();++tIdx)
               {
                    tabuList[tIdx]=tSys.ExecClk()+tabuLength;
               }
            }
            else
            {
                tabuList[tSelcPos]=tSys.ExecClk()+tabuLength;
            }
        }
        cout<< globalMin << endl;

//        for(Idx tIdx=0;tIdx<Protein::p().Length();++tIdx)
//        {
//            cout << c.Coord(tIdx).Comp(0) << " " << c.Coord(tIdx).Comp(1) << " " << c.Coord(tIdx).Comp(2) << endl;
//        }
//        cout << "End of Run" << endl;

}

void lnsbio(int argC, char* argv[])
{


        Int timeOut = parseN(argv[2]);


        Idx timestamp=time(NULL);
		Rnd theRnd(timestamp);
        cout << "Time Stamp: " << timestamp << endl;
		Space::Dimen(3);

		FccLattice fcc;
		Lattice::l(fcc);

		HpEnergy hp;
		Energy::e(hp);
        // PPPHPHHPHHPPPHPHPPPPHPHHPPHPHHHHHPPHHPPHHHHHHPPHPPHHPPHPHPHHHHHPHHPHHHPPPHHHPHHPPHPHPPHPPPHPPHPPHPPHHHPHHHPHPPHPHHPHHHHPHPHHHPHHHPPPPPPHHHHHHPPPPPPPPHHHPPHPHPPPHPHPHPHHPPHHPPPPHHHHHHPPPHHPPPPPHPPPHHPP

        Protein p(2,argv[1]);
		Protein::p(p);
		block1<Point,kmm> allPoints(Protein::p().Length());
        //ReadFile::readTxtFile(argv[3],allPoints);
        //BlockInit bi(allPoints);
        //RandomSphere rs(theRnd);

        RandomStructured rs(theRnd);
        //ChainGrowth rs(theRnd);
        LNSMove dm1;

		Conf c(theRnd);

        //HhNeighCns hn(c,true);
		EnergyObj hpe(c,true);


        HCoreDistObj hcd(c,false);
        CentroidDistObj cd(c,false);
        HhDistObj hhd(c,false);


        //bi.compute(c);
        rs.compute(c);

        //c.initialise(bi.FullChange());
        c.initialise(rs.FullChange());

        //c.writePDB(argv[4]);

        c.writeHPToFile("data.cml");

		//Sys::setVarTabuLimit(c.SysHdl(), p.Length());

		Sys const & tSys = Sys::refc(c.SysHdl());

        Int globalMin= hpe.ExecMetric();
        cout << globalMin << endl;
        exit(0);
        // variables for tabu
        block1<Idx,kmm> tabuList(Protein::p().Length());
        Idx tabuLength = uniform(theRnd,4,(Z)Protein::p().Span()/8);
        // initialize tabu
        for(Idx tIdx=0;tIdx < Protein::p().Length(); ++tIdx)
            tabuList[tIdx]=0;

        block1<platypus::Pos,xmm> curPositions;
        TimeUtl timeUtl;
        R timeInit= timeUtl.getTime();
        Idx nonImp =0;
        Idx initMaxStable = 1000;
        Idx maxStable=initMaxStable;
        Idx initWindowSize = 1;
        Idx windowSize=initWindowSize;
        R mFactor=1.2;

        //int fileCounter=0;

        typedef tuple2a<Int,Idx> mytuple;
        priority_queue<mytuple,vector<mytuple>,greater_equal<mytuple> > pQueue;
        while(timeUtl.getTime()-timeInit <= timeOut)
        {
            //windowSize=uniform(theRnd,4,(Idx)10);
            //c.writeHPToFile("data.cml");
            //Int tMinEng; Idx tMinIdx;

            curPositions.clear();

            dm1.reset();

            // add positions to the move
            platypus::Pos tSelcPos=0;

            Idx windowType = uniform(theRnd,1,2);

            if(windowType==1)
            {


            // consecutive window

                tSelcPos = uniform(theRnd, 1+windowSize/2 , Protein::p().Span()-1-windowSize/2);
                B exitFor = false;
                for(platypus::Pos tPos = 0; tPos < windowSize ; ++tPos)
                {
                    if(tabuList[tSelcPos-windowSize/2+tPos] < tSys.ExecClk())
                    {
                        dm1.addPosition(tSelcPos-windowSize/2+tPos);
                        curPositions.insertMem(tSelcPos-windowSize/2+tPos);
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
                    if(!dm1.alreadyAdded(tSelcPos) && (tabuList[tSelcPos] < tSys.ExecClk()))
                    {
                        dm1.addPosition(tSelcPos);
                        curPositions.insertMem(tSelcPos);
                        //cout << tSelcPos << " " ;
                    }
                    else --tPos;

                }
            }


            dm1.compute(c, tSelcPos);
            //cout << " window size:" << windowSize << " part Size : " << dm1.PartChanges().size()<<endl;
            if(dm1.PartChanges().itemCount()==0) continue;

            //tMinEng = MaxInt;//, tSolCount = 0;
			//tMinIdx = MaxIdx;
            Idx heuristics = uniform(theRnd,1,3);
			for(Idx tIdx = 0; tIdx < dm1.PartChanges().itemCount(); ++tIdx)
			{
                c.simulate(dm1.PartChanges()[tIdx]);

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
                    mytuple t(hcd.SimulMetric(),tIdx);
                    pQueue.push(t);
                    /*if (tMinIdx == MaxIdx || tMinEng >= hhd.SimulMetric())
                    {
                        tMinEng =hhd.SimulMetric();
                        tMinIdx = tIdx;
                    }*/
                }
                else
                {
                    Term::performSimulIncr(Func::ptrm(c.SysHdl(), hpe.FuncHdl()));
                    mytuple t(hcd.SimulMetric(),tIdx);
                    pQueue.push(t);
                    /*if (tMinIdx == MaxIdx || tMinEng >= hpe.SimulMetric())
                    {
                        tMinEng =hpe.SimulMetric();
                        tMinIdx = tIdx;
                    }*/
                }



            }

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
            continue;
        }
            //if (tMinIdx == MaxIdx) continue;
            c.execute(dm1.PartChanges()[selectedId]);
            if(hpe.ExecMetric() < globalMin)
            {
                globalMin = hpe.ExecMetric();
                cout << "global min: "<<globalMin<<endl;
                windowSize = initWindowSize;
                nonImp = 0;
                maxStable = initMaxStable;
                for(Idx tIdx=0;tIdx < Protein::p().Length(); ++tIdx)
                    tabuList[tIdx]=0;
                //std::stringstream strS;
                //strS<<fileCounter;
                //string tt="file"+strS.str()+".pdb";
                //c.writePDB(tt,globalMin/1000.0);
                //sfileCounter++;
            }
            else
            {
                nonImp++;
            }

            if(nonImp>=maxStable)
            {
                //cout << "Stagnation!" << endl;
                nonImp=0;
                windowSize++;
                maxStable=maxStable*mFactor;

			//if(tSys.ExecClk()%1000==0)
			{
			    //cout << " window size:" << windowSize << " part Size : " << dm1.PartChanges().size()<<endl;

                cout << tSys.ExecClk() <<" "<< timeUtl.getTime()-timeInit<< " " <<hhd.ExecMetric() << " " << hpe.ExecMetric() << " "<< globalMin <<endl;
			}
			}
            // Update tabu list
            for(Idx tIdx =0; tIdx < curPositions.itemCount();++tIdx)
            {
                tabuList[tIdx]=tSys.ExecClk()+tabuLength;
            }

        }
        cout<< globalMin << endl;



        /*while(tSys.ExecClk() < tClk)
        {
            c.writeHPToFile("data.cml");
            //Sys::runSelc(c.SysHdl(), VarSelcHdl, theRnd);
			//platypus::Pos tSelcPos = Selc::ref(c.SysHdl(), VarSelcHdl).SelcVars()[0] / Space::Dimen();
            platypus::Pos tSelcPos = uniform(theRnd, 1+windowSize/2, Protein::p().Span()-windowSize/2-1);
            Int tMinEng; Idx tMinIdx;



            while(true)
            {
                cout<< "selected!" << tSelcPos <<endl;
                dm1.reset();
                // add positions to the move
                for(platypus::Pos tPos = 0; tPos < windowSize ; ++tPos)
                {
                    dm1.addPosition(tSelcPos-windowSize/2 + tPos);
                }
                dm1.compute(c, tSelcPos);
                cout << "part Size : " << dm1.PartChanges().size()<<endl;
                //exit(0);
                tMinEng = MaxInt;//, tSolCount = 0;
				tMinIdx = MaxIdx;
				for(Idx tIdx = 0; tIdx < dm1.PartChanges().size(); ++tIdx)
				{
					c.simulate(dm1.PartChanges()[tIdx]);
					Sys::simulIncr(c.SysHdl(), hpe.FuncHdl());
					if (tMinIdx == MaxIdx || tMinEng >= hpe.SimulMetric())
					{
						tMinEng = hpe.SimulMetric();
						tMinIdx = tIdx;
					}
				}
				if (tMinIdx != MaxIdx) break;
				tSelcPos = uniform(theRnd, 1+windowSize/2, Protein::p().Span()-windowSize/2-1);
				//tSelcPos = uniform(theRnd, Protein::p().Length());
				//while(tSelcPos-windowSize/2 > 0||tSelcPos+windowSize/2 >= Protein::p().Span())
                    //tSelcPos = uniform(theRnd, Protein::p().Length());
            }
            c.execute(dm1.PartChanges()[tMinIdx]);
            if(hpe.ExecMetric() < globalMin) globalMin = hpe.ExecMetric();
			cout << tSys.ExecClk() << " " << hn.ExecMetric() << " " << hpe.ExecMetric() << " "<< globalMin <<endl;

        }*/
}
