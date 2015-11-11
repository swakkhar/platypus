using namespace platypus;
#include "pspl/idx.hh"
void lnstrend(int argC, char* argv[])
{
    //try
    {

        if(argC<4)
        {
            cout << "Not Enough parameters provided!" << endl;
            cout << "Usage:" << argv[0] << " <protein_seq> <timeout> <input_file>" << endl;
            return;
        }
        Int timeOut = parseN(argv[2]);

		Rnd theRnd;//(timestamp);
        //cout << "Time Stamp: " << theRnd.Seed() << endl;
		Space::Dimen(3);

		FccLattice fcc;
		Lattice::l(fcc);



		BarerraEnergy br;
		Energy::e(br);
        // PPPHPHHPHHPPPHPHPPPPHPHHPPHPHHHHHPPHHPPHHHHHHPPHPPHHPPHPHPHHHHHPHHPHHHPPPHHHPHHPPHPHPPHPPPHPPHPPHPPHHHPHHHPHPPHPHHPHHHHPHPHHHPHHHPPPPPPHHHHHHPPPPPPPPHHHPPHPHPPPHPHPHPHHPPHHPPPPHHHHHHPPPHHPPPPPHPPPHHPP

        Protein p(20,argv[1]);
		Protein::p(p);

		Dim proteinLength=Protein::p().Length();

		block1<RPoint,xmm> allPoints(Protein::p().Length());
		int start=0, end=0;
		PDBReader::readPDBFile(argv[3],allPoints,start,end);

        block1<Point,kmm> Points(Protein::p().Length());
		for (Idx tIdx=0;tIdx<Protein::p().Length();++tIdx)
		{
		    //cout << allPoints[tIdx].x
		    //<< " " << allPoints[tIdx].y
            //<< " " << allPoints[tIdx].z<<endl;
            Points[tIdx].Comp(0)=allPoints[tIdx].x/2.69;
            Points[tIdx].Comp(1)=allPoints[tIdx].y/2.69;
            Points[tIdx].Comp(2)=allPoints[tIdx].z/2.69;

            //cout << Points[tIdx].Comp(0) << " "<<
		//Points[tIdx].Comp(1) << " "<< Points[tIdx].Comp(2)<<endl;
		}
        //cout << proteinLength<<endl;
		//exit(0);
        //ReadFile::readTxtFile(argv[3],allPoints);
        //BlockInit bi(allPoints);




        LNSMove lns;
        DoublePullMove pull;

		Conf c;

		EnergyObj bre(c,true);
        //bi.compute(c);


        ContactTrendObj cto(c,true);
        RMSDObj rmsd(c,true,allPoints);


        //ACoreDistObj ado(c,true);

        //ContactOrderObj coo(c,true);
        //ContactCountObj cco(c,true);
        //MinusDistObj mdo(c,true);
        //PlusDistObj pdo(c,true);

        //CentroidDistObj cdo(c,true);
        //OriginDistObj odo(c,true);

        RGyrationObj rgo(c,false);


        //RandomSphere rs(theRnd);
        //cout<<"here"<<endl;

		BlockInit rs(Points);
		rs.compute(c);
		c.initialise(rs.FullChange());
		//cout<<"also here"<<endl;
        //c.initialise(bi.FullChange());
   		//Sys::refc(c.SysHdl())setVarTabuLimit(c.SysHdl(), p.Length());

		Sys const & tSys = Sys::refc(c.SysHdl());

        Int globalMin= cto.ExecMetric();
        cout << globalMin << endl;
        exit(0);
        // variables for tabu
        block1<Idx,kmm> tabuList(Protein::p().Length());
        Idx tabuLength = uniform(theRnd,4,(int)Protein::p().Span()/8);
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



        //exit(0);
        //Bln heuRG=true;
        while(timeUtl.getTime()-timeInit <= timeOut)
        {
            //windowSize=uniform(theRnd,4,(Idx)10);
            //c.writeHPToFile("data.cml");
            Int tMinEng; Idx tMinIdx;
            platypus::Pos tSelcPos=0;

            lnsSelected=uniform(theRnd,0,2);

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
                        tSelcPos = uniform(theRnd, 1 , (int)Protein::p().Span()-1);
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
                    tSelcPos = uniform(theRnd, 1 , (int)Protein::p().Span()-1);
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
            tMinEng = MaxInt;//, tSolCount = 0;
			tMinIdx = MaxIdx;
            Idx heuristics = uniform(theRnd,0,5);
			for(Idx tIdx = 0; tIdx < partChanges.itemCount(); ++tIdx)
			{
                c.simulate(partChanges[tIdx]);


                /*if(heuRG==true)
                {
                    Term::performSimulIncr(Func::ptrm(c.SysHdl(), rgo.FuncHdl()));
                    //Sys::simulIncr(c.SysHdl(), bre.FuncHdl());
                    if (tMinIdx == MaxIdx || tMinEng >= rgo.SimulMetric())
                    {
                        tMinEng =rgo.SimulMetric();
                        tMinIdx = tIdx;
                    }
                }
                else
                {
                    Term::performSimulIncr(Func::ptrm(c.SysHdl(), cto.FuncHdl()));
                    //Sys::simulIncr(c.SysHdl(), bre.FuncHdl());
                    if (tMinIdx == MaxIdx || tMinEng >= cto.SimulMetric())
                    {
                        tMinEng =cto.SimulMetric();
                        tMinIdx = tIdx;
                    }
                }*/

                if(heuristics<6)
                {

                    Term::performSimulIncr(Func::ptrm(c.SysHdl(), rmsd.FuncHdl()));
                    //Sys::simulIncr(c.SysHdl(), bre.FuncHdl());
                    if (tMinIdx == MaxIdx || tMinEng >= rmsd.SimulMetric())
                    {
                        tMinEng =rmsd.SimulMetric();
                        tMinIdx = tIdx;
                    }
                }
                /*else if(heuristics<6)
                {

                    Term::performSimulIncr(Func::ptrm(c.SysHdl(), rgo.FuncHdl()));
                    //Sys::simulIncr(c.SysHdl(), bre.FuncHdl());
                    if (tMinIdx == MaxIdx || tMinEng >= rgo.SimulMetric())
                    {
                        tMinEng =rgo.SimulMetric();
                        tMinIdx = tIdx;
                    }
                }
                else if(heuristics<0)
                {

                    Term::performSimulIncr(Func::ptrm(c.SysHdl(), rmsd.FuncHdl()));
                    //Sys::simulIncr(c.SysHdl(), bre.FuncHdl());
                    if (tMinIdx == MaxIdx || tMinEng >= rmsd.SimulMetric())
                    {
                        tMinEng =rmsd.SimulMetric();
                        tMinIdx = tIdx;
                    }
                }*/

            }
            if (tMinIdx == MaxIdx) continue;

            c.execute(partChanges[tMinIdx]);
            /*if(heuRG==true && (rgo.ExecMetric())<50)
                {
                    heuRG=false;
                    cout << "toggle to CTO"<<endl;
                }
                else if(heuRG==false && abs(rgo.ExecMetric())>200)
                {
                    heuRG=true;
                    cout << "toggle to RGO"<<endl;
                }*/


            if(cto.ExecMetric() < globalMin)
            {
                globalMin = cto.ExecMetric();
                windowSize = initWindowSize;
                nonImp = 0;
                maxStable = initMaxStable;
                for(Idx tIdx=0;tIdx < Protein::p().Length(); ++tIdx)
                    tabuList[tIdx]=0;

                // change heuristics select
                std::stringstream strS;
                strS<<fileCounter;
                string tt="file"+strS.str()+".pdb";
                //c.writePDB(tt,globalMin/1000.0);
                fileCounter++;
                c.writeHPToFile("data.cml");
                c.writePDB("data.pdb",bre.ExecMetric());
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
            //if(tSys.ExecClk()%100==0)
            {
			    //cout << " window size:" << windowSize << " part Size : " << dm1.PartChanges().size()<<endl;

                cout << tSys.ExecClk() <<" "<<
                timeUtl.getTime()-timeInit<< " "<<
                0.1*sqrt(2*rmsd.ExecMetric()/(proteinLength*(proteinLength-1))) << " "
                <<"rgo:"<<rgo.ExecMetric()<<" "
                //<<cto.ExecMetric()<<" "<<cdo.ExecMetric()<<" "<< globalMin <<" "
                <<bre.ExecMetric()<<endl;
			}

        }
        c.writePDB("data.pdb",bre.ExecMetric());
        cout<< globalMin << endl;

//        for(Idx tIdx=0;tIdx<Protein::p().Length();++tIdx)
//        {
//            cout << c.Coord(tIdx).Comp(0) << " " << c.Coord(tIdx).Comp(1) << " " << c.Coord(tIdx).Comp(2) << endl;
//        }
//        cout << "End of Run" << endl;
    }
    //catch(Err &e)
    {
      //  cout << e()<<endl;
    }
}
