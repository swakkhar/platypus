
#include "pspl/idx.hh"

/*!
4PTI.pdb 8.06412 7.92843
4RXN.pdb 6.38357 6.36632
2IGD.pdb 9.37923 9.37923
1R69.pdb 6.72978 6.67083
1ENH.pdb 7.04557 6.71863
1YPA.pdb 8.09074 7.98686
swakkhar@swakkhar-laptop:~/platypus$ sh run.sh
4PTI.pdb 7.16868 7.07672
4RXN.pdb 6.37417 6.29444
2IGD.pdb 9.33649 9.33649
1R69.pdb 6.4846 6.47225
1ENH.pdb 6.65132 6.61211
1YPA.pdb 7.80705 7.53591
swakkhar@swakkhar-laptop:~/platypus$


*/



using namespace platypus;

void pclf(int argC, char* argv[])
{
    //try
    {

        if(argC<3)
        {
            cout << "Not Enough parameters provided!" << endl;
            cout << "Usage:" << argv[0] << " <input_file> <timeout>" << endl;
            return;
        }
        Int timeOut = parseN(argv[2]);

		Rnd theRnd;//(timestamp);
        //cout << "Time Stamp: " << theRnd.Seed() << endl;
		Space::Dimen(3);

		CCLattice fcc;
		Lattice::l(fcc);


        block1<RPoint,xmm> allPoints;
        int start=0;
        int end=0;

		BarerraEnergy br;
		Energy::e(br);
        // PPPHPHHPHHPPPHPHPPPPHPHHPPHPHHHHHPPHHPPHHHHHHPPHPPHHPPHPHPHHHHHPHHPHHHPPPHHHPHHPPHPHPPHPPPHPPHPPHPPHHHPHHHPHPPHPHHPHHHHPHPHHHPHHHPPPPPPHHHHHHPPPPPPPPHHHPPHPHPPPHPHPHPHHPPHHPPPPHHHHHHPPPHHPPPPPHPPPHHPP
        PDBReader::readPDBFile(argv[1],allPoints,start,end);

        //cout << start << " "<<end<<endl;

        if(start<=0)
        {
            cout << "ATOM co-ordinates starts from <=0, cant process it: " << argv[1] <<endl;
            exit(0);
        }
        char pro[1000];
		PDBReader::readSequence(argv[1],pro);

		string protein_seq(pro);

        //cout <<pro<<endl;
        Protein p(20,(protein_seq.substr(start-1,end-start)).c_str());
		Protein::p(p);



		Dim proteinLength=Protein::p().Length();





        std::stringstream strS;
        strS<<argv[1];
        string filename=strS.str();

        //cout << filename << " ";

        //std::ofstream out(filename.c_str());

		/*for (Idx tIdx=0;tIdx<Protein::p().Length();++tIdx)
		{
		    cout << allPoints[tIdx].x
		    << " " << allPoints[tIdx].y
            << " " << allPoints[tIdx].z<<endl;
		}*/

        //ReadFile::readTxtFile(argv[3],allPoints);
        //BlockInit bi(allPoints);


        //cout << "here1" <<endl;

        LNSMove lns;
        SinglePullMove pull;

		Conf c(theRnd);

		EnergyObj bre(c,true);
        //bi.compute(c);

        //cout << "here2" <<endl;


        //cout << "here3" <<endl;
        RMSDObj rmsd(c,true,allPoints);
        //cout << "here4" <<endl;

        ContactTrendObj cto(c,true);
        ACoreDistObj ado(c,true);
        AffineDistObj afo(c,true);
        ContactOrderObj coo(c,true);
        ContactCountObj cco(c,true);
        CentroidDistObj cdo(c,true);

        HpObj hp(c,true);
        MjObj mj(c,true);
        HPNXObj hpnx(c,true);




        RandomSphere rs(theRnd);
        //PCLFInit pc(allPoints);
        //cout<<"here"<<endl;
        //cout << "here5" <<endl;
		rs.compute(c);
		//pc.compute(c);
		//cout << "here5.5" <<endl;
		c.initialise(rs.FullChange());
		//cout << "here6" <<endl;
		//cout<<"also here"<<endl;
        //c.initialise(bi.FullChange());
   		//Sys::refc(c.SysHdl())setVarTabuLimit(c.SysHdl(), p.Length());

		Sys const & tSys = Sys::refc(c.SysHdl());

        Int globalMin= bre.ExecMetric();

        //cout << globalMin << endl;
        //cout << proteinLength << " ";
        //cout<< 0.1*sqrt(globalMin/(proteinLength*(2*proteinLength-1))) << " ";
        //cout << 0.1*sqrt(2*globalMin/(proteinLength*(proteinLength-1))) <<" ";
        //exit(0);
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
        Idx initMaxStable = proteinLength;
        Idx maxStable=initMaxStable;
        Idx initWindowSize = 1;
        Idx windowSize=initWindowSize;
        R mFactor=1.1;


        N lnsSelected=0;
        block1<Change,xmm> partChanges;
        //int fileCounter=0;



        //exit(0);
        //Bln heuRG=true;
        //cout << "here7" <<endl;

        int iteration=0;
        while(timeUtl.getTime()-timeInit <= timeOut)
        {
            //windowSize=uniform(theRnd,4,(Idx)10);
            //c.writeHPToFile("data.cml");
            Int tMinEng; Idx tMinIdx;
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

                    Term::performSimulIncr(Func::ptrm(c.SysHdl(),bre.FuncHdl()));
                    //Sys::simulIncr(c.SysHdl(), bre.FuncHdl());
                    if (tMinIdx == MaxIdx || tMinEng >=bre.SimulMetric())
                    {
                        tMinEng =bre.SimulMetric();
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
            iteration++;
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


            if(bre.ExecMetric() < globalMin)
            {
                globalMin = bre.ExecMetric();
                windowSize = initWindowSize;
                nonImp = 0;
                maxStable = initMaxStable;
                for(Idx tIdx=0;tIdx < Protein::p().Length(); ++tIdx)
                    tabuList[tIdx]=0;

                // change heuristics select
                //std::stringstream strS;
                //strS<<fileCounter;
                //string tt="file"+strS.str()+".pdb";
                //c.writePDB(tt,globalMin/1000.0);
                //fileCounter++;
                //c.writeHPToFile("data.cml");
                //c.writePDB("data.pdb",bre.ExecMetric());
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
                cout <<0.1*sqrt(2*rmsd.ExecMetric()/(proteinLength*(proteinLength-1))) << " " << bre.ExecMetric() << " " << mj.ExecMetric() << " "
                << hp.ExecMetric() << " " << hpnx.ExecMetric() << " " << cdo.ExecMetric() << " "
                << ado.ExecMetric() << " " << afo.ExecMetric() << " " << cco.ExecMetric()<< " " << coo.ExecMetric() << " "
                << cto.ExecMetric() << endl;

                //cout << iteration << " " << 0.1*sqrt(2*rmsd.ExecMetric()/(proteinLength*(proteinLength-1))) << endl;
			    //cout << " window size:" << windowSize << " part Size : " << dm1.PartChanges().size()<<endl;

                /*cout << iteration <<" "<<
                timeUtl.getTime()-timeInit<< " "<<
                0.1*sqrt(2*rmsd.ExecMetric()/(proteinLength*(proteinLength-1))) << " "
                <<rgo.ExecMetric()<<" "
                <<cto.ExecMetric()<<" "<<cdo.ExecMetric()<<" "<< globalMin <<" "
                <<bre.ExecMetric()<<endl;*/
			}

        }
        //out.close();
        //c.writePDB("data.pdb",bre.ExecMetric());
        //cout<< 0.1*sqrt(globalMin/(proteinLength*(2*proteinLength-1))) << " ";
        //cout<< 0.1*sqrt(2*globalMin/(proteinLength*(proteinLength-1))) << endl;

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

