#include "pspl/idx.hh"
#include <queue>
#include <vector>

using namespace platypus;


int mj(int argC, char* argv[])
{

        Clk tClk = parseN(argv[1]);
        //Int timeOut = parseN(argv[2]);


        Idx timestamp=time(NULL);//1346289815
		Rnd theRnd(timestamp);
        cout << "Time Stamp: " << timestamp << endl;
		Space::Dimen(3);

		CCLattice fcc;
		Lattice::l(fcc);

		//HpEnergy hp;
		//Energy::e(hp);
        MjEnergy mj;
        Energy::e(mj);


        Protein p(20,"FRTRPLNHDFYNYKIWEPFKPADFPKAWDRMLDHVWDSMASWGHQHCS");
        //Protein p(20,"MKKYTCTVCGYIYDPEDGDPDDGVNPGTDFKDIPDDWVCPLCGVGKDEFEEVEE");


		Protein::p(p);
		//kblock<Point> allPoints(Protein::p().Length());
        //ReadFile::readTxtFile(argv[3],allPoints);
        //BlockInit bi(allPoints);
//	Checked
		//DiagonalMove dm;
		DoublePullMove dm;



		Conf c;

        EnergyObj hpe(c);


        MjMinusDistObj cd1(c,true);
        MjPlusDistObj cd2(c,false);

        ContactCountObj ccc(c,true);
        ContactOrderObj cco(c,true);

		//Hdl Heap = TabuMinHeapHiFrHi::def(c.SysHdl(), Prm(cd.FuncHdl(), Hd));
		//Hdl VarSelcHdl = RankVarSp::def(c.SysHdl(), Heap);

        RandomSphere rs(theRnd);
		rs.compute(c);
		c.initialise(rs.FullChange());
        c.writeHPToFile("data.cml");
        //exit(0);
//		Sys::execAnewRndVars(c.SysHdl(), theRnd);

		Sys::setVarTabuLimit(c.SysHdl(), p.Length());

//		RefHdl(tTop, FcMiHn, c.SysHdl(), hn.FuncHdl());
		Sys const & tSys = Sys::ref(c.SysHdl());



        //cout << tSys.ExecClk() << " " << hn.ExecMetric() << " " << hpe.ExecMetric() << endl;
        // variables for tabu
        kblock<Idx> tabuList(Protein::p().Length());
        Idx tabuLength = uniform(theRnd,4,Protein::p().Span()/4);
        // initialize tabu
        for(Idx tIdx=0;tIdx < Protein::p().Length(); ++tIdx)
            tabuList[tIdx]=0;

        Int globalMin=hpe.ExecMetric();
		cout << "start---!" << endl;
        Change bestChange;
        N nonImp =0;
        N maxStable =20;
		std::priority_queue<Change,std::vector<Change>,std::greater<Change> > changeSet;
		while(tSys.ExecClk() < tClk)
		{
		    //c.writeHPToFile("data.cml");
		    //c.writePDB("data.pdb");
			//Sys::runSelc(c.SysHdl(), VarSelcHdl, theRnd);
			//platypus::Pos tSelcPos = Selc::ref(c.SysHdl(), VarSelcHdl).SelcVars()[0] / Space::Dimen();

            for(platypus::Pos tPos=0;tPos < Protein::p().Length();++tPos)
            {

                if(tabuList[tPos] > tSys.ExecClk() /*&& hpe.SimulMetric() >= hpe.ExecMetric()*/) continue;

                dm.compute(c,tPos);
                for(Idx tIdx = 0; tIdx < dm.PartChanges().size(); ++tIdx)
                {
                    c.simulate(dm.PartChanges()[tIdx]);
					Sys::simulIncr(c.SysHdl(), cd1.FuncHdl());
                    Change tempChange=dm.PartChanges()[tIdx];
                    tempChange.addFVal(cd1.SimulMetric());

                    //Sys::simulIncr(c.SysHdl(), cd2.FuncHdl());
                    //tempChange.addFVal(-cd2.SimulMetric());



                    //cout << hpe.SimulMetric() << " ";
                    tempChange.addFVal(tPos);
					changeSet.push(tempChange);
                }
            }
            //cout << endl;
            if (changeSet.size()==0)
            {
                cout << "No moves possible! "<<tSys.ExecClk()<< endl;
                break;
            }
            // else tiebreak from top elements

            std::vector<Change> topList;
            Change tChange = changeSet.top();
            //cout <<"Selected Metric:" <<tChange.getFVal(0) <<endl;
            topList.push_back(tChange);
            changeSet.pop();
            while(changeSet.size() > 0 && tChange == changeSet.top())
            {
                tChange = changeSet.top();
                topList.push_back(tChange);
                changeSet.pop();
            }
            tChange = topList.at(uniform(theRnd,0,(N)topList.size()-1));
            //cout << tChange.getFVal(0) << endl;
            c.execute(tChange);
            while(!changeSet.empty())changeSet.pop();

            tabuList[tChange.getFVal(tChange.sizeFVal()-1)]=tSys.ExecClk()+tabuLength;
            if(globalMin>hpe.ExecMetric())
            {
                globalMin=hpe.ExecMetric();
                bestChange.reset();
                for(Idx tIdx =0; tIdx < Protein::p().Length();++tIdx)
                {
                    bestChange.add(tIdx,c.Coord(tIdx));
                }
            }
            else
            {
                nonImp++;
            }

            if(tSys.ExecClk()%10==0)
			cout << tSys.ExecClk() << " pos:"
			<< " order: "<<cco.ExecMetric() << " "
			<< " Contacts:" << ccc.ExecMetric() << " "
			<< " co:" << cco.ExecMetric() / (R)ccc.ExecMetric() << " "
			<<tChange.getFVal(tChange.sizeFVal()-1)
			<< " minDis:" << cd1.ExecMetric()
			<<" " << hpe.ExecMetric()<<" "<<globalMin/1000.0 << endl;
            if(nonImp>=maxStable)
            {
                c.execute(bestChange);
                nonImp=0;
            }
		}
		return 0;

}
int mjVarSelector(int argC, char* argv[])
{

        Clk tClk = parseN(argv[1]);
        //Int timeOut = parseN(argv[2]);


        Idx timestamp=time(NULL);//1346289815
		Rnd theRnd(timestamp);
        cout << "Time Stamp: " << timestamp << endl;
		Space::Dimen(3);

		CCLattice fcc;
		Lattice::l(fcc);

		//HpEnergy hp;
		//Energy::e(hp);
        MjEnergy mj;
        Energy::e(mj);


        Protein p(20,"FRTRPLNHDFYNYKIWEPFKPADFPKAWDRMLDHVWDSMASWGHQHCS");
		Protein::p(p);
		//kblock<Point> allPoints(Protein::p().Length());
        //ReadFile::readTxtFile(argv[3],allPoints);
        //BlockInit bi(allPoints);
//	Checked
		DoublePullMove dm;

		RandomValid rs(theRnd);

		Conf c;

        EnergyObj hpe(c);


        MjMinusDistObj cd(c,true);



		Hdl Heap = TabuMinHeapHiFrHi::def(c.SysHdl(), Prm(cd.FuncHdl(), Hd));
		Hdl VarSelcHdl = RankVarSp::def(c.SysHdl(), Heap);


		rs.compute(c);
		c.initialise(rs.FullChange());
        c.writeHPToFile("data.cml");
        //exit(0);
//		Sys::execAnewRndVars(c.SysHdl(), theRnd);

		Sys::setVarTabuLimit(c.SysHdl(), p.Length());

//		RefHdl(tTop, FcMiHn, c.SysHdl(), hn.FuncHdl());
		Sys const & tSys = Sys::ref(c.SysHdl());



        //cout << tSys.ExecClk() << " " << hn.ExecMetric() << " " << hpe.ExecMetric() << endl;


        Int globalMin=hpe.ExecMetric();
		cout << "start---" << endl;
		while(tSys.ExecClk() < tClk)
		{
		    //c.writeHPToFile("data.cml");
		    //c.writePDB("data.pdb");
			Sys::runSelc(c.SysHdl(), VarSelcHdl, theRnd);
			platypus::Pos tSelcPos = Selc::ref(c.SysHdl(), VarSelcHdl).SelcVars()[0] / Space::Dimen();

			Int tMinEng; Idx tMinIdx;
			while(true)
			{
				dm.compute(c, tSelcPos);

				tMinEng = MaxInt;//, tSolCount = 0;
				tMinIdx = MaxIdx;
				for(Idx tIdx = 0; tIdx < dm.PartChanges().size(); ++tIdx)
				{
					c.simulate(dm.PartChanges()[tIdx]);
					Sys::simulIncr(c.SysHdl(), cd.FuncHdl());
					if (tMinIdx == MaxIdx || tMinEng > cd.SimulMetric())
					{
						tMinEng = cd.SimulMetric();
						tMinIdx = tIdx;
					}
				}
				if (tMinIdx != MaxIdx) break;
				tSelcPos = uniform(theRnd, Protein::p().Length());
			}

			c.execute(dm.PartChanges()[tMinIdx]);
            if(globalMin>hpe.ExecMetric())
                globalMin=hpe.ExecMetric();
			cout << tSys.ExecClk() << " " << tSelcPos << " "
			<< cd.ExecMetric() << " " << hpe.ExecMetric()<<" "<<globalMin/1000.0 << endl;
		}
		return 0;

}
