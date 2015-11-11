

#include "pspl/idx.hh"
using namespace platypus;

int tests(int argC, char* argv[])
{

		Clk tClk = parseN(argv[1]);
        //Int timeOut = parseN(argv[2]);


        Idx timestamp=time(NULL);
		Rnd theRnd(timestamp);
        cout << "Time Stamp: " << timestamp << endl;
		Space::Dimen(3);

		FccLattice fcc;
		Lattice::l(fcc);

		HpEnergy hp;
		Energy::e(hp);
        // PPPHPHHPHHPPPHPHPPPPHPHHPPHPHHHHHPPHHPPHHHHHHPPHPPHHPPHPHPHHHHHPHHPHHHPPPHHHPHHPPHPHPPHPPPHPPHPPHPPHHHPHHHPHPPHPHHPHHHHPHPHHHPHHHPPPPPPHHHHHHPPPPPPPPHHHPPHPHPPPHPHPHPHHPPHHPPPPHHHHHHPPPHHPPPPPHPPPHHPP

        Protein p(2,"PPPHPHHPHHPPPHPHPPPPHPHHPPHPHHHHHPPHHPPHHHHHHPPHPPHHPPHPHPHHHHHPHHPHHHPPPHHHPHHPPHPHPPHPPPHPPHPPHPPHHHPHHHPHPPHPHHPHHHHPHPHHHPHHHPPPPPPHHHHHHPPPPPPPPHHHPPHPHPPPHPHPHPHHPPHHPPPPHHHHHHPPPHHPPPPPHPPPHHPP");
		Protein::p(p);
		kblock<Point> allPoints(Protein::p().Length());
        //ReadFile::readTxtFile(argv[3],allPoints);
        //BlockInit bi(allPoints);
//	Checked
		//DiagonalMove dm;
		TriplePullMove dm;
//		PushMove dm;
//		SinglePullMove dm;
//		CrankShaftMove dm;

//	Unchecked
        //TiltMove dm;
        LNSMove dm1;
//        DoublePullMove dm;

//        TwistMove dm;
		RandomValid rv(theRnd);
		RandomValid rs(theRnd);
		RandomSphere rsp(theRnd);

		Conf c;
//		SelfAvoidCns sa(c,true);
//		ConnectedCns cn(c,true);

        HhNeighCns hn(c,true);
		EnergyObj hpe(c);


        HCoreDistObj hcd(c);
        CentroidDistObj cd(c,false);
        HhDistObj hhd(c);


//		Hdl tTopComb = SumXiHFcMiHn::def(Xm|Hd, c.SysHdl(), tCombParts.items(), tCombParts.size());
//
//		Hdl Heap = TabuMaxHeapHiFrHi::def(c.SysHdl(), Prm(tTopComb, Hd));
//		Hdl VarSelcHdl = RankVarSp::def(c.SysHdl(), Heap);
//		Hdl ValSelcHdl = MinValSdXi::def(Xm, c.SysHdl(), Prm(tTopComb, Mm));


		Hdl Heap = TabuMaxHeapHiFrHi::def(c.SysHdl(), Prm(hn.FuncHdl(), Hd));
		Hdl VarSelcHdl = RankVarSp::def(c.SysHdl(), Heap);


		rsp.compute(c);
		rv.compute(c);
		//c.initialise(rs.FullChange());
        //bi.compute(c);


        c.initialise(rsp.FullChange());
        c.writeHPToFile("data.cml");
        //exit(0);
//		Sys::execAnewRndVars(c.SysHdl(), theRnd);

		Sys::setVarTabuLimit(c.SysHdl(), p.Length());

//		RefHdl(tTop, FcMiHn, c.SysHdl(), hn.FuncHdl());
		Sys const & tSys = Sys::ref(c.SysHdl());



        //cout << tSys.ExecClk() << " " << hn.ExecMetric() << " " << hpe.ExecMetric() << endl;



		while(tSys.ExecClk() < tClk)
		{
			c.writeHPToFile("data.cml");
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
					Sys::simulIncr(c.SysHdl(), hpe.FuncHdl());
					if (tMinIdx == MaxIdx || tMinEng > hpe.SimulMetric())
					{
						tMinEng = hpe.SimulMetric();
						tMinIdx = tIdx;
					}
				}
				if (tMinIdx != MaxIdx) break;
				tSelcPos = uniform(theRnd, Protein::p().Length());
			}

			c.execute(dm.PartChanges()[tMinIdx]);

			cout << tSys.ExecClk() << " " << tSelcPos << " " << hn.ExecMetric() << " " << hpe.ExecMetric() << endl;
		}
		return 0;

}
