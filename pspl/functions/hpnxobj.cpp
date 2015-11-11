/*!
	@file hpmodel/HPNXObj.cpp
	@brief The implementation file for HPNXObj class.
	@details This is the implementation file for HPNXObj class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#include "pspl/functions/hpnxobj.hpp"



openPlatypusSpace

const int toHPNX[]={3,0,0,0,0,0,0,0,0,0,3,3,3,3,2,2,1,1,1,0};
const int HPmatrix[4][4]=
{
{-4,0,0,0},
{0,1,-1,0},
{0,-1,1,0},
{0,0,0,0}
};


/*!
	The constructor
*/
HPNXObj::HPNXObj(Conf const & theConf,B const NeedHint) : Function(theConf.SysHdl(), NeedHint)
{
	WatchError


	Warn(!(Energy::e().ModelId() != EnergyCls<HpEnergy>::ModelId ||
      Energy::e().ModelId() != EnergyCls<MjEnergy>::ModelId ||
      Energy::e().ModelId() != EnergyCls<BarerraEnergy>::ModelId), eEnergyMismatch);

	Dim tDimen = Space::Dimen();
	Hdl tSysHdl = theConf.SysHdl();
	Dim tSpan = Protein::p().Span();
	Dim tLength = Protein::p().Length();
	block1<Prm,nmm> tResult(tSpan * tLength / 2);

	Int tSqrNeighDist = Lattice::l().SqrNeighDist();
	for(Pos tPos1 = 0; tPos1 < tSpan; ++tPos1)
	{
		Typ tTyp1 = Protein::p().Type(tPos1);

		for(Pos tPos2 = tPos1 + 2; tPos2 < tLength; ++tPos2)
		{
			Typ tTyp2 = Protein::p().Type(tPos2);

			Eng tEng = HPmatrix[toHPNX[tTyp1]][toHPNX[tTyp2]];  // reads from the array

            if(!tEng) continue;

            Prm tDiff[tDimen], tSqrd[tDimen];
			for(Cmp tCmp = 0; tCmp < tDimen; ++tCmp)
			{
				Prm tVar1 = Prm(TermVar,tPos1 * tDimen + tCmp);
				Prm tVar2 = Prm(TermVar,tPos2 * tDimen + tCmp);
				tDiff[tCmp] = Prm(TermFunc,BsubXiFeVi::def(Xv, tSysHdl, tVar1, tVar2));
				tSqrd[tCmp] = Prm(TermFunc,UsqrXiFeVi::def(Xv, tSysHdl, tDiff[tCmp]));
			}
			Prm tDist = Prm(TermFunc,SumXiFeVi::def(Xv, tSysHdl, tSqrd, tDimen));
			Prm tEqu = Prm(TermFunc,UequXiFcMi::def(Xv, tSysHdl, tDist, UequXiFcMi::bind(tSqrNeighDist)));
			tResult.insert(Prm(TermFunc,UifXiKiFcMi::def(Xm, tSysHdl, tEqu,
                                       UifXiKiFcMi::bind(castInt(0), tEng)),MetricAsEvalMin));
		}
	}
    if(NeedHint)
    {
        mFuncHdl = SumXiEFcMiHn::def(Xm | EvalMin, tSysHdl, tResult.items(), tResult.itemCount());
		mMetricRec = &SumXiEFcMiHn::refc(tSysHdl, mFuncHdl).MetricRec();
    }
    else
    {

        mFuncHdl = SumXiFcMi::def(Xm, tSysHdl, tResult.items(), tResult.itemCount());
        mMetricRec = &SumXiFcMi::refc(tSysHdl, mFuncHdl).MetricRec();
    }

	CatchError
}



/*!
	The destructor.
*/
HPNXObj::~HPNXObj()
{
	WatchError
	// nothing to be done.
	CatchError
}



/*!
	The duplicator.
*/
HPNXObj::HPNXObj(HPNXObj const & that) : Function(that)
{
	WatchError
	Alert(&that, eUndefDuplicator);
	CatchError
}



/*!
	The assigner.
*/
HPNXObj const & HPNXObj::operator = (HPNXObj const & that)
{
	WatchError
	Alert(&that, eUndefAssigner);
	return *this;
	CatchError
}


closePlatypusSpace
