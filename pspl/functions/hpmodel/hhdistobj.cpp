/*!
	@file hpmodel/hhdistobj.cpp
	@brief The implementation file for HhDistObj class.
	@details This is the implementation file for HhDistObj class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#include "pspl/functions/hpmodel/hhdistobj.hpp"



openPlatypusSpace



/*!
	The constructor
*/
HhDistObj::HhDistObj(Conf const & theConf, B const NeedHint) :
		Function(theConf.SysHdl(), NeedHint)
{
	WatchError
	Alert(Energy::e().ModelId() != EnergyCls<HpEnergy>::ModelId, eEnergyMismatch);

	Dim tDimen = Space::Dimen();
	Hdl tSysHdl = theConf.SysHdl();
	Dim tSpan = Protein::p().Span();
	Dim tLength = Protein::p().Length();
	block1<Prm,nmm> tResult(tSpan * tLength / 2);

	for(Pos tPos1 = 0; tPos1 < tSpan; ++tPos1)
	{
		Typ tTyp1 = Protein::p().Type(tPos1);
		if (tTyp1 != Energy::e().Type(2, 'H')) continue;
		for(Pos tPos2 = tPos1 + 2; tPos2 < tLength; ++tPos2)
		{
			Typ tTyp2 = Protein::p().Type(tPos2);
			if (tTyp2 != Energy::e().Type(2, 'H')) continue;

            Prm tDiff[tDimen], tSqrd[tDimen];
			for(Cmp tCmp = 0; tCmp < tDimen; ++tCmp)
			{
				Prm tVar1 = Prm(TermVar,tPos1 * tDimen + tCmp);
				Prm tVar2 = Prm(TermVar,tPos2 * tDimen + tCmp);
				tDiff[tCmp] = Prm(TermFunc,BsubXiFeVi::def(Xv, tSysHdl, tVar1, tVar2));
				tSqrd[tCmp] = Prm(TermFunc,UsqrXiFeVi::def(Xv, tSysHdl, tDiff[tCmp]));
			}
			tResult.insert(Prm(TermFunc,SumXiFeVi::def(Xv, tSysHdl, tSqrd, tDimen),ValueAsEvalMin));
		}
	}
	//mFuncHdl = SumXiFcMi::def(Xv, tSysHdl, tResult.items(), tResult.size());
	//mMetricRec = &SumXiFcMi::ref(tSysHdl, mFuncHdl).MetricRec();
	if (NeedHint)
	{
		mFuncHdl = SumXiEFcMiHn::def(Xv | EvalMin, tSysHdl, tResult.items(), tResult.itemCount());
		mMetricRec = &SumXiEFcMiHn::refc(tSysHdl, mFuncHdl).MetricRec();
	}
	else
	{
		mFuncHdl = SumXiFcMi::def(Xv, tSysHdl, tResult.items(), tResult.itemCount());
		mMetricRec = &SumXiFcMi::refc(tSysHdl, mFuncHdl).MetricRec();
	}
	CatchError
}


/*!
	The destructor.
*/
HhDistObj::~HhDistObj()
{
	WatchError
	// nothing to be done.
	CatchError
}



/*!
	The duplicator.
*/
HhDistObj::HhDistObj(HhDistObj const & that) : Function(that)
{
	WatchError
	Alert(&that, eUndefDuplicator);
	CatchError
}



/*!
	The assigner.
*/
HhDistObj const & HhDistObj::operator = (HhDistObj const & that)
{
	WatchError
	Alert(&that, eUndefAssigner);
	return *this;
	CatchError
}



closePlatypusSpace
