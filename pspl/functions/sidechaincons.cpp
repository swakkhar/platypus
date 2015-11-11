/*!
	@file functions/SideChainCns.cpp
	@brief The implementation file for SideChainCns class.
	@details This is the implementation file for SideChainCns class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#include "pspl/functions/sidechaincons.hpp"



openPlatypusSpace



/*!
	The constructor
*/
SideChainCns::SideChainCns(Conf const & theConf, B const NeedHint) :
		Function(theConf.SysHdl(), NeedHint)
{
	WatchError

	Dim tDimen = Space::Dimen();
	Hdl tSysHdl = theConf.SysHdl();
	Dim tSpan = Protein::p().Span();
	Dim tLength = Protein::p().Length();
	block1<Prm,nmm> tResult(tSpan * tLength / 2);

	Int tSqrNeighDist = Lattice::l().SqrNeighDist();
	for(Pos tPos1 = 0; tPos1 < tLength; ++tPos1)
	{
	    Pos  tPos2;
	    if(tPos1<tLength/2)
            tPos2 = tPos1 + 1;
        else
            tPos2 = tPos1 - tLength/2;
		Prm tDiff[tDimen], tSqrd[tDimen];
		for(Cmp tCmp = 0; tCmp < tDimen; ++tCmp)
		{
			Prm tVar1 = Prm(TermVar,tPos1 * tDimen + tCmp);
			Prm tVar2 = Prm(TermVar,tPos2 * tDimen + tCmp);
			tDiff[tCmp] = Prm(TermFunc,BsubXiFeVi::def(Xv, tSysHdl, tVar1, tVar2));
			tSqrd[tCmp] = Prm(TermFunc,UsqrXiFeVi::def(Xv, tSysHdl, tDiff[tCmp]));
		}
		Prm tDist = Prm(TermFunc,SumXiFeVi::def(Xv, tSysHdl, tSqrd, tDimen));
		tResult.insert(Prm(TermFunc,UequXiFcMi::def(Xv, tSysHdl, tDist, UequXiFcMi::bind(tSqrNeighDist)), MetricAsEvalMin));
	}
	cout << tResult.itemCount() <<endl;
	if (NeedHint)
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
SideChainCns::~SideChainCns()
{
	WatchError
	// nothing to be done.
	CatchError
}



/*!
	The duplicator.
*/
SideChainCns::SideChainCns(SideChainCns const & that) : Function(that)
{
	WatchError
	Alert(&that, eUndefDuplicator);
	CatchError
}



/*!
	The assigner.
*/
SideChainCns const & SideChainCns::operator = (SideChainCns const & that)
{
	WatchError
	Alert(&that, eUndefAssigner);
	return *this;
	CatchError
}



closePlatypusSpace
