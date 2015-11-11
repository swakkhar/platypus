/*!
	@file hpmodel/FreeFtCns.cpp
	@brief The implementation file for FreeFtCns class.
	@details This is the implementation file for FreeFtCns class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#include "pspl/functions/freeftcns.hpp"



openPlatypusSpace



/*!
	The constructor
*/
FreeFtCns::FreeFtCns(Conf const & theConf, B const NeedHint) :
		Function(theConf.SysHdl(), NeedHint)
{
	WatchError
	//Warn(Energy::e().ModelId() != EnergyCls<HpEnergy>::ModelId, eEnergyMismatch);

	Dim tDimen = Space::Dimen();
	Hdl tSysHdl = theConf.SysHdl();
	Dim tSpan = Protein::p().Span();
	Dim tLength = Protein::p().Length();
	Dim tNeighborCount = Lattice::l().NeighCount();

	Int tSqrNeighDist = Lattice::l().SqrNeighDist();

	block1<Prm,nmm> tSubs(tLength);

	for(Pos tPos1 = 0; tPos1 < tLength; ++tPos1)
	{
        block1<Prm,nmm> tResult(tSpan);

		for(Pos tPos2 = 0; tPos2 < tLength; ++tPos2)
		{
		    if(tPos1==tPos2) continue;

            Prm tDiff[tDimen], tSqrd[tDimen];
			for(Cmp tCmp = 0; tCmp < tDimen; ++tCmp)
			{
				Prm tVar1 = Prm(TermVar,tPos1 * tDimen + tCmp);
				Prm tVar2 = Prm(TermVar,tPos2 * tDimen + tCmp);
				tDiff[tCmp] = Prm(TermFunc,BsubXiFeVi::def(Xv, tSysHdl, tVar1, tVar2));
				tSqrd[tCmp] = Prm(TermFunc,UsqrXiFeVi::def(Xv, tSysHdl, tDiff[tCmp]));
			}
			Prm tDist = Prm(TermFunc,SumXiFeVi::def(Xv, tSysHdl, tSqrd, tDimen));
			tResult.insert(Prm(TermFunc,UequXiFeVi::def(Xv, tSysHdl, tDist, UequXiFcMi::bind(tSqrNeighDist))));
		}
        Prm sumTemp = Prm(TermFunc,SumXiFeVi::def(Xv, tSysHdl, tResult.items(), tResult.itemCount()));
		// now subtract this sum of tResults from 10/11 which will give a freecount
        Prm minus = Prm(TermFunc,UmultXiFeVi::def(Xv,tSysHdl,sumTemp,UmultXiFeVi::bind(-1)));

        tSubs.insert(Prm(TermFunc,UaddXiFcMi::def(Xv,tSysHdl,minus,UaddXiFcMi::bind(tNeighborCount)), MetricAsEvalMin));
	}
	if (NeedHint)
	{
		mFuncHdl = SumXiEFcMiHn::def(Xm | EvalMin, tSysHdl, tSubs.items(), tSubs.itemCount());
		mMetricRec = &SumXiEFcMiHn::refc(tSysHdl, mFuncHdl).MetricRec();
	}
	else
	{
		mFuncHdl = SumXiFcMi::def(Xm, tSysHdl, tSubs.items(), tSubs.itemCount());
		mMetricRec = &SumXiFcMi::refc(tSysHdl, mFuncHdl).MetricRec();
	}
	CatchError
}



/*!
	The destructor.
*/
FreeFtCns::~FreeFtCns()
{
	WatchError
	// nothing to be done.
	CatchError
}



/*!
	The duplicator.
*/
FreeFtCns::FreeFtCns(FreeFtCns const & that) : Function(that)
{
	WatchError
	Alert(&that, eUndefDuplicator);
	CatchError
}



/*!
	The assigner.
*/
FreeFtCns const & FreeFtCns::operator = (FreeFtCns const & that)
{
	WatchError
	Alert(&that, eUndefAssigner);
	return *this;
	CatchError
}



closePlatypusSpace
