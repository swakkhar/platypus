/*!
	@file hpmodel/RepCntCns.cpp
	@brief The implementation file for RepCntCns class.
	@details This is the implementation file for RepCntCns class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#include "pspl/functions/realmodel/repcntcns.hpp"



openPlatypusSpace



/*!
	The constructor
*/
RepCntCns::RepCntCns(Conf const & theConf, B const NeedHint) :
		Function(theConf.SysHdl(), NeedHint)
{
	WatchError
	Warn(Energy::e().ModelId() != EnergyCls<HpEnergy>::ModelId, eEnergyMismatch);

	Dim tDimen = Space::Dimen();
	Hdl tSysHdl = theConf.SysHdl();
	Dim tSpan = Protein::p().Span();
	Dim tLength = Protein::p().Length();
	//Dim tNeighborCount = Lattice::l().NeighCount();

	Int tSqrNeighDist = Lattice::l().SqrNeighDist();

	block1<Prm,nmm> tSums(tLength);

	for(Pos tPos1 = 0; tPos1 < tLength; ++tPos1)
	{
        block1<Prm,nmm> tResult(tSpan);
        Typ tTyp1 = Protein::p().Type(tPos1);
		for(Pos tPos2 = 0; tPos2 < tLength; ++tPos2)
		{
		    if(tPos1==tPos2) continue;
            Typ tTyp2 = Protein::p().Type(tPos2);
            if(Energy::e().Level(tTyp1,tTyp2)<=0)continue;
            // if the contact type is negative then its alright

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
        tSums.insert(Prm(TermFunc,SumXiFcMi::def(Xv, tSysHdl, tResult.items(), tResult.itemCount()),MetricAsEvalMin));
		// now subtract this sum of tResults from 10/11 which will give a freecount
        //Prm minus = Prm(TermFunc,UmultXiFeVi::def(Xv,tSysHdl,sumTemp,UmultXiFeVi::bind(-1)));

        //tSubs.insert(Prm(TermFunc,UaddXiFcMi::def(Xv,tSysHdl,minus,UaddXiFcMi::bind(tNeighborCount)), MetricAsEvalMin));
	}
	if (NeedHint)
	{
		mFuncHdl = SumXiEFcMiHn::def(Xm | EvalMin, tSysHdl, tSums.items(), tSums.itemCount());
		mMetricRec = &SumXiEFcMiHn::refc(tSysHdl, mFuncHdl).MetricRec();
	}
	else
	{
		mFuncHdl = SumXiFcMi::def(Xm, tSysHdl, tSums.items(), tSums.itemCount());
		mMetricRec = &SumXiFcMi::refc(tSysHdl, mFuncHdl).MetricRec();
	}
	CatchError
}



/*!
	The destructor.
*/
RepCntCns::~RepCntCns()
{
	WatchError
	// nothing to be done.
	CatchError
}



/*!
	The duplicator.
*/
RepCntCns::RepCntCns(RepCntCns const & that) : Function(that)
{
	WatchError
	Alert(&that, eUndefDuplicator);
	CatchError
}



/*!
	The assigner.
*/
RepCntCns const & RepCntCns::operator = (RepCntCns const & that)
{
	WatchError
	Alert(&that, eUndefAssigner);
	return *this;
	CatchError
}



closePlatypusSpace
