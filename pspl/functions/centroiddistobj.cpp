/*!
	@file hpmodel/centroidistobj.cpp
	@brief The implementation file for CentroidDistObj class.
	@details This is the implementation file for CentroidDistObj class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#include "pspl/functions/centroiddistobj.hpp"



openPlatypusSpace



/*!
	The constructor
*/
CentroidDistObj::CentroidDistObj(Conf const & theConf, B const NeedHint) :
		Function(theConf.SysHdl(), NeedHint)
{
	WatchError
	//Alert(Energy::e().ModelId() != EnergyCls<HpEnergy>::ModelId, eEnergyMismatch);

	Dim tDimen = Space::Dimen();
	Hdl tSysHdl = theConf.SysHdl();
	Dim tSpan = Protein::p().Span();
	Dim tLength = Protein::p().Length();
	block1<Prm,nmm> tResult(tSpan * tLength / 2);


    block1<Prm,nmm> centroid(tDimen);

    for(Cmp tCmp = 0; tCmp < tDimen; ++tCmp)
    {
        block1<Prm,xmm> corePrms;
        for(Pos tPos = 0; tPos < tSpan; ++tPos)
        {
            corePrms.insertMem(Prm(TermVar,tPos * tDimen + tCmp));
        }
        // now sum all the points
        Prm SumComp = Prm(TermFunc,SumXiFeVi::def(Xv, tSysHdl, corePrms.items(), corePrms.itemCount()));
        Prm CompVal = Prm(TermFunc,UdivXiFeVi::def(Xv,tSysHdl, SumComp,UdivXiFeVi::bind(corePrms.itemCount())));
        centroid.insert(CompVal);
    }

    for(Pos tPos = 0; tPos < tSpan; ++tPos)
	{
		Prm tDiff[tDimen], tSqrd[tDimen];
        for(Cmp tCmp = 0; tCmp < tDimen; ++tCmp)
        {
            Prm tVar1 = Prm(TermVar,tPos * tDimen + tCmp);
            Prm tVar2 = Prm(centroid[tCmp]);
            tDiff[tCmp] = Prm(TermFunc,BsubXiFeVi::def(Xv, tSysHdl, tVar1, tVar2));
            tSqrd[tCmp] = Prm(TermFunc,UsqrXiFeVi::def(Xv, tSysHdl, tDiff[tCmp]));
        }
        tResult.insert(Prm(TermFunc,SumXiFeVi::def(Xv,
                    tSysHdl, tSqrd, tDimen),ValueAsEvalMin));
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
CentroidDistObj::~CentroidDistObj()
{
	WatchError
	// nothing to be done.
	CatchError
}



/*!
	The duplicator.
*/
CentroidDistObj::CentroidDistObj(CentroidDistObj const & that) : Function(that)
{
	WatchError
	Alert(&that, eUndefDuplicator);
	CatchError
}



/*!
	The assigner.
*/
CentroidDistObj const & CentroidDistObj::operator = (CentroidDistObj const & that)
{
	WatchError
	Alert(&that, eUndefAssigner);
	return *this;
	CatchError
}
closePlatypusSpace
