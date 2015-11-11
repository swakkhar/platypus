/*!
	@file OriginDistObj.cpp
	@brief The implementation file for OriginDistObj class.
	@details This is the implementation file for OriginDistObj class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#include "pspl/functions/origindistobj.hpp"



openPlatypusSpace



/*!
	The constructor
*/
OriginDistObj::OriginDistObj(Conf const & theConf, B const NeedHint) :
		Function(theConf.SysHdl(), NeedHint)
{
	WatchError
	//Alert(Energy::e().ModelId() != EnergyCls<HpEnergy>::ModelId, eEnergyMismatch);

	Dim tDimen = Space::Dimen();
	Hdl tSysHdl = theConf.SysHdl();
	Dim tSpan = Protein::p().Span();
	Dim tLength = Protein::p().Length();
	block1<Prm,nmm> tResult(tSpan * tLength / 2);


    /*nblock<Prm> centroid(tDimen);

    for(Cmp tCmp = 0; tCmp < tDimen; ++tCmp)
    {
        xblock<Prm> corePrms;
        for(Pos tPos = 0; tPos < tSpan; ++tPos)
        {
            //Typ tTyp = Protein::p().Type(tPos);
            //if (tTyp != Energy::e().Type(2, 'H')) continue;
            corePrms.annex(tPos * tDimen + tCmp);

        }
        // now sum all the points
        Prm SumComp = SumXiFeVi::def(Xv, tSysHdl, corePrms.items(), corePrms.size());
        Prm CompVal = UdivXiFeVi::def(Xv,tSysHdl, SumComp,UdivXiFeVi::bind(corePrms.size()));
        centroid.append(CompVal);
    }*/

    for(Pos tPos = 0; tPos < tSpan; ++tPos)
	{
		//Typ tTyp = Protein::p().Type(tPos);
		//if (tTyp != Energy::e().Type(2, 'H')) continue;

        //Prm tDiff[tDimen], tSqrd[tDimen];
        Prm tSqrd[tDimen];
        for(Cmp tCmp = 0; tCmp < tDimen; ++tCmp)
        {
            Prm tVar1 = Prm(TermVar,tPos * tDimen + tCmp);
            //Prm tVar2 = Prm(centroid[tCmp]);
            //tDiff[tCmp] = BsubXiFeVi::def(Xv, tSysHdl, tVar1, tVar2);
            //tSqrd[tCmp] = UsqrXiFeVi::def(Xv, tSysHdl, tDiff[tCmp]);
            tSqrd[tCmp] = Prm(TermFunc,UsqrXiFeVi::def(Xv, tSysHdl, tVar1));


        }
        tResult.insert(Prm(TermFunc,SumXiFeVi::def(Xv, tSysHdl, tSqrd, tDimen),ValueAsEvalMin));
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
OriginDistObj::~OriginDistObj()
{
	WatchError
	// nothing to be done.
	CatchError
}



/*!
	The duplicator.
*/
OriginDistObj::OriginDistObj(OriginDistObj const & that) : Function(that)
{
	WatchError
	Alert(&that, eUndefDuplicator);
	CatchError
}



/*!
	The assigner.
*/
OriginDistObj const & OriginDistObj::operator = (OriginDistObj const & that)
{
	WatchError
	Alert(&that, eUndefAssigner);
	return *this;
	CatchError
}



closePlatypusSpace
