/*!
	@file hpmodel/rgyrationobj.cpp
	@brief The implementation file for RGyrationObj class.
	@details This is the implementation file for RGyrationObj class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#include "pspl/functions/rgyrationobj.hpp"



openPlatypusSpace



/*!
	The constructor
*/
RGyrationObj::RGyrationObj(Conf const & theConf, B const NeedHint) :
		Function(theConf.SysHdl(), NeedHint)
{
	WatchError
	//Alert(Energy::e().ModelId() != EnergyCls<HpEnergy>::ModelId, eEnergyMismatch);

	Dim tDimen = Space::Dimen();
	Hdl tSysHdl = theConf.SysHdl();
	Dim tSpan = Protein::p().Span();
	Dim tLength = Protein::p().Length();
	block1<Prm,nmm> tResult(tSpan * tLength / 2);


    block1<Prm,nmm> cMass(tDimen);

    Int exptRG=(Int)pow(exp(log(Protein::p().Length())*0.33+1.11)/3.8,2)*Protein::p().Length();
    //cout << "rgyration:"<<exptRG <<endl;


    for(Cmp tCmp = 0; tCmp < tDimen; ++tCmp)
    {
        block1<Prm,xmm> corePrms;
        Dim totalMass = 0;
        for(Pos tPos = 0; tPos < tSpan; ++tPos)
        {
            Prm temp = Prm(TermVar,tPos * tDimen + tCmp);
            corePrms.insertMem(Prm(TermFunc,
            UmultXiFeVi::def(Xv,tSysHdl,temp,UmultXiFeVi::bind(Protein::p().Mass(tPos)))));
            totalMass+=Protein::p().Mass(tPos);
            //cout<<Protein::p().Mass(tPos) <<endl;
        }
        // now sum all the points
        Prm SumComp = Prm(TermFunc,SumXiFeVi::def(Xv, tSysHdl, corePrms.items(), corePrms.itemCount()));
        Prm CompVal = Prm(TermFunc,UdivXiFeVi::def(Xv,tSysHdl, SumComp,UdivXiFeVi::bind(totalMass)));
        cMass.insert(CompVal);
    }

    for(Pos tPos = 0; tPos < tSpan; ++tPos)
	{
		Prm tDiff[tDimen], tSqrd[tDimen];
        for(Cmp tCmp = 0; tCmp < tDimen; ++tCmp)
        {
            Prm tVar1 = Prm(TermVar,tPos * tDimen + tCmp);
            Prm tVar2 = Prm(cMass[tCmp]);
            tDiff[tCmp] = Prm(TermFunc,BsubXiFeVi::def(Xv, tSysHdl, tVar1, tVar2));
            tSqrd[tCmp] = Prm(TermFunc,UsqrXiFeVi::def(Xv, tSysHdl, tDiff[tCmp]));
        }
        tResult.insert(Prm(TermFunc,SumXiFeVi::def(Xv, tSysHdl, tSqrd, tDimen),ValueAsEvalMin));
	}
    //mFuncHdl = SumXiFcMi::def(Xv, tSysHdl, tResult.items(), tResult.size());
	//mMetricRec = &SumXiFcMi::ref(tSysHdl, mFuncHdl).MetricRec();

    Prm sumPrm=Prm(TermFunc,SumXiFeVi::def(Xv , tSysHdl, tResult.items(),
                                              tResult.itemCount()),ValueAsEvalMin);



    /*if(NeedHint)
    {
        mFuncHdl =UdiffXiEFcMiHn::def(Xv|EvalMin,tSysHdl,sumPrm,UDiffXiEFcMiHn::bind(exptRG));
        mMetricRec = &UdiffXiEFcMiHn::refc(tSysHdl,mFuncHdl).MetricRec();
    }
    else*/
    {
        mFuncHdl = UdiffXiFcMi::def(Xv,tSysHdl,sumPrm,UdiffXiFcMi::bind(exptRG));
        mMetricRec = &UdiffXiFcMi::refc(tSysHdl,mFuncHdl).MetricRec();
    }
    /*if (NeedHint)
	{
		mFuncHdl = SumXiEFcMiHn::def(Xv | EvalMin, tSysHdl, tResult.items(), tResult.itemCount());
		mMetricRec = &SumXiEFcMiHn::refc(tSysHdl, mFuncHdl).MetricRec();
	}
	else
	{
		mFuncHdl = SumXiFcMi::def(Xv, tSysHdl, tResult.items(), tResult.itemCount());
		mMetricRec = &SumXiFcMi::refc(tSysHdl, mFuncHdl).MetricRec();
	}*/
	CatchError
}


/*!
	The destructor.
*/
RGyrationObj::~RGyrationObj()
{
	WatchError
	// nothing to be done.
	CatchError
}



/*!
	The duplicator.
*/
RGyrationObj::RGyrationObj(RGyrationObj const & that) : Function(that)
{
	WatchError
	Alert(&that, eUndefDuplicator);
	CatchError
}



/*!
	The assigner.
*/
RGyrationObj const & RGyrationObj::operator = (RGyrationObj const & that)
{
	WatchError
	Alert(&that, eUndefAssigner);
	return *this;
	CatchError
}



closePlatypusSpace

