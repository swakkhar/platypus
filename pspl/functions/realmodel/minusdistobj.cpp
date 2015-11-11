/*!
	@file real/minusdistobj.cpp
	@brief The implementation file for MjMinusDistObj class.
	@details This is the implementation file for MinusDistObj class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#include "pspl/functions/realmodel/minusdistobj.hpp"



openPlatypusSpace



/*!
	The constructor
*/
MinusDistObj::MinusDistObj(Conf const & theConf, B const NeedHint) :
		Function(theConf.SysHdl(), NeedHint)
{
	WatchError
	Warn(!(Energy::e().ModelId() != EnergyCls<MjEnergy>::ModelId ||
      Energy::e().ModelId() != EnergyCls<BarerraEnergy>::ModelId), eEnergyMismatch);

    Z modelCutoff=(Energy::e().ModelId() == EnergyCls<MjEnergy>::ModelId)?0:-1000;

	Dim tDimen = Space::Dimen();
	Hdl tSysHdl = theConf.SysHdl();
	Dim tSpan = Protein::p().Span();
	Dim tLength = Protein::p().Length();
	block1<Prm,nmm> tResult(tSpan * tLength / 2);

	for(Pos tPos1 = 0; tPos1 < tSpan; ++tPos1)
	{
		Typ tTyp1 = Protein::p().Type(tPos1);

		//if (tTyp1 != Energy::e().Type(2, 'H')) continue;
		for(Pos tPos2 = tPos1 + 2; tPos2 < tLength; ++tPos2)
		{
			Typ tTyp2 = Protein::p().Type(tPos2);

			//if (tTyp2 != Energy::e().Type(2, 'H')) continue;
            if(Energy::e().Level(tTyp1,tTyp2)>=modelCutoff) continue;
            //cout << tPos1 << " "<<tPos2<<" "<<tTyp1 << " "<< tTyp2<<endl;
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
	if (NeedHint)
	{
		mFuncHdl = SumXiEFcMiHn::def(Xv | EvalMin, tSysHdl, tResult.items(), tResult.itemCount());
		mMetricRec = &SumXiEFcMiHn::refc(tSysHdl, mFuncHdl).MetricRec();
	}
	else
	{
		mFuncHdl = SumXiFcMi::def(Xv, tSysHdl, tResult.items(), tResult.itemCount());
		mMetricRec = &SumXiFcMi::refc(tSysHdl, mFuncHdl).MetricRec();
	}	CatchError
}


/*!
	The destructor.
*/
MinusDistObj::~MinusDistObj()
{
	WatchError
	// nothing to be done.
	CatchError
}



/*!
	The duplicator.
*/
MinusDistObj::MinusDistObj(MinusDistObj const & that) : Function(that)
{
	WatchError
	Alert(&that, eUndefDuplicator);
	CatchError
}



/*!
	The assigner.
*/
MinusDistObj const & MinusDistObj::operator = (MinusDistObj const & that)
{
	WatchError
	Alert(&that, eUndefAssigner);
	return *this;
	CatchError
}



closePlatypusSpace
