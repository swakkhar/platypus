/*!
	@file ContactCountObj.cpp
	@brief The implementation file for ContactCountObj class.
	@details This is the implementation file for ContactCountObj class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#include "pspl/functions/contactcountobj.hpp"



openPlatypusSpace



/*!
	The constructor
*/
ContactCountObj::ContactCountObj(Conf const & theConf, B const NeedHint) :
		Function(theConf.SysHdl(), NeedHint)
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
		//Typ tTyp1 = Protein::p().Type(tPos1);

		for(Pos tPos2 = tPos1 + 2; tPos2 < tLength; ++tPos2)
		{
			//Typ tTyp2 = Protein::p().Type(tPos2);

			//Eng tEng = Energy::e().Level(tTyp1, tTyp2);
            //Pos tEng=tPos2-tPos1;
            //if(!tEng) continue;

            Prm tDiff[tDimen], tSqrd[tDimen];
			for(Cmp tCmp = 0; tCmp < tDimen; ++tCmp)
			{
				Prm tVar1 = Prm(TermVar,tPos1 * tDimen + tCmp);
				Prm tVar2 = Prm(TermVar,tPos2 * tDimen + tCmp);
				tDiff[tCmp] = Prm(TermFunc, BsubXiFeVi::def(Xv, tSysHdl, tVar1, tVar2));
				tSqrd[tCmp] = Prm(TermFunc, UsqrXiFeVi::def(Xv, tSysHdl, tDiff[tCmp]));
			}
			Prm tDist = Prm(TermFunc,SumXiFeVi::def(Xv, tSysHdl, tSqrd, tDimen));
			Prm tEqu = Prm(TermFunc,UequXiFcMi::def(Xv, tSysHdl, tDist, UequXiFcMi::bind(tSqrNeighDist)));
			tResult.insert(Prm(TermFunc,UifXiKiFcMi::def(Xm, tSysHdl, tEqu,
                                       UifXiKiFcMi::bind(castInt(0), castInt(1))),MetricAsEvalMin));
		}
	}
	//mFuncHdl = SumXiFcMi::def(Xm, tSysHdl, tResult.items(), tResult.size());
	//mMetricRec = &SumXiFcMi::ref(tSysHdl, mFuncHdl).MetricRec();
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
ContactCountObj::~ContactCountObj()
{
	WatchError
	// nothing to be done.
	CatchError
}



/*!
	The duplicator.
*/
ContactCountObj::ContactCountObj(ContactCountObj const & that) : Function(that)
{
	WatchError
	Alert(&that, eUndefDuplicator);
	CatchError
}



/*!
	The assigner.
*/
ContactCountObj const & ContactCountObj::operator = (ContactCountObj const & that)
{
	WatchError
	Alert(&that, eUndefAssigner);
	return *this;
	CatchError
}



closePlatypusSpace
