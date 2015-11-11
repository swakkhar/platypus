/*!
	@file RMSDObj.cpp
	@brief The implementation file for RMSDObj class.
	@details This is the implementation file for RMSDObj class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#include "pspl/functions/rmsdobj.hpp"



openPlatypusSpace



/*!
	The constructor
*/
RMSDObj::RMSDObj(Conf const & theConf, B const NeedHint,block1<RPoint,xmm> &arr) :
		Function(theConf.SysHdl(), NeedHint),mNativePoints(arr)
{
		WatchError


	Warn(!(Energy::e().ModelId() != EnergyCls<HpEnergy>::ModelId ||
      Energy::e().ModelId() != EnergyCls<MjEnergy>::ModelId ||
      Energy::e().ModelId() != EnergyCls<BarerraEnergy>::ModelId), eEnergyMismatch);


	Dim tDimen = Space::Dimen();
	Hdl tSysHdl = theConf.SysHdl();
	Dim tSpan = Protein::p().Span();
	Dim tLength = Protein::p().Length();
	block1<Prm,xmm> tResult;

	for(Pos tPos1 = 0; tPos1 < tSpan; ++tPos1)
	{
		for(Pos tPos2 = tPos1 + 1; tPos2 < tLength; ++tPos2)
		{
		    //cerr << "native dist before: " << errno <<endl;
            //Int realDist = Tuple::dist(mNativePoints[tPos1],mNativePoints[tPos2]);
            //cerr << "native dist after: "<< errno <<endl;
            //if(errno==33) exit(0);

            double sqrDist=(mNativePoints[tPos1].x-mNativePoints[tPos2].x)*(mNativePoints[tPos1].x-mNativePoints[tPos2].x)
                            +(mNativePoints[tPos1].y-mNativePoints[tPos2].y)*(mNativePoints[tPos1].y-mNativePoints[tPos2].y)
                            +(mNativePoints[tPos1].z-mNativePoints[tPos2].z)*(mNativePoints[tPos1].z-mNativePoints[tPos2].z);

            double rDist = sqrt(sqrDist)*10;
            Int realDist = (Int)rDist;

           // cout << realDist<<endl;

            Prm tDiff[tDimen], tSqrd[tDimen];
			for(Cmp tCmp = 0; tCmp < tDimen; ++tCmp)
			{
				Prm tVar1 = Prm(TermVar,tPos1 * tDimen + tCmp);
				Prm tVar2 = Prm(TermVar,tPos2 * tDimen + tCmp);
				tDiff[tCmp] = Prm(TermFunc,BsubXiFeVi::def(Xv, tSysHdl, tVar1, tVar2));
				tSqrd[tCmp] = Prm(TermFunc,UsqrXiFeVi::def(Xv, tSysHdl, tDiff[tCmp]));
			}
			Prm tDist = Prm(TermFunc,SumXiFeVi::def(Xv, tSysHdl, tSqrd, tDimen));
			// tDist contains the square of the model distance
			Prm mrDist = Prm(TermFunc,UsqrtXiFeVi::def(Xv,tSysHdl,tDist));
            // multiply it by 2.69*1000 = 2690
            Prm modelDist = Prm(TermFunc,UmultXiFeVi::def(Xv, tSysHdl, mrDist, UmultXiFeVi::bind(10*Lattice::l().latticeConverter())));

            Prm tSub = Prm(TermFunc,UsubXiFeVi::def(Xv, tSysHdl, modelDist, UsubXiFeVi::bind(realDist)));

			tResult.insertMem(Prm(TermFunc,UsqrXiFeVi::def(Xv, tSysHdl, tSub),ValueAsEvalMin));
		}
	}
	if(NeedHint)
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
RMSDObj::~RMSDObj()
{
	WatchError
	// nothing to be done.
	CatchError
}



/*!
	The duplicator.
*/
RMSDObj::RMSDObj(RMSDObj const & that) : Function(that),mNativePoints(that.mNativePoints)
{
	WatchError
	Alert(&that, eUndefDuplicator);
	CatchError
}



/*!
	The assigner.
*/
RMSDObj const & RMSDObj::operator = (RMSDObj const & that)
{
	WatchError
	Alert(&that, eUndefAssigner);
	return *this;
	CatchError
}
closePlatypusSpace
