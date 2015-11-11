/*!
	@file hpmodel/MjObj.cpp
	@brief The implementation file for MjObj class.
	@details This is the implementation file for MjObj class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#include "pspl/functions/mjobj.hpp"



openPlatypusSpace


const Eng MJmatrix[20][20]=
{
{ -1060, 190, -230, 160, -80, 60, 80, 40, 0, -80, 190, -20, 50, 130, 690, 30, -190, 240, 710, 0 },
{ 190, 40, -420, -280, -200, -140, -670, -130, 250, 190, 190, 140, 460, 80, 440, 650, 990, 310, 0, -340 },
{ -230, -420, -440, -190, -300, -220, -160, 0, 30, 380, 310, 290, 490, 180, 270, 390, -160, 410, 440, 200 },
{ 160, -280, -190, -220, -410, -250, 20, 110, -220, 250, 140, 210, 360, 530, 350, 590, 490, 420, 360, 250 },
{ -80, -200, -300, -410, -270, -290, -90, 240, -10, 230, 200, 250, 260, 300, 430, 670, 160, 350, 190, 420 },
{ 60, -140, -220, -250, -290, -290, -170, 20, -100, 160, 250, 180, 240, 500, 340, 580, 190, 300, 440, 90 },
{ 80, -670, -160, 20, -90, -170, -120, -40, -90, 180, 220, 340, 80, 60, 290, 240, -120, -160, 220, -280 },
{ 40, -130, 0, 110, 240, 20, -40, -60, 90, 140, 130, 90, -200, -200, -100, 0, -340, -250, -210, -330 },
{ 0, 250, 30, -220, -10, -100, -90, 90, -130, -70, -90, -60, 80, 280, 260, 120, 340, 430, 140, 100 },
{ -80, 190, 380, 250, 230, 160, 180, 140, -70, -380, -260, -160, -60, -140, 250, -220, 200, -40, 110, -110 },
{ 190, 190, 310, 140, 200, 250, 220, 130, -90, -260, 30, -80, -140, -110, 0, -290, -190, -350, -90, -70 },
{ -20, 140, 290, 210, 250, 180, 340, 90, -60, -160, -80, 200, -140, -140, -260, -310, -50, 170, -130, 10 },
{ 50, 460, 490, 360, 260, 240, 80, -200, 80, -60, -140, -140, 290, -250, -170, -170, -20, -520, -380, -420 },
{ 130, 80, 180, 530, 300, 500, 60, -200, 280, -140, -110, -140, -250, -530, -320, -300, -240, -140, -330, -180 },
{ 690, 440, 270, 350, 430, 340, 290, -100, 260, 250, 0, -260, -170, -320, -30, -150, -450, -740, -970, -100 },
{ 30, 650, 390, 590, 670, 580, 240, 0, 120, -220, -290, -310, -170, -300, -150, 40, -390, -720, -760, 40 },
{ -190, 990, -160, 490, 160, 190, -120, -340, 340, 200, -190, -50, -20, -240, -450, -390, -290, -120, 220, -210 },
{ 240, 310, 410, 420, 350, 300, -160, -250, 430, -40, -350, 170, -520, -140, -740, -720, -120, 110, 750, -380 },
{ 710, 0, 440, 360, 190, 440, 220, -210, 140, 110, -90, -130, -380, -330, -970, -760, 220, 750, 250, 110 },
{ 0, -340, 200, 250, 420, 90, -280, -330, 100, -110, -70, 10, -420, -180, -100, 40, -210, -380, 110, 260 }
};


/*!
	The constructor
*/
MjObj::MjObj(Conf const & theConf,B const NeedHint) : Function(theConf.SysHdl(), NeedHint)
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
		Typ tTyp1 = Protein::p().Type(tPos1);

		for(Pos tPos2 = tPos1 + 2; tPos2 < tLength; ++tPos2)
		{
			Typ tTyp2 = Protein::p().Type(tPos2);

			Eng tEng = MJmatrix[tTyp1][tTyp2];  // reads from the array

            if(!tEng) continue;

            Prm tDiff[tDimen], tSqrd[tDimen];
			for(Cmp tCmp = 0; tCmp < tDimen; ++tCmp)
			{
				Prm tVar1 = Prm(TermVar,tPos1 * tDimen + tCmp);
				Prm tVar2 = Prm(TermVar,tPos2 * tDimen + tCmp);
				tDiff[tCmp] = Prm(TermFunc,BsubXiFeVi::def(Xv, tSysHdl, tVar1, tVar2));
				tSqrd[tCmp] = Prm(TermFunc,UsqrXiFeVi::def(Xv, tSysHdl, tDiff[tCmp]));
			}
			Prm tDist = Prm(TermFunc,SumXiFeVi::def(Xv, tSysHdl, tSqrd, tDimen));
			Prm tEqu = Prm(TermFunc,UequXiFcMi::def(Xv, tSysHdl, tDist, UequXiFcMi::bind(tSqrNeighDist)));
			tResult.insert(Prm(TermFunc,UifXiKiFcMi::def(Xm, tSysHdl, tEqu,
                                       UifXiKiFcMi::bind(castInt(0), tEng)),MetricAsEvalMin));
		}
	}
    if(NeedHint)
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
MjObj::~MjObj()
{
	WatchError
	// nothing to be done.
	CatchError
}



/*!
	The duplicator.
*/
MjObj::MjObj(MjObj const & that) : Function(that)
{
	WatchError
	Alert(&that, eUndefDuplicator);
	CatchError
}



/*!
	The assigner.
*/
MjObj const & MjObj::operator = (MjObj const & that)
{
	WatchError
	Alert(&that, eUndefAssigner);
	return *this;
	CatchError
}


closePlatypusSpace
