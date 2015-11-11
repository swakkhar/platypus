/*!
	@file inits/MotifInit.cpp
	@brief The implementation file for MotifInit class.
	@details This is the implementation file for MotifInit class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#include "pspl/inits/motifinit.hpp"
#include "pspl/globals/conf.hpp"
#include "pspl/globals/change.hpp"
#include "pspl/globals/lattice.hpp"



openPlatypusSpace

Z helix[]={2,8,7,1};
Z sheet[]={0,6};

Z motifDirs[][4]={
    {0,0,0,0},
    {0,6,0,6},
    {2,8,7,1}
};

/*!
	The costructor.
*/
MotifInit::MotifInit(Rnd & theRnd) : mRnd(theRnd)
{
    WatchError
    // do nothing
    CatchError
}


/*!
	The Destructor
*/
MotifInit::~MotifInit()
{
    WatchError
    // do nothing
    CatchError
}


/*!
	The Dupllicator
*/
MotifInit::MotifInit(MotifInit const & that) :
	Init(that), mRnd(that.mRnd)
{
    WatchError
    Alert(&that, eUndefDuplicator);
    CatchError
}


/*!
	The assigner.
*/

MotifInit const & MotifInit::operator= (MotifInit const & that)
{
    WatchError
	Alert(&that, eUndefAssigner);
    return *this;
    CatchError
}


/*!
	Compute and store a possible initialisation change.
*/
void MotifInit::compute(Conf & theConf)
{
    WatchError

	//Dim tNeighCount = Lattice::l().NeighCount();
	Dim tLength = Protein::p().Length();
    hset<Point,xmmh> tPoints(tLength);

	Idx motifIdx=0;
	Idx motifCode=0;


	mChange.reset();
	tPoints.clear();

	Point tPoint(0), sPoint;
	mChange.add(0,tPoint);
	tPoints.insertBll(tPoint);

    motifCode=theConf.getMotifAtPos(0);

	Pos tPos = 1;
	do
	{
	    Dir tDir=0;
        Idx tCode = theConf.getMotifAtPos(tPos);

        if(motifCode!=tCode) motifIdx=0;
        else motifIdx++;

        motifCode=tCode;


        if(motifIdx==0) tDir=0;
        else tDir=motifDirs[motifCode][(motifIdx-1)%4];

        sPoint=tPoint+Lattice::l().DirVec(tDir);

        tPoint=sPoint;

        mChange.add(tPos,tPoint);
        tPoints.insertBll(tPoint);

        ++tPos;


	}
	while (tPos < tLength);
    CatchError
}



closePlatypusSpace
