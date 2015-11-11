/*!
	@file inits/randomvalid.cpp
	@brief The implementation file for RandomValidSC class.
	@details This is the implementation file for RandomValidSC class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#include "pspl/inits/randomvalidsc.hpp"
#include "pspl/globals/conf.hpp"
#include "pspl/globals/change.hpp"
#include "pspl/globals/lattice.hpp"



openPlatypusSpace



/*!
	The costructor.
*/
RandomValidSC::RandomValidSC(Rnd & theRnd) : mRnd(theRnd)
{
    WatchError
    // do nothing
    CatchError
}


/*!
	The Destructor
*/
RandomValidSC::~RandomValidSC()
{
    WatchError
    // do nothing
    CatchError
}


/*!
	The Dupllicator
*/
RandomValidSC::RandomValidSC(RandomValidSC const & that) :
	Init(that), mRnd(that.mRnd)
{
    WatchError
    Alert(&that, eUndefDuplicator);
    CatchError
}


/*!
	The assigner.
*/

RandomValidSC const & RandomValidSC::operator= (RandomValidSC const & that)
{
    WatchError
	Alert(&that, eUndefAssigner);
    return *this;
    CatchError
}


/*!
	Compute and store a possible initialisation change.
*/
void RandomValidSC::compute(Conf & theConf)
{
    WatchError

	Dim tNeighCount = Lattice::l().NeighCount();
	Dim tLength = Protein::p().Length();
    hset<Point,xmmh> tPoints(tLength);

	mChange.reset();
	tPoints.clear();

	Point tPoint(0), sPoint;
	mChange.add(0,tPoint);
	tPoints.insertBll(tPoint);

	Pos tPos = 1;
	do
	{
		Idx tIdx = 0;
		do
		{
			++tIdx;
			// if tPos is a side chain take the point from backbone
			if(theConf.isPosSideChain(tPos))
			{
			    sPoint = mChange.Destination(tPos-Protein::p().Length()/2);
			}

            else // the previous point
                sPoint = tPoint;
			Dir tDir = uniform(mRnd, tNeighCount);
			sPoint = sPoint + Lattice::l().DirVec(tDir);
		}
		while (!tPoints.insertBll(sPoint) && tIdx < tNeighCount);
		if (tIdx <= tNeighCount)
		{
			mChange.add(tPos, sPoint);
			tPoint = sPoint;
			++tPos;
		}
		else
		{
			tPoints.removeItr(tPoint);
			tPoint = mChange.Destination(mChange.size() - 1);
			mChange.remove();
			--tPos;
		}
	}
	while (tPos < tLength);
    CatchError
}



closePlatypusSpace
