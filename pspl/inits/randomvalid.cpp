/*!
	@file inits/randomvalid.cpp
	@brief The implementation file for RandomValid class.
	@details This is the implementation file for RandomValid class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#include "pspl/inits/randomvalid.hpp"
#include "pspl/globals/conf.hpp"
#include "pspl/globals/change.hpp"
#include "pspl/globals/lattice.hpp"



openPlatypusSpace



/*!
	The costructor.
*/
RandomValid::RandomValid(Rnd & theRnd) : mRnd(theRnd)
{
    WatchError
    // do nothing
    CatchError
}


/*!
	The Destructor
*/
RandomValid::~RandomValid()
{
    WatchError
    // do nothing
    CatchError
}


/*!
	The Dupllicator
*/
RandomValid::RandomValid(RandomValid const & that) :
	Init(that), mRnd(that.mRnd)
{
    WatchError
    Alert(&that, eUndefDuplicator);
    CatchError
}


/*!
	The assigner.
*/

RandomValid const & RandomValid::operator= (RandomValid const & that)
{
    WatchError
	Alert(&that, eUndefAssigner);
    return *this;
    CatchError
}


/*!
	Compute and store a possible initialisation change.
*/
void RandomValid::compute(Conf & theConf)
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
