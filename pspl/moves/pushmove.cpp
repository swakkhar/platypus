/*!
	@file moves/PushMove.cpp
	@brief The implementation file for PushMove class.
	@details This is the implementation file for PushMove class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/

#include "pspl/moves/pushmove.hpp"
#include "pspl/globals/conf.hpp"
#include "pspl/globals/change.hpp"
#include "pspl/globals/lattice.hpp"


openPlatypusSpace


/*!
The costructor
*/
PushMove::PushMove()
{
    WatchError

    // do nothing
    CatchError
}


/*!

The Destructor

*/
PushMove::~PushMove()
{
    WatchError

    // do nothing

    CatchError
}


/*!

The Dupllicator
*/

PushMove::PushMove(PushMove const & that):Move(that)
{
    WatchError

    Alert(&that, eUndefDuplicator);

    CatchError
}


/*!

Assignment Operator
*/

PushMove const & PushMove::operator= (PushMove const & that)
{
    WatchError

    Move::operator=(that);
    return *this;


    CatchError
}


/*!

The compute method
save all possible move results into mChanges
*/
void PushMove::compute(Conf & theConf, Pos const thePos)
{
    WatchError
    mChanges.clear();
	if (!thePos || thePos >= Protein::p().Span())
		return;

	Point const & b = theConf.Coord(thePos - 1);
	Point const & a = theConf.Coord(thePos + 1);

	if (!Lattice::l().areNeighbors(a,b))
		return;

    Change tChange;
    block1<Point,xmm> validPoints;
    for(Int pushDir = -1; pushDir <= 1; pushDir += 2)
    {
        tChange.reset();

		validPoints.clear();

		Point f = (pushDir == -1 ? b : a);
		tChange.add(thePos, f);

		Pos tPos = thePos + pushDir;
        while(tPos && tPos < Protein::p().Span())
        {
			Point n = theConf.Coord(tPos + pushDir);
			for(Dir tDir = 0; tDir < Lattice::l().NeighCount(); ++tDir)
			{
				Point t = n + Lattice::l().DirVec(tDir);
				if(theConf.occupied(t)) continue;
				if(Lattice::l().areNeighbors(t,f))
					validPoints.insertMem(t);
			}
			if (validPoints.itemCount())
				break;
			tChange.add(tPos,n);
			f = n;
			tPos += pushDir;
        }
		if (!validPoints.itemCount())// && (!tPos || tPos == Protein::p().Span()))
			for(Dir tDir = 0; tDir < Lattice::l().NeighCount(); ++tDir)
			{
				Point t = f + Lattice::l().DirVec(tDir);
				if(!theConf.occupied(t))
					validPoints.insertMem(t);
			}

        if(validPoints.itemCount())
			for(Idx tIdx = 0;tIdx < validPoints.itemCount(); ++tIdx)
			{
				Change temp = tChange;
				temp.add(tPos, validPoints[tIdx]);
				mChanges.insertMem(temp);
			}
    }
    CatchError
}
closePlatypusSpace
