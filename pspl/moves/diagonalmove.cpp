/*!
	@file moves/diagonalmove.cpp
	@brief The implementation file for DiagonalMove class.
	@details This is the implementation file for DiagonalMove class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#include "pspl/moves/diagonalmove.hpp"
#include "pspl/globals/conf.hpp"
#include "pspl/globals/change.hpp"
#include "pspl/globals/lattice.hpp"



openPlatypusSpace



/*!
	The costructor.
*/
DiagonalMove::DiagonalMove()
{
    WatchError
    // do nothing
    CatchError
}


/*!
	The Destructor
*/
DiagonalMove::~DiagonalMove()
{
    WatchError
    // do nothing
    CatchError
}


/*!
	The Dupllicator
*/
DiagonalMove::DiagonalMove(DiagonalMove const & that) :
	Move(that)
{
    WatchError
    Alert(&that, eUndefDuplicator);
    CatchError
}


/*!
	The assigner.
*/

DiagonalMove const & DiagonalMove::operator= (DiagonalMove const & that)
{
    WatchError
	Alert(&that, eUndefAssigner);
    return *this;
    CatchError
}


/*!
	The compute method
save all possible move results into mChanges

i-------i+1
|       .
|       .
i-1. . .tPlace

*/
void DiagonalMove::compute(Conf & theConf, Pos const thePos)
{
    WatchError
    Warn(thePos >= Protein::p().Length(), eInvalidPos);
    mChanges.clear(); Change tChange;
	if (!thePos || thePos == Protein::p().Span())
	{
		Point const & sPoint = theConf.Coord(!thePos ? thePos + 1 : thePos -1);
		for(Dir tDir = 0; tDir < Lattice::l().NeighCount(); ++tDir)
		{
			Point tPoint = sPoint + Lattice::l().DirVec(tDir);
			if(!theConf.occupied(tPoint))
			{
				tChange.add(thePos, tPoint);
				mChanges.insertMem(tChange);
				tChange.reset();
			}
		}
	}
	else
	{
		Point const & sPoint = theConf.Coord(thePos -1);
		for(Dir tDir = 0; tDir < Lattice::l().NeighCount(); ++tDir)
		{
			Point tPoint = sPoint + Lattice::l().DirVec(tDir);
			if(theConf.occupied(tPoint)) continue; // not feasible move
			if(!Lattice::l().areNeighbors(theConf.Coord(thePos + 1), tPoint))
					continue;
            tChange.add(thePos, tPoint);
			mChanges.insertMem(tChange);
			tChange.reset();
		}
	}
    CatchError
}
closePlatypusSpace
