/*!
	@file moves/SinglePullMove.cpp
	@brief The implementation file for SinglePullMove class.
	@details This is the implementation file for SinglePullMove class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/

#include "pspl/moves/singlepullmove.hpp"
#include "pspl/globals/conf.hpp"
#include "pspl/globals/change.hpp"
#include "pspl/globals/lattice.hpp"


openPlatypusSpace


/*!
The costructor
*/
SinglePullMove::SinglePullMove()
{
    WatchError

    // do nothing
    CatchError
}


/*!

The Destructor

*/
SinglePullMove::~SinglePullMove()
{
    WatchError

    // do nothing

    CatchError
}


/*!

The Dupllicator
*/

SinglePullMove::SinglePullMove(SinglePullMove const & that):Move(that)
{
    WatchError

    Alert(&that, eUndefDuplicator);

    CatchError
}


/*!

Assignment Operator
*/

SinglePullMove const & SinglePullMove::operator= (SinglePullMove const & that)
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
void SinglePullMove::compute(Conf & theConf, Pos const thePos)
{
    WatchError

    mChanges.clear(); // reset previous info

    if(!thePos || thePos >= Protein::p().Span()) return;
    Change tChange;

	for(Int pullDir = -1; pullDir <= 1; pullDir += 2)
 	{
		Point f = theConf.Coord(thePos + pullDir);
		for(Dir tDir = 0; tDir < Lattice::l().NeighCount(); ++tDir)
		{
			Point t = f + Lattice::l().DirVec(tDir);
			if(theConf.occupied(t)) continue;
			if (!Lattice::l().areNeighbors(t, theConf.Coord(thePos)))
				continue;

			Pos tPos = thePos;
			while(tPos < Protein::p().Length())
			{
				tChange.add(tPos,t);
				Pos sPos = tPos - pullDir;
				if(sPos <= Protein::p().Span() &&
						Lattice::l().areNeighbors(t, theConf.Coord(sPos))) // keeps the move local
					break;
				t = theConf.Coord(tPos);
				tPos = sPos;
			}
			mChanges.insertMem(tChange);
			tChange.reset();
		}
	}
    CatchError
}
closePlatypusSpace
