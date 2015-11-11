/*!
	@file moves/CrankShaftMove.cpp
	@brief The implementation file for CrankShaftMove class.
	@details This is the implementation file for CrankShaftMove class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/

#include "pspl/moves/crankshaftmove.hpp"
#include "pspl/globals/conf.hpp"
#include "pspl/globals/change.hpp"
#include "pspl/globals/lattice.hpp"


openPlatypusSpace


/*!
The costructor
*/
CrankShaftMove::CrankShaftMove()
{
    WatchError

    // do nothing
    CatchError
}


/*!

The Destructor

*/
CrankShaftMove::~CrankShaftMove()
{
    WatchError

    // do nothing

    CatchError
}


/*!

The Dupllicator
*/

CrankShaftMove::CrankShaftMove(CrankShaftMove const & that):Move(that)
{
    WatchError

    Alert(&that, eUndefDuplicator);

    CatchError
}


/*!

Assignment Operator
*/

CrankShaftMove const & CrankShaftMove::operator= (CrankShaftMove const & that)
{
    WatchError

    Move::operator=(that);
    return *this;


    CatchError
}


/*!

The compute method
save all possible move results into mChanges


i-------i+1
|       |
|       |
i-1     i+2
-       +
-       +
-       +
l-------c


*/
void CrankShaftMove::compute(Conf & theConf, Pos const thePos)
{
    WatchError
    mChanges.clear();
    if(thePos >= Protein::p().Span()) return;
    Change tChange;

	if (!thePos || thePos == Protein::p().Span() - 1)
	{
		Pos const tPosA = !thePos ? thePos + 2 : thePos - 1;
		Pos const tPosB = !thePos ? thePos + 1 : thePos;
		Pos const tPosC = !thePos ? thePos : thePos + 1;
		Point const & a = theConf.Coord(tPosA);
		for(Dir tDir1 = 0; tDir1 < Lattice::l().NeighCount(); ++tDir1)
		{
            Point c = a + Lattice::l().DirVec(tDir1);
            if(theConf.occupied(c)) continue;
			for(Dir tDir2 = 0; tDir2 < Lattice::l().NeighCount(); ++tDir2)
			{
				Point d = c + Lattice::l().DirVec(tDir2);
				if(theConf.occupied(d)) continue;
				tChange.add(tPosB, c);
				tChange.add(tPosC, d);
				mChanges.insertMem(tChange);
				tChange.reset();
			}
		}
	}
	else
	{
		Pos const tPosA = thePos - 1;
		Pos const tPosB = thePos;
		Pos const tPosC = thePos + 1;
		Pos const tPosD = thePos + 2;
		Point const & a = theConf.Coord(tPosA);
		for(Dir tDir1 = 0; tDir1 < Lattice::l().NeighCount(); ++tDir1)
		{
            Point c = a + Lattice::l().DirVec(tDir1);
            if(theConf.occupied(c)) continue;
			for(Dir tDir2 = 0; tDir2 < Lattice::l().NeighCount(); ++tDir2)
			{
				Point d = c + Lattice::l().DirVec(tDir2);
				if(theConf.occupied(d)) continue;
				if (!Lattice::l().areNeighbors(d, theConf.Coord(tPosD)))
					continue;
				tChange.add(tPosB, c);
				tChange.add(tPosC, d);
				mChanges.insertMem(tChange);
				tChange.reset();
			}
		}
	}
	CatchError
}
closePlatypusSpace
