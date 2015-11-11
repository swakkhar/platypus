/*!
	@file moves/DoublePullMove.cpp
	@brief The implementation file for DoublePullMove class.
	@details This is the implementation file for DoublePullMove class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/

#include "pspl/moves/doublepullmove.hpp"
#include "pspl/globals/conf.hpp"
#include "pspl/globals/change.hpp"
#include "pspl/globals/lattice.hpp"


openPlatypusSpace


/*!
The costructor
*/
DoublePullMove::DoublePullMove()
{
    WatchError

    // do nothing
    CatchError
}


/*!

The Destructor

*/
DoublePullMove::~DoublePullMove()
{
    WatchError

    // do nothing

    CatchError
}


/*!

The Dupllicator
*/

DoublePullMove::DoublePullMove(DoublePullMove const & that):Move(that)
{
    WatchError

    Alert(&that, eUndefDuplicator);

    CatchError
}


/*!

Assignment Operator
*/

DoublePullMove const & DoublePullMove::operator= (DoublePullMove const & that)
{
    WatchError

    Move::operator=(that);
    return *this;


    CatchError
}


/*!

The compute method
save all possible move results into mChanges

i-1
|
|
|
i.......c
|       .
|       .
|       .
i+1.....l

i-1 --> c
i -- > l
*/
void DoublePullMove::compute(Conf & theConf, Pos const thePos)
{
    WatchError
    mChanges.clear(); // reset previous info
    if (!thePos || thePos >= Protein::p().Span()) return;
	Change tChange;
    for(Int pullDir = -1; pullDir <= 1; pullDir+=2)
    {
		for(Dir tDir1=0; tDir1 < Lattice::l().NeighCount(); ++tDir1)
		{
			for(Dir tDir2=0; tDir2 < Lattice::l().NeighCount(); ++tDir2)
			{
				tChange.reset();
				Point l = theConf.Coord(thePos + pullDir) + Lattice::l().DirVec(tDir1);
				if(theConf.occupied(l))continue;
				Point c = theConf.Coord(thePos) + Lattice::l().DirVec(tDir2);
				if(theConf.occupied(c))continue;
				if(!Lattice::l().areNeighbors(l,c)) continue;

				Pos tPos = thePos;

				while(tPos <= Protein::p().Span()&& l!=theConf.Coord(tPos))
				{
					tChange.add(tPos,l);
					if((Idx)(tPos-pullDir) <=Protein::p().Span() && Lattice::l().areNeighbors(l,theConf.Coord(tPos-pullDir))) // keeps the move local
						break;
					l=c;
					c=theConf.Coord(tPos);
					tPos-=pullDir;
				}
				mChanges.insertMem(tChange);
			}
		}
    }
    CatchError
}
closePlatypusSpace
