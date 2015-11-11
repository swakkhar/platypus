/*!
	@file moves/TriplePullMove.cpp
	@brief The implementation file for TriplePullMove class.
	@details This is the implementation file for TriplePullMove class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/

#include "pspl/moves/triplepullmove.hpp"
#include "pspl/globals/conf.hpp"
#include "pspl/globals/change.hpp"
#include "pspl/globals/lattice.hpp"


openPlatypusSpace


/*!
The costructor
*/
TriplePullMove::TriplePullMove()
{
    WatchError

    // do nothing
    CatchError
}


/*!

The Destructor

*/
TriplePullMove::~TriplePullMove()
{
    WatchError

    // do nothing

    CatchError
}


/*!

The Dupllicator
*/

TriplePullMove::TriplePullMove(TriplePullMove const & that):Move(that)
{
    WatchError

    Alert(&that, eUndefDuplicator);

    CatchError
}


/*!

Assignment Operator
*/

TriplePullMove const & TriplePullMove::operator= (TriplePullMove const & that)
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
void TriplePullMove::compute(Conf & theConf, Pos const thePos)
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

				for(Dir tDir3=0; tDir3 < Lattice::l().NeighCount(); ++tDir3)
				{
				    tChange.reset();
                    Point l = theConf.Coord(thePos + pullDir) + Lattice::l().DirVec(tDir1);
                    if(theConf.occupied(l))continue;

                    Point c = theConf.Coord(thePos) + Lattice::l().DirVec(tDir2);
                    if(theConf.occupied(c))continue;

                    if(l==c) continue;


                    Point n = l + Lattice::l().DirVec(tDir3);
                    if(theConf.occupied(n))continue;
                    if(!Lattice::l().areNeighbors(n,c)) continue;
                    Pos tPos = thePos;

                    while (tPos <= Protein::p().Span()&& l!=theConf.Coord(tPos))
                    {

                        tChange.add(tPos,l);
                        if((Idx)(tPos-pullDir) <=Protein::p().Span() &&
                           Lattice::l().areNeighbors(l,theConf.Coord(tPos-pullDir))) // keeps the move local
                            break;
                        l=n;
                        n=c;
                        c=theConf.Coord(tPos);
                        tPos-=pullDir;

                    }
                    mChanges.insertMem(tChange);
				/*if(!Lattice::l().areNeighbors(l,c)) continue;

				Pos tPos = thePos;

				while(tPos <= Protein::p().Span()&& l!=theConf.Coord(tPos))
				{
					tChange.add(tPos,l);
					if((tPos-pullDir) <=Protein::p().Span() && Lattice::l().areNeighbors(l,theConf.Coord(tPos-pullDir))) // keeps the move local
						break;
					l=c;
					c=theConf.Coord(tPos);
					tPos-=pullDir;
				}
				mChanges.annex(tChange);*/
				}
			}
		}
    }


    CatchError
}
closePlatypusSpace
