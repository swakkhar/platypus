/*!
	@file moves/TiltMove.cpp
	@brief The implementation file for TiltMove class.
	@details This is the implementation file for TiltMove class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/

#include "pspl/moves/tiltmove.hpp"
#include "pspl/globals/conf.hpp"
#include "pspl/globals/change.hpp"
#include "pspl/globals/lattice.hpp"


openPlatypusSpace


/*!
The costructor
*/
TiltMove::TiltMove()
{
    WatchError

    // do nothing
    CatchError
}


/*!

The Destructor

*/
TiltMove::~TiltMove()
{
    WatchError

    // do nothing

    CatchError
}


/*!

The Dupllicator
*/

TiltMove::TiltMove(TiltMove const & that):Move(that)
{
    WatchError

    Alert(&that, eUndefDuplicator);

    CatchError
}


/*!

Assignment Operator
*/

TiltMove const & TiltMove::operator= (TiltMove const & that)
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
void TiltMove::compute(Conf & theConf, Pos const thePos)
{
    WatchError
    mChanges.clear();
    if(!thePos || thePos >= Protein::p().Span()) return;

    Point i = theConf.Coord(thePos);
    Point j = theConf.Coord(thePos+1);



    Change tChange;

    for(Dir tDir1 = 0; tDir1 < Lattice::l().NeighCount(); ++tDir1)
    {
        for(Dir tDir2 = 0; tDir2 < Lattice::l().NeighCount(); ++tDir2)
        {
            Point a = i + Lattice::l().DirVec(tDir1);
            if(theConf.occupied(a)) continue;

            Point b = j + Lattice::l().DirVec(tDir2);
            if(theConf.occupied(b)) continue;
            if(!Lattice::l().areNeighbors(a,b)) continue;

            // now start pull
            Pos tPos = thePos;
            while(tPos <=Protein::p().Span())
            {
                tChange.add(tPos, a);
                if(tPos&&Lattice::l().areNeighbors(a,theConf.Coord(tPos-1)))
                    break;
                a = theConf.Coord(tPos);
                --tPos;
            }
            tPos = thePos+1;
            while(tPos <=Protein::p().Span())
            {
                tChange.add(tPos, b);
                if(tPos !=Protein::p().Span()&&
                   Lattice::l().areNeighbors(b,theConf.Coord(tPos+1)))
                    break;
                b = theConf.Coord(tPos);
                ++tPos;
            }

            mChanges.insertMem(tChange);
            tChange.reset();
        }
    }
	CatchError
}
closePlatypusSpace
