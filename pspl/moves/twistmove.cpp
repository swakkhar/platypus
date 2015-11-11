/*!
	@file moves/TwistMove.cpp
	@brief The implementation file for TwistMove class.
	@details This is the implementation file for TwistMove class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/

#include "pspl/moves/twistmove.hpp"
#include "pspl/globals/conf.hpp"
#include "pspl/globals/change.hpp"
#include "pspl/globals/lattice.hpp"


openPlatypusSpace


/*!
The costructor
*/
TwistMove::TwistMove()
{
    WatchError

    // do nothing
    CatchError
}


/*!

The Destructor

*/
TwistMove::~TwistMove()
{
    WatchError

    // do nothing

    CatchError
}


/*!

The Dupllicator
*/

TwistMove::TwistMove(TwistMove const & that):Move(that)
{
    WatchError

    Alert(&that, eUndefDuplicator);

    CatchError
}


/*!

Assignment Operator
*/

TwistMove const & TwistMove::operator= (TwistMove const & that)
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
void TwistMove::compute(Conf & theConf, Pos const thePos)
{
    //WatchError
    mChanges.clear();
    if(thePos > Protein::p().Span() - 3) return;


    for(Pos tPos=thePos+2; tPos <=Protein::p().Span()-1;++tPos)
    {
        Point a = theConf.Coord(thePos);
        Point b = theConf.Coord(thePos+1);
        Point c = theConf.Coord(tPos);
        Point d = theConf.Coord(tPos+1);
        if(!Lattice::l().areNeighbors(a,c)) continue;
        if(!Lattice::l().areNeighbors(b,d)) continue;
        // else twist possible
        Change tChange;
        for(Cmp tCmp = thePos+1; tCmp <= tPos; ++tCmp)
        {
            tChange.add(tCmp,theConf.Coord(tPos-tCmp+thePos+1));
        }
        mChanges.insertMem(tChange);

    }
    //CatchError
}
closePlatypusSpace
