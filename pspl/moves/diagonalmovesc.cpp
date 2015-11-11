/*!
	@file moves/diagonalmove.cpp
	@brief The implementation file for DiagonalMoveSC class.
	@details This is the implementation file for DiagonalMoveSC class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#include "pspl/moves/diagonalmovesc.hpp"
#include "pspl/globals/conf.hpp"
#include "pspl/globals/change.hpp"
#include "pspl/globals/lattice.hpp"



openPlatypusSpace



/*!
	The costructor.
*/
DiagonalMoveSC::DiagonalMoveSC()
{
    WatchError
    // do nothing
    CatchError
}


/*!
	The Destructor
*/
DiagonalMoveSC::~DiagonalMoveSC()
{
    WatchError
    // do nothing
    CatchError
}


/*!
	The Dupllicator
*/
DiagonalMoveSC::DiagonalMoveSC(DiagonalMoveSC const & that) :
	Move(that)
{
    WatchError
    Alert(&that, eUndefDuplicator);
    CatchError
}


/*!
	The assigner.
*/

DiagonalMoveSC const & DiagonalMoveSC::operator= (DiagonalMoveSC const & that)
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
void DiagonalMoveSC::compute(Conf & theConf, Pos const thePos)
{

    WatchError
    Warn(thePos >= Protein::p().Length(), eInvalidPos);
    mChanges.clear(); Change tChange;

    Z chainLength = Protein::p().Length();

    if(theConf.isPosSideChain(thePos))
    {
        // side chain move
        Pos bPos=thePos-chainLength/2;

        Point bPoint=theConf.Coord(bPos);
        for(Dir tDir = 0; tDir < Lattice::l().NeighCount(); ++tDir)
        {
            Point tPoint=bPoint+Lattice::l().DirVec(tDir);
            if(!theConf.occupied(tPoint))
            {
                tChange.add(thePos,tPoint);
                mChanges.insertMem(tChange);
                tChange.reset();
            }
        }
    }
    else
    {
            // backbone move
        Point curBPoint=theConf.Coord(thePos);
        Point curSPoint=theConf.Coord(thePos+chainLength/2);

        // now this position can be 0 or the last backbone point, in both cases
        // we need only one check, else we also check the neighbors for i+1
        Pos pivotPos=0;
        if(thePos==0)pivotPos=thePos+1;
        else pivotPos=thePos-1;

        // now need to get prevoius backbone point

        Point prevBPoint=theConf.Coord(pivotPos);

        // now we need to seek a free point arround this prevBPoint
        for(Dir tDir = 0; tDir < Lattice::l().NeighCount(); ++tDir)
        {
            // temporary Backbone point
            Point tempBPoint=prevBPoint+Lattice::l().DirVec(tDir);

            // this has to be free or it can be the current Side chain point
            // since that side-chain point is going to be removed

            if(!theConf.occupied(tempBPoint)||tempBPoint==curSPoint)
            {

                // i+1 check

                if(thePos!=0||thePos!=(chainLength/2-1))
                {
                    // check here
                    Point nextBPoint=theConf.Coord(thePos+1);
                    if(!Lattice::l().areNeighbors(nextBPoint, tempBPoint))
                        continue;
                }


                // backbone alright
                // we need to look for a sidechain arround this point
                for(Dir tDir2 = 0; tDir2 < Lattice::l().NeighCount(); ++tDir2)
                {
                    Point tempSPoint=tempBPoint+Lattice::l().DirVec(tDir);
                    // this is required to be free or current backbone point
                    // since that has to be removed
                    if(!theConf.occupied(tempBPoint)||tempBPoint==curSPoint)
                    {
                        // now we have got two points
                        // tempBPoint && tempSPoint
                        // they are neighbors to each other and potentially free

                        // we can add them to the changes

                        tChange.add(thePos,tempBPoint);
                        tChange.add(thePos+chainLength/2,tempSPoint);
                        mChanges.insertMem(tChange);
                        tChange.reset();
                    }
                }

            }


        }



    }
    CatchError
}
closePlatypusSpace
