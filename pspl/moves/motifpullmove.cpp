/*!
	@file moves/MotifPullMove.cpp
	@brief The implementation file for MotifPullMove class.
	@details This is the implementation file for MotifPullMove class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/

#include "pspl/moves/motifpullmove.hpp"
#include "pspl/globals/conf.hpp"
#include "pspl/globals/change.hpp"
#include "pspl/globals/lattice.hpp"


openPlatypusSpace


/*!
The costructor
*/
MotifPullMove::MotifPullMove()
{
    WatchError

    // do nothing
    CatchError
}


/*!

The Destructor

*/
MotifPullMove::~MotifPullMove()
{
    WatchError

    // do nothing

    CatchError
}


/*!

The Dupllicator
*/

MotifPullMove::MotifPullMove(MotifPullMove const & that):Move(that)
{
    WatchError

    Alert(&that, eUndefDuplicator);

    CatchError
}


/*!

Assignment Operator
*/

MotifPullMove const & MotifPullMove::operator= (MotifPullMove const & that)
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
void MotifPullMove::compute(Conf & theConf, Pos const thePos)
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
			B motifBreak=false;
			while(tPos < Protein::p().Length())
			{
			    if(!theConf.isPosMotif(tPos))
			    {
                    tChange.add(tPos,t);
                    Pos sPos = tPos - pullDir;
                    if(sPos <= Protein::p().Span() &&
                            Lattice::l().areNeighbors(t, theConf.Coord(sPos))) // keeps the move local
                        break;
                    t = theConf.Coord(tPos);
                    tPos = sPos;
			    }
			    else
                {

                    Dir tDir1=0,tDir2=1;
                    if(theConf.getMotifAtPos(tPos)==1) // beta sheet
                    {
                        // check if possible to move
                        Pos tPos1=tPos-pullDir;
                        Pos tPos2=tPos-2*pullDir;
                        Point p0=theConf.Coord(tPos);
                        Point p1=theConf.Coord(tPos1);
                        Point p2=theConf.Coord(tPos2);

                        tDir1=Lattice::l().Direction(p0-t);
                        tDir2=Lattice::l().Direction(p2-p1);
                    }
                    else if(theConf.getMotifAtPos(tPos)==2) // alpha helix
                    {
                        // check if possible to move
                        Pos tPos3=tPos-3*pullDir;
                        Pos tPos4=tPos-4*pullDir;
                        Point p0=theConf.Coord(tPos);
                        Point p3=theConf.Coord(tPos3);
                        Point p4=theConf.Coord(tPos4);

                        tDir1=Lattice::l().Direction(p0-t);
                        tDir2=Lattice::l().Direction(p4-p3);
                    }
                    if(tDir1==tDir2)
                    {
                        // possible
                        // for the rest of the rigid body / motif perform a pull
                        Dir code=theConf.getMotifAtPos(tPos);
                        while(tPos<Protein::p().Length()&&theConf.getMotifAtPos(tPos)==code)
                        {
                            // to the last of the subsequece of this current motif
                            tChange.add(tPos,t);
                            Pos sPos = tPos - pullDir;
                            t = theConf.Coord(tPos);
                            tPos = sPos;

                        }

                    }
                    else
                    {
                        motifBreak=true;
                        break;
                    }

			    }

			}
			if(motifBreak==true) // if the motif pull was not possible then discard current changes
			{
			    tChange.reset();
			    continue;
			}
			/*for(int i=0;i<tChange.size();++i)
			{
                cout << "("<<tChange.Destination(i).Comp(0) << " "<< tChange.Destination(i).Comp(1) << " "
                << tChange.Destination(i).Comp(2) << " -" <<tChange.Position(i)<<") ";
			}
			cout << endl;*/
			mChanges.insertMem(tChange);
			tChange.reset();
		}
	}
    CatchError
}
closePlatypusSpace
