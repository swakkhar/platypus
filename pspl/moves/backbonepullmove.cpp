/*!
	@file moves/BackBonePullMove.cpp
	@brief The implementation file for BackBonePullMove class.
	@details This is the implementation file for BackBonePullMove class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/

#include "pspl/moves/backbonepullmove.hpp"
#include "pspl/globals/conf.hpp"
#include "pspl/globals/change.hpp"
#include "pspl/globals/lattice.hpp"


openPlatypusSpace


/*!
The costructor
*/
BackBonePullMove::BackBonePullMove()
{
    WatchError

    // do nothing
    CatchError
}


/*!

The Destructor

*/
BackBonePullMove::~BackBonePullMove()
{
    WatchError

    // do nothing

    CatchError
}


/*!

The Dupllicator
*/

BackBonePullMove::BackBonePullMove(BackBonePullMove const & that):Move(that)
{
    WatchError

    Alert(&that, eUndefDuplicator);

    CatchError
}


/*!

Assignment Operator
*/

BackBonePullMove const & BackBonePullMove::operator= (BackBonePullMove const & that)
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


void BackBonePullMove::compute(Conf & theConf, Pos const thePos)
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





                for(Dir tDir3=0;tDir3<Lattice::l().NeighCount(); ++tDir3)
                {
                    for(Dir tDir4=0;tDir4<Lattice::l().NeighCount(); ++tDir4)
                    {
                        Point l = theConf.Coord(thePos + pullDir) + Lattice::l().DirVec(tDir1);
                        if(theConf.occupied(l))continue;
                        Point c = theConf.Coord(thePos) + Lattice::l().DirVec(tDir2);
                        if(theConf.occupied(c))continue;
                        if(!Lattice::l().areNeighbors(l,c)) continue;

                        // at this situation we also need two more positions
                        // lS and cS that are free and not equal to c and l respectively
                        // then we can start the pull, otherwise not

                        // check for side chains

                        tChange.reset();
                        Point lS=l+Lattice::l().DirVec(tDir3);
                        if(theConf.occupied(lS)||lS==c)continue;
                        Point cS=c+Lattice::l().DirVec(tDir4);
                        if(theConf.occupied(cS)||cS==l)continue;

                        // lS and CS cant be the same position
                        if(lS==cS) continue;

                        // else move possible
                        Pos tPos = thePos;

                        while(tPos>=0 && tPos <= (Protein::p().Length()/2-1) && l!=theConf.Coord(tPos))
                        {
                            tChange.add(tPos,l);
                            tChange.add(tPos+Protein::p().Length()/2,lS);
                            if((Idx)(tPos-pullDir) <=(Protein::p().Length()/2-1) &&
                                (tPos-pullDir) >=0 &&
                                Lattice::l().areNeighbors(l,theConf.Coord(tPos-pullDir))) // keeps the move local
                                break;
                            l=c;
                            lS=cS;
                            c=theConf.Coord(tPos);
                            cS=theConf.Coord(tPos+Protein::p().Length()/2);
                            tPos-=pullDir;

                        }
                        mChanges.insertMem(tChange);

                    }

                }




			}
		}
    }
    CatchError
}





/*!********************************************
COMMENTED ONE IS A SINGLEPULL VERSION
THAT DOES NOT WORK FOR CUBIC LATTICES
I NEED TO CREATE A DOUBLE PULL VERSION
********************************************


*/


/*void BackBonePullMove::compute(Conf & theConf, Pos const thePos)
{
    WatchError
    mChanges.clear();

    if(!thePos||thePos>=Protein::p().Length()/2-1) return;
    Change tChange;

    for(Int pullDir = -1; pullDir <= 1; pullDir += 2)
    {
        Point fB = theConf.Coord(thePos + pullDir);
        Point fS = theConf.Coord(thePos + pullDir+Protein::p().Length()/2);

        for(Dir tDir = 0; tDir < Lattice::l().NeighCount(); ++tDir)
        {
            Point tB = fB + Lattice::l().DirVec(tDir);
            if(theConf.occupied(tB)) continue;
            if (!Lattice::l().areNeighbors(tB, theConf.Coord(thePos)))
				continue;
            // now we need to get another neighbor for the side chain
            for(Dir tDir2 = 0; tDir2 < Lattice::l().NeighCount(); ++tDir2)
            {
                Point tS=tB+Lattice::l().DirVec(tDir);
                if(theConf.occupied(tB)) continue;

                // else pull the loop
                Pos tPos = thePos;
                while(tPos<Protein::p().Length()/2&&tPos>=0)
                {
                    tChange.add(tPos,tB);
                    tChange.add(tPos+Protein::p().Length()/2,tS);
                    Pos sPos = tPos - pullDir;

                    if(sPos>=0 && sPos<=Protein::p().Length()/2-1 &&
                                    Lattice::l().areNeighbors(tB,theConf.Coord(sPos)))
                                break;
                    tB=theConf.Coord(tPos);
                    tS=theConf.Coord(tPos+Protein::p().Length()/2);
                    tPos=sPos;

                }
                mChanges.insertMem(tChange);
                tChange.reset();


            }


        }
    }
    CatchError
}*/
closePlatypusSpace
