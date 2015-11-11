/*!
	@file inits/BlockInit.cpp
	@brief The implementation file for BlockInit class.
	@details This is the implementation file for BlockInit class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#include "pspl/inits/blockinit.hpp"
#include "pspl/globals/conf.hpp"
#include "pspl/globals/change.hpp"
#include "pspl/globals/lattice.hpp"



openPlatypusSpace



/*!
	The costructor.
*/
BlockInit::BlockInit(block1<Point,kmm> &arr) : mPoints(arr)
{
    WatchError
    // do nothing
    CatchError
}


/*!
	The Destructor
*/
BlockInit::~BlockInit()
{
    WatchError
    // do nothing
    CatchError
}


/*!
	The Dupllicator
*/
BlockInit::BlockInit(BlockInit const & that) :
	Init(that), mPoints(that.mPoints)
{
    WatchError
    Alert(&that, eUndefDuplicator);
    CatchError
}


/*!
	The assigner.
*/

BlockInit const & BlockInit::operator= (BlockInit const & that)
{
    WatchError
	Alert(&that, eUndefAssigner);
    return *this;
    CatchError
}


/*!
	Compute and store a possible initialisation change.
*/
void BlockInit::compute(Conf & theConf)
{
    WatchError
    for(Pos tPos=0;tPos<Protein::p().Length();++tPos)
	mChange.add(tPos, mPoints[tPos]);
	CatchError
}



closePlatypusSpace
