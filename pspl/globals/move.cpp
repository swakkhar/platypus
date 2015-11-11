/*!
	@file globals/move.hpp
	@brief The implementation file for Move class.
	@details This is the implementation file for Move class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#include "pspl/globals/move.hpp"


openPlatypusSpace


/*!
	The constructor.
*/
Move::Move() : mChanges()
{
	WatchError
	//	nothing to be done.
	CatchError
}



/*!
	The destructor.
*/
Move::~Move()
{
	WatchError
	//	nothing to be done.
	CatchError
}


/*!
	The duplicator
*/
Move::Move(Move const & that)
{
	WatchError
	Alert(&that, eUndefDuplicator);
	CatchError
}


/*!
	The assigner.
*/
Move const & Move::operator= (Move const & that)
{
	WatchError
	Alert(&that, eUndefAssigner);
    return *this;
	CatchError
}



closePlatypusSpace
