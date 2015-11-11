/*!
	@file globals/init.hpp
	@brief The implementation file for Init class.
	@details This is the implementation file for Init class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#include "pspl/globals/init.hpp"


openPlatypusSpace



/*!
	The constructor.
*/
Init::Init() : mChange()
{
	WatchError
	//	nothing to be done.
	CatchError
}



/*!
	The destructor.
*/
Init::~Init()
{
	WatchError
	//	nothing to be done.
	CatchError
}


/*!
	The duplicator
*/
Init::Init(Init const & that)
{
	WatchError
	Alert(&that, eUndefDuplicator);
	CatchError
}


/*!
	The assigner.
*/
Init const & Init::operator= (Init const & that)
{
	WatchError
	Alert(&that, eUndefAssigner);
    return *this;
	CatchError
}



closePlatypusSpace
