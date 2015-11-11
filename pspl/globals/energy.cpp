/*!
	@file globals/energy.cpp
	@brief The implementation file for Energy class.
	@details This is the implementation file for Energy class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/


#include "pspl/globals/energy.hpp"


openPlatypusSpace

/*!
	Global array of energy models.
*/
Energy const * Energy::mEnergy(Null);



/*!
	The destructor.
*/
Energy::~Energy()
{
	WatchError
	//nothing to be done.
	CatchError
}


closePlatypusSpace
