/*!
	@file energys/hpenergy.cpp
	@brief The implementation file for HpEnergy class.
	@details This is the implementation file for HpEnergy class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/


#include "pspl/energys/hpenergy.hpp"


openPlatypusSpace


Typ HpEnergy::Type(Dim const AlphaSize, Acd const theAcid) const
{
	WatchError

	switch(AlphaSize)
	{
		case 2:
			switch(theAcid)
			{
				case 'P' : return 0;
				case 'H' : return 1;
				default : Throw(eUnknownAcid);
			}
		default:
			Throw(eUnknownAlpha);
	}
	CatchError
}

Eng HpEnergy::Level(Typ const Typ1, Typ const Typ2) const
{
	WatchError
	Warn(Typ1 >= mTypeCount || Typ2 >= mTypeCount, eInvalidType);
	return -(Typ1 & Typ2);
	CatchError
}



closePlatypusSpace
