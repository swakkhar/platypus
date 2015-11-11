/*!
	@file globals/space.hpp
	@brief The prototype file for Space class.
	@details This is the prototype file for Space class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/


#ifndef SpaceHppIncluded
#define SpaceHppIncluded



#include "pspl/globals/data.hpp"



openPlatypusSpace


/*!
	@class Space
	@brief A class to represent the coordinate space.
	@details This class represents the coordinate space.
*/
class Space
{
	public:

		static Dim Dimen();						//!< Return the dimension of the space.
		static void Dimen(Dim const theDimen);	//!< Set the dimension of the space.

	private:
		static Dim mDimen;						//!< The dimension of the space.
};



/*!
	Return the dimension of the coordinate space.
*/
inline Dim Space::Dimen()
{
	WatchError
	return mDimen;
	CatchError
}



/*!
	Set the dimension of the coordinate space.
*/
inline void Space::Dimen(Dim const theDimen)
{
	WatchError
	Alert(mDimen, eNonEmptyDimension);
	Alert(!theDimen, eEmptyDimension);

	mDimen = theDimen;
	CatchError
}



closePlatypusSpace



#endif //SpaceHppIncluded
