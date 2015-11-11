/*!
	@file globals/lattice.hpp
	@brief The implementation file for Lattice class.
	@details This is the implementation file for Lattice class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#include "pspl/globals/lattice.hpp"


openPlatypusSpace


/*!
	Initialiser global lattice.
*/
Lattice const * Lattice::mLattice(Null);



/*!
	The initialiser.
*/
Lattice::Lattice(Mid const theModelId, Dim const theNeighCount, Dim const theSqrNeighDist) :
		mModelId(theModelId), mNeighCount(theNeighCount), mSqrNeighDist(theSqrNeighDist)
{
	WatchError
	Warn(!theNeighCount, eEmptyDimension);
	Warn(!theSqrNeighDist, eEmptyDimension);
	CatchError
}



/*!
	The duplicator.
*/
Lattice::Lattice(Lattice const & that) :
		mModelId(that.mModelId), mNeighCount(that.mNeighCount), mSqrNeighDist(that.mSqrNeighDist)
{
	WatchError
	//	nothing to be done.
	CatchError
}



/*!
	The assigner.
*/
Lattice const & Lattice::operator = (Lattice const & that)
{
	WatchError
	Warn(mModelId != that.mModelId, eLatticeMismatch);
	if (this != &that)
	{
		mNeighCount = that.mNeighCount;
		mSqrNeighDist = that.mSqrNeighDist;
	}
	return * this;
	CatchError
}




/*!
	The destructor.
*/
Lattice::~Lattice()
{
	WatchError
	//	nothing to be done.
	CatchError
}



closePlatypusSpace
