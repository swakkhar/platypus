/*!
	@file globals/protein.cpp
	@brief The implementation file for Protein class.
	@details This is the implementation file for Protein class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/

#include "pspl/globals/protein.hpp"



openPlatypusSpace


/*!
	Static variable.
*/
Protein Protein::mProtein;



/*!
	The initialiser.
*/
Protein::Protein(Dim const AlphaSize, Acd const * theAcids)
{
	WatchError
	mLength = 0;
	while(*theAcids)
	{
		mAcids[mLength] = *theAcids;
		mTypes[mLength] = Energy::e().Type(AlphaSize, *theAcids);
		++mLength; ++theAcids;
	}
	Alert(!mLength, eEmptyDimension);
	mSpan = mLength - 1;
	CatchError
}



/*!
	The duplicator.
*/
Protein::Protein(Protein const & that) : mLength(that.mLength),
		mSpan(that.mSpan)
{
	WatchError
	for(Pos tPos = 0; tPos < mLength; ++tPos)
	{
		mAcids[tPos] = that.mAcids[tPos];
		mTypes[tPos] = that.mTypes[tPos];
	}
	CatchError
}



/*!
	The assigner.
*/
Protein const & Protein::operator = (Protein const & that)
{
	WatchError
	Warn(mLength, eNonEmptyProtein);

	if (this == &that)
		return *this;
	mLength = that.mLength;
	mSpan = that.mSpan;
	for(Pos tPos = 0; tPos < mLength; ++tPos)
	{
		mAcids[tPos] = that.mAcids[tPos];
		mTypes[tPos] = that.mTypes[tPos];
	}
	return *this;
	CatchError
}




closePlatypusSpace
