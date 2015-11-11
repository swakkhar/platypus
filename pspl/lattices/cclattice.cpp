/*!
	@file lattices/cclattice.hpp
	@brief The prototype file for CCLattice class.
	@details This is the prototype file for CCLattice class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/


#include "pspl/lattices/cclattice.hpp"


openPlatypusSpace


/*!
	@struct CCData
	@brief A structure to represent vital data for cc-lattice.
	@details The structure represents vital data for cc-lattice.
*/
struct CCData
{
			char name[3]; 	//!< name of direction
			Spc vec[3];     //!< absolute direction vector
			Spc mat[3][3];  //!< matrix for applying relative direction
			Spc inv[3][3]; //!< inverse of matrix (we are lazy ;-))
};



/*!
	@var ccData
	@brief Vital data for cc-lattice.
	@details The vital data for cc-lattice.
*/
CCData const ccData[] = {
			{"F",{ 1, 0, 0},
			 {{1,0,0},{0,1,0},{0,0,1}},  {{1,0,0},{0,1,0},{0,0,1}}},
			{"L",{ 0, 1, 0},
			 {{0,-1,0},{1,0,0},{0,0,1}}, {{0,1,0},{-1,0,0},{0,0,1}}},
			{"U",{ 0, 0, 1},
			 {{0,0,-1},{0,1,0},{1,0,0}}, {{0,0,1},{0,1,0},{-1,0,0}}},
			{"B",{-1, 0, 0},
			 {{-1,0,0},{0,-1,0},{0,0,1}},{{-1,0,0},{0,-1,0},{0,0,1}}},
			{"R",{ 0,-1, 0},
			 {{0,1,0},{-1,0,0},{0,0,1}}, {{0,-1,0},{1,0,0},{0,0,1}}},
			{"D",{ 0, 0,-1},
			 {{0,0,1},{0,1,0},{-1,0,0}}, {{0,0,-1},{0,1,0},{1,0,0}}}
		};

/*!
	Return the opposite direction of a given direction.
*/
Dir CCLattice::Opposite(Dir const theDir) const
{	//	should be simplified.
	WatchError
	switch(theDir)
	{
		case ccF : return ccB;
		case ccL : return ccR;
		case ccU : return ccD;
		case ccB : return ccF;
		case ccR : return ccL;
		case ccD : return ccU;
	}
	return CCLattice::mDirCount;
	CatchError
}


/*!
	Return the direction for a given vector.
*/
Dir CCLattice::Direction(Vector const & theVector) const
{
	WatchError
    if(Tuple::sqrdist(theVector,0)!=1) return mDirCount;
	Int res = theVector.Comp(0)*4+theVector.Comp(1)*2+theVector.Comp(2)*1;
	switch(res)
	{
	    case -1 :   return ccD;
	    case 1  :   return ccU;
	    case -2 :   return ccR;
	    case 2  :   return ccL;
	    case -4 :   return ccB;
	    case 4 :   return ccF;
	    default: return CCLattice::mDirCount;
	}
	CatchError
}


/*!
	The destructor.
*/
CCLattice::~CCLattice()
{
	WatchError
	CatchError
}


/*!
	The constructor
*/
CCLattice::CCLattice() :
	Lattice(LatticeCls<CCLattice>::ModelId, mDirCount, mSqrNeighDist)
{
	WatchError
	Alert(mSpaceDim != Space::Dimen(), eDimensionMismatch);
	for(Dir tDir = 0; tDir < mDirCount; ++tDir)
	{
		for(Idx tIdx = 0; tIdx <= mNameSize; ++tIdx)
			mDirNames[tDir][tIdx] = ccData[tDir].name[tIdx];
		for(Cmp tCmp1 = 0; tCmp1 < mSpaceDim; ++tCmp1)
		{
			mDirVecs[tDir].Comp(tCmp1) = ccData[tDir].vec[tCmp1];
			for(Cmp tCmp2 = 0; tCmp2 < mSpaceDim; ++tCmp2)
			{
				mDirMats[tDir].Comp(tCmp1,tCmp2) = ccData[tDir].mat[tCmp1][tCmp2];
				mDirInvs[tDir].Comp(tCmp1,tCmp2) = ccData[tDir].inv[tCmp1][tCmp2];
			}
		}
	}
	CatchError
}


/*!
	The duplicator.
*/
CCLattice::CCLattice(CCLattice const & that) : Lattice(that)
{
	WatchError
	for(Dir tDir = 0; tDir < mDirCount; ++tDir)
	{
		for(Idx tIdx = 0; tIdx <= mNameSize; ++tIdx)
			mDirNames[tDir][tIdx] = that.mDirNames[tDir][tIdx];
		for(Cmp tCmp1 = 0; tCmp1 < mSpaceDim; ++tCmp1)
		{
			mDirVecs[tDir].Comp(tCmp1) = that.mDirVecs[tDir].Comp(tCmp1);
			for(Cmp tCmp2 = 0; tCmp2 < mSpaceDim; ++tCmp2)
			{
				mDirMats[tDir].Comp(tCmp1,tCmp2) = that.mDirMats[tDir].Comp(tCmp1,tCmp2);
				mDirInvs[tDir].Comp(tCmp1,tCmp2) = that.mDirInvs[tDir].Comp(tCmp1,tCmp2);
			}
		}
	}
	CatchError
}



/*!
	The assigner.
*/
CCLattice const & CCLattice::operator = (CCLattice const & that)
{
	WatchError
	if (this == &that) return *this;
	Lattice::operator=(that);
	for(Dir tDir = 0; tDir < mDirCount; ++tDir)
	{
		for(Idx tIdx = 0; tIdx <= mNameSize; ++tIdx)
			mDirNames[tDir][tIdx] = that.mDirNames[tDir][tIdx];
		for(Cmp tCmp1 = 0; tCmp1 < mSpaceDim; ++tCmp1)
		{
			mDirVecs[tDir].Comp(tCmp1) = that.mDirVecs[tDir].Comp(tCmp1);
			for(Cmp tCmp2 = 0; tCmp2 < mSpaceDim; ++tCmp2)
			{
				mDirMats[tDir].Comp(tCmp1,tCmp2) = that.mDirMats[tDir].Comp(tCmp1,tCmp2);
				mDirInvs[tDir].Comp(tCmp1,tCmp2) = that.mDirInvs[tDir].Comp(tCmp1,tCmp2);
			}
		}
	}
	return *this;
	CatchError
}



/*!
	Return the direction names.
*/
Drn CCLattice::DirName(Dir const theDir) const
{
	WatchError
	Warn(theDir >= mDirCount, eInvalidDir);
	return mDirNames[theDir];
	CatchError
}



/*!
	Return the direction vectors.
*/
Tuple const & CCLattice::DirVec(Dir const theDir) const
{
	WatchError
	Warn(theDir >= mDirCount, eInvalidDir);
	return mDirVecs[theDir];
	CatchError
}



/*!
	Return the direction matrices.
*/
Matrix const & CCLattice::DirMat(Dir const theDir) const
{
	WatchError
	Warn(theDir >= mDirCount, eInvalidDir);
	return mDirMats[theDir];
	CatchError
}



/*!
	Return the direction inverse matrices.
*/
Matrix const & CCLattice::DirInv(Dir const theDir) const
{
	WatchError
	Warn(theDir >= mDirCount, eInvalidDir);
	return mDirInvs[theDir];
	CatchError
}




/*!
	Return whether a vector is a valid lattice direction.
*/
B CCLattice::ValidDir(Vector const & theVector) const
{
	WatchError
	return Direction(theVector) != mDirCount;
	CatchError
}



/*!
	Return the absolute direction for a given relative direction.
*/
Dir CCLattice::RelToAbsDir(Dir const relDir, Matrix const & baseMatrix)const
{
	WatchError
	Warn(relDir >= mDirCount, eInvalidDir);
	return Direction(mDirVecs[relDir] * baseMatrix);
	CatchError
}




/*!
	Return the relative direction for a given absolute direction.
*/
Dir CCLattice::AbsToRelDir(Matrix const & baseMatrix, Dir const absDir)const
{
	WatchError
	Warn(absDir >= mDirCount, eInvalidDir);
	return Direction(baseMatrix * mDirVecs[absDir]);
	CatchError
}



/*!
	Update the base matrix for relative to absolute direction conversion.
*/
void CCLattice::updtRelToAbsBase(Dir const theDir,Matrix & baseMatrix)const
{
	WatchError
	Matrix tMatrix = baseMatrix;
	baseMatrix = mDirMats[theDir]*tMatrix ;
	CatchError
}



/*!
	Update the base matrix for abolute to relative direction conversion.
*/
void CCLattice::updtAbsToRelBase(Matrix & baseMatrix,Dir const theDir )const
{
	WatchError
	Matrix tMatrix = baseMatrix;
	baseMatrix = tMatrix*mDirInvs[theDir];
	CatchError
}


/*!
	Return whether tow given points are neighbours on the lattice.
*/
B CCLattice::areNeighbors(const Point & point1,const Point & point2) const
{
	WatchError
    return ValidDir(point1-point2);
	CatchError
}

/*!
    Return the lattice Converter

*/
R CCLattice::latticeConverter() const
{
        WatchError
        return 3.8;
        CatchError
}

closePlatypusSpace
