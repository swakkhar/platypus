/*!
	@file lattices/fcclattice.hpp
	@brief The prototype file for FccLattice class.
	@details This is the prototype file for FccLattice class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/


#include "pspl/lattices/fcclattice.hpp"


openPlatypusSpace


/*!
	@struct FccData
	@brief A structure to represent vital data for fcc-lattice.
	@details The structure represents vital data for fcc-lattice.
*/
struct FccData
{
			char name[3]; 	//!< name of direction
			Spc vec[3];     //!< absolute direction vector
			Spc mat[3][3];  //!< matrix for applying relative direction
			Spc inv[3][3]; //!< inverse of matrix (we are lazy ;-))
};



/*!
	@var fccData
	@brief Vital data for fcc-lattice.
	@details The vital data for fcc-lattice.
*/
FccData const fccData[] = {
			{"FL",{1,1,0}, {{1,0,0},{0,1,0},{0,0,1}},   {{1,0,0},{0,1,0},{0,0,1}}},
			{"LU",{0,1,1}, {{0,0,-1},{0,1,0},{1,0,0}},  {{0,0,1},{0,1,0},{-1,0,0}}},
			{"FU",{1,0,1}, {{1,0,0},{0,0,-1},{0,1,0}},  {{1,0,0},{0,0,1},{0,-1,0}}},
			{"BL",{-1,1,0}, {{0,-1,0},{1,0,0},{0,0,1}},  {{0,1,0},{-1,0,0},{0,0,1}}},
			{"RU",{0,-1,1}, {{0,0,-1},{-1,0,0},{0,1,0}}, {{0,-1,0},{0,0,1},{-1,0,0}}},
			{"BU",{-1,0,1}, {{0,-1,0},{0,0,-1},{1,0,0}}, {{0,0,1},{-1,0,0},{0,-1,0}}},
			{"FR",{1,-1,0}, {{0,1,0},{-1,0,0},{0,0,1}},  {{0,-1,0},{1,0,0},{0,0,1}}},
			{"LD",{0,1,-1}, {{0,0,1},{0,1,0},{-1,0,0}},  {{0,0,-1},{0,1,0},{1,0,0}}},
			{"FD",{1,0,-1}, {{1,0,0},{0,0,1},{0,-1,0}},  {{1,0,0},{0,0,-1},{0,1,0}}},
			{"BR",{-1,-1,0}, {{-1,0,0},{0,-1,0},{0,0,1}}, {{-1,0,0},{0,-1,0},{0,0,1}}},
			{"RD",{0,-1,-1}, {{0,0,1},{-1,0,0},{0,-1,0}}, {{0,-1,0},{0,0,-1},{1,0,0}}},
			{"BD",{-1,0,-1}, {{0,-1,0},{0,0,1},{-1,0,0}}, {{0,0,-1},{-1,0,0},{0,1,0}}}
			};


/*!
	Return the opposite direction of a given direction.
*/
Dir FccLattice::Opposite(Dir const theDir) const
{	//	should be simplified.
	WatchError
	switch(theDir)
	{
		case fccFL : return fccBR;
		case fccLU : return fccRD;
		case fccFU : return fccBD;
		case fccBL : return fccFR;
		case fccRU : return fccLD;
		case fccBU : return fccFD;
		case fccBR : return fccFL;
		case fccRD : return fccLU;
		case fccBD : return fccFU;
	}
	return FccLattice::mDirCount;
	CatchError
}


/*!
	Return the direction for a given vector.
*/
Dir FccLattice::Direction(Vector const & theVector) const
{
	WatchError
	switch(theVector.Comp(0))
	{
		case -1 :
			switch(theVector.Comp(1))
			{
				case -1 :
					if (theVector.Comp(2) == 0) return fccBR;
					else return mDirCount;
				case 0 :
					switch(theVector.Comp(2))
					{
						case -1 : return fccBD;
						case 1  : return fccBU;
						default : return mDirCount;
					}
				case 1 :
					if (theVector.Comp(2) == 0) return fccBL;
					else return mDirCount;
				default:
					return mDirCount;
			}

		case 0 :
			switch(theVector.Comp(1))
			{
				case -1 :
					switch(theVector.Comp(2))
					{
						case -1 : return fccRD;
						case 1 : return fccRU;
						default: return mDirCount;
					}
				case 1:
					switch(theVector.Comp(2))
					{
						case -1 : return fccLD;
						case 1 : return fccLU;
						default: return mDirCount;
					}
				default: return mDirCount;
			}

		case 1 :
			switch(theVector.Comp(1))
			{
				case -1 :
					if (theVector.Comp(2) == 0) return fccFR;
					else return mDirCount;
				case 0 :
					switch(theVector.Comp(2))
					{
						case -1 : return fccFD;
						case 1  : return fccFU;
						default : return mDirCount;
					}
				case 1 :
					if (theVector.Comp(2) == 0) return fccFL;
					else return mDirCount;
				default:
					return mDirCount;
			}

		default:
			return mDirCount;
	}
	CatchError
}


/*!
	The destructor.
*/
FccLattice::~FccLattice()
{
	WatchError
	CatchError
}


/*!
	The constructor
*/
FccLattice::FccLattice() :
	Lattice(LatticeCls<FccLattice>::ModelId, mDirCount, mSqrNeighDist)
{
	WatchError
	Alert(mSpaceDim != Space::Dimen(), eDimensionMismatch);
	for(Dir tDir = 0; tDir < mDirCount; ++tDir)
	{
		for(Idx tIdx = 0; tIdx <= mNameSize; ++tIdx)
			mDirNames[tDir][tIdx] = fccData[tDir].name[tIdx];
		for(Cmp tCmp1 = 0; tCmp1 < mSpaceDim; ++tCmp1)
		{
			mDirVecs[tDir].Comp(tCmp1) = fccData[tDir].vec[tCmp1];
			for(Cmp tCmp2 = 0; tCmp2 < mSpaceDim; ++tCmp2)
			{
				mDirMats[tDir].Comp(tCmp1,tCmp2) = fccData[tDir].mat[tCmp1][tCmp2];
				mDirInvs[tDir].Comp(tCmp1,tCmp2) = fccData[tDir].inv[tCmp1][tCmp2];
			}
		}
	}
	CatchError
}


/*!
	The duplicator.
*/
FccLattice::FccLattice(FccLattice const & that) : Lattice(that)
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
FccLattice const & FccLattice::operator = (FccLattice const & that)
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
Drn FccLattice::DirName(Dir const theDir) const
{
	WatchError
	Warn(theDir >= mDirCount, eInvalidDir);
	return mDirNames[theDir];
	CatchError
}



/*!
	Return the direction vectors.
*/
Tuple const & FccLattice::DirVec(Dir const theDir) const
{
	WatchError
	Warn(theDir >= mDirCount, eInvalidDir);
	return mDirVecs[theDir];
	CatchError
}



/*!
	Return the direction matrices.
*/
Matrix const & FccLattice::DirMat(Dir const theDir) const
{
	WatchError
	Warn(theDir >= mDirCount, eInvalidDir);
	return mDirMats[theDir];
	CatchError
}



/*!
	Return the direction inverse matrices.
*/
Matrix const & FccLattice::DirInv(Dir const theDir) const
{
	WatchError
	Warn(theDir >= mDirCount, eInvalidDir);
	return mDirInvs[theDir];
	CatchError
}




/*!
	Return whether a vector is a valid lattice direction.
*/
B FccLattice::ValidDir(Vector const & theVector) const
{
	WatchError
	return Direction(theVector) != mDirCount;
	CatchError
}



/*!
	Return the absolute direction for a given relative direction.
*/
Dir FccLattice::RelToAbsDir(Dir const relDir, Matrix const & baseMatrix)const
{
	WatchError
	Warn(relDir >= mDirCount, eInvalidDir);
	return Direction(mDirVecs[relDir] * baseMatrix);
	CatchError
}




/*!
	Return the relative direction for a given absolute direction.
*/
Dir FccLattice::AbsToRelDir(Matrix const & baseMatrix, Dir const absDir)const
{
	WatchError
	Warn(absDir >= mDirCount, eInvalidDir);
	return Direction(baseMatrix * mDirVecs[absDir]);
	CatchError
}



/*!
	Update the base matrix for relative to absolute direction conversion.
*/
void FccLattice::updtRelToAbsBase(Dir const theDir,Matrix & baseMatrix)const
{
	WatchError
	Matrix tMatrix = baseMatrix;
	baseMatrix = mDirMats[theDir]*tMatrix ;
	CatchError
}



/*!
	Update the base matrix for abolute to relative direction conversion.
*/
void FccLattice::updtAbsToRelBase(Matrix & baseMatrix,Dir const theDir )const
{
	WatchError
	Matrix tMatrix = baseMatrix;
	baseMatrix = tMatrix*mDirInvs[theDir];
	CatchError
}


/*!
	Return whether tow given points are neighbours on the lattice.
*/
B FccLattice::areNeighbors(const Point & point1,const Point & point2) const
{
	WatchError
    return ValidDir(point1-point2);
	CatchError
}
/*!
    Return lattice unit converter
*/

R FccLattice::latticeConverter() const
{
        WatchError
        return 2.67;
        CatchError
}
closePlatypusSpace
