/*!
	@file lattices/KwLattice.hpp
	@brief The prototype file for KwLattice class.
	@details This is the prototype file for KwLattice class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/


#include "pspl/lattices/kwlattice.hpp"


openPlatypusSpace


/*!
	@struct FccData
	@brief A structure to represent vital data for fcc-lattice.
	@details The structure represents vital data for fcc-lattice.
*/
struct KwData
{
			char name[4]; 	//!< name of direction
			Spc vec[3];     //!< absolute direction vector
			Spc mat[3][3];  //!< matrix for applying relative direction
			Spc inv[3][3]; //!< inverse of matrix (we are lazy ;-))
};



/*!
	@var fccData
	@brief Vital data for fcc-lattice.
	@details The vital data for fcc-lattice.
*/
KwData const kwData[] = {
			{"FFR",{2,1,0},{{0,0,0},{0,0,0},{0,0,0}}, {{0,0,0},{0,0,0},{0,0,0}}},
			{"FFL",{2,-1,0},{{0,0,0},{0,0,0},{0,0,0}}, {{0,0,0},{0,0,0},{0,0,0}}},
			{"FFU",{2,0,1},{{0,0,0},{0,0,0},{0,0,0}}, {{0,0,0},{0,0,0},{0,0,0}}},
			{"FFD",{2,0,-1},{{0,0,0},{0,0,0},{0,0,0}}, {{0,0,0},{0,0,0},{0,0,0}}},
			{"BBR",{-2,1,0},{{0,0,0},{0,0,0},{0,0,0}}, {{0,0,0},{0,0,0},{0,0,0}}},
			{"BBL",{-2,-1,0},{{0,0,0},{0,0,0},{0,0,0}}, {{0,0,0},{0,0,0},{0,0,0}}},
			{"BBU",{-2,0,1},{{0,0,0},{0,0,0},{0,0,0}}, {{0,0,0},{0,0,0},{0,0,0}}},
			{"BBD",{-2,0,-1},{{0,0,0},{0,0,0},{0,0,0}}, {{0,0,0},{0,0,0},{0,0,0}}},
			{"RRF",{1,2,0},{{0,0,0},{0,0,0},{0,0,0}}, {{0,0,0},{0,0,0},{0,0,0}}},
			{"RRB",{-1,2,0},{{0,0,0},{0,0,0},{0,0,0}}, {{0,0,0},{0,0,0},{0,0,0}}},
			{"RRU",{0,2,1},{{0,0,0},{0,0,0},{0,0,0}}, {{0,0,0},{0,0,0},{0,0,0}}},
			{"RRD",{0,2,-1},{{0,0,0},{0,0,0},{0,0,0}}, {{0,0,0},{0,0,0},{0,0,0}}},
			{"LLF",{1,-2,0},{{0,0,0},{0,0,0},{0,0,0}}, {{0,0,0},{0,0,0},{0,0,0}}},
			{"LLB",{-1,-2,0},{{0,0,0},{0,0,0},{0,0,0}}, {{0,0,0},{0,0,0},{0,0,0}}},
			{"LLU",{0,-2,1},{{0,0,0},{0,0,0},{0,0,0}}, {{0,0,0},{0,0,0},{0,0,0}}},
			{"LLD",{0,-2,-1},{{0,0,0},{0,0,0},{0,0,0}}, {{0,0,0},{0,0,0},{0,0,0}}},
			{"UUF",{1,0,2},{{0,0,0},{0,0,0},{0,0,0}}, {{0,0,0},{0,0,0},{0,0,0}}},
			{"UUB",{-1,0,2},{{0,0,0},{0,0,0},{0,0,0}}, {{0,0,0},{0,0,0},{0,0,0}}},
			{"UUR",{0,1,2},{{0,0,0},{0,0,0},{0,0,0}}, {{0,0,0},{0,0,0},{0,0,0}}},
			{"UUL",{0,-1,2},{{0,0,0},{0,0,0},{0,0,0}}, {{0,0,0},{0,0,0},{0,0,0}}},
			{"DDF",{1,0,-2},{{0,0,0},{0,0,0},{0,0,0}}, {{0,0,0},{0,0,0},{0,0,0}}},
			{"DDB",{-1,0,-2},{{0,0,0},{0,0,0},{0,0,0}}, {{0,0,0},{0,0,0},{0,0,0}}},
			{"DDR",{0,1,-2},{{0,0,0},{0,0,0},{0,0,0}}, {{0,0,0},{0,0,0},{0,0,0}}},
			{"DDL",{0,-1,-2},{{0,0,0},{0,0,0},{0,0,0}}, {{0,0,0},{0,0,0},{0,0,0}}},
        };


/*!
	Return the opposite direction of a given direction.
*/
Dir KwLattice::Opposite(Dir const theDir) const
{	//	should be simplified.
	WatchError
	switch(theDir)
	{
        case kwFFR : return kwFFL;
        case kwFFL : return kwFFR;
        case kwFFU : return kwFFD;
        case kwFFD : return kwFFU;
        case kwBBR : return kwBBL;
        case kwBBL : return kwBBR;
        case kwBBU : return kwBBD;
        case kwBBD : return kwBBU;
        case kwRRF : return kwRRB;
        case kwRRB : return kwRRF;
        case kwRRU : return kwRRD;
        case kwRRD : return kwRRU;
        case kwLLF : return kwLLB;
        case kwLLB : return kwLLF;
        case kwLLU : return kwLLD;
        case kwLLD : return kwLLU;
        case kwUUF : return kwUUB;
        case kwUUB : return kwUUF;
        case kwUUR : return kwUUL;
        case kwUUL : return kwUUR;
        case kwDDF : return kwDDB;
        case kwDDB : return kwDDF;
        case kwDDR : return kwDDL;
        case kwDDL : return kwDDR;
	}
	return KwLattice::mDirCount;
	CatchError
}


/*!
	Return the direction for a given vector.
*/
Dir KwLattice::Direction(Vector const & theVector) const
{
	WatchError

	switch(theVector.Comp(0))
	{
		case 2 :
			switch(theVector.Comp(1))
			{
				case 1 :
					if (theVector.Comp(2) == 0) return kwFFR;
					else return mDirCount;
				case 0 :
					switch(theVector.Comp(2))
					{
						case -1 : return kwFFD;
						case 1  : return kwFFU;
						default : return mDirCount;
					}
				case -1 :
					if (theVector.Comp(2) == 0) return kwFFL;
					else return mDirCount;
				default:
					return mDirCount;
			}
        case -2:
            switch(theVector.Comp(1))
            {
                case 1:
                    if(theVector.Comp(2) == 0) return kwBBR;
                    else return mDirCount;
                case 0:
                    switch(theVector.Comp(2))
                    {
                        case -1 : return kwBBD;
                        case 1 : return kwBBU;
                        default : return mDirCount;
                    }
                case -1 :
                    if(theVector.Comp(2) == 0) return kwBBL;
                    else return mDirCount;
                default:
                    return mDirCount;

            }
        case -1:
            switch(theVector.Comp(1))
            {
                case 2:
                    if(theVector.Comp(2) == 0) return kwRRB;
                    else return mDirCount;
                case 0:
                    switch(theVector.Comp(2))
                    {
                        case -2 : return kwDDB;
                        case 2 : return kwUUB;
                        default : return mDirCount;
                    }
                case -2 :
                    if(theVector.Comp(2) == 0) return kwLLB;
                    else return mDirCount;
                default:
                    return mDirCount;

            }
        case 1:
            switch(theVector.Comp(1))
            {
                case 2:
                    if(theVector.Comp(2) == 0) return kwRRF;
                    else return mDirCount;
                case 0:
                    switch(theVector.Comp(2))
                    {
                        case -2 : return kwDDF;
                        case 2 : return kwUUF;
                        default : return mDirCount;
                    }
                case -2 :
                    if(theVector.Comp(2) == 0) return kwLLF;
                    else return mDirCount;
                default:
                    return mDirCount;

            }
        case 0:
            switch(theVector.Comp(1))
            {
                case 1:
                    switch(theVector.Comp(2))
                    {
                        case 2: return kwUUR;
                        case -2: return kwDDR;
                        default: return mDirCount;
                    }
                case -1:
                    switch(theVector.Comp(2))
                    {
                        case 2: return kwUUL;
                        case -2: return kwDDL;
                        default: return mDirCount;
                    }
                case 2:
                    switch(theVector.Comp(2))
                    {
                        case 1: return kwRRU;
                        case -1: return kwRRD;
                        default: return mDirCount;
                    }
                case -2:
                    switch(theVector.Comp(2))
                    {
                        case 1: return kwLLU;
                        case -1: return kwLLD;
                        default: return mDirCount;
                    }
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
KwLattice::~KwLattice()
{
	WatchError
	CatchError
}


/*!
	The constructor
*/
KwLattice::KwLattice() :
	Lattice(LatticeCls<KwLattice>::ModelId, mDirCount, mSqrNeighDist)
{
	WatchError
	Alert(mSpaceDim != Space::Dimen(), eDimensionMismatch);
	for(Dir tDir = 0; tDir < mDirCount; ++tDir)
	{
		for(Idx tIdx = 0; tIdx <= mNameSize; ++tIdx)
			mDirNames[tDir][tIdx] = kwData[tDir].name[tIdx];
		for(Cmp tCmp1 = 0; tCmp1 < mSpaceDim; ++tCmp1)
		{
			mDirVecs[tDir].Comp(tCmp1) = kwData[tDir].vec[tCmp1];
			for(Cmp tCmp2 = 0; tCmp2 < mSpaceDim; ++tCmp2)
			{
				mDirMats[tDir].Comp(tCmp1,tCmp2) = kwData[tDir].mat[tCmp1][tCmp2];
				mDirInvs[tDir].Comp(tCmp1,tCmp2) = kwData[tDir].inv[tCmp1][tCmp2];
			}
		}
	}
	CatchError
}


/*!
	The duplicator.
*/
KwLattice::KwLattice(KwLattice const & that) : Lattice(that)
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
KwLattice const & KwLattice::operator = (KwLattice const & that)
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
Drn KwLattice::DirName(Dir const theDir) const
{
	WatchError
	Warn(theDir >= mDirCount, eInvalidDir);
	return mDirNames[theDir];
	CatchError
}



/*!
	Return the direction vectors.
*/
Tuple const & KwLattice::DirVec(Dir const theDir) const
{
	WatchError
	Warn(theDir >= mDirCount, eInvalidDir);
	return mDirVecs[theDir];
	CatchError
}



/*!
	Return the direction matrices.
*/
Matrix const & KwLattice::DirMat(Dir const theDir) const
{
	WatchError
	Warn(theDir >= mDirCount, eInvalidDir);
	return mDirMats[theDir];
	CatchError
}



/*!
	Return the direction inverse matrices.
*/
Matrix const & KwLattice::DirInv(Dir const theDir) const
{
	WatchError
	Warn(theDir >= mDirCount, eInvalidDir);
	return mDirInvs[theDir];
	CatchError
}




/*!
	Return whether a vector is a valid lattice direction.
*/
B KwLattice::ValidDir(Vector const & theVector) const
{
	WatchError
	return Direction(theVector) != mDirCount;
	CatchError
}



/*!
	Return the absolute direction for a given relative direction.
*/
Dir KwLattice::RelToAbsDir(Dir const relDir, Matrix const & baseMatrix)const
{
	WatchError
	Warn(relDir >= mDirCount, eInvalidDir);
	return Direction(mDirVecs[relDir] * baseMatrix);
	CatchError
}




/*!
	Return the relative direction for a given absolute direction.
*/
Dir KwLattice::AbsToRelDir(Matrix const & baseMatrix, Dir const absDir)const
{
	WatchError
	Warn(absDir >= mDirCount, eInvalidDir);
	return Direction(baseMatrix * mDirVecs[absDir]);
	CatchError
}



/*!
	Update the base matrix for relative to absolute direction conversion.
*/
void KwLattice::updtRelToAbsBase(Dir const theDir,Matrix & baseMatrix)const
{
	WatchError
	Matrix tMatrix = baseMatrix;
	baseMatrix = mDirMats[theDir]*tMatrix ;
	CatchError
}



/*!
	Update the base matrix for abolute to relative direction conversion.
*/
void KwLattice::updtAbsToRelBase(Matrix & baseMatrix,Dir const theDir )const
{
	WatchError
	Matrix tMatrix = baseMatrix;
	baseMatrix = tMatrix*mDirInvs[theDir];
	CatchError
}


/*!
	Return whether tow given points are neighbours on the lattice.
*/
B KwLattice::areNeighbors(const Point & point1,const Point & point2) const
{
	WatchError
    return ValidDir(point1-point2);
	CatchError
}
/*!
    Return lattice unit converter
*/

R KwLattice::latticeConverter() const
{
        WatchError
        return 1.7;
        CatchError
}
closePlatypusSpace
