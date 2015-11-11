/*!
	@file lattices/TLattice.hpp
	@brief The prototype file for TLattice class.
	@details This is the prototype file for TLattice class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#ifndef TLatticeHhIncluded
#define TLatticeHhIncluded



#include "pspl/globals/lattice.hpp"



openPlatypusSpace


/*!
	Register Tlattice lattice class.
*/
regLattice(TLattice, castN(4));


/*!
	@enum CcDir
	@brief A enum to enumerate of CC basis directions.
	@details This enum enumerates the CC basis directions.
*/
enum KwDir	//changing values here requires chaging the arrays in the cpp file.
{
	kwFFR=0,
	kwFFL=1,
	kwFFU=2,
	kwFFD=3,
	kwBBR=4,
	kwBBL=5,
	kwBBU=6,
	kwBBD=7,
	kwRRF=8,
	kwRRB=9,
	kwRRU=10,
	kwRRD=11,
	kwLLF=12,
	kwLLB=13,
	kwLLU=14,
	kwLLD=15,
	kwUUF=16,
	kwUUB=17,
	kwUUR=18,
	kwUUL=19,
	kwDDF=20,
	kwDDB=21,
	kwDDR=22,
	kwDDL=23,
};



/*!
	@class TLattice
	@brief A class to represent cc-lattices.
	@details This class represent cc-lattices.
*/
class TLattice : public Lattice
{
	public:
		virtual Drn DirName(Dir const theDir) const;					//!< Return a direction name.
		virtual Vector const & DirVec(Dir const theDir) const;			//!< Return a direction vector.
		virtual Matrix const & DirMat(Dir const theDir) const;			//!< Return a direction matrix.
		virtual Matrix const & DirInv(Dir const theDir) const;			//!< Return a direction inverse matrix.

	private:
		/*!
			@enum Settings
			@brief Internal settings.
			@details
		*/
		enum Settings
		{
			mSpaceDim = 3,		//!< Space dimension.
			mDirCount = 24,		//!< Direction count.
			mNameSize = 3,		//!< Direction name size.
			mSqrNeighDist = 1	//!< Square of unit distance.
		};

	private:
		C mDirNames[mDirCount][mNameSize + 1];		//!< Direction names.
		Vector mDirVecs[mDirCount];					//!< Direction vectors.
		Matrix mDirMats[mDirCount];					//!< Direction matrices.
		Matrix mDirInvs[mDirCount];					//!< Direction inv matrices.

	public:

		TLattice();									 			//!< Constructor.
		TLattice(TLattice const & that);						//!< Duplicator.
		TLattice const & operator = (TLattice const & that);	//!< Assignment.
		~TLattice();												//!< Destructor.

		virtual Dir Opposite(Dir const theDir) const;				//!< Opposite direction.
		virtual Dir Direction(Vector const & theVector) const;		//!< Direction of a vector.

		virtual B ValidDir(Vector const & theVector) const;							//!< Whether direction is a valid basis.
		virtual Dir RelToAbsDir(Dir const relDir, Matrix const & baseMatrix) const;	//!< Convert from relative to absolute direction.
		virtual Dir AbsToRelDir(Matrix const & baseMatrix, Dir const absDir) const;	//!< Convert from absolute to relative direction.
		virtual void updtRelToAbsBase(Dir const theDir,Matrix & baseMatrix) const;	//!< Update base matric for relative to absolute direction conversion.
		virtual void updtAbsToRelBase(Matrix & baseMatrix,Dir const theDir) const;	//!< Update base matric for absolute to relative direction conversion.

        virtual B areNeighbors(const Point & point1,const Point & ponit2) const;	//!< Whether two points are neighbour to each other on the lattice

        virtual R latticeConverter() const;
};



closePlatypusSpace



#endif //TLatticeHhIncluded


