/*!
	@file lattices/fcclattice.hpp
	@brief The prototype file for FccLattice class.
	@details This is the prototype file for FccLattice class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#ifndef FccLatticeHhIncluded
#define FccLatticeHhIncluded



#include "pspl/globals/lattice.hpp"



openPlatypusSpace


/*!
	Register FCC lattice class.
*/
regLattice(FccLattice, castN(1));


/*!
	@enum FccDir
	@brief A enum to enumerate of FCC basis directions.
	@details This enum enumerates the FCC basis directions.
*/
enum FccDir	//changing values here requires chaging the arrays in the cpp file.
{
	fccFL = 0,
	fccLU = 1,
	fccFU = 2,
	fccBL = 3,
	fccRU = 4,
	fccBU = 5,
	fccFR = 6,
	fccLD = 7,
	fccFD = 8,
	fccBR = 9,
	fccRD = 10,
	fccBD = 11
};



/*!
	@class FccLattice
	@brief A class to represent fcc-lattices.
	@details This class represent fcc-lattices.
*/
class FccLattice : public Lattice
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
			mDirCount = 12,		//!< Direction count.
			mNameSize = 2,		//!< Direction name size.
			mSqrNeighDist = 2	//!< Square of unit distance.
		};

	private:
		C mDirNames[mDirCount][mNameSize + 1];		//!< Direction names.
		Vector mDirVecs[mDirCount];					//!< Direction vectors.
		Matrix mDirMats[mDirCount];					//!< Direction matrices.
		Matrix mDirInvs[mDirCount];					//!< Direction inv matrices.

	public:

		FccLattice();									 			//!< Constructor.
		FccLattice(FccLattice const & that);						//!< Duplicator.
		FccLattice const & operator = (FccLattice const & that);	//!< Assignment.
		~FccLattice();												//!< Destructor.

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



#endif //FccLatticeHhIncluded
