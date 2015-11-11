/*!
	@file globals/lattice.hpp
	@brief The prototype file for Lattice class.
	@details This is the prototype file for Lattice class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#ifndef LatticeHppIncluded
#define LatticeHppIncluded



#include "pspl/globals/space.hpp"
#include "pspl/globals/tuple.hpp"
#include "pspl/globals/matrix.hpp"



openPlatypusSpace


/*!
	@class Lattice
	@brief A class to represent lattices.
	@details This class represents lattices.
*/
class Lattice
{
	private:
		static Lattice const * mLattice;			//!< Points to the global lattice.

	public:
		static Lattice const & l();					//!< Return the global lattice.
		static void l(Lattice const & theLattice);	//!< Set the global lattice.

	protected:
		Mid mModelId;				//!< Lattice model id.
		Dim mNeighCount;			//!< Neighbour count of a lattice point.
		Dst mSqrNeighDist;			//!< Square distance between neighbours.

	public:
		Mid ModelId() const;			//!< Return the model id for the lattice.
		Dim NeighCount() const;		//!< Return the neighbour count of a latice point
		Dst SqrNeighDist() const;	//!< Return square distance between two neighbours.

		virtual Drn DirName(Dir const theDir) const = 0;				//!< Return direction name.
		virtual Tuple const & DirVec(Dir const theDir) const = 0;		//!< Return direction vector.
		virtual Matrix const & DirMat(Dir const theDir) const = 0;		//!< Return transformation matrix.
		virtual Matrix const & DirInv(Dir const theDir) const = 0;		//!< Return inverse transformation matrix.
		virtual Dir Opposite(Dir const theDir) const = 0;				//!< Return opposite direction.

		virtual Dir Direction(Vector const & theVector) const = 0;		//!< Return direction of the vector.
		virtual B ValidDir(Vector const & theVector) const = 0;			//!< Return a vector is a valid direction?

		virtual Dir AbsToRelDir(Matrix const & baseMatrix, Dir const absDir) const = 0;	//!< Convert an absolute direction into relative direction.
        virtual void updtAbsToRelBase(Matrix & baseMatrix, Dir const theDir) const = 0;	//!< Update base matrix for absolute to relative transformation.

        virtual Dir RelToAbsDir(Dir const relDir, Matrix const & baseMatrix) const = 0;	//!< Convert a relative direction into absolute direction.
		virtual void updtRelToAbsBase(Dir const theDir, Matrix & baseMatrix) const = 0;	//!< Update base matrix for relative to absolute transfomration.

        virtual B areNeighbors(const Point & point1,const Point & point2) const = 0;	//!< Return whether two points are neighbours to each other on the lattice.


        virtual R latticeConverter() const = 0;

	protected:

		virtual ~Lattice(); 															//!< Destructor.
		Lattice(Lattice const & that);													//!< Duplicator.
		Lattice const & operator = (Lattice const & that); 								//!< Assigner.
		Lattice(Mid const ModelId, Dim const theNeighCount, Dim const theSqrNeighDist); 	//!< Initialiser.
};


/*!
	@class LatticeId
	@brief A class to represent lattice id numbers.
	@details This class represents lattice id numbers.
*/
template <Mid ModelId> class LatticeId {};


/*!
	@class LatticeCls
	@brief A class to represent lattice class meta data.
	@details This class represents lattice class meta data.
*/
template <class ClassName> class LatticeCls{};


/*!
	@def regLattice(Name, Model)
	@brief Register a lattice class.
	@details Register a lattice class.
*/
#define regLattice(Name, Model) class Name; \
	template<> class LatticeId<Model> {public: typedef Name Class;}; \
	template<> class LatticeCls<Name> {public: enum {ModelId = Model};};


/*!
	Return the global lattice.
*/
inline Lattice const & Lattice::l()
{
	WatchError
	Warn(!mLattice, eEmptyLattice);
	return *mLattice;
	CatchError
}



/*!
	Set the global lattice.
*/
inline void Lattice::l(Lattice const & theLattice)
{
	WatchError
	Warn(mLattice, eNonEmptyLattice);
	mLattice = &theLattice;
	CatchError
}




/*!
	Return the neighbour count of the lattice.
*/
inline Mid Lattice::ModelId() const
{
	WatchError
	return mModelId;
	CatchError
}



/*!
	Return the neighbour count of the lattice.
*/
inline Dim Lattice::NeighCount() const
{
	WatchError
	return mNeighCount;
	CatchError
}



/*!
	Return the neighbour count of the lattice.
*/
inline Dst Lattice::SqrNeighDist() const
{
	WatchError
	return mSqrNeighDist;
	CatchError
}



closePlatypusSpace



#endif //LatticeHppIncluded
