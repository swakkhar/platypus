/*!
	@file globals/tuple.hpp
	@brief The prototype file for Tuple class.
	@details This is the prototype file for Tuple class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#ifndef TupleHppIncluded
#define TupleHppIncluded



#include "pspl/globals/space.hpp"



openPlatypusSpace



/*!
	Forward declaration.
*/
class Matrix;



/*!
	@class Tuple
	@brief A class to represent tuples.
	@details This is a class to represent tuples.
*/
class Tuple
{
	private:
		Spc * mComps;	//!< Components.

	public:
		Tuple();										//!< Constructor.
		~Tuple();										//!< Destructor.
		Tuple(Spc const theSpc);						//!< Initialiser.
		Tuple(Tuple const & that);						//!< Duplicator.

		Tuple const & operator = (Spc const theSpc);	//!< Assigner.
		Tuple const & operator = (Tuple const & that);	//!< Assigner.

		Spc const & Comp(Cmp const theCmp) const;		//!< Component.
		Spc & Comp(Cmp const theCmp);					//!< Component.

		B operator == (Tuple const & that) const;		//!< Equality.
		B operator != (Tuple const & that) const;		//!< Not-equality.
		B operator < (Tuple const & that) const;		//!< Less-than.
		B operator <= (Tuple const & that) const;		//!< Less-equal.
		B operator > (Tuple const & that) const;		//!< Greater-than.
		B operator >= (Tuple const & that) const;		//!< Greater-equal.

		Tuple const & operator += (Tuple const & that);			//!< Addition.
		Tuple const & operator -= (Tuple const & that);			//!< Subtraction.
		Tuple const & operator *= (Tuple const & that);			//!< Multiplication.
		Tuple const & operator *= (Matrix const & theMatrix);	//!< Multiplication.

		Tuple const & operator += (Spc const & theBias);		//!< Addition.
		Tuple const & operator -= (Spc const & theBias);		//!< Subtraction.
		Tuple const & operator *= (Spc const & theFactor);		//!< Multiplication.

		Tuple operator - () const;								//!< Subtraction.

		Tuple operator + (Tuple const & that) const;			//!< Addition.
		Tuple operator - (Tuple const & that) const;			//!< Subtraction.
		Tuple operator * (Tuple const & that) const;			//!< Multiplication.
		Tuple operator * (Matrix const & theMatrix) const;		//!< Multiplication.

		Tuple operator + (Spc const & theBias) const;			//!< Addition.
		Tuple operator - (Spc const & theBias) const;			//!< Subtraction.
		Tuple operator * (Spc const & theFactor) const;			//!< Multiplication.

		Dst length() const;										//!< Length.
		Tuple unit() const;										//!< Unit tuple.

	public:
		static Dst manh(Tuple const & Tuple1, Tuple const & Tuple);		//!< Manhattan distance.
		static Dst dist(Tuple const & Tuple1, Tuple const & Tuple2);	//!< Cartesian distance.
		static Dst sqrdist(Tuple const & Tuple1, Tuple const & Tuple2);	//!< Square distance.
		static Spc dot(Tuple const & Tubple1, Tuple const & Tuple2);	//!< Dot product.
		static Tuple cross(Tuple const & Tubple1, Tuple const & Tuple2);//!< Cross product.
};


/*!
	@typedef Point.
	@brief Cartesian coordinates.
	@brief Cartesian coordinates.
*/
typedef Tuple Point;


/*!
	@typedef Vector.
	@brief Cartesian vectors.
	@brief Cartesian vectors.
*/
typedef Tuple Vector;



/*!
	@brief Calculate hash value.
	@details Calculate hash value.
*/
inline Hvl calcHash(Tuple const & theTuple)
{
	WatchError
	Hvl tHvl = 0;
	for(Cmp tCmp = 0; tCmp < Space::Dimen(); ++tCmp)
		tHvl = mixHash(tHvl, emu::calcHash(theTuple.Comp(tCmp)));
	return tHvl;
	CatchError
}



/*!
	Return a component.
*/
inline Spc const & Tuple::Comp(Cmp const theCmp) const
{
	WatchError
	Warn(theCmp >= Space::Dimen(), eInvalidIndex);
	return mComps[theCmp];
	CatchError
}



/*!
	Return a component.
*/
inline Spc & Tuple::Comp(Cmp const theCmp)
{
	WatchError
	Warn(theCmp >= Space::Dimen(), eInvalidIndex);
	return mComps[theCmp];
	CatchError
}



/*!
	The constructor.
*/
inline Tuple::Tuple()
{
	WatchError
	mComps = new Spc[Space::Dimen()];
	CatchError
}



/*!
	The destructor.
*/
inline Tuple::~Tuple()
{
	WatchError
	delete [] mComps;
	CatchError
}


closePlatypusSpace



#endif //TupleHhIncluded
