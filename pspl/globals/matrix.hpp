/*!
	@file globals/matrix.hpp
	@brief The prototype file for Matrix class.
	@details This is the prototype file for Matrix class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/


#ifndef MatrixHppIncluded
#define MatrixHppIncluded


#include "pspl/globals/space.hpp"



openPlatypusSpace



/*!
	Forward declaration.
*/
class Tuple;



/*!
	@class Matrix
	@brief A class to represent matrices.
	@details This class represents matrices.
*/
class Matrix
{
	private:

		Spc ** mComps;	//!< Components.

	public:

		Spc const & Comp(Cmp const theCmp1, Cmp const theCmp2) const;	//!< Components.
		Spc & Comp(Cmp const theCmp1, Cmp const theCmp2);				//!< Components.

		Matrix();												//!< Constrcuctor.
		~Matrix();												//!< Destructor.
		Matrix(Spc const theSpc);								//!< Initialiser.
		Matrix(Spc const ** theComps);							//!< Initialiser.
		Matrix(Matrix const & that);							//!< Duplicator.
		Matrix const & operator = (Matrix const & that);		//!< Assigner.

		Matrix const & operator += (Matrix const & that);		//!< Addition.
		Matrix const & operator -= (Matrix const & that);		//!< Subtraction.
		Matrix const & operator *= (Matrix const & that);		//!< Multiplication.

		Matrix const & operator *= (Spc const theFactor);		//!< Multiplication.
		Matrix const & operator += (Spc const theBias);			//!< Addition.
		Matrix const & operator -= (Spc const theBias);			//!< Subtraction.

		Matrix const operator - () const;						//!< Negation.
		Matrix const operator + (Matrix const & that) const;	//!< Addition.
		Matrix const operator - (Matrix const & that) const;	//!< Subtraction.
		Matrix const operator * (Matrix const & that) const;	//!< Multiplication.

		Matrix const operator + (Spc const theBias) const;		//!< Addition.
		Matrix const operator - (Spc const theBias) const;		//!< Subtraction.
		Matrix const operator * (Spc const theFactor) const;	//!< Multiplication.
		Tuple const operator * (Tuple const & theTuple) const;	//!< Multiplication.
};



/*!
	Return a component.
*/
inline Spc const & Matrix::Comp(Cmp const theCmp1, Cmp const theCmp2) const
{
	WatchError
	Warn(theCmp1 >= Space::Dimen(), eInvalidIndex);
	Warn(theCmp2 >= Space::Dimen(), eInvalidIndex);
	return mComps[theCmp1][theCmp2];
	CatchError
}



/*!
	Return a component.
*/
inline Spc & Matrix::Comp(Cmp const theCmp1, Cmp const theCmp2)
{
	WatchError
	Warn(theCmp1 >= Space::Dimen(), eInvalidIndex);
	Warn(theCmp2 >= Space::Dimen(), eInvalidIndex);
	return mComps[theCmp1][theCmp2];
	CatchError
}


closePlatypusSpace


#endif //MatrixHppIncluded
