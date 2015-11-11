/*!
	@file globals/matrix.cpp
	@brief The implementation file for Matrix class.
	@details This is the implementation file for Matrix class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/


#include "pspl/globals/matrix.hpp"
#include "pspl/globals/tuple.hpp"



openPlatypusSpace



/*!
	The constructor.
*/
Matrix::Matrix()
{
	WatchError
	mComps = new Spc * [Space::Dimen()];
	for(Cmp tCmp = 0; tCmp < Space::Dimen(); ++tCmp)
		mComps[tCmp] = new Spc [Space::Dimen()];
	CatchError
}


/*!
	The destructor.
*/
Matrix::~Matrix()
{
	WatchError
	for(Cmp tCmp = 0; tCmp < Space::Dimen(); ++tCmp)
		delete [] mComps[tCmp];
	delete [] mComps;
	CatchError
}


/*!
	The initialiser.
*/
Matrix::Matrix(Spc const theSpc)
{
	WatchError
	mComps = new Spc * [Space::Dimen()];
	for(Cmp tCmp1 = 0; tCmp1 < Space::Dimen(); ++tCmp1)
	{
		mComps[tCmp1] = new Spc [Space::Dimen()];
		for(Cmp tCmp2 = 0; tCmp2 < Space::Dimen(); ++tCmp2)
			mComps[tCmp1][tCmp2] = (tCmp1 == tCmp2) ? theSpc : 0;
	}
	CatchError
}



/*!
	The initialiser.
*/
Matrix::Matrix(Spc const ** theComps)
{
	WatchError
	mComps = new Spc * [Space::Dimen()];
	for(Cmp tCmp1 = 0; tCmp1 < Space::Dimen(); ++tCmp1)
	{
		mComps[tCmp1] = new Spc [Space::Dimen()];
		for(Cmp tCmp2 = 0; tCmp2 < Space::Dimen(); ++tCmp2)
			mComps[tCmp1][tCmp2] = theComps[tCmp1][tCmp2];
	}
	CatchError
}


/*!
	The duplicator.
*/
Matrix::Matrix(Matrix const & that)
{
	WatchError
	mComps = new Spc * [Space::Dimen()];
	for(Cmp tCmp1 = 0; tCmp1 < Space::Dimen(); ++tCmp1)
	{
		mComps[tCmp1] = new Spc [Space::Dimen()];
		for(Cmp tCmp2 = 0; tCmp2 < Space::Dimen(); ++tCmp2)
			mComps[tCmp1][tCmp2] = that.mComps[tCmp1][tCmp2];
	}
	CatchError
}



/*!
	The assigner.
*/
Matrix const & Matrix::operator = (Matrix const & that)
{
	WatchError
	if (&that != this)
		for(Cmp tCmp1 = 0; tCmp1 < Space::Dimen(); ++tCmp1)
			for(Cmp tCmp2 = 0; tCmp2 < Space::Dimen(); ++tCmp2)
				mComps[tCmp1][tCmp2] = that.mComps[tCmp1][tCmp2];
	return *this;
	CatchError
}



/*!
	Addition of matrices.
*/
Matrix const & Matrix::operator += (Matrix const & that)
{
	WatchError
	for(Cmp tCmp1 = 0; tCmp1 < Space::Dimen(); ++tCmp1)
		for(Cmp tCmp2 = 0; tCmp2 < Space::Dimen(); ++tCmp2)
			mComps[tCmp1][tCmp2] += that.mComps[tCmp1][tCmp2];
	return *this;
	CatchError
}


/*!
	Subtraction of matrices
*/
Matrix const & Matrix::operator -= (Matrix const & that)
{
	WatchError
	for(Cmp tCmp1 = 0; tCmp1 < Space::Dimen(); ++tCmp1)
		for(Cmp tCmp2 = 0; tCmp2 < Space::Dimen(); ++tCmp2)
			mComps[tCmp1][tCmp2] -= that.mComps[tCmp1][tCmp2];
	return *this;
	CatchError
}

/*!
	Multiplication of matrics.
*/
Matrix const & Matrix::operator *= (Matrix const & that)
{
	WatchError
	for(Cmp tCmp1 = 0; tCmp1 < Space::Dimen(); ++tCmp1)
	{
		Tuple tTuple;
		for (Cmp tCmp3 = 0; tCmp3 < Space::Dimen(); ++tCmp3)
			tTuple.Comp(tCmp3) = mComps[tCmp1][tCmp3];

		for(Cmp tCmp2 = 0; tCmp2 < Space::Dimen(); ++tCmp2)
		{
			mComps[tCmp1][tCmp2] = 0;
			for (Cmp tCmp3 = 0; tCmp3 < Space::Dimen(); ++tCmp3)
				mComps[tCmp1][tCmp2] += tTuple.Comp(tCmp3) * that.mComps[tCmp3][tCmp2];
		}
	}
	return *this;
	CatchError
}



/*!
	Multiplication with a scalar factor.
*/
Matrix const & Matrix::operator *= (Spc const theFactor)
{
	WatchError
	for(Cmp tCmp1 = 0; tCmp1 < Space::Dimen(); ++tCmp1)
		for(Cmp tCmp2 = 0; tCmp2 < Space::Dimen(); ++tCmp2)
			mComps[tCmp1][tCmp2] *= theFactor;
	return *this;
	CatchError
}



/*!
	Addition with a scalar factor.
*/
Matrix const & Matrix::operator += (Spc const theBias)
{
	WatchError
	for(Cmp tCmp1 = 0; tCmp1 < Space::Dimen(); ++tCmp1)
		for(Cmp tCmp2 = 0; tCmp2 < Space::Dimen(); ++tCmp2)
			mComps[tCmp1][tCmp2] += theBias;
	return *this;
	CatchError
}



/*!
	Subtraction with a scalar factor.
*/
Matrix const & Matrix::operator -= (Spc const theBias)
{
	WatchError
	for(Cmp tCmp1 = 0; tCmp1 < Space::Dimen(); ++tCmp1)
		for(Cmp tCmp2 = 0; tCmp2 < Space::Dimen(); ++tCmp2)
			mComps[tCmp1][tCmp2] -= theBias;
	return *this;
	CatchError
}



/*!
	Negation of the matrix.
*/
Matrix const Matrix::operator - () const
{
	WatchError
	Matrix tMatrix;
	for(Cmp tCmp1 = 0; tCmp1 < Space::Dimen(); ++tCmp1)
		for(Cmp tCmp2 = 0; tCmp2 < Space::Dimen(); ++tCmp2)
			tMatrix.mComps[tCmp1][tCmp2] = -mComps[tCmp1][tCmp2];
	return tMatrix;
	CatchError
}



/*!
	Addition of two matrices.
*/
Matrix const Matrix::operator + (Matrix const & that) const
{
	WatchError
	Matrix tMatrix;
	for(Cmp tCmp1 = 0; tCmp1 < Space::Dimen(); ++tCmp1)
		for(Cmp tCmp2 = 0; tCmp2 < Space::Dimen(); ++tCmp2)
			tMatrix.mComps[tCmp1][tCmp2] = mComps[tCmp1][tCmp2] + that.mComps[tCmp1][tCmp2];
	return tMatrix;
	CatchError
}


/*!
	Subtraction with a scalar factor.
*/
Matrix const Matrix::operator - (Matrix const & that) const
{
	WatchError
	Matrix tMatrix;
	for(Cmp tCmp1 = 0; tCmp1 < Space::Dimen(); ++tCmp1)
		for(Cmp tCmp2 = 0; tCmp2 < Space::Dimen(); ++tCmp2)
			tMatrix.mComps[tCmp1][tCmp2] = mComps[tCmp1][tCmp2] - that.mComps[tCmp1][tCmp2];
	return *this;
	CatchError
}



/*!
	Multiplication of two matrices.
*/
Matrix const Matrix::operator * (Matrix const & that) const
{
	WatchError
	Matrix tMatrix;
	for(Cmp tCmp1 = 0; tCmp1 < Space::Dimen(); ++tCmp1)
		for(Cmp tCmp2 = 0; tCmp2 < Space::Dimen(); ++tCmp2)
		{
			tMatrix.mComps[tCmp1][tCmp2] = 0;
			for (Cmp tCmp3 = 0; tCmp3 < Space::Dimen(); ++tCmp3)
				tMatrix.mComps[tCmp1][tCmp2] += mComps[tCmp1][tCmp3] * that.mComps[tCmp3][tCmp2];
		}
	return tMatrix;
	CatchError
}



/*!
	Addition with a scalar.
*/
Matrix const Matrix::operator + (Spc const theBias) const
{
	WatchError
	Matrix tMatrix;
	for(Cmp tCmp1 = 0; tCmp1 < Space::Dimen(); ++tCmp1)
		for(Cmp tCmp2 = 0; tCmp2 < Space::Dimen(); ++tCmp2)
			tMatrix.mComps[tCmp1][tCmp2] = mComps[tCmp1][tCmp2] + theBias;
	return tMatrix;
	CatchError
}



/*!
	Subtraction with a scalar.
*/
Matrix const Matrix::operator - (Spc const theBias) const
{
	WatchError
	Matrix tMatrix;
	for(Cmp tCmp1 = 0; tCmp1 < Space::Dimen(); ++tCmp1)
		for(Cmp tCmp2 = 0; tCmp2 < Space::Dimen(); ++tCmp2)
			tMatrix.mComps[tCmp1][tCmp2] = mComps[tCmp1][tCmp2] - theBias;
	return *this;
	CatchError
}



/*!
	Multiplication with a scalar.
*/
Matrix const Matrix::operator * (Spc const theFactor) const
{
	WatchError
	Matrix tMatrix;
	for(Cmp tCmp1 = 0; tCmp1 < Space::Dimen(); ++tCmp1)
		for(Cmp tCmp2 = 0; tCmp2 < Space::Dimen(); ++tCmp2)
			tMatrix.mComps[tCmp1][tCmp2] = mComps[tCmp1][tCmp2] * theFactor;
	return tMatrix;
	CatchError
}



/*!
	Subtraction with a scalar.
*/
Tuple const Matrix::operator * (Tuple const & theTuple) const
{
	WatchError
	Tuple tTuple;
	for(Cmp tCmp1 = 0; tCmp1 < Space::Dimen(); ++tCmp1)
	{
		tTuple.Comp(tCmp1) = 0;
		for (Cmp tCmp2 = 0; tCmp2 < Space::Dimen(); ++tCmp2)
			tTuple.Comp(tCmp1) += mComps[tCmp1][tCmp2] * theTuple.Comp(tCmp2);
	}
	return tTuple;
	CatchError
}


closePlatypusSpace

