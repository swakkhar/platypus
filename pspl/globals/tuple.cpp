/*!
	@file globals/tuple.cpp
	@brief The implementation file for Tuple class.
	@details This is the implementation file for Tuple class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/


#include "pspl/globals/tuple.hpp"
#include "pspl/globals/matrix.hpp"



openPlatypusSpace



/*!
	The initialiser.
*/
Tuple::Tuple(Spc const theSpc)
{
	WatchError
	mComps = new Spc [Space::Dimen()];
	for (Cmp tCmp = 0; tCmp < Space::Dimen(); ++tCmp)
		mComps[tCmp] = theSpc;
	CatchError
}




/*!
	The duplicator.
*/
Tuple::Tuple(Tuple const & that)
{
	WatchError
	mComps = new Spc [Space::Dimen()];
	for (Cmp tCmp = 0; tCmp < Space::Dimen(); ++tCmp)
		mComps[tCmp] = that.mComps[tCmp];
	CatchError
}



/*!
	The assigner.
*/
Tuple const & Tuple::operator = (Tuple const & that)
{
	WatchError
	if (&that != this)
		for (Cmp tCmp = 0; tCmp < Space::Dimen(); ++tCmp)
			mComps[tCmp] = that.mComps[tCmp];
	return *this;
	CatchError
}



/*!
	The assigner.
*/
Tuple const & Tuple::operator = (Spc const theSpc)
{
	WatchError
	for (Cmp tCmp = 0; tCmp < Space::Dimen(); ++tCmp)
		mComps[tCmp] = theSpc;
	return *this;
	CatchError
}



/*!
	Check equal.
*/
B Tuple::operator == (Tuple const & that) const
{
	WatchError
	if (&that == this) return true;
	for (Cmp tCmp = 0; tCmp < Space::Dimen(); ++tCmp)
		if (mComps[tCmp] != that.mComps[tCmp])
			return false;
	return true;
	CatchError
}



/*!
	Check not-equal.
*/
B Tuple::operator != (Tuple const & that) const
{
	WatchError
	if (&that == this) return false;
	for (Cmp tCmp = 0; tCmp < Space::Dimen(); ++tCmp)
		if (mComps[tCmp] != that.mComps[tCmp])
			return true;
	return false;
	CatchError
}



/*!
	Check less-than.
*/
B Tuple::operator < (Tuple const & that) const
{
	WatchError
	if (&that == this) return false;
	for (Cmp tCmp = 0; tCmp < Space::Dimen(); ++tCmp)
		if (mComps[tCmp] < that.mComps[tCmp])
			return true;
		else if (mComps[tCmp] > that.mComps[tCmp])
			return false;
	return false;
	CatchError
}



/*!
	Check greater-than.
*/
B Tuple::operator > (Tuple const & that) const
{
	WatchError
	if (&that == this) return false;
	for (Cmp tCmp = 0; tCmp < Space::Dimen(); ++tCmp)
		if (mComps[tCmp] > that.mComps[tCmp])
			return true;
		else if (mComps[tCmp] < that.mComps[tCmp])
			return false;
	return false;
	CatchError
}


/*!
	Check less-equal.
*/
B Tuple::operator <= (Tuple const & that) const
{
	WatchError
	if (&that == this) return false;
	for (Cmp tCmp = 0; tCmp < Space::Dimen(); ++tCmp)
		if (mComps[tCmp] < that.mComps[tCmp])
			return true;
		else if (mComps[tCmp] > that.mComps[tCmp])
			return false;
	return true;
	CatchError
}



/*!
	Check greater-equal.
*/
B Tuple::operator >= (Tuple const & that) const
{
	WatchError
	if (&that == this) return false;
	for (Cmp tCmp = 0; tCmp < Space::Dimen(); ++tCmp)
		if (mComps[tCmp] > that.mComps[tCmp])
			return true;
		else if (mComps[tCmp] < that.mComps[tCmp])
			return false;
	return true;
	CatchError
}



/*!
	Check plus-assignment.
*/
Tuple const & Tuple::operator += (Tuple const & that)
{
	WatchError
	for(Cmp tCmp = 0; tCmp < Space::Dimen(); ++tCmp)
		mComps[tCmp] += that.mComps[tCmp];
	return *this;
	CatchError
}



/*!
	Check minus-assignment.
*/
Tuple const & Tuple::operator -= (Tuple const & that)
{
	WatchError
	for(Cmp tCmp = 0; tCmp < Space::Dimen(); ++tCmp)
		mComps[tCmp] -= that.mComps[tCmp];
	return *this;
	CatchError
}



/*!
	Check multiplication-assignment.
*/
Tuple const & Tuple::operator *= (Tuple const & that)
{
	WatchError
	for(Cmp tCmp = 0; tCmp < Space::Dimen(); ++tCmp)
		mComps[tCmp] *= that.mComps[tCmp];
	return *this;
	CatchError
}



/*!
	Check multiplication-assignment.
*/
Tuple const & Tuple::operator *= (Matrix const & theMatrix)
{
	WatchError
	Tuple tTuple(*this);

	for(Cmp tCmp1 = 0; tCmp1 < Space::Dimen(); ++tCmp1)
	{
		mComps[tCmp1] = 0;
		for(Cmp tCmp2 = 0; tCmp2 < Space::Dimen(); ++tCmp2)
			mComps[tCmp1] += tTuple.mComps[tCmp2] * theMatrix.Comp(tCmp2, tCmp1);
	}
	return *this;
	CatchError
}



/*!
	Check multiplication-assignment.
*/
Tuple const & Tuple::operator *= (Spc const & theFactor)
{
	WatchError
	for(Cmp tCmp = 0; tCmp < Space::Dimen(); ++tCmp)
		mComps[tCmp] *= theFactor;
	return *this;
	CatchError
}



/*!
	Check plus-assignment.
*/
Tuple const & Tuple::operator += (Spc const & theBias)
{
	WatchError
	for(Cmp tCmp = 0; tCmp < Space::Dimen(); ++tCmp)
		mComps[tCmp] += theBias;
	return *this;
	CatchError
}



/*!
	Check minus-assignment.
*/
Tuple const & Tuple::operator -= (Spc const & theBias)
{
	WatchError
	for(Cmp tCmp = 0; tCmp < Space::Dimen(); ++tCmp)
		mComps[tCmp] -= theBias;
	return *this;
	CatchError
}



/*!
	Negation.
*/
Tuple Tuple::operator - () const
{
	WatchError
	Tuple tTuple;
	for (Cmp tCmp = 0; tCmp < Space::Dimen(); ++tCmp)
		tTuple.mComps[tCmp] = -mComps[tCmp];
	return tTuple;
	CatchError
}



/*!
	Addition.
*/
Tuple Tuple::operator + (Tuple const & that) const
{
	WatchError
	Tuple tTuple;
	for(Cmp tCmp = 0; tCmp < Space::Dimen(); ++tCmp)
		tTuple.mComps[tCmp] = mComps[tCmp] + that.mComps[tCmp];
	return tTuple;
	CatchError
}



/*!
	Subtraction.
*/
Tuple Tuple::operator - (Tuple const & that) const
{
	WatchError
	Tuple tTuple;
	for(Cmp tCmp = 0; tCmp < Space::Dimen(); ++tCmp)
		tTuple.mComps[tCmp] = mComps[tCmp] - that.mComps[tCmp];
	return tTuple;
	CatchError
}



/*!
	Multiplication.
*/
Tuple Tuple::operator * (Tuple const & that) const
{
	WatchError
	Tuple tTuple;
	for(Cmp tCmp = 0; tCmp < Space::Dimen(); ++tCmp)
		tTuple.mComps[tCmp] = mComps[tCmp] * that.mComps[tCmp];
	return tTuple;
	CatchError
}



/*!
	Multiplication.
*/
Tuple Tuple::operator * (Matrix const & theMatrix) const
{
	WatchError
	Tuple tTuple;

	for(Cmp tCmp1 = 0; tCmp1 < Space::Dimen(); ++tCmp1)
	{
		tTuple.mComps[tCmp1] = 0;
		for(Cmp tCmp2 = 0; tCmp2 < Space::Dimen(); ++tCmp2)
			tTuple.mComps[tCmp1] += mComps[tCmp2] * theMatrix.Comp(tCmp2, tCmp1);
	}
	return tTuple;
	CatchError
}



/*!
	Addition.
*/
Tuple Tuple::operator + (Spc const & theBias) const
{
	WatchError
	Tuple tTuple;
	for(Cmp tCmp = 0; tCmp < Space::Dimen(); ++tCmp)
		tTuple.mComps[tCmp] = mComps[tCmp] + theBias;
	return tTuple;
	CatchError
}



/*!
	Subtraction.
*/
Tuple Tuple::operator - (Spc const & theBias) const
{
	WatchError
	Tuple tTuple;
	for(Cmp tCmp = 0; tCmp < Space::Dimen(); ++tCmp)
		tTuple.mComps[tCmp] = mComps[tCmp] - theBias;
	return tTuple;
	CatchError
}



/*!
	Multiplication.
*/
Tuple Tuple::operator * (Spc const & theFactor) const
{
	WatchError
	Tuple tTuple;
	for(Cmp tCmp = 0; tCmp < Space::Dimen(); ++tCmp)
		tTuple.mComps[tCmp] = mComps[tCmp] * theFactor;
	return tTuple;
	CatchError
}



/*!
	Return length.
*/
Dst Tuple::length() const
{
	WatchError
	Dst tDst = 0;
	for(Cmp tCmp = 0; tCmp < Space::Dimen(); ++tCmp)
		tDst += mComps[tCmp] * mComps[tCmp];
	return sqrt(tDst);
	CatchError
}


/*!
	Unit tuple.
*/
Tuple Tuple::unit() const
{
	WatchError
	Dst tDst = 0; Tuple tTuple(*this);
	for(Cmp tCmp = 0; tCmp < Space::Dimen(); ++tCmp)
		tDst += mComps[tCmp] * mComps[tCmp];
	tDst = sqrt(tDst);
	for(Cmp tCmp = 0; tCmp < Space::Dimen(); ++tCmp)
		tTuple.mComps[tCmp] /= tDst;
	return tTuple;
	CatchError
}



/*!
	Return dot product.
*/
Spc Tuple::dot(Tuple const & Tuple1, Tuple const & Tuple2)
{
	WatchError
	Spc tSpc = 0;
	for(Cmp tCmp = 0; tCmp < Space::Dimen(); ++tCmp)
		tSpc += Tuple1.mComps[tCmp] * Tuple2.mComps[tCmp];
	return tSpc;
	CatchError
}




/*!
	Return cross product.
*/
Tuple Tuple::cross(Tuple const & Tuple1, Tuple const & Tuple2)
{
	WatchError
	Warn(Space::Dimen() != 3, eCrossProduct);

	Tuple tTuple;
	tTuple.mComps[0] = Tuple1.mComps[1] * Tuple2.mComps[2] - Tuple1.mComps[2] * Tuple2.mComps[1];
	tTuple.mComps[1] = Tuple1.mComps[2] * Tuple2.mComps[0] - Tuple1.mComps[0] * Tuple2.mComps[2];
	tTuple.mComps[2] = Tuple1.mComps[0] * Tuple2.mComps[1] - Tuple1.mComps[1] * Tuple2.mComps[0];
	return tTuple;
	CatchError
}


/*!
	Return manhattan distance.
*/
Dst Tuple::manh(Tuple const & Tuple1, Tuple const & Tuple2)
{
	WatchError
	Dst tDst = 0;
	for(Cmp tCmp = 0; tCmp < Space::Dimen(); ++tCmp)
		tDst += abs(Tuple1.mComps[tCmp] - Tuple2.mComps[tCmp]);
	return tDst;
	CatchError
}


/*!
	Return cartesian distance.
*/
Dst Tuple::dist(Tuple const & Tuple1, Tuple const & Tuple2)
{
	WatchError
	Dst tDst = 0;
	//cerr << "cerr no,dist1:" << errno << endl;
	for(Cmp tCmp = 0; tCmp < Space::Dimen(); ++tCmp)
	{
		Spc tSpc = Tuple1.mComps[tCmp] - Tuple2.mComps[tCmp];

		tDst += tSpc * tSpc;
        //cout << tSpc * tSpc<<" " <<tDst <<endl;
	}
	//cerr << "cerr no,dist2:" << tDst << " "<< errno << endl;
	return sqrt(tDst);
   	CatchError
}


/*!
	Return square distance.
*/
Dst Tuple::sqrdist(Tuple const & Tuple1, Tuple const & Tuple2)
{
	WatchError
	Dst tDst = 0;
	for(Cmp tCmp = 0; tCmp < Space::Dimen(); ++tCmp)
	{
		Spc tSpc = Tuple1.mComps[tCmp] - Tuple2.mComps[tCmp];
		tDst += tSpc * tSpc;
	}
	return tDst;
	CatchError
}



closePlatypusSpace
