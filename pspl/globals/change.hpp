/*!
	@file globals/change.hpp
	@brief The prototype file for Change class.
	@details This is the prototype file for Change class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#ifndef ChangeHppIncluded
#define ChangeHppIncluded


#include "pspl/globals/tuple.hpp"


openPlatypusSpace


/*!
	@class Change
	@brief A class to represent changes to be made in conformations.
	@details This is a class to represent changes to be made in conformations.
*/
class Change
{
	private:
		block1<tuple2<Pos, Point>,xmm > mPosDsts;				//!< Positions and destinations.
        //xblock<R> fVars;
	public:

		Change();											//!< Constructor.
		~Change();											//!< Destructor.
		Change(Change const & that);						//!< Duplicator.
		Change const & operator = (Change const & that);	//!< Assigner.

		Dim size() const;									//!< Change size.
		Pos Position(Idx const theIdx) const;				//!< Return a position.
		Point const & Destination(Idx const theIdx) const;	//!< Return a destination.

        Point const & DestinationByPos(Pos const tPos) const;

		void reset();										//!< Reset the positions and destinations.
		void remove();										//!< Remove the position-destination pair.
		void add(Pos const thePos, Point const & theDst);	//!< Add a position-destination pair.

        //void addFVal(R tVal); //! add fVal to the end of the lsit
        //R getFVal(Idx tIdx) const; //! get fVal from the position tIdx
        //N sizeFVal() const;     //! number of fVals
        //void resetFVals(); //! clear the list

        //B operator == (Change const & that) const ;
        //B operator > (Change const & that) const;
        //B operator < (Change const & that) const;
};



/*inline void Change::addFVal(R tVal)
{
    fVars.annex(tVal);
}

inline void Change::resetFVals()
{
    fVars.reset();
}
inline R Change::getFVal(Idx tIdx) const
{
    Alert(tIdx>=fVars.size(),eInvalidIndex);
    return fVars[tIdx];
}
inline N Change::sizeFVal() const
{
    return fVars.size();
}*/

/*!
	The constructor.
*/
inline Change::Change()
{
	WatchError
	//	nothing to be done.
	CatchError
}



/*!
	The destructor.
*/
inline Change::~Change()
{
	WatchError
	//	nothing to be done.
	CatchError
}



/*!
	The duplicator.
*/
inline Change::Change(Change const & that) : mPosDsts(that.mPosDsts)//,fVars(that.fVars)
{
	WatchError
	//	nothing to be done.
	CatchError
}



/*!
	The assigner.
*/
inline Change const & Change::operator=(Change const & that)
{
	WatchError
	if (this != &that)
    {
        mPosDsts = that.mPosDsts;
        //fVars = that.fVars;
	}

	return *this;
	CatchError
}


/*!
	Return the size of the change.
*/
inline Dim Change::size() const
{
	WatchError
	return mPosDsts.itemCount();
	CatchError
}



/*!
	Return a position in the change.
*/
inline Pos Change::Position(Idx const theIdx) const
{
	WatchError
	return mPosDsts[theIdx].First;
	CatchError
}



/*!
	Return a destination in the change.
*/
inline Point const & Change::Destination(Idx const theIdx) const
{
	WatchError
	return mPosDsts[theIdx].Second;
	CatchError
}

inline Point const & Change::DestinationByPos(Pos const tPos) const
{
    WatchError
    for(Idx tIdx=0;tIdx<mPosDsts.itemCount();++tIdx)
    {
        if(mPosDsts[tIdx].First==tPos) return mPosDsts[tIdx].Second;
    }
    cout<<"This is a invalid Destination call!" << endl;
    exit(0);
    CatchError
    }



/*!
	Reset the change.
*/
inline void Change::reset()
{
	WatchError
	mPosDsts.clear();
	CatchError
}



/*!
	Remove a position-point pair from the change.
*/
inline void Change::remove()
{
	WatchError
	mPosDsts.remove();
	CatchError
}



/*!
	Add a position-point pair to the change.
*/
inline void Change::add(Pos const thePos, Point const & theDst)
{
	WatchError
	mPosDsts.insertMem(tuple2<Pos,Point>(thePos, theDst));
	CatchError
}



closePlatypusSpace



#endif //ChangeHppIncluded

