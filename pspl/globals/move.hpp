/*!
	@file globals/move.hpp
	@brief The prototype file for Move class.
	@details This is the prototype file for Move class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/


#ifndef MoveHppIncluded
#define MoveHppIncluded


#include "pspl/globals/conf.hpp"
#include "pspl/globals/change.hpp"


openPlatypusSpace


/*!
	@class Move
	@brief A class represent moves.
	@details This class represents moves.
*/
class Move
{
	public:
		block1<Change,xmm> const & PartChanges() const;
		virtual void compute(Conf & theConf, Pos const thePos) = 0;

	protected:

		Move();
		virtual ~Move();
		Move(Move const & that);
		Move const & operator= (Move const & that);

		block1<Change,xmm> mChanges;
};


inline block1<Change,xmm> const & Move::PartChanges() const
{
	WatchError
	return mChanges;
	CatchError
}



closePlatypusSpace


#endif	//	MoveHppIncluded
