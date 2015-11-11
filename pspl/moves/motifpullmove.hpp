
/*!
	@file moves/MotifPullMove.hpp
	@brief The prototype file for MotifPullMove class.
	@details This is the prototype file for MotifPullMove class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/
#ifndef MotifPullHppIncluded
#define MotifPullHppIncluded


#include "pspl/globals/move.hpp"

openPlatypusSpace

class MotifPullMove:public Move
{
    public:

		MotifPullMove();
		~MotifPullMove();
		MotifPullMove(MotifPullMove const & that);
		MotifPullMove const & operator= (MotifPullMove const & that);

        virtual void compute(Conf & theConf, Pos const thePos);
};



closePlatypusSpace

#endif
