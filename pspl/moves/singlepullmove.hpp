/*!
	@file moves/SinglePullMove.hpp
	@brief The prototype file for SinglePullMove class.
	@details This is the prototype file for SinglePullMove class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/
#ifndef SinglePullHppIncluded
#define SinglePullHppIncluded


#include "pspl/globals/move.hpp"

openPlatypusSpace

class SinglePullMove:public Move
{
    public:

		SinglePullMove();
		~SinglePullMove();
		SinglePullMove(SinglePullMove const & that);
		SinglePullMove const & operator= (SinglePullMove const & that);

        virtual void compute(Conf & theConf, Pos const thePos);
};



closePlatypusSpace

#endif
