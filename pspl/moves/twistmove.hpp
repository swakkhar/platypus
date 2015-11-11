
/*!
	@file moves/TwistMove.hpp
	@brief The prototype file for TwistMove class.
	@details This is the prototype file for TwistMove class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/
#ifndef TwistHppIncluded
#define TwistHppIncluded


#include "pspl/globals/move.hpp"

openPlatypusSpace

class TwistMove:public Move
{
    public:

		TwistMove();
		~TwistMove();
		TwistMove(TwistMove const & that);
		TwistMove const & operator= (TwistMove const & that);

        virtual void compute(Conf & theConf, Pos const thePos);
};



closePlatypusSpace

#endif


