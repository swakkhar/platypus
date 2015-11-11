/*!
	@file moves/PushMove.hpp
	@brief The prototype file for PushMove class.
	@details This is the prototype file for PushMove class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/
#ifndef PushHppIncluded
#define PushHppIncluded


#include "pspl/globals/move.hpp"

openPlatypusSpace

class PushMove:public Move
{
    public:

		PushMove();
		~PushMove();
		PushMove(PushMove const & that);
		PushMove const & operator= (PushMove const & that);

        virtual void compute(Conf & theConf, Pos const thePos);
};



closePlatypusSpace

#endif



