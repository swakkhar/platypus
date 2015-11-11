/*!
	@file moves/DoublePull.hpp
	@brief The prototype file for DoublePull class.
	@details This is the prototype file for DoublePull class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/
#ifndef DoublePullHppIncluded
#define DoublePullHppIncluded


#include "pspl/globals/move.hpp"

openPlatypusSpace

class DoublePullMove:public Move
{
    public:

		DoublePullMove();
		~DoublePullMove();
		DoublePullMove(DoublePullMove const & that);
		DoublePullMove const & operator= (DoublePullMove const & that);

        virtual void compute(Conf & theConf, Pos const thePos);

};



closePlatypusSpace

#endif

