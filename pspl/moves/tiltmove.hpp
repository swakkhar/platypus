/*!
	@file moves/TiltMove.hpp
	@brief The prototype file for TiltMove class.
	@details This is the prototype file for TiltMove class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/
#ifndef TiltHppIncluded
#define TiltHppIncluded


#include "pspl/globals/move.hpp"

openPlatypusSpace

class TiltMove:public Move
{
    public:

		TiltMove();
		~TiltMove();
		TiltMove(TiltMove const & that);
		TiltMove const & operator= (TiltMove const & that);

        virtual void compute(Conf & theConf, Pos const thePos);

};



closePlatypusSpace

#endif


