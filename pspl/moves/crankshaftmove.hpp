
/*!
	@file moves/CrankShaftMove.hpp
	@brief The prototype file for CrankShaftMove class.
	@details This is the prototype file for CrankShaftMove class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/
#ifndef CrankShaftHppIncluded
#define CrankShaftHppIncluded


#include "pspl/globals/move.hpp"

openPlatypusSpace

class CrankShaftMove:public Move
{
    public:

		CrankShaftMove();
		~CrankShaftMove();
		CrankShaftMove(CrankShaftMove const & that);
		CrankShaftMove const & operator= (CrankShaftMove const & that);

        virtual void compute(Conf & theConf, Pos const thePos);

};



closePlatypusSpace

#endif


