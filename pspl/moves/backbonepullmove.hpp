/*!
	@file moves/BackBonePullMove.hpp
	@brief The prototype file for BackBonePullMove class.
	@details This is the prototype file for BackBonePullMove class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/
#ifndef BakBonePullHppIncluded
#define BakBonePullHppIncluded


#include "pspl/globals/move.hpp"

openPlatypusSpace

class BackBonePullMove:public Move
{
    public:

		BackBonePullMove();
		~BackBonePullMove();
		BackBonePullMove(BackBonePullMove const & that);
		BackBonePullMove const & operator= (BackBonePullMove const & that);

        virtual void compute(Conf & theConf, Pos const thePos);
};



closePlatypusSpace

#endif

