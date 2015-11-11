/*!
	@file moves/TriplePullMove.hpp
	@brief The prototype file for TriplePullMove class.
	@details This is the prototype file for TriplePullMove class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/
#ifndef TriplePullHppIncluded
#define TriplePullHppIncluded


#include "pspl/globals/move.hpp"

openPlatypusSpace

class TriplePullMove:public Move
{
    public:

		TriplePullMove();
		~TriplePullMove();
		TriplePullMove(TriplePullMove const & that);
		TriplePullMove const & operator= (TriplePullMove const & that);

        virtual void compute(Conf & theConf, Pos const thePos);
};



closePlatypusSpace

#endif



