
/*!
	@file moves/LNSMove.hpp
	@brief The prototype file for LNSMove class.
	@details This is the prototype file for LNSMove class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 05.08.2012 QRL NICTA www.nicta.com.au
*/
#ifndef LNSMoveHppIncluded
#define LNSMoveHppIncluded


#include "pspl/globals/move.hpp"

openPlatypusSpace

class LNSMove:public Move
{
    public:

		LNSMove();
		~LNSMove();
		LNSMove(LNSMove const & that);
		LNSMove const & operator= (LNSMove const & that);

        virtual void compute(Conf & theConf, Pos const thePos);

        void addPosition(Pos tPos);
        void reset();
        bool alreadyAdded(Pos tPos);
        N size();

    private:
        block1<Pos,xmm> aPositions; // it holds the points that are to be altered
        bool next(block1<Dir,kmm> &generator);
        bool skip(block1<Dir,kmm> &generator,Idx pos);
        bool equal(block1<Dir,kmm> &generator, block1<Dir,kmm> &prior);

};



closePlatypusSpace

#endif



