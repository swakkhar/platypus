
/*!
	@file moves/RotationMove.hpp
	@brief The prototype file for RotationMove class.
	@details This is the prototype file for RotationMove class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 05.08.2012 QRL NICTA www.nicta.com.au
*/
#ifndef RotationMoveHppIncluded
#define RotationMoveHppIncluded


#include "pspl/globals/move.hpp"

openPlatypusSpace

class RotationMove:public Move
{
    public:

		RotationMove();
		~RotationMove();
		RotationMove(RotationMove const & that);
		RotationMove const & operator= (RotationMove const & that);

        virtual void compute(Conf & theConf, Pos const thePos);

        /*void addPosition(Pos tPos);
        void reset();
        bool alreadyAdded(Pos tPos);*/

    private:
        /*block1<Pos,xmm> aPositions; // it holds the points that are to be altered
        bool next(block1<Dir,kmm> &generator);
        bool skip(block1<Dir,kmm> &generator,Idx pos);
        bool equal(block1<Dir,kmm> &generator, block1<Dir,kmm> &prior);*/

};



closePlatypusSpace

#endif



