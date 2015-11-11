
/*!
	@file moves/LNSMoveSC.hpp
	@brief The prototype file for LNSMoveSC class.
	@details This is the prototype file for LNSMoveSC class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 05.08.2012 QRL NICTA www.nicta.com.au
*/
#ifndef LNSMoveSCHppIncluded
#define LNSMoveSCHppIncluded


#include "pspl/globals/move.hpp"

openPlatypusSpace

class LNSMoveSC:public Move
{
    public:

		LNSMoveSC();
		~LNSMoveSC();
		LNSMoveSC(LNSMoveSC const & that);
		LNSMoveSC const & operator= (LNSMoveSC const & that);

        virtual void compute(Conf & theConf, Pos const thePos);

        void addPosition(Pos tPos);
        void reset();
        bool alreadyAdded(Pos tPos);

    private:
        block1<Pos,xmm> aPositions; // it holds the points that are to be altered
        bool next(block1<Dir,kmm> &generator);
        bool skip(block1<Dir,kmm> &generator,Idx pos);
        bool equal(block1<Dir,kmm> &generator, block1<Dir,kmm> &prior);

};



closePlatypusSpace

#endif




