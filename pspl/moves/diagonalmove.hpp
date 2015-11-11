/*!
	@file moves/diagonalmove.hpp
	@brief The prototype file for DiagonalMove class.
	@details This is the prototype file for DiagonalMove class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/


#ifndef DiagonalMoveHppIncluded
#define DiagonalMoveHppIncluded



#include "pspl/globals/move.hpp"



openPlatypusSpace


/*!
	@class DiagonalMove
	@brief A class to represent diagonal moves.
	@details This class represents diagonal moves.
*/
class DiagonalMove: public Move
{
    public:
		DiagonalMove();													//!< Constructor.
		~DiagonalMove();												//!< Destructor.
		DiagonalMove(DiagonalMove const & that);						//!< Duplicator.
		DiagonalMove const & operator= (DiagonalMove const & that);		//!< Assigner.

        virtual void compute(Conf & theConf, Pos const thePos);			//!< Compute.

};



closePlatypusSpace

#endif
