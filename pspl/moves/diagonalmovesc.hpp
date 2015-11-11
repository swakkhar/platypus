/*!
	@file moves/diagonalmove.hpp
	@brief The prototype file for DiagonalMoveSC class.
	@details This is the prototype file for DiagonalMoveSC class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/


#ifndef DiagonalMoveSCHppIncluded
#define DiagonalMoveSCHppIncluded



#include "pspl/globals/move.hpp"



openPlatypusSpace


/*!
	@class DiagonalMoveSC
	@brief A class to represent diagonal moves.
	@details This class represents diagonal moves.
*/
class DiagonalMoveSC: public Move
{
    public:
		DiagonalMoveSC();													//!< Constructor.
		~DiagonalMoveSC();												//!< Destructor.
		DiagonalMoveSC(DiagonalMoveSC const & that);						//!< Duplicator.
		DiagonalMoveSC const & operator= (DiagonalMoveSC const & that);		//!< Assigner.

        virtual void compute(Conf & theConf, Pos const thePos);			//!< Compute.

};



closePlatypusSpace

#endif

