/*!
	@file moves/BlockInit.hpp
	@brief The prototype file for BlockInit class.
	@details This is the prototype file for BlockInit class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/


#ifndef BlockInitHppIncluded
#define BlockInitHppIncluded



#include "pspl/globals/init.hpp"



openPlatypusSpace



/*!
	@class BlockInit
	@brief A class to represent BlockInit inits.
	@details This class represents BlockInit inits.
*/
class BlockInit: public Init
{
    public:
		BlockInit(block1<Point,kmm> &arr);										//!< Constructor.
		~BlockInit();													//!< Destructor.
		BlockInit(BlockInit const & that);							//!< Duplicator.
		BlockInit const & operator= (BlockInit const & that);		//!< Assigner.

        virtual void compute(Conf & theConf);							//!< Compute.

	private:
		block1<Point,kmm> &mPoints;														//!< block of points.
};



closePlatypusSpace

#endif

