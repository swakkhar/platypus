/*!
	@file moves/chaingrowthinit.hpp
	@brief The prototype file for ChainGrowth class.
	@details This is the prototype file for ChainGrowth class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/


#ifndef ChainGrowthHppIncluded
#define ChainGrowthHppIncluded



#include "pspl/globals/init.hpp"



openPlatypusSpace



/*!
	@class ChainGrowth
	@brief A class to represent random-sphere inits.
	@details This class represents random-sphere inits.
*/
class ChainGrowth: public Init
{
    public:
		ChainGrowth(Rnd & theRnd);										//!< Constructor.
		~ChainGrowth();													//!< Destructor.
		ChainGrowth(ChainGrowth const & that);							//!< Duplicator.
		ChainGrowth const & operator= (ChainGrowth const & that);		//!< Assigner.

        virtual void compute(Conf & theConf);							//!< Compute.

	private:
		Rnd & mRnd;														//!< Random number generator.
};



closePlatypusSpace

#endif


