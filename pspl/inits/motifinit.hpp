/*!
	@file moves/MotifInit.hpp
	@brief The prototype file for MotifInit class.
	@details This is the prototype file for MotifInit class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/


#ifndef MotifInitHppIncluded
#define MotifInitHppIncluded



#include "pspl/globals/init.hpp"



openPlatypusSpace



/*!
	@class MotifInit
	@brief A class to represent random-valid inits.
	@details This class represents random-valid inits.
*/
class MotifInit: public Init
{
    public:
		MotifInit(Rnd & theRnd);										//!< Constructor.
		~MotifInit();													//!< Destructor.
		MotifInit(MotifInit const & that);							//!< Duplicator.
		MotifInit const & operator= (MotifInit const & that);		//!< Assigner.

        virtual void compute(Conf & theConf);							//!< Compute.

	private:
		Rnd & mRnd;														//!< Random number generator.
};



closePlatypusSpace

#endif

