/*!
	@file moves/randomvalid.hpp
	@brief The prototype file for RandomSphere class.
	@details This is the prototype file for RandomSphere class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/


#ifndef RandomSphereHppIncluded
#define RandomSphereHppIncluded



#include "pspl/globals/init.hpp"



openPlatypusSpace



/*!
	@class RandomSphere
	@brief A class to represent random-sphere inits.
	@details This class represents random-sphere inits.
*/
class RandomSphere: public Init
{
    public:
		RandomSphere(Rnd & theRnd);										//!< Constructor.
		~RandomSphere();													//!< Destructor.
		RandomSphere(RandomSphere const & that);							//!< Duplicator.
		RandomSphere const & operator= (RandomSphere const & that);		//!< Assigner.

        virtual void compute(Conf & theConf);							//!< Compute.

	private:
		Rnd & mRnd;														//!< Random number generator.
};



closePlatypusSpace

#endif

