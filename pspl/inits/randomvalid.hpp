/*!
	@file moves/randomvalid.hpp
	@brief The prototype file for RandomValid class.
	@details This is the prototype file for RandomValid class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/


#ifndef RandomValidHppIncluded
#define RandomValidHppIncluded



#include "pspl/globals/init.hpp"



openPlatypusSpace



/*!
	@class RandomValid
	@brief A class to represent random-valid inits.
	@details This class represents random-valid inits.
*/
class RandomValid: public Init
{
    public:
		RandomValid(Rnd & theRnd);										//!< Constructor.
		~RandomValid();													//!< Destructor.
		RandomValid(RandomValid const & that);							//!< Duplicator.
		RandomValid const & operator= (RandomValid const & that);		//!< Assigner.

        virtual void compute(Conf & theConf);							//!< Compute.

	private:
		Rnd & mRnd;														//!< Random number generator.
};



closePlatypusSpace

#endif
