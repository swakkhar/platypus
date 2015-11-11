/*!
	@file moves/randomvalid.hpp
	@brief The prototype file for RandomValidSC class.
	@details This is the prototype file for RandomValidSC class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/


#ifndef RandomValidSCHppIncluded
#define RandomValidSCHppIncluded



#include "pspl/globals/init.hpp"



openPlatypusSpace



/*!
	@class RandomValidSC
	@brief A class to represent random-valid inits.
	@details This class represents random-valid inits.
*/
class RandomValidSC: public Init
{
    public:
		RandomValidSC(Rnd & theRnd);										//!< Constructor.
		~RandomValidSC();													//!< Destructor.
		RandomValidSC(RandomValidSC const & that);							//!< Duplicator.
		RandomValidSC const & operator= (RandomValidSC const & that);		//!< Assigner.

        virtual void compute(Conf & theConf);							//!< Compute.

	private:
		Rnd & mRnd;														//!< Random number generator.
};



closePlatypusSpace

#endif

