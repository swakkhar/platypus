/*!
	@file moves/randomvalid.hpp
	@brief The prototype file for RandomStructured class.
	@details This is the prototype file for RandomStructured class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/


#ifndef RandomStructuredHppIncluded
#define RandomStructuredHppIncluded



#include "pspl/globals/init.hpp"



openPlatypusSpace



/*!
	@class RandomStructured
	@brief A class to represent random-structured inits.
	@details This class represents random-structured inits.
*/

const int transformMoves[]={-1,0,9,6,3,2,11,8,5,1,10,7,4};
class RandomStructured: public Init
{
    public:
		RandomStructured(Rnd & theRnd);										//!< Constructor.
		~RandomStructured();													//!< Destructor.
		RandomStructured(RandomStructured const & that);							//!< Duplicator.
		RandomStructured const & operator= (RandomStructured const & that);		//!< Assigner.

        virtual void compute(Conf & theConf);							//!< Compute.

	private:
		Rnd & mRnd;														//!< Random number generator.
        void updateForbiden(Z *f, Z m);

};



closePlatypusSpace

#endif

