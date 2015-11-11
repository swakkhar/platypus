/*!
	@file moves/PCLFInit.hpp
	@brief The prototype file for PCLFInit class.
	@details This is the prototype file for PCLFInit class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/


#ifndef PCLFInitHppIncluded
#define PCLFInitHppIncluded



#include "pspl/globals/init.hpp"
#include "pspl/globals/rpoint.hpp"



openPlatypusSpace



/*!
	@class PCLFInit
	@brief A class to represent PCLFInit inits.
	@details This class represents PCLFInit inits.
*/
class PCLFInit: public Init
{
    public:
		PCLFInit(block1<RPoint,xmm> &arr);										//!< Constructor.
		~PCLFInit();													//!< Destructor.
		PCLFInit(PCLFInit const & that);							//!< Duplicator.
		PCLFInit const & operator= (PCLFInit const & that);		//!< Assigner.

        virtual void compute(Conf & theConf);							//!< Compute.

	private:
		block1<RPoint,xmm> &mPoints;														//!< block of points.
};



closePlatypusSpace

#endif


