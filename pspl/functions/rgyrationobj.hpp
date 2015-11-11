
/*!
	@file functions/rgyrationobj.hpp
	@brief The prototype file for RGyrationObj class.
	@details This is the prototype file for RGyrationObj class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#ifndef RGyrationObjHhIncluded
#define RGyrationObjHhIncluded



#include "pspl/globals/function.hpp"
#include "pspl/energys/hpenergy.hpp"
#include "pspl/energys/mjenergy.hpp"
#include "pspl/energys/barerraenergy.hpp"


openPlatypusSpace



/*!
	@class RGyrationObj
	@brief A class to represent distance to the centroid.
	@details This class represent distance to the centroid.
*/
class RGyrationObj : public Function
{
	public:
		RGyrationObj(Conf const & theConf,B const NeedHint);		//!< Constructor.
		RGyrationObj(RGyrationObj const & that);						//!< Duplicator.
		RGyrationObj const & operator = (RGyrationObj const & that);	//!< Assignment.
		~RGyrationObj();											//!< Destructor.
};



closePlatypusSpace



#endif //RGyrationObjHhIncluded


