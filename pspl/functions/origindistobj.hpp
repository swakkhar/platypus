
/*!
	@file functions/origindistobj.hpp
	@brief The prototype file for OriginDistObj class.
	@details This is the prototype file for OriginDistObj class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#ifndef OriginDistObjHhIncluded
#define OriginDistObjHhIncluded



#include "pspl/globals/function.hpp"
#include "pspl/energys/hpenergy.hpp"


openPlatypusSpace



/*!
	@class OriginDistObj
	@brief A class to represent distance to the centroid.
	@details This class represent distance to the centroid.
*/
class OriginDistObj : public Function
{
	public:
		OriginDistObj(Conf const & theConf,B const NeedHint);		//!< Constructor.
		OriginDistObj(OriginDistObj const & that);						//!< Duplicator.
		OriginDistObj const & operator = (OriginDistObj const & that);	//!< Assignment.
		~OriginDistObj();											//!< Destructor.
};



closePlatypusSpace



#endif //OriginDistObjHhIncluded


