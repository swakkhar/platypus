
/*!
	@file hpmodel/hcoredistobj.hpp
	@brief The prototype file for HCoreDistObj class.
	@details This is the prototype file for HCoreDistObj class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#ifndef HCoreDistObjHhIncluded
#define HCoreDistObjHhIncluded



#include "pspl/globals/function.hpp"
#include "pspl/energys/hpenergy.hpp"


openPlatypusSpace



/*!
	@class HCoreDistObj
	@brief A class to represent HP fitness constraint function.
	@details This class represent HP fitness constraint function.
*/
class HCoreDistObj : public Function
{
	public:
		HCoreDistObj(Conf const & theConf, B const NeedHint);						//!< Constructor.
		HCoreDistObj(HCoreDistObj const & that);						//!< Duplicator.
		HCoreDistObj const & operator = (HCoreDistObj const & that);	//!< Assignment.
		~HCoreDistObj();											//!< Destructor.
};



closePlatypusSpace



#endif //HCoreDistObjHhIncluded
