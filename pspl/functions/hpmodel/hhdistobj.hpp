/*!
	@file hpmodel/hhdistobj.hpp
	@brief The prototype file for HhDistObj class.
	@details This is the prototype file for HhDistObj class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#ifndef HhDistObjHhIncluded
#define HhDistObjHhIncluded



#include "pspl/globals/function.hpp"
#include "pspl/energys/hpenergy.hpp"


openPlatypusSpace



/*!
	@class HhDistObj
	@brief A class to represent HP fitness constraint function.
	@details This class represent HP fitness constraint function.
*/
class HhDistObj : public Function
{
	public:
		HhDistObj(Conf const & theConf, B const NeedHint);						//!< Constructor.
		HhDistObj(HhDistObj const & that);						//!< Duplicator.
		HhDistObj const & operator = (HhDistObj const & that);	//!< Assignment.
		~HhDistObj();											//!< Destructor.
};



closePlatypusSpace



#endif //HhDistObjHhIncluded
