/*!
	@file HpObj.hpp
	@brief The prototype file for HpObj class.
	@details This is the prototype file for HpObj class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#ifndef HpObjHhIncluded
#define HpObjHhIncluded



#include "pspl/globals/function.hpp"

openPlatypusSpace



/*!
	@class HpObj
	@brief A class to represent HP energy objective function.
	@details This class represent HP energy objective function.
*/
class HpObj : public Function
{
	public:
		HpObj(Conf const & theConf,B const NeedHint);							//!< Constructor.
		HpObj(HpObj const & that);						//!< Duplicator.
		HpObj const & operator = (HpObj const & that);	//!< Assignment.
		~HpObj();											//!< Destructor.
};



closePlatypusSpace



#endif //HpObjHhIncluded


