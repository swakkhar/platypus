/*!
	@file HPNXObj.hpp
	@brief The prototype file for HPNXObj class.
	@details This is the prototype file for HPNXObj class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#ifndef HPNXObjHhIncluded
#define HPNXObjHhIncluded



#include "pspl/globals/function.hpp"

openPlatypusSpace



/*!
	@class HPNXObj
	@brief A class to represent HP energy objective function.
	@details This class represent HP energy objective function.
*/
class HPNXObj : public Function
{
	public:
		HPNXObj(Conf const & theConf,B const NeedHint);							//!< Constructor.
		HPNXObj(HPNXObj const & that);						//!< Duplicator.
		HPNXObj const & operator = (HPNXObj const & that);	//!< Assignment.
		~HPNXObj();											//!< Destructor.
};



closePlatypusSpace



#endif //HPNXObjHhIncluded



