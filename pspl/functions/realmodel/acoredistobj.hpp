
/*!
	@file realmodel/acoredistobj.hpp
	@brief The prototype file for ACoreDistObj class.
	@details This is the prototype file for ACoreDistObj class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#ifndef ACoreDistObjHhIncluded
#define ACoreDistObjHhIncluded



#include "pspl/globals/function.hpp"


openPlatypusSpace



/*!
	@class ACoreDistObj
	@brief A class to represent HP fitness constraint function.
	@details This class represent HP fitness constraint function.
*/
class ACoreDistObj : public Function
{
	public:
		ACoreDistObj(Conf const & theConf, B const NeedHint);						//!< Constructor.
		ACoreDistObj(ACoreDistObj const & that);						//!< Duplicator.
		ACoreDistObj const & operator = (ACoreDistObj const & that);	//!< Assignment.
		~ACoreDistObj();											//!< Destructor.
};



closePlatypusSpace



#endif //ACoreDistObjHhIncluded

