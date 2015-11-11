/*!
	@file mj/mjplusdistobj.hpp
	@brief The prototype file for PlusDistObj class.
	@details This is the prototype file for PlusDistObj class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#ifndef PlusDistObjHhIncluded
#define PlusDistObjHhIncluded



#include "pspl/globals/function.hpp"
#include "pspl/energys/mjenergy.hpp"
#include "pspl/energys/barerraenergy.hpp"


openPlatypusSpace



/*!
	@class PlusDistObj
	@brief A class to represent HP fitness constraint function.
	@details This class represent HP fitness constraint function.
*/
class PlusDistObj : public Function
{
	public:
		PlusDistObj(Conf const & theConf, B const NeedHint);						//!< Constructor.
		PlusDistObj(PlusDistObj const & that);						//!< Duplicator.
		PlusDistObj const & operator = (PlusDistObj const & that);	//!< Assignment.
		~PlusDistObj();											//!< Destructor.
};



closePlatypusSpace



#endif //PlusDistObjHhIncluded

