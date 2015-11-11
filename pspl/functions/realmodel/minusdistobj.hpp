/*!
	@file mj/mjplusdistobj.hpp
	@brief The prototype file for MinusDistObj class.
	@details This is the prototype file for MinusDistObj class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#ifndef MinusDistObjHhIncluded
#define MinusDistObjHhIncluded



#include "pspl/globals/function.hpp"
#include "pspl/energys/mjenergy.hpp"
#include "pspl/energys/barerraenergy.hpp"


openPlatypusSpace



/*!
	@class MinusDistObj
	@brief A class to represent HP fitness constraint function.
	@details This class represent HP fitness constraint function.
*/
class MinusDistObj : public Function
{
	public:
		MinusDistObj(Conf const & theConf, B const NeedHint);						//!< Constructor.
		MinusDistObj(MinusDistObj const & that);						//!< Duplicator.
		MinusDistObj const & operator = (MinusDistObj const & that);	//!< Assignment.
		~MinusDistObj();											//!< Destructor.
};



closePlatypusSpace



#endif //MinusDistObjHhIncluded


