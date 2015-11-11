/*!
	@file mjobj.hpp
	@brief The prototype file for MjObj class.
	@details This is the prototype file for MjObj class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#ifndef MjObjHhIncluded
#define MjObjHhIncluded



#include "pspl/globals/function.hpp"
#include "pspl/energys/mjenergy.hpp"


openPlatypusSpace



/*!
	@class MjObj
	@brief A class to represent HP energy objective function.
	@details This class represent HP energy objective function.
*/
class MjObj : public Function
{
	public:
		MjObj(Conf const & theConf,B const NeedHint);							//!< Constructor.
		MjObj(MjObj const & that);						//!< Duplicator.
		MjObj const & operator = (MjObj const & that);	//!< Assignment.
		~MjObj();											//!< Destructor.
};



closePlatypusSpace



#endif //MjObjHhIncluded

