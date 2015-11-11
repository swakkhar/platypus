/*!
	@file functions/selfavoidcns.hpp
	@brief The prototype file for SelfAvoidCns class.
	@details This is the prototype file for SelfAvoidCns class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#ifndef SelfAvoidCnsHhIncluded
#define SelfAvoidCnsHhIncluded



#include "pspl/globals/function.hpp"


openPlatypusSpace



/*!
	@class SelfAvoidCns
	@brief A class to represent self avoiding constraint.
	@details This class represent self avoiding constraint.
*/
class SelfAvoidCns : public Function
{
	public:
		SelfAvoidCns(Conf const & theConf, B const NeedHint);		//!< Constructor.
		SelfAvoidCns(SelfAvoidCns const & that);					//!< Duplicator.
		SelfAvoidCns const & operator = (SelfAvoidCns const & that);//!< Assignment.
		~SelfAvoidCns();											//!< Destructor.
};



closePlatypusSpace



#endif //SelfAvoidCnsHhIncluded
