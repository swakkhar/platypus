/*!
	@file functions/FreeFtCns.hpp
	@brief The prototype file for FreeFtCns class.
	@details This is the prototype file for FreeFtCns class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#ifndef FreeFtCnsHhIncluded
#define FreeFtCnsHhIncluded



#include "pspl/globals/function.hpp"


openPlatypusSpace



/*!
	@class FreeFtCns
	@brief A class to represent connectedness constraint.
	@details This class represent connectedness constraint.
*/
class FreeFtCns : public Function
{
	public:
		FreeFtCns(Conf const & theConf, B const NeedHint);		//!< Constructor.
		FreeFtCns(FreeFtCns const & that);					//!< Duplicator.
		FreeFtCns const & operator = (FreeFtCns const & that);//!< Assignment.
		~FreeFtCns();											//!< Destructor.
};



closePlatypusSpace



#endif //FreeFtCnsHhIncluded

