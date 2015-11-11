/*!
	@file functions/SideChainCns.hpp
	@brief The prototype file for SideChainCns class.
	@details This is the prototype file for SideChainCns class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#ifndef SideChainCnsHhIncluded
#define SideChainCnsHhIncluded



#include "pspl/globals/function.hpp"


openPlatypusSpace



/*!
	@class SideChainCns
	@brief A class to represent self avoiding constraint.
	@details This class represent self avoiding constraint.
*/
class SideChainCns : public Function
{
	public:
		SideChainCns(Conf const & theConf, B const NeedHint);		//!< Constructor.
		SideChainCns(SideChainCns const & that);					//!< Duplicator.
		SideChainCns const & operator = (SideChainCns const & that);//!< Assignment.
		~SideChainCns();											//!< Destructor.
};



closePlatypusSpace



#endif //SideChainCnsHhIncluded

