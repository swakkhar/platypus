
/*!
	@file realmodel/RepCntCns.hpp
	@brief The prototype file for RepCntCns class.
	@details This is the prototype file for RepCntCns class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#ifndef RepCntCnsHhIncluded
#define RepCntCnsHhIncluded



#include "pspl/globals/function.hpp"


openPlatypusSpace



/*!
	@class RepCntCns
	@brief A class to represent HP fitness constraint function.
	@details This class represent HP fitness constraint function.
*/
class RepCntCns : public Function
{
	public:
		RepCntCns(Conf const & theConf, B const NeedHint);						//!< Constructor.
		RepCntCns(RepCntCns const & that);						//!< Duplicator.
		RepCntCns const & operator = (RepCntCns const & that);	//!< Assignment.
		~RepCntCns();											//!< Destructor.
};



closePlatypusSpace



#endif //RepCntCnsHhIncluded


