/*!
	@file functions/connectedcns.hpp
	@brief The prototype file for ConnectedCns class.
	@details This is the prototype file for ConnectedCns class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#ifndef ConnectedCnsHhIncluded
#define ConnectedCnsHhIncluded



#include "pspl/globals/function.hpp"


openPlatypusSpace



/*!
	@class ConnectedCns
	@brief A class to represent connectedness constraint.
	@details This class represent connectedness constraint.
*/
class ConnectedCns : public Function
{
	public:
		ConnectedCns(Conf const & theConf, B const NeedHint);		//!< Constructor.
		ConnectedCns(ConnectedCns const & that);					//!< Duplicator.
		ConnectedCns const & operator = (ConnectedCns const & that);//!< Assignment.
		~ConnectedCns();											//!< Destructor.
};



closePlatypusSpace



#endif //ConnectedCnsHhIncluded
