/*!
	@file hpmodel/hhneighcns.hpp
	@brief The prototype file for HhNeighCns class.
	@details This is the prototype file for HhNeighCns class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#ifndef HhNeighCnsHhIncluded
#define HhNeighCnsHhIncluded



#include "pspl/globals/function.hpp"
#include "pspl/energys/hpenergy.hpp"


openPlatypusSpace



/*!
	@class HhNeighCns
	@brief A class to represent HP neighbour constraint.
	@details This class represent HP neighbour constraint.
*/
class HhNeighCns : public Function
{
	public:
		HhNeighCns(Conf const & theConf, B const NeedHint);		//!< Constructor.
		HhNeighCns(HhNeighCns const & that);					//!< Duplicator.
		HhNeighCns const & operator = (HhNeighCns const & that);//!< Assignment.
		~HhNeighCns();											//!< Destructor.
};



closePlatypusSpace



#endif //HhNeighCnsHhIncluded
