
/*!
	@file functions/contactorderobj.hpp
	@brief The prototype file for ContactOrderObj class.
	@details This is the prototype file for ContactOrderObj class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#ifndef ContactOrderObjHhIncluded
#define ContactOrderObjHhIncluded

#include "pspl/globals/function.hpp"
#include "pspl/energys/hpenergy.hpp"

#include "pspl/energys/barerraenergy.hpp"
#include "pspl/energys/mjenergy.hpp"


openPlatypusSpace



/*!
	@class ContactOrderObj
	@brief A class to represent distance to the centroid.
	@details This class represent distance to the centroid.
*/
class ContactOrderObj : public Function
{
	public:
		ContactOrderObj(Conf const & theConf,B const NeedHint);		//!< Constructor.
		ContactOrderObj(ContactOrderObj const & that);						//!< Duplicator.
		ContactOrderObj const & operator = (ContactOrderObj const & that);	//!< Assignment.
		~ContactOrderObj();											//!< Destructor.
};



closePlatypusSpace



#endif //ContactOrderObjHhIncluded


