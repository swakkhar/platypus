
/*!
	@file functions/ContactCountObj.hpp
	@brief The prototype file for ContactCountObj class.
	@details This is the prototype file for ContactCountObj class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#ifndef ContactCountObjHhIncluded
#define ContactCountObjHhIncluded



#include "pspl/globals/function.hpp"
#include "pspl/energys/hpenergy.hpp"
#include "pspl/energys/mjenergy.hpp"
#include "pspl/energys/barerraenergy.hpp"


openPlatypusSpace



/*!
	@class ContactCountObj
	@brief A class to represent distance to the centroid.
	@details This class represent distance to the centroid.
*/
class ContactCountObj : public Function
{
	public:
		ContactCountObj(Conf const & theConf,B const NeedHint);		//!< Constructor.
		ContactCountObj(ContactCountObj const & that);						//!< Duplicator.
		ContactCountObj const & operator = (ContactCountObj const & that);	//!< Assignment.
		~ContactCountObj();											//!< Destructor.
};



closePlatypusSpace



#endif //ContactCountObjHhIncluded



