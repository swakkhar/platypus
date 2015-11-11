/*!
	@file realmodel/AffineDistObj.hpp
	@brief The prototype file for AffineDistObj class.
	@details This is the prototype file for AffineDistObj class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#ifndef AffineDistObjHhIncluded
#define AffineDistObjHhIncluded



#include "pspl/globals/function.hpp"


openPlatypusSpace



/*!
	@class AffineDistObj
	@brief A class to represent HP fitness constraint function for MJ
	@details This class represent HP fitness constraint function for MJ
*/
class AffineDistObj : public Function
{
	public:
		AffineDistObj(Conf const & theConf, B const NeedHint);						//!< Constructor.
		AffineDistObj(AffineDistObj const & that);						//!< Duplicator.
		AffineDistObj const & operator = (AffineDistObj const & that);	//!< Assignment.
		~AffineDistObj();											//!< Destructor.
};



closePlatypusSpace



#endif //AffineDistObjHhIncluded

