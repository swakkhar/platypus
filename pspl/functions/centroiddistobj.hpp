
/*!
	@file functions/centroiddistobj.hpp
	@brief The prototype file for CentroidDistObj class.
	@details This is the prototype file for CentroidDistObj class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#ifndef CentroidDistObjHhIncluded
#define CentroidDistObjHhIncluded



#include "pspl/globals/function.hpp"
#include "pspl/energys/hpenergy.hpp"
#include "pspl/energys/mjenergy.hpp"
#include "pspl/energys/barerraenergy.hpp"


openPlatypusSpace



/*!
	@class CentroidDistObj
	@brief A class to represent distance to the centroid.
	@details This class represent distance to the centroid.
*/
class CentroidDistObj : public Function
{
	public:
		CentroidDistObj(Conf const & theConf,B const NeedHint);		//!< Constructor.
		CentroidDistObj(CentroidDistObj const & that);						//!< Duplicator.
		CentroidDistObj const & operator = (CentroidDistObj const & that);	//!< Assignment.
		~CentroidDistObj();											//!< Destructor.
};



closePlatypusSpace



#endif //CentroidDistObjHhIncluded

