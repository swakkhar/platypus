
/*!
	@file functions/RMSDobj.hpp
	@brief The prototype file for RMSDObj class.
	@details This is the prototype file for RMSDObj class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#ifndef RMSDObjHhIncluded
#define RMSDObjHhIncluded



#include "pspl/globals/rpoint.hpp"
#include "pspl/globals/function.hpp"
#include "pspl/energys/hpenergy.hpp"
#include "pspl/energys/mjenergy.hpp"
#include "pspl/energys/barerraenergy.hpp"


openPlatypusSpace



/*!
	@class RMSDObj
	@brief A class to represent distance to the centroid.
	@details This class represent distance to the centroid.
*/
class RMSDObj : public Function
{
	public:
		RMSDObj(Conf const & theConf,B const NeedHint,block1<RPoint,xmm> &arr);		//!< Constructor.
		RMSDObj(RMSDObj const & that);						//!< Duplicator.
		RMSDObj const & operator = (RMSDObj const & that);	//!< Assignment.
		~RMSDObj();											//!< Destructor.
    private:
        block1<RPoint,xmm> &mNativePoints;
};



closePlatypusSpace



#endif //RMSDObjHhIncluded


