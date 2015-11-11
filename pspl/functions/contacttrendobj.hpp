

/*!
	@file functions/ContactTrendObj.hpp
	@brief The prototype file for ContactTrendObj class.
	@details This is the prototype file for ContactTrendObj class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#ifndef ContactTrendObjHhIncluded
#define ContactTrendObjHhIncluded



#include "pspl/globals/function.hpp"
#include "pspl/energys/hpenergy.hpp"
#include "pspl/energys/mjenergy.hpp"
#include "pspl/energys/barerraenergy.hpp"


openPlatypusSpace



/*!
	@class ContactTrendObj
	@brief A class to represent distance to the centroid.
	@details This class represent distance to the centroid.
*/
class ContactTrendObj : public Function
{
	public:
		ContactTrendObj(Conf const & theConf,B const NeedHint);		//!< Constructor.
		ContactTrendObj(ContactTrendObj const & that);						//!< Duplicator.
		ContactTrendObj const & operator = (ContactTrendObj const & that);	//!< Assignment.
		~ContactTrendObj();											//!< Destructor.
};

//const Dim mVal[] = {0,4,80,446,1442,3466,6917,12143,19422,28941,40800,133398,264754,406959,539329,651745,741779,811147,863217,901609,929569};
const Dim mVal[]=  {0,0,1,4,14,35,69,121,194,289,408,1334,2648,4070,5393,6517,7418,8111,8632,9016,9296,9498,9642,9746,9820,9872,9909,9936,9955};

const Dim cutOff[] =  {0,1,2,3,4,5,6,7,8,9,10,15,20,25,30,35,40,45,50,55,60};

const Dim cutOffCount = 21;
closePlatypusSpace



#endif //ContactTrendObjHhIncluded

