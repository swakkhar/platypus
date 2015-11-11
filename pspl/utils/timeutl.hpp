/*!
	@file utils/timeutl.hpp
	@brief The prototype file for TimeUtl class.
	@details This is the prototype file for TimeUtl class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/


#ifndef TimeUtlHppIncluded
#define TimeUtlHppIncluded

#include "pspl/globals/data.hpp"
openPlatypusSpace

class TimeUtl
{
    public:
        TimeUtl();
        ~TimeUtl();
        R getTime();
};
closePlatypusSpace
#endif
