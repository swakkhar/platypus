/*!
	@file inits/TimeUtl.cpp
	@brief The implementation file for TimeUtl class.
	@details This is the implementation file for TimeUtl class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#include "pspl/utils/timeutl.hpp"

openPlatypusSpace

TimeUtl::TimeUtl()
{

}

TimeUtl::~TimeUtl()
{

}

R TimeUtl::getTime()
{
    struct tms sTime;
    times(&sTime);
    return castR(sTime.tms_utime + sTime.tms_stime + sTime.tms_cutime
				+ sTime.tms_cstime) / castR(sysconf(_SC_CLK_TCK));
}
closePlatypusSpace
