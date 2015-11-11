/*!
	@file globals/rpoint.hpp
	@brief The prototype file for RPoint class.
	@details This is the prototype file for RPoint class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#ifndef RPointHppIncluded
#define RPointHppIncluded

#include "pspl/globals/data.hpp"

openPlatypusSpace
class RPoint
{

    public:
        double x;
        double y;
        double z;
        RPoint();
        RPoint(double a, double b, double c);
        ~RPoint();
        RPoint(RPoint const & that);
        RPoint const & operator = (RPoint const & that);
};
closePlatypusSpace
#endif
