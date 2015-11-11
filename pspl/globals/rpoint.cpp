/*!
	@file globals/rpoint.cpp
	@brief The implementation file for RPoint class.
	@details This is the prototype file for RPoint class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/

#include "pspl/globals/rpoint.hpp"

openPlatypusSpace
/*!

default constructor

*/
RPoint::RPoint()
{

}

/*!

Constructor

*/
RPoint::RPoint(double a, double b, double c):x(a),y(b),z(c)
{

}

/*!

Destrcutor

*/
RPoint::~RPoint()
{

}

/*!

    Copy Constructor

*/
RPoint::RPoint(RPoint const & that)
{
    WatchError
        x=that.x;
        y=that.y;
        z=that.z;
    CatchError

}


RPoint const & RPoint::operator = (RPoint const & that)
{
    WatchError
	if (&that != this)
    {
        x=that.x;
        y=that.y;
        z=that.z;
    }
	return *this;
	CatchError
}

closePlatypusSpace
