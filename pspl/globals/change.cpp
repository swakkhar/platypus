
/*!
	@file globals/change.cpp
	@brief The implementation file for Change class.
	@details This is the implementation file for Change class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#include "pspl/globals/change.hpp"



openPlatypusSpace

/*B Change::operator == (Change const & that) const
{
    if(fVars.size()!=that.sizeFVal()) return false;
    else if(fVars.size())
    {
        for(Idx tIdx=0;tIdx < fVars.size()-1;++tIdx)
        {
            if(fVars[tIdx]!=that.getFVal(tIdx)) return false;
        }
        return true;
    }
    else return false;
}
B Change::operator > (Change const & that) const
{
    if(fVars.size()!=that.sizeFVal()) return false;
    else if(fVars.size())
    {
        for(Idx tIdx=0;tIdx < fVars.size()-1;++tIdx)
        {
            if(fVars[tIdx]<=that.getFVal(tIdx)) return false;
        }
        return true;
    }
    else return false;
}
B Change::operator < (Change const & that) const
{
    if(fVars.size()!=that.sizeFVal()) return false;
    else if(fVars.size())
    {
        for(Idx tIdx=0;tIdx < fVars.size()-1;++tIdx)
        {
            if(fVars[tIdx]>=that.getFVal(tIdx)) return false;
        }
        return true;
    }
    else return false;
}*/

closePlatypusSpace
