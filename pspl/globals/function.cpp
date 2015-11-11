/*!
	@file globals/function.cpp
	@brief The implementation file for Function class.
	@details This is the implementation file for Function class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/


#include "pspl/globals/function.hpp"



openPlatypusSpace



/*!
	The constructor.
*/
Function::Function(Hdl const SysHdl, B const NeedHint) :
	mSysHdl(SysHdl), mFuncHdl(InvHdl), mHaveHint(NeedHint), mMetricRec(Null)
{
	WatchError
	Warn(mSysHdl == InvHdl, eInvalidHandle);
	//	nothing to be done.
	CatchError
}



/*!
	The destructor.
*/
Function::~Function()
{
	WatchError
	//	nothing to be done.
	CatchError
}



/*!
	The assigner.
*/
Function::Function(const Function & that)
{
	WatchError
	Alert(&that, eUndefDuplicator);
	CatchError
}



/*!
	The assigner.
*/
Function const & Function::operator=(Function const & that)
{
	WatchError
	Alert(&that, eUndefAssigner);
	return *this;
	CatchError
}



closePlatypusSpace
