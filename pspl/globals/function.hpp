/*!
	@file globals/function.hpp
	@brief The prototype file for Function class.
	@details This is the prototype file for Function class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/




#ifndef FunctionHppIncluded
#define FunctionHppIncluded



#include "pspl/globals/conf.hpp"



openPlatypusSpace



/*!
	@class Function
	@brief A class represent functions.
	@details This class represent functions.
*/
class Function
{
	protected:
		virtual ~Function();								//!< Destructor.
		Function(Function const & that);					//!< Duplicator.
		Function(Hdl const SysHdl, B const NeedHint);		//!< Constructor.
		Function const & operator=(Function const & that);	//!< Assigner.

	public:

		Hdl SysHdl() const;			//!< System handle.
		Hdl FuncHdl() const;		//!< Function handle.
		B HaveHint() const;			//!< Function has hints.
		Int ExecMetric() const;		//!< Execution metric.
		Int SimulMetric() const;	//!< Simulation data.

	protected:

		Hdl mSysHdl;					//!< The kangaroo system handle.
		Hdl mFuncHdl;					//!< The function handle.
		B  mHaveHint;					//!< The function has hints.
		EvalRecInt const * mMetricRec;	//!< The metric record.

};



/*!
	Return system handle of kangaroo.
*/
inline Hdl Function::SysHdl() const
{
	WatchError
	Warn(mSysHdl == InvHdl, eInvalidHandle);
	return mSysHdl;
	CatchError
}



/*!
	Return system handle of kangaroo.
*/
inline Hdl Function::FuncHdl() const
{
	WatchError
	Warn(mFuncHdl == InvHdl, eInvalidHandle);
	return mFuncHdl;
	CatchError
}



/*!
	Whether the function has fint.
*/
inline B Function::HaveHint() const
{
	WatchError
	return mHaveHint;
	CatchError
}



/*!
	Return execution metric.
*/
inline Int Function::ExecMetric() const
{
	WatchError
	Warn(!mMetricRec, eNullPointer);
	return mMetricRec->CurrData();
	CatchError
}



/*!
	Return simulation metric.
*/
inline Int Function::SimulMetric() const
{
	WatchError
	Warn(!mMetricRec, eNullPointer);
	return mMetricRec->NextData(Sys::refc(mSysHdl).SimulClk());
	CatchError
}



closePlatypusSpace



#endif	//	FunctionHppIncluded

