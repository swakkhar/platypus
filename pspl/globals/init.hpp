/*!
	@file globals/init.hpp
	@brief The prototype file for Init class.
	@details This is the prototype file for Init class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/


#ifndef InitHppIncluded
#define InitHppIncluded


#include "pspl/globals/conf.hpp"
#include "pspl/globals/change.hpp"


openPlatypusSpace


/*!
	@class Init
	@brief A class represent initialisation.
	@details This class represents initiasation.
*/
class Init
{
    protected:
		Change mChange;

	public:

		Change const & FullChange() const;
		virtual void compute(Conf & theConf) = 0;

	protected:

		Init();
		virtual ~Init();
		Init(Init const & that);
		Init const & operator= (Init const & that);
};


inline Change const & Init::FullChange() const
{
	WatchError
	return mChange;
	CatchError
}


closePlatypusSpace


#endif	//	InitHppIncluded
