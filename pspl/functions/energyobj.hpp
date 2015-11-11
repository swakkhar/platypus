/*!
	@file hpmodel/hpengobj.hpp
	@brief The prototype file for EnergyObj class.
	@details This is the prototype file for EnergyObj class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#ifndef EnergyObjHhIncluded
#define EnergyObjHhIncluded



#include "pspl/globals/function.hpp"
#include "pspl/energys/hpenergy.hpp"
#include "pspl/energys/mjenergy.hpp"
#include "pspl/energys/barerraenergy.hpp"


openPlatypusSpace



/*!
	@class EnergyObj
	@brief A class to represent HP energy objective function.
	@details This class represent HP energy objective function.
*/
class EnergyObj : public Function
{
	public:
		EnergyObj(Conf const & theConf,B const NeedHint);							//!< Constructor.
		EnergyObj(EnergyObj const & that);						//!< Duplicator.
		EnergyObj const & operator = (EnergyObj const & that);	//!< Assignment.
		~EnergyObj();											//!< Destructor.
};



closePlatypusSpace



#endif //EnergyObjHhIncluded
