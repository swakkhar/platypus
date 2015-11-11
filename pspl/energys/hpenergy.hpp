/*!
	@file energys/hpenergy.hpp
	@brief The prototype file for HpEnergy class.
	@details This is the prototype file for HpEnergy class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/


#ifndef HpEnergyHppIncluded
#define HpEnergyHppIncluded


#include "pspl/globals/energy.hpp"


openPlatypusSpace


/*!
	Register HP energy class.
*/
regEnergy(HpEnergy, castN(1));



/*!
	@class HpEnergy
	@brief A class to represent HP energy models.
	@details This class represents HP energy models.
*/
class HpEnergy: public Energy
{
	public:
		virtual Typ Type(Dim const AlphaSize, Acd const theAcd) const;	//!< Return the type of an amino acid.
		virtual Eng Level(Typ const Typ1, Typ const Typ2) const;		//!< Return the energy level between two amino acids.

		HpEnergy();														//!< The constructor.
		HpEnergy(HpEnergy const & that);								//!< The duplicator.
		HpEnergy const & operator = (HpEnergy const & that);			//!< The assigner.
		~HpEnergy();													//!< The destructor.

	private:
		/*!
			@enum Settings
			@brief Internal settings.
			@details Internal settings.
		*/
		enum Settings
		{
			mTypeSize = 2		//!<	Number of types.
		};
};




/*!
	The constructor.
*/
inline HpEnergy::HpEnergy() :
	Energy(EnergyCls<HpEnergy>::ModelId, mTypeSize)
{
	WatchError
	//	nothing to be done.
	CatchError
}



/*!
	The duplicator.
*/
inline HpEnergy::HpEnergy(HpEnergy const & that) : Energy(that)
{
	WatchError
	//	nothing to be done.
	CatchError
}



/*!
	The assigner.
*/
inline HpEnergy const & HpEnergy::operator = (HpEnergy const & that)
{
	WatchError
	if (this != &that)
		Energy::operator=(that);
	return *this;
	CatchError
}



/*!
	The destructor.
*/
inline HpEnergy::~HpEnergy()
{
	//	nothing to be done.
}



closePlatypusSpace


#endif //HpEnergyHhIncluded
