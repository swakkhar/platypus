/*!
	@file energys/barerraenergy.hpp
	@brief The prototype file for BarerraEnergy class.
	@details This is the prototype file for BarerraEnergy class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/


#ifndef BarerraEnergyHppIncluded
#define BarerraEnergyHppIncluded


#include "pspl/globals/energy.hpp"


openPlatypusSpace


/*!
	Register HP energy class.
*/
regEnergy(BarerraEnergy, castN(3));



/*!
	@class BarerraEnergy
	@brief A class to represent HP energy models.
	@details This class represents HP energy models.
*/
class BarerraEnergy: public Energy
{
	public:
		virtual Typ Type(Dim const AlphaSize, Acd const theAcd) const;	//!< Return the type of an amino acid.
		virtual Eng Level(Typ const Typ1, Typ const Typ2) const;		//!< Return the energy level between two amino acids.

		BarerraEnergy();														//!< The constructor.
		BarerraEnergy(BarerraEnergy const & that);								//!< The duplicator.
		BarerraEnergy const & operator = (BarerraEnergy const & that);			//!< The assigner.
		~BarerraEnergy();													//!< The destructor.

	private:
		/*!
			@enum Settings
			@brief Internal settings.
			@details Internal settings.
		*/
		enum Settings
		{
			mTypeSize = 20		//!<	Number of types.
		};
};




/*!
	The constructor.
*/
inline BarerraEnergy::BarerraEnergy() :
	Energy(EnergyCls<BarerraEnergy>::ModelId, mTypeSize)
{
	WatchError
	//	nothing to be done.
	CatchError
}



/*!
	The duplicator.
*/
inline BarerraEnergy::BarerraEnergy(BarerraEnergy const & that) : Energy(that)
{
	WatchError
	//	nothing to be done.
	CatchError
}



/*!
	The assigner.
*/
inline BarerraEnergy const & BarerraEnergy::operator = (BarerraEnergy const & that)
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
inline BarerraEnergy::~BarerraEnergy()
{
	//	nothing to be done.
}



closePlatypusSpace


#endif //BarerraEnergyHhIncluded


