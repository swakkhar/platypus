/*!
	@file energys/mjenergy.hpp
	@brief The prototype file for MjEnergy class.
	@details This is the prototype file for MjEnergy class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/


#ifndef MjEnergyHppIncluded
#define MjEnergyHppIncluded


#include "pspl/globals/energy.hpp"


openPlatypusSpace


/*!
	Register HP energy class.
*/
regEnergy(MjEnergy, castN(2));



/*!
	@class MjEnergy
	@brief A class to represent HP energy models.
	@details This class represents HP energy models.
*/
class MjEnergy: public Energy
{
	public:
		virtual Typ Type(Dim const AlphaSize, Acd const theAcd) const;	//!< Return the type of an amino acid.
		virtual Eng Level(Typ const Typ1, Typ const Typ2) const;		//!< Return the energy level between two amino acids.

		MjEnergy();														//!< The constructor.
		MjEnergy(MjEnergy const & that);								//!< The duplicator.
		MjEnergy const & operator = (MjEnergy const & that);			//!< The assigner.
		~MjEnergy();													//!< The destructor.

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
inline MjEnergy::MjEnergy() :
	Energy(EnergyCls<MjEnergy>::ModelId, mTypeSize)
{
	WatchError
	//	nothing to be done.
	CatchError
}



/*!
	The duplicator.
*/
inline MjEnergy::MjEnergy(MjEnergy const & that) : Energy(that)
{
	WatchError
	//	nothing to be done.
	CatchError
}



/*!
	The assigner.
*/
inline MjEnergy const & MjEnergy::operator = (MjEnergy const & that)
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
inline MjEnergy::~MjEnergy()
{
	//	nothing to be done.
}



closePlatypusSpace


#endif //MjEnergyHhIncluded

