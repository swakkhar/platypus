/*!
	@file globals/energy.hpp
	@brief The prototype file for Energy class.
	@details This is the prototype file for Energy class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#ifndef EnergyHhIncluded
#define EnergyHhIncluded



#include "pspl/globals/space.hpp"
#include "pspl/globals/tuple.hpp"


openPlatypusSpace


/*!
	Forward declarations.
*/
class Conf;



/*!
	@class Energy
	@brief A class to represent energy models.
	@details This class represents energy models.
*/
class Energy
{
	private:
		static Energy const * mEnergy;				//!< Points to the global energy.

	public:
        static Energy const & e(); 					//!< Return the default energy model.
		static void e(Energy & theEnergy); 			//!< Set the default energy model.

	public:
		Mid ModelId() const;						//!< Return model id of the energy.
		Dim TypeCount() const;						//!< Return the number types in the model.

		virtual Typ Type(Dim const AlphaSize, Acd const theAcd) const = 0;	//!< The type value for the acid.
		virtual Eng Level(Typ const Typ1, Typ const Typ2) const = 0;		//!< The energy level for the types.

	protected:

		virtual ~Energy();										//!< Destructor.
		Energy(Energy const & that);							//!< Duplicator.
		Energy const & operator = (Energy const & that);		//!< Assigner.
		Energy(Mid const theModelId, Dim const theTypeCount);	//!< Initialiser.

		Mid mModelId;											//!< Model id.
		Dim mTypeCount;											//!< Number of types.
};



/*!
	@class EnergyId
	@brief A class to represent energy id numbers.
	@details This class represents energy id numbers.
*/
template <Mid ModelId> class EnergyId {};


/*!
	@class EnergyCls
	@brief A class to represent energy class meta data.
	@details This class represents energy class meta data.
*/
template <class ClassName> class EnergyCls{};


/*!
	@def regEnergy(Name, Model)
	@brief Register a energy class.
	@details Register a energy class.
*/
#define regEnergy(Name, Model) class Name; \
	template<> class EnergyId<Model> {public: typedef Name Class;}; \
	template<> class EnergyCls<Name> {public: enum {ModelId = Model};}



/*!
	Return the energy model.
*/
inline Energy const & Energy::e()
{
	WatchError
	Warn(!mEnergy, eEmptyEnergy);
	return *mEnergy;
	CatchError
}



/*!
	User another energy model.
*/
inline void Energy::e(Energy & theEnergy)
{
	WatchError
	Warn(mEnergy, eNonEmptyEnergy);
	mEnergy = &theEnergy;
	CatchError
}



/*!
	Return energy model.
*/
inline Mid Energy::ModelId() const
{
	WatchError
    return mModelId;
	CatchError
}



/*!
	Return the type count of the energy model.
*/
inline Dim Energy::TypeCount() const
{
	WatchError
    return mTypeCount;
	CatchError
}



/*!
	The initialiser.
*/
inline Energy::Energy(Mid const theModelId, Dim const theTypeCount) :
	mModelId(theModelId), mTypeCount(theTypeCount)
{
	WatchError
	Warn(!theTypeCount, eEmptyDimension);
	CatchError
}



/*!
	The duplicator.
*/
inline Energy::Energy(Energy const & that) :
	mModelId(that.mModelId), mTypeCount(that.mTypeCount)
{
	WatchError
	//	nothing to be done.
	CatchError
}



/*!
	The assigner.
*/
inline Energy const & Energy::operator = (Energy const & that)
{
	WatchError
	Warn(this->mModelId != that.mModelId, eEnergyMismatch);
	if (this != &that)
		mTypeCount = that.mTypeCount;
	return *this;
	CatchError
}


closePlatypusSpace


#endif // EnergyHhIncluded
