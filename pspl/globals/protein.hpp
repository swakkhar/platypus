/*!
	@file globals/protein.hpp
	@brief The prototype file for Protein class.
	@details This is the prototype file for Protein class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/


#ifndef ProteinHppIncluded
#define ProteinHppIncluded



#include "pspl/globals/data.hpp"
#include "pspl/globals/energy.hpp"



openPlatypusSpace



/*!
	@class Protien
	@brief A class to represent proteins.
	@details This class represents proteins.
*/
class Protein
{
	private:
		enum { mMaxLength = 1025 };

		static Protein mProtein; 				//!< The global protein.

		Dim mLength;							//!< Length of the protein.
		Dim mSpan;								//!< Span of the protien = Length - 1.
		Acd mAcids[mMaxLength];					//!< Amino acids of the protein.
		Typ mTypes[mMaxLength];					//!< Amino acid types of the protein.

	public:

		Protein();											//!< Constructor.
		Protein(Protein const & that);						//!< Duplicator.
		Protein(Dim const AlphaSize, Acd const * theAcids);	//!< Initialiser.
		Protein const & operator = (Protein const & that); 	//!< Assigner.

		Dim Span() const;							//!< Return protein span.
		Dim Length() const;							//!< Return protein length.

		Acd const * Acids() const;					//!< Return protein acids.
		Typ const * Types() const;					//!< Return acid types.

        Typ Type(Pos const thePos) const;			//!< Return type at a given position.
        Acd Acid(Pos const thePos) const;			//!< Return acid at a given position.

		static void p(Protein const & theProtein);	//!< Set the global protein.
		static Protein const & p();					//!< Get the global protein.

        Dim Mass(Pos const thePos) const;           //! Return the mass of the particular position

        const char * acid3LCode(Pos const thePos) const;


};
const char listAcidTypes[][4]={"CYS","MET","PHE","ILE","LEU","VAL","TRP","TYR","ALA",
"GLY","THR","SER","GLN","ASN","GLU","ASP","HIS","ARG","LYS","PRO"};
const int allMass[]={71,0,103,115,129,147,57,137,113,0,113,113,131,114,0,97,128,156,87,101,0,99,186,0,163,0};



/*!
    return the mass at the position
*/
inline Dim Protein::Mass(Pos const thePos) const
{
    WatchError
    Warn(thePos >= mLength, eInvalidIndex);
    return allMass[mTypes[thePos]];
    CatchError
}

/*!
	Set the global protein.
*/
inline void Protein::p(Protein const & theProtein)
{
	WatchError
	mProtein = theProtein;
	CatchError
}


/*!
	Return the global protien.
*/
inline Protein const & Protein::p()
{
	WatchError
	Warn(!mProtein.mLength, eEmptyProtein);
	return mProtein;
	CatchError
}



/*!
	Return the amino acid type at a given position.
*/
inline Typ Protein::Type(Pos const thePos) const
{
	WatchError
    Warn(thePos >= mLength, eInvalidIndex);
    return mTypes[thePos];
	CatchError
}



/*!
	Return the amino acid at a given position.
*/
inline Acd Protein::Acid(Pos const thePos) const
{
	WatchError
    Warn(thePos >= mLength, eInvalidIndex);
    return mAcids[thePos];
	CatchError
}
/*!
 return a three letter code
*/
inline const char * Protein::acid3LCode(Pos const thePos) const
{
    WatchError
    Warn(thePos >= mLength, eInvalidIndex);
    return listAcidTypes[mTypes[thePos]];
	CatchError
}

/*!
	Return the span of the protein.
*/
inline Dim Protein::Span() const
{
	WatchError
	Warn(!mSpan, eEmptyDimension);
	return mSpan;
	CatchError
}


/*!
	Return the length of the protein.
*/
inline Dim Protein::Length() const
{
	WatchError
	Warn(!mLength, eEmptyDimension);
	return mLength;
	CatchError
}


/*!
	Return the amino acids of the protein.
*/
inline Acd const * Protein::Acids() const
{
	WatchError
	Warn(!mLength, eEmptyDimension);
	return mAcids;
	CatchError
}



/*!
	Return the amino acid types of the protein.
*/
inline Typ const * Protein::Types() const
{
	WatchError
	Warn(!mLength, eEmptyDimension);
	return mTypes;
	CatchError
}



/*!
	The constructor.
*/
inline Protein::Protein() : mLength(0), mSpan(0)
{
	WatchError
	mAcids[0] = 0;
	CatchError
}


closePlatypusSpace



#endif // ProteinHppIncluded
