/*!
	@file globals/data.hpp
	@brief The prototype file for data types.
	@details This is the prototype file for data types.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/


#ifndef DataHppIncluded
#define DataHppIncluded



#include "pspl/prj.hh"
#include "pspl/lib.hh"
#include "pspl/opt.hh"
#include "pspl/msg.hh"


openPlatypusSpace


/*!
	@var typedef Spc
	@brief Space type.
	@details Space type.
*/
typedef Z Spc;



/*!
	@typedef Dst
	@brief Distance type.
	@details Distance type.
*/
typedef Z  Dst;



/*!
	@typedef Eng
	@brief Energy value.
	@details Energy value.
*/
typedef Z  Eng;


/*!
	@typedef Mid
	@brief Model id.
	@details Model id.
*/
typedef S Mid;



/*!
	@typedef Acd
	@brief Amino acid monomer.
	@details Amino acid monomer.
*/
typedef C Acd;



/*!
	@typedef Typ
	@brief Amino acid type.
	@details Amino acid type.
*/
typedef D Typ;



/*!
	@typedef Drn
	@brief Direction name.
	@details Direction name.
*/
typedef C const * Drn;




/*!
	@typedef Dir
	@brief Basis directions.
	@details Basis directions.
*/
typedef S  Dir;



/*!
	@typedef Pos
	@brief Amino acid positions.
	@details Amino acid positions.
*/
typedef D Pos;



/*!
	@typedef Cmp
	@brief Vector or matrix component.
	@details Vector or matrix component.
*/
typedef D Cmp;



/*!
	@var typedef AbsDirs
	@brief Absolute Directions
	@details Array of Dirs.
*/
typedef block1<Dir,xmm> AbsDirs;


/*!
	@var typedef Packed
	@brief Packed Type.
	@details Array of N.
*/
typedef block1<N,xmm> Packed;


closePlatypusSpace


#endif //DataHppIncluded
