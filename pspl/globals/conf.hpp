/*!
	@file globals/conf.hpp
	@brief The prototype file for Conf class.
	@details This is the prototype file for Conf class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/




#ifndef ConfHppIncluded
#define ConfHppIncluded


#include "pspl/globals/protein.hpp"
#include "pspl/globals/lattice.hpp"
#include "pspl/globals/data.hpp"
#include "pspl/globals/tuple.hpp"
#include "pspl/globals/change.hpp"



openPlatypusSpace



/*!
	@class Conf
	@brief A class to represent conformations.
	@details This class represents conformations.
*/
class Conf
{
	private:
		Hdl mSysHdl;					//!< The kangaroo system handle.
		Dim mVarCount;					//!< The number of variable.
		hmapd<Point,Dim,nmmh> mPointCounts;	//!< The conformation points and their counts.

		Dim tabuLength ;    //! <tabulaength for variables>

        B sideChain;
        B motif;
        block1<Z,xmm> mMotifs; // 0=coil 1=beta 2=helix
        block1<Pos,xmm> matchPos;

	public:

		Conf(Rnd theRnd);												//!< Constructor.
		~Conf();											//!< Destructor.
		Conf(Conf const & that);							//!< Duplicator.
		Conf const & operator = (Conf const & that);		//!< Assigner.

        B isMotifEnabled() const;
        void setMotif(B tMotif);
        void insertMotifs(char* tMotifs);
        B isPosMotif(Pos tPos) const;
        Z getMotifAtPos(Pos tPos) const;

        Z getMatchPos(Pos tPos) const;

        B isSideChainEnabled() const;
        void setSideChain(B tBool);

        B isPosSideChain(Pos tPos) const;

		Hdl SysHdl() const;									//!< Return sys handle.
		Dim VarCount() const;								//!< Return the number of variables.
		B occupied(Point const & thePoint) const;			//!< Whether a point is occupied.
		Point Coord(Pos const thePos) const;				//!< Coordinate of a position.
		void initialise(Change const & theChange);			//!< Initialise the conformation.

		void execute(Change const & theChange);				//!< Execute with the changes.
		void simulate(Change const & theChange) const;		//!< Simulate with the changes.

		void NonIsoDirs(AbsDirs & absdirs);                 //! Gets NonIso Directions
        void PackIt(AbsDirs const absDirs, Packed & thePck);//! Gets Packed Dirs
		void writeHPToFile(const string filename) const;	//!< Write the conformation to a file.

		void writeMotifToFile(const string filename) const;	//!< Write the conformation to a file.
        void writeHPSCToFile(const string filename) const;

		void writePDB(const string filename,R energy) const;

		void savePoints(Change & theChange)const; //! Return a Change containing all the points
};


/*!
    Return motif boolean

*/
inline B Conf::isMotifEnabled()const
{
    return motif;
}

/*!
    Set/reset boolean motif
*/
inline void Conf::setMotif(B tMotif)
{
    motif=tMotif;
}

/*!
    set all motif positions
*/
inline void Conf::insertMotifs(char *tMotifs)
{
    Idx motifIdx=0;
    cout << "REMARK ";
    for(Idx tIdx=0;tMotifs[tIdx];++tIdx)
    {

        Z code=0;
        if(tMotifs[tIdx]=='E')
        {
            code=1;
        }
        else if(tMotifs[tIdx]=='H')
        {
            code=2;
        }
        if(tMotifs[tIdx]==tMotifs[tIdx-1]&&tMotifs[tIdx-1]!='C'&&tMotifs[tIdx+1]==tMotifs[tIdx])
        {
            motifIdx++;
        }
        else
            motifIdx=0;
        if(code==1 && motifIdx >=2)
        {
            matchPos.insertMem(tIdx-2);
        }
        else if(code==2 && motifIdx >=4 )
        {
            matchPos.insertMem(tIdx-4);
        }
        else
        {
            matchPos.insertMem(tIdx);
        }
        mMotifs.insertMem(code);
        cout << matchPos[tIdx] << " ";
    }
    cout << endl;

}
/*!
    Return whether a position is part of a motif or not
*/
inline B Conf::isPosMotif(Pos tPos) const
{
    WatchError
    if(motif==true&&mMotifs[tPos]>0)
        return true;
    else
        return false;
    CatchError
}



/*!

    Return the motif code for Position
*/
inline Z Conf::getMotifAtPos(Pos tPos)const
{
    return mMotifs[tPos];
}

/*!

    Return the match position for Position
*/
inline Z Conf::getMatchPos(Pos tPos)const
{
    return matchPos[tPos];
}

/*!
    Return sideChain boolean
*/
inline B Conf::isSideChainEnabled() const
{
    return sideChain;
}

/*!
    Return whether a position is part of side chain or not
*/

inline B Conf::isPosSideChain(Pos tPos) const
{
    WatchError
    if(sideChain==true&&tPos>=Protein::p().Length()/2)
        return true;
    else
        return false;
    CatchError
}

/*!
    Set/Reset side chain
*/
inline void Conf::setSideChain(B tBool)
{
    sideChain=tBool;
}

/*!
	Return the number of variable.
*/
inline Dim Conf::VarCount() const
{
	WatchError
	return Protein::p().Length() * Space::Dimen();
	CatchError
}




closePlatypusSpace



#endif //ConfHppIncluded

