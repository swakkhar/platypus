/*!
	@file globals/conf.cpp
	@brief The implementation file for Conf class.
	@details This is the implementation file for Conf class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/


#include "pspl/globals/conf.hpp"


openPlatypusSpace


/*!
	The constructor.
*/
Conf::Conf(Rnd theRnd) : mPointCounts(Protein::p().Length() * 2)
{
	WatchError
	mSysHdl = Sys::def();
	//EvalTf::def(mSysHdl);
	//HintTf::def(mSysHdl);

	EvalTi::def(mSysHdl);
	HintTi::def(mSysHdl);

    if(Protein::p().Span()/2 <=12)
        tabuLength=uniform(theRnd,1,4);
    else
    tabuLength = uniform(theRnd,(Dim)12,Protein::p().Span()/2);
	QcSv2Tabu::def(mSysHdl, tabuLength);

	Sys::refm(mSysHdl).setMultiCandExec();
	Sys::refm(mSysHdl).setMultiFixedSimul();
	Sys::refm(mSysHdl).setMultiFlexiSimul();

    sideChain=false;


	Dim const tDimen = Space::Dimen();
	Dim const tLength = Protein::p().Length();
	for(Pos tPos = 0; tPos < tLength; ++tPos)
		for(Cmp tCmp = 0; tCmp < tDimen; ++tCmp)
			StatRangeVarVi::def(mSysHdl, -tLength, tLength);
	CatchError
}



/*!
	The destructor.
*/
Conf::~Conf()
{
	WatchError
	CatchError
}



/*!
	The duplicator.
*/
Conf::Conf(Conf const & that)
{
	WatchError
	Alert(&that, eUndefDuplicator);
	CatchError
}



/*!
	The assigner.
*/
Conf const & Conf::operator = (Conf const & that)
{
	WatchError
	Alert(&that, eUndefAssigner);
	return *this;
	CatchError
}


/*!
	Return a handle to the kangaroo system.
*/
Hdl Conf::SysHdl() const
{
	WatchError
	return mSysHdl;
	CatchError
}



/*!
	Return whether a point is occupied.
*/
B Conf::occupied(Point const & thePoint) const
{
	WatchError
	Idx tIdx = mPointCounts.findItr(thePoint);
	return tIdx != InvIdx && mPointCounts.item(tIdx);
	CatchError
}



/*!
	Return the point at the position.
*/
Point Conf::Coord(Pos const thePos) const
{
	WatchError
	Point tPoint;
	for (Cmp tCmp = 0; tCmp < Space::Dimen(); ++tCmp)
		tPoint.Comp(tCmp) = StatRangeVarVi::refc(mSysHdl,thePos * Space::Dimen() + tCmp).ValueRec().CurrData();
	return tPoint;
	CatchError
}


/*!
	Initialise the conformation.
*/
void Conf::initialise(Change const & theChange)
{
	WatchError
	mPointCounts.clear();
	Warn(theChange.size() != Protein::p().Length(), eDimensionMismatch);
	Hdl tVarHdls[theChange.size() * Space::Dimen()]; // array of Hdls
    Wrp tVarWrps[theChange.size() * Space::Dimen()]; // array of wraps

//	Hdl tValItrs[theChange.size() * Space::Dimen()];
	for(Idx tIdx = 0; tIdx < theChange.size(); ++tIdx)
	{
		Warn(theChange.Position(tIdx) != tIdx, eValNotAssignable);
		for(Cmp tCmp = 0; tCmp < Space::Dimen(); ++tCmp)
        {
            Idx ttIdx = tIdx * Space::Dimen() + tCmp;
            tVarHdls[ttIdx] = theChange.Position(tIdx) * Space::Dimen() + tCmp;
            tVarWrps[ttIdx] = theChange.Destination(tIdx).Comp(tCmp);
			//cerr << theChange.Destination(tIdx).Comp(tCmp) << " ";
        }
		//cerr << "| ";
//			tValItrs[tIdx * Space::Dimen() + tCmp] =
	//				theChange.Destination(tIdx).Comp(tCmp) + Protein::p().Length();
		Idx sIdx = mPointCounts.insertItrIfOld(theChange.Destination(tIdx), 1);
		if (sIdx != InvIdx) ++mPointCounts.item(sIdx);
	}
	//cerr << endl;
	Sys::refm(mSysHdl).initialiseVarsWrap(tVarHdls,tVarWrps);
	CatchError
}



/*!
	Execute with the given changes.
*/
void Conf::execute(Change const & theChange)
{
	WatchError
	Hdl tVarHdls[theChange.size() * Space::Dimen()];
    Wrp tVarWrps[theChange.size() * Space::Dimen()];


	for(Idx tIdx = 0; tIdx < theChange.size(); ++tIdx)
	{
	    for(Cmp tCmp = 0; tCmp < Space::Dimen(); ++tCmp)
		{
		    Idx ttIdx = tIdx * Space::Dimen() + tCmp;
		    tVarHdls[ttIdx] = theChange.Position(tIdx) * Space::Dimen() + tCmp;
            tVarWrps[ttIdx] = theChange.Destination(tIdx).Comp(tCmp);
        }
		Idx sIdx = mPointCounts.insertItrIfOld(theChange.Destination(tIdx), 1);
		if (sIdx != InvIdx) { ++mPointCounts.item(sIdx); }
	}
	for(Idx tIdx = 0; tIdx < theChange.size(); ++tIdx)
	{
		Idx sIdx = mPointCounts.findItr(Coord(theChange.Position(tIdx)));

		Warn(sIdx == InvIdx, eKeyNotFound);
		if (--mPointCounts.item(sIdx) == 0) mPointCounts.removeWithItr(sIdx);
	}
	Sys& tSys=Sys::refm(mSysHdl);
    Dim tVarCount = theChange.size() * Space::Dimen();
	if(tVarCount < tSys.varCount()/4)
        tSys.execIncrDiffVarsWrap(tVarHdls, tVarWrps, theChange.size() * Space::Dimen());
	else
        tSys.execAnewVarsWrap(tVarHdls, tVarWrps, theChange.size() * Space::Dimen());

	CatchError
}


/*!
	Simulate with the given changes.
*/
void Conf::simulate(Change const & theChange) const
{
	WatchError
	Hdl tVarHdls[theChange.size() * Space::Dimen()];
    Wrp tVarWrps[theChange.size() * Space::Dimen()];

	//Hdl tValItrs[theChange.size() * Space::Dimen()];
	for(Idx tIdx = 0; tIdx < theChange.size(); ++tIdx)
		for(Cmp tCmp = 0; tCmp < Space::Dimen(); ++tCmp)
		{
		    Idx ttIdx = tIdx * Space::Dimen() + tCmp;
		    tVarHdls[ttIdx] = theChange.Position(tIdx) * Space::Dimen() + tCmp;

			//tVarHdls[tIdx * Space::Dimen() + tCmp] = theChange.Position(tIdx) * Space::Dimen() + tCmp;
			//tValItrs[tIdx * Space::Dimen() + tCmp] =
				//theChange.Destination(tIdx).Comp(tCmp) + Protein::p().Length();
            tVarWrps[ttIdx] = theChange.Destination(tIdx).Comp(tCmp);
        }

    Sys& tSys=Sys::refm(mSysHdl);
    //Dim tVarCount = theChange.size() * Space::Dimen();
//	if(tVarCount < tSys.varCount()/4)
        tSys.simulIncrDiffVarsWrap(tVarHdls, tVarWrps, theChange.size() * Space::Dimen());
//    else
//       tSys.simulAnewVarsWrap(tVarHdls, tVarWrps, theChange.size() * Space::Dimen());
	CatchError
}

void Conf::savePoints(Change & theChange)const
{
    WatchError
    theChange.reset();
    for(Pos tPos=0;tPos<Protein::p().Length();++tPos)
    {
        theChange.add(tPos,Coord(tPos));
    }
    CatchError
}


/*!
    this function returns isometric directions
    from the co-ordinates
*/
void Conf::NonIsoDirs(AbsDirs & absdirs)
{
        N contactNum = Lattice::l().NeighCount();
        block1<Int,xmm> tMap(contactNum);

        for(Idx tIdx = 0; tIdx < contactNum ; ++tIdx)
            tMap[tIdx]=-1;

        N countDir=0;
        Idx tIdx = 0;
        for(  ; tIdx < Protein::p().Span(); ++tIdx)
        {
            Dir dirVal = Lattice::l().Direction(Coord(tIdx+1)-Coord(tIdx)); //conf.absdir(tCmp);
            if(countDir < contactNum && tMap[dirVal] == -1)
            {
                tMap[dirVal]=countDir;
                ++countDir;
            }
            absdirs[tIdx]=tMap[dirVal];
        }
}
/*!
    This funciton packs a conformation in to Packed var
*/
void Conf::PackIt(AbsDirs const absDirs, Packed & thePck)
{
    thePck.clear();
    N t=0;
    Idx tIdx = 0;
    for(  ; tIdx < Protein::p().Span(); ++tIdx)
    {

        N mapVal= absDirs[tIdx];
        t = t | mapVal;

        if ( tIdx % 4 == 3)
        {
                thePck[tIdx/4] = t;
                t = 0;
        }
        else
            t = t << 4;
    }
    if( tIdx %4 != 0) thePck[tIdx/4] = t;

}


void Conf::writePDB(const string filename, R energy) const
{
    //std::ofstream out(/*filename.c_str()*/stdout);
    cout << "REMARK " << energy<<endl;
    /*!
    secondary structure information
    */
    if(isMotifEnabled())
    {
            int motifIDx=0;
            int helixId=0;
            int strandId=0;
            for(Idx tIdx=0; tIdx<Protein::p().Length(); ++ tIdx)
            {
                if(isPosMotif(tIdx))
                {
                    if(motifIDx==0)
                    {
                        /// start of sheet or helix
                        if(getMotifAtPos(tIdx)==1) /// sheet
                        {
                            strandId++;
                            /// 1-5 SHEET
                            cout << "SHEET";
                            /// 7-10 strand nbumber
                            cout.width(5); cout << std::right << strandId;
                            /// 12-14 sheet ID - treat all as different for the time being
                            char hId[]=" HA#";
                            hId[3]='A'+strandId;
                            cout.width(4); cout << std::left <<hId;
                            /// 15-16 number of strands - should be always one since all are diff sheets
                            cout.width(2); cout << std::right <<"1";
                            /// 18-20 init residue name 17 18 19 20
                            cout.width(4); cout << std::right << Protein::p().acid3LCode(tIdx);
                            /// 22 init chain id always A
                            cout.width(2); cout << std::right << "A";
                            /// 23-26 initi SeqNumber 23 24 25 26
                            cout.width(4); cout << std::right << tIdx+1;
                        }
                        else /// helix
                        {
                            helixId++;
                            /// 1-5 HELIX
                            cout << "HELIX";
                            ///7-10 serial number of the helix
                            cout.width(5); cout << std::right << helixId;
                            /// 12-14	helix id
                            char hId[]=" HA#";
                            hId[3]='A'+helixId;
                            cout.width(4); cout << std::left <<hId;
                            /// 16-18 init residue name
                            cout.width(4); cout << std::right << Protein::p().acid3LCode(tIdx);
                            /// 20 chain ID always set to A
                            cout.width(2); cout << std::right << "A";
                            /// 22-25 sequence number 21 22 23 24 25
                            cout.width(5); cout << std::right << tIdx+1;

                        }
                    }
                    motifIDx++;
                }
                else
                {
                    if(motifIDx>0) /// end of sheet or helix
                    {
                        if(getMotifAtPos(tIdx)==1) /// sheet
                        {
                            /// 27 init I code - leave optional
                            /// 28-31 end residue name 27 28 29 30 31
                            cout.width(5); cout << std::right << Protein::p().acid3LCode(tIdx);
                            /// 33 end chain ID always A
                            cout.width(2); cout << std::right << "A";
                            /// 34-37 end eq number 34 35 36 37
                            cout.width(4); cout << std::right << tIdx+1;
                            /// 38 end I Code - leave optional
                            /// 39-40 sense always 0
                            cout.width(3); cout << std::right << "0";
                            cout << endl;
                        }
                        else /// helix
                        {
                            /// 26 initICode - leave optional
                            /// 28-30 end residure name
                            /// 26 27 28 29 30
                            cout.width(5); cout << std::right << Protein::p().acid3LCode(tIdx);
                            /// 32 chain ID always set to A
                            cout.width(2); cout << std::right << "A";
                            /// 34-37 sequence number 33 34 35 36 37
                            cout.width(5); cout << std::right << tIdx+1;
                            /// 38 enICode - leave optional
                            /// 39-40 helix class always 1
                            cout.width(3); cout << std::right << "1";
                            cout << endl;
                        }
                    }
                    motifIDx=0;
                }

            }

    }
    for(Idx tIdx=0; tIdx<Protein::p().Length(); ++ tIdx)
    {
        // 1-5	"HELIX"
        cout << "ATOM";
        // 7-11	Atom serial number	right	integer
        cout.width(7); cout << std::right << tIdx+1;
        // 13-16	Atom name	left*	character
        // 17	Alternate location indicator		character
        cout.width(6); cout << std::left << "  CA";
        // 18-20§	Residue name	right	character
        cout.width(3); cout << std::right << Protein::p().acid3LCode(tIdx);
	// 22	Chain identifier		character
        // We give the name A to the chain identifier, but this may be wrong
	cout.width(2); cout << std::right << "A";
	// 23-26	Residue sequence number	right	integer
	cout.width(4); cout << std::right << tIdx+1;
	// 27	Code for insertions of residues		character
	cout.width(1);
	// 31-38 	X orthogonal  coordinate	right	floating
	cout.width(9); cout << std::right << (R)Coord(tIdx).Comp(0)*2.69;
	// 39-46	Y orthogonal Å coordinate	right	floating
	cout.width(8); cout << std::right << (R)Coord(tIdx).Comp(1)*2.69;
	// 47-54	Z orthogonal Å coordinate	right	floating
	cout.width(8); cout << std::right << (R)Coord(tIdx).Comp(2)*2.69;
	// 55-60	Occupancy	right	floating
	cout.width(6); cout << std::right;
	// 61-66	Temperature factor	right	floating
	cout.width(6); cout << std::right;
	// 77-78	Element symbol	right	character
	// In the RNA we specify the C1' carbon
	cout.width(12); cout << std::right << "C";
	// 79-80	Charge		character
	cout.width(2); cout << endl;
 }
 cout << "TER" << endl;
/*!
for(Idx tIdx=1;tIdx<Protein::p().Length();++tIdx)
	{
	    cout << "CONNECT ";
	    cout.width(4);
	    cout << std::right << tIdx;
	    cout.width(5);
	    cout << std::right << tIdx+1 << endl;

	}*/
}

/*!

Write the conformation into a file with motif colors
*/

void Conf::writeMotifToFile(const string filename) const
{
    	std::ofstream out(filename.c_str());
	 /*
	for(int i=0;i<protein->sequenceLength;i++)
	{
		out<<latticePoints[i][0]<<" "<<latticePoints[i][1]<<" "<<latticePoints[i][2]<<" "<<(char)('A'+protein->getProteinType(i))<<endl;;
	}
	out.close();*/
	Dim const tLength = Protein::p().Length();
	int distance=3;
		int shift = tLength*distance;	// shift output to display with positive dimensions only

				out	<<"<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>"
					<<"\n\n<!-- FOR BEST VIEWING e.g. IN JMOL USE AFTER LOADING THE SCRIPT -->"
					<<"\n<!-- select all; color bonds grey; -->"
					<<"\n\n<!-- created with CPSP-package : http://www.bioinf.uni-freiburg.de/sw/cpsp/ -->"
					<<"\n"
					<<"\n<list dictRef=\"cdk:model\" xmlns=\"http://www.xml-cml.org/schema\">"
					<<"\n  <list>"
					<<"\n    <molecule id=\"lattice-protein\">"
					<<"\n    <!-- model    : HP-model -->"
					<<"\n    <!-- sequence : " <<Protein::p().Acids()<<" -->"
					<<"\n    <!-- encoding : H = Cl-atom, P = C-atom  -->"
					<<"\n    <!-- lattice  : " <<"fcc"<<" -->"
					<<"\n      <atomArray>";
				//cout << "---------------"<<endl;
				for (Idx i = 0; i < tLength; ++i) {
				    //cout << i<<": "<<Coord(i).Comp(0) << " "<<Coord(i).Comp(1) << " "<<Coord(i).Comp(2) <<endl;
				    out <<"\n        <atom id=\""<<i+1<<"_"<<Protein::p().Acid(i)<<"\""
						<<" elementType=\"C"<<(getMotifAtPos(i)==1?"o":getMotifAtPos(i)==2?"d":"")<<"\""
						<<" x3=\""<<(shift + (distance * Coord(i).Comp(0)))<<"\""
						<<" y3=\""<<(shift + (distance * Coord(i).Comp(1)))<<"\""
						<<" z3=\""<<(shift + (distance * Coord(i).Comp(2)))<<"\""
						<<" />";
				}
				out	<<"\n      </atomArray>"
					<<"\n      <bondArray>";
				for (Idx i = 1; i < tLength; i++) {
					//if(!isPosSideChain(i))
					out	<<"\n        <bond id=\"b"<<i<<"\" atomRefs2=\""<<i<<"_"<<Protein::p().Acid(i-1)<<" "<<i+1<<"_"<<Protein::p().Acid(i)<<"\" order=\"S\"/>";
                    //else
                    //out	<<"\n        <bond id=\"b"<<i<<"\" atomRefs2=\""<<i-Protein::p().Length()/2+1<<"_"<<Protein::p().Acid(i-Protein::p().Length()/2)<<" "<<i+1<<"_"<<Protein::p().Acid(i)<<"\" order=\"S\"/>";
				}

				out	<<"\n      </bondArray>"
					<<"\n    </molecule>"
					<<"\n  </list>"
					<<"\n</list>"
					<<"\n\n<!-- FOR BEST VIEWING e.g. IN JMOL USE THE SCRIPT -->"
					<<"\n<!-- select all; color bonds grey; -->"
					<<std::endl;
				out.flush();
				out.close();
}

/*!
	Write the conformation into a file.
*/
void Conf::writeHPToFile(const string filename) const
{

	std::ofstream out(filename.c_str());
	 /*
	for(int i=0;i<protein->sequenceLength;i++)
	{
		out<<latticePoints[i][0]<<" "<<latticePoints[i][1]<<" "<<latticePoints[i][2]<<" "<<(char)('A'+protein->getProteinType(i))<<endl;;
	}
	out.close();*/
	Dim const tLength = Protein::p().Length();
	int distance=3;
		int shift = tLength*distance;	// shift output to display with positive dimensions only

				out	<<"<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>"
					<<"\n\n<!-- FOR BEST VIEWING e.g. IN JMOL USE AFTER LOADING THE SCRIPT -->"
					<<"\n<!-- select all; color bonds grey; -->"
					<<"\n\n<!-- created with CPSP-package : http://www.bioinf.uni-freiburg.de/sw/cpsp/ -->"
					<<"\n"
					<<"\n<list dictRef=\"cdk:model\" xmlns=\"http://www.xml-cml.org/schema\">"
					<<"\n  <list>"
					<<"\n    <molecule id=\"lattice-protein\">"
					<<"\n    <!-- model    : HP-model -->"
					<<"\n    <!-- sequence : " <<Protein::p().Acids()<<" -->"
					<<"\n    <!-- encoding : H = Cl-atom, P = C-atom  -->"
					<<"\n    <!-- lattice  : " <<"fcc"<<" -->"
					<<"\n      <atomArray>";
				//cout << "---------------"<<endl;
				for (Idx i = 0; i < tLength; ++i) {
				    //cout << i<<": "<<Coord(i).Comp(0) << " "<<Coord(i).Comp(1) << " "<<Coord(i).Comp(2) <<endl;
				    out <<"\n        <atom id=\""<<i+1<<"_"<<Protein::p().Acid(i)<<"\""
						<<" elementType=\"C"<<(Protein::p().Type(i)==1?"l":"")<<"\""
						<<" x3=\""<<(shift + (distance * Coord(i).Comp(0)))<<"\""
						<<" y3=\""<<(shift + (distance * Coord(i).Comp(1)))<<"\""
						<<" z3=\""<<(shift + (distance * Coord(i).Comp(2)))<<"\""
						<<" />";
				}
				out	<<"\n      </atomArray>"
					<<"\n      <bondArray>";
				for (Idx i = 1; i < tLength; i++) {
					if(!isPosSideChain(i))
					out	<<"\n        <bond id=\"b"<<i<<"\" atomRefs2=\""<<i<<"_"<<Protein::p().Acid(i-1)<<" "<<i+1<<"_"<<Protein::p().Acid(i)<<"\" order=\"S\"/>";
                    else
                    out	<<"\n        <bond id=\"b"<<i<<"\" atomRefs2=\""<<i-Protein::p().Length()/2+1<<"_"<<Protein::p().Acid(i-Protein::p().Length()/2)<<" "<<i+1<<"_"<<Protein::p().Acid(i)<<"\" order=\"S\"/>";
				}

				out	<<"\n      </bondArray>"
					<<"\n    </molecule>"
					<<"\n  </list>"
					<<"\n</list>"
					<<"\n\n<!-- FOR BEST VIEWING e.g. IN JMOL USE THE SCRIPT -->"
					<<"\n<!-- select all; color bonds grey; -->"
					<<std::endl;
				out.flush();
				out.close();
}


/*!
	Write the conformation into a file.
*/
void Conf::writeHPSCToFile(const string filename) const
{

	std::ofstream out(filename.c_str());
	 /*
	for(int i=0;i<protein->sequenceLength;i++)
	{
		out<<latticePoints[i][0]<<" "<<latticePoints[i][1]<<" "<<latticePoints[i][2]<<" "<<(char)('A'+protein->getProteinType(i))<<endl;;
	}
	out.close();*/
	Dim const tLength = Protein::p().Length();
	int distance=3;
		int shift = tLength*distance;	// shift output to display with positive dimensions only

				out	<<"<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>"
					<<"\n\n<!-- FOR BEST VIEWING e.g. IN JMOL USE AFTER LOADING THE SCRIPT -->"
					<<"\n<!-- select all; color bonds grey; -->"
					<<"\n\n<!-- created with CPSP-package : http://www.bioinf.uni-freiburg.de/sw/cpsp/ -->"
					<<"\n"
					<<"\n<list dictRef=\"cdk:model\" xmlns=\"http://www.xml-cml.org/schema\">"
					<<"\n  <list>"
					<<"\n    <molecule id=\"lattice-protein\">"
					<<"\n    <!-- model    : HP-model -->"
					<<"\n    <!-- sequence : " <<Protein::p().Acids()<<" -->"
					<<"\n    <!-- encoding : H = Cl-atom, P = C-atom  -->"
					<<"\n    <!-- lattice  : " <<"fcc"<<" -->"
					<<"\n      <atomArray>";
				//cout << "---------------"<<endl;
				for (Idx i = 0; i < tLength; ++i) {
				    //cout << i<<": "<<Coord(i).Comp(0) << " "<<Coord(i).Comp(1) << " "<<Coord(i).Comp(2) <<endl;
				    if(!isPosSideChain(i))
				    out <<"\n        <atom id=\""<<i+1<<"_"<<Protein::p().Acid(i)<<"\""
						<<" elementType=\"C"<<(Protein::p().Type(i)==1?"l":"")<<"\""
						<<" x3=\""<<(shift + (distance * Coord(i).Comp(0)))<<"\""
						<<" y3=\""<<(shift + (distance * Coord(i).Comp(1)))<<"\""
						<<" z3=\""<<(shift + (distance * Coord(i).Comp(2)))<<"\""
						<<" />";
				    else

				    out <<"\n        <atom id=\""<<i+1<<"_"<<Protein::p().Acid(i)<<"\""
						<<" elementType=\"C"<<(Protein::p().Type(i)==1?"l":"o")<<"\""
						<<" x3=\""<<(shift + (distance * Coord(i).Comp(0)))<<"\""
						<<" y3=\""<<(shift + (distance * Coord(i).Comp(1)))<<"\""
						<<" z3=\""<<(shift + (distance * Coord(i).Comp(2)))<<"\""
						<<" />";
				}
				out	<<"\n      </atomArray>"
					<<"\n      <bondArray>";
				for (Idx i = 1; i < tLength; i++) {
					if(!isPosSideChain(i))
					out	<<"\n        <bond id=\"b"<<i<<"\" atomRefs2=\""<<i<<"_"<<Protein::p().Acid(i-1)<<" "<<i+1<<"_"<<Protein::p().Acid(i)<<"\" order=\"S\"/>";
                    else
                    out	<<"\n        <bond id=\"b"<<i<<"\" atomRefs2=\""<<i-Protein::p().Length()/2+1<<"_"<<Protein::p().Acid(i-Protein::p().Length()/2)<<" "<<i+1<<"_"<<Protein::p().Acid(i)<<"\" order=\"S\"/>";
				}

				out	<<"\n      </bondArray>"
					<<"\n    </molecule>"
					<<"\n  </list>"
					<<"\n</list>"
					<<"\n\n<!-- FOR BEST VIEWING e.g. IN JMOL USE THE SCRIPT -->"
					<<"\n<!-- select all; color bonds grey; -->"
					<<std::endl;
				out.flush();
				out.close();
}


closePlatypusSpace
