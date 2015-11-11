/*!
	@file utils/PDBReader.hh
	@brief The header file for PDBReader.
	@details This the header file for PDBReader.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#ifndef PDBReaderHhIncluded
#define PDBReaderHhIncluded

#include "pspl/globals/rpoint.hpp"
#include "pspl/globals/protein.hpp"

openPlatypusSpace
class PDBReader
{
  public:
    static void readPDBFile(string filename,block1<RPoint,xmm> & points,int &start,int &end); // read the contents of the file into the kblock
    static void readSequence(string filename,char* protein);
    PDBReader();
    ~PDBReader();

    static int Type(string acdThree);
};


const string threeAA[] ={"CYS","MET","PHE","ILE","LEU","VAL","TRP","TYR","ALA","GLY","THR","SER","GLN","ASN","GLU","ASP","HIS","ARG","LYS","PRO"};
const char chAA[]={'C','M','F','I','L','V','W','Y','A','G','T','S','Q','N','E','D','H','R','K','P'};
closePlatypusSpace
#endif // PDBReaderHhIncluded


