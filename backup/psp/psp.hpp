/*!
	@file psp.hpp
	@brief The header file for psp.
	@details This is the header file for psp.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@date 08.05.2012 QRL NICTA www.nicta.com.au
*/

#include "cbls/idx.hh"

#include "psp/global.hh"

#ifndef PspHhIncluded
#define PspHhIncluded


//using namespace kangaroo;
//using namespace koala;

openKangarooSpace

int connectedPSP();
int proteinSP();
int proteinSPWithAbsDirModel();
void writeHPToFile(const string filename, char* Protein, const RngVi** PointsX, const  RngVi** PointsY,const  RngVi** PointsZ);
void writeHPToFile(const string filename, char* Protein, kblock2<int> & latticePoints);
void inputFromFile(string filename, Itr **latticePoints, int length);
void new_init(char* protein, kblock2<Int> & latticePoints, Rnd & theRnd);

void updateForbiden(int *f, int m);
static int const transformMoves[]={-1,0,9,6,3,2,11,8,5,1,10,7,4};

struct NeighborData {
			const char *name;       //!< name of move
			int vec[3];       //!< absolute move vector
			int mat[3][3];    //!< matrix for applying rel. move
			int invmat[3][3]; //!< inverse of matrix (we are lazy ;-))
		};
static const NeighborData fccNeighborData[] = {
			{"FL",{1,1,0},
			 {{1,0,0},{0,1,0},{0,0,1}},   {{1,0,0},{0,1,0},{0,0,1}}},
			{"LU",{0,1,1},
			 {{0,0,-1},{0,1,0},{1,0,0}},  {{0,0,1},{0,1,0},{-1,0,0}}},
			{"FU",{1,0,1},
			 {{1,0,0},{0,0,-1},{0,1,0}},  {{1,0,0},{0,0,1},{0,-1,0}}},
			{"BL",{-1,1,0},
			 {{0,-1,0},{1,0,0},{0,0,1}},  {{0,1,0},{-1,0,0},{0,0,1}}},
			{"RU",{0,-1,1},
			 {{0,0,-1},{-1,0,0},{0,1,0}}, {{0,-1,0},{0,0,1},{-1,0,0}}},
			{"BU",{-1,0,1},
			 {{0,-1,0},{0,0,-1},{1,0,0}}, {{0,0,1},{-1,0,0},{0,-1,0}}},
			{"FR",{1,-1,0},
			 {{0,1,0},{-1,0,0},{0,0,1}},  {{0,-1,0},{1,0,0},{0,0,1}}},
			{"LD",{0,1,-1},
			 {{0,0,1},{0,1,0},{-1,0,0}},  {{0,0,-1},{0,1,0},{1,0,0}}},
			{"FD",{1,0,-1},
			 {{1,0,0},{0,0,1},{0,-1,0}},  {{1,0,0},{0,0,-1},{0,1,0}}},
			{"BR",{-1,-1,0},
			 {{-1,0,0},{0,-1,0},{0,0,1}}, {{-1,0,0},{0,-1,0},{0,0,1}}},
			{"RD",{0,-1,-1},
			 {{0,0,1},{-1,0,0},{0,-1,0}}, {{0,-1,0},{0,0,-1},{1,0,0}}},
			{"BD",{-1,0,-1},
			 {{0,-1,0},{0,0,1},{-1,0,0}}, {{0,0,-1},{-1,0,0},{0,1,0}}}
			};

closeKangarooSpace

#endif //TestsHhIncluded
