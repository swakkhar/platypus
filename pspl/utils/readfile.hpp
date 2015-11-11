/*!
	@file utils/readfile.hh
	@brief The header file for ReadFile.
	@details This the header file for ReadFile.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#ifndef ReadFileHhIncluded
#define ReadFileHhIncluded

#include "pspl/globals/tuple.hpp"

openPlatypusSpace
class ReadFile
{
  public:
    static void readTxtFile(string filename,block1<Point,kmm> & points); // read the contents of the file into the kblock
    ReadFile();
    ~ReadFile();
};

closePlatypusSpace
#endif // ReadFileHhIncluded

