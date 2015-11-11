/*!
	@file inits/ReadFile.cpp
	@brief The implementation file for ReadFile class.
	@details This is the implementation file for ReadFile class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#include "pspl/utils/readfile.hpp"
#include <fstream>
#include <sstream>

openPlatypusSpace



/*!
	The costructor.
*/
ReadFile::ReadFile()
{
    WatchError
    // do nothing
    CatchError
}


/*!
	The Destructor
*/
ReadFile::~ReadFile()
{
    WatchError
    // do nothing
    CatchError
}

/*!

    This funciton reads from a text file
*/

void ReadFile::readTxtFile(string filename,block1<Point,kmm> &allPoints)
{
    //cout << "In here " << filename<< endl;
    std::ifstream inFile(filename.c_str());
    if(!inFile)
        cout<<"Stream is Null"<<endl;
    Spc a,b,c;
    char line[1000];
    Idx tIdx=0;
    while(inFile.getline(line,1000,'\n'))
    {
        // do nothing
        std::stringstream sin(line);
        sin>>a;
        sin>>b;
        sin>>c;
        //cout<< "table" << a<<" "<<b<<" "<<c<<endl;
        Point p;
        p.Comp(0)=a;
        p.Comp(1)=b;
        p.Comp(2)=c;
        allPoints[tIdx++]=p;
    }
}



closePlatypusSpace
