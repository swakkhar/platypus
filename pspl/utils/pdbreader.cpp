/*!
	@file inits/PDBReader.cpp
	@brief The implementation file for PDBReader class.
	@details This is the implementation file for PDBReader class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#include "pspl/utils/pdbreader.hpp"
#include <fstream>
#include <sstream>

openPlatypusSpace



/*!
	The costructor.
*/
PDBReader::PDBReader()
{
    WatchError
    // do nothing
    CatchError
}


/*!
	The Destructor
*/
PDBReader::~PDBReader()
{
    WatchError
    // do nothing
    CatchError
}

/*!

    This funciton reads from a text file
*/

void PDBReader::readPDBFile(string filename,block1<RPoint,xmm> &allPoints,int &start,int &end)
{
    //cout << "In here " << filename<< endl;
    std::ifstream file(filename.c_str());
    if(!file)
    {
        cout<<"REMARK Stream is Null! Can't read PDB file! "<<filename<<endl;
        exit(0);
    }


    char buffer[2048];
    Idx i=0;
    bool chainstart=false;
    while(!file.eof())
    {
        file.getline(buffer,2048);
        //cout << buffer << endl;
        string bufString(buffer);

        if(bufString.find("ATOM")==0)
        {
            if(bufString.find("CA")!=string::npos)
            {
                //cout << buffer << endl;
                std::stringstream strIn(std::stringstream::in | std::stringstream::out);
                strIn << buffer;
                string atm;
                string srl;
                string name;
                string typ1;
                string typ2;
                string seq;
                double x,y,z;

                seq=bufString.substr(22,4);

                std::istringstream(bufString.substr(30,8)) >> x;
                std::istringstream(bufString.substr(38,8)) >> y;
                std::istringstream(bufString.substr(46,8)) >> z;
                //strIn >> atm>>srl >>name>>typ1>>typ2>>seq>>x>>y>>z;
                //cout << x<<y <<z<<endl;
                //Point p;
                //p.Comp(0)=(Spc)(x*100);
                //p.Comp(1)=(Spc)(y*100);
                //p.Comp(2)=(Spc)(z*100);
                //cout << i << endl;
                RPoint rp(x,y,z);
                allPoints.insertMem(rp);
                // i is just a counter
                i++;

                if(chainstart==false)
                {
                    chainstart=true;
                    std::istringstream(seq) >> start;
                }
                else
                {
                    std::istringstream(seq)>>end;
                }

            }
        }
        else if(bufString.find("TER")==0)
            break;

        //if(Protein::p().Length()==i) break;


    }
    //cout << filename << " " << i<<" ";
    //if(i!=Protein::p().Length())
    //{
      //  cout << "Error in ATOM co-ordinates reading!"<<endl;
        //exit(0);
    //}
}

int PDBReader::Type(string acdThree)
{
    Idx tIdx =0;
    while (tIdx<20)
    {
        if(acdThree==threeAA[tIdx])
            return tIdx;
        tIdx++;
    }
    return -1;
}

void PDBReader::readSequence(string filename,char* protein)
{
    //cout << "In here " << filename<< endl;
    std::ifstream file(filename.c_str());
    if(!file)
    {
        cout<<"Stream is Null! Can't read PDB file!"<<endl;
        exit(0);
    }


    char buffer[2048];
    Int i=0;
    while(!file.eof())
    {
        file.getline(buffer,2048);
        //cout << buffer << endl;
        string bufString(buffer);
        if(bufString.find("SEQRES")==0)
        {
            std::stringstream strIn(std::stringstream::in | std::stringstream::out);
            strIn << buffer;
            string seqres;
            string srl;
            string chain;
            string noRes;
            string name;

            strIn>> seqres>>srl>>chain>>noRes;
            //cout << noRes<<endl;
            while(strIn>>name)
            {
                if(Type(name)!=-1)
                {
                    //cout <<ch<<endl;
                    protein[i]=chAA[Type(name)];
                }
                else
                {
                    cout << "Error in reading sequence:" << filename << " seq:"<<buffer <<" i:"<<i<< endl;
                    exit(0);
                }
                i++;
            }
            int j=-1;
            std::istringstream(noRes) >> j;
            if(i==j) break;
        }



    }
    protein[i]=0;
    //cout<<i<<endl;
}



closePlatypusSpace
