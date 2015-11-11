
#include "pspl/idx.hh"
#include "psps/operator.h"
#include "psps/lnsbio.h"
#include "psps/pclf.h"
#include "psps/motif.h"
#include "psps/sidechain.h"

using namespace platypus;


int main(int argC, char* argv[])
{
    //cout << kangaroo::usqrt<int,float>::iof(2) << endl;
    //exit(0);
    //cout << "Hello World!" << endl;
    //operatorSearch(argC,argv);
    lnsbio(argC,argv);
    //sidechainSearch(argC,argv);
    //lnsmj(argC,argv);
    //motifSearch(argC,argv);
    //memeticSearch(argC,argv);
    return 0;
}




