#ifndef PascalHppIncluded
#define PascalHppIncluded

#include "global.hh"
#include "conf.hpp"
#include "tmove.hpp"
#include "move.hpp"
#include <queue>
#include <limits>
#include <vector>
#include "utility.h"
#include "fcache.hpp"

class Pascal
{
    private:
    N mTabuTenure;
    R mDistSum;
    N mUnitCost;
    Conf mConf;
    Rnd mRnd;

    kblock2<N> mTabuList;

    N mMaxIt;
    N mMaxStable;
    N mMaxContacts;
    N mLimit;
    N mIter;

    N mRestarts;


    FitnessCache distanceCache;
    FitnessCache energyCache;


    N mNeighborH;
    N mNeighborP;
    N mLockedNP;
    N mLockedNH;
    N mNonImp;


    N mBestUnitCost;
    N mGlobalBestUnitCost;


    Conf bestConf;
    Conf globalBestConf;

    R mBestDistCost;
    R mGlobalBestDistCost;
    R mLastCost;


    R myDist;

    Z selectedPos;
    Dir selectedDir;

	R selectedDist;
	Bln isSelected;

    N mHasMoved; // if there has been a movement we continue the search, otherwise we re-start
    N mHasMovedNH;
    N mHasMovedNP;



    priority_queue<tMove,vector<tMove>,greater_equal<tMove> > candidatesH;
    vector<tMove> candidatesP;

    void initializeCosts();
    void generateHMoves();
    void generatePMoves();
    void selectHMove();
    void selectPMove();
    void updateTabuList();
    void updateCosts();


    public:
    Pascal(Rnd aRnd, N maxIt, N th, N limit);
    void searchBasicPascal();
    ~Pascal();
};
#endif
