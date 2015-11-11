#include "pascal.hpp"


Pascal::Pascal(Rnd aRnd, N maxIt, N th, N limit):mRnd(aRnd),mTabuList(Protein::Length(),Lmodel::lm().Count()),mMaxIt(maxIt), mMaxStable(th), mLimit(limit),mRestarts(0)
{
    mTabuTenure=aRnd((N)4, 4 < Protein::Length()/4?Protein::Length()/4:4);


    mDistSum=0; // sum of all distances
    mUnitCost=0; // energy according to the hp Model

    mGlobalBestUnitCost=0;
    mGlobalBestDistCost=numeric_limits<R>::max();

    mIter = 0;
    mMaxContacts = 10000;

}
Pascal::~Pascal()
{
    // do nothing
}
void Pascal::initializeCosts()
{
    mDistSum=mConf.fitness(1);
    mUnitCost=mConf.fitness(0);



    mNeighborH=1;
    mNeighborP=0;
    mLockedNP=0;
    mLockedNH=0;
    mNonImp=0;
    mHasMoved=0; // if there has been a movement we continue the search, otherwise we re-start
    mHasMovedNH=0;
    mHasMovedNP=0;



    mBestUnitCost=mUnitCost;
    mBestDistCost=mDistSum;
    mLastCost=mDistSum;

    for(Cmp tCmp1=0;tCmp1<Protein::Length();++tCmp1)
     for(Cmp tCmp2=0;tCmp2 < Lmodel::lm().Count(); ++tCmp2)
        mTabuList[tCmp1][tCmp2]=0;


}
void Pascal::updateTabuList()
{
    if( Protein::getType(selectedPos)==1)
        for(Cmp tCmpj=0;tCmpj < Lmodel::lm().Count();++tCmpj)
                mTabuList[selectedPos][tCmpj]=mIter+mTabuTenure;
}
void Pascal::updateCosts()
{
    mDistSum=mConf.fitness(1);

    mUnitCost=mConf.fitness(0);

    if(mDistSum<mBestDistCost)
        mBestDistCost=mDistSum;
    if(mDistSum<mGlobalBestDistCost)
            mGlobalBestDistCost=mDistSum;

    if(mUnitCost>mBestUnitCost)
    {
        mBestUnitCost=mUnitCost;
        mNonImp=0;
        bestConf = mConf;
    }
    if(mUnitCost>mGlobalBestUnitCost)
    {
        mGlobalBestUnitCost=mUnitCost;
        globalBestConf = mConf;

    }
}
void Pascal::generatePMoves()
{

    if(!candidatesP.empty())candidatesP.clear();

    for(Cmp mi=0;mi<Protein::Length();++mi)
    {
        if(Protein::getType(mi)==1)continue;
        for(Cmp mv = 0; mv < Lmodel::lm().Count(); ++mv)
        {
            if( Move::mv(0).feasible(mConf,mi,mv)==false)continue;
            if(mTabuList[mi][mv] > mIter) continue;
                tMove tm(mi,mv,0);
                candidatesP.push_back(tm);
        }
    }
}
void Pascal::selectPMove()
{
    while(!candidatesP.empty())
    {
        Idx iTemp=mRnd(0,(int)candidatesP.size()-1);
        tMove selectedMove=(tMove)candidatesP.at(iTemp);
        selectedPos=selectedMove.getPos();
        selectedDir=selectedMove.getDir();
		selectedDist=selectedMove.getDist();
        candidatesP.clear();
        //performTrialPascalMove(selectedPos,selectedDir);
    }

}
void Pascal::generateHMoves()
{
    while(!candidatesH.empty())candidatesH.pop();
    for(N mi=0;mi<Protein::Length();++mi)
    {
		if(Protein::getType(mi)!=1)continue;
		for(Cmp mv=0;mv<Lmodel::lm().Count();++mv)
        {
            if(Move::mv(0).feasible(mConf,mi,mv)==false)continue;
            R dist= Move::mv(0).simulate(mConf,mi,mv,1); // 1 means all pair only //treeDelta(mi);
            if(mTabuList[mi][mv]>mIter&&mLastCost+dist>=mBestDistCost)continue;
            tMove tm(mi,mv,dist);
            candidatesH.push(tm);

        }

    }
}
void Pascal::selectHMove()
{

    if(candidatesH.size()>0)
    {
        vector<tMove> topList;

        tMove tm=candidatesH.top();

        topList.push_back(tm);
        candidatesH.pop();

        while(candidatesH.size()>0 && tm == candidatesH.top())
        {
            tm=candidatesH.top();
            topList.push_back(tm);
            candidatesH.pop();
        }

		tMove selectedMove=topList.at(mRnd(0,(int)topList.size()-1));
        selectedDist=selectedMove.getDist();
        selectedPos=selectedMove.getPos();
        selectedDir=selectedMove.getDir();
        isSelected=true;
        topList.clear();
		while(!candidatesH.empty())
        			candidatesH.pop();
    }
}
void Pascal::searchBasicPascal()
{
    Utility utl;
    Flt timeInit = utl.getTime();

    mConf.randomStructuredFCC(mRnd); // initialize
    mConf.setCache(distanceCache,1);
    mConf.setCache(energyCache,0);
    mConf.initCache();
    initializeCosts();


    while(mIter++<=mMaxIt && mBestUnitCost<mMaxContacts)
	{
		mHasMovedNP = 0;
        myDist = 0;
		selectedPos = -1;
		selectedDir = -1;
		selectedDist = 1;
		isSelected = false;

		//(*mConf).printSoln();

		if(mNeighborH==1)
		{
			mHasMovedNH=0;
			generateHMoves(); // all H moves are generated into HQUEUE
            selectHMove();

			if(isSelected==true)
			{
				if((selectedDist<=0)||(mNeighborP==0))
				{
					//performTrialPascalMove(selectedPos,selectedDir);
					if(mNeighborP==0)
					{
						mLockedNP=mLockedNP+1;
						if(mLockedNP==2*mLimit)
						{
							mLockedNP=0;
							mLockedNH=0;
							mNeighborP=1;
						}
					}
				}
				mHasMovedNH=1;
				myDist=selectedDist;
			}
		}

		if((mHasMovedNH==0)&&(mNeighborH==1))
		{
			myDist=1;
			mNeighborP=1;
			mLockedNP=0;
		}
		else
			mHasMoved=1;

		if((mNeighborH==0)||((myDist>0)&&(mNeighborP==1)))
		{
			mNeighborH=0;
			selectedPos=-1;

			generatePMoves();
            selectPMove();

			mLockedNH=mLockedNH+1;
			if(mLockedNH==mLimit)
			{
				mLockedNH=0;
				mLockedNP=0;
				mNeighborH=1;
				mNeighborP=0;
			}
			mHasMoved=1;
			mHasMovedNP=1;
		}
		if(selectedPos >= 0)
		{
		    Move::mv(0).execute(mConf,selectedPos,selectedDir);
		    //performPascalMove(selectedPos,selectedDir);
			updateTabuList();

            updateCosts();
			mHasMoved=1;
		}
		else
			mHasMoved=0;

		if(mHasMoved==1)
		{
			mHasMoved=0;
		}
		else
		{
            cout << "No moves possible:" <<mIter<< endl;
			mNonImp=mMaxStable+1;
		}
		if(mHasMovedNP==1) ++mNonImp;
		if(mNonImp>mMaxStable)
		{

			// restart
			++mRestarts;
			mConf=globalBestConf;
            initializeCosts();
		}
		mLastCost=mDistSum;
        if(mIter%1000==0)
        {


            cout<<mIter << " " << utl.getTime() - timeInit << " " << mConf.fitness(0)<<" " <<mDistSum<<" "<<mUnitCost<<" "<<mGlobalBestDistCost<<" "<<mGlobalBestUnitCost<< " "<< mRestarts << endl;
        }
	}
//	mConf.writeHPToFile("data.cml");
//    system("java -jar /home/swakkhar/Downloads/jmol-12.2.2/Jmol.jar -s data.cml.script data.cml");


}

