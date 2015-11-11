/*!
	@file inits/chaingrowthinit.cpp
	@brief The implementation file for ChainGrowth class.
	@details This is the implementation file for ChainGrowth class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#include "pspl/inits/chaingrowthinit.hpp"
#include "pspl/globals/conf.hpp"
#include "pspl/globals/change.hpp"
#include "pspl/globals/lattice.hpp"



openPlatypusSpace



/*!
	The costructor.
*/
ChainGrowth::ChainGrowth(Rnd & theRnd) : mRnd(theRnd)
{
    WatchError
    // do nothing
    CatchError
}


/*!
	The Destructor
*/
ChainGrowth::~ChainGrowth()
{
    WatchError
    // do nothing
    CatchError
}


/*!
	The Dupllicator
*/
ChainGrowth::ChainGrowth(ChainGrowth const & that) :
	Init(that), mRnd(that.mRnd)
{
    WatchError
    Alert(&that, eUndefDuplicator);
    CatchError
}


/*!
	The assigner.
*/

ChainGrowth const & ChainGrowth::operator= (ChainGrowth const & that)
{
    WatchError
	Alert(&that, eUndefAssigner);
    return *this;
    CatchError
}


/*!
	Compute and store a possible initialisation change.
*/
void ChainGrowth::compute(Conf & theConf)
{
    WatchError

    block1<R,kmm> avgEnergy(Protein::p().Length());

    for(Idx i=2;i<Protein::p().Length();++i)
    {
        R sum=0;
        Idx bonds=0;
        for(Idx j=i+2;j<Protein::p().Length();++j)
        {
            sum+=Energy::e().Level(Protein::p().Type(i),Protein::p().Type(j));
            ++bonds;
        }
        avgEnergy[i]=sum/bonds;
    }


    rackl<block1<Point,xmm>,kmm> rackNodes(Protein::p().Length());
    Bln isVisited[1000];

    for(Idx tIdx=2;tIdx<Protein::p().Length();++tIdx)
    {
        isVisited[tIdx]=false;
        //selected=0;
    }

    hset<Point,xmmh> tPoints(Protein::p().Length());
    //add first two points
    // - 1
    Point tPoint=0;
    mChange.add(0,tPoint);
    tPoints.insertBll(tPoint);

    // - 2
    tPoint=tPoint+Lattice::l().DirVec(0);
    mChange.add(1,tPoint);
    tPoints.insertBll(tPoint);

    Idx tIdx=2;

    while(tIdx<Protein::p().Length()&&tIdx>=2)
    {
        tPoint=mChange.Destination(tIdx-1);

        if(isVisited[tIdx]==false)
        {
            // rank all possible moves
            block1<double,xmm> rmsdArr;
            for(Idx tIdx2=0;tIdx2<Lattice::l().NeighCount();++tIdx2)
            {
                //cout << tPoint.Comp(0) <<" "<<tPoint.Comp(1) << " "<<tPoint.Comp(2)<<" "<<endl;
                Point nPoint = tPoint + Lattice::l().DirVec(tIdx2);
                if(tPoints.findBll(nPoint))
                    continue;
                else
                {
                    // record the number and add tPoint to the rack
                    double tempVal=0;
                    // in case of pclf init we used to calculate
                    // RMSD distances with already added poitns
                    // here we would calculate energy based on the
                    // already added points
                    Idx freeCount=(tIdx==Protein::p().Span())?11:10;
                    for(Idx ttIdx=0;ttIdx<tIdx-1;++ttIdx)
                    {
                        // check each of the elements in the change array
                        // if they are neighbors with tIdx / nPoint
                        // add energy value
                        Point dPoint = mChange.Destination(ttIdx);
                        Pos dPos =mChange.Position(ttIdx);
                        if(Lattice::l().areNeighbors(dPoint,nPoint))
                        {
                            tempVal+=Energy::e().Level(Protein::p().Type(dPos),Protein::p().Type(tIdx));
                            freeCount--;
                        }
                    }
                    // now calculate free energy

                    tempVal+=(freeCount*avgEnergy[tIdx]);


                    /*for(Idx ttIdx=0;ttIdx<tIdx-1;++ttIdx)
                    {
                        double distM=Tuple::dist(mChange.Destination(ttIdx),nPoint);
                        double distR=sqrt(
                        (mPoints[ttIdx].x-mPoints[tIdx].x)*(mPoints[ttIdx].x-mPoints[tIdx].x)+
                        (mPoints[ttIdx].y-mPoints[tIdx].y)*(mPoints[ttIdx].y-mPoints[tIdx].y)+
                        (mPoints[ttIdx].z-mPoints[tIdx].z)*(mPoints[ttIdx].z-mPoints[tIdx].z));

                        tempVal+=((distR-distM*Lattice::l().latticeConverter())*(distR-distM*Lattice::l().latticeConverter()));
                    }*/

                    // add this value and point to the arrays then sort
                    rmsdArr.insertMem(tempVal);
                    rackNodes[tIdx].insertMem(nPoint);

                }
            }

            // now all the points and values are added
            // if the list is empty then no positions possible; go back to previous position
            if(rmsdArr.itemCount()==0)
            {
                tIdx--;
                continue;
            }
            // sort the arrays in a descending order, so that last item is the preferred
            for(Idx tii=0;tii<rmsdArr.itemCount();++tii)
            {
                for(Idx tjj=tii+1;tjj<rmsdArr.itemCount();++tjj)
                {
                    if(rmsdArr[tii]<rmsdArr[tjj])
                    {
                        // swap
                        double temp=rmsdArr[tii];
                        rmsdArr[tii]=rmsdArr[tjj];
                        rmsdArr[tjj]=temp;

                        Point tempP=rackNodes[tIdx][tii];
                        rackNodes[tIdx][tii]=rackNodes[tIdx][tjj];
                        rackNodes[tIdx][tjj]=tempP;

                    }
                }
            }

            // now sorted, flag as visited and select the top element and remove it
            isVisited[tIdx]=true;
            mChange.add(tIdx,rackNodes[tIdx][rackNodes[tIdx].itemCount()-1]);
            tPoints.insertBll(rackNodes[tIdx][rackNodes[tIdx].itemCount()-1]);
            rackNodes[tIdx].remove();
            tIdx++;
            continue;




        }
        else  // if isVisited == true
        {
            // remove the previous assignment
            Point toRemove=mChange.Destination(mChange.size()-1);
            tPoints.removeItr(toRemove);
            mChange.remove();

            // now select a new point from the racks

            // if rack is empty then fold back
            if(rackNodes[tIdx].itemCount()==0)
            {
                isVisited[tIdx]=false;
                tIdx--;
                continue;
            }
            // else select and remove and continue
            else
            {
                mChange.add(tIdx,rackNodes[tIdx][rackNodes[tIdx].itemCount()-1]);
                tPoints.insertBll(rackNodes[tIdx][rackNodes[tIdx].itemCount()-1]);
                rackNodes[tIdx].remove();
                tIdx++;
                continue;
            }
        }




    }

    CatchError
}



closePlatypusSpace
