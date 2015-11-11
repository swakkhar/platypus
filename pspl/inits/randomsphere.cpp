/*!
	@file inits/randomvalid.cpp
	@brief The implementation file for RandomSphere class.
	@details This is the implementation file for RandomSphere class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#include "pspl/inits/randomsphere.hpp"
#include "pspl/globals/conf.hpp"
#include "pspl/globals/change.hpp"
#include "pspl/globals/lattice.hpp"



openPlatypusSpace



/*!
	The costructor.
*/
RandomSphere::RandomSphere(Rnd & theRnd) : mRnd(theRnd)
{
    WatchError
    // do nothing
    CatchError
}


/*!
	The Destructor
*/
RandomSphere::~RandomSphere()
{
    WatchError
    // do nothing
    CatchError
}


/*!
	The Dupllicator
*/
RandomSphere::RandomSphere(RandomSphere const & that) :
	Init(that), mRnd(that.mRnd)
{
    WatchError
    Alert(&that, eUndefDuplicator);
    CatchError
}


/*!
	The assigner.
*/

RandomSphere const & RandomSphere::operator= (RandomSphere const & that)
{
    WatchError
	Alert(&that, eUndefAssigner);
    return *this;
    CatchError
}


/*!
	Compute and store a possible initialisation change.
*/
void RandomSphere::compute(Conf & theConf)
{
    WatchError
    Dim tNeighCount = Lattice::l().NeighCount();
	Dim tLength = Protein::p().Length();
    hset<Point,xmmh> tPoints(tLength);

	mChange.reset();
	tPoints.clear();

	Point tPoint(0), sPoint;
	mChange.add(0,tPoint);
	tPoints.insertBll(tPoint);

    Z radius=pow(tLength*tNeighCount/4.0,1.0/3.0);
    Idx tIdx=1;
    N retry=0;
    Dir basis=uniform(mRnd,Lattice::l().NeighCount());
	while(tIdx<tLength)
	{
	    if(retry>tNeighCount)
	    {
            // remove points
            tPoints.removeItr(tPoint);
            mChange.remove();
            tPoint = mChange.Destination(mChange.size() - 1);
            --tIdx;
            retry=0;
            basis=uniform(mRnd,tNeighCount);
            continue;
	    }


		sPoint = tPoint + Lattice::l().DirVec(basis);
		if(tPoints.findBll(sPoint))
		{
			basis=uniform(mRnd,tNeighCount);
			retry++;
			continue;
		}

		R dist=Tuple::sqrdist(sPoint,0); //tempAxes[0]*tempAxes[0]+tempAxes[1]*tempAxes[1]+tempAxes[2]*tempAxes[2];
		N effectiveR=uniform(mRnd,radius/2,radius);
		if(dist>=effectiveR*effectiveR)
		{
			basis=uniform(mRnd,tNeighCount);
			retry++;
			continue;
		}
        tPoint = sPoint;
		tPoints.insertBllMem(sPoint);
		mChange.add(tIdx,sPoint);
		tIdx++;
	}
	CatchError
}



closePlatypusSpace
