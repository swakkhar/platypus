/*!
	@file moves/LNSMove.cpp
	@brief The implementation file for LNSMove class.
	@details This is the implementation file for LNSMove class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 05.08.2012 QRL NICTA www.nicta.com.au
*/

#include "pspl/moves/lnsmove.hpp"
#include "pspl/globals/conf.hpp"
#include "pspl/globals/change.hpp"
#include "pspl/globals/lattice.hpp"


openPlatypusSpace


/*!
The costructor
*/
LNSMove::LNSMove()
{
    WatchError
    // do nothing
    CatchError
}


/*!

The Destructor

*/
LNSMove::~LNSMove()
{
    WatchError

    // do nothing

    CatchError
}


/*!

The Dupllicator
*/

LNSMove::LNSMove(LNSMove const & that):Move(that)
{
    WatchError

    Alert(&that, eUndefDuplicator);

    CatchError
}


/*!

Assignment Operator
*/

LNSMove const & LNSMove::operator= (LNSMove const & that)
{
    WatchError

    Move::operator=(that);
    return *this;


    CatchError
}


/*!

The compute method
save all possible move results into mChanges

In the case of the LNS move thePos is a dummy
variable, it contains another internal array
aPoints that contains all the variables; it is
assumed that the array is sorted if the
array contains any element that is not valid or
any element twice it will crash, for the time
being we can trust the caller


dont allow the call with initial position


need to add shortest path and self check

*/
void LNSMove::compute(Conf & theConf, Pos const thePos)
{
    //WatchError
    mChanges.clear();

    //cout << "In Compute:" << endl;


    // create a dictionary using all the positions
    hset<Pos,xmmh> posDict;
    for(Idx tIdx=0;tIdx<aPositions.itemCount();++tIdx)
    {
        posDict.insertBllMem(aPositions[tIdx]);
    }
    // now create another dictionary by taking only the co-ordinates
    // leaving out the positions only
    hset<Point,xmmh> pointsDict;
    for(Idx tIdx = 0; tIdx < Protein::p().Length(); ++tIdx)
    {
        if(!posDict.findBll(tIdx))
        {
            pointsDict.insertBllMem(theConf.Coord(tIdx));
        }
    }
    block1<Dir,kmm> prior(aPositions.itemCount());

    // generator
    block1<Dir,kmm> generator(aPositions.itemCount());
    for(Idx tIdx=0;tIdx<aPositions.itemCount();++tIdx)
        generator[tIdx]=0; // initialize by zero

    // it is also assumed that there are no duplicates

    // sort the positions array please
    // for the time being its just and insertion sort
    for(Idx tIdx1=0;tIdx1<aPositions.itemCount();++tIdx1)
    {
        for(Idx tIdx2=tIdx1+1;tIdx2<aPositions.itemCount();++tIdx2)
        {
            if(aPositions[tIdx1]>=aPositions[tIdx2])
            {
                // swap
                Pos temp = aPositions[tIdx1];
                aPositions[tIdx1]=aPositions[tIdx2];
                aPositions[tIdx2]=temp;
            }
        }
    }
    for(Idx tIdx1=0;tIdx1<aPositions.itemCount();++tIdx1)
    {
        prior[tIdx1]=Lattice::l().Direction(theConf.Coord(aPositions[tIdx1])-theConf.Coord(aPositions[tIdx1]-1));
    }
    // if sorted then we can start doing it
    bool tSkip=true;
    while(1)
    {
        //for each generated string try to build the chain please
        Change tChange;
        Point lastPoint;
        hset<Point,xmmh> tempDict;
        tempDict.clear();
        /*cout<<"genrator: ";
        for(Idx ttIdx=0;ttIdx<generator.size();++ttIdx)
        {
          cout<<(int)generator[ttIdx];
        }
        cout<<endl;*/
        Idx tIdx =0;
        tSkip=true;

        //if(equal(generator,prior))
        {
            //if(!next(generator)) break;
        }


        for(;tIdx < aPositions.itemCount(); ++tIdx)
        {
            Pos thisPos = aPositions[tIdx];
            Pos lastPos = thisPos - 1;
            Point newPoint;
            //cout<<"lastPos:" << lastPos<< endl;
            if(posDict.findBll(lastPos))
            {
                    newPoint = lastPoint + Lattice::l().DirVec(generator[tIdx]);
              //      cout << "last " <<lastPoint.Comp(0) << " " << lastPoint.Comp(1) << " " << lastPoint.Comp(2) << " " << endl;
            }
            else
            {
                newPoint = theConf.Coord(lastPos) + Lattice::l().DirVec(generator[tIdx]);
                //cout << "conf "<<theConf.Coord(lastPos).Comp(0) << " " << theConf.Coord(lastPos).Comp(1) << " " << theConf.Coord(lastPos).Comp(2) << " " << endl;
            }


            //cout << "vector "<<Lattice::l().DirVec(generator[tIdx]).Comp(0) << " " << Lattice::l().DirVec(generator[tIdx]).Comp(1) << " " << Lattice::l().DirVec(generator[tIdx]).Comp(2) << " " << endl;
            //cout << newPoint.Comp(0) << " " << newPoint.Comp(1) << " " << newPoint.Comp(2) << " " << endl;
            // now check the new point, if it is already occupied then end

            if(pointsDict.findBll(newPoint)||tempDict.findBll(newPoint))
            {
                // skip the rest of the upddate generator string and
                // break
                tSkip=skip(generator,tIdx);
                break;
            }
            else
            {
                // need some other checks
                // now if this point is the end of some subchain and not the end of the
                // sequence, we must check whether it is in contact with the
                // next point
                if(thisPos!=Protein::p().Span()&&
                   (tIdx==generator.itemCount()-1 /* last item in the list*/
                    || aPositions[tIdx+1]!=thisPos+1 /*last item in the subsequence*/))
                {
                    if(!Lattice::l().areNeighbors(newPoint,theConf.Coord(thisPos+1)))
                    {
                        tSkip=skip(generator,tIdx);
                        break;
                    }
                }
                // add the new point
                // update the generator string
                // add the newPoint into the tempDict
                // update the lastPoint
                tChange.add(thisPos,newPoint);
                tempDict.insertBllMem(newPoint);
                lastPoint=newPoint;


            }

        }
        //exit(0);

        if(!tSkip)
        {
            break;
        }
        if(tIdx==aPositions.itemCount())  // means the chain ended successfully
        {
            /*for(int i=0;i<tChange.size();++i)
			{
                cout << "("<<tChange.Destination(i).Comp(0) << " "<< tChange.Destination(i).Comp(1) << " "
                << tChange.Destination(i).Comp(2) << " -" <<tChange.Position(i)<<") ";
			}
			cout << endl;*/
            //cout<<tChange.Destination(0).Comp(0);
            mChanges.insertMem(tChange);
            if(!next(generator))
                break; // all the string expired
        }

    }


	//CatchError
}

bool LNSMove::next(block1<Dir,kmm> &generator)
{
    Dir carry=1;
    Int tIdx=generator.itemCount()-1;
    while(carry!=0&&tIdx>=0)
    {
        Dir tDir= generator[tIdx] + carry;
        generator[tIdx] = (tDir) % (Lattice::l().NeighCount());
        carry = (tDir) / Lattice::l().NeighCount();
        --tIdx;
    }
    //if(carry)cout << "Break!" <<endl;
    return (carry==0);
}
bool LNSMove::skip(block1<Dir,kmm> &generator, Idx pos)
{
    for(Idx tIdx=pos+1;tIdx<generator.itemCount();++tIdx)
    {
        generator[tIdx] = 0;
    }
    N carry=1;
    Int tIdx=pos;
    while(carry!=0&&tIdx>=0)
    {
        //cout<< tIdx << endl;
       // cout <<"In skip:"<< tIdx << " " << (int)generator[tIdx]<<endl;
       Dir tDir= generator[tIdx] + carry;
        generator[tIdx] = (tDir) % Lattice::l().NeighCount();
        carry = (tDir) / Lattice::l().NeighCount();
        //cout << "In next line:" << endl;
        --tIdx;
    }
    //if(carry)cout << "Break!" <<endl;
    return (carry==0);
}

N LNSMove::size()
{
    return aPositions.itemCount();
}
void LNSMove::addPosition(Pos tPos)
{
    Warn(tPos>=Protein::p().Length(), eInvalidIndex);
    //cout << tPos << endl;
    aPositions.insertMem(tPos);
}
bool LNSMove::alreadyAdded(Pos tPos)
{
    for(Idx tIdx=0;tIdx<aPositions.itemCount();++tIdx)
    {
        if(aPositions[tIdx]==tPos) return true;
    }
    return false;
}
void LNSMove::reset()
{
    aPositions.clear();
}
bool LNSMove::equal(block1<Dir,kmm> &generator, block1<Dir,kmm> &prior)
{
    if(generator.itemCount()!=prior.itemCount()) return false;

    for(Idx tIdx =0; tIdx < generator.itemCount();++tIdx)
    {
        if(generator[tIdx]!=prior[tIdx]) return false;
    }
    return true;

}
closePlatypusSpace
