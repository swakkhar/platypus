/*!
	@file moves/RotationMove.cpp
	@brief The implementation file for RotationMove class.
	@details This is the implementation file for RotationMove class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 05.08.2012 QRL NICTA www.nicta.com.au
*/

#include "pspl/moves/rotationmove.hpp"
#include "pspl/globals/conf.hpp"
#include "pspl/globals/change.hpp"
#include "pspl/globals/lattice.hpp"


openPlatypusSpace


/*!
The costructor
*/
RotationMove::RotationMove()
{
    WatchError

    // do nothing
    CatchError
}


/*!

The Destructor

*/
RotationMove::~RotationMove()
{
    WatchError

    // do nothing

    CatchError
}


/*!

The Dupllicator
*/

RotationMove::RotationMove(RotationMove const & that):Move(that)
{
    WatchError

    Alert(&that, eUndefDuplicator);

    CatchError
}


/*!

Assignment Operator
*/

RotationMove const & RotationMove::operator= (RotationMove const & that)
{
    WatchError

    Move::operator=(that);
    return *this;


    CatchError
}


/*!

The compute method
save all possible move results into mChanges

In the case of the rotation move thePos is a dummy
variable, it contains another internal array
aPoints that contains all the variables; it is
assumed that the array is sorted if the
array contains any element that is not valid or
any element twice it will crash, for the time
being we can trust the caller


dont allow the call with initial position


step 1: calculates the relative string / relative direction array
step 2: for each of the positions, alter and generate a new string
step 3: for each of the strings if the string is a valid one then add the points to mChanges

*/
void RotationMove::compute(Conf & theConf, Pos const thePos)
{
    WatchError
    mChanges.clear();

    /// this version of rotation move will have only one position as rotation
    /// uses the absolute vectors and the encoding with partial maps
    block1<Dir,kmm> mAbsDirs(Protein::p().Span());

    for(Cmp tCmp=0;tCmp<Protein::p().Span();++tCmp)
    {
        Vector absVec= theConf.Coord(tCmp+1)-theConf.Coord(tCmp);
        // now we have got the absVec
        mAbsDirs[tCmp]=Lattice::l().Direction(absVec);
        //cout << (int)mAbsDirs[tCmp] << " #"<<tCmp<<" ";
    }
    /// now we have got the absolute vectors

    for(Dir tDir=0;tDir<Lattice::l().NeighCount();++tDir)
    {
        if(mAbsDirs[thePos]==tDir) continue;
        Vector replaceVector=Lattice::l().DirVec(tDir);
        Vector absVec = Lattice::l().DirVec(mAbsDirs[thePos]);

        /// we need a two dimensional array for maps 4 x 3
        typedef block1<int,kmm> list;
        rackl<list,kmm> rackMap(4);
        for(Idx tIdx=0;tIdx<4;++tIdx)
        rackMap[tIdx]=list(3);

        Idx nonzeroCount=0;
        /// prepare map
        for(Cmp tCmp=0;tCmp<Space::Dimen();++tCmp)
        {
            /// for each of the components of replace Vector
            if(absVec.Comp(tCmp)==0) /// if zero then find zero
            {
                for(Cmp tCmp2=0;tCmp2<Space::Dimen();++tCmp2)
                {
                    if(replaceVector.Comp(tCmp2)==0)
                    {
                        /// assign this position
                        rackMap[0][tCmp]=tCmp2+1;
                        rackMap[1][tCmp]=tCmp2+1;
                        rackMap[2][tCmp]=-(tCmp2+1);
                        rackMap[3][tCmp]=-(tCmp2+1);
                        break;
                    }
                }
            }
            else /// if non-zero
            {
                Idx srl=0;
                for(Cmp tCmp2=0;tCmp2<Space::Dimen();++tCmp2)
                {
                    if(replaceVector.Comp(tCmp2)!=0)
                    {
                        Z sign = replaceVector.Comp(tCmp2) * absVec.Comp(tCmp);
                        rackMap[(srl+nonzeroCount)%4][tCmp]=(tCmp2+1)*sign;
                        rackMap[(srl+2+nonzeroCount)%4][tCmp]=(tCmp2+1)*sign;
                        srl++;
                    }

                }
                nonzeroCount++;
            }
        }
        /// print map???
        ///for(Idx tIdx=0;tIdx<4;++tIdx)
         ///  cout << rackMap[tIdx][0]<< " "<<rackMap[tIdx][1]<< " "<<rackMap[tIdx][2]<< endl;


        for(Idx tIdx=0;tIdx<4;++tIdx) /// for all maps
        {
            hset<Point,xmmh> tPointsDict(Protein::p().Length());
            tPointsDict.clear();
            Pos tPos=0;
            while(tPos<=thePos)
            {
                tPointsDict.insertBllMem(theConf.Coord(tPos));
                tPos++;
            }
            /// already added all the points prior to the pivot position
            /// now change all the points or the vectors according to
            /// the rotation angle / matrix saved in the map
            Point tPoint = theConf.Coord(tPos-1);
            Change tChange;
            tChange.reset();
            while(tPos<Protein::p().Length())
            {
                Vector tAbsVec = Lattice::l().DirVec(mAbsDirs[tPos-1]);
                /// transform the absolute vector using the map
                Vector newVec;
                //cout <<tPos<<" "<< tAbsVec.Comp(0)<<" "<<tAbsVec.Comp(1)<<" "<<tAbsVec.Comp(2);

                for(Cmp tCmp=0;tCmp<Space::Dimen();++tCmp)
                {
                    Z sign = abs(rackMap[tIdx][tCmp])/rackMap[tIdx][tCmp];
                    Cmp tCmpMap= abs(rackMap[tIdx][tCmp])-1;
                    newVec.Comp(tCmpMap)=sign*tAbsVec.Comp(tCmp);
                }
                //cout << newVec.Comp(0)<<" "<<newVec.Comp(1)<<" "<<newVec.Comp(2);
                tPoint=tPoint+newVec;
                //cout <<tPoint.Comp(0)<<" "<<tPoint.Comp(1)<<" "<<tPoint.Comp(2)<<endl;
                if(tPointsDict.findBll(tPoint))
                {
                    ///
                    break;
                }
                tPointsDict.insertBllMem(tPoint);
                tChange.add(tPos,tPoint);
                ++tPos;
            }
            /// if tPos == span() then add change to the set
            if(tPos == Protein::p().Length())
            {
                mChanges.insertMem(tChange);
            }

        }


    }

    CatchError
}
/*{
    WatchError
    mChanges.clear();



    // first calculate absolute dirs

    block1<Dir,kmm> mAbsDirs(Protein::p().Span());
    block1<Dir,kmm> mRelDirs(Protein::p().Span());

    Vector absVec;
    for(Cmp tCmp=0;tCmp<Protein::p().Span();++tCmp)
    {
        absVec= theConf.Coord(tCmp+1)-theConf.Coord(tCmp);

        // now we have got the absVec
        mAbsDirs[tCmp]=Lattice::l().Direction(absVec);
        //cout << (int)mAbsDirs[tCmp] << " #"<<tCmp<<" ";
    }


    // declare a base Matrix dimension == dimension of the space
    Spc matrixDim=Space::Dimen();
    Matrix base(matrixDim);

    // initialize the base from the first element from the lattice model

    base=Lattice::l().DirMat(0);

    for(Cmp tCmp=0;tCmp<Protein::p().Span();++tCmp)
    {
        Dir relMove=Lattice::l().AbsToRelDir(base,mAbsDirs[tCmp]);
        mRelDirs[tCmp]=relMove;
        //cout << (int)relMove << " #"<<tCmp<<" ";
        Lattice::l().updtAbsToRelBase(base,relMove);
    }
    // now the relDirs is ready
    block1<Dir,kmm> generator(aPositions.itemCount());
    for(Idx tIdx=0;tIdx<aPositions.itemCount();++tIdx)
        generator[tIdx]=0; // initialize by zero

    //this is the generator string - needed for multimutation
    //cout <<"rotation size:"<<aPositions.itemCount()<<endl;
    // create a dictionary using all the positions
    hset<Pos,xmmh> posDict;
    for(Idx tIdx=0;tIdx<aPositions.itemCount();++tIdx)
    {
        posDict.insertBllMem(aPositions[tIdx]);
    }
    cout << "pos:" << aPositions[0] << " "<<aPositions.itemCount() ;

    //
    bool tSkip=true;
    while(1)
    {
        // we need to skip the all zero pos
        if(!next(generator))
            break;

        // in this loop we use the generator to create new strings
        // and keep track of the failures
        hset<Point,xmmh> tempDict;
        tempDict.clear();

        // declare a base Matrix dimension == dimension of the space
        Spc matrixDim=Space::Dimen();
        Matrix base(matrixDim);

        // initialize the base from the first element from the lattice model
        Dir absMove=0;//Lmodel::lm().DirVec(0);
        Change tChange;
        tChange.reset();
        Point cur=0;
        tempDict.insertBllMem(cur);
        tChange.add(0,cur);
        base=Lattice::l().DirInv(0);
        Idx genIdx=0;
        tSkip=true;
        //cout << endl;
        for(Cmp tCmp=0;tCmp<Protein::p().Span();++tCmp)
        {
            Vector absVec;
            if(posDict.findBll(tCmp))
            {

                Lattice::l().updtRelToAbsBase((mRelDirs[tCmp]+
                                               generator[genIdx])%Lattice::l().NeighCount(),base);
                ///++genIdx; /// should it be here?

            }
            else
                Lattice::l().updtRelToAbsBase(mRelDirs[tCmp],base);
            //absVec=Lattice::l().RelToAbsDir(absMove,base);
            //cout << (int)Lattice::l().RelToAbsDir(absMove,base) << " #" <<tCmp << " ";
            Dir absDir=Lattice::l().RelToAbsDir(absMove,base);
            absVec=Lattice::l().DirVec(absDir);
            mAbsDirs[tCmp]=absDir;
            if(absDir!=mAbsDirs[theConf.getMatchPos(tCmp)])
            {
                // error
                tSkip=skip(generator,genIdx);
                break;
            }
            if(posDict.findBll(tCmp)) genIdx++;
            cur = cur + absVec;
            // if cur already exists then error
            if(tempDict.findBll(cur))
            {
                // error
                tSkip=skip(generator,genIdx);
                break;
            }

            tChange.add(tCmp+1,cur);
            tempDict.insertBllMem(cur);


        }
        if(!tSkip)
        {
            break;
        }

        if(tChange.size()==Protein::p().Length())
        {

            mChanges.insertMem(tChange);
        }
    }
    cout << " changesize:"<<mChanges.itemCount();
     cout << endl;
	CatchError
}

bool RotationMove::next(block1<Dir,kmm> &generator)
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
bool RotationMove::skip(block1<Dir,kmm> &generator, Idx pos)
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


void RotationMove::addPosition(Pos tPos)
{
    Warn(tPos>=Protein::p().Length(), eInvalidIndex);
    //cout << tPos << endl;
    aPositions.insertMem(tPos);
}
bool RotationMove::alreadyAdded(Pos tPos)
{
    for(Idx tIdx=0;tIdx<aPositions.itemCount();++tIdx)
    {
        if(aPositions[tIdx]==tPos) return true;
    }
    return false;
}
void RotationMove::reset()
{
    aPositions.clear();
}
bool RotationMove::equal(block1<Dir,kmm> &generator, block1<Dir,kmm> &prior)
{
    if(generator.itemCount()!=prior.itemCount()) return false;

    for(Idx tIdx =0; tIdx < generator.itemCount();++tIdx)
    {
        if(generator[tIdx]!=prior[tIdx]) return false;
    }
    return true;

}*/
closePlatypusSpace
