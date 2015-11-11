
#include "tmove.hpp"
#include "emodel.hpp"

tMove::tMove():mPos(0),mDir(0),mDists(Emodel::ModelCount())
{

}

tMove::tMove(Z pos, Dir dir, R dist):mPos(pos),mDir(dir),mDists(Emodel::ModelCount())
{
    mDists[0] = dist;
}
tMove::~tMove()
{

}

tMove::tMove(const tMove & that):mPos(that.mPos),mDir(that.mDir),mDists(that.mDists)
{

}
const tMove & tMove::operator = (const tMove & that)
{
        if(this!=&that)
        {

            mPos=that.mPos;
            mDir=that.mDir;
            mDists=that.mDists;
        }
        return *this;
}
Bln tMove::operator==(const tMove & that) const
{
        return mPos==that.mPos&&mDir==that.mDir&&mDists[0]==that.mDists[0];
}
Bln tMove::operator < (const tMove & that) const
{

        return (mDists[0]<that.mDists[0]);
}
Bln tMove::operator <= (const tMove & that) const
{

        return (mDists[0] <= that.mDists[0]);
}
Bln tMove::operator > (const tMove & that) const
{

        return (mDists[0]>that.mDists[0]);
}
Bln tMove::operator >= (const tMove & that) const
{

        return (mDists[0]>=that.mDists[0]);
}
Z tMove::getPos()
{
        return mPos;
}
Dir tMove::getDir()
{
        return mDir;
}
R tMove::getDist()
{
    Warn(Emodel::ModelCount() == 0, eInvalidComp);
    return mDists[0];
}
R tMove::getDist(Idx idx)
{
    Warn(idx>=Emodel::ModelCount(),eInvalidComp);
    return mDists[idx];
}
void tMove::setPos(Z pos)
{
    mPos=pos;
}
void tMove::setDir(Dir dir)
{
    mDir=dir;
}

void tMove::setDist(R dist)
{
    Warn(Emodel::ModelCount() == 0, eInvalidComp);
    mDists[0]=dist;
}

void tMove::setDist(R dist, Idx idx)
{
    Warn(idx>=Emodel::ModelCount(),eInvalidComp);
    mDists[idx]=dist;
}
