#include "dmove.hpp"


DMove::DMove()
{
    // do nothing
}
DMove::~DMove()
{

}

DMove::DMove(DMove const & that):Move(that)
{
    // do nothing
}
DMove const & DMove::operator=(DMove const & that)
{
    Move::operator=(that);
    return *this;
}


void DMove::execute(Conf & pConf, Pos const pos, Dir const dir)const
{

    Pos tPlace = 1;
    if(pos == 0) tPlace = pos +1;
    else tPlace = pos -1;

    Point mTrial = pConf[tPlace] + Lmodel::lm().DirVec(dir);

    for(Idx tIdx =0; tIdx < Emodel::ModelCount();++tIdx)
    {
        R tDelta = 0;
        R oldDist =0;
        R newDist = 0;
        for(Idx tCmp = 0 ; tCmp < Protein::Length(); ++tCmp)
        {
            if(pConf.cacheEnabled(tIdx))
                oldDist = pConf.cacheFitness(tCmp,pos,tIdx);
            else
                oldDist = Emodel::em(tIdx).Level(pConf,tCmp,pos);
            newDist = Emodel::em(tIdx).Level(tCmp,pConf[tCmp],pos,mTrial);
            if(pConf.cacheEnabled(tIdx))
                pConf.cacheFitness(tCmp,pos,tIdx) = newDist;
            tDelta = tDelta + newDist - oldDist;
        }

        pConf.fitness(tIdx)+=tDelta;
    }
    pConf.remove(pConf[pos]);
    pConf[pos]=mTrial;
    pConf.add(mTrial);
}

void DMove::simulate(Conf const & pConf, Pos const pos, Dir const dir, kblock<R> & deltaFitness)const
{
    Pos tPlace = 1;
	if(pos==0) tPlace = pos + 1;
	else tPlace = pos - 1;

    Point mTrial= pConf[tPlace]+ Lmodel::lm().DirVec(dir);


    for(Idx tIdx= 0; tIdx < Emodel::ModelCount(); ++tIdx)
    {
         R newDist = 0;
         R oldDist = 0;
         R tDelta = 0;
        for(Idx tCmp = 0; tCmp < Protein::Length() ; ++tCmp)
        {
            if(pConf.cacheEnabled(tIdx))
                oldDist = pConf.cacheFitness(tCmp,pos,tIdx);
            else
                oldDist = Emodel::em(tIdx).Level(pConf,tCmp,pos);
            newDist = Emodel::em(tIdx).Level(tCmp,pConf[tCmp],pos,mTrial);
            tDelta = tDelta + newDist - oldDist;
        }
        deltaFitness[tIdx]=tDelta;
    }

}

R DMove::simulate(Conf const & pConf, Pos const pos, Dir const dir,Idx const modelIdx)const
{

    Pos tPlace = 1;
	if(pos==0) tPlace = pos + 1;
	else tPlace = pos - 1;

    Point mTrial= pConf[tPlace]+ Lmodel::lm().DirVec(dir);

    R newDist = 0;
    R oldDist = 0;
	R tDelta = 0;

    for(Cmp tCmp = 0; tCmp < Protein::Length() ; ++tCmp)
    {
        if(pConf.cacheEnabled(modelIdx))
                oldDist = pConf.cacheFitness(tCmp,pos,modelIdx);
            else
                oldDist = Emodel::em(modelIdx).Level(pConf,tCmp,pos);
        newDist = Emodel::em(modelIdx).Level(tCmp,pConf[tCmp],pos,mTrial);
        tDelta = tDelta + newDist - oldDist;

    }
    return tDelta;

}

Bln DMove::feasible(Conf const & pConf, Pos const pos, Dir const dir)const
{
    Warn(pos > Protein::Span(), eInvalidComp);

    Pos tPos = 0;
    if(pos == 0) tPos = pos + 1;
    else tPos = pos - 1;

    Point tPoint = pConf[tPos]+Lmodel::lm().DirVec(dir);

    if(pConf.find(tPoint) == true )
        return false;
    else if(pos != Protein::Span() && pos != 0)
    {
        if(Lmodel::lm().areNeighbors(tPoint,pConf[pos + 1]))
            return true;
        else
            return false;
    }
    else
        return true;
}
