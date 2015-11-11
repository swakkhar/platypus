#include "pullmove.hpp"


PullMove::PullMove()
{
    pullDir = 1;
}
PullMove::PullMove(Z pDir)
{
    pullDir = pDir;
}
PullMove::~PullMove()
{

}

PullMove::PullMove(PullMove const & that):Move(that), pullDir(that.pullDir)
{
    // do nothing
}
PullMove const & PullMove::operator=(PullMove const & that)
{
    Move::operator=(that);
    pullDir = that.pullDir;
    return *this;
}


void PullMove::execute(Conf & pConf, Pos const pos, Dir const dir)const
{
    Point l = pConf[pos + pullDir] + Lmodel::lm().DirVec(dir);
    Idx idx = 0;
    while(idx < Lmodel::lm().Count())
    {
        Point t = l + Lmodel::lm().DirVec(idx);
        if(pConf.find(t)==false && Lmodel::lm().areNeighbors(t,pConf[pos]))
                break;
        ++idx;
    }
    Point c = l + Lmodel::lm().DirVec(idx);

    Int start = pos;
    Int itr = pos;
    xblock<Point> pointsT; // temporary point holder
    xset<Int> pullWindow;

    kblock<R> before(2);


    while(itr >=0 && itr < (int)Protein::Length()&& l!=pConf[itr])
    {
        pointsT.annex(l);
        pullWindow.annex(itr);
        l=c;
        c=pConf[itr];
        itr+=-pullDir;
    }

    Int end = itr;


    for(Idx modelIdx = 0 ; modelIdx < Emodel::ModelCount(); ++modelIdx)
    {

        R tDelta = 0;
        Int tStart = start;
        Int xblockItr = 0;
        while(tStart !=end)
        {
            R oldDist=0;
            R newDist=0;
            for(Int tIdx = 0; tIdx < (int)Protein::Length(); ++tIdx )
            {

                oldDist = Emodel::em(modelIdx).Level(pConf,tIdx,tStart);
                if(pullWindow.find(tIdx))
                {
                    if(tIdx<=tStart) newDist=oldDist;
                    else newDist = Emodel::em(modelIdx).Level(tIdx,pointsT[(start-tIdx)*pullDir],tStart,pointsT[xblockItr]);
                }
                else
                {
                    newDist = Emodel::em(modelIdx).Level(tIdx,pConf[tIdx],tStart,pointsT[xblockItr]);
                }
                tDelta = tDelta + newDist - oldDist;
            }

            ++xblockItr;
            tStart+=-pullDir;
        }

        pConf.fitness(modelIdx)+=tDelta;
    }


    // now make the changes
    Int tStart = start;
    Int xblockItr = 0;
    while(tStart != end)
    {
        pConf.remove(pConf[tStart]);
        pConf[tStart] = pointsT[xblockItr];
        pConf.add(pointsT[xblockItr]);
        ++xblockItr;
        tStart+=-pullDir;
    }
    for(Idx modelIdx = 0 ; modelIdx < Emodel::ModelCount(); ++modelIdx)
    {
        if(pConf.calcFitness(modelIdx)!=pConf.fitness(modelIdx))
        {
            exit(0);
        }
    }
}

void PullMove::simulate(Conf const & pConf, Pos const pos, Dir const dir, kblock<R> & deltaFitness)const
{
    Point l = pConf[pos + pullDir] + Lmodel::lm().DirVec(dir);
    Idx idx = 0;
    while(idx < Lmodel::lm().Count())
    {
        Point t = l + Lmodel::lm().DirVec(idx);
        if(pConf.find(t)==false && Lmodel::lm().areNeighbors(t,pConf[pos]))
                break;
        ++idx;
    }
    Point c = l + Lmodel::lm().DirVec(idx);\

    Int start = pos;
    Int itr = pos;
    xblock<Point> pointsT; // temporary point holder
    xset<Int> pullWindow;

    while(itr >=0 && itr < (int)Protein::Length() && l!=pConf[itr])
    {
        pointsT.annex(l);
        pullWindow.annex(itr);
        l=c;
        c=pConf[itr];
        itr+=-pullDir;
    }
    Int end =itr;


    for(Idx modelIdx = 0 ; modelIdx < Emodel::ModelCount(); ++modelIdx)
    {
        R tDelta = 0;
        Int tStart = start;
        Int xblockItr = 0;
        while(tStart !=end)
        {
            R oldDist=0;
            R newDist=0;
            for(Int tIdx = 0; tIdx < (int)Protein::Length(); ++tIdx )
            {

                oldDist = Emodel::em(modelIdx).Level(pConf,tIdx,tStart);
                if(pullWindow.find(tIdx))
                {
                    if(tIdx<=tStart) newDist=oldDist;
                    else newDist = Emodel::em(modelIdx).Level(tIdx,pointsT[(start-tIdx)*pullDir],tStart,pointsT[xblockItr]);
                }
                else
                {
                    newDist = Emodel::em(modelIdx).Level(tIdx,pConf[tIdx],tStart,pointsT[xblockItr]);
                }
                tDelta = tDelta + newDist - oldDist;
            }

            ++xblockItr;
            tStart+=-pullDir;
        }
        deltaFitness[modelIdx] = tDelta;
    }


}

R PullMove::simulate(Conf const & pConf, Pos const pos, Dir const dir,Idx const modelIdx)const
{


    Point l = pConf[pos + pullDir] + Lmodel::lm().DirVec(dir);
    Idx idx = 0;
    while(idx < Lmodel::lm().Count())
    {
        Point t = l + Lmodel::lm().DirVec(idx);
        if(pConf.find(t)==false && Lmodel::lm().areNeighbors(t,pConf[pos]))
                break;
        ++idx;
    }
    Point c = l + Lmodel::lm().DirVec(idx);\

    Int start = pos;
    Int itr = pos;
    xblock<Point> pointsT; // temporary point holder
    xset<Int> pullWindow;


    while(itr >=0 && itr < (int)Protein::Length()&& l!=pConf[itr])
    {
        pointsT.annex(l);
        pullWindow.annex(itr);
        l=c;
        c=pConf[itr];
        itr+=-pullDir;
    }
    Int end =itr;
    Int xblockItr = 0;
    R tDelta = 0;
    Int tStart = start;
    while(tStart !=end)
    {
        R oldDist=0;
        R newDist=0;
        for(Int tIdx = 0; tIdx < (int)Protein::Length(); ++tIdx )
        {

            oldDist = Emodel::em(modelIdx).Level(pConf,tIdx,tStart);
            if(pullWindow.find(tIdx))
            {
                if(tIdx<=tStart) newDist=oldDist;
                else newDist = Emodel::em(modelIdx).Level(tIdx,pointsT[(start-tIdx)*pullDir],tStart,pointsT[xblockItr]);
            }
            else
            {
                newDist = Emodel::em(modelIdx).Level(tIdx,pConf[tIdx],tStart,pointsT[xblockItr]);
            }
            tDelta = tDelta + newDist - oldDist;
        }

        ++xblockItr;
        tStart+=-pullDir;
    }
    return tDelta;




}

Bln PullMove::feasible(Conf const & pConf, Pos const pos, Dir const dir)const
{
    if(pos == 0|| pos == Protein::Span())
        return false;

    Point l = pConf[pos + pullDir] + Lmodel::lm().DirVec(dir);

    if(pConf.find(l)) return false;
    else
    {
        //find a c that is neighbor to pos
        Idx idx=0;
        while(idx < Lmodel::lm().Count())
        {
            Point c = l + Lmodel::lm().DirVec(idx);
            if(pConf.find(c)==false && Lmodel::lm().areNeighbors(c,pConf[pos]))
                return true;
            ++idx;
        }
        return false;
    }
}

