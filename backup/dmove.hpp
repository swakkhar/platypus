#ifndef DmoveHppIncluded
#define DmoveHppIncluded

#include "move.hpp"
class DMove:public Move
{
    public:

    DMove();
    ~DMove();
    DMove(DMove const &);
    DMove const & operator=(DMove const & that);


    virtual void execute(Conf & conf, Pos const pos, Dir const dir)const;

    virtual void simulate(Conf const & conf, Pos const pos, Dir const dir, kblock<R> & deltaFitness)const; // returns new fitness
    virtual R simulate(Conf const & conf, Pos const pos, Dir const dir,Idx const modelIdx=0)const; // returns new fitness

    virtual Bln feasible(Conf const & conf, Pos const pos, Dir const dir)const;

};

#endif
