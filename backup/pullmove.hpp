#ifndef PullmoveHppIncluded
#define PullmoveHppIncluded

#include "move.hpp"
class PullMove:public Move
{
    private:
        Z pullDir;

    public:

    PullMove();
    PullMove(Z pDir);
    ~PullMove();
    PullMove(PullMove const &);
    PullMove const & operator=(PullMove const & that);


    virtual void execute(Conf & conf, Pos const pos, Dir const dir)const;

    virtual void simulate(Conf const & conf, Pos const pos, Dir const dir, kblock<R> & deltaFitness)const; // returns new fitness
    virtual R simulate(Conf const & conf, Pos const pos, Dir const dir,Idx const modelIdx=0)const; // returns new fitness

    virtual Bln feasible(Conf const & conf, Pos const pos, Dir const dir)const;

};

#endif

