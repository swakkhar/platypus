
#ifndef TmoveHppIncluded
#define TmoveHppIncluded

#include "global.hh"
class tMove
{
        private:
            Pos mPos; // it is Z since we need -1 for invali
            Dir mDir;
            kblock<R> mDists;

        public:
            tMove();
            ~tMove();
            tMove(const tMove & that);
            tMove(Z pos, Dir dir, R dist);

            /*!
                getters and setters
            */
            Z getPos();
            Dir getDir();

            R getDist(); // at zero
            R getDist(Idx idx);


            void setPos(Z pos);
            void setDir(Dir dir);

            void setDist(R dist); // at zero
            void setDist(R dist, Idx idx);

            const tMove & operator =( const tMove & that);

            Bln operator ==(const tMove & other) const;
            Bln operator <(const tMove & other) const;
            Bln operator <=(const tMove & other) const;
            Bln operator >(const tMove & other) const;
            Bln operator >=(const tMove & other) const;
};

#endif
