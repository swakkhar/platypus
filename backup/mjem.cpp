#include "conf.hpp"
#include "protein.hpp"
#include "lmodel.hpp"
#include "mjem.hpp"


Typ MjEm::Type(Ltr const theLtr) const
{
	switch(mAlpha)
	{
		case 20:
			if(ALL_PAIR_TRANSLATION_TABLE[theLtr-'A'])
                return ALL_PAIR_TRANSLATION_TABLE[theLtr-'A'];
            else
                Throw(eUnknownAlpha);
		default:
			Throw(eUnknownAlpha);
	}
}

Eng MjEm::Level(Typ const Typ1, Typ const Typ2) const
{
	Warn(Typ1 >= mSize || Typ2 >= mSize, eInvalidType);
	if(Typ1 < Typ2)  return MJmatrix[Typ1][Typ2];
    else return MJmatrix[Typ1][Typ2];
}


Eng MjEm::Level(Conf const & conf, Pos const pos1, Pos const pos2) const
{
    if((pos1 >= pos2 && pos1 <=pos2+1 )||(pos2 >= pos1 && pos2 <=pos1+1)) return 0;

    if(Lmodel::lm().areNeighbors(conf[pos1],conf[pos2]))
        return Level(Protein::getType(pos1),Protein::getType(pos2));
    else return 0;
}
Eng MjEm::Level(Pos const pos1, Point const point1, Pos const pos2, Point const point2) const
{
    if((pos1 >= pos2 && pos1 <=pos2+1 )||(pos2 >= pos1 && pos2 <=pos1+1)) return 0;

    if(Lmodel::lm().areNeighbors(point1,point2))
        return Level(Protein::getType(pos1),Protein::getType(pos2));
    else return 0;
}

