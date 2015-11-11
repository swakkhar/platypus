#include "conf.hpp"
#include "protein.hpp"
#include "lmodel.hpp"
#include "realbem.hpp"


Typ RealBEm::Type(Ltr const theLtr) const
{
	switch(mAlpha)
	{
		case 20:
			if(realType[theLtr-'A'])
                return realType[theLtr-'A'];
            else
                Throw(eUnknownAlpha);
		default:
			Throw(eUnknownAlpha);
	}
}

Eng RealBEm::Level(Typ const Typ1, Typ const Typ2) const
{
	Warn(Typ1 >= mSize || Typ2 >= mSize, eInvalidType);
	if(Typ1 < Typ2)  return RealB[Typ1][Typ2];
    else return RealB[Typ1][Typ2];
}


Eng RealBEm::Level(Conf const & conf, Pos const pos1, Pos const pos2) const
{
    if((pos1 >= pos2 && pos1 <=pos2+1 )||(pos2 >= pos1 && pos2 <=pos1+1)) return 0;

    if(Lmodel::lm().areNeighbors(conf[pos1],conf[pos2]))
        return Level(Protein::getType(pos1),Protein::getType(pos2));
    else return 0;
}
Eng RealBEm::Level(Pos const pos1, Point const point1, Pos const pos2, Point const point2) const
{
    if((pos1 >= pos2 && pos1 <=pos2+1 )||(pos2 >= pos1 && pos2 <=pos1+1)) return 0;

    if(Lmodel::lm().areNeighbors(point1,point2))
        return Level(Protein::getType(pos1),Protein::getType(pos2));
    else return 0;
}


