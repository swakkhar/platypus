#include "conf.hpp"
#include "protein.hpp"
#include "lmodel.hpp"
#include "allpairhem.hpp"


Typ AllPairHEm::Type(Ltr const theLtr) const
{
	switch(mAlpha)
	{
		case 2:
			switch(theLtr)
			{
				case 'P' : return 0;
				case 'H' : return 1;
				default : Throw(eUnknownLetter);
			}
		default:
			Throw(eUnknownAlpha);
	}
}

Eng AllPairHEm::Level(Typ const Typ1, Typ const Typ2) const
{
	Warn(Typ1 >= mSize || Typ2 >= mSize, eInvalidType);
	return (Typ1 & Typ2);
}

Eng AllPairHEm::Level(Conf const & conf, Pos const pos1, Pos const pos2) const
{
    if((pos1 >= pos2 && pos1 <=pos2+1 )||(pos2 >= pos1 && pos2 <=pos1+1)) return 0;

    if(Protein::getType(pos1) ==1 && Protein::getType(pos2) == 1)
    {
        R dist = Tuple::sqrdist(conf[pos1],conf[pos2]);
        return (dist-2)*(dist-2);
    }
    else return 0;
}
Eng AllPairHEm::Level(Pos const pos1, Point const point1, Pos const pos2, Point const point2) const
{
    if((pos1 >= pos2 && pos1 <=pos2+1 )||(pos2 >= pos1 && pos2 <=pos1+1)) return 0;

    if(Protein::getType(pos1) ==1 && Protein::getType(pos2) == 1)
    {
        R dist = Tuple::sqrdist(point1,point2);
        return (dist-2)*(dist-2);
    }
    else return 0;
}

