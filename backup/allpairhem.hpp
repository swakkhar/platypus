#ifndef AllpairHEmHhIncluded
#define AllpairHEmHhIncluded


#include "emodel.hpp"


class AllPairHEm: public Emodel
{
	private:
		Dim mAlpha;

	public:
		virtual Typ Type(Ltr const theLtr) const;
		virtual Eng Level(Typ const Typ1, Typ const Typ2) const;
		virtual Eng Level(Conf const & conf, Pos const pos1, Pos const pos2) const;
        virtual Eng Level(Pos const pos1, Point const point1, Pos const pos2, Point const point2) const;

		AllPairHEm();
		AllPairHEm(Dim const Alpha);
		AllPairHEm(AllPairHEm const & that);
		AllPairHEm const & operator = (AllPairHEm const & that);
		~AllPairHEm();
};


inline AllPairHEm::AllPairHEm() : Emodel(2), mAlpha(2)
{
	//	nothing to be done.
}


inline AllPairHEm::AllPairHEm(Dim const Alpha) : Emodel(2), mAlpha(Alpha)
{
	Warn(!mAlpha, eEmptyAlpha);
}


inline AllPairHEm::AllPairHEm(AllPairHEm const & that) : Emodel(that), mAlpha(that.mAlpha)
{
	//	nothing to be done.
}


inline AllPairHEm const & AllPairHEm::operator = (AllPairHEm const & that)
{
	Emodel::operator=(that);
	mAlpha = that.mAlpha;
	return *this;
}


inline AllPairHEm::~AllPairHEm()
{
	//	nothing to be done.
}


#endif //HpEmHhIncluded

