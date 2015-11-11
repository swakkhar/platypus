#ifndef MjEmHhIncluded
#define MjEmHhIncluded


#include "emodel.hpp"


class MjEm: public Emodel
{
	private:
		Dim mAlpha;

	public:
		virtual Typ Type(Ltr const theLtr) const;
		virtual Eng Level(Typ const Typ1, Typ const Typ2) const;
		virtual Eng Level(Conf const & conf, Pos const pos1, Pos const pos2) const;
        virtual Eng Level(Pos const pos1, Point const point1, Pos const pos2, Point const point2) const;

		MjEm();
		MjEm(Dim const Alpha);
		MjEm(MjEm const & that);
		MjEm const & operator = (MjEm const & that);
		~MjEm();
};


inline MjEm::MjEm() : Emodel(20), mAlpha(20)
{
	//	nothing to be done.
}


inline MjEm::MjEm(Dim const Alpha) : Emodel(20), mAlpha(Alpha)
{
	Warn(!mAlpha, eEmptyAlpha);
}


inline MjEm::MjEm(MjEm const & that) : Emodel(that), mAlpha(that.mAlpha)
{
	//	nothing to be done.
}


inline MjEm const & MjEm::operator = (MjEm const & that)
{
	Emodel::operator=(that);
	mAlpha = that.mAlpha;
	return *this;
}


inline MjEm::~MjEm()
{
	//	nothing to be done.
}

const R MJmatrix[20][20]=
{
{ -1060, 190, -230, 160, -80, 60, 80, 40, 0, -80, 190, -20, 50, 130, 690, 30, -190, 240, 710, 0 },
{ 190, 40, -420, -280, -200, -140, -670, -130, 250, 190, 190, 140, 460, 80, 440, 650, 990, 310, 0, -340 },
{ -230, -420, -440, -190, -300, -220, -160, 0, 30, 380, 310, 290, 490, 180, 270, 390, -160, 410, 440, 200 },
{ 160, -280, -190, -220, -410, -250, 20, 110, -220, 250, 140, 210, 360, 530, 350, 590, 490, 420, 360, 250 },
{ -80, -200, -300, -410, -270, -290, -90, 240, -10, 230, 200, 250, 260, 300, 430, 670, 160, 350, 190, 420 },
{ 60, -140, -220, -250, -290, -290, -170, 20, -100, 160, 250, 180, 240, 500, 340, 580, 190, 300, 440, 90 },
{ 80, -670, -160, 20, -90, -170, -120, -40, -90, 180, 220, 340, 80, 60, 290, 240, -120, -160, 220, -280 },
{ 40, -130, 0, 110, 240, 20, -40, -60, 90, 140, 130, 90, -200, -200, -100, 0, -340, -250, -210, -330 },
{ 0, 250, 30, -220, -10, -100, -90, 90, -130, -70, -90, -60, 80, 280, 260, 120, 340, 430, 140, 100 },
{ -80, 190, 380, 250, 230, 160, 180, 140, -70, -380, -260, -160, -60, -140, 250, -220, 200, -40, 110, -110 },
{ 190, 190, 310, 140, 200, 250, 220, 130, -90, -260, 30, -80, -140, -110, 0, -290, -190, -350, -90, -70 },
{ -20, 140, 290, 210, 250, 180, 340, 90, -60, -160, -80, 200, -140, -140, -260, -310, -50, 170, -130, 10 },
{ 50, 460, 490, 360, 260, 240, 80, -200, 80, -60, -140, -140, 290, -250, -170, -170, -20, -520, -380, -420 },
{ 130, 80, 180, 530, 300, 500, 60, -200, 280, -140, -110, -140, -250, -530, -320, -300, -240, -140, -330, -180 },
{ 690, 440, 270, 350, 430, 340, 290, -100, 260, 250, 0, -260, -170, -320, -30, -150, -450, -740, -970, -100 },
{ 30, 650, 390, 590, 670, 580, 240, 0, 120, -220, -290, -310, -170, -300, -150, 40, -390, -720, -760, 40 },
{ -190, 990, -160, 490, 160, 190, -120, -340, 340, 200, -190, -50, -20, -240, -450, -390, -290, -120, 220, -210 },
{ 240, 310, 410, 420, 350, 300, -160, -250, 430, -40, -350, 170, -520, -140, -740, -720, -120, 110, 750, -380 },
{ 710, 0, 440, 360, 190, 440, 220, -210, 140, 110, -90, -130, -380, -330, -970, -760, 220, 750, 250, 110 },
{ 0, -340, 200, 250, 420, 90, -280, -330, 100, -110, -70, 10, -420, -180, -100, 40, -210, -380, 110, 260 }
};
/*{
{-1.06,  0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00},
{ 0.19,	 0.04,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00},
{-0.23,	-0.42,	-0.44,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00},
{ 0.16,	-0.28,	-0.19,	-0.22,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00},
{-0.08,	-0.20,	-0.30,	-0.41,	-0.27,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00},
{ 0.06,	-0.14,	-0.22,	-0.25,	-0.29,	-0.29,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00},
{ 0.08,	-0.67,	-0.16,	 0.02,	-0.09,	-0.17,	-0.12,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00},
{ 0.04,	-0.13,	 0.00,	 0.11,	 0.24,	 0.02,	-0.04,	-0.06,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00},
{ 0.00,	 0.25,	 0.03,	-0.22,	-0.01,	-0.10,	-0.09,	 0.09,	-0.13,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00},
{-0.08,	 0.19,	 0.38,	 0.25,	 0.23,	 0.16,	 0.18,	 0.14,	-0.07,	-0.38,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00},
{ 0.19,	 0.19,	 0.31,	 0.14,	 0.20,	 0.25,	 0.22,	 0.13,	-0.09,	-0.26,	 0.03,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00},
{-0.02,	 0.14,	 0.29,	 0.21,	 0.25,	 0.18,	 0.34,	 0.09,	-0.06,	-0.16,	-0.08,	 0.20,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00},
{ 0.05,	 0.46,	 0.49,	 0.36,	 0.26,	 0.24,	 0.08,	-0.20,	 0.08,	-0.06,	-0.14,	-0.14,	 0.29,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00},
{ 0.13,	 0.08,	 0.18,	 0.53,	 0.30,	 0.50,	 0.06,	-0.20,	 0.28,	-0.14,	-0.11,	-0.14,	-0.25,	-0.53,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00},
{ 0.69,	 0.44,	 0.27,	 0.35,	 0.43,	 0.34,	 0.29,	-0.10,	 0.26,	 0.25,	 0.00,	-0.26,	-0.17,	-0.32,	-0.03,	 0.00,	 0.00,	 0.00,	 0.00,	 0.00},
{ 0.03,	 0.65,	 0.39,	 0.59,	 0.67,	 0.58,	 0.24,	 0.00,	 0.12,	-0.22,	-0.29,	-0.31,	-0.17,	-0.30,	-0.15,	 0.04,	 0.00,	 0.00,	 0.00,	 0.00},
{-0.19,	 0.99,	-0.16,	 0.49, 	 0.16, 	 0.19,	-0.12,	-0.34,	 0.34,	 0.20,	-0.19, 	-0.05,	-0.02,	-0.24,	-0.45,	-0.39,	-0.29,	 0.00,	 0.00,	 0.00},
{ 0.24,	 0.31,	 0.41,	 0.42, 	 0.35,	 0.30,	-0.16,	-0.25,	 0.43,	-0.04,	-0.35,	 0.17,	-0.52,	-0.14,	-0.74,	-0.72,	-0.12,	 0.11,	 0.00,	 0.00},
{ 0.71,	 0.00,	 0.44,	 0.36,	 0.19,	 0.44,	 0.22,	-0.21,	 0.14,	 0.11,	-0.09,	-0.13,	-0.38,	-0.33,	-0.97,	-0.76,	 0.22,	 0.75,	 0.25,	 0.00},
{ 0.00,	-0.34,	 0.20,	 0.25,	 0.42,	 0.09,	-0.28,	-0.33,	 0.10,	-0.11,	-0.07,	 0.01,	-0.42,	-0.18,	-0.10,	 0.04,	-0.21, 	-0.38,	 0.11,	 0.26},
};*/

const int ALL_PAIR_TRANSLATION_TABLE[]={8,-1,0,15,14,2,9,16,3,-1,18,4,1,13,-1,19,12,17,11,10,-1,5,6,-1,7,-1};
const int HP_TRANSLATION_TABLE[]={0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,1};


#endif //MjEmHhIncluded

