/*!
	@file inits/randomvalid.cpp
	@brief The implementation file for RandomStructured class.
	@details This is the implementation file for RandomStructured class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#include "pspl/inits/randomstructured.hpp"
#include "pspl/globals/conf.hpp"
#include "pspl/globals/change.hpp"
#include "pspl/globals/lattice.hpp"



openPlatypusSpace



/*!
	The costructor.
*/
RandomStructured::RandomStructured(Rnd & theRnd) : mRnd(theRnd)
{
    WatchError
    // do nothing
    CatchError
}


/*!
	The Destructor
*/
RandomStructured::~RandomStructured()
{
    WatchError
    // do nothing
    CatchError
}


/*!
	The Dupllicator
*/
RandomStructured::RandomStructured(RandomStructured const & that) :
	Init(that), mRnd(that.mRnd)
{
    WatchError
    Alert(&that, eUndefDuplicator);
    CatchError
}


/*!
	The assigner.
*/

RandomStructured const & RandomStructured::operator= (RandomStructured const & that)
{
    WatchError
	Alert(&that, eUndefAssigner);
    return *this;
    CatchError
}


/*!
	Compute and store a possible initialisation change.
*/
void RandomStructured::compute(Conf & theConf)
{
    //WatchError
    Z ran,ranf=0,pTurn=0,pUp=0,up,turn,mv,mv1,mvUp,mvP,mvP_,mv1_,mvS,pUpRate = 20,pTuRate = 15;

    Z neighbors= Lattice::l().NeighCount();
    Idx dimension = Space::Dimen();

	Z *forbiden=new Z[Lattice::l().NeighCount()+1];


	// WE FIRSTLY INITIALIZE THE COORDINATES
	for (Cmp i=1;i<=neighbors;i++)
		if(((i-neighbors/dimension-1)%(neighbors/dimension)==0 && i>=(neighbors/dimension-1)))
		{
			forbiden[i]=1;
			forbiden[i+1]=1;
		}


    Point tPoint(0);

    if(Protein::p().Length()%2!=0) tPoint.Comp(dimension-1)=Protein::p().Length()+1;

    mChange.add(0,tPoint);

	do{
		ran=uniform(mRnd,1,neighbors);
	}while(forbiden[ran]==1);

	mv=ran;
	mv1=ran;

	if(mv1%2!=0)
		mv1_=mv1+1;
	else
		mv1_=mv1-1;

	updateForbiden(forbiden,mv1);

	do{
		//		ran=(rand()%12)+1;
		ran=uniform(mRnd,1,neighbors);
	}while(forbiden[ran]==1);
	mvP=ran;

	if(mvP%2!=0) mvP_=mvP+1;
	else mvP_=mvP-1;

	forbiden[mvP]=1;

	forbiden[mvP_]=1;

	do{
		ran=uniform(mRnd,1,neighbors);
	}while(forbiden[ran]==1);
	mvUp=ran;



	for(Cmp i=1;i<Protein::p().Length();++i)
	{

		ranf=uniform(mRnd,(N)0,(N)100);
		turn=0;
		up=0;
		if(pTurn>ranf){
                ranf=uniform(mRnd,(N)0,(N)100);
                if(pUp>ranf){
				mv=mvUp;
				pUp=0;
				up=1;
				pTurn=0;
			}else{
				pUp=pUp+pUpRate;
				mv=mvP;
				pTurn=0;
				turn=1;
			}
		}else{
			pTurn+=pTuRate;
		}

		// transform mv into our own coding for vectors
		Point mTrial = tPoint+ Lattice::l().DirVec(transformMoves[mv]);

        mChange.add(i,mTrial);
        tPoint = mTrial;
        if(up==1)
  		{
  			mv=mv1_;
  			//swapp moves
  			mvS=mv1;
  			mv1=mv1_;
  			mv1_=mvS;
  			mvS=mvP;
  			mvP=mvP_;
  			mvP_=mvS;
  		}
  		else
  		{
  			if(turn==1)
  			{
  				mv=mv1_;
  				//swapp moves
  				mvS=mv1;
  				mv1=mv1_;
  				mv1_=mvS;
  			}

  		}

  }

    delete [] forbiden;


    //CatchError
}



void RandomStructured::updateForbiden(Z *f, Z m)
{

	if((m-1)%4<2)
	{
		for(int i=1;i<=12;i++)
			if(((i-3)%4==0 && i>=3))
			{

				f[i]=1;
				f[i+1]=1;
			}
	}
	else
	{
		for(int i=1;i<=12;i++)
			if(((i)%4==1))
			{
				f[i]=1;
				f[i+1]=1;
			}
	}

	if(m<5)
	{
		for(int i=1;i<=4;i++)
			f[i]=1;
	}
	else if(m<9)
	{
		for(int i=4;i<=8;i++)
    		  f[i]=1;
	}
	else if(m<13)
	{
		for(int i=8;i<=12;i++)
			f[i]=1;
	}

}




closePlatypusSpace
