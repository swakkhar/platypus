/*!
	@file ContactTrendObj.cpp
	@brief The implementation file for ContactTrendObj class.
	@details This is the implementation file for ContactTrendObj class.
	@author Swakkhar Shatabda swakkhar.shatabda@nicta.com.au
	@author M.A.Hakim Newton hakim.newton@nicta.com.au
	@date 01.12.2011 QRL NICTA www.nicta.com.au
*/



#include "pspl/functions/contacttrendobj.hpp"



openPlatypusSpace

//int weightArr[]={0,1,1,2,3,2,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0};
Int weightArr[]={1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

/*!
	The constructor
*/
ContactTrendObj::ContactTrendObj(Conf const & theConf, B const NeedHint) :
		Function(theConf.SysHdl(), NeedHint)
{
		WatchError

    //cout << "NeedHint: "<< NeedHint << endl;
	Warn(!(Energy::e().ModelId() != EnergyCls<HpEnergy>::ModelId ||
      Energy::e().ModelId() != EnergyCls<MjEnergy>::ModelId ||
      Energy::e().ModelId() != EnergyCls<BarerraEnergy>::ModelId), eEnergyMismatch);


    Dim tDimen = Space::Dimen();
	Hdl tSysHdl = theConf.SysHdl();
	Dim tSpan = Protein::p().Span();
	Dim tLength = Protein::p().Length();
	// make a toplist from the acid types present in the sequence

	block1<Dim, kmm> types(20);
	for(Idx tIdx=0;tIdx<20;++tIdx)
	{
        types[tIdx]=0;
	}
	for(Idx tIdx=0;tIdx<tLength;++tIdx)
	{
	    ++types[Protein::p().Type(tIdx)];
	}
	/*cout << "the types: ";
	for(Idx tIdx=0;tIdx<20;++tIdx)
	{
	    cout << types[tIdx]<< " ";
	}
	cout << endl;*/
    // now take top
    block1<Dim, kmm> sort(20);
    for(Idx tIdx=0;tIdx<20;++tIdx)
	{
        sort[tIdx]=tIdx;
	}
	for(Idx tIdx1=0;tIdx1<20;++tIdx1)
        for(Idx tIdx2=tIdx1+1;tIdx2<20;++tIdx2)
        {
            if(types[sort[tIdx1]]<types[sort[tIdx2]])
            {
                //swap
                Dim tDim = sort[tIdx1];
                sort[tIdx1]=sort[tIdx2];
                sort[tIdx2]=tDim;
            }
        }

    // now sorted, take top 5 and corresponding pairs
    N numberOfRes=5;
    N numberOfPairs=15;
    block1<Dim,kmm> topList(numberOfRes);
    for(Idx tIdx=0;tIdx<20;++tIdx)
	{
	    if(tIdx<numberOfRes)
            topList[tIdx]=sort[tIdx];
	    else
            types[sort[tIdx]]=0;
	}
	// now types got non zero only in the toplist elements
    Dim tLatticeUnitDist=10*Lattice::l().latticeConverter();

    block1<Prm,nmm> tResults(numberOfPairs); // 5*4/2
	for(Idx tIdx1=0;tIdx1<numberOfRes;++tIdx1)
    {
        //cout << topList[tIdx1] <<endl;
        for(Idx tIdx2=tIdx1;tIdx2<numberOfRes;++tIdx2)
        {
            // for each pairs
            rackl<block1<Prm,xmm>,kmm> rackNodes(cutOffCount);

            Idx yMax=0;

            for(Pos tPos1 = 0; tPos1 < tSpan; ++tPos1)
            {
                Typ tTyp1 = Protein::p().Type(tPos1);

                for(Pos tPos2 = tPos1 + 2; tPos2 < tLength; ++tPos2)
                {
                    Typ tTyp2 = Protein::p().Type(tPos2);
                    if(!((tTyp1==topList[tIdx1]&&tTyp2==topList[tIdx2])||
                       (tTyp1==topList[tIdx2]&&tTyp2==topList[tIdx1])))
                        continue;
                    // else add this 20 different nodes!
                    ++yMax;

                    Prm tDiff[tDimen], tSqrd[tDimen];
                    for(Cmp tCmp = 0; tCmp < tDimen; ++tCmp)
                    {
                        Prm tVar1 = Prm(TermVar,tPos1 * tDimen + tCmp);
                        Prm tVar2 = Prm(TermVar,tPos2 * tDimen + tCmp);
                        tDiff[tCmp] = Prm(TermFunc,BsubXiFeVi::def(Xv, tSysHdl, tVar1, tVar2));
                        tSqrd[tCmp] = Prm(TermFunc,UsqrXiFeVi::def(Xv, tSysHdl, tDiff[tCmp]));
                    }
                    Prm tDist = Prm(TermFunc,SumXiFeVi::def(Xv, tSysHdl, tSqrd, tDimen));
                    Prm modelDist = Prm(TermFunc,UsqrtXiFeVi::def(Xv,tSysHdl,tDist));
                    Prm rDist = Prm(TermFunc,
                        UmultXiFeVi::def(Xv,tSysHdl,
                        modelDist,UmultXiFeVi::bind(tLatticeUnitDist)));

                    // now compare with this rDist making a loop
                    // add add all the resulting nodes
                    for(Idx tIdxDist=0;tIdxDist<cutOffCount;++tIdxDist)
                    {
                        rackNodes[tIdxDist].insertMem(Prm(TermFunc, UlesXiFeVi::def(Xv,tSysHdl,rDist,UlesXiFeVi::bind(cutOff[tIdxDist]*10))));                    }
                    // so now the rack contains all the comparing nodes
                    // now go to the outher loop of the points and add all these nodes
                }
            }
            // here add all the nodes from rack with a loop
            block1<Prm,kmm> tCount(cutOffCount);
            for(Idx tIdxDist=0;tIdxDist<cutOffCount;++tIdxDist)
            {
                Prm tNode = Prm (TermFunc,SumXiFeVi::def(Xv,tSysHdl,rackNodes[tIdxDist].items(),rackNodes[tIdxDist].itemCount()));
                tCount[tIdxDist] = Prm (TermFunc,UmultXiFeVi::def(Xv,tSysHdl,tNode,UmultXiFeVi::bind(10000)));
            }


            /////// PREVIOUS tCount[maxCutoff] == yMAX, now this is the max TYPE COUNT
            block1<Prm,kmm> toSub(cutOffCount); // multiply tCount[20] by the array value
            for(Idx tIdxDist=0;tIdxDist<cutOffCount;++tIdxDist)
            {
                toSub[tIdxDist] = Prm (TermFunc,UmultXiFeVi::def(Xv,tSysHdl,tCount[cutOffCount-1],UmultXiFeVi::bind(mVal[tIdxDist])));
            }

            block1<Prm,kmm> tPairFinal(cutOffCount);
            block1<Prm,kmm> tNorm(cutOffCount);
            for(Idx tIdxDist=0;tIdxDist<cutOffCount;++tIdxDist)
            {
                //tNorm[tIdxDist] = Prm (TermFunc,BdiffXiFeVi::def(Xv,tSysHdl,tCount[tIdxDist],toSub[tIdxDist]));
                tNorm[tIdxDist] = Prm(TermFunc,UdiffXiFeVi::def(Xv,tSysHdl,tCount[tIdxDist],UdiffXiFeVi::bind(mVal[tIdxDist]*yMax)));
                tPairFinal[tIdxDist] = Prm (TermFunc,
                    UdivXiFeVi::def(Xv,tSysHdl,tNorm[tIdxDist],UdivXiFeVi::bind(10000)));
            }
            // now put the sum the of all these PairFinalscore to results
            //tResults.insert(Prm(TermFunc,SumXiFeVi::def(Xv,tSysHdl,tPairFinal.items(),tPairFinal.itemCount()),ValueAsEvalMin));

            // or weighted sum
            tResults.insert(Prm(TermFunc,WghtSumKiXiFeVi::def(Xv,tSysHdl,weightArr,tPairFinal.items(),tPairFinal.itemCount()),ValueAsEvalMin));
        }
    }
    //cerr << "NeedHint:"<< NeedHint << endl;
    if (NeedHint)
	{
		mFuncHdl = SumXiEFcMiHn::def(Xv | EvalMin, tSysHdl, tResults.items(), tResults.itemCount());
		mMetricRec = &SumXiEFcMiHn::refc(tSysHdl, mFuncHdl).MetricRec();
	}
	else
	{
		mFuncHdl = SumXiFcMi::def(Xv, tSysHdl, tResults.items(), tResults.itemCount());
		mMetricRec = &SumXiFcMi::refc(tSysHdl, mFuncHdl).MetricRec();
	}


	CatchError
}


/*!
	The destructor.
*/
ContactTrendObj::~ContactTrendObj()
{
	WatchError
	// nothing to be done.
	CatchError
}



/*!
	The duplicator.
*/
ContactTrendObj::ContactTrendObj(ContactTrendObj const & that) : Function(that)
{
	WatchError
	Alert(&that, eUndefDuplicator);
	CatchError
}



/*!
	The assigner.
*/
ContactTrendObj const & ContactTrendObj::operator = (ContactTrendObj const & that)
{
	WatchError
	Alert(&that, eUndefAssigner);
	return *this;
	CatchError
}
closePlatypusSpace
