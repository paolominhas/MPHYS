#include "ScintAna.hh"


ScintAna::ScintAna() : nHits(0)
{
	rndm = new TRandom3();
}

ScintAna::ScintAna(const ScintAna& ana)
{
	nHits = ana.nHits;
	side = ana.side;
	layer = ana.layer;
	bar = ana.bar;
	nPh1 = ana.nPh1;
	nPh2 = ana.nPh2;
	nPh3 = ana.nPh3;
	nPh4 = ana.nPh4;
	nPh = ana.nPh;
	Pos = ana.Pos;
	
	rndm = new TRandom3();
}

const ScintAna& ScintAna::operator=(const ScintAna& ana)
{
	nHits = ana.nHits;
	side = ana.side;
	layer = ana.layer;
	bar = ana.bar;
	nPh1 = ana.nPh1;
	nPh2 = ana.nPh2;
	nPh3 = ana.nPh3;
	nPh4 = ana.nPh4;
	nPh = ana.nPh;
	Pos = ana.Pos;

	return *this;
}

void ScintAna::response(ScintHit *sc)
{
	nHits=sc->nHits;
	side = sc->side;
	layer = sc->layer;
	bar = sc->bar;
	for(int i=0; i<nHits; i++){
		double x1 = length/2.-sc->Pos_X[i];
		double x2 = length/2.+sc->Pos_X[i];
		double y1 = width/4.-sc->Pos_Y[i];
		double y2 = width/4.+sc->Pos_Y[i];
		double r1 = sqrt(pow(y1,2)+pow(sc->Pos_Z[i],2));
		r1 = (r1>width/2.) ? width/2. : r1;
		double r2 = sqrt(pow(y2,2)+pow(sc->Pos_Z[i],2));
		r2 = (r2>width/2.) ? width/2. : r2;

		double fr1 = A[0]*exp(-r1/A[1])+A[2];
		double fr2 = A[0]*exp(-r2/A[1])+A[2];
		fr1=rndm->Gaus(fr1,A[3]);
		fr2=rndm->Gaus(fr2,A[3]);

		double fx1 = B[0]*exp(-x1/B[1]);
		double fx2 = B[0]*exp(-x2/B[1]);
		fx1=rndm->Gaus(fx1,B[2]);
		fx2=rndm->Gaus(fx2,B[2]);

		double nph=sc->Edep[i]*phYield;
		nph=rndm->Gaus(nph,sqrt(nph));

		nPh1.push_back(nph*fr1*fx1*rndm->Gaus(C[0],C[1]));
		nPh2.push_back(nph*fr2*fx1*rndm->Gaus(C[0],C[1]));
		nPh3.push_back(nph*fr1*fx2*rndm->Gaus(C[0],C[1]));
		nPh4.push_back(nph*fr2*fx2*rndm->Gaus(C[0],C[1]));

		}
}

void ScintAna::add_up()
{
	for(int i=0; i<nHits-1; i++){
		for (int j=i+1; j<nHits; j++){
			if(side[i]==side[j] && layer[i]==layer[j] && bar[i]==bar[j]){
				nPh1[i] += nPh1[j];
				nPh2[i] += nPh2[j];
				nPh3[i] += nPh3[j];
				nPh4[i] += nPh4[j];
				erase(j);
				j--;
			}
		}
	}
	for(int i=0; i<nHits; i++){
		nPh.push_back(nPh1[i]+nPh2[i]+nPh3[i]+nPh4[i]);
		Pos.push_back(log((nPh1[i]+nPh2[i])/(nPh3[i]+nPh4[i]))*A[1]/2.);
	}
}

void ScintAna::sort()
{
	for(int i=0; i<nHits-1; i++){
		for (int j=i+1; j<nHits; j++){
			if(nPh[i]<nPh[j]){
				swap(i,j);
			}
		}
	}
}


void ScintAna::sumE()
{
	for(int i=0; i<nSides; i++){
		mult.push_back(0);
		dE.push_back(0);
		Esum1.push_back(0);
		Esum2.push_back(0);
		Pos1.push_back(0);
	}
	std::vector<std::vector<bool>> multArr;
	std::vector<bool> multVec;
	for(int i=0; i<nSides; i++){
		for(int i=0; i<nLayers; i++){
			multVec.push_back(0);
		}
		multArr.push_back(multVec);
	}
	for(int i=0; i<nHits; i++){
		Esum1[side[i]-1]+=nPh[i];
		if(layer[i]==0){
			dE[side[i]-1]+=nPh[i];
			Pos1[side[i]-1] = (Pos1[side[i]-1]*(Esum1[side[i]-1]-nPh[i])+Pos[i]*nPh[i])/(Esum1[side[i]-1]);
		}
		else{
			Esum2[side[i]-1]+=nPh[i];
		}
		if(multArr[side[i]-1][layer[i]]==0 && nPh[i]>threshold*phYield){
			mult[side[i]-1]++;
			multArr[side[i]-1][layer[i]]=1;
		}
	}
	multVec.clear();
	multArr.clear();
}


void ScintAna::erase(int i)
{
	side.erase(side.begin()+i);
	layer.erase(layer.begin()+i);
	bar.erase(bar.begin()+i);
	nPh1.erase(nPh1.begin()+i);
	nPh2.erase(nPh2.begin()+i);
	nPh3.erase(nPh3.begin()+i);
	nPh4.erase(nPh4.begin()+i);
	//E.erase(E.begin()+i);
	//Pos.erase(Pos.begin()+i);
	nHits--;
}

void ScintAna::swap(int i, int j)
{
	iter_swap(side.begin()+i,side.begin()+j);
	iter_swap(layer.begin()+i,layer.begin()+j);
	iter_swap(bar.begin()+i,bar.begin()+j);
	iter_swap(nPh1.begin()+i,nPh1.begin()+j);
	iter_swap(nPh2.begin()+i,nPh2.begin()+j);
	iter_swap(nPh3.begin()+i,nPh3.begin()+j);
	iter_swap(nPh4.begin()+i,nPh4.begin()+j);
	iter_swap(nPh.begin()+i,nPh.begin()+j);
	iter_swap(Pos.begin()+i,Pos.begin()+j);
}


void ScintAna::clear()
{
	nHits=0;
	side.clear();
	layer.clear();
	bar.clear();
	nPh1.clear();
	nPh2.clear();
	nPh3.clear();
	nPh4.clear();
	nPh.clear();
	Pos.clear();
	mult.clear();
	dE.clear();
	Esum1.clear();
	Esum2.clear();
	Pos1.clear();
}

ClassImp(ScintAna)
