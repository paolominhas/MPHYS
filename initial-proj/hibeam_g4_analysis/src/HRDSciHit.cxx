#include "HRDSciHit.hh"

HRDSciHit::HRDSciHit() : nHits(0)
{
	rndm = new TRandom3();
}

HRDSciHit::HRDSciHit(const HRDSciHit& hit)
{
	nHits = hit.nHits;
	Pos_X = hit.Pos_X;
	Pos_Y = hit.Pos_Y;
	Pos_Z = hit.Pos_Z;
	Edep = hit.Edep;
	Time = hit.Time;
	side = hit.side;
	layer = hit.layer;
	bar = hit.bar;
	track = hit.track;
	
	rndm = new TRandom3();
}

const HRDSciHit& HRDSciHit::operator=(const HRDSciHit& hit)
{
	nHits = hit.nHits;
	Pos_X = hit.Pos_X;
	Pos_Y = hit.Pos_Y;
	Pos_Z = hit.Pos_Z;
	Edep = hit.Edep;
	Time = hit.Time;
	side = hit.side;
	layer = hit.layer;
	bar = hit.bar;
	track = hit.track;

	return *this;
}

bool HRDSciHit::operator==(const HRDSciHit& hit) const { return (track == hit.track && side == hit.side && layer == hit.layer && bar == hit.bar); }

void HRDSciHit::push_back(int tr, int b, int l, int s, double e, double t, double x, double y, double z)
{
	nHits++;
	track.push_back(tr);
	side.push_back(s);
	layer.push_back(l);
	bar.push_back(b);
	Edep.push_back(e);
	Time.push_back(t);
	Pos_X.push_back(x);
	Pos_Y.push_back(y);
	Pos_Z.push_back(z);
}

void HRDSciHit::process()
{
	sum_up();
	response();
	sort();
}

void HRDSciHit::sum_up()
{
	for(int i=0; i<nHits-1; i++){
		for (int j=i+1; j<nHits; j++){
			if(track[i]==track[j] && side[i]==side[j] && layer[i]==layer[j] && bar[i]==bar[j]){
				Pos_X[i] = (Pos_X[i]*Edep[i]+Pos_X[j]*Edep[j])/(Edep[i]+Edep[j]);
				Pos_Y[i] = (Pos_Y[i]*Edep[i]+Pos_Y[j]*Edep[j])/(Edep[i]+Edep[j]);
				Pos_Z[i] = (Pos_Z[i]*Edep[i]+Pos_Z[j]*Edep[j])/(Edep[i]+Edep[j]);
				Edep[i] += Edep[j];
				Time[i] = (Time[i]<Time[j]) ? Time[i] : Time[j];
				erase(j);
				j--;
			}
		}
	}
}

void HRDSciHit::response()
{
	for(int i=0; i<nHits; i++){
		double x1 = length/2.-Pos_X[i];
		double x2 = length/2.+Pos_X[i];
		double y1 = width/4.-Pos_Y[i];
		double y2 = width/4.+Pos_Y[i];
		double r1 = sqrt(pow(y1,2)+pow(Pos_Z[i],2));
		r1 = (r1>width/2.) ? width/2. : r1;
		double r2 = sqrt(pow(y2,2)+pow(Pos_Z[i],2));
		r2 = (r2>width/2.) ? width/2. : r2;

		double fr1 = A[0]*exp(-r1/A[1])+A[2];
		double fr2 = A[0]*exp(-r2/A[1])+A[2];
		fr1=rndm->Gaus(fr1,A[3]);
		fr2=rndm->Gaus(fr2,A[3]);

		double fx1 = B[0]*exp(-x1/B[1]);
		double fx2 = B[0]*exp(-x2/B[1]);
		fx1=rndm->Gaus(fx1,B[2]);
		fx2=rndm->Gaus(fx2,B[2]);

		double nph=Edep[i]*phYield;
		nph=rndm->Gaus(nph,sqrt(nph));

		nPh1.push_back(nph*fr1*fx1*rndm->Gaus(C[0],C[1]));
		nPh2.push_back(nph*fr2*fx1*rndm->Gaus(C[0],C[1]));
		nPh3.push_back(nph*fr1*fx2*rndm->Gaus(C[0],C[1]));
		nPh4.push_back(nph*fr2*fx2*rndm->Gaus(C[0],C[1]));

		nPh.push_back(nPh1[i]+nPh2[i]+nPh3[i]+nPh4[i]);
		Pos.push_back(log((nPh1[i]+nPh2[i])/(nPh3[i]+nPh4[i]))*2*length);

	}
}

void HRDSciHit::sort()
{
	for(int i=0; i<nHits-1; i++){
		for (int j=i+1; j<nHits; j++){
			if(nPh[i]<nPh[j]){
				swap(i,j);
			}
		}
	}
}

void HRDSciHit::erase(int i)
{
	track.erase(track.begin()+i);
	side.erase(side.begin()+i);
	layer.erase(layer.begin()+i);
	bar.erase(bar.begin()+i);
	Edep.erase(Edep.begin()+i);
	Time.erase(Time.begin()+i);
	Pos_X.erase(Pos_X.begin()+i);
	Pos_Y.erase(Pos_Y.begin()+i);
	Pos_Z.erase(Pos_Z.begin()+i);
	nHits--;
}

void HRDSciHit::swap(int i, int j)
{
	iter_swap(track.begin()+i,track.begin()+j);
	iter_swap(side.begin()+i,side.begin()+j);
	iter_swap(layer.begin()+i,layer.begin()+j);
	iter_swap(bar.begin()+i,bar.begin()+j);
	iter_swap(Edep.begin()+i,Edep.begin()+j);
	iter_swap(Time.begin()+i,Time.begin()+j);
	iter_swap(Pos_X.begin()+i,Pos_X.begin()+j);
	iter_swap(Pos_Y.begin()+i,Pos_Y.begin()+j);
	iter_swap(Pos_Z.begin()+i,Pos_Z.begin()+j);
	iter_swap(nPh1.begin()+i,nPh1.begin()+j);
	iter_swap(nPh2.begin()+i,nPh2.begin()+j);
	iter_swap(nPh3.begin()+i,nPh3.begin()+j);
	iter_swap(nPh4.begin()+i,nPh4.begin()+j);
	iter_swap(nPh.begin()+i,nPh.begin()+j);
	iter_swap(Pos.begin()+i,Pos.begin()+j);
}


void HRDSciHit::clear()
{
	nHits=0;
	track.clear();
	side.clear();
	layer.clear();
	bar.clear();
	Edep.clear();
	Time.clear();
	Pos_X.clear();
	Pos_Y.clear();
	Pos_Z.clear();
	nPh1.clear();
	nPh2.clear();
	nPh3.clear();
	nPh4.clear();
	nPh.clear();
	Pos.clear();
}

ClassImp(HRDSciHit)
