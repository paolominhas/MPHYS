#include "GenSciHit.hh"

GenSciHit::GenSciHit() : nHits(0)
{}

GenSciHit::GenSciHit(const GenSciHit& hit)
{
	nHits = hit.nHits;
	Edep = hit.Edep;
	Time = hit.Time;
	Pos_X = hit.Pos_X;
	Pos_Y = hit.Pos_Y;
	Pos_Z = hit.Pos_Z;
}

const GenSciHit& GenSciHit::operator=(const GenSciHit& hit)
{
	nHits = hit.nHits;
	Edep = hit.Edep;
	Time = hit.Time;
	Pos_X = hit.Pos_X;
	Pos_Y = hit.Pos_Y;
	Pos_Z = hit.Pos_Z;

	return *this;
}

bool GenSciHit::operator==(const GenSciHit& hit) const { return detID == hit.detID && Time == hit.Time; }

void GenSciHit::push_back(int i, double e, double t, double x, double y, double z)
{
	nHits++;
	detID.push_back(i);
	Edep.push_back(e);
	Time.push_back(t);
	Pos_X.push_back(x);
	Pos_Y.push_back(y);
	Pos_Z.push_back(z);
}

void GenSciHit::process()
{
	sum_up();
	sort();
}

void GenSciHit::sum_up()
{
	for(int i=0; i<nHits-1; i++){
		for (int j=i+1; j<nHits; j++){
			if(detID[i]==detID[j]){
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

void GenSciHit::erase(int i)
{
	detID.erase(detID.begin()+i);
	Edep.erase(Edep.begin()+i);
	Time.erase(Time.begin()+i);
	Pos_X.erase(Pos_X.begin()+i);
	Pos_Y.erase(Pos_Y.begin()+i);
	Pos_Z.erase(Pos_Z.begin()+i);
	nHits--;
}

void GenSciHit::sort()
{
	for(int i=0; i<nHits-1; i++){
		for (int j=i+1; j<nHits; j++){
			if(detID[i]>detID[j]){
				swap(i,j);
			}
		}
	}
}

void GenSciHit::swap(int i, int j)
{
	iter_swap(detID.begin()+i,detID.begin()+j);
	iter_swap(Edep.begin()+i,Edep.begin()+j);
	iter_swap(Time.begin()+i,Time.begin()+j);
	iter_swap(Pos_X.begin()+i,Pos_X.begin()+j);
	iter_swap(Pos_Y.begin()+i,Pos_Y.begin()+j);
	iter_swap(Pos_Z.begin()+i,Pos_Z.begin()+j);
}

void GenSciHit::clear()
{
	nHits=0;
	detID.clear();
	Edep.clear();
	Time.clear();
	Pos_X.clear();
	Pos_Y.clear();
	Pos_Z.clear();
}

ClassImp(GenSciHit)
