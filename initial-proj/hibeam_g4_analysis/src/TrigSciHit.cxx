#include "TrigSciHit.hh"

TrigSciHit::TrigSciHit() : nHits(0)
{}

TrigSciHit::TrigSciHit(const TrigSciHit& hit)
{
	nHits = hit.nHits;
	detID = hit.detID;
	Edep = hit.Edep;
	Time = hit.Time;
}

const TrigSciHit& TrigSciHit::operator=(const TrigSciHit& hit)
{
	nHits = hit.nHits;
	detID = hit.detID;
	Edep = hit.Edep;
	Time = hit.Time;

	return *this;
}

bool TrigSciHit::operator==(const TrigSciHit& hit) const { return detID == hit.detID && Time == hit.Time; }

void TrigSciHit::push_back(int i, double e, double t)
{
	nHits++;
	detID.push_back(i);
	Edep.push_back(e);
	Time.push_back(t);
}

void TrigSciHit::process()
{
	sum_up();
	sort();
}

void TrigSciHit::sum_up()
{
	for(int i=0; i<nHits-1; i++){
		for (int j=i+1; j<nHits; j++){
			if(detID[i]==detID[j] && abs(Time[i]-Time[j])<100.){
				Edep[i] += Edep[j];
				Time[i] = (Time[i]<Time[j]) ? Time[i] : Time[j];
				erase(j);
				j--;
			}
		}
	}
}

void TrigSciHit::erase(int i)
{
	detID.erase(detID.begin()+i);
	Edep.erase(Edep.begin()+i);
	Time.erase(Time.begin()+i);
	nHits--;
}

void TrigSciHit::sort()
{
	for(int i=0; i<nHits-1; i++){
		for (int j=i+1; j<nHits; j++){
			if(detID[i]>detID[j]){
				swap(i,j);
			}
		}
	}
}

void TrigSciHit::swap(int i, int j)
{
	iter_swap(detID.begin()+i,detID.begin()+j);
	iter_swap(Edep.begin()+i,Edep.begin()+j);
	iter_swap(Time.begin()+i,Time.begin()+j);
}

void TrigSciHit::clear()
{
	nHits=0;
	detID.clear();
	Edep.clear();
	Time.clear();
}

ClassImp(TrigSciHit)
