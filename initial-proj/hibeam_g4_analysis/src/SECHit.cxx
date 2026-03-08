#include "SECHit.hh"

SECHit::SECHit() : nHits(0)
{}

SECHit::SECHit(const SECHit& hit)
{
	nHits = hit.nHits;
	Pos_X = hit.Pos_X;
	Pos_Y = hit.Pos_Y;
	Pos_Z = hit.Pos_Z;
	Edep = hit.Edep;
	Time = hit.Time;
	ringNo = hit.ringNo;
	crystalNo = hit.crystalNo;
}

const SECHit& SECHit::operator=(const SECHit& hit)
{
	nHits = hit.nHits;
	Pos_X = hit.Pos_X;
	Pos_Y = hit.Pos_Y;
	Pos_Z = hit.Pos_Z;
	Edep = hit.Edep;
	Time = hit.Time;
	ringNo = hit.ringNo;
	crystalNo = hit.crystalNo;

	return *this;
}

bool SECHit::operator==(const SECHit& hit) const { return ringNo == hit.ringNo && crystalNo == hit.crystalNo; }

void SECHit::push_back(int ring, int cr, double e, double t, double x, double y, double z)
{
	nHits++;
	ringNo.push_back(ring);
	crystalNo.push_back(cr);
	Edep.push_back(e);
	Time.push_back(t);
	Pos_X.push_back(x);
	Pos_Y.push_back(y);
	Pos_Z.push_back(z);
}

void SECHit::add_up()
{
	for(int i=0; i<nHits-1; i++){
		for (int j=i+1; j<nHits; j++){
			if(ringNo[i]==ringNo[j] && crystalNo[i]==crystalNo[j]){
				Edep[i] += Edep[j];
				Time[i] = (Time[i]<Time[j]) ? Time[i] : Time[j];
				erase(j);
				j--;
			}
		}
	}
}

void SECHit::sort()
{
	for(int i=0; i<nHits-1; i++){
		for (int j=i+1; j<nHits; j++){
			if(Edep[i]<Edep[j]){
				swap(i,j);
			}
		}
	}
}


void SECHit::erase(int i)
{
	ringNo.erase(ringNo.begin()+i);
	crystalNo.erase(crystalNo.begin()+i);
	Edep.erase(Edep.begin()+i);
	Time.erase(Time.begin()+i);
	Pos_X.erase(Pos_X.begin()+i);
	Pos_Y.erase(Pos_Y.begin()+i);
	Pos_Z.erase(Pos_Z.begin()+i);
	nHits--;
}

void SECHit::swap(int i, int j)
{
	iter_swap(ringNo.begin()+i,ringNo.begin()+j);
	iter_swap(crystalNo.begin()+i,crystalNo.begin()+j);
	iter_swap(Edep.begin()+i,Edep.begin()+j);
	iter_swap(Time.begin()+i,Time.begin()+j);
	iter_swap(Pos_X.begin()+i,Pos_X.begin()+j);
	iter_swap(Pos_Y.begin()+i,Pos_Y.begin()+j);
	iter_swap(Pos_Z.begin()+i,Pos_Z.begin()+j);
}


void SECHit::clear()
{
	nHits=0;
	ringNo.clear();
	crystalNo.clear();
	Edep.clear();
	Time.clear();
	Pos_X.clear();
	Pos_Y.clear();
	Pos_Z.clear();
}

ClassImp(SECHit)
