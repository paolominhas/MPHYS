#include "CVHit.hh"

CVHit::CVHit() : nHits(0)
{}

CVHit::CVHit(const CVHit& hit)
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
}

const CVHit& CVHit::operator=(const CVHit& hit)
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

	return *this;
}

bool CVHit::operator==(const CVHit& hit) const { return (side == hit.side && layer == hit.layer && bar == hit.bar); }

void CVHit::push_back(int b, int l, int s, double e, double t, double x, double y, double z)
{
	nHits++;
	side.push_back(s);
	layer.push_back(l);
	bar.push_back(b);
	Edep.push_back(e);
	Time.push_back(t);
	Pos_X.push_back(x);
	Pos_Y.push_back(y);
	Pos_Z.push_back(z);
}

void CVHit::add_up()
{
	for(int i=0; i<nHits-1; i++){
		for (int j=i+1; j<nHits; j++){
			if(side[i]==side[j] && layer[i]==layer[j] && bar[i]==bar[j]){
				Edep[i] += Edep[j];
				Time[i] = (Time[i]<Time[j]) ? Time[i] : Time[j];
				erase(j);
				j--;
			}
		}
	}
}

void CVHit::sort()
{
	for(int i=0; i<nHits-1; i++){
		for (int j=i+1; j<nHits; j++){
			if(Edep[i]<Edep[j]){
				swap(i,j);
			}
		}
	}
}


void CVHit::erase(int i)
{
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

void CVHit::swap(int i, int j)
{
	iter_swap(side.begin()+i,side.begin()+j);
	iter_swap(layer.begin()+i,layer.begin()+j);
	iter_swap(bar.begin()+i,bar.begin()+j);
	iter_swap(Edep.begin()+i,Edep.begin()+j);
	iter_swap(Time.begin()+i,Time.begin()+j);
	iter_swap(Pos_X.begin()+i,Pos_X.begin()+j);
	iter_swap(Pos_Y.begin()+i,Pos_Y.begin()+j);
	iter_swap(Pos_Z.begin()+i,Pos_Z.begin()+j);
}


void CVHit::clear()
{
	nHits=0;
	side.clear();
	layer.clear();
	bar.clear();
	Edep.clear();
	Time.clear();
	Pos_X.clear();
	Pos_Y.clear();
	Pos_Z.clear();
}

ClassImp(CVHit)
