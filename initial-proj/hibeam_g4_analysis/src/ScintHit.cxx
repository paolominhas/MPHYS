#include "ScintHit.hh"

ScintHit::ScintHit() : nHits(0)
{}

ScintHit::ScintHit(const ScintHit& hit)
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
}

const ScintHit& ScintHit::operator=(const ScintHit& hit)
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

bool ScintHit::operator==(const ScintHit& hit) const { return (track == hit.track && side == hit.side && layer == hit.layer && bar == hit.bar); }

void ScintHit::push_back(int tr, int b, int l, int s, double e, double t, double x, double y, double z)
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

void ScintHit::add_up()
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

void ScintHit::sort()
{
	for(int i=0; i<nHits-1; i++){
		for (int j=i+1; j<nHits; j++){
			if(Edep[i]<Edep[j]){
				swap(i,j);
			}
		}
	}
}

void ScintHit::erase(int i)
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

void ScintHit::swap(int i, int j)
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
}


void ScintHit::clear()
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
}

ClassImp(ScintHit)
