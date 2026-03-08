#include "TPCHit.hh"

TPCHit::TPCHit() : nHits(0)
{
	rndm = new TRandom3();
}

TPCHit::TPCHit(const TPCHit& hit)
{
	nHits = hit.nHits;
	Edep = hit.Edep;
	nEl = hit.nEl;
	timestamp = hit.timestamp;
	padRow = hit.padRow;
	padColumn = hit.padColumn;
	
	rndm = new TRandom3();
}

const TPCHit& TPCHit::operator=(const TPCHit& hit)
{
	nHits = hit.nHits;
	Edep = hit.Edep;
	nEl = hit.nEl;
	timestamp = hit.timestamp;
	padRow = hit.padRow;
	padColumn = hit.padColumn;

	return *this;
}

bool TPCHit::operator==(const TPCHit& hit) const { return padRow == hit.padRow && padColumn == hit.padColumn; }

void TPCHit::push_back(double e, double t, double x, double y, double z)
{
	nHits++;
	Edep.push_back(e);
	Time.push_back(t);
	Pos_X.push_back(x);
	Pos_Y.push_back(y);
	Pos_Z.push_back(z);
}


void TPCHit::assign_pads(double e, double t, double x, double y, double z)
{
	//if(abs(x)<GEMsize/2.&&abs(y)<GEMsize/2.){
		int ne = int(rndm->Gaus(e/ionE,sqrt(e/ionE)));
		ne = (ne<0) ? -ne : ne;
		for(int i=0; i<ne; i++){
			double padC = (rndm->Gaus(x,cloudR)+GEMsize/2.)/padWidth;
			double padR = (rndm->Gaus(y,cloudR)+GEMsize/2.)/padHeight;
			double padT = (t + (rndm->Gaus(z,cloudR)+length/2.)/driftV)/timeBin;
			padRow.push_back(int(padR));
			padColumn.push_back(int(padC));
			timestamp.push_back(int(padT));
		}
		Edep.push_back(e);
		nEl.push_back(ne);
		nHits++;
	//}
}

void TPCHit::add_up()
{
//	for(int i=0; i<nHits; i++){
//		if(Edep[i]<0.003){
//			erase(i);
//			i--;
//		}
//	}
	for(int i=0; i<nHits-1; i++){
		for (int j=i+1; j<nHits; j++){
			if(padRow[i]==padRow[j] && padColumn[i]==padColumn[j] && timestamp[i]==timestamp[j]){
				Edep[i] += Edep[j];
				nEl[i] += nEl[j];
				erase(j);
				j--;
			}
		}
	}
}

void TPCHit::erase(int i)
{
	padRow.erase(padRow.begin()+i);
	padColumn.erase(padColumn.begin()+i);
	timestamp.erase(timestamp.begin()+i);
	Edep.erase(Edep.begin()+i);
	nEl.erase(nEl.begin()+i);
	nHits--;
}

void TPCHit::clear()
{
	nHits=0;
	padRow.clear();
	padColumn.clear();
	Edep.clear();
	nEl.clear();
	timestamp.clear();
}

ClassImp(TPCHit)
