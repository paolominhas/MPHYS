#include "ProtoTPCHit.hh"

ProtoTPCHit::ProtoTPCHit() : nHits(0)
{
	rndm = new TRandom3();
}

ProtoTPCHit::ProtoTPCHit(const ProtoTPCHit& hit)
{
	nHits = hit.nHits;
	Edep = hit.Edep;
	signal = hit.signal;
	timestamp = hit.timestamp;
	padRow = hit.padRow;
	padColumn = hit.padColumn;
	
	rndm = new TRandom3();
}

const ProtoTPCHit& ProtoTPCHit::operator=(const ProtoTPCHit& hit)
{
	nHits = hit.nHits;
	Edep = hit.Edep;
	signal = hit.signal;
	timestamp = hit.timestamp;
	padRow = hit.padRow;
	padColumn = hit.padColumn;

	return *this;
}

bool ProtoTPCHit::operator==(const ProtoTPCHit& hit) const { return (padRow == hit.padRow && padColumn == hit.padColumn && timestamp == hit.timestamp); }

void ProtoTPCHit::push_back(double e, double t, double x, double y, double z, int l)
{
	nHits++;
	Edep.push_back(e);
	Time.push_back(t);
	Pos_X.push_back(x);
	Pos_Y.push_back(y);
	Pos_Z.push_back(z);
	Layer.push_back(l);
}

void ProtoTPCHit::process()
{
	swap_xy();
	shift();
	assign_pads();
	sum_up();
}

void ProtoTPCHit::add_pad(int pr, int pc, int ts, double s)
{
	padRow.push_back(pr);
	padColumn.push_back(pc);
	timestamp.push_back(ts);
	int sig = int(rndm->Gaus(s,sqrt(s)));
	signal.push_back(sig);
	nSignals++;
}

void ProtoTPCHit::assign_pads()
{
	int row, col, ts;
	double sig;
	t0=44;
	for(int i=0; i<nHits; i++){
		row = int(Pos_Y[i]/padHeight);
		col = int(Pos_X[i]/padWidth)-1;
		ts = int(Pos_Z[i]/timeBin/driftV)+t0-5;
		sig = Edep[i]/ionE;
		for(int j=0; j<3; j++){
			for(int k=0; k<40; k++){
				add_pad(row,col+j,ts+k,sig*xgaus[j]*zlandau[k]);
			}
		}
	}
}

void ProtoTPCHit::swap_xy()
{
	swap(Pos_X,Pos_Y);
}

void ProtoTPCHit::shift()
{
	for(int i=0; i<nHits; i++){
		Pos_X[i] += GEMsize/2.;
		Pos_Y[i] += double(Layer[i]);
		Pos_Z[i] += length/2.;
	}
}

void ProtoTPCHit::sum_up()
{
	for(int i=0; i<nSignals-1; i++){
		for (int j=i+1; j<nSignals; j++){
			if(padRow[i]==padRow[j] && padColumn[i]==padColumn[j] && timestamp[i]==timestamp[j]){
				signal[i] += signal[j];
				erase(j);
				j--;
			}
		}
	}

	for(int i=0; i<nSignals; i++){
		if(signal[i]<threshold){
			erase(i);
			i--;
			continue;
		}
		if(signal[i]>1023) signal[i]=1023;
	}
}

void ProtoTPCHit::erase(int i)
{
	padRow.erase(padRow.begin()+i);
	padColumn.erase(padColumn.begin()+i);
	timestamp.erase(timestamp.begin()+i);
	signal.erase(signal.begin()+i);
	nSignals--;
}

void ProtoTPCHit::clear()
{
	nHits=0;
	Edep.clear();
	Time.clear();
	Pos_X.clear();
	Pos_Y.clear();
	Pos_Z.clear();
	Layer.clear();

	nSignals=0;
	padColumn.clear();
	padRow.clear();
	timestamp.clear();
	signal.clear();
	t0=0;

	SumEdep=0.;
}

ClassImp(ProtoTPCHit)
