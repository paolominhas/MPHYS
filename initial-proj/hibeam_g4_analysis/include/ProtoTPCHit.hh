#ifndef PROTOTPCHit_h
#define PROTOTPCHit_h

#include <algorithm>
#include <vector>
#include "TObject.h"
#include "TClass.h"
#include "TRandom3.h"

class ProtoTPCHit : public TObject{
	private :
		void add_pad(int pr, int pc, int ts, double s);
		void erase(int i);

		const double length = 23.;		//! cm
		const double GEMsize = 10.;		//!.cm
		const double padHeight = 1.;	//! cm
		const double padWidth = 0.45;	//! cm
		const double timeBin = 50.;		//! ns
		const double driftV = 0.002;	//! cm/ns
		const double ionE = 1.5e-5;		//! MeV
		const double cloudR = 0.12;		//! cm
		const int threshold = 250;		//! electrons
		const double xgaus[3] = { 0.16, 0.68, 0.16 };
		const double zlandau[40] = { 0.03159, 0.1735, 0.4436, 0.7146, 0.8716, 0.9004, 0.8435, 0.7471, 
			0.6417, 0.5429, 0.4567, 0.3842, 0.3242, 0.2751, 0.2347, 0.2017,
			0.1744, 0.1518, 0.1330, 0.1172, 0.1038, 0.0925, 0.08282, 0.07449,
			0.06729, 0.06103, 0.05557, 0.05077, 0.04655, 0.04281, 0.03948, 0.03651,
			0.03386, 0.03147, 0.02932, 0.02737, 0.02561, 0.02401, 0.02254, 0.02121 };

	public :
		ProtoTPCHit();
		ProtoTPCHit(const ProtoTPCHit& hit);
		virtual ~ProtoTPCHit() {};
  
		const ProtoTPCHit& operator=(const ProtoTPCHit& right);
  		bool operator==(const ProtoTPCHit& right) const;

		void push_back(double e, double t, double x, double y, double z, int l);
		void swap_xy();
		void shift();
		void assign_pads();
		void sum_up();
		void clean_up();
		void clear();

		int nHits;
		double SumEdep;
		std::vector<double> Edep;
		std::vector<double> Time;
		std::vector<double> Pos_X;
		std::vector<double> Pos_Y;
		std::vector<double> Pos_Z;
		std::vector<int> Layer;

		int nSignals;
		std::vector<int> padRow;
		std::vector<int> padColumn;
		std::vector<int> timestamp;
		std::vector<int> signal;

		TRandom3 *rndm;	//!
		ClassDef(ProtoTPCHit,1);
};
#endif
