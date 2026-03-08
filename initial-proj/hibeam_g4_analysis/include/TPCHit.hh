#ifndef TPCHit_h
#define TPCHit_h

#include <algorithm>
#include <vector>
#include "TObject.h"
#include "TClass.h"
#include "TRandom3.h"

class TPCHit : public TObject{
	public :
		TPCHit();
  		TPCHit(const TPCHit& hit);
		virtual ~TPCHit() {};
  
		const TPCHit& operator=(const TPCHit& right);
  		bool operator==(const TPCHit& right) const;

		void push_back(double e, double t, double x, double y, double z);
		void assign_pads(double e, double t, double x, double y, double z);
		void add_up();
		void erase(int i);
		void clear();

		const double length = 30.;		//! cm
		const double GEMsize = 10.;		//!.cm
		const double padHeight = 1.;	//! cm
		const double padWidth = 0.45;	//! cm
		const double timeBin = 50.;		//! ns
		const double driftV = 0.0007;	//! cm/ns
		const double ionE = 1.5e-5;		//! MeV
		const double cloudR = 0.12;	//! cm
		
		int nHits;
		
		std::vector<double> Edep;
		std::vector<double> Time;
		std::vector<double> Pos_X;
		std::vector<double> Pos_Y;
		std::vector<double> Pos_Z;

		std::vector<int> padRow;
		std::vector<int> padColumn;
		std::vector<int> timestamp;
		//std::vector<double> Edep;
		std::vector<int> nEl;

		TRandom3 *rndm;	//!
		ClassDef(TPCHit,1);
};
#endif
