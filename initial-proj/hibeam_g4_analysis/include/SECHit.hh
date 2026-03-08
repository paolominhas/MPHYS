#ifndef SECHit_h
#define SECHit_h

#include <algorithm>
#include <vector>
#include "TObject.h"
#include "TClass.h"

class SECHit : public TObject{
	public :
		SECHit();
  		SECHit(const SECHit& hit);
		virtual ~SECHit() {};
  
		const SECHit& operator=(const SECHit& right);
  		bool operator==(const SECHit& right) const;

		void push_back(int ring, int cr, double e, double t, double x, double y, double z);
		void add_up();
		void sort();
		void erase(int i);
		void swap(int i, int j);
		void clear();

		int nHits;
		std::vector<int> ringNo;
		std::vector<int> crystalNo;
		std::vector<double> Edep;
		std::vector<double> Time;
		std::vector<double> Pos_X;
		std::vector<double> Pos_Y;
		std::vector<double> Pos_Z;

		ClassDef(SECHit,1);
};
#endif
