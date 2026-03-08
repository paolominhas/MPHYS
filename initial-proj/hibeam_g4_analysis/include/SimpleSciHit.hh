#ifndef SimpleSciHit_h
#define SimpleSciHit_h

#include <algorithm>
#include <vector>
#include "TObject.h"
#include "TClass.h"

class SimpleSciHit : public TObject{
	public :
		SimpleSciHit();
  		SimpleSciHit(const SimpleSciHit& hit);
		virtual ~SimpleSciHit() {};
  
		const SimpleSciHit& operator=(const SimpleSciHit& right);
  		bool operator==(const SimpleSciHit& right) const;

		void push_back(int i, double e, double t, double x, double y, double z);
		void add_up();
		void erase(int i);
		void sort();
		void swap(int i, int j);
		void clear();
		
		int nHits;
		std::vector<int> detID;
		std::vector<double> Edep;
		std::vector<double> Time;
		std::vector<double> Pos_X;
		std::vector<double> Pos_Y;
		std::vector<double> Pos_Z;

		ClassDef(SimpleSciHit,1);
};
#endif
