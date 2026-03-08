#ifndef CVHit_h
#define CVHit_h

#include <algorithm>
#include <vector>
#include "TObject.h"
#include "TClass.h"

class CVHit : public TObject{
	public :
		CVHit();
  		CVHit(const CVHit& hit);
		virtual ~CVHit() {};
  
		const CVHit& operator=(const CVHit& right);
  		bool operator==(const CVHit& right) const;

		void push_back(int b, int l, int s, double e, double t, double x, double y, double z);
		void add_up();
		void sort();
		void erase(int i);
		void swap(int i, int j);
		void clear();

		int nHits;
		std::vector<int> side;
		std::vector<int> layer;
		std::vector<int> bar;
		std::vector<double> Edep;
		std::vector<double> Time;
		std::vector<double> Pos_X;
		std::vector<double> Pos_Y;
		std::vector<double> Pos_Z;

		ClassDef(CVHit,1);
};
#endif
