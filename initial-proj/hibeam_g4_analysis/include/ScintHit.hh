#ifndef ScintHit_h
#define ScintHit_h

#include <algorithm>
#include <vector>
#include "TObject.h"
#include "TClass.h"

class ScintHit : public TObject{
	public :
		ScintHit();
  		ScintHit(const ScintHit& hit);
		virtual ~ScintHit() {};
  
		const ScintHit& operator=(const ScintHit& right);
  		bool operator==(const ScintHit& right) const;

		void push_back(int tr, int b, int l, int s, double e, double t, double x, double y, double z);
		void add_up();
		void sort();
		void erase(int i);
		void swap(int i, int j);
		void clear();

		int nHits;
		std::vector<int> track;	//! not written to output file
		std::vector<int> side;
		std::vector<int> layer;
		std::vector<int> bar;
		std::vector<double> Edep;
		std::vector<double> Time;
		std::vector<double> Pos_X;
		std::vector<double> Pos_Y;
		std::vector<double> Pos_Z;

		ClassDef(ScintHit,1);

};
#endif
