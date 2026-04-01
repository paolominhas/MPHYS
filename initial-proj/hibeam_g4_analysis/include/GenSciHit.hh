#ifndef GenSciHit_h
#define GenSciHit_h

#include <algorithm>
#include <vector>
#include "TObject.h"
#include "TClass.h"

class GenSciHit : public TObject{
	private :
		void sum_up();
		void erase(int i);
		void sort();
		void swap(int i, int j);

	public :
		GenSciHit();
  		GenSciHit(const GenSciHit& hit);
		virtual ~GenSciHit() {};
  
		const GenSciHit& operator=(const GenSciHit& right);
  		bool operator==(const GenSciHit& right) const;

		void push_back(int i, double e, double t, double x, double y, double z);
		void process();
		void clear();
		
		int nHits;
		std::vector<int> detID;
		std::vector<double> Edep;
		std::vector<double> Time;
		std::vector<double> Pos_X;
		std::vector<double> Pos_Y;
		std::vector<double> Pos_Z;

		ClassDef(GenSciHit,1);
};
#endif
