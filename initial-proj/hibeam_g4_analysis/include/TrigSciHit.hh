#ifndef TrigSciHit_h
#define TrigSciHit_h

#include <algorithm>
#include <vector>
#include "TObject.h"
#include "TClass.h"

class TrigSciHit : public TObject{
	private :
		void sum_up();
		void erase(int i);
		void sort();
		void swap(int i, int j);

	public :
		TrigSciHit();
  		TrigSciHit(const TrigSciHit& hit);
		virtual ~TrigSciHit() {};
  
		const TrigSciHit& operator=(const TrigSciHit& right);
  		bool operator==(const TrigSciHit& right) const;

		void push_back(int i, double e, double t);
		void process();
		void clear();
		
		int nHits;
		std::vector<int> detID;
		std::vector<double> Edep;
		std::vector<double> Time;

		ClassDef(TrigSciHit,1);
};
#endif
