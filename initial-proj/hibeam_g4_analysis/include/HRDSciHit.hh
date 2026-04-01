#ifndef HRDSciHit_h
#define HRDSciHit_h

#include <algorithm>
#include <vector>
#include "TObject.h"
#include "TClass.h"
#include "TRandom3.h"

class HRDSciHit : public TObject{
	private :
		void sum_up();
		void response();
		void sort();
		void erase(int i);
		void swap(int i, int j);

	public :
		HRDSciHit();
  		HRDSciHit(const HRDSciHit& hit);
		virtual ~HRDSciHit() {};
  
		const HRDSciHit& operator=(const HRDSciHit& right);
  		bool operator==(const HRDSciHit& right) const;

		void push_back(int tr, int b, int l, int s, double e, double t, double x, double y, double z);
		void process();
		void clear();

		const int nSides=2;	//! not written to output file
		const int nLayers=10;	//! not written to output file
		const double phYield=17400*0.64;	//! not written to output file
		const double width=5.;	//! not written to output file
		const double length=50.;	//! not written to output file
		const double thickness=2.;	//! not written to output file
		const double threshold=0.1;	//! not written to output file
		const double A[4]={0.328,18.34,0.0042,0.015};	//! not written to output file
		const double B[3]={0.2616,180.395,0.022};	//! not written to output file
		const double C[2]={0.3,0.003};	//! not written to output file
	
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
	
		std::vector<double> nPh1;
		std::vector<double> nPh2;
		std::vector<double> nPh3;
		std::vector<double> nPh4;
		std::vector<double> nPh;
		std::vector<double> Pos;

		TRandom3 *rndm;
		ClassDef(HRDSciHit,2);

};
#endif
