#ifndef ScintAna_h
#define ScintAna_h

#include <algorithm>
#include <vector>
#include "TObject.h"
#include "TClass.h"
#include "TRandom3.h"

#include "ScintHit.hh"

class ScintAna : public TObject{
	public :
		ScintAna();
  		ScintAna(const ScintAna& ana);
		virtual ~ScintAna() {};
  
		const ScintAna& operator=(const ScintAna& right);

		void response(ScintHit *sc);
		void add_up();
		void sort();
		void sumE();
		void erase(int i);
		void swap(int i, int j);
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
		std::vector<int> side;
		std::vector<int> layer;
		std::vector<int> bar;
		std::vector<double> nPh1;
		std::vector<double> nPh2;
		std::vector<double> nPh3;
		std::vector<double> nPh4;
		std::vector<double> nPh;
		std::vector<double> Pos;
		std::vector<int> mult;
		std::vector<double> dE;
		std::vector<double> Esum1;
		std::vector<double> Esum2;
		std::vector<double> Pos1;

		TRandom3 *rndm;
		ClassDef(ScintAna,1);

};
#endif
