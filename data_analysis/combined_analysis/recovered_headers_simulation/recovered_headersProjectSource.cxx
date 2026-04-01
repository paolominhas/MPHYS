namespace std {}
using namespace std;
#include "recovered_headersProjectHeaders.h"

#include "recovered_headersLinkDef.h"

#include "recovered_headersProjectDict.cxx"

struct DeleteObjectFunctor {
   template <typename T>
   void operator()(const T *ptr) const {
      delete ptr;
   }
   template <typename T, typename Q>
   void operator()(const std::pair<T,Q> &) const {
      // Do nothing
   }
   template <typename T, typename Q>
   void operator()(const std::pair<T,Q*> &ptr) const {
      delete ptr.second;
   }
   template <typename T, typename Q>
   void operator()(const std::pair<T*,Q> &ptr) const {
      delete ptr.first;
   }
   template <typename T, typename Q>
   void operator()(const std::pair<T*,Q*> &ptr) const {
      delete ptr.first;
      delete ptr.second;
   }
};

#ifndef ProtoTPCHit_cxx
#define ProtoTPCHit_cxx
ProtoTPCHit::ProtoTPCHit() {
}
ProtoTPCHit &ProtoTPCHit::operator=(const ProtoTPCHit & rhs)
{
   // This is NOT a copy operator=. This is actually a move operator= (for stl container's sake).
   // Use at your own risk!
   (void)rhs; // avoid warning about unused parameter
   TObject::operator=(const_cast<ProtoTPCHit &>( rhs ));
   nHits = (const_cast<ProtoTPCHit &>( rhs ).nHits);
   SumEdep = (const_cast<ProtoTPCHit &>( rhs ).SumEdep);
   Edep = (const_cast<ProtoTPCHit &>( rhs ).Edep);
   Time = (const_cast<ProtoTPCHit &>( rhs ).Time);
   Pos_X = (const_cast<ProtoTPCHit &>( rhs ).Pos_X);
   Pos_Y = (const_cast<ProtoTPCHit &>( rhs ).Pos_Y);
   Pos_Z = (const_cast<ProtoTPCHit &>( rhs ).Pos_Z);
   Layer = (const_cast<ProtoTPCHit &>( rhs ).Layer);
   nSignals = (const_cast<ProtoTPCHit &>( rhs ).nSignals);
   padRow = (const_cast<ProtoTPCHit &>( rhs ).padRow);
   padColumn = (const_cast<ProtoTPCHit &>( rhs ).padColumn);
   timestamp = (const_cast<ProtoTPCHit &>( rhs ).timestamp);
   signal = (const_cast<ProtoTPCHit &>( rhs ).signal);
   for (Int_t i=0;i<3;i++) xgaus[i] = rhs.xgaus[i];
   for (Int_t i=0;i<40;i++) zlandau[i] = rhs.zlandau[i];
   ProtoTPCHit &modrhs = const_cast<ProtoTPCHit &>( rhs );
   modrhs.Edep.clear();
   modrhs.Time.clear();
   modrhs.Pos_X.clear();
   modrhs.Pos_Y.clear();
   modrhs.Pos_Z.clear();
   modrhs.Layer.clear();
   modrhs.padRow.clear();
   modrhs.padColumn.clear();
   modrhs.timestamp.clear();
   modrhs.signal.clear();
   return *this;
}
ProtoTPCHit::ProtoTPCHit(const ProtoTPCHit & rhs)
   : TObject(const_cast<ProtoTPCHit &>( rhs ))
   , nHits(const_cast<ProtoTPCHit &>( rhs ).nHits)
   , SumEdep(const_cast<ProtoTPCHit &>( rhs ).SumEdep)
   , Edep(const_cast<ProtoTPCHit &>( rhs ).Edep)
   , Time(const_cast<ProtoTPCHit &>( rhs ).Time)
   , Pos_X(const_cast<ProtoTPCHit &>( rhs ).Pos_X)
   , Pos_Y(const_cast<ProtoTPCHit &>( rhs ).Pos_Y)
   , Pos_Z(const_cast<ProtoTPCHit &>( rhs ).Pos_Z)
   , Layer(const_cast<ProtoTPCHit &>( rhs ).Layer)
   , nSignals(const_cast<ProtoTPCHit &>( rhs ).nSignals)
   , padRow(const_cast<ProtoTPCHit &>( rhs ).padRow)
   , padColumn(const_cast<ProtoTPCHit &>( rhs ).padColumn)
   , timestamp(const_cast<ProtoTPCHit &>( rhs ).timestamp)
   , signal(const_cast<ProtoTPCHit &>( rhs ).signal)
{
   // This is NOT a copy constructor. This is actually a move constructor (for stl container's sake).
   // Use at your own risk!
   (void)rhs; // avoid warning about unused parameter
   for (Int_t i=0;i<3;i++) xgaus[i] = rhs.xgaus[i];
   for (Int_t i=0;i<40;i++) zlandau[i] = rhs.zlandau[i];
   ProtoTPCHit &modrhs = const_cast<ProtoTPCHit &>( rhs );
   modrhs.Edep.clear();
   modrhs.Time.clear();
   modrhs.Pos_X.clear();
   modrhs.Pos_Y.clear();
   modrhs.Pos_Z.clear();
   modrhs.Layer.clear();
   modrhs.padRow.clear();
   modrhs.padColumn.clear();
   modrhs.timestamp.clear();
   modrhs.signal.clear();
}
ProtoTPCHit::~ProtoTPCHit() {
}
#endif // ProtoTPCHit_cxx

#ifndef ScintHit_cxx
#define ScintHit_cxx
ScintHit::ScintHit() {
}
ScintHit &ScintHit::operator=(const ScintHit & rhs)
{
   // This is NOT a copy operator=. This is actually a move operator= (for stl container's sake).
   // Use at your own risk!
   (void)rhs; // avoid warning about unused parameter
   TObject::operator=(const_cast<ScintHit &>( rhs ));
   nHits = (const_cast<ScintHit &>( rhs ).nHits);
   side = (const_cast<ScintHit &>( rhs ).side);
   layer = (const_cast<ScintHit &>( rhs ).layer);
   bar = (const_cast<ScintHit &>( rhs ).bar);
   Edep = (const_cast<ScintHit &>( rhs ).Edep);
   Time = (const_cast<ScintHit &>( rhs ).Time);
   Pos_X = (const_cast<ScintHit &>( rhs ).Pos_X);
   Pos_Y = (const_cast<ScintHit &>( rhs ).Pos_Y);
   Pos_Z = (const_cast<ScintHit &>( rhs ).Pos_Z);
   ScintHit &modrhs = const_cast<ScintHit &>( rhs );
   modrhs.side.clear();
   modrhs.layer.clear();
   modrhs.bar.clear();
   modrhs.Edep.clear();
   modrhs.Time.clear();
   modrhs.Pos_X.clear();
   modrhs.Pos_Y.clear();
   modrhs.Pos_Z.clear();
   return *this;
}
ScintHit::ScintHit(const ScintHit & rhs)
   : TObject(const_cast<ScintHit &>( rhs ))
   , nHits(const_cast<ScintHit &>( rhs ).nHits)
   , side(const_cast<ScintHit &>( rhs ).side)
   , layer(const_cast<ScintHit &>( rhs ).layer)
   , bar(const_cast<ScintHit &>( rhs ).bar)
   , Edep(const_cast<ScintHit &>( rhs ).Edep)
   , Time(const_cast<ScintHit &>( rhs ).Time)
   , Pos_X(const_cast<ScintHit &>( rhs ).Pos_X)
   , Pos_Y(const_cast<ScintHit &>( rhs ).Pos_Y)
   , Pos_Z(const_cast<ScintHit &>( rhs ).Pos_Z)
{
   // This is NOT a copy constructor. This is actually a move constructor (for stl container's sake).
   // Use at your own risk!
   (void)rhs; // avoid warning about unused parameter
   ScintHit &modrhs = const_cast<ScintHit &>( rhs );
   modrhs.side.clear();
   modrhs.layer.clear();
   modrhs.bar.clear();
   modrhs.Edep.clear();
   modrhs.Time.clear();
   modrhs.Pos_X.clear();
   modrhs.Pos_Y.clear();
   modrhs.Pos_Z.clear();
}
ScintHit::~ScintHit() {
}
#endif // ScintHit_cxx

#ifndef ScintAna_cxx
#define ScintAna_cxx
ScintAna::ScintAna() {
   rndm = 0;
}
ScintAna &ScintAna::operator=(const ScintAna & rhs)
{
   // This is NOT a copy operator=. This is actually a move operator= (for stl container's sake).
   // Use at your own risk!
   (void)rhs; // avoid warning about unused parameter
   TObject::operator=(const_cast<ScintAna &>( rhs ));
   nHits = (const_cast<ScintAna &>( rhs ).nHits);
   side = (const_cast<ScintAna &>( rhs ).side);
   layer = (const_cast<ScintAna &>( rhs ).layer);
   bar = (const_cast<ScintAna &>( rhs ).bar);
   nPh1 = (const_cast<ScintAna &>( rhs ).nPh1);
   nPh2 = (const_cast<ScintAna &>( rhs ).nPh2);
   nPh3 = (const_cast<ScintAna &>( rhs ).nPh3);
   nPh4 = (const_cast<ScintAna &>( rhs ).nPh4);
   nPh = (const_cast<ScintAna &>( rhs ).nPh);
   Pos = (const_cast<ScintAna &>( rhs ).Pos);
   mult = (const_cast<ScintAna &>( rhs ).mult);
   dE = (const_cast<ScintAna &>( rhs ).dE);
   Esum1 = (const_cast<ScintAna &>( rhs ).Esum1);
   Esum2 = (const_cast<ScintAna &>( rhs ).Esum2);
   Pos1 = (const_cast<ScintAna &>( rhs ).Pos1);
   rndm = (const_cast<ScintAna &>( rhs ).rndm);
   ScintAna &modrhs = const_cast<ScintAna &>( rhs );
   modrhs.side.clear();
   modrhs.layer.clear();
   modrhs.bar.clear();
   modrhs.nPh1.clear();
   modrhs.nPh2.clear();
   modrhs.nPh3.clear();
   modrhs.nPh4.clear();
   modrhs.nPh.clear();
   modrhs.Pos.clear();
   modrhs.mult.clear();
   modrhs.dE.clear();
   modrhs.Esum1.clear();
   modrhs.Esum2.clear();
   modrhs.Pos1.clear();
   modrhs.rndm = 0;
   return *this;
}
ScintAna::ScintAna(const ScintAna & rhs)
   : TObject(const_cast<ScintAna &>( rhs ))
   , nHits(const_cast<ScintAna &>( rhs ).nHits)
   , side(const_cast<ScintAna &>( rhs ).side)
   , layer(const_cast<ScintAna &>( rhs ).layer)
   , bar(const_cast<ScintAna &>( rhs ).bar)
   , nPh1(const_cast<ScintAna &>( rhs ).nPh1)
   , nPh2(const_cast<ScintAna &>( rhs ).nPh2)
   , nPh3(const_cast<ScintAna &>( rhs ).nPh3)
   , nPh4(const_cast<ScintAna &>( rhs ).nPh4)
   , nPh(const_cast<ScintAna &>( rhs ).nPh)
   , Pos(const_cast<ScintAna &>( rhs ).Pos)
   , mult(const_cast<ScintAna &>( rhs ).mult)
   , dE(const_cast<ScintAna &>( rhs ).dE)
   , Esum1(const_cast<ScintAna &>( rhs ).Esum1)
   , Esum2(const_cast<ScintAna &>( rhs ).Esum2)
   , Pos1(const_cast<ScintAna &>( rhs ).Pos1)
   , rndm(const_cast<ScintAna &>( rhs ).rndm)
{
   // This is NOT a copy constructor. This is actually a move constructor (for stl container's sake).
   // Use at your own risk!
   (void)rhs; // avoid warning about unused parameter
   ScintAna &modrhs = const_cast<ScintAna &>( rhs );
   modrhs.side.clear();
   modrhs.layer.clear();
   modrhs.bar.clear();
   modrhs.nPh1.clear();
   modrhs.nPh2.clear();
   modrhs.nPh3.clear();
   modrhs.nPh4.clear();
   modrhs.nPh.clear();
   modrhs.Pos.clear();
   modrhs.mult.clear();
   modrhs.dE.clear();
   modrhs.Esum1.clear();
   modrhs.Esum2.clear();
   modrhs.Pos1.clear();
   modrhs.rndm = 0;
}
ScintAna::~ScintAna() {
   delete rndm;   rndm = 0;
}
#endif // ScintAna_cxx

