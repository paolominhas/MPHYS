// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME recovered_headersProjectDict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "ROOT/RConfig.hxx"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Header files passed as explicit arguments
#include "recovered_headersProjectHeaders.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void *new_ProtoTPCHit(void *p = nullptr);
   static void *newArray_ProtoTPCHit(Long_t size, void *p);
   static void delete_ProtoTPCHit(void *p);
   static void deleteArray_ProtoTPCHit(void *p);
   static void destruct_ProtoTPCHit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ProtoTPCHit*)
   {
      ::ProtoTPCHit *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::ProtoTPCHit >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("ProtoTPCHit", ::ProtoTPCHit::Class_Version(), "ProtoTPCHit.h", 17,
                  typeid(::ProtoTPCHit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::ProtoTPCHit::Dictionary, isa_proxy, 4,
                  sizeof(::ProtoTPCHit) );
      instance.SetNew(&new_ProtoTPCHit);
      instance.SetNewArray(&newArray_ProtoTPCHit);
      instance.SetDelete(&delete_ProtoTPCHit);
      instance.SetDeleteArray(&deleteArray_ProtoTPCHit);
      instance.SetDestructor(&destruct_ProtoTPCHit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ProtoTPCHit*)
   {
      return GenerateInitInstanceLocal(static_cast<::ProtoTPCHit*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::ProtoTPCHit*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_ScintHit(void *p = nullptr);
   static void *newArray_ScintHit(Long_t size, void *p);
   static void delete_ScintHit(void *p);
   static void deleteArray_ScintHit(void *p);
   static void destruct_ScintHit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ScintHit*)
   {
      ::ScintHit *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::ScintHit >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("ScintHit", ::ScintHit::Class_Version(), "ScintHit.h", 17,
                  typeid(::ScintHit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::ScintHit::Dictionary, isa_proxy, 4,
                  sizeof(::ScintHit) );
      instance.SetNew(&new_ScintHit);
      instance.SetNewArray(&newArray_ScintHit);
      instance.SetDelete(&delete_ScintHit);
      instance.SetDeleteArray(&deleteArray_ScintHit);
      instance.SetDestructor(&destruct_ScintHit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ScintHit*)
   {
      return GenerateInitInstanceLocal(static_cast<::ScintHit*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::ScintHit*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_ScintAna(void *p = nullptr);
   static void *newArray_ScintAna(Long_t size, void *p);
   static void delete_ScintAna(void *p);
   static void deleteArray_ScintAna(void *p);
   static void destruct_ScintAna(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ScintAna*)
   {
      ::ScintAna *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::ScintAna >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("ScintAna", ::ScintAna::Class_Version(), "ScintAna.h", 18,
                  typeid(::ScintAna), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::ScintAna::Dictionary, isa_proxy, 4,
                  sizeof(::ScintAna) );
      instance.SetNew(&new_ScintAna);
      instance.SetNewArray(&newArray_ScintAna);
      instance.SetDelete(&delete_ScintAna);
      instance.SetDeleteArray(&deleteArray_ScintAna);
      instance.SetDestructor(&destruct_ScintAna);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ScintAna*)
   {
      return GenerateInitInstanceLocal(static_cast<::ScintAna*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::ScintAna*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr ProtoTPCHit::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *ProtoTPCHit::Class_Name()
{
   return "ProtoTPCHit";
}

//______________________________________________________________________________
const char *ProtoTPCHit::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ProtoTPCHit*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int ProtoTPCHit::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ProtoTPCHit*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *ProtoTPCHit::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ProtoTPCHit*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *ProtoTPCHit::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ProtoTPCHit*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr ScintHit::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *ScintHit::Class_Name()
{
   return "ScintHit";
}

//______________________________________________________________________________
const char *ScintHit::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ScintHit*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int ScintHit::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ScintHit*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *ScintHit::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ScintHit*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *ScintHit::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ScintHit*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr ScintAna::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *ScintAna::Class_Name()
{
   return "ScintAna";
}

//______________________________________________________________________________
const char *ScintAna::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ScintAna*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int ScintAna::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ScintAna*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *ScintAna::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ScintAna*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *ScintAna::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ScintAna*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void ProtoTPCHit::Streamer(TBuffer &R__b)
{
   // Stream an object of class ProtoTPCHit.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(ProtoTPCHit::Class(),this);
   } else {
      R__b.WriteClassBuffer(ProtoTPCHit::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_ProtoTPCHit(void *p) {
      return  p ? new(p) ::ProtoTPCHit : new ::ProtoTPCHit;
   }
   static void *newArray_ProtoTPCHit(Long_t nElements, void *p) {
      return p ? new(p) ::ProtoTPCHit[nElements] : new ::ProtoTPCHit[nElements];
   }
   // Wrapper around operator delete
   static void delete_ProtoTPCHit(void *p) {
      delete (static_cast<::ProtoTPCHit*>(p));
   }
   static void deleteArray_ProtoTPCHit(void *p) {
      delete [] (static_cast<::ProtoTPCHit*>(p));
   }
   static void destruct_ProtoTPCHit(void *p) {
      typedef ::ProtoTPCHit current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::ProtoTPCHit

//______________________________________________________________________________
void ScintHit::Streamer(TBuffer &R__b)
{
   // Stream an object of class ScintHit.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(ScintHit::Class(),this);
   } else {
      R__b.WriteClassBuffer(ScintHit::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_ScintHit(void *p) {
      return  p ? new(p) ::ScintHit : new ::ScintHit;
   }
   static void *newArray_ScintHit(Long_t nElements, void *p) {
      return p ? new(p) ::ScintHit[nElements] : new ::ScintHit[nElements];
   }
   // Wrapper around operator delete
   static void delete_ScintHit(void *p) {
      delete (static_cast<::ScintHit*>(p));
   }
   static void deleteArray_ScintHit(void *p) {
      delete [] (static_cast<::ScintHit*>(p));
   }
   static void destruct_ScintHit(void *p) {
      typedef ::ScintHit current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::ScintHit

//______________________________________________________________________________
void ScintAna::Streamer(TBuffer &R__b)
{
   // Stream an object of class ScintAna.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(ScintAna::Class(),this);
   } else {
      R__b.WriteClassBuffer(ScintAna::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_ScintAna(void *p) {
      return  p ? new(p) ::ScintAna : new ::ScintAna;
   }
   static void *newArray_ScintAna(Long_t nElements, void *p) {
      return p ? new(p) ::ScintAna[nElements] : new ::ScintAna[nElements];
   }
   // Wrapper around operator delete
   static void delete_ScintAna(void *p) {
      delete (static_cast<::ScintAna*>(p));
   }
   static void deleteArray_ScintAna(void *p) {
      delete [] (static_cast<::ScintAna*>(p));
   }
   static void destruct_ScintAna(void *p) {
      typedef ::ScintAna current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::ScintAna

namespace ROOT {
   static TClass *vectorlEintgR_Dictionary();
   static void vectorlEintgR_TClassManip(TClass*);
   static void *new_vectorlEintgR(void *p = nullptr);
   static void *newArray_vectorlEintgR(Long_t size, void *p);
   static void delete_vectorlEintgR(void *p);
   static void deleteArray_vectorlEintgR(void *p);
   static void destruct_vectorlEintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<int>*)
   {
      vector<int> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<int>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<int>", -2, "vector", 423,
                  typeid(vector<int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEintgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<int>) );
      instance.SetNew(&new_vectorlEintgR);
      instance.SetNewArray(&newArray_vectorlEintgR);
      instance.SetDelete(&delete_vectorlEintgR);
      instance.SetDeleteArray(&deleteArray_vectorlEintgR);
      instance.SetDestructor(&destruct_vectorlEintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<int> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<int>","std::vector<int, std::allocator<int> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<int>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<int>*>(nullptr))->GetClass();
      vectorlEintgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEintgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<int> : new vector<int>;
   }
   static void *newArray_vectorlEintgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<int>[nElements] : new vector<int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEintgR(void *p) {
      delete (static_cast<vector<int>*>(p));
   }
   static void deleteArray_vectorlEintgR(void *p) {
      delete [] (static_cast<vector<int>*>(p));
   }
   static void destruct_vectorlEintgR(void *p) {
      typedef vector<int> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<int>

namespace ROOT {
   static TClass *vectorlEdoublegR_Dictionary();
   static void vectorlEdoublegR_TClassManip(TClass*);
   static void *new_vectorlEdoublegR(void *p = nullptr);
   static void *newArray_vectorlEdoublegR(Long_t size, void *p);
   static void delete_vectorlEdoublegR(void *p);
   static void deleteArray_vectorlEdoublegR(void *p);
   static void destruct_vectorlEdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<double>*)
   {
      vector<double> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<double>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<double>", -2, "vector", 423,
                  typeid(vector<double>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEdoublegR_Dictionary, isa_proxy, 0,
                  sizeof(vector<double>) );
      instance.SetNew(&new_vectorlEdoublegR);
      instance.SetNewArray(&newArray_vectorlEdoublegR);
      instance.SetDelete(&delete_vectorlEdoublegR);
      instance.SetDeleteArray(&deleteArray_vectorlEdoublegR);
      instance.SetDestructor(&destruct_vectorlEdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<double> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<double>","std::vector<double, std::allocator<double> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<double>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<double>*>(nullptr))->GetClass();
      vectorlEdoublegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEdoublegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEdoublegR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<double> : new vector<double>;
   }
   static void *newArray_vectorlEdoublegR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<double>[nElements] : new vector<double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEdoublegR(void *p) {
      delete (static_cast<vector<double>*>(p));
   }
   static void deleteArray_vectorlEdoublegR(void *p) {
      delete [] (static_cast<vector<double>*>(p));
   }
   static void destruct_vectorlEdoublegR(void *p) {
      typedef vector<double> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<double>

namespace {
  void TriggerDictionaryInitialization_recovered_headersProjectDict_Impl() {
    static const char* headers[] = {
"recovered_headersProjectHeaders.h",
nullptr
    };
    static const char* includePaths[] = {
"/cvmfs/sft.cern.ch/lcg/releases/ROOT/6.30.02-fb5be/x86_64-el9-gcc12-opt/include",
"/cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc12-opt/include/Geant4",
"/cvmfs/sft.cern.ch/lcg/releases/jsonmcpp/3.10.5-f26c3/x86_64-el9-gcc12-opt/include",
"/cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc12-opt/src/cpp",
"/cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc12-opt/include",
"/cvmfs/sft.cern.ch/lcg/releases/Python/3.9.12-9a1bc/x86_64-el9-gcc12-opt/include/python3.9",
"/cvmfs/sft.cern.ch/lcg/releases/R/4.3.0-2a3db/x86_64-el9-gcc12-opt/lib64/R/include",
"/cvmfs/sft.cern.ch/lcg/releases/R/4.3.0-2a3db/x86_64-el9-gcc12-opt/lib64/R/library/RInside/include",
"/cvmfs/sft.cern.ch/lcg/releases/R/4.3.0-2a3db/x86_64-el9-gcc12-opt/lib64/R/library/Rcpp/include",
"/cvmfs/sft.cern.ch/lcg/releases/ROOT/6.30.02-fb5be/x86_64-el9-gcc12-opt/etc/",
"/cvmfs/sft.cern.ch/lcg/releases/ROOT/6.30.02-fb5be/x86_64-el9-gcc12-opt/etc//cling",
"/cvmfs/sft.cern.ch/lcg/releases/ROOT/6.30.02-fb5be/x86_64-el9-gcc12-opt/etc//cling/plugins/include",
"/cvmfs/sft.cern.ch/lcg/releases/ROOT/6.30.02-fb5be/x86_64-el9-gcc12-opt/include/",
"/cvmfs/sft.cern.ch/lcg/releases/Python/3.9.12-9a1bc/x86_64-el9-gcc12-opt/include",
"/cvmfs/sft.cern.ch/lcg/releases/ROOT/6.30.02-fb5be/x86_64-el9-gcc12-opt/include/",
"/home/s2289940/ESS/HIBEAM/data_analysis/combined_analysis/recovered_headers/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "recovered_headersProjectDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$ProtoTPCHit.h")))  __attribute__((annotate("$clingAutoload$recovered_headersProjectHeaders.h")))  ProtoTPCHit;
class __attribute__((annotate("$clingAutoload$ScintHit.h")))  __attribute__((annotate("$clingAutoload$recovered_headersProjectHeaders.h")))  ScintHit;
class __attribute__((annotate("$clingAutoload$ScintAna.h")))  __attribute__((annotate("$clingAutoload$recovered_headersProjectHeaders.h")))  ScintAna;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "recovered_headersProjectDict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "recovered_headersProjectHeaders.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"ProtoTPCHit", payloadCode, "@",
"ScintAna", payloadCode, "@",
"ScintHit", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("recovered_headersProjectDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_recovered_headersProjectDict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_recovered_headersProjectDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_recovered_headersProjectDict() {
  TriggerDictionaryInitialization_recovered_headersProjectDict_Impl();
}
