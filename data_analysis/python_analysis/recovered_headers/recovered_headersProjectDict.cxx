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
   static TClass *tpcData_Dictionary();
   static void tpcData_TClassManip(TClass*);
   static void *new_tpcData(void *p = nullptr);
   static void *newArray_tpcData(Long_t size, void *p);
   static void delete_tpcData(void *p);
   static void deleteArray_tpcData(void *p);
   static void destruct_tpcData(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::tpcData*)
   {
      ::tpcData *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::tpcData));
      static ::ROOT::TGenericClassInfo 
         instance("tpcData", "tpcData.h", 14,
                  typeid(::tpcData), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &tpcData_Dictionary, isa_proxy, 4,
                  sizeof(::tpcData) );
      instance.SetNew(&new_tpcData);
      instance.SetNewArray(&newArray_tpcData);
      instance.SetDelete(&delete_tpcData);
      instance.SetDeleteArray(&deleteArray_tpcData);
      instance.SetDestructor(&destruct_tpcData);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::tpcData*)
   {
      return GenerateInitInstanceLocal(static_cast<::tpcData*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::tpcData*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *tpcData_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::tpcData*>(nullptr))->GetClass();
      tpcData_TClassManip(theClass);
   return theClass;
   }

   static void tpcData_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *CentroidData_Dictionary();
   static void CentroidData_TClassManip(TClass*);
   static void *new_CentroidData(void *p = nullptr);
   static void *newArray_CentroidData(Long_t size, void *p);
   static void delete_CentroidData(void *p);
   static void deleteArray_CentroidData(void *p);
   static void destruct_CentroidData(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::CentroidData*)
   {
      ::CentroidData *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::CentroidData));
      static ::ROOT::TGenericClassInfo 
         instance("CentroidData", "CentroidData.h", 14,
                  typeid(::CentroidData), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &CentroidData_Dictionary, isa_proxy, 4,
                  sizeof(::CentroidData) );
      instance.SetNew(&new_CentroidData);
      instance.SetNewArray(&newArray_CentroidData);
      instance.SetDelete(&delete_CentroidData);
      instance.SetDeleteArray(&deleteArray_CentroidData);
      instance.SetDestructor(&destruct_CentroidData);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::CentroidData*)
   {
      return GenerateInitInstanceLocal(static_cast<::CentroidData*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::CentroidData*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *CentroidData_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::CentroidData*>(nullptr))->GetClass();
      CentroidData_TClassManip(theClass);
   return theClass;
   }

   static void CentroidData_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *TrackData_Dictionary();
   static void TrackData_TClassManip(TClass*);
   static void *new_TrackData(void *p = nullptr);
   static void *newArray_TrackData(Long_t size, void *p);
   static void delete_TrackData(void *p);
   static void deleteArray_TrackData(void *p);
   static void destruct_TrackData(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TrackData*)
   {
      ::TrackData *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::TrackData));
      static ::ROOT::TGenericClassInfo 
         instance("TrackData", "TrackData.h", 16,
                  typeid(::TrackData), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &TrackData_Dictionary, isa_proxy, 4,
                  sizeof(::TrackData) );
      instance.SetNew(&new_TrackData);
      instance.SetNewArray(&newArray_TrackData);
      instance.SetDelete(&delete_TrackData);
      instance.SetDeleteArray(&deleteArray_TrackData);
      instance.SetDestructor(&destruct_TrackData);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TrackData*)
   {
      return GenerateInitInstanceLocal(static_cast<::TrackData*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::TrackData*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *TrackData_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::TrackData*>(nullptr))->GetClass();
      TrackData_TClassManip(theClass);
   return theClass;
   }

   static void TrackData_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_tpcData(void *p) {
      return  p ? new(p) ::tpcData : new ::tpcData;
   }
   static void *newArray_tpcData(Long_t nElements, void *p) {
      return p ? new(p) ::tpcData[nElements] : new ::tpcData[nElements];
   }
   // Wrapper around operator delete
   static void delete_tpcData(void *p) {
      delete (static_cast<::tpcData*>(p));
   }
   static void deleteArray_tpcData(void *p) {
      delete [] (static_cast<::tpcData*>(p));
   }
   static void destruct_tpcData(void *p) {
      typedef ::tpcData current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::tpcData

namespace ROOT {
   // Wrappers around operator new
   static void *new_CentroidData(void *p) {
      return  p ? new(p) ::CentroidData : new ::CentroidData;
   }
   static void *newArray_CentroidData(Long_t nElements, void *p) {
      return p ? new(p) ::CentroidData[nElements] : new ::CentroidData[nElements];
   }
   // Wrapper around operator delete
   static void delete_CentroidData(void *p) {
      delete (static_cast<::CentroidData*>(p));
   }
   static void deleteArray_CentroidData(void *p) {
      delete [] (static_cast<::CentroidData*>(p));
   }
   static void destruct_CentroidData(void *p) {
      typedef ::CentroidData current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::CentroidData

namespace ROOT {
   // Wrappers around operator new
   static void *new_TrackData(void *p) {
      return  p ? new(p) ::TrackData : new ::TrackData;
   }
   static void *newArray_TrackData(Long_t nElements, void *p) {
      return p ? new(p) ::TrackData[nElements] : new ::TrackData[nElements];
   }
   // Wrapper around operator delete
   static void delete_TrackData(void *p) {
      delete (static_cast<::TrackData*>(p));
   }
   static void deleteArray_TrackData(void *p) {
      delete [] (static_cast<::TrackData*>(p));
   }
   static void destruct_TrackData(void *p) {
      typedef ::TrackData current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::TrackData

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
   static TClass *vectorlEfloatgR_Dictionary();
   static void vectorlEfloatgR_TClassManip(TClass*);
   static void *new_vectorlEfloatgR(void *p = nullptr);
   static void *newArray_vectorlEfloatgR(Long_t size, void *p);
   static void delete_vectorlEfloatgR(void *p);
   static void deleteArray_vectorlEfloatgR(void *p);
   static void destruct_vectorlEfloatgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<float>*)
   {
      vector<float> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<float>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<float>", -2, "vector", 423,
                  typeid(vector<float>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEfloatgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<float>) );
      instance.SetNew(&new_vectorlEfloatgR);
      instance.SetNewArray(&newArray_vectorlEfloatgR);
      instance.SetDelete(&delete_vectorlEfloatgR);
      instance.SetDeleteArray(&deleteArray_vectorlEfloatgR);
      instance.SetDestructor(&destruct_vectorlEfloatgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<float> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<float>","std::vector<float, std::allocator<float> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<float>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEfloatgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<float>*>(nullptr))->GetClass();
      vectorlEfloatgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEfloatgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEfloatgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<float> : new vector<float>;
   }
   static void *newArray_vectorlEfloatgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<float>[nElements] : new vector<float>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEfloatgR(void *p) {
      delete (static_cast<vector<float>*>(p));
   }
   static void deleteArray_vectorlEfloatgR(void *p) {
      delete [] (static_cast<vector<float>*>(p));
   }
   static void destruct_vectorlEfloatgR(void *p) {
      typedef vector<float> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<float>

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
"/home/s2289940/ESS/HIBEAM/data_analysis/python_analysis/recovered_headers/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "recovered_headersProjectDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$tpcData.h")))  __attribute__((annotate("$clingAutoload$recovered_headersProjectHeaders.h")))  tpcData;
class __attribute__((annotate("$clingAutoload$CentroidData.h")))  __attribute__((annotate("$clingAutoload$recovered_headersProjectHeaders.h")))  CentroidData;
class __attribute__((annotate("$clingAutoload$TrackData.h")))  __attribute__((annotate("$clingAutoload$recovered_headersProjectHeaders.h")))  TrackData;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "recovered_headersProjectDict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "recovered_headersProjectHeaders.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"CentroidData", payloadCode, "@",
"TrackData", payloadCode, "@",
"tpcData", payloadCode, "@",
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
