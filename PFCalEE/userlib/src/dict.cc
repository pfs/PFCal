// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME srcdIdict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
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

// The generated code does not explicitly qualifies STL entities
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "HGCSSInfo.hh"
#include "HGCSSEvent.hh"
#include "HGCSSSamplingSection.hh"
#include "HGCSSSimHit.hh"
#include "HGCSSGenParticle.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSRecoJet.hh"
#include "HGCSSCluster.hh"
#include "HGCSSMipHit.hh"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_HGCSSInfo(void *p = 0);
   static void *newArray_HGCSSInfo(Long_t size, void *p);
   static void delete_HGCSSInfo(void *p);
   static void deleteArray_HGCSSInfo(void *p);
   static void destruct_HGCSSInfo(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::HGCSSInfo*)
   {
      ::HGCSSInfo *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::HGCSSInfo >(0);
      static ::ROOT::TGenericClassInfo 
         instance("HGCSSInfo", ::HGCSSInfo::Class_Version(), "HGCSSInfo.hh", 9,
                  typeid(::HGCSSInfo), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::HGCSSInfo::Dictionary, isa_proxy, 4,
                  sizeof(::HGCSSInfo) );
      instance.SetNew(&new_HGCSSInfo);
      instance.SetNewArray(&newArray_HGCSSInfo);
      instance.SetDelete(&delete_HGCSSInfo);
      instance.SetDeleteArray(&deleteArray_HGCSSInfo);
      instance.SetDestructor(&destruct_HGCSSInfo);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::HGCSSInfo*)
   {
      return GenerateInitInstanceLocal((::HGCSSInfo*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::HGCSSInfo*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_HGCSSEvent(void *p = 0);
   static void *newArray_HGCSSEvent(Long_t size, void *p);
   static void delete_HGCSSEvent(void *p);
   static void deleteArray_HGCSSEvent(void *p);
   static void destruct_HGCSSEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::HGCSSEvent*)
   {
      ::HGCSSEvent *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::HGCSSEvent >(0);
      static ::ROOT::TGenericClassInfo 
         instance("HGCSSEvent", ::HGCSSEvent::Class_Version(), "HGCSSEvent.hh", 9,
                  typeid(::HGCSSEvent), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::HGCSSEvent::Dictionary, isa_proxy, 4,
                  sizeof(::HGCSSEvent) );
      instance.SetNew(&new_HGCSSEvent);
      instance.SetNewArray(&newArray_HGCSSEvent);
      instance.SetDelete(&delete_HGCSSEvent);
      instance.SetDeleteArray(&deleteArray_HGCSSEvent);
      instance.SetDestructor(&destruct_HGCSSEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::HGCSSEvent*)
   {
      return GenerateInitInstanceLocal((::HGCSSEvent*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::HGCSSEvent*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_HGCSSSamplingSection(void *p = 0);
   static void *newArray_HGCSSSamplingSection(Long_t size, void *p);
   static void delete_HGCSSSamplingSection(void *p);
   static void deleteArray_HGCSSSamplingSection(void *p);
   static void destruct_HGCSSSamplingSection(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::HGCSSSamplingSection*)
   {
      ::HGCSSSamplingSection *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::HGCSSSamplingSection >(0);
      static ::ROOT::TGenericClassInfo 
         instance("HGCSSSamplingSection", ::HGCSSSamplingSection::Class_Version(), "HGCSSSamplingSection.hh", 9,
                  typeid(::HGCSSSamplingSection), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::HGCSSSamplingSection::Dictionary, isa_proxy, 4,
                  sizeof(::HGCSSSamplingSection) );
      instance.SetNew(&new_HGCSSSamplingSection);
      instance.SetNewArray(&newArray_HGCSSSamplingSection);
      instance.SetDelete(&delete_HGCSSSamplingSection);
      instance.SetDeleteArray(&deleteArray_HGCSSSamplingSection);
      instance.SetDestructor(&destruct_HGCSSSamplingSection);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::HGCSSSamplingSection*)
   {
      return GenerateInitInstanceLocal((::HGCSSSamplingSection*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::HGCSSSamplingSection*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_HGCSSSimHit(void *p = 0);
   static void *newArray_HGCSSSimHit(Long_t size, void *p);
   static void delete_HGCSSSimHit(void *p);
   static void deleteArray_HGCSSSimHit(void *p);
   static void destruct_HGCSSSimHit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::HGCSSSimHit*)
   {
      ::HGCSSSimHit *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::HGCSSSimHit >(0);
      static ::ROOT::TGenericClassInfo 
         instance("HGCSSSimHit", ::HGCSSSimHit::Class_Version(), "HGCSSSimHit.hh", 27,
                  typeid(::HGCSSSimHit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::HGCSSSimHit::Dictionary, isa_proxy, 4,
                  sizeof(::HGCSSSimHit) );
      instance.SetNew(&new_HGCSSSimHit);
      instance.SetNewArray(&newArray_HGCSSSimHit);
      instance.SetDelete(&delete_HGCSSSimHit);
      instance.SetDeleteArray(&deleteArray_HGCSSSimHit);
      instance.SetDestructor(&destruct_HGCSSSimHit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::HGCSSSimHit*)
   {
      return GenerateInitInstanceLocal((::HGCSSSimHit*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::HGCSSSimHit*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_HGCSSGenParticle(void *p = 0);
   static void *newArray_HGCSSGenParticle(Long_t size, void *p);
   static void delete_HGCSSGenParticle(void *p);
   static void deleteArray_HGCSSGenParticle(void *p);
   static void destruct_HGCSSGenParticle(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::HGCSSGenParticle*)
   {
      ::HGCSSGenParticle *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::HGCSSGenParticle >(0);
      static ::ROOT::TGenericClassInfo 
         instance("HGCSSGenParticle", ::HGCSSGenParticle::Class_Version(), "HGCSSGenParticle.hh", 12,
                  typeid(::HGCSSGenParticle), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::HGCSSGenParticle::Dictionary, isa_proxy, 4,
                  sizeof(::HGCSSGenParticle) );
      instance.SetNew(&new_HGCSSGenParticle);
      instance.SetNewArray(&newArray_HGCSSGenParticle);
      instance.SetDelete(&delete_HGCSSGenParticle);
      instance.SetDeleteArray(&deleteArray_HGCSSGenParticle);
      instance.SetDestructor(&destruct_HGCSSGenParticle);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::HGCSSGenParticle*)
   {
      return GenerateInitInstanceLocal((::HGCSSGenParticle*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::HGCSSGenParticle*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_HGCSSRecoHit(void *p = 0);
   static void *newArray_HGCSSRecoHit(Long_t size, void *p);
   static void delete_HGCSSRecoHit(void *p);
   static void deleteArray_HGCSSRecoHit(void *p);
   static void destruct_HGCSSRecoHit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::HGCSSRecoHit*)
   {
      ::HGCSSRecoHit *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::HGCSSRecoHit >(0);
      static ::ROOT::TGenericClassInfo 
         instance("HGCSSRecoHit", ::HGCSSRecoHit::Class_Version(), "HGCSSRecoHit.hh", 16,
                  typeid(::HGCSSRecoHit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::HGCSSRecoHit::Dictionary, isa_proxy, 4,
                  sizeof(::HGCSSRecoHit) );
      instance.SetNew(&new_HGCSSRecoHit);
      instance.SetNewArray(&newArray_HGCSSRecoHit);
      instance.SetDelete(&delete_HGCSSRecoHit);
      instance.SetDeleteArray(&deleteArray_HGCSSRecoHit);
      instance.SetDestructor(&destruct_HGCSSRecoHit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::HGCSSRecoHit*)
   {
      return GenerateInitInstanceLocal((::HGCSSRecoHit*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::HGCSSRecoHit*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_HGCSSRecoJet(void *p = 0);
   static void *newArray_HGCSSRecoJet(Long_t size, void *p);
   static void delete_HGCSSRecoJet(void *p);
   static void deleteArray_HGCSSRecoJet(void *p);
   static void destruct_HGCSSRecoJet(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::HGCSSRecoJet*)
   {
      ::HGCSSRecoJet *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::HGCSSRecoJet >(0);
      static ::ROOT::TGenericClassInfo 
         instance("HGCSSRecoJet", ::HGCSSRecoJet::Class_Version(), "HGCSSRecoJet.hh", 9,
                  typeid(::HGCSSRecoJet), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::HGCSSRecoJet::Dictionary, isa_proxy, 4,
                  sizeof(::HGCSSRecoJet) );
      instance.SetNew(&new_HGCSSRecoJet);
      instance.SetNewArray(&newArray_HGCSSRecoJet);
      instance.SetDelete(&delete_HGCSSRecoJet);
      instance.SetDeleteArray(&deleteArray_HGCSSRecoJet);
      instance.SetDestructor(&destruct_HGCSSRecoJet);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::HGCSSRecoJet*)
   {
      return GenerateInitInstanceLocal((::HGCSSRecoJet*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::HGCSSRecoJet*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_HGCSSCluster(void *p = 0);
   static void *newArray_HGCSSCluster(Long_t size, void *p);
   static void delete_HGCSSCluster(void *p);
   static void deleteArray_HGCSSCluster(void *p);
   static void destruct_HGCSSCluster(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::HGCSSCluster*)
   {
      ::HGCSSCluster *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::HGCSSCluster >(0);
      static ::ROOT::TGenericClassInfo 
         instance("HGCSSCluster", ::HGCSSCluster::Class_Version(), "HGCSSCluster.hh", 21,
                  typeid(::HGCSSCluster), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::HGCSSCluster::Dictionary, isa_proxy, 4,
                  sizeof(::HGCSSCluster) );
      instance.SetNew(&new_HGCSSCluster);
      instance.SetNewArray(&newArray_HGCSSCluster);
      instance.SetDelete(&delete_HGCSSCluster);
      instance.SetDeleteArray(&deleteArray_HGCSSCluster);
      instance.SetDestructor(&destruct_HGCSSCluster);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::HGCSSCluster*)
   {
      return GenerateInitInstanceLocal((::HGCSSCluster*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::HGCSSCluster*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_HGCSSMipHit(void *p = 0);
   static void *newArray_HGCSSMipHit(Long_t size, void *p);
   static void delete_HGCSSMipHit(void *p);
   static void deleteArray_HGCSSMipHit(void *p);
   static void destruct_HGCSSMipHit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::HGCSSMipHit*)
   {
      ::HGCSSMipHit *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::HGCSSMipHit >(0);
      static ::ROOT::TGenericClassInfo 
         instance("HGCSSMipHit", ::HGCSSMipHit::Class_Version(), "HGCSSMipHit.hh", 11,
                  typeid(::HGCSSMipHit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::HGCSSMipHit::Dictionary, isa_proxy, 4,
                  sizeof(::HGCSSMipHit) );
      instance.SetNew(&new_HGCSSMipHit);
      instance.SetNewArray(&newArray_HGCSSMipHit);
      instance.SetDelete(&delete_HGCSSMipHit);
      instance.SetDeleteArray(&deleteArray_HGCSSMipHit);
      instance.SetDestructor(&destruct_HGCSSMipHit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::HGCSSMipHit*)
   {
      return GenerateInitInstanceLocal((::HGCSSMipHit*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::HGCSSMipHit*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr HGCSSInfo::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *HGCSSInfo::Class_Name()
{
   return "HGCSSInfo";
}

//______________________________________________________________________________
const char *HGCSSInfo::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::HGCSSInfo*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int HGCSSInfo::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::HGCSSInfo*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *HGCSSInfo::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::HGCSSInfo*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *HGCSSInfo::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::HGCSSInfo*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr HGCSSEvent::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *HGCSSEvent::Class_Name()
{
   return "HGCSSEvent";
}

//______________________________________________________________________________
const char *HGCSSEvent::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::HGCSSEvent*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int HGCSSEvent::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::HGCSSEvent*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *HGCSSEvent::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::HGCSSEvent*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *HGCSSEvent::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::HGCSSEvent*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr HGCSSSamplingSection::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *HGCSSSamplingSection::Class_Name()
{
   return "HGCSSSamplingSection";
}

//______________________________________________________________________________
const char *HGCSSSamplingSection::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::HGCSSSamplingSection*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int HGCSSSamplingSection::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::HGCSSSamplingSection*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *HGCSSSamplingSection::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::HGCSSSamplingSection*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *HGCSSSamplingSection::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::HGCSSSamplingSection*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr HGCSSSimHit::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *HGCSSSimHit::Class_Name()
{
   return "HGCSSSimHit";
}

//______________________________________________________________________________
const char *HGCSSSimHit::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::HGCSSSimHit*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int HGCSSSimHit::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::HGCSSSimHit*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *HGCSSSimHit::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::HGCSSSimHit*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *HGCSSSimHit::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::HGCSSSimHit*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr HGCSSGenParticle::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *HGCSSGenParticle::Class_Name()
{
   return "HGCSSGenParticle";
}

//______________________________________________________________________________
const char *HGCSSGenParticle::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::HGCSSGenParticle*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int HGCSSGenParticle::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::HGCSSGenParticle*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *HGCSSGenParticle::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::HGCSSGenParticle*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *HGCSSGenParticle::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::HGCSSGenParticle*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr HGCSSRecoHit::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *HGCSSRecoHit::Class_Name()
{
   return "HGCSSRecoHit";
}

//______________________________________________________________________________
const char *HGCSSRecoHit::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::HGCSSRecoHit*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int HGCSSRecoHit::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::HGCSSRecoHit*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *HGCSSRecoHit::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::HGCSSRecoHit*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *HGCSSRecoHit::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::HGCSSRecoHit*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr HGCSSRecoJet::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *HGCSSRecoJet::Class_Name()
{
   return "HGCSSRecoJet";
}

//______________________________________________________________________________
const char *HGCSSRecoJet::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::HGCSSRecoJet*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int HGCSSRecoJet::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::HGCSSRecoJet*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *HGCSSRecoJet::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::HGCSSRecoJet*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *HGCSSRecoJet::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::HGCSSRecoJet*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr HGCSSCluster::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *HGCSSCluster::Class_Name()
{
   return "HGCSSCluster";
}

//______________________________________________________________________________
const char *HGCSSCluster::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::HGCSSCluster*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int HGCSSCluster::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::HGCSSCluster*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *HGCSSCluster::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::HGCSSCluster*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *HGCSSCluster::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::HGCSSCluster*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr HGCSSMipHit::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *HGCSSMipHit::Class_Name()
{
   return "HGCSSMipHit";
}

//______________________________________________________________________________
const char *HGCSSMipHit::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::HGCSSMipHit*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int HGCSSMipHit::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::HGCSSMipHit*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *HGCSSMipHit::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::HGCSSMipHit*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *HGCSSMipHit::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::HGCSSMipHit*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void HGCSSInfo::Streamer(TBuffer &R__b)
{
   // Stream an object of class HGCSSInfo.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(HGCSSInfo::Class(),this);
   } else {
      R__b.WriteClassBuffer(HGCSSInfo::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_HGCSSInfo(void *p) {
      return  p ? new(p) ::HGCSSInfo : new ::HGCSSInfo;
   }
   static void *newArray_HGCSSInfo(Long_t nElements, void *p) {
      return p ? new(p) ::HGCSSInfo[nElements] : new ::HGCSSInfo[nElements];
   }
   // Wrapper around operator delete
   static void delete_HGCSSInfo(void *p) {
      delete ((::HGCSSInfo*)p);
   }
   static void deleteArray_HGCSSInfo(void *p) {
      delete [] ((::HGCSSInfo*)p);
   }
   static void destruct_HGCSSInfo(void *p) {
      typedef ::HGCSSInfo current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::HGCSSInfo

//______________________________________________________________________________
void HGCSSEvent::Streamer(TBuffer &R__b)
{
   // Stream an object of class HGCSSEvent.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(HGCSSEvent::Class(),this);
   } else {
      R__b.WriteClassBuffer(HGCSSEvent::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_HGCSSEvent(void *p) {
      return  p ? new(p) ::HGCSSEvent : new ::HGCSSEvent;
   }
   static void *newArray_HGCSSEvent(Long_t nElements, void *p) {
      return p ? new(p) ::HGCSSEvent[nElements] : new ::HGCSSEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_HGCSSEvent(void *p) {
      delete ((::HGCSSEvent*)p);
   }
   static void deleteArray_HGCSSEvent(void *p) {
      delete [] ((::HGCSSEvent*)p);
   }
   static void destruct_HGCSSEvent(void *p) {
      typedef ::HGCSSEvent current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::HGCSSEvent

//______________________________________________________________________________
void HGCSSSamplingSection::Streamer(TBuffer &R__b)
{
   // Stream an object of class HGCSSSamplingSection.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(HGCSSSamplingSection::Class(),this);
   } else {
      R__b.WriteClassBuffer(HGCSSSamplingSection::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_HGCSSSamplingSection(void *p) {
      return  p ? new(p) ::HGCSSSamplingSection : new ::HGCSSSamplingSection;
   }
   static void *newArray_HGCSSSamplingSection(Long_t nElements, void *p) {
      return p ? new(p) ::HGCSSSamplingSection[nElements] : new ::HGCSSSamplingSection[nElements];
   }
   // Wrapper around operator delete
   static void delete_HGCSSSamplingSection(void *p) {
      delete ((::HGCSSSamplingSection*)p);
   }
   static void deleteArray_HGCSSSamplingSection(void *p) {
      delete [] ((::HGCSSSamplingSection*)p);
   }
   static void destruct_HGCSSSamplingSection(void *p) {
      typedef ::HGCSSSamplingSection current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::HGCSSSamplingSection

//______________________________________________________________________________
void HGCSSSimHit::Streamer(TBuffer &R__b)
{
   // Stream an object of class HGCSSSimHit.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(HGCSSSimHit::Class(),this);
   } else {
      R__b.WriteClassBuffer(HGCSSSimHit::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_HGCSSSimHit(void *p) {
      return  p ? new(p) ::HGCSSSimHit : new ::HGCSSSimHit;
   }
   static void *newArray_HGCSSSimHit(Long_t nElements, void *p) {
      return p ? new(p) ::HGCSSSimHit[nElements] : new ::HGCSSSimHit[nElements];
   }
   // Wrapper around operator delete
   static void delete_HGCSSSimHit(void *p) {
      delete ((::HGCSSSimHit*)p);
   }
   static void deleteArray_HGCSSSimHit(void *p) {
      delete [] ((::HGCSSSimHit*)p);
   }
   static void destruct_HGCSSSimHit(void *p) {
      typedef ::HGCSSSimHit current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::HGCSSSimHit

//______________________________________________________________________________
void HGCSSGenParticle::Streamer(TBuffer &R__b)
{
   // Stream an object of class HGCSSGenParticle.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(HGCSSGenParticle::Class(),this);
   } else {
      R__b.WriteClassBuffer(HGCSSGenParticle::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_HGCSSGenParticle(void *p) {
      return  p ? new(p) ::HGCSSGenParticle : new ::HGCSSGenParticle;
   }
   static void *newArray_HGCSSGenParticle(Long_t nElements, void *p) {
      return p ? new(p) ::HGCSSGenParticle[nElements] : new ::HGCSSGenParticle[nElements];
   }
   // Wrapper around operator delete
   static void delete_HGCSSGenParticle(void *p) {
      delete ((::HGCSSGenParticle*)p);
   }
   static void deleteArray_HGCSSGenParticle(void *p) {
      delete [] ((::HGCSSGenParticle*)p);
   }
   static void destruct_HGCSSGenParticle(void *p) {
      typedef ::HGCSSGenParticle current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::HGCSSGenParticle

//______________________________________________________________________________
void HGCSSRecoHit::Streamer(TBuffer &R__b)
{
   // Stream an object of class HGCSSRecoHit.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(HGCSSRecoHit::Class(),this);
   } else {
      R__b.WriteClassBuffer(HGCSSRecoHit::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_HGCSSRecoHit(void *p) {
      return  p ? new(p) ::HGCSSRecoHit : new ::HGCSSRecoHit;
   }
   static void *newArray_HGCSSRecoHit(Long_t nElements, void *p) {
      return p ? new(p) ::HGCSSRecoHit[nElements] : new ::HGCSSRecoHit[nElements];
   }
   // Wrapper around operator delete
   static void delete_HGCSSRecoHit(void *p) {
      delete ((::HGCSSRecoHit*)p);
   }
   static void deleteArray_HGCSSRecoHit(void *p) {
      delete [] ((::HGCSSRecoHit*)p);
   }
   static void destruct_HGCSSRecoHit(void *p) {
      typedef ::HGCSSRecoHit current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::HGCSSRecoHit

//______________________________________________________________________________
void HGCSSRecoJet::Streamer(TBuffer &R__b)
{
   // Stream an object of class HGCSSRecoJet.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(HGCSSRecoJet::Class(),this);
   } else {
      R__b.WriteClassBuffer(HGCSSRecoJet::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_HGCSSRecoJet(void *p) {
      return  p ? new(p) ::HGCSSRecoJet : new ::HGCSSRecoJet;
   }
   static void *newArray_HGCSSRecoJet(Long_t nElements, void *p) {
      return p ? new(p) ::HGCSSRecoJet[nElements] : new ::HGCSSRecoJet[nElements];
   }
   // Wrapper around operator delete
   static void delete_HGCSSRecoJet(void *p) {
      delete ((::HGCSSRecoJet*)p);
   }
   static void deleteArray_HGCSSRecoJet(void *p) {
      delete [] ((::HGCSSRecoJet*)p);
   }
   static void destruct_HGCSSRecoJet(void *p) {
      typedef ::HGCSSRecoJet current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::HGCSSRecoJet

//______________________________________________________________________________
void HGCSSCluster::Streamer(TBuffer &R__b)
{
   // Stream an object of class HGCSSCluster.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(HGCSSCluster::Class(),this);
   } else {
      R__b.WriteClassBuffer(HGCSSCluster::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_HGCSSCluster(void *p) {
      return  p ? new(p) ::HGCSSCluster : new ::HGCSSCluster;
   }
   static void *newArray_HGCSSCluster(Long_t nElements, void *p) {
      return p ? new(p) ::HGCSSCluster[nElements] : new ::HGCSSCluster[nElements];
   }
   // Wrapper around operator delete
   static void delete_HGCSSCluster(void *p) {
      delete ((::HGCSSCluster*)p);
   }
   static void deleteArray_HGCSSCluster(void *p) {
      delete [] ((::HGCSSCluster*)p);
   }
   static void destruct_HGCSSCluster(void *p) {
      typedef ::HGCSSCluster current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::HGCSSCluster

//______________________________________________________________________________
void HGCSSMipHit::Streamer(TBuffer &R__b)
{
   // Stream an object of class HGCSSMipHit.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(HGCSSMipHit::Class(),this);
   } else {
      R__b.WriteClassBuffer(HGCSSMipHit::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_HGCSSMipHit(void *p) {
      return  p ? new(p) ::HGCSSMipHit : new ::HGCSSMipHit;
   }
   static void *newArray_HGCSSMipHit(Long_t nElements, void *p) {
      return p ? new(p) ::HGCSSMipHit[nElements] : new ::HGCSSMipHit[nElements];
   }
   // Wrapper around operator delete
   static void delete_HGCSSMipHit(void *p) {
      delete ((::HGCSSMipHit*)p);
   }
   static void deleteArray_HGCSSMipHit(void *p) {
      delete [] ((::HGCSSMipHit*)p);
   }
   static void destruct_HGCSSMipHit(void *p) {
      typedef ::HGCSSMipHit current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::HGCSSMipHit

namespace ROOT {
   static TClass *vectorlEHGCSSSimHitgR_Dictionary();
   static void vectorlEHGCSSSimHitgR_TClassManip(TClass*);
   static void *new_vectorlEHGCSSSimHitgR(void *p = 0);
   static void *newArray_vectorlEHGCSSSimHitgR(Long_t size, void *p);
   static void delete_vectorlEHGCSSSimHitgR(void *p);
   static void deleteArray_vectorlEHGCSSSimHitgR(void *p);
   static void destruct_vectorlEHGCSSSimHitgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<HGCSSSimHit>*)
   {
      vector<HGCSSSimHit> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<HGCSSSimHit>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<HGCSSSimHit>", -2, "vector", 339,
                  typeid(vector<HGCSSSimHit>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEHGCSSSimHitgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<HGCSSSimHit>) );
      instance.SetNew(&new_vectorlEHGCSSSimHitgR);
      instance.SetNewArray(&newArray_vectorlEHGCSSSimHitgR);
      instance.SetDelete(&delete_vectorlEHGCSSSimHitgR);
      instance.SetDeleteArray(&deleteArray_vectorlEHGCSSSimHitgR);
      instance.SetDestructor(&destruct_vectorlEHGCSSSimHitgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<HGCSSSimHit> >()));

      ::ROOT::AddClassAlternate("vector<HGCSSSimHit>","std::vector<HGCSSSimHit, std::allocator<HGCSSSimHit> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<HGCSSSimHit>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEHGCSSSimHitgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<HGCSSSimHit>*)0x0)->GetClass();
      vectorlEHGCSSSimHitgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEHGCSSSimHitgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEHGCSSSimHitgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<HGCSSSimHit> : new vector<HGCSSSimHit>;
   }
   static void *newArray_vectorlEHGCSSSimHitgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<HGCSSSimHit>[nElements] : new vector<HGCSSSimHit>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEHGCSSSimHitgR(void *p) {
      delete ((vector<HGCSSSimHit>*)p);
   }
   static void deleteArray_vectorlEHGCSSSimHitgR(void *p) {
      delete [] ((vector<HGCSSSimHit>*)p);
   }
   static void destruct_vectorlEHGCSSSimHitgR(void *p) {
      typedef vector<HGCSSSimHit> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<HGCSSSimHit>

namespace ROOT {
   static TClass *vectorlEHGCSSSamplingSectiongR_Dictionary();
   static void vectorlEHGCSSSamplingSectiongR_TClassManip(TClass*);
   static void *new_vectorlEHGCSSSamplingSectiongR(void *p = 0);
   static void *newArray_vectorlEHGCSSSamplingSectiongR(Long_t size, void *p);
   static void delete_vectorlEHGCSSSamplingSectiongR(void *p);
   static void deleteArray_vectorlEHGCSSSamplingSectiongR(void *p);
   static void destruct_vectorlEHGCSSSamplingSectiongR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<HGCSSSamplingSection>*)
   {
      vector<HGCSSSamplingSection> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<HGCSSSamplingSection>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<HGCSSSamplingSection>", -2, "vector", 339,
                  typeid(vector<HGCSSSamplingSection>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEHGCSSSamplingSectiongR_Dictionary, isa_proxy, 4,
                  sizeof(vector<HGCSSSamplingSection>) );
      instance.SetNew(&new_vectorlEHGCSSSamplingSectiongR);
      instance.SetNewArray(&newArray_vectorlEHGCSSSamplingSectiongR);
      instance.SetDelete(&delete_vectorlEHGCSSSamplingSectiongR);
      instance.SetDeleteArray(&deleteArray_vectorlEHGCSSSamplingSectiongR);
      instance.SetDestructor(&destruct_vectorlEHGCSSSamplingSectiongR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<HGCSSSamplingSection> >()));

      ::ROOT::AddClassAlternate("vector<HGCSSSamplingSection>","std::vector<HGCSSSamplingSection, std::allocator<HGCSSSamplingSection> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<HGCSSSamplingSection>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEHGCSSSamplingSectiongR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<HGCSSSamplingSection>*)0x0)->GetClass();
      vectorlEHGCSSSamplingSectiongR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEHGCSSSamplingSectiongR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEHGCSSSamplingSectiongR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<HGCSSSamplingSection> : new vector<HGCSSSamplingSection>;
   }
   static void *newArray_vectorlEHGCSSSamplingSectiongR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<HGCSSSamplingSection>[nElements] : new vector<HGCSSSamplingSection>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEHGCSSSamplingSectiongR(void *p) {
      delete ((vector<HGCSSSamplingSection>*)p);
   }
   static void deleteArray_vectorlEHGCSSSamplingSectiongR(void *p) {
      delete [] ((vector<HGCSSSamplingSection>*)p);
   }
   static void destruct_vectorlEHGCSSSamplingSectiongR(void *p) {
      typedef vector<HGCSSSamplingSection> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<HGCSSSamplingSection>

namespace ROOT {
   static TClass *vectorlEHGCSSRecoJetgR_Dictionary();
   static void vectorlEHGCSSRecoJetgR_TClassManip(TClass*);
   static void *new_vectorlEHGCSSRecoJetgR(void *p = 0);
   static void *newArray_vectorlEHGCSSRecoJetgR(Long_t size, void *p);
   static void delete_vectorlEHGCSSRecoJetgR(void *p);
   static void deleteArray_vectorlEHGCSSRecoJetgR(void *p);
   static void destruct_vectorlEHGCSSRecoJetgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<HGCSSRecoJet>*)
   {
      vector<HGCSSRecoJet> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<HGCSSRecoJet>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<HGCSSRecoJet>", -2, "vector", 339,
                  typeid(vector<HGCSSRecoJet>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEHGCSSRecoJetgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<HGCSSRecoJet>) );
      instance.SetNew(&new_vectorlEHGCSSRecoJetgR);
      instance.SetNewArray(&newArray_vectorlEHGCSSRecoJetgR);
      instance.SetDelete(&delete_vectorlEHGCSSRecoJetgR);
      instance.SetDeleteArray(&deleteArray_vectorlEHGCSSRecoJetgR);
      instance.SetDestructor(&destruct_vectorlEHGCSSRecoJetgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<HGCSSRecoJet> >()));

      ::ROOT::AddClassAlternate("vector<HGCSSRecoJet>","std::vector<HGCSSRecoJet, std::allocator<HGCSSRecoJet> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<HGCSSRecoJet>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEHGCSSRecoJetgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<HGCSSRecoJet>*)0x0)->GetClass();
      vectorlEHGCSSRecoJetgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEHGCSSRecoJetgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEHGCSSRecoJetgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<HGCSSRecoJet> : new vector<HGCSSRecoJet>;
   }
   static void *newArray_vectorlEHGCSSRecoJetgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<HGCSSRecoJet>[nElements] : new vector<HGCSSRecoJet>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEHGCSSRecoJetgR(void *p) {
      delete ((vector<HGCSSRecoJet>*)p);
   }
   static void deleteArray_vectorlEHGCSSRecoJetgR(void *p) {
      delete [] ((vector<HGCSSRecoJet>*)p);
   }
   static void destruct_vectorlEHGCSSRecoJetgR(void *p) {
      typedef vector<HGCSSRecoJet> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<HGCSSRecoJet>

namespace ROOT {
   static TClass *vectorlEHGCSSRecoHitgR_Dictionary();
   static void vectorlEHGCSSRecoHitgR_TClassManip(TClass*);
   static void *new_vectorlEHGCSSRecoHitgR(void *p = 0);
   static void *newArray_vectorlEHGCSSRecoHitgR(Long_t size, void *p);
   static void delete_vectorlEHGCSSRecoHitgR(void *p);
   static void deleteArray_vectorlEHGCSSRecoHitgR(void *p);
   static void destruct_vectorlEHGCSSRecoHitgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<HGCSSRecoHit>*)
   {
      vector<HGCSSRecoHit> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<HGCSSRecoHit>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<HGCSSRecoHit>", -2, "vector", 339,
                  typeid(vector<HGCSSRecoHit>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEHGCSSRecoHitgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<HGCSSRecoHit>) );
      instance.SetNew(&new_vectorlEHGCSSRecoHitgR);
      instance.SetNewArray(&newArray_vectorlEHGCSSRecoHitgR);
      instance.SetDelete(&delete_vectorlEHGCSSRecoHitgR);
      instance.SetDeleteArray(&deleteArray_vectorlEHGCSSRecoHitgR);
      instance.SetDestructor(&destruct_vectorlEHGCSSRecoHitgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<HGCSSRecoHit> >()));

      ::ROOT::AddClassAlternate("vector<HGCSSRecoHit>","std::vector<HGCSSRecoHit, std::allocator<HGCSSRecoHit> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<HGCSSRecoHit>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEHGCSSRecoHitgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<HGCSSRecoHit>*)0x0)->GetClass();
      vectorlEHGCSSRecoHitgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEHGCSSRecoHitgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEHGCSSRecoHitgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<HGCSSRecoHit> : new vector<HGCSSRecoHit>;
   }
   static void *newArray_vectorlEHGCSSRecoHitgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<HGCSSRecoHit>[nElements] : new vector<HGCSSRecoHit>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEHGCSSRecoHitgR(void *p) {
      delete ((vector<HGCSSRecoHit>*)p);
   }
   static void deleteArray_vectorlEHGCSSRecoHitgR(void *p) {
      delete [] ((vector<HGCSSRecoHit>*)p);
   }
   static void destruct_vectorlEHGCSSRecoHitgR(void *p) {
      typedef vector<HGCSSRecoHit> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<HGCSSRecoHit>

namespace ROOT {
   static TClass *vectorlEHGCSSMipHitgR_Dictionary();
   static void vectorlEHGCSSMipHitgR_TClassManip(TClass*);
   static void *new_vectorlEHGCSSMipHitgR(void *p = 0);
   static void *newArray_vectorlEHGCSSMipHitgR(Long_t size, void *p);
   static void delete_vectorlEHGCSSMipHitgR(void *p);
   static void deleteArray_vectorlEHGCSSMipHitgR(void *p);
   static void destruct_vectorlEHGCSSMipHitgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<HGCSSMipHit>*)
   {
      vector<HGCSSMipHit> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<HGCSSMipHit>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<HGCSSMipHit>", -2, "vector", 339,
                  typeid(vector<HGCSSMipHit>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEHGCSSMipHitgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<HGCSSMipHit>) );
      instance.SetNew(&new_vectorlEHGCSSMipHitgR);
      instance.SetNewArray(&newArray_vectorlEHGCSSMipHitgR);
      instance.SetDelete(&delete_vectorlEHGCSSMipHitgR);
      instance.SetDeleteArray(&deleteArray_vectorlEHGCSSMipHitgR);
      instance.SetDestructor(&destruct_vectorlEHGCSSMipHitgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<HGCSSMipHit> >()));

      ::ROOT::AddClassAlternate("vector<HGCSSMipHit>","std::vector<HGCSSMipHit, std::allocator<HGCSSMipHit> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<HGCSSMipHit>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEHGCSSMipHitgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<HGCSSMipHit>*)0x0)->GetClass();
      vectorlEHGCSSMipHitgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEHGCSSMipHitgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEHGCSSMipHitgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<HGCSSMipHit> : new vector<HGCSSMipHit>;
   }
   static void *newArray_vectorlEHGCSSMipHitgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<HGCSSMipHit>[nElements] : new vector<HGCSSMipHit>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEHGCSSMipHitgR(void *p) {
      delete ((vector<HGCSSMipHit>*)p);
   }
   static void deleteArray_vectorlEHGCSSMipHitgR(void *p) {
      delete [] ((vector<HGCSSMipHit>*)p);
   }
   static void destruct_vectorlEHGCSSMipHitgR(void *p) {
      typedef vector<HGCSSMipHit> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<HGCSSMipHit>

namespace ROOT {
   static TClass *vectorlEHGCSSGenParticlegR_Dictionary();
   static void vectorlEHGCSSGenParticlegR_TClassManip(TClass*);
   static void *new_vectorlEHGCSSGenParticlegR(void *p = 0);
   static void *newArray_vectorlEHGCSSGenParticlegR(Long_t size, void *p);
   static void delete_vectorlEHGCSSGenParticlegR(void *p);
   static void deleteArray_vectorlEHGCSSGenParticlegR(void *p);
   static void destruct_vectorlEHGCSSGenParticlegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<HGCSSGenParticle>*)
   {
      vector<HGCSSGenParticle> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<HGCSSGenParticle>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<HGCSSGenParticle>", -2, "vector", 339,
                  typeid(vector<HGCSSGenParticle>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEHGCSSGenParticlegR_Dictionary, isa_proxy, 4,
                  sizeof(vector<HGCSSGenParticle>) );
      instance.SetNew(&new_vectorlEHGCSSGenParticlegR);
      instance.SetNewArray(&newArray_vectorlEHGCSSGenParticlegR);
      instance.SetDelete(&delete_vectorlEHGCSSGenParticlegR);
      instance.SetDeleteArray(&deleteArray_vectorlEHGCSSGenParticlegR);
      instance.SetDestructor(&destruct_vectorlEHGCSSGenParticlegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<HGCSSGenParticle> >()));

      ::ROOT::AddClassAlternate("vector<HGCSSGenParticle>","std::vector<HGCSSGenParticle, std::allocator<HGCSSGenParticle> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<HGCSSGenParticle>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEHGCSSGenParticlegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<HGCSSGenParticle>*)0x0)->GetClass();
      vectorlEHGCSSGenParticlegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEHGCSSGenParticlegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEHGCSSGenParticlegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<HGCSSGenParticle> : new vector<HGCSSGenParticle>;
   }
   static void *newArray_vectorlEHGCSSGenParticlegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<HGCSSGenParticle>[nElements] : new vector<HGCSSGenParticle>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEHGCSSGenParticlegR(void *p) {
      delete ((vector<HGCSSGenParticle>*)p);
   }
   static void deleteArray_vectorlEHGCSSGenParticlegR(void *p) {
      delete [] ((vector<HGCSSGenParticle>*)p);
   }
   static void destruct_vectorlEHGCSSGenParticlegR(void *p) {
      typedef vector<HGCSSGenParticle> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<HGCSSGenParticle>

namespace ROOT {
   static TClass *vectorlEHGCSSClustergR_Dictionary();
   static void vectorlEHGCSSClustergR_TClassManip(TClass*);
   static void *new_vectorlEHGCSSClustergR(void *p = 0);
   static void *newArray_vectorlEHGCSSClustergR(Long_t size, void *p);
   static void delete_vectorlEHGCSSClustergR(void *p);
   static void deleteArray_vectorlEHGCSSClustergR(void *p);
   static void destruct_vectorlEHGCSSClustergR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<HGCSSCluster>*)
   {
      vector<HGCSSCluster> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<HGCSSCluster>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<HGCSSCluster>", -2, "vector", 339,
                  typeid(vector<HGCSSCluster>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEHGCSSClustergR_Dictionary, isa_proxy, 4,
                  sizeof(vector<HGCSSCluster>) );
      instance.SetNew(&new_vectorlEHGCSSClustergR);
      instance.SetNewArray(&newArray_vectorlEHGCSSClustergR);
      instance.SetDelete(&delete_vectorlEHGCSSClustergR);
      instance.SetDeleteArray(&deleteArray_vectorlEHGCSSClustergR);
      instance.SetDestructor(&destruct_vectorlEHGCSSClustergR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<HGCSSCluster> >()));

      ::ROOT::AddClassAlternate("vector<HGCSSCluster>","std::vector<HGCSSCluster, std::allocator<HGCSSCluster> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<HGCSSCluster>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEHGCSSClustergR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<HGCSSCluster>*)0x0)->GetClass();
      vectorlEHGCSSClustergR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEHGCSSClustergR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEHGCSSClustergR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<HGCSSCluster> : new vector<HGCSSCluster>;
   }
   static void *newArray_vectorlEHGCSSClustergR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<HGCSSCluster>[nElements] : new vector<HGCSSCluster>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEHGCSSClustergR(void *p) {
      delete ((vector<HGCSSCluster>*)p);
   }
   static void deleteArray_vectorlEHGCSSClustergR(void *p) {
      delete [] ((vector<HGCSSCluster>*)p);
   }
   static void destruct_vectorlEHGCSSClustergR(void *p) {
      typedef vector<HGCSSCluster> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<HGCSSCluster>

namespace ROOT {
   static TClass *maplEHGCSSRecoHitmUcOdoublegR_Dictionary();
   static void maplEHGCSSRecoHitmUcOdoublegR_TClassManip(TClass*);
   static void *new_maplEHGCSSRecoHitmUcOdoublegR(void *p = 0);
   static void *newArray_maplEHGCSSRecoHitmUcOdoublegR(Long_t size, void *p);
   static void delete_maplEHGCSSRecoHitmUcOdoublegR(void *p);
   static void deleteArray_maplEHGCSSRecoHitmUcOdoublegR(void *p);
   static void destruct_maplEHGCSSRecoHitmUcOdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<HGCSSRecoHit*,double>*)
   {
      map<HGCSSRecoHit*,double> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<HGCSSRecoHit*,double>));
      static ::ROOT::TGenericClassInfo 
         instance("map<HGCSSRecoHit*,double>", -2, "map", 100,
                  typeid(map<HGCSSRecoHit*,double>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEHGCSSRecoHitmUcOdoublegR_Dictionary, isa_proxy, 0,
                  sizeof(map<HGCSSRecoHit*,double>) );
      instance.SetNew(&new_maplEHGCSSRecoHitmUcOdoublegR);
      instance.SetNewArray(&newArray_maplEHGCSSRecoHitmUcOdoublegR);
      instance.SetDelete(&delete_maplEHGCSSRecoHitmUcOdoublegR);
      instance.SetDeleteArray(&deleteArray_maplEHGCSSRecoHitmUcOdoublegR);
      instance.SetDestructor(&destruct_maplEHGCSSRecoHitmUcOdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<HGCSSRecoHit*,double> >()));

      ::ROOT::AddClassAlternate("map<HGCSSRecoHit*,double>","std::map<HGCSSRecoHit*, double, std::less<HGCSSRecoHit*>, std::allocator<std::pair<HGCSSRecoHit* const, double> > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<HGCSSRecoHit*,double>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEHGCSSRecoHitmUcOdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<HGCSSRecoHit*,double>*)0x0)->GetClass();
      maplEHGCSSRecoHitmUcOdoublegR_TClassManip(theClass);
   return theClass;
   }

   static void maplEHGCSSRecoHitmUcOdoublegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEHGCSSRecoHitmUcOdoublegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<HGCSSRecoHit*,double> : new map<HGCSSRecoHit*,double>;
   }
   static void *newArray_maplEHGCSSRecoHitmUcOdoublegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<HGCSSRecoHit*,double>[nElements] : new map<HGCSSRecoHit*,double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEHGCSSRecoHitmUcOdoublegR(void *p) {
      delete ((map<HGCSSRecoHit*,double>*)p);
   }
   static void deleteArray_maplEHGCSSRecoHitmUcOdoublegR(void *p) {
      delete [] ((map<HGCSSRecoHit*,double>*)p);
   }
   static void destruct_maplEHGCSSRecoHitmUcOdoublegR(void *p) {
      typedef map<HGCSSRecoHit*,double> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<HGCSSRecoHit*,double>

namespace {
  void TriggerDictionaryInitialization_dict_Impl() {
    static const char* headers[] = {
"HGCSSInfo.hh",
"HGCSSEvent.hh",
"HGCSSSamplingSection.hh",
"HGCSSSimHit.hh",
"HGCSSGenParticle.hh",
"HGCSSRecoHit.hh",
"HGCSSRecoJet.hh",
"HGCSSCluster.hh",
"HGCSSMipHit.hh",
0
    };
    static const char* includePaths[] = {
"/cvmfs/sft.cern.ch/lcg/releases/ROOT/6.22.00-dd202/x86_64-centos7-gcc8-opt/include/",
"/afs/cern.ch/work/p/psilva/HGCal/Standalone/move2cc7/PFCal/PFCalEE/userlib/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$HGCSSSamplingSection.hh")))  HGCSSSamplingSection;
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$bits/allocator.h")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}
class __attribute__((annotate("$clingAutoload$HGCSSSimHit.hh")))  HGCSSSimHit;
class __attribute__((annotate("$clingAutoload$HGCSSGenParticle.hh")))  HGCSSGenParticle;
class __attribute__((annotate("$clingAutoload$HGCSSRecoHit.hh")))  HGCSSRecoHit;
class __attribute__((annotate("$clingAutoload$HGCSSRecoJet.hh")))  HGCSSRecoJet;
class __attribute__((annotate("$clingAutoload$HGCSSCluster.hh")))  HGCSSCluster;
class __attribute__((annotate("$clingAutoload$HGCSSMipHit.hh")))  HGCSSMipHit;
class __attribute__((annotate("$clingAutoload$HGCSSInfo.hh")))  HGCSSInfo;
class __attribute__((annotate("$clingAutoload$HGCSSEvent.hh")))  HGCSSEvent;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "dict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "HGCSSInfo.hh"
#include "HGCSSEvent.hh"
#include "HGCSSSamplingSection.hh"
#include "HGCSSSimHit.hh"
#include "HGCSSGenParticle.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSRecoJet.hh"
#include "HGCSSCluster.hh"
#include "HGCSSMipHit.hh"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"HGCSSCluster", payloadCode, "@",
"HGCSSEvent", payloadCode, "@",
"HGCSSGenParticle", payloadCode, "@",
"HGCSSInfo", payloadCode, "@",
"HGCSSMipHit", payloadCode, "@",
"HGCSSRecoHit", payloadCode, "@",
"HGCSSRecoJet", payloadCode, "@",
"HGCSSSamplingSection", payloadCode, "@",
"HGCSSSimHit", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_dict() {
  TriggerDictionaryInitialization_dict_Impl();
}
