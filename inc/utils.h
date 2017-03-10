#ifndef UTILS_HEADER
#define UTILS_HEADER

#include "inc/RooStatsHead.h"
#include "inc/RooFitHead.h"
#include "inc/CommonHead.h"

using namespace std;
using namespace RooFit;
using namespace RooStats;

namespace utils {

  static void collectEverything( ModelConfig *mc, RooArgSet *set ) {
    set->add(*mc->GetNuisanceParameters());
    set->add(*mc->GetGlobalObservables());
    set->add(*mc->GetParametersOfInterest());
  }

  static void setAllConstant( RooArgSet *set, bool flag ) {
    TIterator *iter = set->createIterator();
    RooRealVar *parg = NULL;
    while ( (parg=(RooRealVar*)iter->Next()) ) parg->setConstant(flag);
    SafeDelete(iter);
  }
  
  static void setAllValue( RooArgSet *set, double val ) {
    TIterator *iter = set->createIterator();
    RooRealVar *parg = NULL;
    while ( (parg=(RooRealVar*)iter->Next()) ) { 
      if (val > parg->getMax()) parg->setRange(parg->getMin(), val);
      if (val < parg->getMin()) parg->setRange(val, parg->getMax());
      parg->setVal(val);
    }
    SafeDelete(iter);
  }

  static void setAllConstant( const RooArgSet *set, bool flag ){
    TIterator *iter = set->createIterator();
    RooRealVar *parg = NULL;
    while ( (parg=(RooRealVar*)iter->Next()) ) parg->setConstant(flag);
    SafeDelete(iter);
  }

}

#endif
