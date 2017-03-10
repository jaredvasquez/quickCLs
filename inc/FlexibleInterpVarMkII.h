// @(#)root/roostats:$Id:  cranmer $
// Author: Kyle Cranmer, Akira Shibata
// Modified by Hongtao Yang for xmlAnaWSBuilder: we need to have both upper and lower uncertainties as functions, not just numbers
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOSTATS_FLEXIBLEINTERPVARMKII
#define ROOSTATS_FLEXIBLEINTERPVARMKII

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include <vector>
#include <ostream>

using namespace std;

class RooRealVar;
class RooArgList ;

class FlexibleInterpVarMkII : public RooAbsReal {
public:

  FlexibleInterpVarMkII() ;
  FlexibleInterpVarMkII(const char *name, const char *title,
			const RooArgList& _paramList, 
			double nominal, const RooArgList& lowList, const RooArgList& highList, std::vector<int> code);

  FlexibleInterpVarMkII(const FlexibleInterpVarMkII&, const char*);

  void setInterpCode(RooAbsReal& param, int code);
  void setAllInterpCodes(int code);
  void setGlobalBoundary(double boundary) {_interpBoundary = boundary;}

  virtual TObject* clone(const char* newname) const { return new FlexibleInterpVarMkII(*this, newname); }
  virtual ~FlexibleInterpVarMkII() ;

  virtual void printMultiline(ostream& os, Int_t contents, Bool_t verbose = kFALSE, TString indent = "") const;
  virtual void printFlexibleInterpVarMkIIs(ostream& os) const;

private:

  double PolyInterpValue(int i, double x) const;

protected:

  double _nominal;
  
  RooListProxy _paramList ;  
  RooListProxy _lowList;
  RooListProxy _highList;
  std::vector<int> _interpCode;
  Double_t _interpBoundary;

  TIterator* _paramIter ;  //! do not persist
  TIterator* _lowIter ;  //! do not persist
  TIterator* _highIter ;  //! do not persist

  mutable Bool_t         _logInit ;            //! flag used for chaching polynomial coefficients
  mutable std::vector< double>  _polCoeff;     //! cached polynomial coefficients
  
  Double_t evaluate() const;
  
  ClassDef(FlexibleInterpVarMkII,2) // flexible interpolation
};

#endif
