#ifndef ROOT_HggMG5aMCNLOLineShapePdf
#define ROOT_HggMG5aMCNLOLineShapePdf

#include <math.h>
#include "Math/ProbFuncMathCore.h"
#include "Riostream.h"
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooFit.h"
#include "RooMath.h"
#include "RooRealProxy.h"
#include "RooRealVar.h"
#include "TMath.h"
 
// ==================================================================
// Implemented by Hongtao Yang <Hongtao.Yang@cern.ch> on May. 25, 2016
// With theory inputs provided by Martin White <martin.white@coepp.org.au> et. al.
// More details in https://indico.cern.ch/event/522134/contributions/2176494/attachments/1277773/1896607/20160524-YeeY-SignalModel.pdf
// GeV should be used as unit
// Can be used for both graivton (alpha=1.5) and scalar (alpha=0) signal
// ==================================================================
class RooRealVar;

class HggMG5aMCNLOLineShapePdf : public RooAbsPdf {
  
 public:
  
  HggMG5aMCNLOLineShapePdf();
  HggMG5aMCNLOLineShapePdf(const char *name, const char *title, RooAbsReal& x,
			   RooAbsReal& mean, RooAbsReal& width, RooAbsReal& alpha, int cme=13);
  
  HggMG5aMCNLOLineShapePdf(const HggMG5aMCNLOLineShapePdf& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { return new HggMG5aMCNLOLineShapePdf(*this,newname); }
  inline virtual ~HggMG5aMCNLOLineShapePdf() { }

 protected:

  RooRealProxy _x;		// Invariant mass
  RooRealProxy _mean;		// Mean of Briet-Wigner
  RooRealProxy _width;		// Width of Briet-Wigner
  RooRealProxy _alpha;		// Fraction of qqbar production
  int _cme;			// Center-of-mass energy
  
  Double_t ME(Double_t mgg) const;		// Matrix element part
  Double_t PL(Double_t mgg) const;		// Parton luminosity part
  Double_t BW(Double_t mgg) const;		// Breit-Wigner part
  Double_t evaluate() const;

 private:
  
  ClassDef(HggMG5aMCNLOLineShapePdf,1); // Crystal Ball lineshape PDF
};
  
#endif
