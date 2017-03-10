// class declaration include file below retrieved from workspace code storage
// class declaration include file below retrieved from workspace code storage
#include "HggTwoSidedCBPdf.h"

ClassImp(HggTwoSidedCBPdf);

//_____________________________________________________________________________
HggTwoSidedCBPdf:: HggTwoSidedCBPdf() {
}

//_____________________________________________________________________________
HggTwoSidedCBPdf::HggTwoSidedCBPdf(const char *name, const char *title,
				   RooAbsReal& _m, RooAbsReal& _m0,
				   RooAbsReal& _sigma, RooAbsReal& _alphaLo,
				   RooAbsReal& _nLo, RooAbsReal& _alphaHi,
				   RooAbsReal& _nHi) :
  RooAbsPdf(name, title),
  m("m", "Dependent", this, _m),
  m0("m0", "M0", this, _m0),
  sigma("sigma", "Sigma", this, _sigma),
  alphaLo("alphaLo", "Low-side Alpha", this, _alphaLo),
  nLo("nLo", "Low-side Order", this, _nLo),
  alphaHi("alphaHi", "High-side Alpha", this, _alphaHi),
  nHi("nHi", "Hig-side Order", this, _nHi)
{
}


//_____________________________________________________________________________
HggTwoSidedCBPdf::HggTwoSidedCBPdf(const HggTwoSidedCBPdf& other, const char* name) :
  RooAbsPdf(other, name), m("m", this, other.m), m0("m0", this, other.m0),
  sigma("sigma", this, other.sigma), 
  alphaLo("alphaLo", this, other.alphaLo), nLo("nLo", this, other.nLo),
  alphaHi("alphaHi", this, other.alphaHi), nHi("nHi", this, other.nHi)
{
}


//_____________________________________________________________________________
Double_t HggTwoSidedCBPdf::evaluate() const {

  Double_t t = (m-m0)/sigma;

  if (t < -alphaLo) {
    Double_t a = exp(-0.5*alphaLo*alphaLo);
    Double_t b = nLo/alphaLo - alphaLo; 
    return a/TMath::Power(alphaLo/nLo*(b - t), nLo);
  }
  else if (t > alphaHi) {
    Double_t a = exp(-0.5*alphaHi*alphaHi);
    Double_t b = nHi/alphaHi - alphaHi; 
    return a/TMath::Power(alphaHi/nHi*(b + t), nHi);
  }
  return exp(-0.5*t*t);
}


//_____________________________________________________________________________
Int_t HggTwoSidedCBPdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
{
  if( matchArgs(allVars,analVars,m) )
    return 1;
  
  return 0;
}


//_____________________________________________________________________________
Double_t HggTwoSidedCBPdf::analyticalIntegral(Int_t code, const char* rangeName) const
{
  assert(code==1);
  double result = 0;
    
  double sig = fabs((Double_t)sigma);
  double tmin = (m.min(rangeName)-m0)/sig;
  double tmax = (m.max(rangeName)-m0)/sig;
  
  if (tmin < -alphaLo) 
    result += powerLawIntegral(tmin, TMath::Min(tmax, -alphaLo), alphaLo, nLo);
  if (tmin < alphaHi && tmax > -alphaLo)
    result += gaussianIntegral(TMath::Max(tmin, -alphaLo), TMath::Min(tmax, alphaHi));
  if (tmax > alphaHi)
    result += powerLawIntegral(-tmax, TMath::Min(-tmin, -alphaHi), alphaHi, nHi);

  return sig*result;
}

//_____________________________________________________________________________
double HggTwoSidedCBPdf::gaussianIntegral(double tmin, double tmax) const
{
  return sqrt(TMath::TwoPi())*(ROOT::Math::gaussian_cdf(tmax) - ROOT::Math::gaussian_cdf(tmin));
}

//_____________________________________________________________________________
double HggTwoSidedCBPdf::powerLawIntegral(double tmin, double tmax, double alpha, double n) const
{
  double a = exp(-0.5*alpha*alpha);
  double b = n/alpha - alpha;
  return a/(1 - n)*( (b - tmin)/(TMath::Power(alpha/n*(b - tmin), n)) - (b - tmax)/(TMath::Power(alpha/n*(b - tmax), n)) );
}
