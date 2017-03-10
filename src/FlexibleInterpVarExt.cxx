// @(#)root/roostats:$Id: FlexibleInterpVarExt.cxx 157152 2015-05-19 05:14:58Z yanght $
// Author: Kyle Cranmer, Akira Shibata
// Author: Giovanni Petrucciani (UCSD) (log-interpolation)
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//_________________________________________________
/*
BEGIN_HTML
<p>
</p>
END_HTML
*/
//

#include "RooFit.h"

#include "Riostream.h"
#include <math.h>
#include "TMath.h"

#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooArgList.h"
#include "RooMsgService.h"
#include "TMath.h"

#include "FlexibleInterpVarExt.h"

using namespace std;

ClassImp(RooStats::HistFactory::FlexibleInterpVarExt)

using namespace RooStats;
using namespace HistFactory;

//_____________________________________________________________________________
FlexibleInterpVarExt::FlexibleInterpVarExt()
{
    // Default constructor
    _paramIter = _paramList.createIterator() ;
    _nominal = 0;
    _interpBoundary=1.;
}

//_____________________________________________________________________________
FlexibleInterpVarExt::FlexibleInterpVarExt(const FlexibleInterpVar& other, const char* name) :
    // RooAbsReal(other, name),
    // _paramList("paramList",this,other._paramList),
    // _nominal(other._nominal), _low(other._low), _high(other._high), _interpCode(other._interpCode), _interpBoundary(other._interpBoundary)
    FlexibleInterpVar(other, name)
{
    // Copy constructor
    _paramIter = _paramList.createIterator() ;

}


//_____________________________________________________________________________
FlexibleInterpVarExt::~FlexibleInterpVarExt()
{
    // Destructor
    // delete _paramIter ;
}


//_____________________________________________________________________________
void FlexibleInterpVarExt::setInterpCode(RooAbsReal& param, int code){

    int index = _paramList.index(&param);
    if(index<0){
        coutE(InputArguments) << "FlexibleInterpVarExt::setInterpCode ERROR:  " << param.GetName()
            << " is not in list" << endl ;
    } else {
        coutW(InputArguments) << "FlexibleInterpVarExt::setInterpCode :  " << param.GetName()
            << " is now " << code << endl ;
        _interpCode.at(index) = code;
    }
}

//_____________________________________________________________________________
void FlexibleInterpVarExt::setAllInterpCodes(int code){

    for(unsigned int i=0; i<_interpCode.size(); ++i){
        _interpCode.at(i) = code;
    }
}

//_____________________________________________________________________________
void FlexibleInterpVarExt::printAllInterpCodes(){

    for(unsigned int i=0; i<_interpCode.size(); ++i){
        coutI(InputArguments) <<"interp code for " << _paramList.at(i)->GetName() << " = " << _interpCode.at(i) <<endl;
    }

}

//_____________________________________________________________________________
Double_t FlexibleInterpVarExt::evaluate() const
{
    // Calculate and return value of polynomial

    Double_t total(_nominal) ;
    _paramIter->Reset() ;

    RooAbsReal* param ;
    //const RooArgSet* nset = _paramList.nset() ;
    int i=0;

    while((param=(RooAbsReal*)_paramIter->Next())) {
        //    param->Print("v");

        if(_interpCode.at(i)==0){
            // piece-wise linear
            if(param->getVal()>0)
                total +=  param->getVal()*(_high.at(i) - _nominal );
            else
                total += param->getVal()*(_nominal - _low.at(i));
        } else if(_interpCode.at(i)==1){
            // pice-wise log
            if(param->getVal()>=0)
                total *= pow(_high.at(i)/_nominal, +param->getVal());
            else
                total *= pow(_low.at(i)/_nominal,  -param->getVal());
        } else if(_interpCode.at(i)==2){
            // parabolic with linear
            double a = 0.5*(_high.at(i)+_low.at(i))-_nominal;
            double b = 0.5*(_high.at(i)-_low.at(i));
            double c = 0;
            if(param->getVal()>1 ){
                total += (2*a+b)*(param->getVal()-1)+_high.at(i)-_nominal;
            } else if(param->getVal()<-1 ) {
                total += -1*(2*a-b)*(param->getVal()+1)+_low.at(i)-_nominal;
            } else {
                total +=  a*pow(param->getVal(),2) + b*param->getVal()+c;
            }
        } else if(_interpCode.at(i)==3){
            //parabolic version of log-normal
            double a = 0.5*(_high.at(i)+_low.at(i))-_nominal;
            double b = 0.5*(_high.at(i)-_low.at(i));
            double c = 0;
            if(param->getVal()>1 ){
                total += (2*a+b)*(param->getVal()-1)+_high.at(i)-_nominal;
            } else if(param->getVal()<-1 ) {
                total += -1*(2*a-b)*(param->getVal()+1)+_low.at(i)-_nominal;
            } else {
                total +=  a*pow(param->getVal(),2) + b*param->getVal()+c;
            }

        } else if(_interpCode.at(i)==4){ // Aaron Armbruster - exponential extrapolation, polynomial interpolation
            double boundary = _interpBoundary;
            // piece-wise log + parabolic
            if(param->getVal()>=boundary)
            {
                total *= pow(_high.at(i)/_nominal, +param->getVal());
            }
            else if (param->getVal() < boundary && param->getVal() > -boundary && boundary != 0)
            {
                double x0 = boundary;
                double x  = param->getVal();

                // 	double pow_up       = pow(_high.at(i)/_nominal, x0);
                // 	double pow_down     = pow(_low.at(i)/_nominal,  x0);
                // 	double pow_up_log   = pow_up*TMath::Log(_high.at(i));
                // 	double pow_down_log =-pow_down*TMath::Log(_low.at(i));
                // 	double pow_up_log2  = pow_up_log*TMath::Log(_high.at(i));
                // 	double pow_down_log2=-pow_down*TMath::Log(_low.at(i));


                //fcns+der are eq at bd
                // 	  double a =  1./(4*pow(x0, 1))*(3*A0 - x0*S1);
                // 	double b =  1./(4*pow(x0, 2))*(4*S0 - x0*A1 - 8);
                // 	double c = -1./(4*pow(x0, 3))*(  A0 - x0*S1);
                // 	double d = -1./(4*pow(x0, 4))*(2*S0 - x0*A1 - 4);
                // 	total *= 1 + a*x + b*pow(x, 2) + c*pow(x, 3) + d*pow(x, 4);


                //fcns+der+2nd_der are eq at bd

                double pow_up       = pow(_high.at(i)/_nominal, x0);
                double pow_down     = pow(_low.at(i)/_nominal,  x0);
                double pow_up_log   = pow_up*TMath::Log(_high.at(i));
                double pow_down_log = -pow_down*TMath::Log(_low.at(i));
                double pow_up_log2  = pow_up_log*TMath::Log(_high.at(i));
                double pow_down_log2= pow_down_log*TMath::Log(_low.at(i));

                double S0 = (pow_up+pow_down)/2;
                double A0 = (pow_up-pow_down)/2;
                double S1 = (pow_up_log+pow_down_log)/2;
                double A1 = (pow_up_log-pow_down_log)/2;
                double S2 = (pow_up_log2+pow_down_log2)/2;
                double A2 = (pow_up_log2-pow_down_log2)/2;

                //fcns+der+2nd_der are eq at bd
                double a = 1./(8*pow(x0, 1))*(      15*A0 -  7*x0*S1 + x0*x0*A2);
                double b = 1./(8*pow(x0, 2))*(-24 + 24*S0 -  9*x0*A1 + x0*x0*S2);
                double c = 1./(4*pow(x0, 3))*(    -  5*A0 +  5*x0*S1 - x0*x0*A2);
                double d = 1./(4*pow(x0, 4))*( 12 - 12*S0 +  7*x0*A1 - x0*x0*S2);
                double e = 1./(8*pow(x0, 5))*(    +  3*A0 -  3*x0*S1 + x0*x0*A2);
                double f = 1./(8*pow(x0, 6))*( -8 +  8*S0 -  5*x0*A1 + x0*x0*S2);

                total *= 1 + a*x + b*pow(x, 2) + c*pow(x, 3) + d*pow(x, 4) + e*pow(x, 5) + f*pow(x, 6);
            }
            else if (param->getVal()<=-boundary)
            {
                total *= pow(_low.at(i)/_nominal,  -param->getVal());
            }
        } else {
            coutE(InputArguments) << "FlexibleInterpVarExt::evaluate ERROR:  " << param->GetName()
                << " with unknown interpolation code" << endl ;
        }
        ++i;
    }

    if(total<=0) {
        total=1E-9;
    }

    return total;
}



