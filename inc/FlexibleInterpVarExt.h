// @(#)root/roostats:$id:  cranmer $
// author: kyle cranmer, akira shibata
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOSTATS_FLEXIBLEINTERPVAREXT
#define ROOSTATS_FLEXIBLEINTERPVAREXT

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "RooStats/HistFactory/FlexibleInterpVar.h"

#include <iostream>
#include <iomanip>

using namespace RooFit;
using namespace RooStats;
using namespace HistFactory;

class RooRealVar;
class RooArgList ;

namespace RooStats{
    namespace HistFactory{

        class FlexibleInterpVarExt : public FlexibleInterpVar {
            public:

                FlexibleInterpVarExt() ;
                FlexibleInterpVarExt(const FlexibleInterpVar&, const char*);
                virtual ~FlexibleInterpVarExt();

                void setInterpCode(RooAbsReal& param, int code);
                void setAllInterpCodes(int code);
                void setGlobalBoundary(double boundary) {_interpBoundary = boundary;}

                void printAllInterpCodes();

                virtual TObject* clone(const char* newname) const { return new FlexibleInterpVarExt(*this, newname); }

                void printUncerts()
                {
                  int num  =  _low.size();
                  int varSize  =  _paramList.getSize();
                  std::cout.setf(std::ios::left);
                  for ( int  i= 0; i < num; i++ ) {
                    if ( i >= varSize ) {
                      std::cout << "WARNING: variation size: " << num << ", RooArgSet size: " << varSize << std::endl;
                      continue;
                    }
                    // if(
                    // 		std::string(_paramList[i].GetName()).find("pdf_")!=std::string::npos ||
                    // 		std::string(_paramList[i].GetName()).find("QCDscale_")!=std::string::npos
                    // 		)
                    {
                      std::cout << "\tnp name: " << std::setw(30) << _paramList[i].GetName() << ", high: " << std::setw(10) << _high[i]-1 << ", low: " << _low[i]-1 << std::endl;
                    }
                  }
                }

            protected:

                // RooListProxy _paramList ;
                // Double_t _nominal;
                // std::vector<double> _low;
                // std::vector<double> _high;
                // std::vector<int> _interpCode;
                // Double_t _interpBoundary;

                // TIterator* _paramIter ;  //! do not persist

                Double_t evaluate() const;

                ClassDef(RooStats::HistFactory::FlexibleInterpVarExt,3) // flexible interpolation
        };
    }
}

#endif
