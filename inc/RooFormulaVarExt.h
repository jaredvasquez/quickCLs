/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 *    File: $Id: RooFormulaVarExt.h 154634 2013-10-15 19:33:09Z yanght $
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
#ifndef ROO_FORMULA_VAREXT
#define ROO_FORMULA_VAREXT

#include "RooAbsReal.h"
#include "RooFormula.h"
#include "RooFormulaVar.h"
#include "RooArgList.h"
#include "RooListProxy.h"

class RooArgSet ;
class RooFormulaVar ;

class RooFormulaVarExt : public RooFormulaVar {
    public:
        // Constructors, assignment etc
        inline RooFormulaVarExt() : _formula(0), _nset(0) { }
        // RooFormulaVarExt(const RooFormulaVarExt& other, const char* name=0);
        RooFormulaVarExt(const RooFormulaVar& other, const char* name=0);
        virtual TObject* clone(const char* newname) const { return new RooFormulaVarExt(*this,newname); }
        virtual ~RooFormulaVarExt();

        inline Bool_t ok() const { return formula().ok() ; }

        inline RooAbsArg* getParameter(const char* name) const {
            // Return pointer to parameter with given name
            return _actualVars.find(name) ;
        }
        inline RooAbsArg* getParameter(Int_t index) const {
            // Return pointer to parameter at given index
            return _actualVars.at(index) ;
        }

        // I/O streaming interface (machine readable)
        virtual Bool_t readFromStream(std::istream& is, Bool_t compact, Bool_t verbose=kFALSE) ;
        virtual void writeToStream(std::ostream& os, Bool_t compact) const ;

        // Printing interface (human readable)
        virtual void printMultiline(std::ostream& os, Int_t contents, Bool_t verbose=kFALSE, TString indent= "") const ;
        void printMetaArgs(std::ostream& os) const ;

        // Debugging
        void dumpFormula() { formula().dump() ; }

        virtual Double_t defaultErrorLevel() const ;

        TString formulaStr(){
            return _formExpr;
        }
        void Rebuild(RooFormulaVar* & newVar, std::string newName, bool hardcode=false);

        TString RebuildStr(RooFormulaVar* & newVar, std::string newName, bool hardcode=false);
        RooArgList dependVars(){
            return _actualVars;
        }
    protected:

        // Function evaluation
        virtual Double_t evaluate() const ;
        RooFormula& formula() const ;

        // Post-processing of server redirection
        virtual Bool_t redirectServersHook(const RooAbsCollection& newServerList, Bool_t mustReplaceAll, Bool_t nameChange, Bool_t isRecursive) ;

        virtual Bool_t isValidReal(Double_t value, Bool_t printError) const ;

        // RooListProxy _actualVars ;     // Actual parameters used by formula engine
        mutable RooFormula* _formula ; //! Formula engine
        mutable RooArgSet* _nset ;     //! Normalization set to be passed along to contents
        // TString _formExpr ;            // Formula expression string

        ClassDef(RooFormulaVarExt,1) // Real-valued function of other RooAbsArgs calculated by a TFormula expression
};

#endif
