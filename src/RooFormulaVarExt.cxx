/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id: RooFormulaVarExt.cxx 157152 2015-05-19 05:14:58Z yanght $
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

//////////////////////////////////////////////////////////////////////////////
//
// RooRealVar is a generic implementation of a real valued object
// which takes a RooArgList of servers and a C++ expression string defining how
// its value should be calculated from the given list of servers.
// RooRealVar uses a RooFormula object to perform the expression evaluation.
//
// If RooAbsPdf objects are supplied to RooRealVar as servers, their
// raw (unnormalized) values will be evaluated. Use RooGenericPdf, which
// constructs generic PDF functions, to access their properly normalized
// values.
//
// The string expression can be any valid TFormula expression referring to the
// listed servers either by name or by their ordinal list position:
//
//   RooRealVar("gen","x*y",RooArgList(x,y))  or
//   RooRealVar("gen","@0*@1",RooArgList(x,y))
//
// The latter form, while slightly less readable, is more versatile because it
// doesn't hardcode any of the variable names it expects
//


#include "RooFit.h"
#include "Riostream.h"

#include "RooStreamParser.h"
#include "RooNLLVar.h"
#include "RooChi2Var.h"
#include "RooMsgService.h"

#include "TMath.h"

#include "RooFormulaVarExt.h"

ClassImp(RooFormulaVarExt)



//_____________________________________________________________________________
RooFormulaVarExt::RooFormulaVarExt(const RooFormulaVar& other, const char* name) :
    // RooAbsReal(other, name),
    // _actualVars("actualVars",this,other._actualVars),
    // _formula(0),
    // _formExpr(other._formExpr)
    RooFormulaVar(other, name)
{
    // Copy constructor
}




//_____________________________________________________________________________
RooFormulaVarExt::~RooFormulaVarExt()
{
    // Destructor

    // if (_formula) delete _formula ;
}

TString RooFormulaVarExt::RebuildStr(RooFormulaVar* & newVar,
        std::string newVarName,
        bool hardCoded)
{
    TString ret = "";
    if(!hardCoded)
    {
        if ( _formExpr.Contains("@") ) {
            std::cout << "\tIt's already using index, formula: " << _formExpr << std::endl;
        }
        else{
            TString newFormExpr = _formExpr;
            std::map<std::string, std::string> reNameMap;
            /* want to replace the long ones first, so that there is no mis-replacement */
            std::vector<std::string> allOldNames;
            std::vector<int> indice;
            std::vector<int> nameLength;
            int num = _actualVars.getSize();
            for ( int i= 0; i < num; i++ ) {
                TString oldName = _actualVars.at(i)->GetName();
                TString newName = TString::Format("@%d", i);
                indice.push_back(i);
                nameLength.push_back(oldName.Length());
                allOldNames.push_back(oldName.Data());
                reNameMap[oldName.Data()] = newName.Data();
            }
            TMath::Sort(num, &nameLength[0], &indice[0], true);
            for ( int i= 0; i < num; i++ ) {
                int index = indice[i];
                TString oldName = allOldNames[index];
                TString newName = reNameMap[oldName.Data()].c_str();
                newFormExpr = newFormExpr.ReplaceAll(oldName, newName);
                std::cout << "\tWhat??? " << newFormExpr << std::endl;
            }
            ret = newFormExpr;
            // newVar = new RooFormulaVar(newVarName.c_str(), newFormExpr.Data(), _actualVars);
        }
    }
    else
    {
        if ( !_formExpr.Contains("@") ) {
            // std::cout << "\tIt's already using hardCoded, formula: " << _formExpr << std::endl;
        }
        else{
            TString newFormExpr = _formExpr;
            std::map<std::string, std::string> reNameMap;
            /* want to replace the long ones first, so that there is no mis-replacement */
            std::vector<std::string> allOldNames;
            std::vector<int> indice;
            std::vector<int> nameLength;
            int num = _actualVars.getSize();
            for ( int i= 0; i < num; i++ ) {
                TString oldName = _actualVars.at(i)->GetName();
                TString newName = TString::Format("@%d", i);
                indice.push_back(i);
                // nameLength.push_back(oldName.Length());
                // allOldNames.push_back(oldName.Data());
                // reNameMap[oldName.Data()] = newName.Data();
                nameLength.push_back(newName.Length());
                allOldNames.push_back(newName.Data());
                reNameMap[newName.Data()] = oldName.Data();
            }
            TMath::Sort(num, &nameLength[0], &indice[0], true);
            for ( int i= 0; i < num; i++ ) {
                int index = indice[i];
                TString oldName = allOldNames[index];
                TString newName = reNameMap[oldName.Data()].c_str();
                newFormExpr = newFormExpr.ReplaceAll(oldName, newName);
                std::cout << "\tWhat??? " << newFormExpr << std::endl;
            }
            ret = newFormExpr;
        }
    }
    return ret;
}


void RooFormulaVarExt::Rebuild(RooFormulaVar* & newVar,
        std::string newVarName,
        bool hardCoded)
{
    if(!hardCoded)
    {
        if ( _formExpr.Contains("@") ) {
            // std::cout << "\tIt's already using index, formula: " << _formExpr << std::endl;
            return;
        }
        else{
            TString newFormExpr = _formExpr;
            std::map<std::string, std::string> reNameMap;
            /* want to replace the long ones first, so that there is no mis-replacement */
            std::vector<std::string> allOldNames;
            std::vector<int> indice;
            std::vector<int> nameLength;
            int num = _actualVars.getSize();
            for ( int i= 0; i < num; i++ ) {
                TString oldName = _actualVars.at(i)->GetName();
                TString newName = TString::Format("@%d", i);
                indice.push_back(i);
                nameLength.push_back(oldName.Length());
                allOldNames.push_back(oldName.Data());
                reNameMap[oldName.Data()] = newName.Data();
            }
            TMath::Sort(num, &nameLength[0], &indice[0], true);
            for ( int i= 0; i < num; i++ ) {
                int index = indice[i];
                TString oldName = allOldNames[index];
                TString newName = reNameMap[oldName.Data()].c_str();
                newFormExpr = newFormExpr.ReplaceAll(oldName, newName);
                std::cout << "\tWhat??? " << newFormExpr << std::endl;
            }
            newVar = new RooFormulaVar(newVarName.c_str(), newFormExpr.Data(), _actualVars);
        }
    }
    else
    {
        if ( !_formExpr.Contains("@") ) {
            std::cout << "\tIt's already using hardCoded, formula: " << _formExpr << std::endl;
            return;
        }
        else{
            TString newFormExpr = _formExpr;
            std::map<std::string, std::string> reNameMap;
            /* want to replace the long ones first, so that there is no mis-replacement */
            std::vector<std::string> allOldNames;
            std::vector<int> indice;
            std::vector<int> nameLength;
            int num = _actualVars.getSize();
            for ( int i= 0; i < num; i++ ) {
                TString oldName = _actualVars.at(i)->GetName();
                TString newName = TString::Format("@%d", i);
                indice.push_back(i);
                // nameLength.push_back(oldName.Length());
                // allOldNames.push_back(oldName.Data());
                // reNameMap[oldName.Data()] = newName.Data();
                nameLength.push_back(newName.Length());
                allOldNames.push_back(newName.Data());
                reNameMap[newName.Data()] = oldName.Data();
            }
            TMath::Sort(num, &nameLength[0], &indice[0], true);
            for ( int i= 0; i < num; i++ ) {
                int index = indice[i];
                TString oldName = allOldNames[index];
                TString newName = reNameMap[oldName.Data()].c_str();
                newFormExpr = newFormExpr.ReplaceAll(oldName, newName);
                std::cout << "\tWhat??? " << newFormExpr << std::endl;
            }
            newVar = new RooFormulaVar(newVarName.c_str(), newFormExpr.Data(), _actualVars);
        }

    }
}


//_____________________________________________________________________________
RooFormula& RooFormulaVarExt::formula() const
{
    // Return reference to internal RooFormula object

    if (!_formula) {
        _formula = new RooFormula(GetName(),_formExpr,_actualVars) ;
    }
    return *_formula ;
}



//_____________________________________________________________________________
Double_t RooFormulaVarExt::evaluate() const
{
    // Calculate current value of object from internal formula
    return formula().eval(_lastNSet) ;
}



//_____________________________________________________________________________
Bool_t RooFormulaVarExt::isValidReal(Double_t /*value*/, Bool_t /*printError*/) const
{
    // Check if given value is valid
    return kTRUE ;
}



//_____________________________________________________________________________
Bool_t RooFormulaVarExt::redirectServersHook(const RooAbsCollection& newServerList, Bool_t mustReplaceAll, Bool_t nameChange, Bool_t /*isRecursive*/)
{
    // Propagate server change information to embedded RooFormula object
    return _formula ? _formula->changeDependents(newServerList,mustReplaceAll,nameChange) : kFALSE ;
}



//_____________________________________________________________________________
void RooFormulaVarExt::printMultiline(std::ostream& os, Int_t contents, Bool_t verbose, TString indent) const
{
    // Print info about this object to the specified stream.

    // RooAbsReal::printMultiline(os,contents,verbose,indent);
    // if(verbose) {
    //   indent.Append("  ");
    //   os << indent;
    //   formula().printMultiline(os,contents,verbose,indent);
    // }
}



//_____________________________________________________________________________
void RooFormulaVarExt::printMetaArgs(std::ostream& os) const
{
    // Add formula expression as meta argument in printing interface
    // os << "formula=\"" << _formExpr << "\" " ;
}




//_____________________________________________________________________________
Bool_t RooFormulaVarExt::readFromStream(std::istream& /*is*/, Bool_t /*compact*/, Bool_t /*verbose*/)
{
    return kTRUE ;
}



//_____________________________________________________________________________
void RooFormulaVarExt::writeToStream(std::ostream& os, Bool_t compact) const
{
    // Write object contents to given stream
}


//_____________________________________________________________________________
Double_t RooFormulaVarExt::defaultErrorLevel() const
{
    return 1.0 ;
}




