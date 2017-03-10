#include "fitTool.h"

fitTool::fitTool() {
  // Default values for the minimizer
  _minAlgo = "Minuit2";
  _minTolerance = 1E-03;
  _minStrat = 1;
  _nllOffset = true;
  _nCPU = 1;
  _optConst = 2;
  _printLevel = 2;
  _useHESSE = true;
  _useMINOS = true;
  _useSIMPLEX = false;
  _saveWS= false;
}

using namespace std;
using namespace RooFit;
using namespace RooStats;

bool fitTool::checkModel(const RooStats::ModelConfig &model, bool throwOnFail) {
// ----------------------------------------------------------------------------------------------------- 
  bool ok = true; 
  std::ostringstream errors;
  std::auto_ptr<TIterator> iter;
  RooAbsPdf *pdf = model.GetPdf(); 
  if (pdf == 0) throw std::invalid_argument("Model without Pdf");
   
  RooArgSet allowedToFloat;
  
  // Check model observables
  if (model.GetObservables() == 0) {
    ok = false; errors << "ERROR: model does not define observables.\n";
    std::cout << errors.str() << std::endl;
    if (throwOnFail) throw std::invalid_argument(errors.str()); else return false;
  } else {
    allowedToFloat.add(*model.GetObservables());
  }

  // Check model parameters of interset
  if (model.GetParametersOfInterest() == 0) {
    ok = false; errors << "ERROR: model does not define parameters of interest.\n";
  } else {
    iter.reset(model.GetParametersOfInterest()->createIterator());
    for (RooAbsArg *a = (RooAbsArg *) iter->Next(); a != 0; a = (RooAbsArg *) iter->Next()) {
      RooRealVar *v = dynamic_cast<RooRealVar *>(a);
      if (!v) { 
        errors << "ERROR: parameter of interest " << a->GetName() << " is a " << a->ClassName() << " and not a RooRealVar\n"; 
        ok = false; continue; 
      }
      if (!pdf->dependsOn(*v)) { 
        errors << "ERROR: pdf does not depend on parameter of interest " << a->GetName() << "\n"; 
        ok = false; continue; 
      }
      allowedToFloat.add(*v);
    }
  }
  
  // Check model nuisance parameters 
  if (model.GetNuisanceParameters() != 0) {
    iter.reset(model.GetNuisanceParameters()->createIterator());
    for (RooAbsArg *a = (RooAbsArg *) iter->Next(); a != 0; a = (RooAbsArg *) iter->Next()) {
      RooRealVar *v = dynamic_cast<RooRealVar *>(a);
      if (!v) { 
        errors << "ERROR: nuisance parameter " << a->GetName() << " is a " << a->ClassName() << " and not a RooRealVar\n"; 
        ok = false; continue; 
      }
      if (v->isConstant()) { 
        errors << "ERROR: nuisance parameter " << a->GetName() << " is constant\n"; 
        ok = false; continue; 
      }
      if (!pdf->dependsOn(*v))
      {
          errors << "WARNING: pdf does not depend on nuisance parameter, removing " << a->GetName() << "\n";
          const_cast<RooArgSet*>(model.GetNuisanceParameters())->remove(*a);
          continue;
      }
      allowedToFloat.add(*v);
    }
  }

  // check model global observables 
  if (model.GetGlobalObservables() != 0) {
    iter.reset(model.GetGlobalObservables()->createIterator());
    for (RooAbsArg *a = (RooAbsArg *) iter->Next(); a != 0; a = (RooAbsArg *) iter->Next()) {
      RooRealVar *v = dynamic_cast<RooRealVar *>(a);
      if (!v) { ok = false; 
        errors << "ERROR: global observable " << a->GetName() << " is a " << a->ClassName() << " and not a RooRealVar\n"; continue; }
      if (!v->isConstant()) { ok = false; errors << "ERROR: global observable " << a->GetName() << " is not constant\n"; continue; }
      if (!pdf->dependsOn(*v)) { errors << "WARNING: pdf does not depend on global observable " << a->GetName() << "\n"; continue; }
    }
  }

  // check the rest of the pdf
  std::auto_ptr<RooArgSet> params(pdf->getParameters(*model.GetObservables()));
  iter.reset(params->createIterator());
  for (RooAbsArg *a = (RooAbsArg *) iter->Next(); a != 0; a = (RooAbsArg *) iter->Next()) {
    if (a->isConstant() || allowedToFloat.contains(*a)) continue;
    if (a->getAttribute("flatParam")) continue;
    errors << "WARNING: pdf parameter " << a->GetName() << " (type " << a->ClassName() << ")"
            << " is not allowed to float (it's not nuisance, poi, observable or global observable)\n";
  }
  iter.reset();
  std::cout << errors.str() << std::endl;
  if (!ok && throwOnFail) throw std::invalid_argument(errors.str());
  return ok;
}



int fitTool::profileToData(ModelConfig *mc, RooAbsData *data){
// ----------------------------------------------------------------------------------------------------- 
  RooAbsPdf *pdf=mc->GetPdf();
  RooWorkspace *w=mc->GetWS();

  RooArgSet funcs = w->allPdfs();
  std::auto_ptr<TIterator> iter(funcs.createIterator());
  for ( RooAbsPdf* v = (RooAbsPdf*)iter->Next(); v!=0; v = (RooAbsPdf*)iter->Next() ) {
    std::string name = v->GetName();
    if (v->IsA() == RooRealSumPdf::Class()) {
      std::cout << "\tset binned likelihood for: " << v->GetName() << std::endl;
      v->setAttribute("BinnedLikelihood", true);
    }
  }
  
  TStopwatch timer1;
  std::cout << "   Building NLL..." << std::endl;
  //RooAbsReal *nll = pdf->createNLL(*data, NumCPU(_nCPU,3), 
  RooAbsReal *nll = pdf->createNLL(*data, NumCPU(_nCPU), 
      Constrain(*mc->GetNuisanceParameters()), GlobalObservables(*mc->GetGlobalObservables()));
  nll->enableOffsetting(1);
  timer1.Stop();
  double t_cpu_ = timer1.CpuTime()/60.;
  double t_real_ = timer1.RealTime()/60.;
  printf("   NLL built in %.2f min (cpu), %.2f min (real)\n", t_cpu_, t_real_);

  //ROOT::Math::MinimizerOptions::SetDefaultTolerance( _minTolerance / 0.001 );
  
  RooMinimizer minim(*nll);
  minim.setStrategy( _minStrat );
  minim.setPrintLevel( _printLevel-1 );
  minim.setProfile(); /* print out time */
  minim.setEps( _minTolerance / 0.001 );
  minim.setOffsetting( _nllOffset );
  if (_optConst > 0) minim.optimizeConst( _optConst );

  int status = 0;
  
  if ( _useSIMPLEX ) {
    cout << endl << "Starting fit with SIMPLEX..." << endl;
    status += minim.simplex();
  }

  // Perform fit with MIGRAD
  status += minim.minimize( _minAlgo );

  if ( _useHESSE ) {
    cout << endl << "Starting fit with HESSE..." << endl;
    status += minim.hesse();
  }

  if ( _useMINOS ) {
    cout << endl << "Starting fit with MINOS..." << endl;
    status += minim.minos( *mc->GetParametersOfInterest() );
  }

  if ( _outputFile != "" ) { 
    cout << endl << "Saving results to " << _outputFile << endl;
    // Create output file and save fit results
    TFile *tout = new TFile( _outputFile, "RECREATE" );
    tout->cd();
    RooFitResult *result = minim.save("fitResult","Fit Results");
    result->Write();
    
    // Get important values to save
    double nllVal = nll->getVal();
    std::map<std::string, double> muMap;
    for (RooLinkedListIter it = mc->GetParametersOfInterest()->iterator(); RooRealVar* POI = dynamic_cast<RooRealVar*>(it.Next());) {
      muMap[POI->GetName()] = POI->getVal();
    }
    
    // Save values to TTree
    TTree *nllTree = new TTree("nllscan", "nllscan");
    nllTree->Branch( "status", &status);
    nllTree->Branch( "nll", &nllVal);
    for (RooLinkedListIter it = mc->GetParametersOfInterest()->iterator(); RooRealVar* POI = dynamic_cast<RooRealVar*>(it.Next());) {
      nllTree->Branch( POI->GetName(), &(muMap[POI->GetName()]) );
    }

    nllTree->Fill();
    nllTree->Write();
    if (_saveWS) {
      RooArgSet everything;
      utils::collectEverything(mc, &everything);
      w->saveSnapshot("postfit", everything);
      w->Write();
    }
  }

  return status;
}


