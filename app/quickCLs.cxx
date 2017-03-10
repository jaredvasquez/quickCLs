#include "inc/CommonHead.h"
#include "inc/RooFitHead.h"
#include "inc/RooStatsHead.h"

#include <boost/program_options.hpp>

#include "inc/asymCLsTool.h"
#include "inc/auxUtils.h"
//#include "inc/fitTool.h"

std::string _outputFile = "";
std::string _inputFile = "";

std::string _minAlgo  = "Minuit2";
std::string _dataName = "combData";
std::string _asimovName = "AsimovSB";
std::string _wsName = "combWS";
std::string _mcName = "ModelConfig";

std::string _poiStr = "";
std::string _fixNPStr = "";

//band configuration
bool _betterBands = 1;
bool _betterNegativeBands = 0;
bool _profileNegativeAtZero = 0;

//other configuration
int _maxRetries = 3;
int _printLevel = 0;
int _minStrategy = 1;

bool _doExp = 1;
bool _verbose = 0;
bool _doBlind = 1;
bool _doTilde = 1;
bool _killBelowFatal = 1;
bool _usePredFit = 0;
bool _doObs = 1 && !_doBlind;
bool _conditionalExpected = 1 && !_doBlind;

double _precision = 0.005;
  
string OKGREEN = "\033[92m";
string FAIL = "\033[91m";
string ENDC = "\033[0m";

int main( int argc, char** argv )
{
  namespace po = boost::program_options;
  po::options_description desc( "quickFit options" );
  desc.add_options()
    // IO Options 
    ( "inputFile,f",    po::value<std::string>(&_inputFile),  "Specify the input TFile (REQUIRED)" )
    ( "outputFile,o",   po::value<std::string>(&_outputFile), "Save fit results to output TFile" )
    ( "dataName,d",     po::value<std::string>(&_dataName)->default_value(_dataName),   
                          "Name of the observed dataset" )
    ( "asimovName,a",     po::value<std::string>(&_asimovName)->default_value(_asimovName),   
                          "Name of the Asimov dataset" )
    ( "wsName,w",       po::value<std::string>(&_wsName)->default_value(_wsName),
                          "Name of the workspace" )
    ( "mcName,m",       po::value<std::string>(&_mcName)->default_value(_mcName), 
                          "Name of the model config" )
    // Model Options
    ( "poi,p",          po::value<std::string>(&_poiStr),     "Specify POIs to be used in fit" )
    ( "fixNP,n",        po::value<std::string>(&_fixNPStr),   "Specify NPs to be used in fit" )

    // Band Configuration
    ( "betterBands",    po::value<bool>(&_betterBands)->default_value(_betterBands),
                          "Improve bands by using a more appropriate asimov dataset for those points" )
    ( "betterNegBands", po::value<bool>(&_betterNegativeBands)->default_value(_betterNegativeBands),
                          "Also improve negative bands (not recommended)" )
    ( "setNegAtZero",   po::value<bool>(&_profileNegativeAtZero)->default_value(_profileNegativeAtZero),
                          "Profile Asimov for negative bands at zero (not recommended)" )
    // Minimizer Options
    ( "minStrat",       po::value<int>(&_minStrategy)->default_value(_minStrategy),
                          "Set minimizer strategy" )
    ( "printLevel",     po::value<int>(&_printLevel)->default_value(_printLevel),
                          "Set minimizer print level" )
    ( "maxRetries",     po::value<int>(&_maxRetries)->default_value(_maxRetries),
                          "Number of minimize (fcn) retries before giving up" )
    ( "precision",      po::value<double>(&_precision)->default_value(_precision),
                          "Set \% precision in mu that defines iterative cutoff" )
    ( "verbose",        po::value<bool>(&_verbose)->default_value(_verbose),
                          "Set verbose (very spammy)" )
    // Limit Options
    ( "doExp",          po::value<bool>(&_doExp)->default_value(_doExp),
                          "Compute expected limit" )
    ( "doObs",          po::value<bool>(&_doObs)->default_value(_doObs),
                          "Compute observed limit" )
    ( "doBlind",        po::value<bool>(&_doBlind)->default_value(_doBlind),
                          "Blind analysis from observed limits" )
    ( "doTilde",        po::value<bool>(&_doTilde)->default_value(_doTilde),
                          "Bound mu at zero if true and do the \\tilde{q}_{mu} asymptotics" )
    ( "killBelowFatal", po::value<bool>(&_killBelowFatal)->default_value(_killBelowFatal),
                          "Bound mu at zero if true and do the \\tilde{q}_{mu} asymptotics" )
    ( "usePredFit",     po::value<bool>(&_usePredFit)->default_value(_usePredFit),
                          "(Experimental) extrapolate best fit nuisance parameters based on previous fit results" )
    ( "condExp",        po::value<bool>(&_conditionalExpected)->default_value(_conditionalExpected),
                          "Profiling mode for Asimov data: 0 = conditional MLEs, 1 = nominal MLEs" )
    // Other
    ( "help,h",  "Print help message")
    ;

  po::variables_map vm;
  try
  {
    po::store( po::command_line_parser( argc, argv ).options( desc ).run(), vm );
    po::notify( vm );
  }
  catch ( std::exception& ex )
  {
    std::cerr << "Invalid options: " << ex.what() << std::endl;
    std::cout << "Invalid options: " << ex.what() << std::endl;
    std::cout << "Use manager --help to get a list of all the allowed options"  << std::endl;
    return 999;
  }
  catch ( ... )
  {
    std::cerr << "Unidentified error parsing options." << std::endl;
    return 1000;
  }

  // if help, print help
  if ( !vm.count("inputFile") || vm.count( "help" ) )
  {
    std::cout << "Usage: manager [options]\n";
    std::cout << desc;
    return 0;
  }

  RooMsgService::instance().getStream(1).removeTopic(RooFit::NumIntegration) ;
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Fitting) ;
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Minimization) ;
  RooMsgService::instance().getStream(1).removeTopic(RooFit::InputArguments) ;
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Eval) ;
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  
  // Set fit options
  //fitTool *fitter = new fitTool();
  asymCLsTool *limTool = new asymCLsTool();
  //fitter->setMinAlgo( (TString) _minAlgo );
  //fitter->useHESSE( _useHESSE );
  //fitter->useMINOS( _useMINOS );
  //fitter->useSIMPLEX( _useSIMPLEX );
  //fitter->setNLLOffset( _nllOffset );
  //fitter->setTolerance( _minTolerance );
  //fitter->setStrategy( _minStrategy );
  //fitter->setOptConst( _optConst );
  //fitter->setPrintLevel( _printLevel );
  //fitter->setNCPU( _nCPU );
  //fitter->setOutputFile( (TString) _outputFile );
  //fitter->saveWorkspace( _saveWS );

  // Get workspace, model, and data from file
  TFile *tf = new TFile( (TString) _inputFile );
  RooWorkspace *ws = (RooWorkspace*)tf->Get( (TString) _wsName );
  RooStats::ModelConfig *mc = (RooStats::ModelConfig*)ws->obj( (TString) _mcName );
  RooDataSet *data   = (RooDataSet*)ws->data((TString) _dataName);
  RooDataSet *asimov = (RooDataSet*)ws->data((TString) _asimovName);

  // Prepare model as expected
  utils::setAllConstant( mc->GetGlobalObservables(), true );
  utils::setAllConstant( mc->GetNuisanceParameters(), false );
  utils::setAllConstant( mc->GetParametersOfInterest(), true );

  // Sanity checks on model 
  cout << "Performing sanity checks on model..." << endl;
  bool validModel = limTool->checkModel( *mc, true );
  cout << "Sanity checks on the model: " << (validModel ? "OK" : "FAIL") << endl;

  // Fix nuisance narameters
  if ( vm.count("fixNP") ) {
    cout << endl << "Fixing nuisance parameters : " << endl;
    std::vector<std::string> fixNPStrs = auxUtils::Tokenize( _fixNPStr, "," );
    for( unsigned int inp(0); inp < fixNPStrs.size(); inp++ ) {
      RooAbsCollection *fixNPs = mc->GetNuisanceParameters()->selectByName( (TString) fixNPStrs[inp]);
      for (RooLinkedListIter it = fixNPs->iterator(); RooRealVar* NP = dynamic_cast<RooRealVar*>(it.Next());) {
        cout << "   Fixing nuisance parameter " << NP->GetName() << endl;
        NP->setConstant( kTRUE );
      }
    }
  }

  // Prepare parameters of interest
  RooArgSet fitPOIs;
  if ( vm.count("poi") ) {
    cout << endl << "Preparing parameters of interest :" << endl;
    std::vector<std::string> poiStrs = auxUtils::Tokenize( _poiStr, "," );
    for( unsigned int ipoi(0); ipoi < poiStrs.size(); ipoi++ ) {
      std::vector<std::string> poiTerms = auxUtils::Tokenize( poiStrs[ipoi], "=" );
      TString poiName = (TString) poiTerms[0];

      // check if variable is in workspace
      if (not ws->var(poiName))  {
        cout << FAIL << "Variable " << poiName << " not in workspace. Skipping." << ENDC << endl;
        continue;
      }

      // set variable for fit
      fitPOIs.add( *(ws->var(poiName)) );
      if (poiTerms.size() > 1) {
        std::vector<std::string> poiVals = auxUtils::Tokenize( poiTerms[1], "_" );
        if (poiVals.size() == 3) {
          ws->var(poiName)->setConstant( kFALSE );
          ws->var(poiName)->setVal( std::stof(poiVals[0]) );
          ws->var(poiName)->setRange( std::stof(poiVals[1]), std::stof(poiVals[2]) );
        } else {
          ws->var(poiName)->setVal( std::stof(poiVals[0]) );
        }
        cout << "   ";
        ws->var(poiName)->Print();
      } else {
        ws->var(poiName)->setConstant( kFALSE );
        cout << "   ";
        ws->var(poiName)->Print();
      }
    }
  } else {
    RooRealVar *firstPOI = (RooRealVar*)mc->GetParametersOfInterest()->first();
    cout << endl << "No POIs specified. Will only float the first POI " << firstPOI->GetName() << endl;
    firstPOI->setConstant(kFALSE);
    cout << "   ";
    firstPOI->Print();
    fitPOIs.add( *firstPOI );
  }

  // Fitting 
  //TStopwatch timer;
  //cout << endl << "Starting fit..." << endl;
  //int status = 0; //fitter->profileToData( mc, data ); // Perform fit

  //limTool->runAsymptoticsCLs(infile.Data(),wname.Data(),mname.Data(),dname.Data(),aname.Data(),pname.Data(),mass,CL,option);

  cout << endl << "Start limit setting:" << endl;
  limTool->runAsymptoticsCLs( ws, mc, data, asimov, "folder", "125.09", 0.95, "" );

  //timer.Stop();
  //double t_cpu_ = timer.CpuTime()/60.;
  //double t_real_ = timer.RealTime()/60.;
  //printf("\nAll fits done in %.2f min (cpu), %.2f min (real)\n", t_cpu_, t_real_);

  //string STATMSG = (status) ? "\033[91m STATUS FAILED \033[0m" : "\033[92m STATUS OK \033[0m" ;

  // Print summary 
  //cout << endl << "  Fit Summary of POIs (" << STATMSG << ")" << endl;
  //cout << "------------------------------------------------" << endl;
  //for (RooLinkedListIter it = fitPOIs.iterator(); RooRealVar* POI = dynamic_cast<RooRealVar*>(it.Next());) {
  //  if (POI->isConstant()) continue;
  //  POI->Print();
  //}

  //if (status) {
  //  cout << FAIL << endl;
  //  cout << "   *****************************************" << endl;
  //  cout << "          WARNING: Fit status failed.       " << endl;
  //  cout << "   *****************************************" << ENDC << endl;
  //}

  cout << endl;
  return 1;
}
