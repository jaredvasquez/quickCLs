#ifndef FITUTIL_HEADER
#define FITUTIL_HEADER

#include "CommonHead.h"
#include "RooFitHead.h"
#include "RooStatsHead.h"
#include "utils.h"

using namespace std;
using namespace RooFit;
using namespace RooStats;

class fitTool : public TObject{
  private:
    TString _minAlgo, _outputFile;
    float _minTolerance;
    bool _nllOffset, _useHESSE, _useMINOS, _useSIMPLEX, _saveWS;
    int _minStrat, _optConst, _printLevel, _nCPU; 

    //band configuration
    bool betterBands           = 1; // (recommendation = 1) improve bands by using a more appropriate asimov dataset for those points
    bool betterNegativeBands   = 0; // (recommendation = 0) improve also the negative bands
    bool profileNegativeAtZero = 0; // (recommendation = 0) profile asimov for negative bands at zero

    //other configuration
    string defaultMinimizer    = "Minuit2";     // or "Minuit"
    int defaultPrintLevel      = 0;            // Minuit print level
    int defaultStrategy        = 0;             // Minimization strategy. 0-2. 0 = fastest, least robust. 2 = slowest, most robust
    bool killBelowFatal        = 1;             // In case you want to suppress RooFit warnings further, set to 1
    bool doBlind               = 1;             // in case your analysis is blinded
    bool conditionalExpected   = 1 && !doBlind; // Profiling mode for Asimov data: 0 = conditional MLEs, 1 = nominal MLEs
    bool doTilde               = 1;             // bound mu at zero if true and do the \tilde{q}_{mu} asymptotics
    bool doExp                 = 1;             // compute expected limit
    bool doObs                 = 1 && !doBlind; // compute observed limit
    double precision           = 0.005;         // % precision in mu that defines iterative cutoff
    bool verbose               = 0;             // 1 = very spammy
    bool usePredictiveFit      = 0;             // experimental, extrapolate best fit nuisance parameters based on previous fit results
    bool extrapolateSigma      = 0;             // experimantal, extrapolate sigma based on previous fits
    int maxRetries             = 3;             // number of minimize(fcn) retries before giving up

    TString Option="";

    //don't touch!
    map<RooNLLVar*, double> map_nll_muhat;
    map<RooNLLVar*, double> map_muhat;
    map<RooDataSet*, RooNLLVar*> map_data_nll;
    map<RooNLLVar*, string> map_snapshots;
    map<RooNLLVar*, map<double, double> > map_nll_mu_sigma;
    RooWorkspace* w = NULL;
    ModelConfig* mc = NULL;
    RooDataSet* data = NULL;
    RooRealVar* firstPOI = NULL;
    RooNLLVar* asimov_0_nll = NULL;
    RooNLLVar* obs_nll = NULL;
    int nrMinimize=0;
    int direction=1;
    int global_status=0;
    double target_CLs=0.05;

    double Hmass=-1;
    double CT=-998;

    double getLimit(RooNLLVar* nll, double initial_guess = 0);
    double getSigma(RooNLLVar* nll, double mu, double muhat, double& qmu);
    double getQmu(RooNLLVar* nll, double mu);
    void saveSnapshot(RooNLLVar* nll, double mu);
    void loadSnapshot(RooNLLVar* nll, double mu);
    void doPredictiveFit(RooNLLVar* nll, double mu1, double m2, double mu);
    RooNLLVar* createNLL(RooDataSet* _data);
    double getNLL(RooNLLVar* nll);
    double findCrossing(double sigma_obs, double sigma, double muhat);
    void setMu(double mu);
    double getQmu95_brute(double sigma, double mu);
    double getQmu95(double sigma, double mu);
    double calcCLs(double qmu_tilde, double sigma, double mu);
    double calcPmu(double qmu_tilde, double sigma, double mu);
    double calcPb(double qmu_tilde, double sigma, double mu);
    double calcDerCLs(double qmu, double sigma, double mu);
    int minimize(RooNLLVar* nll);
    int minimize(RooAbsReal* nll);
    //RooDataSet* makeAsimovData2(RooDataSet* conditioningData, double mu_true, double mu_prof = -999, string* mu_str = NULL, string* mu_prof_str = NULL);
    //RooDataSet* makeAsimovData2(RooNLLVar* conditioningNLL, double mu_true, double mu_prof = -999, string* mu_str = NULL, string* mu_prof_str = NULL);

    void unfoldConstraints(RooArgSet& initial, RooArgSet& final, RooArgSet& obs, RooArgSet& nuis, int& counter);
    RooDataSet* makeAsimovData(bool doConditional, RooNLLVar* conditioning_nll, double mu_val, string* mu_str = NULL, string* mu_prof_str = NULL, double mu_val_profile = -999, bool doFit = true);


  public:
    fitTool();
    void useHESSE( bool flag ) { _useHESSE = flag; };
    void useMINOS( bool flag ) { _useMINOS = flag; };
    void useSIMPLEX( bool flag ) { _useSIMPLEX = flag; };
    void setNLLOffset( bool flag ) { _nllOffset = flag; };
    void saveWorkspace( bool flag ) { _saveWS = flag; };
    void setTolerance( float val ) { _minTolerance = val; };
    void setNCPU( int val ) { _nCPU = val; };
    void setStrategy( int val ) { _minStrat = val; };
    void setOptConst( int val ) { _optConst = val; };
    void setPrintLevel( int val ) { _printLevel = val; };
    void setOutputFile( TString str ) { _outputFile = str; };
    void setMinAlgo( TString str ) { _minAlgo = str; };
    
    //main
    void runAsymptoticsCLs(const char* infile,
               const char* workspaceName,
               const char* modelConfigName,
               const char* dataName,
               const char* asimovDataName,
               string folder,
               string mass,
               double CL,
               TString option);

    bool checkModel(const RooStats::ModelConfig &model, bool throwOnFail=false) ;
    int profileToData(ModelConfig *mc, RooAbsData *data);
};

#endif
