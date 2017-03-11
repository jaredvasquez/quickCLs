#ifndef ASYMCLS_HEADER
#define ASYMCLS_HEADER

#include "CommonHead.h"
#include "RooFitHead.h"
#include "RooStatsHead.h"
#include "utils.h"

using namespace std;
using namespace RooFit;
using namespace RooStats;

class asymCLsTool : public TObject {
  private:
    
    //band configuration
    bool betterBands;
    bool betterNegativeBands;
    bool profileNegativeAtZero;

    //other configuration
    string defaultMinimizer;
    int maxRetries;
    int defaultPrintLevel;
    int defaultStrategy;
    int _optConst;

    bool _nllOffset;
    bool killBelowFatal;
    bool doBlind;
    bool conditionalExpected;
    bool doTilde;
    bool doExp;
    bool doObs;
    bool verbose;
    bool usePredictiveFit;
    bool extrapolateSigma;
    double precision;

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
    double CT=-999;

    RooNLLVar* createNLL(RooDataSet* _data);
    int minimize(RooNLLVar* nll);
    int minimize(RooAbsReal* nll);

    double getLimit(RooNLLVar* nll, double initial_guess = 0);
    double getSigma(RooNLLVar* nll, double mu, double muhat, double& qmu);
    double getQmu(RooNLLVar* nll, double mu);

    void setMu(double mu);
    void saveSnapshot(RooNLLVar* nll, double mu);
    void loadSnapshot(RooNLLVar* nll, double mu);
    void doPredictiveFit(RooNLLVar* nll, double mu1, double m2, double mu);

    double getNLL(RooNLLVar* nll);
    double findCrossing(double sigma_obs, double sigma, double muhat);

    double getQmu95_brute(double sigma, double mu);
    double getQmu95(double sigma, double mu);
    double calcCLs(double qmu_tilde, double sigma, double mu);
    double calcPmu(double qmu_tilde, double sigma, double mu);
    double calcPb(double qmu_tilde, double sigma, double mu);
    double calcDerCLs(double qmu, double sigma, double mu);

    //RooDataSet* makeAsimovData2(RooDataSet* conditioningData, double mu_true, 
    //              double mu_prof = -999, string* mu_str = NULL, string* mu_prof_str = NULL);
    //RooDataSet* makeAsimovData2(RooNLLVar* conditioningNLL, double mu_true, 
    //              double mu_prof = -999, string* mu_str = NULL, string* mu_prof_str = NULL);

    void unfoldConstraints( RooArgSet& initial, RooArgSet& final, RooArgSet& obs, RooArgSet& nuis, int& counter);
    RooDataSet* makeAsimovData( bool doConditional, RooNLLVar* conditioning_nll, double mu_val, string* mu_str = NULL, 
          string* mu_prof_str = NULL, double mu_val_profile = -999, bool doFit = true);


  public:
    asymCLsTool();
    
    // set limit options
    void setBetterBands( bool flag ) { betterBands = flag; };
    void setProfileNegAtZero( bool flag ) { profileNegativeAtZero = flag; };
    void setBetterNegativeBands( bool flag ) { betterNegativeBands = flag; };
    
    void setMinAlgo( string str ) { defaultMinimizer = str; };
    void setPrecision( double val ) { precision = val; };

    void setStrategy( int val ) { defaultStrategy = val; };
    void setOptConst( int val ) { _optConst = val; };
    void setMaxRetries( int val ) { maxRetries = val; };
    void setPrintLevel( int val ) { defaultPrintLevel = val; };

    void setDoTilde( bool flag ) { doTilde = flag; };
    void setDoBlind( bool flag ) { doBlind = flag; }
    void setVerbose( bool flag ) { verbose = flag; };
    void setNLLOffset( bool flag ) { _nllOffset = flag; };
    void setDoExpected( bool flag ) { doExp = flag; };
    void setDoObserved( bool flag ) { doObs = flag && !doBlind; };
    void setCondExpected( bool flag ) { conditionalExpected = flag && !doBlind; };
    void setPredictiveFit( bool flag ) { usePredictiveFit = flag; };
    void setKillBelowFatal( bool flag) { killBelowFatal = flag; };
    
    void runAsymptoticsCLs(  RooWorkspace *ws,
               RooStats::ModelConfig *mc,
               RooDataSet *data,
               string snapshotIn,
		           string folder,
		           string mass,
		           double CL,
		           string outputFile );

    bool checkModel(const RooStats::ModelConfig &model, bool throwOnFail=false) ;
};

#endif
