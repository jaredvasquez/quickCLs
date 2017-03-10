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
    
    bool checkModel(const RooStats::ModelConfig &model, bool throwOnFail=false) ;
    int profileToData(ModelConfig *mc, RooAbsData *data);
};

#endif
