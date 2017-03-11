# quickCLs
Lightweight tool for quickly extracting asymptotic CLs limits 

## Installation
```
source setup.sh
make
```

## Fitting your workspace
The script requires you provide a dataset but by default will keep this dataset blinded
and generate it's own asimov based on the nominal pdf's values or a designated snapshot 
( setting the first POI to zero ). 

Setting limit on `mu_tH`
```
quickCLs -f filename.root -d dataset -p mu_tH=1_0_50
```

Setting limit on `mu_ggH` while fixing `mu_VBF=2`
```
quickCLs -f filename.root -d dataset -p mu_ggH=1_-5_5,mu_VBF=2
```

Setting limit on `mu_ttH` while profiling `mu_ggH`, `mu_VBF`, and `mu_VH`
```
quickCLs -f filename.root -d dataset -p mu_ttH=1_0_5,mu_ggH=1_0_5,mu_VBF=1_0_5,mu_VH=1_0_5
```

Fixing all systematics with `ATLAS_*` prefix
```
quickCLs -f filename.root -d dataset -p mu_ZH=1_0_5 -n ATLAS_*
```

Additional features can be discovered by asking for help
```
Usage: manager [options]
quickFit options:
  -f [ --inputFile ] arg                Specify the input TFile (REQUIRED)
  -o [ --outputFile ] arg               Save fit results to output TFile
  -d [ --dataName ] arg                 Name of the observed dataset
  -w [ --wsName ] arg (=combWS)         Name of the workspace
  -m [ --mcName ] arg (=ModelConfig)    Name of the model config
  -s [ --snapshot ] arg                 Load snapshot for generating Asimov
                                        dataset.
  -p [ --poi ] arg                      Specify POIs to be used in fit
  -n [ --fixNP ] arg                    Specify NPs to be used in fit
  --betterBands arg (=1)                Improve bands by using a more
                                        appropriate asimov dataset for those
                                        points
  --betterNegBands arg (=0)             Also improve negative bands (not
                                        recommended)
  --setNegAtZero arg (=0)               Profile Asimov for negative bands at
                                        zero (not recommended)
  --minStrat arg (=0)                   Set minimizer strategy
  --printLevel arg (=-1)                Set minimizer print level
  --maxRetries arg (=3)                 Number of minimize (fcn) retries before
                                        giving up
  --precision arg (=0.005)              Set % precision in mu that defines
                                        iterative cutoff
  --verbose arg (=0)                    Set verbose (very spammy)
  --nllOffset arg (=1)                  Set NLL offset
  --optConst arg (=2)                   Set optimize constant
  --doExp arg (=1)                      Compute expected limit
  --doObs arg (=0)                      Compute observed limit
  --doBlind arg (=1)                    Blind analysis from observed limits
  --doTilde arg (=1)                    Bound mu at zero if true and do the
                                        \tilde{q}_{mu} asymptotics
  --killBelowFatal arg (=1)             Bound mu at zero if true and do the
                                        \tilde{q}_{mu} asymptotics
  --usePredFit arg (=0)                 (Experimental) extrapolate best fit
                                        nuisance parameters based on previous
                                        fit results
  --condExp arg (=0)                    Profiling mode for Asimov data: 0 =
                                        conditional MLEs, 1 = nominal MLEs
  -h [ --help ]                         Print help message
```
