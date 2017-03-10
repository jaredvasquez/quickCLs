# quickFit
Tool for quickly fitting a data in a RooWorkspace

## Installation
```
source setup.sh
make
```

## Fitting your workspace
Simple fit with MIGRAD
```
quickFit -f filename.root -d dataset -p mu_ggH=1_-5_5,mu_VBF=1_-5_5
```

Fitting `mu_ggH` but fixing `mu_VBF=2`
```
quickFit -f filename.root -d dataset -p mu_ggH=1_-5_5,mu_VBF=2
```

Fixing all systematics with `ATLAS_*` prefix
```
quickFit -f filename.root -d dataset -p mu_ggH=1_-5_5,mu_VBF=1_-5_5 -n ATLAS_*
```

Including HESSE + MINOS fits
```
quickFit -f filename.root -d dataset -p mu_ggH=1_-5_5,mu_VBF=1_-5_5 --hesse 1 --minos 1
```

Outputting fit results to `output.root`
```
quickFit -f filename.root -d dataset -p mu_ggH=1_-5_5,mu_VBF=1_-5_5 -o output.root
```
