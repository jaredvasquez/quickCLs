# quickCLs
Lightweight tool for quickly extracting asymptotic CLs limits 

## Installation
```
source setup.sh
make
```

## Fitting your workspace
Simple fit with MIGRAD
```
quickCLs -f filename.root -d dataset -p mu_ggH=1_-5_5,mu_VBF=1_-5_5
```

Fitting `mu_ggH` but fixing `mu_VBF=2`
```
quickCLs -f filename.root -d dataset -p mu_ggH=1_-5_5,mu_VBF=2
```

Fixing all systematics with `ATLAS_*` prefix
```
quickCLs -f filename.root -d dataset -p mu_ggH=1_-5_5,mu_VBF=1_-5_5 -n ATLAS_*
```

Including HESSE + MINOS fits
```
quickCLs -f filename.root -d dataset -p mu_ggH=1_-5_5,mu_VBF=1_-5_5 --hesse 1 --minos 1
```

Outputting fit results to `output.root`
```
quickCLs -f filename.root -d dataset -p mu_ggH=1_-5_5,mu_VBF=1_-5_5 -o output.root
```
