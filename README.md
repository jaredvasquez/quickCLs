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
