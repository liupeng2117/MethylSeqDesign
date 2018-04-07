 20# MethylSeqDesign
An R package for the study design and power calculation issues of Methyl-Seq experiments.
## Installation
You can install this package by copying and paste the following code into R
```
install.packages("devtools") ## In case you didn't install devtools
library(devtools)
devtools::install_github("gitlongor/pi0") # install the dependency pi0
install_github("liupeng2117/MethylSeqDesign")
```
## How to use MethylSeqDesign
* Call the example data. The toy example contains the mouse methylation data from 6 sample in 2 groups, each group has 3 samples.
```
data(Mouse)
```
* Perform differential methylation analysis
```
dmr.result<-DMR.analysis(c(3,3),Mouse$cov.matrix,Mouse$methyl.matrix)
```
* Predict expected discovery rate(EDR) when there are 20 samples in each group, which is the genome wide power we are interested in. EDR=the number of detected DMR/The total number of DMR
```
Estimate.EDR.from.pilot(dmr.result[[1]],model=dmr.result[[2]],N0=3,target.N=20)
```
## How to interpret the results

```
$EDR
[1] 0.9330382

$DeclareDMR
[1] 159.3

$FDR
[1] 0.04767486
```

EDR is the predicted genome wide power, DeclareDMR is the number of declared DMRs. FDR is the type I error level estimated by bootstrap procedure, by default we control type I error at 0.05. 
