# MethylSeqDesign
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
dmr.result<-DMR.analysis(3, Mouse$cov.matrix, Mouse$methyl.matrix)
```
* Predict expected discovery rate(EDR) when there are 20 samples in each group, which is the genome wide power we are interested in. EDR=the number of detected DMR/The total number of DMR
```
set.seed(123)
Estimate.EDR.from.pilot(dmr.result$p.values, model=dmr.result$model, N0=3, target.N=20)
```
## How to interpret the results

```
$`EDR`
[1] 0.9326

$DeclareDMR
[1] 1938

$FDR
[1] 0.05
```

EDR is the predicted genome wide power, DeclareDMR is the number of declared DMRs. FDR is the type I error level estimated by bootstrap procedure, by default we control type I error at 0.05. 
