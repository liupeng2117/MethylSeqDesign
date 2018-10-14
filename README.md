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
## Genome-wide power calculation
* Call the example data. The toy example contains the mouse methylation data from 6 sample in 2 groups, each group has 3 samples.
```
library(MethylSeqDesign)
data(Mouse)
```
* Perform differential methylation analysis. prop specifies the proportion of subsampling. By subsampling, we vary the sequencing depth of pilot data.
```
prop=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)
dmr.result<-DMR.analysis(N0=3,cov.matrix=Mouse$cov.matrix, methyl.matrix=Mouse$methyl.matrix, prop=prop)
```
* Predict expected discovery rate(EDR) when there are N=2,6,12,30,40,50 samples in each group. EDR is the genome wide power we are interested in. EDR=the number of detected DMR/The total number of DMR. This step may take a few minutes
```
N=c(2,6,12,20,30,40,50)
predict.result<-Estimate.EDR.from.pilot(dmr.result, N0=3, target.N=N)
```
* Interpret the results

```
$`EDR`
    2 vs 2 6 vs 6 12 vs 12 20 vs 20 30 vs 30 40 vs 40 50 vs 50
0.1 0.0003 0.1874   0.6601   0.8678   0.9484   0.9763   0.9869
0.2 0.0033 0.4192   0.8175   0.9406   0.9792   0.9902   0.9944
0.3 0.0063 0.5055   0.8542   0.9563   0.9858   0.9935   0.9971
0.4 0.0075 0.5662   0.8787   0.9640   0.9868   0.9937   0.9967
0.5 0.0188 0.5914   0.8888   0.9658   0.9886   0.9947   0.9975
0.6 0.0148 0.6089   0.8947   0.9697   0.9890   0.9956   0.9969
0.7 0.0217 0.6196   0.8982   0.9701   0.9906   0.9953   0.9972
0.8 0.0305 0.6379   0.8996   0.9718   0.9882   0.9952   0.9975
0.9 0.0354 0.6501   0.9079   0.9748   0.9920   0.9960   0.9976
1   0.0350 0.6399   0.9075   0.9716   0.9910   0.9957   0.9974

$DeclareDMR
    2 vs 2 6 vs 6 12 vs 12 20 vs 20 30 vs 30 40 vs 40 50 vs 50
0.1      0    232      820     1078     1178     1213     1226
0.2      4    526     1026     1180     1229     1243     1248
0.3      8    622     1052     1178     1214     1224     1228
0.4      9    708     1099     1206     1235     1243     1247
0.5     24    748     1125     1222     1251     1259     1262
0.6     18    754     1108     1201     1225     1233     1235
0.7     27    769     1115     1204     1230     1236     1238
0.8     38    805     1136     1226     1247     1256     1259
0.9     44    811     1132     1216     1237     1242     1244
1       44    805     1141     1222     1246     1252     1255

$FDR.matrix
    2 vs 2 6 vs 6 12 vs 12 20 vs 20 30 vs 30 40 vs 40 50 vs 50
0.1 0.0000 0.0486   0.0496   0.0499   0.0496   0.0497   0.0497
0.2 0.0000 0.0493   0.0496   0.0497   0.0497   0.0497   0.0497
0.3 0.0000 0.0495   0.0496   0.0498   0.0497   0.0495   0.0498
0.4 0.0033 0.0495   0.0496   0.0497   0.0493   0.0499   0.0498
0.5 0.0328 0.0495   0.0497   0.0499   0.0496   0.0498   0.0498
0.6 0.0269 0.0495   0.0494   0.0495   0.0497   0.0497   0.0499
0.7 0.0286 0.0495   0.0496   0.0497   0.0498   0.0497   0.0496
0.8 0.0306 0.0494   0.0497   0.0495   0.0498   0.0498   0.0498
0.9 0.0334 0.0494   0.0496   0.0496   0.0497   0.0496   0.0493
1   0.0385 0.0497   0.0496   0.0497   0.0496   0.0499   0.0499

attr(,"class")
[1] "list"            "MethylSeqDesign"
```

The result is a MethylSeqDesign object, where the first element is a matrix of predicted genome wide power. DeclareDMR is a matrix of the number of declared differential methylation regions. FDR is the type I error level estimated by bootstrap procedure, by default we control type I error at 0.05. 

## Plot
* Contour plot
```
plotContour(predict.result)
```
![Plot of EDR by change of N and R](https://github.com/liupeng2117/MethylSeqDesign/img/Contour.png)
* 3D plot
```
plot3D(predict.result)
```
![Plot of EDR by change of N and R](https://github.com/liupeng2117/MethylSeqDesign/img/3D.png)
## Study Design
* The pilot data has 40 million reads per sample, and we want to know the minimun budget to achieve 80% power
```
designOptim(predict.result, pilot.R=40, prop, N, targetEDR=0.8)
```

* The pilot data has 40 million reads per sample, and we want to know the maximum power under $20000 budget
```
designOptim(predict.result, pilot.R=40, prop, N, budget=20000)
```