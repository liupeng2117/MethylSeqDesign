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
pilot.R=1/4
R=c(1/15, 1/10, 1/8, 1/6, 1/4)
dmr.result<-DMR.analysis(N0=3,cov.matrix=Mouse$cov.matrix, methyl.matrix=Mouse$methyl.matrix, pilot.R=pilot.R, R=R)
```
* Predict expected discovery rate(EDR) when there are N=2,6,12,30,40,50 samples in each group. EDR is the genome wide power we are interested in. EDR=the number of detected DMR/The total number of DMR. This step may take a few minutes
```
N=seq(4,50,2)
predict.result<-Estimate.EDR.from.pilot(dmr.result, N0=3, target.N=N)
```
* Interpret the results

```
$`EDR`
                   4 vs 4 6 vs 6 8 vs 8 10 vs 10 12 vs 12 14 vs 14 16 vs 16 18 vs 18 20 vs 20 22 vs 22 24 vs 24 26 vs 26 28 vs 28 30 vs 30
0.0666666666666667 0.1654 0.4640 0.6651   0.7705   0.8482   0.8906   0.9183   0.9379   0.9527   0.9622   0.9701   0.9754   0.9795   0.9833
0.1                0.2431 0.5482 0.7195   0.8202   0.8778   0.9113   0.9346   0.9515   0.9627   0.9699   0.9751   0.9797   0.9832   0.9859
0.125              0.3033 0.6133 0.7577   0.8471   0.8977   0.9291   0.9490   0.9630   0.9702   0.9751   0.9804   0.9838   0.9883   0.9909
0.166666666666667  0.3203 0.6144 0.7648   0.8556   0.8997   0.9303   0.9508   0.9643   0.9724   0.9786   0.9828   0.9865   0.9896   0.9913
0.25               0.3643 0.6402 0.7889   0.8630   0.9049   0.9377   0.9537   0.9654   0.9736   0.9796   0.9836   0.9872   0.9904   0.9918
                   32 vs 32 34 vs 34 36 vs 36 38 vs 38 40 vs 40 42 vs 42 44 vs 44 46 vs 46 48 vs 48 50 vs 50
0.0666666666666667   0.9860   0.9881   0.9898   0.9910   0.9925   0.9933   0.9940   0.9948   0.9953   0.9955
0.1                  0.9883   0.9900   0.9914   0.9923   0.9928   0.9935   0.9943   0.9949   0.9954   0.9957
0.125                0.9923   0.9935   0.9944   0.9949   0.9958   0.9965   0.9970   0.9973   0.9975   0.9979
0.166666666666667    0.9925   0.9932   0.9943   0.9953   0.9957   0.9963   0.9970   0.9974   0.9979   0.9981
0.25                 0.9928   0.9933   0.9945   0.9952   0.9959   0.9964   0.9967   0.9972   0.9975   0.9977

$DeclareDMR
                   4 vs 4 6 vs 6 8 vs 8 10 vs 10 12 vs 12 14 vs 14 16 vs 16 18 vs 18 20 vs 20 22 vs 22 24 vs 24 26 vs 26 28 vs 28 30 vs 30
0.0666666666666667    203    573    822      952     1048     1100     1134     1159     1177     1189     1199     1205     1210     1215
0.1                   304    686    902     1028     1100     1142     1171     1192     1206     1215     1222     1227     1232     1235
0.125                 378    765    946     1058     1120     1159     1185     1202     1211     1217     1224     1228     1234     1237
0.166666666666667     398    764    952     1065     1120     1158     1183     1200     1210     1218     1224     1228     1232     1234
0.25                  459    807    994     1087     1140     1182     1202     1216     1227     1234     1239     1244     1248     1250
                   32 vs 32 34 vs 34 36 vs 36 38 vs 38 40 vs 40 42 vs 42 44 vs 44 46 vs 46 48 vs 48 50 vs 50
0.0666666666666667     1218     1221     1223     1225     1226     1227     1228     1229     1230     1230
0.1                    1238     1240     1242     1243     1244     1245     1246     1246     1247     1247
0.125                  1238     1240     1241     1242     1243     1244     1244     1245     1245     1245
0.166666666666667      1235     1236     1238     1239     1239     1240     1241     1242     1242     1243
0.25                   1251     1252     1253     1254     1255     1256     1256     1256     1257     1257

$FDR.matrix
                   4 vs 4 6 vs 6 8 vs 8 10 vs 10 12 vs 12 14 vs 14 16 vs 16 18 vs 18 20 vs 20 22 vs 22 24 vs 24 26 vs 26 28 vs 28 30 vs 30
0.0666666666666667 0.0474 0.0494 0.0496   0.0494   0.0498   0.0497   0.0495   0.0498   0.0498   0.0494   0.0499   0.0497   0.0497   0.0498
0.1                0.0481 0.0493 0.0496   0.0495   0.0496   0.0497   0.0495   0.0497   0.0496   0.0497   0.0497   0.0498   0.0497   0.0496
0.125              0.0486 0.0493 0.0497   0.0498   0.0497   0.0496   0.0497   0.0498   0.0494   0.0496   0.0496   0.0496   0.0497   0.0500
0.166666666666667  0.0493 0.0495 0.0494   0.0497   0.0498   0.0497   0.0496   0.0499   0.0496   0.0497   0.0499   0.0499   0.0496   0.0494
0.25               0.0493 0.0496 0.0498   0.0498   0.0495   0.0498   0.0497   0.0496   0.0497   0.0495   0.0495   0.0495   0.0498   0.0498
                   32 vs 32 34 vs 34 36 vs 36 38 vs 38 40 vs 40 42 vs 42 44 vs 44 46 vs 46 48 vs 48 50 vs 50
0.0666666666666667   0.0497   0.0497   0.0496   0.0495   0.0496   0.0495   0.0495   0.0496   0.0496   0.0496
0.1                  0.0497   0.0498   0.0498   0.0498   0.0498   0.0498   0.0497   0.0496   0.0495   0.0495
0.125                0.0499   0.0499   0.0499   0.0499   0.0498   0.0497   0.0497   0.0497   0.0496   0.0496
0.166666666666667    0.0494   0.0495   0.0495   0.0495   0.0496   0.0495   0.0496   0.0496   0.0498   0.0498
0.25                 0.0499   0.0499   0.0499   0.0499   0.0500   0.0498   0.0498   0.0498   0.0498   0.0497

attr(,"class")
[1] "list"            "MethylSeqDesign"
```

The result is a MethylSeqDesign object, where the first element is a matrix of predicted genome wide power. DeclareDMR is a matrix of the number of declared differential methylation regions. FDR is the type I error level estimated by bootstrap procedure, by default we control type I error at 0.05. 

## Plot
* Contour plot
```
plotContour(predict.result)
```
![Plot of EDR by change of N and R](https://github.com/liupeng2117/MethylSeqDesign/raw/master/img/Contour.png)
* 3D plot
```
plot3d(predict.result)
```
![Plot of EDR by change of N and R](https://github.com/liupeng2117/MethylSeqDesign/raw/master/img/3D.png)
## Study Design
* Goal1: the minimun budget to achieve 80% power
```
design.res1<-designOptim(predict.result, pilot.R=pilot.R, targetEDR=0.8)
design.res1$res
```
```
$`optim.N`
[1] 10

$optim.R
[1] 0.1

$EDR.target
[1] 0.8

$EDR.achive
[1] 0.8202

$min.budget
[1] 9000
```
* N and R plot, the admissible points are in black, and inadmissible points are in grey. The optimal N and R combination is in red.
```
design.res1$plot1
```
![Plot of EDR by change of N and R](https://github.com/liupeng2117/MethylSeqDesign/raw/master/img/Goal1.1.png)

* Cost and EDR plot, the corresponding points for admissible points are in black, and for inadmissible points are in grey. The corresponding point for optimal N and R combination is in red.
```
design.res1$plot2
```
![Plot of EDR by change of N and R](https://github.com/liupeng2117/MethylSeqDesign/raw/master/img/Goal1.2.png)

* Goal2: we want to know the maximum power under $20000 budget
```
design.res2<-designOptim(predict.result, pilot.R=pilot.R, R=R, N=N, budget=20000)
design.res2$res
```
```
$`optim.N`
[1] 20

$optim.R
[1] 0.125

$budget
[1] 20000

$cost
[1] 19500

$max.EDR
[1] 0.9702
```
* N and R plot, the admissible points are in black, and inadmissible points are in grey. The optimal N and R combination is in red.
```
design.res2$plot1
```
![Plot of EDR by change of N and R](https://github.com/liupeng2117/MethylSeqDesign/raw/master/img/Goal2.1.png)

* Cost and EDR plot, the corresponding points for admissible points are in black, and for inadmissible points are in grey. The corresponding point for optimal N and R combination is in red.
```
design.res2$plot2
```
![Plot of EDR by change of N and R](https://github.com/liupeng2117/MethylSeqDesign/raw/master/img/Goal2.2.png)

