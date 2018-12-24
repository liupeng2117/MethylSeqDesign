# MethylSeqDesign
An R package for the study design and power calculation issues of Methyl-Seq experiments.
## Installation
You can install this package by copying and paste the following code into R
```
install.packages("devtools") ## In case you didn't install devtools
devtools::install_github("gitlongor/pi0") # install the dependency pi0
devtools::install_github("liupeng2117/MethylSeqDesign")
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
* Predict expected discovery rate(EDR) when there are N=2,6,12,30,40,50 samples per group. EDR is the genome wide power we are interested in. EDR=the number of detected DMR/The total number of DMR. This step may take a few minutes
```
N=seq(4,50,2)
predict.result<-Estimate.EDR.from.pilot(dmr.result, N0=3, target.N=N)
```
* Interpret the results

```
$`EDR`
      4 vs 4 6 vs 6 8 vs 8 10 vs 10 12 vs 12 14 vs 14 16 vs 16 18 vs 18 20 vs 20 22 vs 22 24 vs 24 26 vs 26 28 vs 28 30 vs 30
0.067 0.1731 0.4776 0.6715   0.7725   0.8455   0.8853   0.9156   0.9361   0.9511   0.9606   0.9682   0.9745   0.9790   0.9832
0.1   0.2602 0.5674 0.7258   0.8279   0.8792   0.9149   0.9374   0.9536   0.9637   0.9702   0.9742   0.9781   0.9827   0.9858
0.125 0.2777 0.5875 0.7377   0.8331   0.8907   0.9233   0.9460   0.9601   0.9676   0.9738   0.9790   0.9830   0.9871   0.9895
0.167 0.3248 0.6164 0.7613   0.8460   0.8927   0.9244   0.9453   0.9581   0.9677   0.9745   0.9798   0.9845   0.9881   0.9898
0.25  0.3709 0.6426 0.7874   0.8641   0.9049   0.9362   0.9535   0.9655   0.9723   0.9793   0.9838   0.9860   0.9891   0.9911
      32 vs 32 34 vs 34 36 vs 36 38 vs 38 40 vs 40 42 vs 42 44 vs 44 46 vs 46 48 vs 48 50 vs 50
0.067   0.9857   0.9878   0.9893   0.9907   0.9919   0.9930   0.9937   0.9946   0.9949   0.9956
0.1     0.9882   0.9900   0.9910   0.9921   0.9937   0.9944   0.9952   0.9955   0.9965   0.9972
0.125   0.9909   0.9920   0.9929   0.9938   0.9950   0.9953   0.9957   0.9963   0.9965   0.9969
0.167   0.9917   0.9928   0.9939   0.9946   0.9956   0.9964   0.9967   0.9971   0.9974   0.9978
0.25    0.9924   0.9937   0.9943   0.9952   0.9954   0.9960   0.9970   0.9971   0.9973   0.9977

$DeclareDMR
      4 vs 4 6 vs 6 8 vs 8 10 vs 10 12 vs 12 14 vs 14 16 vs 16 18 vs 18 20 vs 20 22 vs 22 24 vs 24 26 vs 26 28 vs 28 30 vs 30
0.067    216    597    840      966     1057     1107     1144     1170     1188     1201     1210     1218     1224     1229
0.1      322    701    897     1023     1086     1130     1158     1178     1191     1199     1204     1209     1214     1218
0.125    344    729    916     1034     1106     1146     1174     1192     1202     1209     1216     1221     1226     1229
0.167    405    768    949     1054     1112     1152     1178     1194     1206     1215     1221     1227     1231     1234
0.25     457    791    969     1064     1114     1153     1174     1189     1197     1206     1211     1214     1218     1220
      32 vs 32 34 vs 34 36 vs 36 38 vs 38 40 vs 40 42 vs 42 44 vs 44 46 vs 46 48 vs 48 50 vs 50
0.067     1232     1235     1237     1239     1240     1242     1242     1243     1244     1244
0.1       1221     1223     1225     1226     1228     1229     1230     1230     1232     1233
0.125     1230     1232     1233     1234     1236     1236     1236     1237     1237     1238
0.167     1236     1237     1239     1240     1241     1242     1242     1243     1243     1244
0.25      1222     1223     1224     1226     1226     1226     1228     1228     1228     1229

$FDR.matrix
      4 vs 4 6 vs 6 8 vs 8 10 vs 10 12 vs 12 14 vs 14 16 vs 16 18 vs 18 20 vs 20 22 vs 22 24 vs 24 26 vs 26 28 vs 28 30 vs 30
0.067 0.0474 0.0495 0.0494   0.0497   0.0496   0.0498   0.0496   0.0496   0.0494   0.0497   0.0496   0.0497   0.0498   0.0497
0.1   0.0490 0.0496 0.0495   0.0498   0.0497   0.0498   0.0495   0.0495   0.0497   0.0495   0.0497   0.0498   0.0497   0.0499
0.125 0.0484 0.0496 0.0497   0.0494   0.0496   0.0495   0.0496   0.0496   0.0499   0.0496   0.0498   0.0497   0.0497   0.0496
0.167 0.0491 0.0493 0.0495   0.0495   0.0496   0.0496   0.0498   0.0497   0.0499   0.0495   0.0496   0.0497   0.0497   0.0496
0.25  0.0495 0.0495 0.0494   0.0497   0.0495   0.0498   0.0498   0.0497   0.0496   0.0496   0.0498   0.0500   0.0497   0.0496
      32 vs 32 34 vs 34 36 vs 36 38 vs 38 40 vs 40 42 vs 42 44 vs 44 46 vs 46 48 vs 48 50 vs 50
0.067   0.0496   0.0496   0.0498   0.0499   0.0499   0.0498   0.0498   0.0497   0.0497   0.0496
0.1     0.0495   0.0495   0.0493   0.0494   0.0497   0.0498   0.0498   0.0497   0.0498   0.0498
0.125   0.0495   0.0497   0.0497   0.0497   0.0497   0.0497   0.0497   0.0497   0.0497   0.0497
0.167   0.0498   0.0497   0.0496   0.0496   0.0496   0.0497   0.0497   0.0497   0.0496   0.0497
0.25    0.0497   0.0496   0.0496   0.0497   0.0498   0.0498   0.0497   0.0497   0.0497   0.0497

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
design.res1<-designOptim(predict.result, pilot.R=pilot.R, R=R, N=N, targetEDR=0.8)
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

