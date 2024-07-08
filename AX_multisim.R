## AX multiple-simulations conde
# 2021/05/09 This code is intended to run the AX_model_titration some
# number of times (it has random uncertainties included), then interpreted
# with AX_interpretation, and in the end means and standard deviations are
# reported. I would like mean and std both for estimated AXTotal, and also
# uncertainty for any given point on the AX curve. -MLP
clear all
n           <- 1# number of simulations
SampleType  <- 'SW'# options are NaCl, KCl, and SW
S_or_I      <- 33.5

for (i in 1:n){
    fileNameGen <- AX_model_titration(SampleType,S_or_I,i)
}

Files <- dir(fileNameGen)

for (i in 1:n){
    File <- Files(i).name
    FileName <- char(File)
    [Results(1,i),ResultsText] <- AX_interpretation(FileName)
#     disp(ResultsText)
}


for (i in 1:length(Results)){
    AT_BWD(i) <- Results(1,i).AT_BWD
    [AXcurve{i}] <-  Results(1,i).AXcurve(:)
    Length(i) <- length(AXcurve{i})

}
for (i in 1:min(Length)){
    for (j in 1:n){
        AXcurvematrix(i,j) <- AXcurve{j}(i)
    }
}
for (i in 1:min(Length)){
        AXcurvemean(i) <- mean(AXcurvematrix(i,:))
        AXcurvestd(i)  <- std(AXcurvematrix(i,:))
}

figure(2)
hold on
H <- shadedErrorBar(-log10(Results(1).H2T(1:length(AXcurvemean))),AXcurvemean*1e6,AXcurvestd*1e6)
ylabel('\Delta{\itA}_{X} (\mumol kg^{-1})')
xlabel('pH_{T}')
plot(0:1:14,zeros(15),'-k')

xlim([2 10])

# figure(1)
# hold on
# plot(0:1:14,zeros(15),'-k')

ATmean <- mean(AT_BWD)
ATstd  <- std(AT_BWD)
Summary <- sprintf('For %d samples, average AT is %.2f +/- %.2f umol/kg.\r\n',n,ATmean*1e6,ATstd*1e6)
disp(Summary)
