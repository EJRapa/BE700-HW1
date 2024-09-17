%% Load and separate good/poor prognosis
clear, close all

HW1_dat = readtable("HW1_data _2024.csv");
good_table = table2array(HW1_dat(strcmp(HW1_dat.prognosis,'GOOD'),1:30));
bad_table = table2array(HW1_dat(strcmp(HW1_dat.prognosis,'POOR'),1:30));

%% perform two-sample Student's ttest

[NLL,p_val,CI,stats] = ttest2(good_table,bad_table);

p_log = -log10(p_val);
z_score = norminv(p_val/2);   % Assuming one-tailed

[p_val,idx] = sort(p_val);

%% Sort columns based on results

HW1_dat = HW1_dat(:,[idx,31]);

writetable(HW1_dat,'sorted_HW1_dat.csv');

%% Create Graphs

rowName = HW1_dat.Properties.VariableNames(1:30);

% Visualization of gene significance
bar(rowName,NLL);
ylim([0 2]);
xlabel("Gene Type");
ylabel("Null Hypothesis Rejection | 1 = true | 0 = false");
title("Comparison of Good Prognosis Versus Poor Prognosis Patients Significance");
saveas(gcf,"Significance.png");

figure 

% Raw data comparison of good and poor prognosis
good_mean = mean(good_table);
poor_mean = mean(bad_table);

SE_good = std(good_table)/size(good_table,1);
SE_poor = std(bad_table)/size(bad_table,1);

bar(rowName,good_mean,'g');

hold on 

er_good = errorbar(1:length(good_mean),good_mean,SE_good);
er_good.LineStyle = 'none';

hold on 

bar(rowName,poor_mean,'r');

hold on 

er_bad = errorbar(1:length(poor_mean),poor_mean,SE_poor);
er_bad.LineStyle = 'none';

legend(["Good Prognosis","","Poor Prognosis"]);

title("Average Microarray Expression At Specific Gene Based on Prognosis");
xlabel("Gene");
ylabel("Gene Expression Level");

saveas(gcf,"RawDat.png");

% Plot p-values 

figure
bar(rowName,p_val);

hold on

yline(0.05,"--");
ylim([0 0.055])

legend(["P-value","Failure to Reject Null Hypothesis"]);

title("Visualization of p-values for Each Gene");
ylabel("Probability of Observing Test Statistic or More Extreme");
xlabel("Gene");

saveas(gcf,"pvals.png");

%% kNN Test Success

% rows are number of folds, cols are k (1,3,5,10,20)
kNN_correct = [
    87.36, 86.21, 85.06, 90.80, 87.36 ;  % 5-fold
    87.36, 88.51, 83.91, 87.36, 88.51 ;  % 10-fold
    0    , 0    , 0    , 0    , 90.80 ;  % 15-fold
    0    , 0    , 0    , 0    , 88.51 ;  % 20-fold
    0    , 0    , 0    , 0    , 89.66 ;  % 25-fold
    0    , 0    , 0    , 0    , 89.66 ;  % 30-fold
];

x = [1,3,5,10,20];

figure
bar(x,kNN_correct);
ylim([80 91]);

hold on

[high_val,index] = max(kNN_correct,[],"all");
scatter(x(round(index/size(kNN_correct,2))),high_val,'r',"filled");

legend(["5-fold","10-fold","15-fold","20-fold","25-fold","30-fold","Max Correct Classification"]);

saveas(gcf,"wekaplot.png");