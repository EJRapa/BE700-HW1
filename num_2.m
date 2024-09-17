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

saveas(gcf,"RawDat.png");