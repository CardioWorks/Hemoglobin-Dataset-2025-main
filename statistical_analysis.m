%% SNR Data Statistical Analysis
clear; clc; close all;

%% 1. Data Reading and Preprocessing
data_path = '\Hb-PPG SNR.xlsx';
data_table = readtable(data_path);
snr_data = [data_table.x660nm, data_table.x730nm, data_table.x850nm, data_table.x940nm];
wavelengths = {'660nm', '730nm', '850nm', '940nm'};

[n_subjects, n_wavelengths] = size(snr_data);
fprintf('Data dimensions: %d subjects × %d wavelengths\n', n_subjects, n_wavelengths);

%% 2. Normality Test - Determine Appropriate Statistical Test
use_nonparametric = false;
for i = 1:n_wavelengths
    [h, p] = lillietest(snr_data(:, i));
    if h == 1
        fprintf('%s: Non-normal distribution (p=%.3f)\n', wavelengths{i}, p);
        use_nonparametric = true;
    else
        fprintf('%s: Normal distribution (p=%.3f)\n', wavelengths{i}, p);
    end
end

if use_nonparametric      
    fprintf('\n→ Using non-parametric tests\n\n');
else
    fprintf('\n→ Data follows normal distribution, parametric tests can be used\n\n');
end

%% 3. Calculate median, interquartile range (IQR), and confidence intervals (CI)
rng(42); 
n_boot = 1000;
median_vals = zeros(1, n_wavelengths);
q1_vals = zeros(1, n_wavelengths);  
q3_vals = zeros(1, n_wavelengths);
median_cis = zeros(n_wavelengths, 2);

for i = 1:n_wavelengths
    median_val = median(snr_data(:, i));
    q1_val = quantile(snr_data(:, i), 0.25);
    q3_val = quantile(snr_data(:, i), 0.75);
    
    boot_medians = arrayfun(@(x) median(snr_data(randi(n_subjects, n_subjects, 1), i)), 1:n_boot);
    median_ci = prctile(boot_medians, [2.5, 97.5]);

    median_vals(i) = median_val;
    q1_vals(i) = q1_val;        
    q3_vals(i) = q3_val;
    median_cis(i, :) = median_ci;
    
    fprintf('%s: Median=%.2f, IQR=[%.2f, %.2f], 95%% CI [%.2f, %.2f]\n', ...
        wavelengths{i}, median_val, q1_val, q3_val, median_ci(1), median_ci(2));
end

%% 4. Friedman Overall Test
[p_friedman, tbl] = friedman(snr_data, 1, 'off');
fprintf('χ²(%d)=%.3f, p=%.2e\n', n_wavelengths-1, tbl{2,5}, p_friedman);

if p_friedman < 0.05
    fprintf('Significant differences in SNR between wavelengths\n\n');
else
    fprintf('No significant differences in SNR between wavelengths\n\n');
end

%% 5. Pairwise Comparisons with Effect Size (Wilcoxon Signed-Rank Test)
pairs = nchoosek(1:n_wavelengths, 2);
n_comparisons = size(pairs, 1);
results = zeros(n_comparisons, 4); 

for i = 1:n_comparisons
    idx1 = pairs(i, 1); idx2 = pairs(i, 2);
    [p, ~, stats] = signrank(snr_data(:, idx1), snr_data(:, idx2));
    r_effect = abs(stats.zval) / sqrt(n_subjects);
    results(i, :) = [p, stats.zval, r_effect, stats.signedrank];
    
    fprintf('%s vs %s: ', wavelengths{idx1}, wavelengths{idx2});
    fprintf('p=%.2e, r=%.3f\n', p, r_effect);
end

%% 6. Non-parametric Bootstrap Confidence Intervals for Paired Differences
bootstrap_results = zeros(n_comparisons, 3); 
for i = 1:n_comparisons
    idx1 = pairs(i, 1); idx2 = pairs(i, 2);
    x = snr_data(:, idx1);
    y = snr_data(:, idx2);
    diff_vec = x - y;
    median_diff = median(diff_vec);
    n_boot = 1000;
    boot_medians = zeros(n_boot, 1);
    
    for b = 1:n_boot
        indices = randi(n_subjects, n_subjects, 1);
        boot_sample = diff_vec(indices);
        boot_medians(b) = median(boot_sample);
    end
    ci_lower = prctile(boot_medians, 2.5);
    ci_upper = prctile(boot_medians, 97.5);
    
    bootstrap_results(i, :) = [median_diff, ci_lower, ci_upper];
    fprintf('%s vs %s: ', wavelengths{idx1}, wavelengths{idx2});
    fprintf('Median difference = %.3f, 95%% Bootstrap CI = [%.3f, %.3f]', ...
        median_diff, ci_lower, ci_upper);
    if ci_lower > 0 || ci_upper < 0
        fprintf(' *\n');  
    else
        fprintf('\n');
    end
end

%% 7. Multiple Comparisons Correction (Significance level α=0.05)
p_vals = results(:, 1);
bonf_p = min(p_vals * n_comparisons, 1);
[sorted_p, sort_idx] = sort(p_vals);
fdr_p = zeros(size(p_vals));

for i = 1:length(p_vals)
    fdr_p(sort_idx(i)) = sorted_p(i) * n_comparisons / i;
end
for i = length(fdr_p)-1:-1:1
    if fdr_p(sort_idx(i)) > fdr_p(sort_idx(i+1))
        fdr_p(sort_idx(i)) = fdr_p(sort_idx(i+1));
    end
end

fdr_p = min(fdr_p, 1);
fprintf('\n=== Multiple Comparisons Correction Results ===\n');
fprintf('Comparison Pair\tRaw p-value\tBonferroni\tFDR Corrected\n');
fprintf('──────────────────────────\n');

for i = 1:n_comparisons
    idx1 = pairs(i, 1); idx2 = pairs(i, 2);
    fprintf('%s vs %s\t%.2e\t%.2e\t%.2e\n', ...
        wavelengths{idx1}, wavelengths{idx2}, p_vals(i), bonf_p(i), fdr_p(i));
end

%% 8. Statistical Power Analysis
for i = 1:n_comparisons
    r = results(i, 3);
    ncp = r * sqrt(n_subjects/2);
    z_alpha = norminv(0.975);
    power_val = 1 - normcdf(z_alpha - ncp) + normcdf(-z_alpha - ncp);
    
    idx1 = pairs(i, 1); idx2 = pairs(i, 2);
    fprintf('%s vs %s: Statistical Power=%.3f\n', wavelengths{idx1}, wavelengths{idx2}, power_val);
end

%% 9. Results Summary
fprintf('\n=== Analysis Results Summary ===\n');

fprintf('1. Data Characteristics:\n');
fprintf('   - Sample size: n=%d\n', n_subjects);
fprintf('   - Number of wavelengths: %d wavelengths\n', n_wavelengths);

fprintf('2. Normality Test Results:\n');
normality_count = 0;
for i = 1:n_wavelengths
    [h, p] = lillietest(snr_data(:, i));
    if h == 1
        fprintf('   - %s: Non-normal distribution (p=%.3f)\n', wavelengths{i}, p);
        normality_count = normality_count + 1;
    else
        fprintf('   - %s: Normal distribution (p=%.3f)\n', wavelengths{i}, p);
    end
end
fprintf('   - Summary: %d/%d wavelengths show normal distribution\n', n_wavelengths-normality_count, n_wavelengths);

fprintf('3. Descriptive Statistics (Median[IQR], 95%% Confidence Interval):\n');
for i = 1:n_wavelengths
    fprintf('   - %s: %.2f[%.2f, %.2f], 95%% CI[%.2f, %.2f]\n', ...
        wavelengths{i}, median_vals(i), q1_vals(i), q3_vals(i), median_cis(i, 1), median_cis(i, 2));
end

fprintf('4. Overall Test Results:\n');
fprintf('   - Friedman test: p=%.2e\n', p_friedman);
if p_friedman < 0.05
    fprintf('   - Proceed with pairwise comparisons\n');
end

fprintf('5. Pairwise Difference Analysis (Wilcoxon Signed-Rank Test):\n');
fprintf('   ┌─────────────────────────────────────────────┐\n');
fprintf('   │    Comparison Pair    │    p-value   │  Bonferroni p│   Effect Size r  │  Power │\n');
fprintf('   ├─────────────────────────────────────────────┤\n');
    
for i = 1:n_comparisons
    idx1 = pairs(i, 1); 
    idx2 = pairs(i, 2);

    r = results(i, 3);
    ncp = r * sqrt(n_subjects/2);
    z_alpha = norminv(0.975);
    power_val = 1 - normcdf(z_alpha - ncp) + normcdf(-z_alpha - ncp);

    pair_name = sprintf('%s vs %s', wavelengths{idx1}, wavelengths{idx2});
    bonf_sign = '';
    if bonf_p(i) < 0.05
        bonf_sign = '*';
    end
    
    fprintf('   │     %-15s   │  %-8.2e%-1s   │   %-12.2e │   %-12.3f│ %-5.3f │\n', ...
            pair_name, p_vals(i), bonf_sign, bonf_p(i), r, power_val);
end
fprintf('   ├─────────────────────────────────────────────┤\n');

fprintf('6. Non-parametric Bootstrap Confidence Intervals (Median Differences):\n');
fprintf('   Each comparison shows the median difference (W1 - W2) with 95%% bootstrap confidence interval\n');
fprintf('   ┌─────────────────────────────────────┐\n');
fprintf('   │     Comparison Pair    │    Median Difference (95%% Bootstrap CI)    │\n');
fprintf('   ├─────────────────────────────────────┤\n');

for i = 1:n_comparisons
    idx1 = pairs(i, 1); idx2 = pairs(i, 2);
    median_diff = bootstrap_results(i, 1);
    ci_lower = bootstrap_results(i, 2);
    ci_upper = bootstrap_results(i, 3);
    
    pair_name = sprintf('%s - %s', wavelengths{idx1}, wavelengths{idx2});
   
    if ci_lower > 0 || ci_upper < 0
        significance = '*';
    else
        significance = ' ';
    end
    
    fprintf('   │      %-15s   │       %7.3f [%6.3f, %6.3f] %-1s          │\n', ...
            pair_name, median_diff, ci_lower, ci_upper, significance);
end
fprintf('   ├─────────────────────────────────────┤\n');
fprintf('\nNotes:\n');
fprintf('   - For Wilcoxon test: * indicates p < 0.05 after Bonferroni correction\n');
fprintf('   - For Bootstrap CI: * indicates 95%% CI does not include 0 (significant difference)\n');