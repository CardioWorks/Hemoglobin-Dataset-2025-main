%% SNR Calculation Process

clc, clear;
cd '\Hb-PPG Dataset';

%%%%%%----Read all mat files----%%%%%%
files = dir('*.mat');
fileNames = {files.name}; 
fileIDs = str2double(extractBefore(fileNames, '.mat')); 
[~, sortedIndex] = sort(fileIDs);
files = files(sortedIndex); 

%%%%%%----Create a table to record SNR results----%%%%%%
results = cell(length(files), 5); 
results(:, 1) = num2cell(fileIDs(sortedIndex)); 

%%%%%%----Read and process data----%%%%%%
channelNames = {'nm_660', 'nm_730', 'nm_850', 'nm_940'};
for i = 1:length(files)-251
    filename = files(i).name;
    fprintf('Processing file: %s\n', filename);
    loadedData = load(filename);
    if isfield(loadedData, 'PPGdata')
        data = loadedData.PPGdata;
    else
        error('The file %s does not contain the PPGdata field.', filename);
    end
    snrValues = zeros(1, 4);
    % figure;
    for channel = 1:4
        fieldName = channelNames{channel};
        y = data.(fieldName);
        if any(~isfinite(y))
            warning('Data contains non-finite values. Replacing with NaNs.');
            y(~isfinite(y)) = NaN; 
        end
        y = y(~isnan(y));
        if isempty(y)
            error('Data is empty after removing non-finite values.');
        end
        fprintf('File: %s, Channel: %s\n', filename, fieldName);
        %%%%%%----Preprocessing----%%%%%%
        [b_hp, a_hp] = butter(1, 0.5 / (200 / 2), 'high'); 
        y_baseline_removed = filtfilt(b_hp, a_hp, y);
        %%%%%%----Calculate SNR----%%%%%%
        [b_bp, a_bp] = cheby2(1, 10, [0.5, 10] / (200 / 2), 'bandpass');
        y_pass = filtfilt(b_bp, a_bp, y_baseline_removed);
        [b_sp, a_sp] = cheby2(1, 10, [0.5, 10] / (200 / 2), 'stop');
        y_stop = filtfilt(b_sp, a_sp, y_baseline_removed);
        snrValues(channel) = snr( y_pass , y_stop );
        %%%%%%----Waveform visualization----%%%%%%
        % subplot(2, 2, channel); 
        % plot(y_baseline_removed, 'b'); hold on;
        % plot(y_pass, 'r');
        % title(['Channel ', num2str(channel)]);
        % xlabel('Sample Index');
        % ylabel('Amplitude');
        % legend('Baseline Removed', 'Filtered');
        % grid on;
        % hold off;
    end
    %%%%%%----SNR results retained and exported----%%%%%%
    snrValues = round(snrValues, 2);
    results(i, 2:5) = num2cell(snrValues);
end
resultsTable = cell2table(results, ...
    'VariableNames', {'ID', 'nm_660', 'nm_730', 'nm_850', 'nm_940'});
outputFile = 'SNR.xlsx';
writetable(resultsTable, outputFile);
fprintf('All files have been processed. Results saved to %s\n', outputFile);