%%     Signal length trade-off
clear; clc;
%%%%%%----Read all mat files----%%%%%%
cd '\Hb-PPG Dataset';
files = dir('*.mat');
filename = '1.mat';
data = load(filename);  

%%%%%%----Extract signals from four wavelengths----%%%%%%
nm_660 = data.PPGdata.nm_660;   nm_730 = data.PPGdata.nm_730;
nm_850 = data.PPGdata.nm_850;   nm_940 = data.PPGdata.nm_940;

%%%%%%----Visualize original signals from four wavelengths----%%%%%%
figure(1);
sgtitle('Original Signals'); 
subplot(4, 1, 1); plot(nm_660, 'r'); title('660nm');
subplot(4, 1, 2); plot(nm_730, 'g'); title('730nm');
subplot(4, 1, 3); plot(nm_850, 'b'); title('850nm');
subplot(4, 1, 4); plot(nm_940, 'k'); title('940nm');

fs = 200; 
%%%%%%----Choose to discard from front or end----%%%%%%
cut_position = 'front';  % 'front': discard front, 'end': discard end

%%%%%%----Choose the length of signal segment to discard (seconds)----%%%%%%
cut_seconds = 15;  % Discard 5/10/15/20/25/30 seconds
cut_samples = cut_seconds * fs;

%%%%%%----Check if signal length is sufficient----%%%%%%
min_required_length = 30 * fs;
total_length = length(nm_660);

if total_length > min_required_length
    if strcmp(cut_position, 'front')
        if total_length > (cut_samples + min_required_length)
            nm_660_cut = nm_660(cut_samples+1:end);
            nm_730_cut = nm_730(cut_samples+1:end);
            nm_850_cut = nm_850(cut_samples+1:end);
            nm_940_cut = nm_940(cut_samples+1:end);
            fprintf('Signal cutting completed: Discarded %d seconds from front, kept %d seconds of valid signal\n', ...
                cut_seconds, floor(length(nm_660_cut)/fs));
        else
            error('After discarding %d seconds from front, remaining signal length is less than 30 seconds!', cut_seconds);
        end
    else
        if total_length > (cut_samples + min_required_length)
            nm_660_cut = nm_660(1:end-cut_samples);
            nm_730_cut = nm_730(1:end-cut_samples);
            nm_850_cut = nm_850(1:end-cut_samples);
            nm_940_cut = nm_940(1:end-cut_samples);
            fprintf('Signal cutting completed: Discarded %d seconds from end, kept %d seconds of valid signal\n', ...
                cut_seconds, floor(length(nm_660_cut)/fs));
        else
            error('After discarding %d seconds from end, remaining signal length is less than 30 seconds!', cut_seconds);
        end
    end
else
    error('Total signal length insufficient! Original signal length: %d seconds, need at least %d seconds of valid signal', ...
        floor(total_length/fs), 30);
end

%%%%%%----Visualize valid signals from four wavelengths----%%%%%%
figure(2);
if strcmp(cut_position, 'front')
    position_str = 'front';
else
    position_str = 'end';
end
sgtitle('Valid Signals'); 
subplot(4, 1, 1); plot(nm_660_cut, 'r'); title(['660nm (Discarded ' position_str ' ' num2str(cut_seconds) ' seconds)']);
subplot(4, 1, 2); plot(nm_730_cut, 'g'); title(['730nm (Discarded ' position_str ' ' num2str(cut_seconds) ' seconds)']);
subplot(4, 1, 3); plot(nm_850_cut, 'b'); title(['850nm (Discarded ' position_str ' ' num2str(cut_seconds) ' seconds)']);
subplot(4, 1, 4); plot(nm_940_cut, 'k'); title(['940nm (Discarded ' position_str ' ' num2str(cut_seconds) ' seconds)']);

%%%%%%----Save valid signals and output----%%%%%%
PPGdata = struct('nm_660', nm_660_cut, ...
                 'nm_730', nm_730_cut, ...
                 'nm_850', nm_850_cut, ...
                 'nm_940', nm_940_cut);
[~, name, ~] = fileparts(filename);  
output_filename = [name, 'cut.mat'];  

save(output_filename, 'PPGdata');
fprintf('Valid signals saved to %s\n', output_filename);