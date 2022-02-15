
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FILE = 0;
CHANNEL = 2;
TRIALS = [1 2];
EPOCH = 1;
PSD = 0;
SPECT = 0;
FILTER = 1;
EPOCH_ANALYIS = 1;
EPOCH_PSD = 1;
PSD_NA = 1;

%can define the threshold here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if FILE == 0
 folder = 'C:\Users\Ryan\Desktop\Kiel MEG Data\Nov-12 Noise Measurement with LED';

 file = {'dat_meg_measurement_165544', 'dat_meg_measurement_170159'};
 
 titlename = {"MEG Noise Measurement at 16:55 Nov-12", ...
              "MEG Noise Measurement at 17:01 Nov-12"};

elseif FILE == 1
 folder = 'C:\Users\Ryan\Desktop\Kiel MEG Data\Nov-21-2019_Initial_MEG_Measurement';
    
 file = {'dat_Human_inside_fiber_optics_off_person_1_7_Hz_3min133118', ...
         'dat_Human_inside_fiber_optics_on_person_1_7Hz_3min133833', ...
         'dat_Human_inside_white_led_on_person_1_7_Hz_3min134811', ...
         'dat_Human_inside_green_led_on_person_1_7_Hz_3min135530', ...
         'dat_Human_inside_green_led_on_person_1_3min_10Hz140043', ...
         'dat_noise_measurement_10_Hz_green140838', ...
         'dat_dat_Human_inside_fiber_optics_off_person_2_10Hz_3min143008', ...
         'dat_dat_Human_inside_fiber_optics_on_person_2_10Hz_3min143408', ...
         'dat_dat_Human_inside_fiber_optics_on_noise_measurment144252'};
     
 titlename = {"Fiber Optic Blinking 7Hz & OPM Off Person 1 [13:31] Nov-21", ...
              "Fiber Optic Blinking 7Hz & OPMs On Person 1 [13:38] Nov-21", ...
              "White LED Blinking 7Hz & OPMs On Person 1 [13:48] Nov-21", ...
              "Green LED Blinking 7Hz & OPMs On Person 1 [13:55] Nov-21", ...
              "Green LED Blinking 10Hz & OPMs On Person 1 [14:04] Nov-21", ...
              "Green LED Blinking 10Hz & OPMs Noise Measurement [14:08] Nov-21", ...
              "Fiber Optic Blinking 10Hz & OPMs Off Person 2 [14:30] Nov-21", ...
              "Fiber Optic Blinking 10Hz & OPMs On Person 2 [14:34] Nov-21", ...
              "Fiber Optic Blinking 10Hz & OPMs Off Person 2 [14:42] Nov-21"};
     
elseif FILE == 2
folder =  'C:\Users\Ryan\Desktop\Kiel MEG Data\Nov-28-2019_2nd_MEG_Measurements';
    
 file = {'dat_white_led_subject_r_10_Hz_on_154528', ...
         'dat_white_led_subject_j_10_Hz_on_160414', ...
         'dat_white_led_subject_j_10_Hz_off_160818', ...
         'dat_green_led_subject_a_10_Hz_off164510', ...
         'dat_green_led_subject_a_10_Hz_on164908', ...
         'dat_white_led_subject_a_10_Hz_on_current_source_on', ...
         'dat_white_led_subject_a_7_Hz_on_current_source_on', ...
         'dat_white_led_subject_a_15_Hz_on_current_source_on', ...
         'dat_white_led_subject_a_off'};
     
 titlename = {"White LED Blinking 10Hz & OPMs on Subject R [15:45] Nov-28", ...
              "White LED Blinking 10Hz & OPMs on Subject J [16:04] Nov-28", ...
              "White LED Off & OPMs On Subject J [16:08] Nov-28", ...
              "Green LED Off & OPMs on Subject A [16:45] Nov-28", ...
              "Green LED Blinking 10Hz & OPMs on Subject A [16:49] Nov-28", ...
              "Current Source - White LED Blinking 10Hz & OPMs on Subject A Nov-28", ...
              "Current Source - White LED Blinking 7Hz & OPMs on Subject A Nov-28", ...
              "Current Source - White LED Blinking 15Hz & OPMs on Subject A Nov-28", ...
              "White LED Off & OPMs on Subject A Nov-28"};
          
 legNames = {'Arduino - W LED 10 Hz - Subject R', 'Arduino - W LED 10 Hz - Subject J', ...
             'Arduino - W No Blinking - Subject J', 'Arduino - G No Blinking - Subject A', ...
             'Arduino - G LED 10 Hz', 'Current Source - 10 Hz', 'Current Source - 7 Hz', ...
             'Current Source - 15 Hz', 'Current Source - No Blinking'};
end

ChannelName = {'Y ', 'V ', 'X ', 'W '};




addpath('C:\Users\Ryan\Downloads\MATLAB Functions')
%x = ones(10,10)
x = cell(9,1);
% Trials 2, 3, 4 have no trigger threshold. Trials 3 & 4 the LEDs were off
% so trial 2 is a problem. Must epoch just like trial 1 or 5.
for i = [1:2]
    %: length(file)
    Matfile = fullfile(folder, [file{i}  '.mat']);
    x{i} = load(Matfile);
    if  (i >= 6) && (FILE == 2)
        x{i}.xin.in(:,5) = x{i}.x.';
        x{i}.x = x{i}.xin;
    end
    
    thresh = mean(x{i}.x.in(:,5));
    % Because there is no trigger for these signals
    if ismember(i, [2 3])
        thresh = mean(x{1}.x.in(:,5));
        x{i}.x.in(:,5) = (x{1}.x.in(:,5));
    elseif ismember(i, 4)
        thresh = mean(x{5}.x.in(:,5));
        x{i}.x.in(:,5) = (x{1}.x.in(:,5));
    end
    xnow = x{i};
    [x{i}.x.trigger, x{i}.x.trigger_change, x{i}.x.epoch.times x{i}.x.epoch.epoch_opm, x{i}.x.epoch.standard_opm, ...
     x{i}.x.epoch.epoch_trigger, x{i}.x.epoch.standard_trigger] = epochdetection(xnow, thresh);         

    %Sampling Rate
    Fs = 10000;

    %Spectrogram
    if SPECT == 1;
        
        % s - 8193 x 3581; 
        % f - 8193 x 1; Real # so (nfft(512*32)/2) + 1 = 8193
        % specT - 1 x 3581; ((length(1800000) - window(10000))/(window-overlap(500)) + 1 = 3581
        % p - 8193 x 3581;
        % Concatenate along the spectT with the same size f, t, p, s
        [s, f, specT, p] = spectrogram(x{i}.x.b(:,CHANNEL)*1e9, 10000, 9500, (512*32), 10000, 'yaxis');
        colormap lbmap
    

    % Making the Figure, 
    
        map = flipud(lbmap(100, 'redblue'));
        figure; h = pcolor(specT, f, db(p)/2); set(h, 'edgecolor', 'none'); ylim([2 40]); colorbar; colormap(map);
        ylabel("Frequency (Hz)", 'FontSize', 14, 'FontWeight', 'bold'); hcb=colorbar; title(hcb, 'Power (dB/Hz)', 'FontSize', 14, 'FontWeight', 'bold'); xlabel("Time (sec)", 'FontSize', 14, 'FontWeight', 'bold');
        
        if (FILE == 2) && ismember(i,[6 7 9])
            % Change the caxis
             caxis([-70 -40]);

        else  
             caxis([-50 10]);
        
        end
        
        title(sprintf("Spectrogram - High Pass 1Hz - Scaled 1E12: %s",titlename{i}),'FontSize', 15);
        ax = gca; ax.FontSize = 14;
    end


%

    % epochdetection.m lines 60 to 86
    % Cut out the last epoch, it is shorter than the others
    % Take the length of epochs and save into a vector
    if EPOCH_ANALYIS == 1
         x{i}.x.epoch_analysis.epoch_length = zeros(length(x{i}.x.epoch.epoch_opm)-1,1) ;
         for kk = 1:length(x{i}.x.epoch.epoch_opm)-1
             x{i}.x.epoch_analysis.epoch_length(kk) = length(x{i}.x.epoch.epoch_opm{kk}); 
         end
         % Find the minimum epoch size and cut every value to that size.
         x{i}.x.epoch_analysis.epoch_same = min(x{i}.x.epoch_analysis.epoch_length(2:end));

        % Create a matrix containing all OPM epochs and another
        % containing all standardized OPM epochs

        x{i}.x.epoch_analysis.epoch_matrix_std = cell(4,1);
        x{i}.x.epoch_analysis.epoch_matrix_opm = cell(4,1);

        for j = 1:4
            x{i}.x.epoch_analysis.epoch_matrix_opm{j} = zeros(x{i}.x.epoch_analysis.epoch_same, length(x{i}.x.epoch.epoch_opm)-1);
            x{i}.x.epoch_analysis.epoch_matrix_std{j} = zeros(x{i}.x.epoch_analysis.epoch_same, length(x{i}.x.epoch.epoch_opm)-1);
            for ik = 1:(length(x{i}.x.epoch.epoch_opm)-1)
                %epoch matrix
                x{i}.x.epoch_analysis.epoch_matrix_opm{j}(:, ik) = x{i}.x.epoch.epoch_opm{ik}(1:x{i}.x.epoch_analysis.epoch_same,j);
                x{i}.x.epoch_analysis.epoch_matrix_std{j}(:, ik) = x{i}.x.epoch.standard_opm{ik}(1:x{i}.x.epoch_analysis.epoch_same,j);
            end
        end
    end

    % Time to average in sec
    Time = 60;
    Fs = 10000;
    trigFreq = 15;
    % Average a number of B epoch trials - The window size
    B = Time * trigFreq;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Preprocessing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Parameters for zero-phase butterworth notch, highpass, and lowpass filters
    [hp2, hp1] = butter(3, 1/(Fs/2), 'high');
    [lp2, lp1] = butter(3, 40/(Fs/2), 'low');


    % x{i}.x.b: HP 1Hz -> LP 40 Hz -> Average
    x{i}.x.epoch_analysis.epoch_matrix_vec = cell(4,1);
    x{i}.x.epoch_analysis.hp_epoch_vec = cell(4,1);
    x{i}.x.epoch_analysis.lp_hp_epoch_vec = cell(4,1);
    
    % j is the 4 channels
    for j = 1:4
        % Store epoch_matrix_opm trials as one long vector to remove edges
        x{i}.x.epoch_analysis.epoch_matrix_vec{j} = x{i}.x.epoch_analysis.epoch_matrix_opm{j}(:);

        % epoch_matrix_vec: Notch 50 Hz -> HP 1Hz -> LP 40 Hz -> Average
        x{i}.x.epoch_analysis.hp_epoch_vec{j} = filtfilt(hp2, hp1, (x{i}.x.epoch_analysis.epoch_matrix_vec{j}));
        x{i}.x.epoch_analysis.lp_hp_epoch_vec{j} = filtfilt(lp2, lp1, (x{i}.x.epoch_analysis.hp_epoch_vec{j}));
    end

    % Seperate into epoched trials again
    x{i}.x.epoch_analysis.lp_hp_epoch_opm = cell(4,1);
    % j is the 4 channels
    for j = 1:4
        x{i}.x.epoch_analysis.lp_hp_epoch_opm{j} = zeros(length(x{i}.x.epoch_analysis.epoch_matrix_opm{CHANNEL}(:,1)), length(x{i}.x.epoch_analysis.epoch_matrix_opm{CHANNEL}(1,:)));
        for ce = 1:length(x{i}.x.epoch_analysis.epoch_length)
            x{i}.x.epoch_analysis.lp_hp_epoch_opm{j}(:, ce) = x{i}.x.epoch_analysis.lp_hp_epoch_vec{j}((1: length(x{i}.x.epoch_analysis.epoch_matrix_opm{CHANNEL}(:,1)))+(( length(x{i}.x.epoch_analysis.epoch_matrix_opm{CHANNEL}(:,1)))*(ce-1)));
        end
    end
end   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSD Thompsons Multitaper Method - Not Averaged Dataset - PREPROCESSED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Multitaper PSD
    [x{i}.x.psd.pxx, x{i}.x.psd.fxx] = pmtm(detrend(x{i}.x.epoch_analysis.lp_hp_epoch_vec{CHANNEL}), 3, (512*4096), 10000);

    if PSD == 1;
        figure; p1 = plot(x{i}.x.psd.fxx, db(x{i}.x.psd.pxx(:)));
        xlim([0 1000]); 
        ylabel("Power/Frequency (dB/Hz)", 'FontSize', 14); xlabel("Frequency (Hz)", 'FontSize', 14);
        title(sprintf("Multitaper PSD - High Pass 1Hz - Scaled 1E12: %s",titlename{i}),'FontSize', 15);
   
    elseif PSD == 2;
    elseif PSD == 0; 
    end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time Domain - Average epochs - NO PREPROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i = 8;
% Time to average in sec
Time = 60;
Fs = 10000;
trigFreq = 15;
% Average a number of B epoch trials - The window size
B = Time * trigFreq;

x{i}.x.epoch_analysis.Avg60 = zeros(length(x{i}.x.epoch_analysis.epoch_matrix_opm{1}(:, 1)), 4);

i = 9;
figure(); hold on;

% j is the channel number
for j = CHANNEL
    x{i}.x.epoch_analysis.Avg60(:, j) = mean(x{i}.x.epoch_analysis.epoch_matrix_opm{j}(:, (901:(900+B))), 2);
    x{i}.x.epoch_analysis.stdAvg60(:, j) = zscore(x{i}.x.epoch_analysis.Avg60(:, j));
    plot(x{i}.x.epoch_analysis.stdAvg60(:, j))
end

% Plot the standardardized triggers on top
for kl = 1:length(x{i}.x.epoch.standard_trigger)-1
hold on; plot(x{i}.x.epoch.standard_trigger{kl});
end

xlabel('Time (ms)'); xticks([0 100 200 300 400 500 600 700]);
xticklabels({'0', '1', '2', '3', '4', '5', '6', '7'});
ylabel('Standardized values');
title(sprintf('Average Epochs 1 - 2 Min, Channel %s', ChannelName{j}, legNames{i}));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time Domain - Average epochs - AFTER PREPROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time to average in sec
i = 9;
Time = 60;
Fs = 10000;
trigFreq = 15;
% Average a number of B epoch trials - The window size
B = Time * trigFreq;

x{i}.x.epoch_analysis.Avg60 = zeros(length(x{i}.x.epoch_analysis.lp_hp_epoch_opm{1}(:, 1)), 4);
i = 9;
figure(); hold on;
for j = CHANNEL
    x{i}.x.epoch_analysis.Avg60(:, j) = mean(x{i}.x.epoch_analysis.lp_hp_epoch_opm{j}(:, (901:(900+B))), 2);
    x{i}.x.epoch_analysis.stdAvg60(:, j) = zscore(x{i}.x.epoch_analysis.Avg60(:, j));
    plot(x{i}.x.epoch_analysis.stdAvg60(:, j))
end

for kl = 1:length(x{i}.x.epoch.standard_trigger)-1
hold on; plot(x{i}.x.epoch.standard_trigger{kl});
end

xlabel('Time (ms)'); xticks([0 100 200 300 400 500 600 700]);
xticklabels({'0', '1', '2', '3', '4', '5', '6', '7'});
ylabel('Standardized values');
title('Average Epochs from 1 Minute to 2 Minutes, Channel Y');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calc Spectrogram for each epoch and concatenate together to make single Spectrogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% s - 8193 x 3581; 
% f - 8193 x 1; Real # so (nfft(512*32)/2) + 1 = 8193
% specT - 1 x 3581; Round((length(1799) - window(50))/(window(50)-overlap(45) (5)) + 1 = 350
% p - 8193 x 3581;
% Concatenate along the spectT with the same size f, t, p, s
clear -regexp ^s ^f ^p ^specT;
CHANNEL = 2;
WINDOW = 50;
NOVERLAP = 45;
NFFT = (512*8);
s = cell(length(x{i}.x.epoch_analysis.epoch_matrix_opm{CHANNEL}),1);
f = cell(length(x{i}.x.epoch_analysis.epoch_matrix_opm{CHANNEL}),1);
p = cell(length(x{i}.x.epoch_analysis.epoch_matrix_opm{CHANNEL}),1);
specT = cell(length(x{i}.x.epoch_analysis.epoch_matrix_opm{CHANNEL}),1);


% For  a specified number of epochs, compute a spectrogram on each epoch
for tt = 1:2699
[s{tt}, f{tt}, specT{tt}, p{tt}] = spectrogram(x{i}.x.epoch_analysis.epoch_matrix_opm{1}(:,tt)*1e9, WINDOW, NOVERLAP, NFFT, 10000, 'yaxis');
specT{tt}(1,:) = specT{tt}(1,:) - specT{tt} (1,1);
end

%epoch interval - total length is up to tt length
eInterval = 2000:2699;
deltaT = specT{1} (1,2) - specT{1} (1,1);
Tlength = length(specT{1}) * deltaT;
specTconcat = zeros(1, length(specT{1})*length(eInterval));
pConcat = zeros( (NFFT/2+1) , length(specT{1})*length(eInterval) );
pConcat(:,1) = p{1}(:,1);

sConcat = zeros( (NFFT/2+1) , length(specT{1})*length(eInterval));
sConcat(:,1) = s{1}(:,1);
% pConcat(:, (1:length(p{1}(1,:)))) = p{1};
% sConcat(:, 1:length(s{1}(1,:))) = s{1};
for te = eInterval
    specTconcat(1, (1:length(specT{te}))+(te-1)*(length(specT{te})))  = (Tlength*(te-1) + specT{te});
    pConcat(:, (1:length(specT{te}))+(te-1)*(length(specT{te}))) =  (Tlength*(te-1) + p{te});
    sConcat(:, (1:length(specT{te}))+(te-1)*(length(specT{te}))) =  (Tlength*(te-1) + s{te});
end
    
map = flipud(lbmap(100, 'redblue'));
figure; h = pcolor(specTconcat, f{1}, db(pConcat)/2); set(h, 'edgecolor', 'none'); ylim([2 40]); colorbar; colormap(map); caxis([-50 -20]);
view(0,90);
axis tight; 
axis([0 10 0 40]);
title(sprintf("Spectrogram - High Pass 1Hz - Scaled 1E12: %s",titlename{i}),'FontSize', 15);
ax = gca; ax.FontSize = 14;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSD Current Source - Blinking 15 Hz & No Blinking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:9
[x{i}.x.psd.pxx, x{i}.x.psd.fxx] = pmtm(detrend(x{i}.x.epoch_analysis.lp_hp_epoch_vec{CHANNEL}), 3, (512*4096), 10000);
end

figure; p1 = plot(x{7}.x.psd.fxx, db(x{7}.x.psd.pxx(:,2)))
hold on; p2 = plot(x{9}.x.psd.fxx, db(x{9}.x.psd.pxx(:,2))); xlim([0 100]);
legend(legNames{7}, legNames{9});
ylabel('Power/Frequency (dB/Hz)'); xlabel('Frequency (Hz)');
ax = gca; ax.FontSize = 15;
title('PSD Current Source - Preprocessed - Blinking 7 Hz & No Blinking', 'FontSize', 17);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSD Sliding Average - Blinking 15 Hz & No Blinking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time to average in sec
i = 9; i2 = 8;
Time = 120;
Fs = 10000;
trigFreq = 15;
% Average a number of B epoch trials - The window size
B = Time * trigFreq;
%Make a sliding window - average B trials with a sliding window S
S = 10;
x{i}.x.epoch_analysis.block_matrix_opm = zeros(666, floor(((2699-900)/S)));
x{i2}.x.epoch_analysis.block_matrix_opm = zeros(666, floor(((2699-900)/S)));

for y = 1:floor(((length(x{i}.x.epoch_analysis.epoch_length) - 1800)/S))
    % Preprocessed
    x{i}.x.epoch_analysis.block_matrix_opm(:, y) = mean(x{i}.x.epoch_analysis.lp_hp_epoch_opm{1}(:, ((1:B) + ((y-1)*S))), 2);
    x{i2}.x.epoch_analysis.block_matrix_opm(:, y) = mean(x{i2}.x.epoch_analysis.lp_hp_epoch_opm{1}(:, ((1:B) + ((y-1)*S))), 2);

    % Not Preprocessed
     x{i}.x.epoch_analysis.block_matrix_npp(:, y) = mean(x{i}.x.epoch_analysis.epoch_matrix_opm{1}(:, ((1:B) + ((y-1)*S))), 2);
     x{i2}.x.epoch_analysis.block_matrix_opm(:, y) = mean(x{i2}.x.epoch_analysis.lp_hp_epoch_opm{1}(:, ((1:B) + ((y-1)*S))), 2);
     
end

% Make one long vector of the sliding window averaged epochs
x{i}.x.epoch_analysis.block_matrix_vec = x{i}.x.epoch_analysis.block_matrix_opm(:);
x{i2}.x.epoch_analysis.block_matrix_vec = x{i2}.x.epoch_analysis.block_matrix_opm(:);

% PSD of the sliding window averaged epochs
[x{i}.x.epoch_analysis.pBlk x{i}.x.epoch_analysis.fBlk] = pmtm(detrend(x{i}.x.epoch_analysis.block_matrix_vec), 3, 512*512, 10000);
figure; p1 = plot(x{i}.x.epoch_analysis.fBlk, db(x{i}.x.epoch_analysis.pBlk(:,1)));
xlim([0 100]);
hold on;

[x{i2}.x.epoch_analysis.pBlk x{i2}.x.epoch_analysis.fBlk] = pmtm(detrend(x{i2}.x.epoch_analysis.block_matrix_vec), 3, 512*512, 10000);
p2 = plot(x{8}.x.epoch_analysis.fBlk, db(x{8}.x.epoch_analysis.pBlk(:,1)));
x{i}.x.epoch_analysis.std_block_opm = zscore(x{i}.x.epoch_analysis.block_matrix_opm);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSD of average epoched signal, copied to the number of epochs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i = 9; i2 = 8; CHANNEL = 1;
averag = mean(x{i}.x.epoch_analysis.epoch_matrix_opm{CHANNEL}, 2);
averagCopy = repmat(averag, 1, 2699);
[dp df] = pmtm(detrend(averagCopy(:)), 3, 512*512, 10000);
figure; p3 = plot(df, db(dp));
xlim([0 100]);

% P Welch
i = 9; i2 = 8; CHANNEL = 1;
averag = mean(x{i}.x.epoch_analysis.epoch_matrix_opm{CHANNEL}, 2);
averagCopy = repmat(averag, 1, 2699);
[pxx1, F1] = pwelch(averagCopy(:), [], [],  1024*16*16, 10000);
figure; plot(F1(:), db(pxx1)); xlim([0 25]);

hold on;

averag2 = mean(x{i2}.x.epoch_analysis.epoch_matrix_opm{CHANNEL}, 2);
averagCopy2 = repmat(averag2, 1, 2699);
[pxx2, F2] = pwelch(averagCopy2(:), [], [],  1024*16*16, 10000);
figure; plot(F2(:), db(pxx2)); xlim([0 25]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spectogram of average epoched signal copied
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change i based on trial you want to run
i = 7; CHANNEL = 2;
% Take the average epochs of trials column-wise
averag = mean(x{i}.x.epoch_analysis.epoch_matrix_opm{CHANNEL}, 2);
averagCopy = repmat(averag, 1, 2699);

[s, f, specT, p] = spectrogram(x{i}.x.epoch_analysis.epoch_matrix_opm{CHANNEL}(:)*1e9, 10000, 9500, (512*32), 10000, 'yaxis');

colormap lbmap
map = flipud(lbmap(100, 'redblue'));
figure; h = pcolor(specT, f, db(p)/2); set(h, 'edgecolor', 'none'); ylim([0 40]); colorbar; colormap(map); caxis([-100 10]);
ylabel("Frequency (Hz)", 'FontSize', 14, 'FontWeight', 'bold'); hcb=colorbar; title(hcb, 'Power (dB/Hz)', 'FontSize', 14, 'FontWeight', 'bold'); xlabel("Time (sec)", 'FontSize', 14, 'FontWeight', 'bold');
title(sprintf("Spectrogram - High Pass 1Hz - Scaled 1E12: %s",titlename{i}),'FontSize', 15);
ax = gca; ax.FontSize = 14;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% log(FFT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change i & i2 based on trial you want to run
i = 9; i2 = 8;
% Change the range ra 000001:100000; 100001:200000
ra = 000001:1800000;
sig = x{i}.x.b(ra, 1);

figure(); plot(20*log10(abs(fft(sig))))
hold on;
sig2 = x{8}.x.b(ra, 1);
plot(20*log10(abs(fft(sig2))))
ylabel('Power/Frequency (dB/Hz)')
xlabel('Time (sec)');
xticks([0 0.5*10^5 1.0*10^5 1.5*10^5 2.0*10^5 2.5*10^5 3.0*10^5]);
xticklabels({'0', '5', '10', '15', '20', '25', '30'});
title(sprintf(['FFT comparing:  ', legNames{8}, '  &  ', legNames{9}]))
legend(legNames{i}, legNames{i2});

i=9;
sig = x{i}.x.b(1:end);
figure(); psd(sig,1024*16*16,10000); xlim([0 25]);
hold on;

i2 = 6;
sig = x{i2}.x.in(1:end,5);
psd(sig,1024*16*16,10000), xlim([0 25])
legend(legNames{i}, legNames{i2});


%%




