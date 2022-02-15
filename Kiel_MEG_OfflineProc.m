
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FILE = 2;
CHANNEL = 1;
TRIALS = [8 9];
EPOCH = 1;
PSD = 0;
SPECT = 0;
FILTER = 1;
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
     
 titlename = {"White LED Blinking 5Hz & OPMs on Subject R [15:45] Nov-28", ...
              "White LED Blinking 5Hz & OPMs on Subject J [16:04] Nov-28", ...
              "White LED Off & OPMs On Subject J [16:08] Nov-28", ...
              "Green LED Off & OPMs on Subject A [16:45] Nov-28", ...
              "Green LED Blinking 10Hz & OPMs on Subject A [16:49] Nov-28", ...
              "Current Source - White LED Blinking 10Hz & OPMs on Subject A Nov-28", ...
              "Current Source - White LED Blinking 7Hz & OPMs on Subject A Nov-28", ...
              "Current Source - White LED Blinking 15Hz & OPMs on Subject A Nov-28", ...
              "White LED Off & OPMs on Subject A Nov-28"};
end

ChannelName = {"Y", "V", "X", "W"};



addpath('C:\Users\Ryan\Downloads\MATLAB Functions')
%x = ones(10,10)
for i = (TRIALS)
    %: length(file)
    Matfile = fullfile(folder, [file{i}  '.mat']);
    load(Matfile);
    if  (i >= 6) && (FILE == 2)
        xin.in(:,5) = x.';
        x = xin;
    end
    
    
    thresh = mean(x.in(:,5));
    [ x.trigger, x.trigger_change, x.epoch.times x.epoch.epoch_opm, x.epoch.standard_opm, ...
      x.epoch.epoch_trigger, x.epoch.standard_trigger] = epochdetection(x, thresh);         


    %Sampling Rate
    Fs = 10000;
    
    %High Pass data at 1 Hz
    [b, a] = butter(4, 1/(Fs/2), 'high'); pause(.2);
    x.filter.scale = x.b*1e12;
    x.filter.hp = filtfilt(b, a, (x.b));
    x.filter.scale_hp = filtfilt(b, a, (x.filter.scale));

    %Thompsons Multitaper Method
    if PSD == 1;
        [pxx fxx] = pmtm(detrend(x.filter.hp)*1e12, 3, (512*1024), 10000);
        figure; p1 = plot(fxx, db(pxx(:,CHANNEL)));
        xlim([0 100]); 
        ylabel("Power/Frequency (dB/Hz)", 'FontSize', 14); xlabel("Frequency (Hz)", 'FontSize', 14);
        title(sprintf("Multitaper PSD - High Pass 1Hz - Scaled 1E12: %s",titlename{i}),'FontSize', 15);
    elseif PSD == 2;
    elseif PSD == 0; 
     end

    %Spectrogram
    if SPECT == 1;
        
        % s - 8193 x 3581; 
        % f - 8193 x 1; Real # so (nfft(512*32)/2) + 1 = 8193
        % specT - 1 x 3581; ((length(1800000) - window(10000))/(window-overlap(500)) + 1 = 3581
        % p - 8193 x 3581;
        % Concatenate along the spectT with the same size f, t, p, s
        [s, f, specT, p] = spectrogram(x.filter.hp(:,2)*1e12, 10000, 9500, (512*32), 10000, 'yaxis');
        colormap lbmap
    

    % Making the Figure, 
        if (FILE == 2) && ismember(i,[6 7 9)
            map = flipud(lbmap(100, 'redblue'));

            figure; h = pcolor(specT, f, db(p)/2); set(h, 'edgecolor', 'none'); ylim([2 40]); colorbar; colormap(map); caxis([-70 -40]);
            ylabel("Frequency (Hz)", 'FontSize', 14, 'FontWeight', 'bold'); hcb=colorbar; title(hcb, 'Power (dB/Hz)', 'FontSize', 14, 'FontWeight', 'bold'); xlabel("Time (sec)", 'FontSize', 14, 'FontWeight', 'bold');
            title(sprintf("Spectrogram - High Pass 1Hz - Scaled 1E12: %s",titlename{i}),'FontSize', 15);
            ax = gca; ax.FontSize = 14;
        else
            map = flipud(lbmap(100, 'redblue'));
            figure; h = pcolor(specT, f, db(p)/2); set(h, 'edgecolor', 'none'); ylim([2 40]); colorbar; colormap(map); caxis([-50 10]);
            ylabel("Frequency (Hz)", 'FontSize', 14, 'FontWeight', 'bold'); hcb=colorbar; title(hcb, 'Power (dB/Hz)', 'FontSize', 14, 'FontWeight', 'bold'); xlabel("Time (sec)", 'FontSize', 14, 'FontWeight', 'bold');
            title(sprintf("Spectrogram - High Pass 1Hz - Scaled 1E12: %s",titlename{i}),'FontSize', 15);
            ax = gca; ax.FontSize = 14;
        end
    end


%

% epochdetection.m lines 60 to 86
% Cut out the last epoch, it is shorter than the others
% Take the length of epochs and save into a vector
if EPOCH_ANALYIS == 1
     epoch_length = zeros(length(x.epoch.epoch_opm)-1,1) ;
     for i = 1:length(x.epoch.epoch_opm)-1
         epoch_length(i) = length(x.epoch.epoch_opm{i}); 
     end
     % Find the minimum epoch size and cut every value to that size.
     epoch_same = min(epoch_length(2:end));
        
    % Create a matrix containing all OPM epochs and another
    % containing all standardized OPM epochs

    epoch_matrix_std = cell(4,1);
    epoch_matrix_opm = cell(4, 1);

    for j = 1:4
        epoch_matrix_std{j} = zeros(epoch_same, length(x.epoch.epoch_opm)-1);
        epoch_matrix_opm{j} = zeros(epoch_same, length(x.epoch.epoch_opm)-1);
        for i = 1:length(x.epoch.epoch_opm)-1
            %epoch matrix
            epoch_matrix_std{j}(:, i) = x.epoch.standard_opm{i}(1:epoch_same,j);
            epoch_matrix_opm{j}(:, i) = x.epoch.epoch_opm{i}(1:epoch_same,j);
        end
    end
end

% Average all OPM epochs and create a multitaper PSD
if EPOCH_PSD == 1
%     avg_epoch = mean(epoch_matrix_opm{1},2);
    [pAvg fAvg] = pmtm(epoch_matrix_opm{1}(:),3, 1024*2048, 10000);
    figure; p1 = plot(fAvg, db(pAvg(:,1)));
    xlim([0 100]);
    hold on; 
    for jk  = 2:4
        avg_epoch = mean(epoch_matrix_opm{jk},2);
        [pAvg fAvg] = pmtm(avg_epoch,3, 512*512, 10000);
        plot(fAvg, db(pAvg(:,1)));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time to average in sec
Time = 60;
Fs = 10000;
trigFreq = 15;
% Average a number of B epoch trials - The window size
B = Time * trigFreq;

% Parameters for zero-phase butterworth notch, highpass, and lowpass filters
[nf2, nf1] = butter(3, [49 51]./(Fs/2), 'stop');
[hp2, hp1] = butter(3, 1/(Fs/2), 'high');
[lp2, lp1] = butter(3, 40/(Fs/2), 'low');

% Store epoch_matrix_opm trials as one long vector to remove edges
epoch_matrix_vec = epoch_matrix_opm{CHANNEL}(:);

% epoch_matrix_vec: Notch 50 Hz -> HP 1Hz -> LP 40 Hz -> Average
nf_epoch_vec = filtfilt(nf2, nf1, (epoch_matrix_vec));
hp_nf_epoch_vec = filtfilt(hp2, hp1, (epoch_matrix_vec));
lp_hp_nf_epoch_vec = filtfilt(hp2, hp1, (epoch_matrix_vec));

% Seperate the epoched trials again
lp_hp_nf_epoch_opm = zeros(666,2699);
for ce = 1:2699
    lp_hp_nf_epoch_opm(:, ce) = lp_hp_nf_epoch_vec((1:666)+((666)*(ce-1)));
end

% Average epochs over a set time
figure; hold on;
for j = 1
    Avg60 = mean(lp_hp_nf_epoch_opm(:, (901:(900+B))), 2)
    stdAvg60 = zscore(Avg60);
    plot(stdAvg60)
end
for i = 1:length(x.epoch.standard_trigger)-1
hold on; plot(x.epoch.standard_trigger{i});
end
xlabel('Time (ms)'); xticks([0 100 200 300 400 500 600 700]);
xticklabels({'0', '1', '2', '3', '4', '5', '6', '7'});
ylabel('Standardized values');
title('Average Epochs from 1 Minute to 2 Minutes, Channel Y');
 
%Make a sliding window - average B trials every 50 samples
S = 10;
block_matrix_opm = zeros(666, floor(((2699-900)/S)));
for y = 1:floor(((2699-900)/S))
    block_matrix_opm(:, y) = mean(lp_hp_nf_epoch_opm(:, (1:B) + ((y-1)*S)), 2);
end


block_matrix_vec = block_matrix_opm(:);
[pBlk fBlk] = pmtm(detrend(block_matrix_vec)*1e12,3, 512*512, 10000);
figure; p1 = plot(fBlk, db(pBlk(:,1)));
xlim([0 100]);

std_block_opm = zscore(block_matrix_opm);


block_matrix_opm = zeros(length(lp_hp_nf_epoch_opm(:,1)), floor(length(lp_hp_nf_epoch_opm(1, :))/B));
% 1:(floor(2699/900))
for e = 1:(floor(length(lp_hp_nf_epoch_opm(1, :))/B))
    %i.e. @ e = 1 ==> mean(epoch_opm((1:60)+((0)*1), :))
    block_matrix_opm(:, e) = mean(lp_hp_nf_epoch_opm(:, (1:B)+((e-1)*B)), 2);
end
block_matrix_vec = block_matrix_opm(:);
[pBlk fBlk] = pmtm(block_matrix_vec)*1e12,3, 512*512, 10000);
figure; p1 = plot(fBlk, db(pBlk(:,1)));
xlim([0 100]);

std_block_opm = zscore(block_matrix_opm);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% make a spectrogram of each epoch and concatenate together
% s - 8193 x 3581; 
% f - 8193 x 1; Real # so (nfft(512*32)/2) + 1 = 8193
% specT - 1 x 3581; Round((length(1799) - window(50))/(window(50)-overlap(45) (5)) + 1 = 350
% p - 8193 x 3581;
% Concatenate along the spectT with the same size f, t, p, s
clear -regexp ^s ^f ^p ^specT;
s = cell(length(epoch_matrix),1);
f = cell(length(epoch_matrix),1);
p = cell(length(epoch_matrix),1);
specT = cell(length(epoch_matrix),1);
deltaT = cell(length(epoch_matrix),1);


for i = 1:100
[s{i}, f{i}, specT{i}, p{i}] = spectrogram(epoch_matrix(i,:), 50, 45, (512*32), 10000, 'yaxis');
specT{i} (1,:) = specT{i} (1,:) - specT{i} (1,1);
end


e = 150;
deltaT = specT{1} (1,2) - specT{1} (1,1);
Tlength = length(specT{1}) * deltaT;
specTconcat = specT{1};
pConcat = zeros(8193, e*350);
sConcat = zeros(8193, e*350);
pConcat(:, (1:length(p{1}(1,:)))) = p{1};
sConcat(:, 1:length(s{1}(1,:))) = s{1};
for i = 2:e
    specTconcat = [specTconcat (Tlength*(i-1) + specT{i})];
    pConcat = [pConcat p{i}];
    sConcat = [sConcat s{i}];
end
    
map = flipud(lbmap(100, 'redblue'));
figure; h = pcolor(specTconcat, f{1}, db(pConcat)/2); set(h, 'edgecolor', 'none'); ylim([2 40]); colorbar; colormap(map); caxis([-50 -20]);
view(0,90);
axis tight; 
axis([0 inf 2 40]);
title(sprintf("Spectrogram - High Pass 1Hz - Scaled 1E12: %s",titlename{i}),'FontSize', 15);
ax = gca; ax.FontSize = 14;







%%


S={};
for i = 1:3
    %Choose the index of the subjectt and condition (ie. Subject 1 is "3")
    subjectNo = i;
    %Naming the subject
    subj{i} = ['S' num2str(subjectTag(subjectNo))];
    S.(subj{i}) = [];
    for k = 1:(length(conditionTag)-1)
        %Loading specific subject task data and saving it in S structure
        %S structure is S.(Subject#).(Task Name)
        %i.e Subject 3 Open & Close data is S.S3.OpenClose
        conditionNo = k;
        content = [num2str(subjectTag(subjectNo)) '_' conditionTag{conditionNo}];
        Matfile = fullfile(folder, ['fullrun' content  '.mat']);
        
        %Load Data
        load(Matfile);
        S.(subj{i}).(task{k}) = megDat;
        
        %Scale up data 1e12
        S.(subj{i}).(task{k}).data = S.(subj{i}).(task{k}).data * 1e12;
        
        %High Pass data at .1Hz
        S.(subj{i}).(task{k}).data = filtfilt(b,a,(S.(subj{i}).(task{k}).data)); pause(.2);
        %Mu Bandpass
        S.(subj{i}).(task{k}).mu = filtfilt(d,c,(S.(subj{i}).(task{k}).data)); pause(.2);;
        %Beta Bandpass
        S.(subj{i}).(task{k}).beta = filtfilt(f,e,(S.(subj{i}).(task{k}).data)); pause(.2)
        
        
        S.(subj{i}).(task{k}).data = double(S.(subj{i}).(task{k}).data);
       
        %Band Power
        S.(subj{i}).(task{k}).mu = db((S.(subj{i}).(task{k}).mu).^2);
        S.(subj{i}).(task{k}).beta = db((S.(subj{i}).(task{k}).beta).^2); pause(0.2);
        
%                 S.(subj{i}).(task{k}).data = db((S.(subj{i}).(task{k}).data).^2);
        
        %Read the VMRK file because it has the trial type and time labels
        vmrkfile = fullfile(folder2, ['fullrun' content '.vmrk']); pause(0.2)
        %Remove the extra VMRK text
        %Create v structure to organize task types and time labels based on the overall task
        v.(subj{i}).(task{k}) = fileread(vmrkfile); pause(0.2);
        v.(subj{i}).(task{k}) = v.(subj{i}).(task{k})(464:end); pause(0.2);
        v.(subj{i}).(task{k}) = splitlines(v.(subj{i}).(task{k})); pause(0.2);
        v.(subj{i}).(task{k})(end) = [];
        
        pause(0.4);
        
        mkr = {};
        mkrsplit = {};
        ttype = [];
        timeL = [];
        for m = 1:length(v.(subj{i}).(task{k}))
            %Seperate the task type (ttype) and time labels (timeL) into 2
            %seperate matrices
            mkr(m) = regexp(v.(subj{i}).(task{k}){m}, '(?<=S  ).*(?=,1,0)', ('match'));
            mkrsplit{m} = strsplit(mkr{m}, ',');
            ttype = [ttype str2num(mkrsplit{m}{1})];
            timeL = [timeL str2num(mkrsplit{m}{2})];
        end
        
        events = [ttype; timeL]';
        v.(subj{i}).(task{k}) = events;
        
        %Add latencies to S structure. Add time point 1 for the 1st time stamp
        S.(subj{i}).(task{k}).latency = [1; v.(subj{i}).(task{k})(:,2)]; 
        
        %I want to save the data points in each trial. {[4999x1]},
        %{[3000x1]}...etc. {[Trial 1]}, {[Rest 1]}, ...etc.
        taskTemp = S.(subj{i}).(task{k});
        taskTemp.trials = {};
        taskTemp.muT = {};
        taskTemp.betaT ={};
        for t = 1:length(taskTemp.latency) - 1
            taskTemp.trials{t} = ...
                taskTemp.data([taskTemp.latency(t) : taskTemp.latency(t+1)-1], :);
            
            taskTemp.muT{t} = taskTemp.mu([taskTemp.latency(t) : taskTemp.latency(t+1)-1], :);
            taskTemp.betaT{t} = taskTemp.beta([taskTemp.latency(t) : taskTemp.latency(t+1)-1], :);
        end
        
        taskTemp.trials{t+1} = taskTemp.data(taskTemp.latency(t+1)+1:end, :);
        taskTemp.muT{t+1} = taskTemp.mu(taskTemp.latency(t+1)+1:end, :);
        taskTemp.betaT{t+1} = taskTemp.beta(taskTemp.latency(t+1)+1:end, :);
        
       %To create an average of the task trials
       %1st Make all task trials the same number
%        minCountEv = min(cellfun(@length, taskTemp.trials(1:2:end)));
%        minCountOd = min(cellfun(@length, taskTemp.trials(2:2:end)));
%        
%        %Only keep numbers that are
%        for j = 1:2:80
%            taskTemp.trials{j} = taskTemp.trials{j}(1:minCountEv,:);
%            taskTemp.muT{j} = taskTemp.muT{j}(1:minCountEv,:);
%            taskTemp.betaT{j} = taskTemp.betaT{j}(1:minCountEv,:);
%        end
%        
%        for n = 2:2:80
%            taskTemp.trials{n} = taskTemp.trials{n}(1:minCountOd,:);
%            taskTemp.muT{n} = taskTemp.muT{n}(1:minCountOd,:);
%            taskTemp.betaT{n} = taskTemp.betaT{n}(1:minCountOd,:);
%        end
       S.(subj{i}).(task{k}) = taskTemp; 
    
    %Adding in the Baseline data to the S structure, there are no latencies
    %or eventTypes for Baseline Data
        content = [num2str(subjectTag(subjectNo)) '_' conditionTag{5}];
        Matfile = fullfile(folder, ['fullrun' content  '.mat']);
        load(Matfile);
        S.(subj{i}).(task{5}) = megDat.data;
        
        
    end
end


%Could've also done EEG = pop_loadbv('fullrun5_I2.vhdr')...
%EEG.event.latency for time labels; EEG.event.type for trial type
%
%allfiles = dir(fullfile(folder2,'*.vhdr'))
% for ii = 1:length(allfiles)
% EEG = pop_loadbv(folder2,allfiles(ii).name);
% end


%% Boxplot with stats test (COMPLETE)

%SPECIFY SUBJECT NUMBER & TASK
%subjectTag = [3 4 5];
subject = ['S' num2str(subjectTag(3))];
% task -> 1:'Close';  2: 'OpenClose'; 3: 'imgClose'; 4: 'imgOpenClose'; 5: 'Base';
lookat = task{3};


restM = [];
for o = 1:2:80
restM = [restM; mean(S.(subject).(lookat).muT{o})];
end

TaskM = [];
for j = 2:2:80
TaskM = [TaskM; mean(S.(subject).(lookat).muT{j}(500:end))];
end
restB = [];
for o = 1:2:80
restB = [restB; mean(S.(subject).(lookat).betaT{o})];
end

TaskB = [];
for j = 2:2:80
TaskB = [TaskB; mean(S.(subject).(lookat).betaT{j}(500:end))];
end

figure; boxplot([restM(:,1),  TaskM(:,1)])

figure; boxplot([restB(:,1),  TaskB(:,1)])

%Kruskal-Wallis DOF = 1, Chi^2 = 20.54, Prob = .000005837 SIGNIFICANT
x = [restM(:,1),  TaskM(:,1)];
kru = kruskalwallis(x);
%ANOVA Prob>F = .00000179 SIGNIFICANT
ano = anova1(x);
%% Average PSD Plots

%SPECIFY SUBJECT # & TASK
%subjectTag = [3 4 5];
subject = ['S' num2str(subjectTag(3))];
% task -> 1:'Close';  2: 'OpenClose'; 3: 'imgClose'; 4: 'imgOpenClose'; 5: 'Base';
lookat = task{3};

Fs = 10000;
nfft = 512*16;
nw = 4;

pxx = []; fxx = [];
sumRest = 0;
for j = 2:2:78 
[pxx{j} fxx{j}] = pmtm(detrend(S.(subject).(lookat).trials{j}), 4, nfft, 1000);
sumRest = [sumRest+pxx{j}];
end
AvgRest = sumRest./39;

sumTask = 0;
for o = 1:2:80
    [pxx{o} fxx{o}] = pmtm(detrend(S.(subject).(lookat).trials{o}), 4, nfft, 1000);
    sumTask = [sumTask+pxx{o}];
end
AvgTask = sumTask./40; 
figure; p1 = plot(fxx{1}, db(AvgTask(:,1)), 'r'); hold on; plot(fxx{2}, db(AvgRest(:,1)), 'b');

xlim([0 100]); xlabel("Frequency (Hz)"); ylabel("Power (db)")





