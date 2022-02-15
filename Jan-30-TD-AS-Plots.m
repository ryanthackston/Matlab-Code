
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% alldata takes about 30 seconds to load.
load('C:\Users\Ryan\Desktop\Kiel MEG Data\Data\alldata.mat')

% alldata -> Dates
% alldata.nov12_opm_noise -> Trial -> x -> b
% alldata.nov21_opm -> Trial type, Stim Freq, Subject -> x -> b (Magnetic field data, 5th column is trigger)
% alldata.nov28_opm -> Trial type, Stim Freq, Subject -> xin -> b (Magnetic field data)
% alldata.dec16_eeg -> Subject -> Stim Type (Checkerboard 1.7Hz, BlinkLED EEG Matrix - Each column is 1 trial
% alldata.jan30_eeg -> Subject -> Stim Frequency -> Trial # 
% alldata.feb6_eeg -> Subject -> Stim Frequency -> Trial # and type


% Feb-6 data
% % Subject code - AS; Epoch size - 750ms;
% data = {AS750.F7o5.T1_Frame2_Epoch; AS750.F7o5.T2_Frame2_Epoch; AS750.F7o5.T3_Frame2_Epoch; ...
%         AS750.F7o5.T4_Frame2_Epoch; AS750.F7o5.T5_Frame2_Epoch; AS750.F7o5.T6_Frame2_Epoch; ...
%         AS750.F15.T1_Frame1_Epoch; AS750.F15.T2_Frame1_Epoch; AS750.F15.T3_Frame1_Epoch; ...
%         AS750.F15.T4_Frame1_Epoch; AS750.F15.T5_Frame1_Epoch; AS750.F15.T6_Frame1_Epoch; 
%         AS750.F1o9.T1_Frame2_Avg; AS750.F1o9.T2_Frame2_Avg; ...
%         AS750.F10.T1_Frame2_Avg; AS750.F10.T2_Frame2_Avg; ...
%         AS750.F10.T1_Frame1_Avg; AS750.F10.T2_Frame1_Avg};
data = [];
data = [data struct2cell(alldata.feb6_eeg.AS750.F7o5);...
        struct2cell(alldata.feb6_eeg.AS750.F15);...
        struct2cell(alldata.feb6_eeg.AS750.F1o9);...
        struct2cell(alldata.feb6_eeg.AS750.F10)]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multitaper Method

x = reshape([data{1:6}], [], 1);

% Pmtm(data input, 4th order(usually recommended), NFFT, Fs)
[pxx7, fxx7] = pmtm(x, 4, length(x)*4, 960/0.75);
figure; plot(fxx7, db(pxx7))
ylabel('Power (dB/Hz)');
xlabel('Frequency (Hz)');
title('Multitaper - 7.5Hz Epochs Concatenated');
set(gca,'fontsize', 18)


% Feb 6 data - Averaging Epoch data and doing Multitaper
x1 = [data{1:6}];
x1AVG = mean(x1, 2);
x2 = [data{7:12}];
x2AVG = mean(x2, 2);

[pxx1, fxx1] = pmtm(x1AVG, 4, length(x)*128, 960/0.75);
[pxx2, fxx2] = pmtm(x2AVG, 4, length(x)*128, 960/0.75);

% Multitaper Averaged Data vs Epoched Data
figure; subplot(2,1,1); plot(fxx1, db(pxx1));
hold on; plot(fxx2, db(pxx2));
ylabel('Power (dB/Hz)');
xlabel('Frequency (Hz)');
title('Multitaper - 7.5Hz 2-Frame Avg vs 15 Hz 1-Frame Avg');
set(gca,'fontsize', 21)
leg = legend('7.5 Hz 2-Frame', '15 Hz 1-Frame');
xlim([0 40]);
ylim([-40 40]);

[pxx3, fxx3] = pmtm(x1(:), 4, length(x)*128, 960/0.75);
[pxx4, fxx4] = pmtm(x2(:), 4, length(x)*128, 960/0.75);
subplot(2,1,2); plot(fxx3, db(pxx3));

hold on; plot(fxx4, db(pxx4));
ylabel('Power (dB/Hz)');
xlabel('Frequency (Hz)');
title('Multitaper - 7.5Hz 2-Frame vs 15 Hz 1-Frame Concatenated');
set(gca,'fontsize', 21)
leg = legend('7.5 Hz 2-Frame', '15 Hz 1-Frame');
xlim([0 40]);
ylim([-40 40]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Welch's Method
x = [data{7:12}];
x = x(:);
noise_x = x + (50/70)*randn(length(x),1);
[pxx fxx] = pwelch(x, hann(1920), 1910, 512*512, 960/0.75);
[pxxn fxxn] = pwelch(noise_x, hann(960), 950, 512*512, 3200);

figure; plot(fxx, db(pxx)); hold on;
plot(fxx, db(pxxn));
legend('Trials Concatenated', 'Trials Concatenated + Noise')
ylabel('Power (dB/Hz)');
xlabel('Frequency (Hz)');
title('Pwelch - 7.5Hz Epochs 1-6 - Trials vs Trials with OPM Noise');
set(gca,'fontsize', 18)

% Feb 6 Pwelch - Averaged Data vs Epoched Data
x1 = [data{1:6}];
x1AVG = mean(x1, 2);
x2 = [data{7:12}];
x2AVG = mean(x2, 2);

[pxx1, fxx1] = pwelch(x1AVG, hann(480), 470, 512*512, 960/0.75);
[pxx2, fxx2] = pwelch(x2AVG, hann(480), 470, 512*512, 960/0.75);
figure; subplot(2,1,1);
plot(fxx1, db(pxx1));
hold on; plot(fxx2, db(pxx2));
ylabel('Power (dB/Hz)');
xlabel('Frequency (Hz)');
title('Pwelch 7.5Hz 2-Frame Avg vs 15 Hz 1-Frame Avg');
set(gca,'fontsize', 21)
xlim([0 40]);
ylim([-40 40]);

leg = legend('7.5 Hz - 2-Frame', '15 Hz 1-Frame');

[pxx3, fxx3] = pwelch(x1(:), hann(1920), 1910, 512*512, 960/0.75);
[pxx4, fxx4] = pwelch(x2(:), hann(1920), 1910, 512*512, 960/0.75);

subplot(2,1,2);
plot(fxx3, db(pxx3));
hold on; plot(fxx4, db(pxx4));
ylabel('Power (dB/Hz)');
xlabel('Frequency (Hz)');
title('Pwelch - 7.5Hz 2-Frame Avg vs 15 Hz 1-Frame Concatenated (Larger Window - More Data)');
set(gca,'fontsize', 21)
xlim([0 40]);
ylim([-40 40]);

leg = legend('7.5 Hz 2-Frame', '15 Hz 1-Frame');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPECTROGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% x = [data{10:11}];
% x = x(:);

x = [x1AVG; x2AVG]
[s, f, t, p] = spectrogram(x, hamming(960),950, 256*256, (960/0.75), 'yaxis');

% Making the Figure, 
map = flipud(lbmap(100, 'redblue'));
figure; h = pcolor(t, f, db(p)/2); 
set(h, 'edgecolor', 'none'); ylim([0 100]); hcb = colorbar; colormap(map); caxis([-25 11]);
ylabel("Frequency (Hz)"); xlabel('Time (sec)');
colorTitleHandle = get(hcb,'Title');
titleString = 'Power (dB/Hz)';
set(colorTitleHandle ,'String',titleString);

% title('Spectrogram - AS 15Hz 1-Frame 750ms Epoch - Trial 7-12')
title('Spectrogram - AS 15Hz 1-Frame 750ms Epoch - Trial 10-11')
set(gca,'fontsize', 18)
 
%  figure; spectrogram(x, hamming(750), 300, 256*256, 1280, 'yaxis')
%  title('Spectrogram - AS - 7Hz - All Trials - 50 wind - 25 over - 256*256 nfft');

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPECTRAL AMPLITUDE PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 7.5 Hz Trials Epoch
% x = [AS750.F7o5.T1_Frame2_Epoch; AS750.F7o5.T2_Frame2_Epoch; AS750.F7o5.T3_Frame2_Epoch;...
%     AS750.F7o5.T4_Frame2_Epoch; AS750.F7o5.T5_Frame2_Epoch; AS750.F7o5.T6_Frame2_Epoch];

% 15 Hz Trials Epoch
% x = [AS750.F15.T1_Frame1_Epoch; AS750.F15.T2_Frame1_Epoch; AS750.F15.T3_Frame1_Epoch;...
%      AS750.F15.T4_Frame1_Epoch; AS750.F15.T5_Frame1_Epoch; AS750.F15.T6_Frame1_Epoch];

% 10 Hz Frame 1 Avg
% x = [AS750.F10.T1_Frame1_Avg; AS750.F10.T2_Frame1_Avg];


x = [data{15:16}];
x = x(:);

stdx = (x - mean(x))/std(x);

f = 1600*(0:(length(x)))/length(x)

P2 = abs(fft(stdx,2*length(x))/length(x));
P1 = P2(1:(length(stdx)+1));
P1(2:end-1) = 2*P1(2:end-1);

x2 = [data{17:18}];
x2 = x2(:);
stdz = (x2 - mean(x2))/std(x2);
fz = 1600*(0:(length(x2)))/length(x2)
Pz2 = abs(fft(stdz, 2*length(x2))/length(x2));
Pz1 = Pz2(1:(length(x2)+1));
Pz1(2:end-1) = 2*Pz1(2:end-1);

figure; plot(f, db(P1)); hold on;
plot(fz, db(Pz1));
leg = legend ('10Hz 2-Frame', '10 Hz 1-Frame');
% xlim([0 30]);

xlabel('Frequency (Hz)');
ylabel('Spectral Amplitude (dB)');
title('Spectral Amplitude of AS: Comparing 1-Frame to 2-Frame with Trigger at 10 Hz');

set(gca,'fontsize', 15)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency Domain SNR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    load('C:\Users\Ryan\Downloads\MATLAB Functions\AS_Trials\AS_EEGtrials_Scaled1000')
    
    x = [data{1:6}];
    x = x(:);
%     x = (x-min(x));
    x_n = x + (50/70)*randn(length(x),1);
    
    Fs = length(data{1})/0.75;

    %FFT Power
    jc = abs(fft(x));
    jc_n = abs(fft(x_n));
    jc_n(find(jc == 0)) = [];
    jc(find(jc == 0)) = [];
%     % f = 0 : (Fs/(NFFT/2 - 1)) : 1600
    f = 0:(Fs/(length(jc) - 1)):(Fs);
    df = Fs/(length(jc)*2 - 1);
    % Round up f values
    Ndec = 4;
    rou = 10.^Ndec;
    f = round(rou*f)/rou;
    f = f';
    % Take out every other point that is 0
%     jc([find(jc==0)]) = [];

    snrat = zeros(length(f),1);
    bin = 1.0;
%     snrat(1) = 0;

    % finding the SNR for a specific frequency - length(x) is too low for
    % frequency resolution
    for i = 2:length(f)
        % frequency measured is the index of vector f
        freq = f(i);
        % Make a vector of index values in the frequency bin.
        % Bin*2 gives bandwidth to measure the signal frequency against.
        vec = find(f <= (freq + bin) & f >= (freq - bin));
        if i == 1
            snrat(i) = jc(i) / (sum(jc(vec(find(vec > i)))));
        elseif i == length(f)
            snrat(i) = jc(length(f))/ (sum(jc(vec(find(vec < i)))));
        else
        % SNR equation = (FFT value of specific frequency) / (FFT Values of Vec)
            snrat(i) = jc(find(f==freq))/( (sum(jc_n(vec(find(f==freq) > vec)))) + (sum(jc_n(vec(find(f==freq) < vec)))) );
        end
    end

    snrat(length(f)) = 0;

    figure; plot(f, snrat);
    xlim([0 200]);
    xlabel('Frequency (Hz)');
    ylabel('SNR');
    ylim([-inf inf]);
    title(['SSVEP SNR - AS 15Hz 1-Frame Epoch Trial 7-12']);
    set(gca,'fontsize', 15);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME DOMAIN PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 10 Hz: 2 Frames/Trigger Vs 1 Frame/Trigger
figure; subplot(2,1,1); 

x1 = data{15};
x2 = data{16};

t = 0:((0.75)/(length(x1)-1)):0.75;
plot(t, x1); hold on;
plot(t, x2);

xlim([0 0.75]);
xlabel(['Time (Sec)']);
ylabel('Amplitude (\muV)');
title('Time Domain - AS 10Hz - 2 Frames per Trigger Period');
leg = legend('Trial 1','Trial 2','Trial 3');
set(gca,'fontsize', 15)

subplot(2,1,2); 
x1 = data{17};
x2 = data{18};

t = 0:((0.75)/(length(x1)-1)):0.75;
plot(t, x1); hold on;
plot(t, x2);

xlim([0 0.75]);
xlabel(['Time (Sec)']);
ylabel('Amplitude (\muV)');
title('Time Domain - AS 10Hz - 1 Frame per Trigger Period');
leg = legend('Trial 1','Trial 2','Trial 3');
set(gca,'fontsize', 15)

% 7.5 Hz & 15 Hz 2 Epochs
figure;

x1 = [data{1:6}];
x1AVG = mean(x1, 2);
x2 = [data{7:12}];
x2AVG = mean(x2, 2);

t = 0:((0.75)/(length(x1)-1)):0.75;
plot(t, x1AVG); 
hold on;
plot(t, x2AVG);

xlim([0 0.75]);
xlabel(['Time (Sec)']);
ylabel('Amplitude (\muV)');
title('Time Domain - AS 7.5Hz - 2 Frames per Trigger Period');
leg = legend('Trial 1','Trial 2','Trial 3');
set(gca,'fontsize', 15)



subplot(2,1,2); 
x1 = data{10};
x2 = data{11};

t = 0:((0.75)/(length(x1)-1)):0.75;
plot(t, x1); hold on;
plot(t, x2);

xlim([0 0.75]);
xlabel(['Time (Sec)']);
ylabel('Amplitude (\muV)');
title('Time Domain - AS 15Hz - 1 Frame per Trigger Period');
leg = legend('Trial 1','Trial 2','Trial 3');
set(gca,'fontsize', 15)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Correlation Analysis - ... AS 10 Hz - 1 Frame/Trigger (not finished)
figure(20); title('Spectral Amplitude of AS - Zero-Padded');

[c, lags] = xcorr(data{17}, data{18});
figure; stem(lags, c)
xlim([0 1000])

set(gca, 'XTick',[ 0 0.15 0.30 0.45 0.60 0.75])
set(gca, 'XTickLabel', str2mat(num2str([0 0.15 0.30 0.45 0.60 0.75])))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRIGGER SIGNAL PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Checker board plots
m = 30; n = 48;
% Checkerboard 2 - ON
[C, K] = checkerb(m, n);
h = imshow(C, 'Border', 'tight', 'InitialMagnification', 'fit');
hold on; pause(0.3);
h2 = plot(n/2, m/2, 'r.', 'MarkerSize', 60)
hold off; pause(0.3);
% CHECKERBOARD 1 - OFF
h.CData = K; drawnow('expose')


%Trigger signal at 1Hz - 2 Frame Checkerboard
t1 = linspace(0, 1, 2000000);
% t2 = linspace(1, 2, 10000);
x1 = zeros(1000000*2,1);
x1(1000001:end) = 2.8;

figure; plot(t1, x1);
xlabel('Time (sec)');
ylabel('Amplitude (V)');
title('Trigger Signal 1Hz - 2-Frame');
set(gca,'fontsize', 18)


%Trigger signal at 1Hz - 1 Frame Checkerboard
t1 = linspace(0, 2, 2000000*2);
% t2 = linspace(1, 2, 10000);
x1 = zeros(1000000*2,1);
x1(1000001:2000000) = 2.8;
x1(3000001:4000000) = 2.8;

figure; plot(t1, x1);
xlabel('Time (sec)');
ylabel('Amplitude (V)');
title('Trigger Signal 1Hz - 1 Frame');
set(gca,'fontsize', 18)



































