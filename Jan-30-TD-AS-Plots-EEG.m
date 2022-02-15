%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('C:\Users\Ryan\Downloads\MATLAB Functions\AS_Trials\Feb6_AS_EEGtrials_Scaled1000.mat')
data = {AS750.F7o5.T1_Frame2_Epoch; AS750.F7o5.T2_Frame2_Epoch; AS750.F7o5.T3_Frame2_Epoch; ...
        AS750.F7o5.T4_Frame2_Epoch; AS750.F7o5.T5_Frame2_Epoch; AS750.F7o5.T6_Frame2_Epoch; ...
        AS750.F15.T1_Frame1_Epoch; AS750.F15.T2_Frame1_Epoch; AS750.F15.T3_Frame1_Epoch; ...
        AS750.F15.T4_Frame1_Epoch; AS750.F15.T5_Frame1_Epoch; AS750.F15.T6_Frame1_Epoch; 
        AS750.F1o9.T1_Frame2_Avg; AS750.F1o9.T2_Frame2_Avg; ...
        AS750.F10.T1_Frame2_Avg; AS750.F10.T1_Frame2_Avg; ...
        AS750.F10.T1_Frame1_Avg; AS750.F10.T1_Frame1_Avg};

x = AS.F7.T1 - min(AS.F7.T1);
repx = repmat(x, 100, 1);
[pxx7, fxx7] = pmtm(repx, [], length(repx)*4, 3200);
figure; plot(fxx7, db(pxx7))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = AS.F3o4.T1 - min(AS.F3o4.T1);
repx = repmat(x, 100, 1);
[pxx7, fxx7] = pmtm(repx, 4, length(repx)*4, 3200);
figure; plot(fxx7, db(pxx7))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = AS.F1o7.T1 - min(AS.F1o7.T1);
repx = repmat(x, 100, 1);
[pxx1o7, fxx1o7] = pmtm(repx, 4, length(repx)*4, 3200);
figure; plot(fxx1o7, db(pxx1o7))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = [AS.F7.T1; AS.F7.T2; AS.F7.T3];
[pxx fxx] = pwelch(repx, hann(100), 95, 512*512, 3200);
[pxxn fxxn] = pwelch(noise_x, hann(100), 95, 512*512, 3200);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPECTROGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = [data{10:11}];
x = x(:);
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





x = (x - mean(x))/std(x);

f = 1600*(0:(length(x)))/length(x)

P2 = abs(fft(x,2*length(x))/length(x));
P1 = P2(1:(length(x)+1));
P1(2:end-1) = 2*P1(2:end-1);
figure; plot(f, db(P1)); hold on;

% 1.7 Hz Trials
% z = [AS.F1o7.T1; AS.F1o7.T2; AS.F1o7.T3]
x2 = [data{17:18}];
x2 = x2(:);
z = (x2 - mean(x2))/std(x2);
fz = 1600*(0:(length(z)))/length(z)
Pz2 = abs(fft(z, 2*length(z))/length(z));
Pz1 = Pz2(1:(length(z)+1));
Pz1(2:end-1) = 2*Pz1(2:end-1);

plot(fz, db(Pz1));
leg = legend ('10Hz 2-Frame', '10 Hz 1-Frame');
% xlim([0 30]);

xlabel('Frequency (Hz)');
ylabel('Spectral Amplitude (dB)');
title('Spectral Amplitude of AS: Comparing 1-Frame to 2-Frame with Trigger at 10 Hz');

set(gca,'fontsize', 15)


 [hp2, hp1] = butter(3, 0.5/(Fs/2), 'high');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME DOMAIN PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% d1 = AS750.F7o5.T1_Frame2_Epoch;
% d2 = AS750.F7o5.T2_Frame2_Epoch;
% d3 = AS750.F7o5.T3_Frame2_Epoch;
% d4 = AS750.F7o5.T4_Frame2_Epoch;
% d5 = AS750.F7o5.T5_Frame2_Epoch;
% d6 = AS750.F7o5.T6_Frame2_Epoch;

repx = repmat(x, 100, 1);
t = 0:((0.75)/(length(d1)-1)):0.75;
figure; plot(t, data{13}); hold on;
plot(t, data{14});
% plot(t, d3);
% plot(t, d4);
% plot(t, d5);
% plot(t, d6);

% plot([0.3;0.3], [ylim]')
% plot([0.3*2;0.3*2], [ylim]')
% plot([0.3*3;0.3*3], [ylim]')
% plot([0.3*4;0.3*4], [ylim]')
xlim([0 0.75]);
xlabel(['Time (Sec)']);
ylabel('Amplitude (\muV)');
title('Time Domain - AS 1.7Hz - All Trials');
leg = legend('Trial 1','Trial 2','Trial 3');
set(gca,'fontsize', 15)

% title('Time Domain - AS 7Hz - Edge Artifacts From Copying 1 Trial');

% 961 - 943
p = [0.2944:(1/9600):0.3];
vq1 = interp1(t, repx, p, 'spline')
figure; plot(t, repx, [], p, vq1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(20); title('Spectral Amplitude of AS - Zero-Padded');

[c, lags] = xcorr(data{13}, data{14});
figure; stem(lags, c)

set(gca, 'XTick',[-0.75 -0.6 -0.45 -0.30 -0.15 0 0.15 0.30 0.45 0.60 0.75])
set(gca, 'XTickLabel', str2mat(num2str([-0.75 -0.6 -0.45 -0.30 -0.15 0 0.15 0.30 0.45 0.60 0.75])))

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



































