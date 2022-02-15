
for i = 1:4
    load(['C:\Users\Ryan\Downloads\matlabplotscodechristmasbreak2019\Christin_EEG_Vector' num2str(i)])
end
    load(['C:\Users\Ryan\Downloads\matlabplotscodechristmasbreak2019\ChristinEEGvaluesAll'])


x =[CeegMat1; CeegMat2; CeegMat3; CeegMat4];
Fs = 100;
figure; subplot(2, 2, 1:2); 
t = (0:(length(x)-1))/3200;
p1 = plot(t, x);
title(sprintf('Time Domain Signal'));
hold on;
xlabel('Time (Seconds)');
ylabel('Amplitude (\muV)');

plot([0.3;0.3], [ylim]')
plot([0.3*2;0.3*2], [ylim]')
plot([0.3*3;0.3*3], [ylim]')
plot([0.3*4;0.3*4], [ylim]')


% plot([9.6;9.6], [ylim]')
% plot([9.6*2;9.6*2], [ylim]')
% plot([9.6*3;9.6*3], [ylim]')
% plot([9.6*4;9.6*4], [ylim]')
% plot([9.6*5;9.6*5], [ylim]')
% plot([9.6*6;9.6*6], [ylim]')
% plot([9.6*7;9.6*7], [ylim]')



[pxx fxx] = pwelch(x, 100, 95, 512*512, 3200);
subplot(2, 2, 3); 
p2 = plot(fxx, db(pxx/2));
title(sprintf('PSD of Checkerboard Trials'));
xlim([0 40]);
hold on;
xlabel(['Frequency (Hz)']);
ylabel(['Power (dB/\surd(Hz))']);

[pxx2 fxx2] = pmtm(x, 4, 512*512, 3200);
plot(fxx2, db(pxx2/2));
plot([1.7;1.7], [ylim]')
lgd = legend({'Welch PSD', 'Multitaper PSD'}, 'FontSize', 10)

subplot(2,2,4);
t2 = (0 : length(CeegMat1)-1)/3200
p3 = plot(t2, [CeegMat1 CeegMat2 CeegMat3 CeegMat4]);
hold on;
p4= plot(t2, CeegMatAvgCheck);
p4.LineWidth = 2;
title('Grand Average Time Domain of Checkerboard Trials - Christin');
xlabel('Time (Seconds)');
ylabel('Amplitude');

suptitle('EEG of Oz - Checkerboard Pattern - Christin')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Blinking LED - Christin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 5:8
    load(['C:\Users\Ryan\Downloads\matlabplotscodechristmasbreak2019\Christin_EEG_Vector' num2str(i)])
end
      

x =[CeegMat5; CeegMat6; CeegMat7; CeegMat8];
Fs = 100;
figure; subplot(2, 2, 1:2); 
t = (0:(length(x)-1))/3200;
p1 = plot(t, x);
title(sprintf('Time Domain Signal'));
hold on;
xlabel('Time (Seconds)');
ylabel('Amplitude (\muV)');

plot([0.5;0.5], [ylim]')
plot([0.5*2;0.5*2], [ylim]')
plot([0.5*3;0.5*3], [ylim]')
plot([0.5*4;0.5*4], [ylim]')


% plot([9.6;9.6], [ylim]')
% plot([9.6*2;9.6*2], [ylim]')
% plot([9.6*3;9.6*3], [ylim]')
% plot([9.6*4;9.6*4], [ylim]')
% plot([9.6*5;9.6*5], [ylim]')
% plot([9.6*6;9.6*6], [ylim]')
% plot([9.6*7;9.6*7], [ylim]')



[pxx fxx] = pwelch(x, 100, 95, 512*512, 1920);
subplot(2, 2, 3); 
p2 = plot(fxx, db(pxx/2));
title(sprintf('PSD of Blinking LED Trials'));
xlim([0 40]);
hold on;
xlabel(['Frequency (Hz)']);
ylabel(['Power (dB/\surd(Hz))']);

[pxx2 fxx2] = pmtm(x, 4, 512*512, 1920);
plot(fxx2, db(pxx2/2));
plot([2.1;2.1], [ylim]')
lgd = legend({'Welch PSD', 'Multitaper PSD'}, 'FontSize', 10)

subplot(2,2,4);
t2 = (0 : length(CeegMat1)-1)/1920
p3 = plot(t2, [CeegMat5 CeegMat6 CeegMat7 CeegMat8]);
hold on;
p4= plot(t2, CeegMatAvgLED);
p4.LineWidth = 2;
title('Grand Average Time Domain of Blinking LED Trials - Christin');
xlabel('Time (Seconds)');
ylabel('Amplitude');

suptitle('EEG of Oz - Blinking LED - Christin')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checkerboard Pattern - Ryan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:4
    load(['C:\Users\Ryan\Downloads\matlabplotscodechristmasbreak2019\Ryan_EEG_Vector' num2str(i)])
end
      
load(['C:\Users\Ryan\Downloads\matlabplotscodechristmasbreak2019\RyanEEGvaluesAll'])

x =[ReegMat1; ReegMat2; ReegMat3; ReegMat4];
Fs = 100;
figure; subplot(2, 2, 1:2); 
t = (0:(length(x)-1))/3200;
p1 = plot(t, x);
title(sprintf('Time Domain Signal'));
hold on;
xlabel('Time (Seconds)');
ylabel('Amplitude (\muV)');

plot([0.3;0.3], [ylim]')
plot([0.3*2;0.3*2], [ylim]')
plot([0.3*3;0.3*3], [ylim]')
plot([0.3*4;0.3*4], [ylim]')


% plot([9.6;9.6], [ylim]')
% plot([9.6*2;9.6*2], [ylim]')
% plot([9.6*3;9.6*3], [ylim]')
% plot([9.6*4;9.6*4], [ylim]')
% plot([9.6*5;9.6*5], [ylim]')
% plot([9.6*6;9.6*6], [ylim]')
% plot([9.6*7;9.6*7], [ylim]')



[pxx fxx] = pwelch(x, 100, 95, 512*512, 3200);
subplot(2, 2, 3); 
p2 = plot(fxx, db(pxx/2));
title(sprintf('PSD of Checkerboard Trials'));
xlim([0 40]);
hold on;
xlabel(['Frequency (Hz)']);
ylabel(['Power (dB/\surd(Hz))']);

[pxx2 fxx2] = pmtm(x, 4, 512*512, 3200);
plot(fxx2, db(pxx2/2));
plot([1.7;1.7], [ylim]')
lgd = legend({'Welch PSD', 'Multitaper PSD'}, 'FontSize', 10)

subplot(2,2,4);
t2 = (0 : length(ReegMat1)-1)/3200
p3 = plot(t2, [ReegMat1 ReegMat2 ReegMat3 ReegMat4]);
hold on;
p4= plot(t2, ReegMatAvgCheck);
p4.LineWidth = 2;
title('Grand Average Time Domain of Checkerboard Trials - Ryan');
xlabel('Time (Seconds)');
ylabel('Amplitude');

suptitle('EEG of Oz - Checkerboard Pattern - Ryan')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Blinking LED - Ryan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 5:8
    load(['C:\Users\Ryan\Downloads\matlabplotscodechristmasbreak2019\Ryan_EEG_Vector' num2str(i)])
end
      
% 
% 
x = AS.F7.T1 - min(AS.F7.T1);

Fs = 3200;
figure; subplot(2, 2, 1:2); 
t = (0:(length(x)-1))/3200;
p1 = plot(t, x);
title(sprintf('Time Domain Signal'));
hold on;
xlabel('Time (Seconds)');
ylabel('Amplitude (\muV)');

plot([0.5;0.5], [ylim]')
plot([0.5*2;0.5*2], [ylim]')
plot([0.5*3;0.5*3], [ylim]')
plot([0.5*4;0.5*4], [ylim]')


% plot([9.6;9.6], [ylim]')
% plot([9.6*2;9.6*2], [ylim]')
% plot([9.6*3;9.6*3], [ylim]')
% plot([9.6*4;9.6*4], [ylim]')
% plot([9.6*5;9.6*5], [ylim]')
% plot([9.6*6;9.6*6], [ylim]')
% plot([9.6*7;9.6*7], [ylim]')



[pxx fxx] = pwelch(x, 100, 95, 512*512, 1920);
subplot(2, 2, 3); 
p2 = plot(fxx, db(pxx/2));
title(sprintf('PSD of Blinking LED Trials'));
xlim([0 40]);
hold on;
xlabel(['Frequency (Hz)']);
ylabel(['Power (dB/\surd(Hz))']);

[pxx2 fxx2] = pmtm(x, 4, 512*512, 1920);
plot(fxx2, db(pxx2/2));
plot([2.1;2.1], [ylim]')
lgd = legend({'Welch PSD', 'Multitaper PSD'}, 'FontSize', 10)

subplot(2,2,4);
t2 = (0 : length(ReegMat1)-1)/1920
p3 = plot(t2, [ReegMat5 ReegMat6 ReegMat7 ReegMat8]);
hold on;
p4= plot(t2, ReegMatAvgLED);
p4.LineWidth = 2;
title('Grand Average Time Domain of Blinking LED Trials - Ryan');
xlabel('Time (Seconds)');
ylabel('Amplitude');
leg2 = legend({'Trial 5', 'Trial 6', 'Trial 7', 'Trial 8', 'Grand Average'})

suptitle('EEG of Oz - Blinking LED - Ryan')
