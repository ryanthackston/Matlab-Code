%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MEG Measurement
%% 12.11.2019 (ee)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate measurement name
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = datestr(now,'HHMMSS');
filename = ['meg_measurement_' dt];
% addpath('data')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T.sig = 20;               % Duration of recoding in seconds
f.s   = 10e3;               % Sampling frequency
N.ch  = 4;                  % Number of channels
conv.v_b = (1/2.7)*1e-9;    % Conversation factor
frequency = 3.4;
duration = T.sig;
samples = 1:1:duration*f.s;
% t = (0:length(samples)-1)*1/f.s;
% signal =sin(2*pi*frequency*t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Derive further paramters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N.sig = T.sig * f.s;        % Number of samples of the recording
n.sig = 0:(N.sig-1);        % Sample vector
t.sig = n.sig / f.s;        % Time vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hardware setup (NI-Card) - Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hw.session                      = daq.createSession('ni');
hw.session.Rate                 = f.s;
hw.session.DurationInSeconds    = T.sig;
hw.name                         = 'cDAQ2Mod'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Measurement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ii=0;  
% for iCh = 1:N.ch 
%     if mod(iCh,4) == 1
%         ii=ii+1;
%     end 
%     ch          = addAnalogInputChannel(hw.session,[hw.name num2str(ii)] , ['ai' num2str(mod(iCh-1,4))], 'Voltage');
%     ch.Range    = [-1 1]*10;
%     ch.Coupling = 'DC';
% end

ii = ii+1;
chCou = addCounterOutputChannel(hw.session, 'cDAQ2Mod2', [0 1], 'PulseGeneration');
chCou(1).IdleState='Low';

ii = ii+1;
chCou(1).Frequency = 3.4;
chAn = addAnalogInputChannel(s,'cDAQ2Mod3', [0 1], 'Voltage');

[data,time] = s.startForeground();

figure; plot(time,data);
xlabel('Time (secs)');
ylabel('Current')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Measurement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Start measurement');
x.in = hw.session.startForeground; 

disp('Finished measurement');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Convert data to B-field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iCh=1:N.ch
   x.b(:,iCh) = x.in(:,iCh)*conv.v_b; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add naming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chname = {'Sensor DV - Z','Sensor DY - Z','Sensor DW - Z','Sensor DX - Z'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save measurement to file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save(['data\dat_' filename],'x','f','N','T','n','t','ch','chname');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analysis of measurement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
linew = 1.5;
colors = parula(N.ch);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sensornames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
measurementname = 'Noise-Measurement of QZFM-OPMs';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ff=0;
ff=ff+1; names{ff} ='QZFM-00DV';
ff=ff+1; names{ff} ='QZFM-00DY';
ff=ff+1; names{ff} ='QZFM-00DW';
ff=ff+1; names{ff} ='QZFM-00DX';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot RAW input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('units','normalized','outerposition',[0 0 1 1])
for cc = 1:4
    subplot(4,1,cc)
    plot(t.sig,x.in(:,cc)-mean(x.in(:,cc)),...
        'Color',colors(cc,:),'LineWidth',linew);
    grid on; xlabel('Time / s'); ylabel('Amplitude / V');
    title(names{cc})
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Some preprocessing / detrending
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [b.lp a.lp] = butter(4,100/(f.s/2));
% [b.hp a.hp] = butter(4,10/(f.s/2),'high');
%  x.in = filtfilt(b.hp,a.hp,x.in);
% [x.in(:,1:4),~] = detrending(x.in,f.s);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make conversition to B-field (2.7 V/nT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
conv.t_v = 1/2.7*1e-9; 
for cc = 1:4
    x.Bin(:,cc) = conv.t_v.*x.in(:,cc);    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot input after conversion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('units','normalized','outerposition',[0 0 1 1])
%suptitle(measurementname)
for cc = 1:4
    subplot(4,1,cc)
    plot(t.sig,(x.Bin(:,cc)-mean(x.Bin(:,cc)))*1e12,'Color',...
        colors(cc,:),'LineWidth',linew);
    grid on; xlabel('Time / s'); ylabel('Amplitude / pT');
    title(names{cc})
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate PSD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N.fft = f.s;
ana.window = flattopwin(N.fft);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for cc=1:N.ch
    [P.flux(:,cc) f.pvec] = pwelch(x.Bin(:,cc)-mean(x.Bin(:,cc)),...
        ana.window,[],N.fft,f.s);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot ASD-Overview
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('units','normalized','outerposition',[0 0 1 1])
for cc=1:N.ch
    semilogy(f.pvec,sqrt(P.flux(:,cc)),'Color',colors(cc,:),...
        'DisplayName',names{cc},'LineWidth',linew); 
    hold on;
end
grid on;
set(gca,'XScale','log');
legend('show');
xlabel('Frequency / Hz','Interpreter','latex','FontSize',14);
ylabel('Amplitude density $$T/ \sqrt{Hz}$$','Interpreter','latex','FontSize',14)
legend('show','Location','Best');
title(measurementname)
xlim([1e0 1e3]);
savefig(['data\fig_' filename '.fig'])
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Load old Measurement for comparison
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% uiopen('C:\Users\demo\Desktop\MSR\opm_noise_measurement_psd.fig',1)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%