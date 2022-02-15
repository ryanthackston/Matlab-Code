%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all, clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tic Toc method
frequency = 10;
fs = 10e3;
duration = 20;
samples = 1:1:duration*fs;
t = (0:length(samples)-1)*1/fs;
signal =sin(2*pi*frequency*t);

% Refresh the daq toolbox
daqreset

% To discover a the NI device and save as variable 'devices'
devices = daq.getDevices;
vend = daq.getVendors;

% find what slot # to use in the cDAQ chasis
slot = devices.SlotNumber;
% Create a session
s = daq.createSession('ni')

% Add Counter output channel
addCounterOutputChannel(s, 'cDAQ2Mod1', 2, 'PulseGeneration')

% Set the Session Rate
% ~3.6 Million scans/sec. Need 2.55 Million scans/sec to keep up with 
% Matlabs tic toc function, but it's not exact. This is roughly 1 scan every 392 ns
s.Channels(1).Terminal;
s.Rate = 3600000;
s.NumberOfScans = s.Rate*duration;
s.Channels(1).Terminal;
ch = s.Channels(1);
ch.IdleState='Low';
ch.Frequency = frequency;
ch.InitialDelay = 0.0008;
ch.DutyCycle = 0.5;
s.DurationInSeconds = duration;

% Create checkerboard image
m = 100; n = 160;
[C, K] = checkerb(m, n);
h = imshow(C, 'Border', 'tight', 'InitialMagnification', 'fit');

% f - Length of time for image to be displayed before a blink
f = 1/(frequency*2);
trigger = zeros(8000000, 1);
% k is the number of periods/cycles for fully blinking on then off
k=1;
% n is the index of cycles in the tic toc function
q = 1;
startBackground(s);
while k < (frequency * duration)
        tic
        % d - pause time
        h.CData = K; 
        drawnow('expose')
        while(toc < f)
            trigger(q) = 0.05;
%             queueOutputData(s, trigger(q));
            q = q+1;
        end
        
        tic
        h.CData = C; 
        drawnow('expose')
        while(toc < f)
            trigger(q) = 3.2;
%             queueOutputData(s, trigger(q));
            q = q+1;
        end
        k = k+1;
end

%Delete extra trigger signals
del_trig = find(trigger == 0);
trigger(del_trig(1:end)) = [];

z = trigger;
zs = z(1:(length(z)/(1000*duration)):length(z));
thresh = mean(zs);
tictoc_change = zeros(length(zs),1);
 for c = 1:length(zs)
     if c==1
        tictoc_change(c) = 0;
        continue;
     elseif zs(c-1) == zs(c)
        tictoc_change(c-1) = 0;
     elseif zs(c-1) ~= zs(c)
        tictoc_change(c-1) = 1;
     elseif c == length(zs)
        tictoc_change(c) = 0;
     end
 end

tictoc_times = find(tictoc_change == 1);
for y = 2:length(tictoc_times)
    if y == 1
        tictoc_times(y,2) = 0;
    else
        tictoc_times(y,2) = tictoc_times(y,1) - tictoc_times(y-1,1);
    end
end

clearvars z signal samples rectang










