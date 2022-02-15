%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all, clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m = 30; n = 48;
[C, K] = checkerb(m, n);
h = imshow(C, 'Border', 'tight', 'InitialMagnification', 'fit');

% Tic Toc method
frequency = 1.7;
duration = 20;
fs = 10e3;
samples = 1:1:duration*fs;
t = (0:length(samples)-1)*1/fs;
signal =sin(2*pi*frequency*t);


daqreset
% To discover a the NI device and save as variable 'devices'
d = daq.getDevices;
vend = daq.getVendors;
% find what slot # to use in the cDAQ chasis
slot = d.SlotNumber;

s = daq.createSession('ni');
% Add Counter output channel
ch = addCounterOutputChannel(s, 'cDAQ2Mod1', 1:2, 'PulseGeneration')
chAn = addAnalogInputChannel(s,'cDAQ2Mod2', 1:2, 'Voltage');

lh = addlistener(s,'DataAvailable', @(src,event) plot(event.TimeStamps, event.Data));
% lh = addlistener(src,'EventName',@(src,evnt)obj.callbackMethod(src,evnt,arg1,...argn)
% Set the Session Rate
% ~3.6 Million scans/sec. Need 2.55 Million scans/sec to keep up with 
% Matlabs tic toc function, but it's not exact. This is roughly 1 scan every 392 ns
s.Rate = 10000;
s.NumberOfScans = s.Rate*duration;
ch(1).IdleState='Low';
ch(1).Frequency = frequency;
ch(1).InitialDelay = 0.0008;
ch(2).IdleState='Low';
ch(2).Frequency = frequency;
ch(2).InitialDelay = 0.0008;
s.DurationInSeconds = duration;

% Create checkerboard image
m = 30; n = 48;
[C, K] = checkerb(m, n);
h = imshow(C, 'Border', 'tight', 'InitialMagnification', 'fit');
hold on;
h2 = plot(n/2, m/2, 'r.', 'MarkerSize', 60)
hold off;

% f - Length of time for image to be displayed before a blink
f = 1/(frequency*2);
trigger = zeros(8000000, 1);
% k is the number of periods/cycles for fully blinking on then off
k=1;
% n is the index of cycles in the tic toc function
q = 1;

s.NotifyWhenDataAvailableExceeds = length(samples);
s.startBackground();
while k < (frequency * duration)
        tic
        % d - pause time
        h.CData = K; 
        drawnow('expose')
        while(toc < f)
            trigger(q) = 1;
            q = q+1;
        end
        
        tic
        h.CData = C; 
        drawnow('expose')
        while(toc < f)
            trigger(q) = 2;
            q = q+1;
        end
        k = k+1;
end

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

stop(s);
release(s);
clearvars z signal samples rectang
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
















