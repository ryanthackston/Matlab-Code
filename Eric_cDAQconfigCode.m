%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test Config
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -global
daqreset
% To discover a the NI device and save as variable 'devices'
d = daq.getDevices;
vend = daq.getVendors;
% find what slot # to use in the cDAQ chasis
slot = d.SlotNumber;

global data
s = daq.createSession('ni');
s.addAnalogInputChannel('cDAQ2Mod2', 0,'Voltage')
s.Rate = 10000
s.DurationInSeconds = 10
lh = s.addlistener('DataAvailable',@plotData);
s.startBackground();
% Do something
    s.IsDone % will report 0
    s.wait(10) % rather than while
    s.IsDone % will report 1
    
close(gcf);
plot(data); %  plot global data

plotData(src,event)

% ch.IdleState='Low' Keeps the voltage on the low end (... 1V) when no
% there are no inputs. ch.IdleState='High' would keep it on the high end
% (... 3.5V).

% Add Counter output channel
% chCou = addCounterOutputChannel(s, 'cDAQ2Mod1', [0 1], 'PulseGeneration');
% chCou(1).IdleState='Low';
% 
% chCou(1).Frequency = 3.4;

% chAn = addAnalogInputChannel(s,'cDAQ2Mod2', [0 1], 'Voltage');
% lh = addlistener(s,'DataAvailable', @(src,event) plot(event.TimeStamps, event.Data));

% acquireData
% % outputLow = 1.5;
% % outputHigh = 2.5;
% % outputSingleScan(s,[outputLow outputHigh]);
% 
% s.startBackground();
% 
% xlabel('Time (secs)');
% ylabel('Current')
% 
% % Set the Session Rate
% % ~3.6 Million scans/sec. Need 2.55 Million scans/sec to keep up with 
% % Matlabs tic toc function, but it's not exact. This is roughly 1 scan every 392 ns
% % s.Channels(1).Terminal;
% % s.Rate = 3600000;
% % s.NumberOfScans = s.Rate*duration;
% % s.Channels(1).Terminal;
% % ch = s.Channels(1);
% % ch.IdleState='Low';
% % ch.Frequency = frequency;
% % ch.InitialDelay = 0.0008;
% % ch.DutyCycle = 0.5;
% % s.DurationInSeconds = 10;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % startBackground(s);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%