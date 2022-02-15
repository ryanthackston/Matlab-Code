%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function BlinkLoop(C, K)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
Frequency = 10;
Seconds = 10;
SamplingFreq = 1000;
trigger = zeros (Seconds*SamplingFreq,1);
pauseTime = (1/(Frequency*2));

trigger = zeros (Seconds*SamplingFreq,1);
for i = 1:(length(trigger)/Frequency)
    if mod(i,2) == 1
        trigger(i*Frequency) = 1;
    elseif mod(i,2) == 0
        trigger(i*Frequency) = 2;
    end
end


img1 = imshow(C, 'Border', 'tight', 'InitialMagnification', 'fit'), ...
h1 = getframe;
img2 = imshow(K,  'Border', 'tight', 'InitialMagnification', 'fit'), ...
h2 = getframe;

hFig = figure('Name','APP',...
    'Numbertitle','off',...
    'Position', [0 0 1680 950],...
    'WindowStyle','normal',...
    'Color',[0.5 0.5 0.5],...
    'Toolbar','none');

%    MATLAB:concatenation:integerInteraction
warning('off', 'Images:initSize:adjustingMag');
hold all;
% Every loop is 1/Frequency seconds, i.e freq = 10, loop is .1
for z = 1:length(trigger)
    if trigger(z) == 0
       continue;
    elseif trigger(z) == 1
       imshow(h1.cdata, 'Border', 'tight', 'InitialMagnification', 'fit'),
    elseif trigger(z)== 2
       imshow(h2.cdata, 'Border', 'tight', 'InitialMagnification', 'fit'),
    end
end


th = linspace(1, 1000, length(trigger)/Frequency);
tic
if trigger(th*Frequency) == 1
   imshow(h1.cdata, 'Border', 'tight', 'InitialMagnification', 'fit'), pause(1/(Frequency*2));
elseif trigger(th*Frequency)== 2
   imshow(h2.cdata, 'Border', 'tight', 'InitialMagnification', 'fit'), pause(1/(Frequency*2));
end
toc
% One full loop is 1 period. 

