function BlinkLoop(C, K)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

figure(1), set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% Get rid of tool bar and pulldown menus that are along top of figure.
set(gcf, 'Toolbar', 'none', 'Menu', 'none')

for SomeLoop = 1:10
    
 figure(1), set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]), ...
 % Get rid of tool bar and pulldown menus that are along top of figure.
 set(gcf, 'Toolbar', 'none', 'Menu', 'none'), ...
 imshow(C, 'Border', 'tight', 'InitialMagnification', 1000), ...
 pause(0.5),
 figure(1), set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]), ...
 % Get rid of tool bar and pulldown menus that are along top of figure.
 set(gcf, 'Toolbar', 'none', 'Menu', 'none'), ...
 imshow(K, 'Border', 'tight'), ...
 pause(0.5)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function BlinkLoop(C, K)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
Frequency = 17;
Seconds = 10;
pauseTime = (1/(Frequency*2));
% figure(1), set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% Get rid of tool bar and pulldown menus that are along top of figure.

img1 = imshow(C, 'Border', 'tight', 'InitialMagnification', 'fit'), ...
h1 = getframe;
img2 = imshow(K,  'Border', 'tight', 'InitialMagnification', 'fit'), ...
h2 = getframe;

hFig = figure('Name','APP',...
    'Numbertitle','off',...
    'Position', [0 0 1680 1050],...
    'WindowStyle','modal',...
    'Color',[0.5 0.5 0.5],...
    'Toolbar','none');
% figure('Position',[0 0 1680 1050]),
%    MATLAB:concatenation:integerInteraction
warning('off', 'Images:initSize:adjustingMag');
hold all;
tic
for SomeLoop = 1:(Frequency*Seconds)
 imshow(h1.cdata, 'Border', 'tight', 'InitialMagnification', 'fit'),
 set(gca, 'xlimmode','manual',...
           'ylimmode','manual',...
           'zlimmode','manual',...
           'climmode','manual',...
           'alimmode','manual'), ...
 set(gcf, 'doublebuffer', 'on');
 pause(pauseTime);
 imshow(h2.cdata,  'Border', 'tight', 'InitialMagnification', 'fit'), ...
 pause(pauseTime);
end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function BlinkLoop(C, K)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
Frequency = 10;
Seconds = 10;
SamplingFreq = 1000;
trigger = zeros (Seconds*SamplingFreq,1);
pauseTime = (1/(Frequency*2));

trigger = zeros (Seconds*SamplingFreq,1);
for i = 1:length(trigger)
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
    'Position', [0 0 1680 1050],...
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

% One full loop is 1 period. 

