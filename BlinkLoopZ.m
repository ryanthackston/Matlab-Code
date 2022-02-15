

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
hold on;
% Every loop is 1/Frequency seconds, i.e freq = 10, loop is .1


while tic
    if toc < 0.1
        if toc >= 0.05
            imshow(h1.cdata, 'Border', 'tight', 'InitialMagnification', 'fit')
        end
    elseif toc >= 0.10
        imshow(h2.cdata, 'Border', 'tight', 'InitialMagnification', 'fit')
    end
end
    
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

seconds = 20;
intervals = 10;
trial = seconds/intervals;
    for i = 1:2
    figure('Name','APP',...
    'Numbertitle','off',...
    'Position', [0 0 1680 1050],...
    'WindowStyle','normal',...
    'Color',[0.5 0.5 0.5],...
    'Toolbar','none')

    imshow(A(i).cdata, 'Border', 'tight', 'InitialMagnification', 'fit')
      F(i) = getframe(gcf) ;
      drawnow
    end
  % create the video writer with 1 fps
  writerObj = VideoWriter('C:\Users\Ryan\Desktop\myVideo12-15.avi');
  writerObj.FrameRate = 120;
  framerate = writerObj.FrameRate;
  % set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
trigger = zeros(framerate * seconds, 1);
for i= 1:(framerate * seconds) 
 j = i - currentTime;
 
    if i < 1200
        if j <= 5
            frame = F(1) ;
            writeVideo(writerObj, frame);
        elseif mod(i,5) == 0
            frame = F(2);
            writeVideo(writerObj, frame);
        end
    elseif i < 2400
        if mod(i,4) ~= 0
            frame = F(1) ;
            writeVideo(writerObj, frame);
            trigger(i) = 1;
        elseif mod(i,4) == 0
            frame = F(2);
            writeVideo(writerObj, frame);
            trigger(i) = 2
        end  
    end
    
end
% close the writer object
close(writerObj);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
seconds = 30;
Frequency1 = 10;
Frequency2 = 15;

    for i = 1:2
    figure('Name','APP',...
    'Numbertitle','off',...
    'Position', [0 0 1680 1050],...
    'WindowStyle','normal',...
    'Color',[0.5 0.5 0.5],...
    'Toolbar','none')


    m = 100; n = 160;
    [C, K] = checkerb(m, n);

    imshow(A(i).cdata, 'Border', 'tight', 'InitialMagnification', 'fit')
      F(i) = getframe(gcf) ;
      drawnow
    end
  % create the video writer with 1 fps
  writerObj = VideoWriter('C:\Users\Ryan\Desktop\myVideo10-15.avi');
  writerObj.FrameRate = 60;
  framerate = writerObj.FrameRate;
  % set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
trigger = zeros(framerate * seconds, 1);
currentTime = 0;
for i= 1:(framerate * seconds) 
    % if less than 30 seconds, blink checkerboard at 12 Hz
    if i < (framerate * seconds/2)
        if mod(i, (framerate / Frequency1)) == 0
            J = 0;
            currentTime = i;
            frame = F(1) ;
            writeVideo(writerObj, frame);
            trigger(i) = 1;
        else
            J = i - currentTime;
            if J < (framerate / (Frequency1 * 2) )
                frame = F(1) ;
                writeVideo(writerObj, frame);
                trigger(i) = 1;
            elseif J >= (framerate / (Frequency1 * 2) )
                frame = F(2) ;
                writeVideo(writerObj, frame);
                trigger(i) = 2;  
            end
        end
    % after 30 seconds, switch to 15 Hz    
    elseif i >= (framerate * seconds/2)
        if mod(i, (framerate/Frequency2) ) == 0
            J = 0;
            currentTime = i;
            frame = F(1) ;
            writeVideo(writerObj, frame);
            trigger(i) = 3;
        else
            J = i - currentTime;
            if J < framerate / (Frequency2 * 2)
                frame = F(1) ;
                writeVideo(writerObj, frame);
                trigger(i) = 3;
            elseif J >= framerate / (Frequency2 * 2)
                frame = F(2) ;
                writeVideo(writerObj, frame);
                trigger(i) = 4;  
            end
        end
    end
end

close(writerObj);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
seconds = 30;
Frequency1 = 10;
SamplingFreq = 1000;

    m = 100; n = 160;
    [C, K] = checkerb(m, n);
    imshow(C, 'Border', 'tight', 'InitialMagnification', 'fit'), ...
    h1 = getframe; 
    imshow(K,  'Border', 'tight', 'InitialMagnification', 'fit'), ...
    h2 = getframe;
    A = [h1 h2];
    
    for i = 1:2
    figure('Name','APP',...
    'Numbertitle','off',...
    'Position', [0 0 1680 1050],...
    'WindowStyle','normal',...
    'Color',[0.5 0.5 0.5],...
    'Toolbar','none')

    imshow(A(i).cdata, 'Border', 'tight', 'InitialMagnification', 'fit')
      F(i) = getframe(gcf) ;
      drawnow
    end
  % create the video writer with 1 fps
  writerObj = VideoWriter('C:\Users\Ryan\Desktop\myVideo10.avi');
  writerObj.FrameRate = 60;
  framerate = writerObj.FrameRate;
  % set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
A(1).trigger = [];
A(2).trigger = [];
A(3).trigger = zeros(framerate * seconds, 1);


currentTime = 0;
for i= 1:(framerate * seconds) 
        if mod(i, (framerate / Frequency1)) == 0
            J = 0;
            currentTime = i;
            frame = F(1) ;
            writeVideo(writerObj, frame);
            A(3).trigger(i) = 1;
        else
            J = i - currentTime;
            if J < (framerate / (Frequency1 * 2) )
                frame = F(1) ;
                writeVideo(writerObj, frame);
                A(3).trigger(i) = 1;
            elseif J >= (framerate / (Frequency1 * 2) )
                frame = F(2) ;
                writeVideo(writerObj, frame);
                A(3).trigger(i) = 2;  
            end
        end
end

close(writerObj);


    figure('Name','APP',...
    'Numbertitle','off',...
    'Position', [0 0 1680 1050],...
    'WindowStyle','normal',...
    'Color',[0.5 0.5 0.5],...
    'Toolbar','none')

    v = VideoReader('C:\Users\Ryan\Desktop\myVideo15.avi');
    ax = gca; ax.Position = [0 0 1 1];     v = VideoReader('C:\Users\Ryan\Desktop\myVideo10.avi');
    while hasFrame(v)
        tic
        vidFrame = readFrame(v);
        image(vidFrame, 'Parent', ax);
        ax.Visible = 'off';
        toc = 1/(v.FrameRate*2);
    end

%%
v = VLC();
v.Fullscreen = 'on';
v.play('C:\Users\Ryan\Desktop\myVideo10.avi') 








