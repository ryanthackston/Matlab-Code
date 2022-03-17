function varargout = slidingPLV2(varargin)
% SLIDINGPLV2 MATLAB code for slidingPLV2.fig
%      SLIDINGPLV2, by itself, creates a new SLIDINGPLV2 or raises the existing
%      singleton*.
%
%      H = SLIDINGPLV2 returns the handle to a new SLIDINGPLV2 or the handle to
%      the existing singleton*.
%
%      SLIDINGPLV2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SLIDINGPLV2.M with the given input arguments.
%
%      SLIDINGPLV2('Property','Value',...) creates a new SLIDINGPLV2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before slidingPLV2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to slidingPLV2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help slidingPLV2

% Last Modified by GUIDE v2.5 26-Feb-2022 19:06:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @slidingPLV2_OpeningFcn, ...
                   'gui_OutputFcn',  @slidingPLV2_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before slidingPLV2 is made visible.
function slidingPLV2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to slidingPLV2 (see VARARGIN)

% get the default parameter values for the sliders - (In my code I get
% frequency here.)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD THE DATA
% Load the data
% Motor img - calib
MIcal = {[('D:\OPM data - August 2019 - Colorado\BCI2000 mat files\imag_calibration\S1_Block1_imag_calib.mat');...
              ('D:\OPM data - August 2019 - Colorado\BCI2000 mat files\imag_calibration\S1_Block2_imag_calib.mat')],...
              [('D:\OPM data - August 2019 - Colorado\BCI2000 mat files\imag_calibration\S2_Block1_imag_calib.mat');...
              ('D:\OPM data - August 2019 - Colorado\BCI2000 mat files\imag_calibration\S2_Block2_imag_calib.mat')],...
              [('D:\OPM data - August 2019 - Colorado\BCI2000 mat files\imag_calibration\S3_Block1_imag_calib.mat');...
              ('D:\OPM data - August 2019 - Colorado\BCI2000 mat files\imag_calibration\S3_Block2_imag_calib.mat')],...
              [('D:\OPM data - August 2019 - Colorado\BCI2000 mat files\imag_calibration\S4_Block1_imag_calib.mat');...
              ('D:\OPM data - August 2019 - Colorado\BCI2000 mat files\imag_calibration\S4_Block2_imag_calib.mat')],...
              [('D:\OPM data - August 2019 - Colorado\BCI2000 mat files\imag_calibration\S5_Block1_imag_calib.mat');...
              ('D:\OPM data - August 2019 - Colorado\BCI2000 mat files\imag_calibration\S5_Block2_imag_calib.mat')],...
              [('D:\OPM data - August 2019 - Colorado\BCI2000 mat files\imag_calibration\S6_Block1_imag_calib.mat');...
              ('D:\OPM data - August 2019 - Colorado\BCI2000 mat files\imag_calibration\S6_Block2_imag_calib.mat')],...
              [('D:\OPM data - August 2019 - Colorado\BCI2000 mat files\imag_calibration\S7_Block1_imag_calib.mat');...
              ('D:\OPM data - August 2019 - Colorado\BCI2000 mat files\imag_calibration\S7_Block2_imag_calib.mat')],...
              [('D:\OPM data - August 2019 - Colorado\BCI2000 mat files\imag_calibration\S8_Block1_imag_calib.mat');...
              ('D:\OPM data - August 2019 - Colorado\BCI2000 mat files\imag_calibration\S8_Block2_imag_calib.mat')],...
              [('D:\OPM data - August 2019 - Colorado\BCI2000 mat files\imag_calibration\S9_Block1_imag_calib.mat');...
              ('D:\OPM data - August 2019 - Colorado\BCI2000 mat files\imag_calibration\S9_Block2_imag_calib.mat')]};
          
MIcont = {[('D:\OPM data - August 2019 - Colorado\BCI2000 mat files\imag_control\S1_Block1_imag_control.mat');...
                ('D:\OPM data - August 2019 - Colorado\BCI2000 mat files\imag_control\S1_Block2_imag_control.mat')],...
                [('D:\OPM data - August 2019 - Colorado\BCI2000 mat files\imag_control\S2_Block1_imag_control.mat');...
                ('D:\OPM data - August 2019 - Colorado\BCI2000 mat files\imag_control\S2_Block2_imag_control.mat')],...
                [('D:\OPM data - August 2019 - Colorado\BCI2000 mat files\imag_control\S3_Block1_imag_control.mat');...
                ('D:\OPM data - August 2019 - Colorado\BCI2000 mat files\imag_control\S3_Block2_imag_control.mat')],...
                [('D:\OPM data - August 2019 - Colorado\BCI2000 mat files\imag_control\S4_Block1_imag_control.mat');...
                ('D:\OPM data - August 2019 - Colorado\BCI2000 mat files\imag_control\S4_Block2_imag_control.mat')],...
                [('D:\OPM data - August 2019 - Colorado\BCI2000 mat files\imag_control\S5_Block1_imag_control.mat');...
                ('D:\OPM data - August 2019 - Colorado\BCI2000 mat files\imag_control\S5_Block2_imag_control.mat')],...
                [('D:\OPM data - August 2019 - Colorado\BCI2000 mat files\imag_control\S6_Block1_imag_control.mat');...
                ('D:\OPM data - August 2019 - Colorado\BCI2000 mat files\imag_control\S6_Block2_imag_control.mat')],...
                [('D:\OPM data - August 2019 - Colorado\BCI2000 mat files\imag_control\S7_Block1_imag_control.mat');...
                ('D:\OPM data - August 2019 - Colorado\BCI2000 mat files\imag_control\S7_Block2_imag_control.mat')],...
                [('D:\OPM data - August 2019 - Colorado\BCI2000 mat files\imag_control\S8_Block1_imag_control.mat');...
                ('D:\OPM data - August 2019 - Colorado\BCI2000 mat files\imag_control\S8_Block2_imag_control.mat')],...
                [('D:\OPM data - August 2019 - Colorado\BCI2000 mat files\imag_control\S9_Block1_imag_control.mat');...
                ('D:\OPM data - August 2019 - Colorado\BCI2000 mat files\imag_control\S9_Block2_imag_control.mat')]};
            
questions = { '\fontsize{12} Select Task: [1] Motor Img Calib [2] Motor Img Control';...
                    '\fontsize{12} Which subject number would you like to analyze (1-9)?';...
                    '\fontsize{12} What frequency would you like to analyze first?';...
                    '\fontsize{12} What frequency would like to analyze last?';...
                    '\fontsize{12} How much should the frequency increment by?';...
                    '\fontsize{12} What should be the frequency width for the filter?';...
                    '\fontsize{12} What should be the time window size (seconds)?';...
                    '\fontsize{12} What time do you want to start at?'};

definput = cell(8,1);
definput(:) = [{'1'};{'9'};{'23.1'};{'26'};{'0.5'};{'4'};{'1'};{'1'}];
options.Interpreter = 'tex';           
answers = inputdlg(questions, 'Parameters for PLV', [1 70], definput, options )

            
% Add as input for later
 % stud{1} -> MIcal; stud{2} -> MIcont
stud = {MIcal; MIcont};  
stud = stud{str2double(answers{1})};
subj = str2double(answers{2});
blocks = [load(stud{subj}(1,:)); load(stud{subj}(2,:))]; 
% Motot img - control
% blocks = [load('D:\OPM data - August 2019 - Colorado\BCI2000 mat files\imag_control\S2_Block1_imag_control.mat')...
%           load('D:\OPM data - August 2019 - Colorado\BCI2000 mat files\imag_control\S2_Block2_imag_control.mat')];

              
RestOnset = [blocks(1).RestOnset(:); blocks(2).RestOnset(:)+253000];
MoveOnset = [blocks(1).MoveOnset(:); blocks(2).MoveOnset(:)+253000];

% Delete unused channels
blocks(1).meg(:,[4 6 7 11 15 16 20 23]) = [];
blocks(2).meg(:,[4 6 7 11 15 16 20 23]) = [];

meg = [blocks(1).meg; blocks(2).meg];
meg = meg';

% Rearrange the channel numbers to go sequentially from top-left to bottom-right
megMat = [];
megMat(1,:) = meg(12,:);
megMat(2,:) = meg(7,:);
megMat(3,:) = meg(6,:);
megMat(4,:) = meg(17,:);
megMat(5,:) = meg(14,:);
megMat(6,:) = meg(8,:);
megMat(7,:) = meg(11,:);
megMat(8,:) = meg(16,:);
megMat(9,:) = meg(5,:);
megMat(10,:) = meg(10,:);
megMat(11,:) = meg(4,:);
megMat(12,:) = meg(15,:);
megMat(13,:) = meg(9,:);
megMat(14,:) = meg(3,:);
megMat(15,:) = meg(1,:);
megMat(16,:) = meg(13,:);
megMat(17,:) = meg(2,:);
meg = megMat(:,:);

% Rearrange data to be channels x time points x trials 
% data(:, 1:5000, :) time points is rest task. 
% data(:, 5001:10000, :) time points is move task.
megR3d = zeros(17,5000,size(RestOnset,1));
megM3d = zeros(17,5000,size(MoveOnset,1));
data = zeros(17,10000,size(RestOnset,1));

% motor img - calib
for i = 1:size(RestOnset,1)
    megR3d(:,:,i) = meg(:,(RestOnset(i):MoveOnset(i)-1));
    megM3d(:,:,i) = meg(:,(MoveOnset(i):MoveOnset(i)+4999));
    data(:,:,i) = [megR3d(:,:,i) megM3d(:,:,i)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS

% Beginning frequency to analyze
begF = str2double(answers{3});
% End frequency to analyze
endF = str2double(answers{4});
% Increment steps by
stepf = str2double(answers{5});
fi = [stepf/(endF-begF) 1.0];
%FWHM for filterdat
freqwidt = str2double(answers{6});
fwMax = 2.0;
fwMin = 0.1;
fwi = [ 0.2/(fwMax-fwMin) 1.0];
% time window parameters (in seconds) - leave out 1st sec for each task
tws = str2double(answers{7});
timewin = cell(2,1);
timestart = str2double(answers{8});
timewin{1} = [timestart  timestart+tws];
timewin{2} = [timestart+5 timestart+5+tws];
% time step increment
tstepi = [ 0.5/(5.0-0.1) 1.0];
% sample rate
srate = 1000;
% Frequencies to analyze
% freqs = linspace(begF,endF,incr);


set(handles.slider_freq, 'Max', endF);
set(handles.slider_freq, 'Min', begF);
set(handles.slider_freq, 'Value', begF);
set(handles.slider_freq, 'SliderStep', fi);
centfreq = get(handles.slider_freq, 'Value');
set(handles.text_centfreq, 'String', ['Center Frequency = ' num2str(handles.slider_freq.Value) ' Hz']);

set(handles.slider_width, 'Max', fwMax);
set(handles.slider_width, 'Min', fwMin);
set(handles.slider_width, 'Value', freqwidt);
set(handles.slider_width, 'SliderStep', fwi);
freqwidt = get(handles.slider_width, 'Value');
set(handles.text_freqwidt, 'String', ['Frequency Width = ' num2str(handles.slider_width.Value)]);

set(handles.slider_timewin_size, 'Max', 5.0);
set(handles.slider_timewin_size, 'Min', 0.1);
set(handles.slider_timewin_size, 'Value', tws);
set(handles.slider_timewin_size, 'SliderStep', tstepi);
timeWinSize = get(handles.slider_timewin_size, 'Value');
set(handles.text_timewinSize, 'String', ['Time Window Size = ' num2str(handles.slider_timewin_size.Value)]);

set(handles.slider_time, 'Max', 5.0);
set(handles.slider_time, 'Min', 0.1);
set(handles.slider_time, 'Value', timestart);
set(handles.slider_time, 'SliderStep', tstepi);
timeWin = get(handles.slider_time, 'Value');
set(handles.text_time, 'String', ['Starting Time = ' num2str(handles.slider_time.Value) ' Sec']);


% time vector
timevec = 0+1/srate:(1/(srate)):10;

% convert to time indices (sec)
tidx1 = dsearchn(timevec',timewin{1}');
tidx2 = dsearchn(timevec',timewin{2}');

synchmat = zeros(2,17,17);
filtdat = filterFGx(data,srate, centfreq, freqwidt);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE INITIAL PLOTS HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%

% REST SYNCHRONIZATION INITIAL PLOT
        handles.restTW.XLim = [ 0.5 17.5 ];
        handles.restTW.YLim = [ 0.5 17.5 ];
        colorbar(handles.restTW)
        handles.restTW.Colormap = lbmap(handles.restTW, 'RedBlue')
        handles.restTW.CLim = [0 0.9];
        % Make the axes square
        handles.restTW.PlotBoxAspectRatio = [1 1 1];
        % drawGrid(handles.restTW, 0.75);
        handles.plotR = imagesc( handles.restTW, 'CData', squeeze(synchmat(1,:,:))) 
        handles.plotR.Parent.Title.String= [ 'Synch: ' num2str(timewin{1}(1)) '-' num2str(timewin{1}(2)) ' Sec, '  '? \pm ? Hz'];
        handles.plotR.Parent.Title.Position = [9.0000 17.85 0]
        % Plot grid lines 
        hold (handles.plotR.Parent,'on');
            for j = 1:17
                pRy{j} = plot(handles.plotR.Parent, [0.5, 18.5],[j-0.5, j-0.5],'k-');
                pRx{j} = plot(handles.plotR.Parent,[j-0.5, j-0.5],[0.5, 18.5],'k-');
            end
        hold (handles.plotR.Parent,'off')


% MOVEMENT SYNCHRONIZATION INITIAL PLOT
        handles.moveTW.XLim = [0.5 17.5];
        handles.moveTW.YLim = [0.5 17.5];
        colorbar(handles.moveTW)
        handles.moveTW.Colormap = lbmap(handles.moveTW, 'RedBlue')
        handles.moveTW.CLim = [0 0.9];
        % Make the axes square
        handles.moveTW.PlotBoxAspectRatio = [1 1 1];
        handles.plotM = imagesc( handles.moveTW, 'CData', squeeze(synchmat(2,:,:)));
        handles.plotM.Parent.Title.String = [ 'Synch: ' num2str(timewin{2}(1)) '-' num2str(timewin{2}(2)) ' Sec, '  '? \pm ? Hz'];
        handles.plotM.Parent.Title.Position = [9.0000 17.85 0]
        % Plot grid lines 
        hold (handles.plotM.Parent,'on')
            for j = 1:17
                pMy{j} = plot(handles.moveTW, [0.5,18.5],[j-.5,j-.5],'k-');
                pMx{j} = plot(handles.moveTW,[j-.5,j-.5],[0.5, 18.5],'k-');
            end
        hold (handles.plotM.Parent,'off')


% PHASE DIFFERENCE SYNCHRONIZATION PLOT
        handles.phaseDif.XLim = [0.5 17.5];
        handles.phaseDif.YLim = [0.5 17.5];
        colorbar(handles.phaseDif)
        handles.phaseDif.Colormap = lbmap(handles.phaseDif, 'RedBlue')
        handles.phaseDif.CLim = [-1 1];
        % Make the axes square
        handles.phaseDif.PlotBoxAspectRatio = [1 1 1];
        handles.plotD = imagesc( handles.phaseDif, 'CData', squeeze(diff(synchmat)))
        handles.plotD.Parent.Title.String = ['Phase Synchronization Difference: Move - Rest'];
        handles.plotD.Parent.Title.Position = [9.0000 17.85 0]
        hold (handles.plotD.Parent,'on');
            for j = 1:17
                pDy{j} = plot(handles.phaseDif, [0.5,18.5],[j-.5,j-.5],'k-');
                pDx{j} = plot(handles.phaseDif,[j-.5,j-.5],[0.5, 18.5],'k-');
            end
        hold (handles.plotD.Parent,'off');

% KRUSKAL WALLIS STASTICALLY SIGNIFICANT PHASE SYNCHRONIZATION
        handles.krusWal.XLim = [0.5 17.5];
        handles.krusWal.YLim = [0.5 17.5];
        colorbar(handles.krusWal)
        handles.krusWal.Colormap = lbmap(handles.krusWal, 'RedBlue')
        handles.krusWal.CLim = [0    0.01];
        % Mark Statistically Significant Phase Synchs in Kruskal-Wallis Test
        hori = repmat(1:size(data,1),size(data,1),1);
        vert = hori';
        g = text(hori(:), vert(:), ' ', 'HorizontalAlignment', 'Center', 'Color','k', 'FontWeight', 'bold', 'Parent', handles.krusWal);
        % Make the axes square
        handles.krusWal.PlotBoxAspectRatio = [1 1 1];
        handles.plotSS = imagesc( handles.krusWal, 'CData', squeeze(diff(synchmat)))
        handles.plotSS.Parent.Title.String = ['Kruskal-Wallis Stat. Sig. \alpha = 0.05/17'];
        handles.plotSS.Parent.Title.Position = [9.0000 17.85 0]
        hold (handles.plotSS.Parent,'on');
            for j = 1:17
                pSy{j} = plot(handles.krusWal, [.5, 18.5],[j-.5,j-.5],'k-');
                pSx{j} = plot(handles.krusWal,[j-.5,j-.5],[.5,18.5],'k-');
            end
        hold (handles.plotSS.Parent,'off');









% (In my code I want to compute PLV here)
    % angle time-series
    angst = zeros(size(filtdat));
    
%     How do I take out this for loop?
    for triali=1:size(data,3)
        % Applying hilbert transform to matrix inputs, phase angles should be
        % sawtooth slanted to the right
        angst(:,:,triali) = angle(hilbert(squeeze(filtdat(:,:,triali))').');
        % Plot angle time series and check phase angles
        % figure(5); plot(angst(5,:,10));
    end
    
    nchans = size(data,1);
     for chani=1:nchans
        for chanj=1:nchans
            for tria = 1:size(data,3)

                %%% time window 1
                % extract angles in channel i, all trials
                tmpAi = angst(chani,tidx1(1):tidx1(2),tria);
                % extract angles in channel j all trials
                tmpAj = angst(chanj,tidx1(1):tidx1(2),tria);
                % compute synchronization between the two channels
                % abs( mean( eulers formula( phase angle diff)), 2nd dim: time points)
                trialsynch = abs(mean(exp(1i*( tmpAi-tmpAj )),2));
                % average over trials - synchmat 1st Time Window
                synchmat1(tria,chani,chanj) = mean(trialsynch);

                %%% time window 2
                % extract angles
                tmpAi = angst(chani,tidx2(1):tidx2(2),tria);
                tmpAj = angst(chanj,tidx2(1):tidx2(2),tria);
                % Phase Locking Value on each trial
                trialsynch = abs(mean(exp(1i*( tmpAi-tmpAj )),2));
                % average over trials - synchmat 2nd Time Window

                synchmat2(tria,chani,chanj) = mean(trialsynch); 
            end            
        end
     end
    % synchmat Difference - Time Window x channels x channels 
    synchmatD = [mean(synchmat1,1); mean(synchmat2,1)];
%     synchmatDif(1,chani,chanj) = mean((synchmat2(:,chani,chanj) - synchmat1(:,chani,chanj)),1);
    % diffS = synchmatD(2,:,:) - synchmatD(1,:,:);
    % diffS == diff(synchmatD)
%     subplot(1,4,1)
%     title([ 'Synch: ' num2str(timewin{1}(1)) '-' num2str(timewin{1}(2)) ' Sec, ' num2str(centfreq) ' Hz \pm ' num2str(freqwidt)]);

%% UPDATE THE IMAGES
% Update Rest Trials Plot
        handles.plotR = imagesc( handles.restTW, 'CData', squeeze(mean(synchmat1,1))' ) 
        handles.plotR.Parent.Title.String =  [ 'Rest Phase Synch: ' num2str(timewin{1}(1)) '-' num2str(timewin{1}(2)) ' Sec, ' num2str(centfreq) ' \pm ' num2str(freqwidt) ' Hz'];
        handles.plotR.Parent.Title.FontSize = 13.5;
        handles.plotM.Parent.Title.FontSize = 13.5;
        [p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17] = deal(pRx{:});
        uistack([p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17], 'top');
        [p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17] = deal(pRy{:});
        uistack([p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17], 'top');

% Update Move Trials Plot
        handles.plotM = imagesc( handles.moveTW, 'CData', squeeze(mean(synchmat2,1)) )
        handles.plotM.Parent.Title.String =  [ 'Move Phase Synch: ' num2str(timewin{2}(1)) '-' num2str(timewin{2}(2)) ' Sec, ' num2str(centfreq) ' Hz \pm ' num2str(freqwidt) ' Hz'];
        handles.plotM.Parent.Title.FontSize = 13.5;
        [p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17] = deal(pMx{:});
        uistack([p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17], 'top');
        [p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17] = deal(pMy{:});
        uistack([p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17], 'top');
  
% Update Phase Difference Plot
        handles.plotD = imagesc( handles.phaseDif, 'CData', squeeze(diff(synchmatD)) ) 
%         handles.plotD.Parent.Title.String = ['Phase Synchronization Difference: Move - Rest'];
        handles.plotD.Parent.Title.FontSize = 13.5;
        [p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17] = deal(pDx{:});
        uistack([p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17], 'top');
        [p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17] = deal(pDy{:});
        uistack([p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17], 'top');
    
    % Kruskal-Wallis test on Phase-Locking Value data
    for ii = 1:17
        for jj = 1:17
            [krusP(ii,jj), tbl{ii,jj}, stats{ii,jj}] = kruskalwallis([synchmat1(:,ii,jj) synchmat2(:,ii,jj)], [], 'off');
        end
    end
    % Arrange kruskal-wallis values for labeling in image
    kk = num2cell(krusP);
    kk = cellfun(@num2str, kk, 'UniformOutput', false); % convert to string  

%     ss = krusP <= alphaT;
%     ssi = find(ss);
    
% Update Kruskal-Wallis Statistical Test Plot
        handles.plotSS = imagesc( handles.krusWal, 'CData', krusP )
        handles.plotSS.Parent.Title.FontSize = 13.5;
        [p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17] = deal(pSx{:});
        uistack([p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17], 'top');
        [p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17] = deal(pSy{:});
        uistack([p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17], 'top');
    
    for tt = 1:length(g)
        % Reset previous stat sig labels
        g(tt).String = ' ';
        alphaT = 0.05/size(data,1);
        if str2double(kk{tt})<= alphaT
            % Label images that are stat sig
             g(tt) = text(hori(tt), vert(tt), 'X', 'HorizontalAlignment', 'Center', 'Color','w', 'FontWeight', 'bold', 'FontSize', 12.5, 'Parent', handles.krusWal);
        else
            continue;
        end
    end




% Choose default command line output for slidingPLV2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes slidingPLV2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = slidingPLV2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider_freq_Callback(hObject, eventdata, handles)
% hObject    handle to slider_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updatePhaseSync(handles);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider_freq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_width_Callback(hObject, eventdata, handles)
% hObject    handle to slider_width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updatePhaseSync(handles);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider_width_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_timewin_size_Callback(hObject, eventdata, handles)
% hObject    handle to slider_timewin_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updatePhaseSync(handles);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider_timewin_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_timewin_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_time_Callback(hObject, eventdata, handles)
% hObject    handle to slider_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updatePhaseSync(handles);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function updatePhaseSync(handles)

    blocks = [load('D:\OPM data - August 2019 - Colorado\BCI2000 mat files\imag_calibration\S9_Block1_imag_calib.mat')...
                  load('D:\OPM data - August 2019 - Colorado\BCI2000 mat files\imag_calibration\S9_Block2_imag_calib.mat')];
              
              RestOnset = [blocks(1).RestOnset(:); blocks(2).RestOnset(:)+253000];
MoveOnset = [blocks(1).MoveOnset(:); blocks(2).MoveOnset(:)+253000];

% Delete unused channels
blocks(1).meg(:,[4 6 7 11 15 16 20 23]) = [];
blocks(2).meg(:,[4 6 7 11 15 16 20 23]) = [];

meg = [blocks(1).meg; blocks(2).meg];
meg = meg';

% Rearrange the channel numbers to go sequentially from top-left to bottom-right
megMat = [];
megMat(1,:) = meg(12,:);
megMat(2,:) = meg(7,:);
megMat(3,:) = meg(6,:);
megMat(4,:) = meg(17,:);
megMat(5,:) = meg(14,:);
megMat(6,:) = meg(8,:);
megMat(7,:) = meg(11,:);
megMat(8,:) = meg(16,:);
megMat(9,:) = meg(5,:);
megMat(10,:) = meg(10,:);
megMat(11,:) = meg(4,:);
megMat(12,:) = meg(15,:);
megMat(13,:) = meg(9,:);
megMat(14,:) = meg(3,:);
megMat(15,:) = meg(1,:);
megMat(16,:) = meg(13,:);
megMat(17,:) = meg(2,:);
meg = megMat(:,:);

% Rearrange data to be channels x time points x trials 
% data(:, 1:5000, :) time points is rest task. 
% data(:, 5001:10000, :) time points is move task.
megR3d = zeros(17,5000,size(RestOnset,1));
megM3d = zeros(17,5000,size(MoveOnset,1));
data = zeros(17,10000,size(RestOnset,1));

% motor img - calib
for i = 1:size(RestOnset,1)
    megR3d(:,:,i) = meg(:,(RestOnset(i):MoveOnset(i)-1));
    megM3d(:,:,i) = meg(:,(MoveOnset(i):MoveOnset(i)+4999));
    data(:,:,i) = [megR3d(:,:,i) megM3d(:,:,i)];
end

    % get the new parameter values
    centfreq = get(handles.slider_freq, 'Value');
    set(handles.text_centfreq, 'String', ['Center Frequency = ' num2str(handles.slider_freq.Value) ' Hz']);
    
    freqwidt = get(handles.slider_width, 'Value');
    set(handles.text_freqwidt, 'String', ['Frequency Width = ' num2str(handles.slider_width.Value) ' Hz' ]);


    tws = get(handles.slider_timewin_size, 'Value');
    set(handles.text_timewinSize, 'String', ['Time Window Size = ' num2str(handles.slider_timewin_size.Value) ' Sec']);

    timestart = get(handles.slider_time, 'Value');
    set(handles.text_time, 'String', ['Starting Time = ' num2str(handles.slider_time.Value) ' Sec']);

    timewin{1} = [timestart  timestart+tws];
    timewin{2} = [timestart+5 timestart+5+tws];

    srate = 1000;
    % time vector
    timevec = 0+1/srate:(1/(srate)):10;


    % convert to time indices (sec)
    tidx1 = dsearchn(timevec',timewin{1}');
    tidx2 = dsearchn(timevec',timewin{2}');

    synchmat = zeros(2,17,17);
    filtdat = filterFGx(data,srate, centfreq, freqwidt);
    
    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% % PLOTS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% % REST SYNCHRONIZATION INITIAL PLOT
        handles.restTW.XLim = [ 0.5 17.5 ];
        handles.restTW.YLim = [ 0.5 17.5 ];
        colorbar(handles.restTW)
        handles.restTW.Colormap = lbmap(handles.restTW, 'RedBlue')
        handles.restTW.CLim = [0 0.9];
        % Make the axes square
        handles.restTW.PlotBoxAspectRatio = [1 1 1];
        % drawGrid(handles.restTW, 0.75);
        handles.plotR = imagesc( handles.restTW, 'CData', squeeze(synchmat(1,:,:))) 
        handles.plotR.Parent.Title.String= [ 'Synch: ' num2str(timewin{1}(1)) '-' num2str(timewin{1}(2)) ' Sec, '  '? \pm ? Hz'];
        handles.plotR.Parent.Title.Position = [9.0000 17.85 0]
        % Plot grid lines 
        hold (handles.plotR.Parent,'on');
            for j = 1:17
                pRy{j} = plot(handles.plotR.Parent, [0.5, 18.5],[j-0.5, j-0.5],'k-');
                pRx{j} = plot(handles.plotR.Parent,[j-0.5, j-0.5],[0.5, 18.5],'k-');
            end
        hold (handles.plotR.Parent,'off')


% % MOVEMENT SYNCHRONIZATION INITIAL PLOT
        handles.moveTW.XLim = [0.5 17.5];
        handles.moveTW.YLim = [0.5 17.5];
        colorbar(handles.moveTW)
        handles.moveTW.Colormap = lbmap(handles.moveTW, 'RedBlue')
        handles.moveTW.CLim = [0 0.9];
        % Make the axes square
        handles.moveTW.PlotBoxAspectRatio = [1 1 1];
        handles.plotM = imagesc( handles.moveTW, 'CData', squeeze(synchmat(2,:,:)));
        handles.plotM.Parent.Title.String = [ 'Synch: ' num2str(timewin{2}(1)) '-' num2str(timewin{2}(2)) ' Sec, '  '? \pm ? Hz'];
        handles.plotM.Parent.Title.Position = [9.0000 17.85 0]
        % Plot grid lines 
        hold (handles.plotM.Parent,'on')
            for j = 1:17
                pMy{j} = plot(handles.moveTW, [0.5,18.5],[j-.5,j-.5],'k-');
                pMx{j} = plot(handles.moveTW,[j-.5,j-.5],[0.5, 18.5],'k-');
            end
        hold (handles.plotM.Parent,'off')


% % PHASE DIFFERENCE SYNCHRONIZATION PLOT
        handles.phaseDif.XLim = [0.5 17.5];
        handles.phaseDif.YLim = [0.5 17.5];
        colorbar(handles.phaseDif)
        handles.phaseDif.Colormap = lbmap(handles.phaseDif, 'RedBlue')
        handles.phaseDif.CLim = [-1 1];
        % Make the axes square
        handles.phaseDif.PlotBoxAspectRatio = [1 1 1];
        handles.plotD = imagesc( handles.phaseDif, 'CData', squeeze(diff(synchmat)))
        handles.plotD.Parent.Title.String = ['Phase Synchronization Difference: Move - Rest'];
        handles.plotD.Parent.Title.Position = [9.0000 17.85 0]
        hold (handles.plotD.Parent,'on');
            for j = 1:17
                pDy{j} = plot(handles.phaseDif, [0.5,18.5],[j-.5,j-.5],'k-');
                pDx{j} = plot(handles.phaseDif,[j-.5,j-.5],[0.5, 18.5],'k-');
            end
        hold (handles.plotD.Parent,'off');


% % KRUSKAL WALLIS STASTICALLY SIGNIFICANT PHASE SYNCHRONIZATION

        handles.krusWal.XLim = [0.5 17.5];
        handles.krusWal.YLim = [0.5 17.5];
        colorbar(handles.krusWal)
        handles.krusWal.Colormap = lbmap(handles.krusWal, 'RedBlue')
        handles.krusWal.CLim = [0    0.01];
        % Mark Statistically Significant Phase Synchs in Kruskal-Wallis Test
        hori = repmat(1:size(data,1),size(data,1),1);
        vert = hori';
        g = text(hori(:), vert(:), ' ', 'HorizontalAlignment', 'Center', 'Color','k', 'FontWeight', 'bold', 'Parent', handles.krusWal);
        % Make the axes square
        handles.krusWal.PlotBoxAspectRatio = [1 1 1];
        handles.plotSS = imagesc( handles.krusWal, 'CData', squeeze(diff(synchmat)))
        handles.plotSS.Parent.Title.String = ['Kruskal-Wallis Stat. Sig. \alpha = 0.05/17'];
        handles.plotSS.Parent.Title.Position = [9.0000 17.85 0]
        hold (handles.plotSS.Parent,'on');
        pSy = cell(17,1);
        pSx = cell(1,17);
            for j = 1:17
                pSy{j} = plot(handles.krusWal, [.5, 18.5],[j-.5,j-.5],'k-');
                pSx{j} = plot(handles.krusWal,[j-.5,j-.5],[.5,18.5],'k-');
            end
        hold (handles.plotSS.Parent,'off');
        
    for j = 1:17
%         pRy{j} = plot(handles.plotR, [0.5, 18.5],[j-0.5, j-0.5],'k-');
%         pRx{j} = plot(handles.plotR,[j-0.5, j-0.5],[0.5, 18.5],'k-');
%         pMy{j} = plot(handles.moveTW, [0.5, 18.5],[j-0.5, j-0.5],'k-');
%         pMx{j} = plot(handles.moveTW,[j-0.5, j-0.5],[0.5, 18.5],'k-');
%         pDy{j} = plot(handles.phaseDif, [0.5, 18.5],[j-0.5, j-0.5],'k-');
%         pDx{j} = plot(handles.phaseDif,[j-0.5, j-0.5],[0.5, 18.5],'k-');
%         pSSy{j} = plot(handles.krusWal, [0.5, 18.5],[j-0.5, j-0.5],'k-');
%         pSSx{j} = plot(handles.krusWal,[j-0.5, j-0.5],[0.5, 18.5],'k-');

    end


    nchans = size(data,1);
       for triali=1:size(data,3)
            % Applying hilbert transform to matrix inputs
            angst(:,:,triali) = angle(hilbert(squeeze(filtdat(:,:,triali))').');
       end

         for chani=1:nchans
            for chanj=1:nchans
                for tria = 1:size(data,3)

                    %%% time window 1
                    % extract angles in channel i, all trials
                    tmpAi = angst(chani,tidx1(1):tidx1(2),tria);
                    % extract angles in channel j all trials
                    tmpAj = angst(chanj,tidx1(1):tidx1(2),tria);
                    % compute synchronization between the two channels
                    % abs( mean( eulers formula( phase angle diff)), 2nd dim: time points)
                    trialsynch = abs(mean(exp(1i*( tmpAi-tmpAj )),2));
                    % average over trials - synchmat 1st Time Window
                    synchmat1(tria,chani,chanj) = mean(trialsynch);

                    %%% time window 2
                    % extract angles
                    tmpAi = angst(chani,tidx2(1):tidx2(2),tria);
                    tmpAj = angst(chanj,tidx2(1):tidx2(2),tria);
                    % Phase Locking Value on each trial
                    trialsynch = abs(mean(exp(1i*( tmpAi-tmpAj )),2));
                    % average over trials - synchmat 2nd Time Window

                    synchmat2(tria,chani,chanj) = mean(trialsynch); 
                end            
            end
         end
        % synchmat Difference - Time Window x channels x channels 
        synchmatD = [mean(synchmat1,1); mean(synchmat2,1)];
    %     synchmatDif(1,chani,chanj) = mean((synchmat2(:,chani,chanj) - synchmat1(:,chani,chanj)),1);
        % diffS = synchmatD(2,:,:) - synchmatD(1,:,:);
        % diffS == diff(synchmatD)
    %     subplot(1,4,1)
    %     title([ 'Synch: ' num2str(timewin{1}(1)) '-' num2str(timewin{1}(2)) ' Sec, ' num2str(centfreq) ' Hz \pm ' num2str(freqwidt)]);

    %% UPDATE THE IMAGES
    % Update Rest Trials Plot
%             set(handles.plotR, 'CData', squeeze(mean(synchmat1,1))' );
            handles.restTW = imagesc( handles.restTW, 'CData', squeeze(mean(synchmat1,1))' ) 
            handles.restTW.Parent.Title.String =  [ 'Rest Phase Synch: ' num2str(timewin{1}(1)) '-' num2str(timewin{1}(2)) ' Sec, ' num2str(centfreq) ' \pm ' num2str(freqwidt) ' Hz'];
            handles.restTW.Parent.Title.FontSize = 13.5;
            [p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17] = deal(pRx{:});
            uistack([p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17], 'top');
            [p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17] = deal(pRy{:});
            uistack([p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17], 'top');

    % Update Move Trials Plot
%              set(handles.plotM, 'CData', squeeze(mean(synchmat2,1))' );

            handles.plotM = imagesc( handles.moveTW, 'CData', squeeze(mean(synchmat2,1))' )
            handles.plotM.Parent.Title.String =  [ 'Move Phase Synch: ' num2str(timewin{2}(1)) '-' num2str(timewin{2}(2)) ' Sec, ' num2str(centfreq) ' Hz \pm ' num2str(freqwidt) ' Hz'];
            handles.plotM.Parent.Title.FontSize = 13.5;
            [p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17] = deal(pMx{:});
            uistack([p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17], 'top');
            [p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17] = deal(pMy{:});
            uistack([p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17], 'top');

    % Update Phase Difference Plot
%             set(handles.plotD, 'CData', squeeze(diff(synchmatD)) );
            handles.plotD = imagesc( handles.phaseDif, 'CData', squeeze(diff(synchmatD)) ) 
    %         handles.plotD.Parent.Title.String = ['Phase Synchronization Difference: Move - Rest'];
            handles.plotD.Parent.Title.FontSize = 13.5;
            [p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17] = deal(pDx{:});
            uistack([p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17], 'top');
            [p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17] = deal(pDy{:});
            uistack([p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17], 'top');

        % Kruskal-Wallis test on Phase-Locking Value data
        for ii = 1:17
            for jj = 1:17
                [krusP(ii,jj), tbl{ii,jj}, stats{ii,jj}] = kruskalwallis([synchmat1(:,ii,jj) synchmat2(:,ii,jj)], [], 'off');
            end
        end
        % Arrange kruskal-wallis values for labeling in image
        kk = num2cell(krusP);
        kk = cellfun(@num2str, kk, 'UniformOutput', false); % convert to string  

    %     ss = krusP <= alphaT;
    %     ssi = find(ss);
    
    % Update Kruskal-Wallis Statistical Test Plot
%         set(handles.plotSS, 'CData', krusP );
            handles.plotSS = imagesc( handles.krusWal, 'CData', krusP )
            handles.plotSS.Parent.Title.FontSize = 13.5;
            [p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17] = deal(pSx{:});
            uistack([p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17], 'top');
            [p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17] = deal(pSy{:});
            uistack([p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17], 'top');
        
        for tt = 1:length(g)
            % Reset previous stat sig labels
            g(tt).String = ' ';
            alphaT = 0.05/size(data,1);
            if str2double(kk{tt} ) <= alphaT
                % Label images that are stat sig
                 g(tt) = text(hori(tt), vert(tt), 'X', 'HorizontalAlignment', 'Center', 'Color','w', 'FontWeight', 'bold', 'FontSize', 12.5, 'Parent', handles.krusWal);
            else
                continue;
            end
        end

'Done'


% % (In my code I want to compute PLV here)
% sigmoid = a ./ (1 + exp(-b*(handles.x-c)) );
% 
% % update the sigmoid handle
% set(handles.plots, 'YData', sigmoid);
% 
% % update the lines
% set(handles.plotc, 'XData', [c c])
% set(handles.plota, 'YData', [1 1]*a/2)
