function varargout = slidingPLV_Aug2021(varargin)
% SLIDINGPLV_AUG2021 MATLAB code for slidingPLV_Aug2021.fig
%      SLIDINGPLV_AUG2021, by itself, creates a new SLIDINGPLV_AUG2021 or raises the existing
%      singleton*.
%
%      H = SLIDINGPLV_AUG2021 returns the handle to a new SLIDINGPLV_AUG2021 or the handle to
%      the existing singleton*.
%
%      SLIDINGPLV_AUG2021('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SLIDINGPLV_AUG2021.M with the given input arguments.
%
%      SLIDINGPLV_AUG2021('Property','Value',...) creates a new SLIDINGPLV_AUG2021 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before slidingPLV_Aug2021_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to slidingPLV_Aug2021_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help slidingPLV_Aug2021

% Last Modified by GUIDE v2.5 23-May-2022 16:17:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @slidingPLV_Aug2021_OpeningFcn, ...
                   'gui_OutputFcn',  @slidingPLV_Aug2021_OutputFcn, ...
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


% --- Executes just before slidingPLV_Aug2021 is made visible.
function slidingPLV_Aug2021_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to slidingPLV_Aug2021 (see VARARGIN)

% get the default parameter values for the sliders - (In my code I get
% frequency here.)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD THE DATA
% Load the data
% Motor img - calib
global answers stud subj blocks RestOnset MoveOnset

questions = { '\fontsize{12} Select Task: [1] Image Control [2] Real Control';...
                    '\fontsize{12} Which subject number would you like to analyze (12-18)?';...
                    '\fontsize{12} What frequency would you like to analyze first (Hz)?';...
                    '\fontsize{12} What frequency would like to analyze last (Hz)?';...
                    '\fontsize{12} How much should the frequency increment by?';...
                    '\fontsize{12} What should be the frequency width for the filter?';...
                    '\fontsize{12} What length of time do you want to analyze in each trial (sec)?';...
                    '\fontsize{12} What time do you want to start at (sec)?';...
                    '\fontsize{12} Do you want to arrange filtered data by time points [1] or trials [2]? (Always set 2)';...
                    '\fontsize{12} What will be the size of the arrangement? (Not used)';...
                    '\fontsize{12} Which day of sessions do you want to study from the subject (number 1-5)'};

definput = cell(11,1);
definput(:) = [{'1'};{'13'};{'14'};{'30'};{'0.5'};{'5'};{'7.5'};{'0.5'};{'2'};{'50'};{'2'}];
options.Interpreter = 'tex';
answers = inputdlg(questions, 'Parameters for PLV', [1 85], definput, options );

            
% Add as input for later
 % stud{1} -> MIcal; stud{2} -> MIcont
%  stud_type = {micon_imag_aug2021(); mrcon_real_aug2021() };
stud_type = {mical_imag_aug2021(); mrcal_real_aug2021() };
stud = stud_type{str2double(answers{1})};
subj = str2double(answers{2});
% day = 1;
% blocks = [load(stud{subj-11}{day}{1}(1,:)); load(stud{subj-11}{day}{1}(2,:))]; 

D = str2num(answers{11});

    if D > 0
        for Sess = 1: size(stud{subj-11}{D}{1},1)

            if (Sess==1)
                blocks = load([stud{subj-11}{D}{1}(Sess,:)]);
            else
                blocks = [blocks; load([stud{subj-11}{D}{1}(Sess,:)]) ];
            end

        end


    % if D==0, concatenate all blocks  of all days of the same subject
    elseif D == 0
        z = [];
        for D = 1: size(stud{subj-11}{1}{1},1)
            tmp = stud{subj-11}{D}{1};
            z = char(z,tmp);
        end
        z(1,:) = [];

        for Sess = 1:size(z,1)
            if Sess == 1
                blocks = [load([stud{subj-11}{Sess}{1}(1,:)]); load([stud{subj-11}{Sess}{1}(2,:)])];
            else
                blocks = [blocks; [load([stud{subj-11}{Sess}{1}(1,:)]); load([stud{subj-11}{Sess}{1}(2,:)])] ];
            end
        end
    else
        error('D cannot be negative')
    end



RestOnset = [blocks(1).RestOnset(:); blocks(2).RestOnset(:)+length(blocks(1).meg)];
MoveOnset = [blocks(1).MoveOnset(:); blocks(2).MoveOnset(:)+length(blocks(1).meg)];
        try
            blocks.MoveOffset;
            
            for i = 1:2
                % Make Equal Onset, Offset, & trialEnd indexes
                if length(blocks(i).MoveOnset) ~= length(blocks(i).MoveOffset)
                    blocks(i).trialEnd(length(blocks(i).trialEnd) + 1) = length(blocks(i).meg);                        
                    blocks(i).MoveOffset(length(blocks(i).MoveOffset) + 1) = length(blocks(i).meg);
                elseif length(blocks(i).RestOnset) ~= length(blocks(i).RestOffset)
                      blocks(i).trialEnd(length(blocks(i).trialEnd) + 1) = length(blocks(i).meg);   
                      blocks(i).RestOffset(length(blocks(i).RestOffset) + 1) = length(blocks(i).meg);
                end
            end
            
            % Combine RestOffset, MoveOffset, trialStart, trialEnd,
            % trialLength index time points for 
            RestOffset = [blocks(1).RestOffset(:); blocks(2).RestOffset(:)+length(blocks(1).meg)];
            MoveOffset = [blocks(1).MoveOffset(:); blocks(2).MoveOffset(:)+length(blocks(1).meg)];
            trialStart = [blocks(1).trialStart(:); blocks(2).trialStart(:)+length(blocks(1).meg)];
            trialEnd = [blocks(1).trialEnd(:); blocks(2).trialEnd(:)+length(blocks(1).meg)];       
            trialLength = cell(size(trialStart,1),2);
            
            for i = 1:size(trialStart,1)
                    trialLength{i,1} = trialEnd(i) - trialStart(i);  

                    if any(trialStart(i) == MoveOnset)
                       trialLength{i, 2}= 'MoveOnset';
                    elseif any(trialStart(i) == RestOnset)
                        trialLength{i, 2} = 'RestOnset';
                    end
            end
            % For visuals
%             trials = cat(2, num2cell(trialStart), num2cell(trialEnd), trialLength(:, 1), trialLength(:, 2));
        catch
            disp('There is no MoveOffset or RestOffset in these trial blocks');
            RestOffset = 0;
            MoveOffset = 0;
         end


              


% Delete unused channels
meg = [blocks(1).meg; blocks(2).meg];
meg = meg';
channels = size(meg,1);

labels = [];

for b = 1:size(blocks,1)
    labels = [labels; blocks(b).stimulusCode ];
end

temp_ind = find(labels ~= 0);

% Rearrange the channel numbers to go sequentially from top-left to bottom-right
% megMat = [];
% megMat(1,:) = meg(12,:);
% megMat(2,:) = meg(7,:);
% megMat(3,:) = meg(6,:);
% megMat(4,:) = meg(17,:);
% megMat(5,:) = meg(14,:);
% megMat(6,:) = meg(8,:);
% megMat(7,:) = meg(11,:);
% megMat(8,:) = meg(16,:);
% megMat(9,:) = meg(5,:);
% megMat(10,:) = meg(10,:);
% megMat(11,:) = meg(4,:);
% megMat(12,:) = meg(15,:);
% megMat(13,:) = meg(9,:);
% megMat(14,:) = meg(3,:);
% megMat(15,:) = meg(1,:);
% megMat(16,:) = meg(13,:);
% megMat(17,:) = meg(2,:);
% meg = megMat(:,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% NEED TO CUT OUT EXTRA DATA FOR MOTOR IMG CONTROL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS

global begF endF stepf srate

% Beginning frequency to analyze
begF = str2double(answers{3});
% End frequency to analyze
endF = str2double(answers{4});
% Increment steps by
stepf = str2double(answers{5});
fi = [stepf/(endF-begF) 1.0];
%FWHM for filterdat
freqwidt = str2double(answers{6});
fwMax = 10.0;
fwMin = 0.1;
fwi = [ 0.2/(fwMax-fwMin) 1.0];
% time window parameters (in seconds) - leave out 1st sec for each task
tws = str2double(answers{7});
timewin = cell(2,1);
timestart = str2double(answers{8});

timeMax = 8.0;
timeMin = 0.0;

% time step increment
tstepi = [ tws/(timeMax-timeMin) 1.0];
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

set(handles.slider_timewin_size, 'Max', 8.0);
set(handles.slider_timewin_size, 'Min', 0.0);
set(handles.slider_timewin_size, 'Value', tws);
set(handles.slider_timewin_size, 'SliderStep', tstepi);
timeWinSize = get(handles.slider_timewin_size, 'Value');
set(handles.text_timewinSize, 'String', ['Time Window Size = ' num2str(handles.slider_timewin_size.Value)]);

set(handles.slider_time, 'Max', 8.0);
set(handles.slider_time, 'Min', 0.0);
set(handles.slider_time, 'Value', timestart);
set(handles.slider_time, 'SliderStep', tstepi);
timeWin = get(handles.slider_time, 'Value');
set(handles.text_time, 'String', ['Starting Time = ' num2str(handles.slider_time.Value) ' Sec']);

%% 2 sec blank screen or 1 sec in beginning? 
%% If trial success, 2 or 1 sec pause at the end?
%% Lowest trial length: Rest 3901, Max Length: 13101
%% One Move trial at 12401
%% Current: Cut 1000 points before, 2100 points after
% Start with {Trials} x [ Channels x Time Points ]
% Rearrange data to be channels x time points x trials in for loop
global data megM3d megR3d roundnum

% if RestOffset ~=0
%     megR3d = cell(size(MoveOnset, 1), 1);
%     megM3d = cell(size(MoveOnset, 1), 1);
% 
%     filtdatM = [];
%     filtdatR = [];
%     % Cut out unrelated data in trials
%     % 1 second before, 2 seconds after
%     for i = 1:size(RestOnset, 1)
%         % Rest
%         megR3d{i} = meg( :, (RestOnset(i) + 1000) : RestOffset(i)-2000 );
%         % Move
%         megM3d{i} = meg( :, (MoveOnset(i) + 1000) : MoveOffset(i)-2000 );
%         
%         % filterFGx - Frequency domain narrow-band filter Gaussian
%         filtdatR =  [filtdatR filterFGx(megR3d{i}, srate, handles.slider_freq.Value,handles.slider_width.Value)];
%         filtdatM =  [filtdatM filterFGx(megM3d{i}, srate, handles.slider_freq.Value, handles.slider_width.Value)];
%     end
% 
%     % Make filtered move and rest data the same size, Round down from roundnum on
%     % the smallest matrix
%             roundnum = str2double(answers{10});
%             if size(filtdatR,2) <= size(filtdatM,2)
%                     filtdatR = filtdatR(: , 1:roundnum*(floor(length(filtdatR)/roundnum)));
%                     filtdatM = filtdatM(: , 1:roundnum*(floor(length(filtdatR)/roundnum)));
%             else
%                     filtdatM = filtdatM(: , 1:roundnum*(floor(length(filtdatM)/roundnum)));
%                     filtdatR = filtdatR(: , 1:roundnum*(floor(length(filtdatM)/roundnum)));
%             end
%             
%         if str2double(answers{9}) == 1
%              % Round by time points
%                filtdat = [reshape(filtdatR, size(filtdatR,1), roundnum, [])       reshape(filtdatM, size(filtdatM,1), roundnum, [])];
%         elseif str2double(answers{9}) == 2
%             % Round by trials
%               filtdat = [reshape(filtdatR, size(filtdatR,1), [] , roundnum)      reshape(filtdatM, size(filtdatR,1), [], roundnum)];      
%         end
    
% RestOffset == 0



%     sorted_onsets = sort([MoveOnset; RestOnset]);
%     largest_trial = max(diff(sorted_onsets));
%     
%     megR3d = zeros(channels, largest_trial, length(sorted_onsets) );
%     megM3d = zeros(channels, largest_trial, length(sorted_onsets) );
%     
%     
%     for i = 1:length(sorted_onsets)
%         if i == length(sorted_onsets)
%             if find(RestOnset == sorted_onsets(i))
%                 megR3d(:,:,i) = meg(:, sorted_onsets(i):length(meg) );
%             else find(MoveOnset == sorted_onsets(i))
%                 megM3d(:,:,i) = meg(:, sorted_onsets(i):length(meg) );
%             end
%             
%         else
%             if find(RestOnset == sorted_onsets(i));
%                 megR3d(:,:,i) = meg(:, sorted_onsets(i):(sorted_onsets(i+1)-1) );
%             else find(MoveOnset == sorted_onsets(i)) ;
%                 megM3d(:,:,i) = meg(:, sorted_onsets(i):(sorted_onsets(i+1)-1) ) ;
%             end
%         end
%         
%     end

% Initialize 3D Matrix of Rest Trial and Move Trial meg data
% Channels X Data X Trials
megR3d = zeros(channels, size(RestOnset(1): (MoveOnset(1)-4001) , 2) ,size(RestOnset,1));
    megM3d = zeros(channels, size(MoveOnset(1): RestOnset(2)-4001, 2) ,size(MoveOnset,1));
    data = zeros(channels, size(MoveOnset(1): RestOnset(2)-4001, 2)*2 ,size(RestOnset,1));
    

    for i = 1:size(RestOnset,1)
        % 8 second trials -> subtract out the last 4000 points
        megR3d(:,:,i) = meg(:,(RestOnset(i): RestOnset(i) + 7999 ));
        % Move trials -> subtract out the last 4000 points
        megM3d(:,:,i) = meg(:,(MoveOnset(i):MoveOnset(i)+7999 ) );
        data(:,:,i) = [megR3d(:,:,i) megM3d(:,:,i)];
    end
    % 62 Channels x 16000 Data Points X 30 Trials (First 8000 points are
    % Rest, 2nd 8000 pnts are Move
    filtdat = filterFGx(data, srate, handles.slider_freq.Value, handles.slider_width.Value);

        
    
%     megR3d = zeros(channels, size(RestOnset(1): (MoveOnset(2)-4001) , 2) ,size(RestOnset,1));
%     megM3d = zeros(channels, size(MoveOnset(1): RestOnset(2)-4001, 2) ,size(MoveOnset,1));
%     data = zeros(channels, size(MoveOnset(1): RestOnset(2)-4001, 2)*2 ,size(RestOnset,1));
%     
% 
% 
%     % motor img - calib
%     for i = 1:size(RestOnset,1)
%         % 8 second trials -> subtract out the last 4000 points
%         megR3d(:,:,i) = meg(:,(RestOnset(i): (MoveOnset(i)-4001) ));
%         % Move trials -> subtract out the last 4000 points
%         megM3d(:,:,i) = meg(:,(MoveOnset(i): (MoveOnset(i) + (size(MoveOnset(1): RestOnset(2)-4002, 2)) ) ) );
%         data(:,:,i) = [megR3d(:,:,i) megM3d(:,:,i)];
%     end
%     filtdat = filterFGx(data, srate, handles.slider_freq.Value, handles.slider_width.Value);
% end


% time vector
timevec = 0+1/srate:(1/(srate)): size(filtdat,2)/srate;

% Rest Time Window - default 0 to 0.5 sec
timewin{1} = [timeWin  timeWin+timeWinSize];
% Move Time Window - default 8 to 8.5 sec
timewin{2} = [timeWin+((size(filtdat,2)/2)/srate)    timeWin+timeWinSize+((size(filtdat,2)/2)/srate)];

% convert to time indices
% Rest time indices
tidx1 = dsearchn(timevec',timewin{1}');
% Move time indices
tidx2 = dsearchn(timevec',timewin{2}');

synchmat = zeros(2,channels,channels);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE INITIAL PLOTS HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%
global pRx pRy pMx pMy pDx pDy pSx pSy g hori vert SS_chan_inter

% % REST SYNCHRONIZATION INITIAL PLOT
%         handles.restTW.XLim = [ 0.5 17.5 ];
%         handles.restTW.YLim = [ 0.5 17.5 ];
%         handles.restTW.CLimMode = 'manual'
%         handles.restTW.CLim = [0 0.9];
%         colorbar(handles.restTW)
%         handles.restTW.Colormap = lbmap(500, 'RedBlue');
%         % Make the axes square
%         handles.restTW.PlotBoxAspectRatio = [1 1 1];
%         % drawGrid(handles.restTW, 0.75);
%         handles.plotR = imagesc( handles.restTW, 'CData', squeeze(synchmat(1,:,:)));
%         handles.plotR.Parent.Title.String= [ 'Rest Synch: ' num2str(timewin{1}(1)) '-' num2str(timewin{1}(2)) ' Sec, '  '? \pm ? Hz'];
%         handles.plotR.Parent.Title.Position = [10.5 17.85 0];
%         handles.plotR.Parent.CLim = [0 0.9];
%         % Plot grid lines 
%         hold (handles.plotR.Parent,'on');
%             for j = 1:channels
%                 pRy{j} = plot(handles.plotR.Parent, [0.5, 17.5],[j-0.5, j-0.5],'k-');
%                 pRx{j} = plot(handles.plotR.Parent,[j-0.5, j-0.5],[0.5, 17.5],'k-');
%             end
%         hold (handles.plotR.Parent,'off')
% 
% 
% % MOVEMENT SYNCHRONIZATION INITIAL PLOT
%         handles.moveTW.XLim = [0.5 17.5];
%         handles.moveTW.YLim = [0.5 17.5];
%         handles.moveTW.CLimMode = 'manual'
%         handles.moveTW.CLim = [0 0.9];
%         colorbar(handles.moveTW)
%         handles.moveTW.Colormap = lbmap(500, 'RedBlue');
%         % Make the axes square
%         handles.moveTW.PlotBoxAspectRatio = [1 1 1];
%         handles.plotM = imagesc( handles.moveTW, 'CData', squeeze(synchmat(2,:,:)));
%         handles.plotM.Parent.Title.String = [ 'Move Synch: ' num2str(timewin{2}(1)) '-' num2str(timewin{2}(2)) ' Sec, '  '? \pm ? Hz'];
%         handles.plotM.Parent.Title.Position = [11.5 17.85 0];
%         handles.plotM.Parent.CLim = [0 0.9];
%         % Plot grid lines 
%         hold (handles.plotM.Parent,'on')
%             for j = 1:channels
%                 pMy{j} = plot(handles.moveTW, [0.5,18.5],[j-.5,j-.5],'k-');
%                 pMx{j} = plot(handles.moveTW,[j-.5,j-.5],[0.5, 18.5],'k-');
%             end
%         hold (handles.plotM.Parent,'off')
% 
% 
% % PHASE DIFFERENCE SYNCHRONIZATION PLOT
%         handles.phaseDif.XLim = [0.5 17.5];
%         handles.phaseDif.YLim = [0.5 17.5];
%         handles.phaseDif.CLimMode = 'manual';
%         colorbar(handles.phaseDif)
%         handles.phaseDif.Colormap = lbmap(500, 'RedBlue');
%         handles.phaseDif.CLim = [-.1 .1];
%         % Make the axes square
%         handles.phaseDif.PlotBoxAspectRatio = [1 1 1];
%         handles.plotD = imagesc( handles.phaseDif, 'CData', squeeze(diff(synchmat)));
%         handles.plotD.Parent.Title.String = 'Phase Synchronization Difference: Move - Rest';
%         handles.plotD.Parent.Title.Position = [10.5000 17.85 0];
%         handles.plotR.Parent.CLim = [-.1 .1];        
%         hold (handles.plotD.Parent,'on');
%             for j = 1:channels
%                 pDy{j} = plot(handles.phaseDif, [0.5,18.5],[j-.5,j-.5],'k-');
%                 pDx{j} = plot(handles.phaseDif,[j-.5,j-.5],[0.5, 18.5],'k-');
%             end
%         hold (handles.plotD.Parent,'off');

% KRUSKAL WALLIS STASTICALLY SIGNIFICANT PHASE SYNCHRONIZATION
        % GUI Limits, colors
        handles.krusWal.XLim = [0.5 62.5];
        handles.krusWal.YLim = [0.5 62.5];
        handles.krusWal.CLimMode = 'manual';
        colorbar(handles.krusWal)
        handles.krusWal.Colormap = lbmap(500, 'RedBlue');
        handles.krusWal.CLim = [0    0.01];
        % Mark Statistically Significant Phase Synchs in Kruskal-Wallis Test
        hori = repmat(1:size(filtdat,1),size(filtdat,1),1);
        vert = hori';
        g = text(hori(:), vert(:), ' ', 'HorizontalAlignment', 'Center', 'Color','k', 'FontWeight', 'bold', 'Parent', handles.krusWal);
        % Make the axes square
        handles.krusWal.PlotBoxAspectRatio = [1 1 1];
        handles.plotSS = imagesc( handles.krusWal, 'CData', squeeze(diff(synchmat)));
        handles.plotSS.Parent.Title.String = 'Kruskal-Wallis Stat. Sig. \alpha = 0.05/ 62 Channels';
        handles.plotSS.Parent.Title.Position = [30 62.85 0];
        handles.plotSS.Parent.CLim = [0 0.01];        
        hold (handles.plotSS.Parent,'on');
            for j = 1:channels
                pSy{j} = plot(handles.krusWal, [.5, 63.5],[j-.5,j-.5],'k-');
                pSx{j} = plot(handles.krusWal,[j-.5,j-.5],[.5,63.5],'k-');
            end
        hold (handles.plotSS.Parent,'off');

    % angle time-series
    angts = zeros(size(filtdat));
    
% (In my code I want to compute Phase-Locking Value here)
    %     How do I take out this for loop?
    for triali=1:size(filtdat,3)
        % Applying hilbert transform to matrix inputs, phase angles should be
        % sawtooth slanted to the right
        angts(:,:,triali) = angle(hilbert(squeeze(filtdat(:,:,triali))').');
        % Plot angle time series and check phase angles
        % figure(5); plot(angts(5,:,10));
    end
    
    synchmat1 = zeros( (size(filtdat,3)), (size(filtdat,1)), (size(filtdat,1)) );
    synchmat2 = zeros( (size(filtdat,3)), (size(filtdat,1)), (size(filtdat,1)) );

    % Is it possible to take out these for loops?
    nchans = size(filtdat,1);
     for chani=1:nchans
        for chanj=1:nchans
            for tria = 1:size(filtdat,3)

                %%% time window 1
                % extract angles in channel i, all trials
                tmpAi = angts(chani,tidx1(1):tidx1(2),tria);
                % extract angles in channel j all trials
                tmpAj = angts(chanj,tidx1(1):tidx1(2),tria);
                % compute synchronization between the two channels
                % abs( mean( eulers formula( phase angle diff)), 2nd dim: time points)
                trialsynch = abs(mean(exp(1i*( tmpAi-tmpAj )),2));
                % average over trials - synchmat 1st Time Window
                synchmat1(tria,chani,chanj) = mean(trialsynch);

                %%% time window 2
                % extract angles
                tmpAi = angts(chani,tidx2(1):tidx2(2),tria);
                tmpAj = angts(chanj,tidx2(1):tidx2(2),tria);
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
    diffS = squeeze(abs(synchmatD(2,:,:) - synchmatD(1,:,:)));
%     diffS == diff(synchmatD)
%     subplot(1,4,1)
%     title([ 'Synch: ' num2str(timewin{1}(1)) '-' num2str(timewin{1}(2)) ' Sec, ' num2str(centfreq) ' Hz \pm ' num2str(freqwidt)]);

%% UPDATE THE IMAGES
% % Update Rest Trials Plot
%         handles.plotR = imagesc( handles.restTW, 'CData', squeeze(mean(synchmat1,1))' );
%         handles.plotR.Parent.Title.String =  [ 'Rest Phase Synch: ' num2str(timewin{1}(1)) '-' num2str(timewin{1}(2)) ' Sec, ' num2str(centfreq) ' \pm ' num2str(freqwidt) ' Hz'];
%         handles.plotR.Parent.Title.FontSize = 11;
%         handles.plotR.Parent.CLim = [0 0.9];
%         [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30, p31, p32, p33, p34, p35, p36, p37, p38, p39, p40, p41, p42, p43, p44, p45, p46, p47, p48, p49, p50, p51, p52, p53, p54, p55, p56, p57, p58, p59, p60, p61, p62, p63, p64, p65, p66, p67, p68, p69, p70, p71, p72, p73, p74, p75, p76, p77, p78, p79, p80, p81, p82, p83, p84, p85, p86, p87, p88, p89, p90, p91, p92, p93, p94, p95, p96] = deal(pRx{:});
%         uistack([p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30, p31, p32, p33, p34, p35, p36, p37, p38, p39, p40, p41, p42, p43, p44, p45, p46, p47, p48, p49, p50, p51, p52, p53, p54, p55, p56, p57, p58, p59, p60, p61, p62, p63, p64, p65, p66, p67, p68, p69, p70, p71, p72, p73, p74, p75, p76, p77, p78, p79, p80, p81, p82, p83, p84, p85, p86, p87, p88, p89, p90, p91, p92, p93, p94, p95, p96], 'top');
%         [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30, p31, p32, p33, p34, p35, p36, p37, p38, p39, p40, p41, p42, p43, p44, p45, p46, p47, p48, p49, p50, p51, p52, p53, p54, p55, p56, p57, p58, p59, p60, p61, p62, p63, p64, p65, p66, p67, p68, p69, p70, p71, p72, p73, p74, p75, p76, p77, p78, p79, p80, p81, p82, p83, p84, p85, p86, p87, p88, p89, p90, p91, p92, p93, p94, p95, p96] = deal(pRy{:});
%         uistack([p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30, p31, p32, p33, p34, p35, p36, p37, p38, p39, p40, p41, p42, p43, p44, p45, p46, p47, p48, p49, p50, p51, p52, p53, p54, p55, p56, p57, p58, p59, p60, p61, p62, p63, p64, p65, p66, p67, p68, p69, p70, p71, p72, p73, p74, p75, p76, p77, p78, p79, p80, p81, p82, p83, p84, p85, p86, p87, p88, p89, p90, p91, p92, p93, p94, p95, p96], 'top');
% 
% % Update Move Trials Plot
%         handles.plotM = imagesc( handles.moveTW, 'CData', squeeze(mean(synchmat2,1)) );
%         handles.plotM.Parent.Title.String =  [ 'Move Phase Synch: ' num2str(timewin{2}(1)) '-' num2str(timewin{2}(2)) ' Sec, ' num2str(centfreq) ' Hz \pm ' num2str(freqwidt) ' Hz'];
%         handles.plotM.Parent.Title.FontSize = 11;
%         handles.plotM.Parent.CLim = [0 0.9];
%         [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30, p31, p32, p33, p34, p35, p36, p37, p38, p39, p40, p41, p42, p43, p44, p45, p46, p47, p48, p49, p50, p51, p52, p53, p54, p55, p56, p57, p58, p59, p60, p61, p62, p63, p64, p65, p66, p67, p68, p69, p70, p71, p72, p73, p74, p75, p76, p77, p78, p79, p80, p81, p82, p83, p84, p85, p86, p87, p88, p89, p90, p91, p92, p93, p94, p95, p96] = deal(pMx{:});
%         uistack([p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30, p31, p32, p33, p34, p35, p36, p37, p38, p39, p40, p41, p42, p43, p44, p45, p46, p47, p48, p49, p50, p51, p52, p53, p54, p55, p56, p57, p58, p59, p60, p61, p62, p63, p64, p65, p66, p67, p68, p69, p70, p71, p72, p73, p74, p75, p76, p77, p78, p79, p80, p81, p82, p83, p84, p85, p86, p87, p88, p89, p90, p91, p92, p93, p94, p95, p96], 'top');
%         [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30, p31, p32, p33, p34, p35, p36, p37, p38, p39, p40, p41, p42, p43, p44, p45, p46, p47, p48, p49, p50, p51, p52, p53, p54, p55, p56, p57, p58, p59, p60, p61, p62, p63, p64, p65, p66, p67, p68, p69, p70, p71, p72, p73, p74, p75, p76, p77, p78, p79, p80, p81, p82, p83, p84, p85, p86, p87, p88, p89, p90, p91, p92, p93, p94, p95, p96] = deal(pMy{:});
%         uistack([p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30, p31, p32, p33, p34, p35, p36, p37, p38, p39, p40, p41, p42, p43, p44, p45, p46, p47, p48, p49, p50, p51, p52, p53, p54, p55, p56, p57, p58, p59, p60, p61, p62, p63, p64, p65, p66, p67, p68, p69, p70, p71, p72, p73, p74, p75, p76, p77, p78, p79, p80, p81, p82, p83, p84, p85, p86, p87, p88, p89, p90, p91, p92, p93, p94, p95, p96], 'top');
%   
% % Update Phase Difference Plot
%         handles.plotD = imagesc( handles.phaseDif, 'CData', squeeze(diff(synchmatD)) ); 
% %         handles.plotD.Parent.Title.String = ['Phase Synchronization Difference: Move - Rest'];
%         handles.plotD.Parent.Title.FontSize = 13.5;
%         handles.plotD.Parent.CLim = [-0.1 0.1];
%         [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30, p31, p32, p33, p34, p35, p36, p37, p38, p39, p40, p41, p42, p43, p44, p45, p46, p47, p48, p49, p50, p51, p52, p53, p54, p55, p56, p57, p58, p59, p60, p61, p62, p63, p64, p65, p66, p67, p68, p69, p70, p71, p72, p73, p74, p75, p76, p77, p78, p79, p80, p81, p82, p83, p84, p85, p86, p87, p88, p89, p90, p91, p92, p93, p94, p95, p96] = deal(pDx{:});
%         uistack([p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30, p31, p32, p33, p34, p35, p36, p37, p38, p39, p40, p41, p42, p43, p44, p45, p46, p47, p48, p49, p50, p51, p52, p53, p54, p55, p56, p57, p58, p59, p60, p61, p62, p63, p64, p65, p66, p67, p68, p69, p70, p71, p72, p73, p74, p75, p76, p77, p78, p79, p80, p81, p82, p83, p84, p85, p86, p87, p88, p89, p90, p91, p92, p93, p94, p95, p96], 'top');
%         [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30, p31, p32, p33, p34, p35, p36, p37, p38, p39, p40, p41, p42, p43, p44, p45, p46, p47, p48, p49, p50, p51, p52, p53, p54, p55, p56, p57, p58, p59, p60, p61, p62, p63, p64, p65, p66, p67, p68, p69, p70, p71, p72, p73, p74, p75, p76, p77, p78, p79, p80, p81, p82, p83, p84, p85, p86, p87, p88, p89, p90, p91, p92, p93, p94, p95, p96] = deal(pDy{:});
%         uistack([p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30, p31, p32, p33, p34, p35, p36, p37, p38, p39, p40, p41, p42, p43, p44, p45, p46, p47, p48, p49, p50, p51, p52, p53, p54, p55, p56, p57, p58, p59, p60, p61, p62, p63, p64, p65, p66, p67, p68, p69, p70, p71, p72, p73, p74, p75, p76, p77, p78, p79, p80, p81, p82, p83, p84, p85, p86, p87, p88, p89, p90, p91, p92, p93, p94, p95, p96], 'top');
    
    % Kruskal-Wallis test on Phase-Locking Value data
    krusP = zeros(nchans,nchans);
    for ii = 1:nchans
        for jj = 1:nchans
            % Kruskal-Wallis compares one against many trials and gives
            % significance
            [krusP(ii,jj)] = kruskalwallis([synchmat1(:,ii,jj) synchmat2(:,ii,jj)], [], 'off');
        end
    end
    % Arrange kruskal-wallis values for labeling in image
   chans_krus_wall_vals  = num2cell(krusP);
    chans_krus_wall_vals = cellfun(@num2str, chans_krus_wall_vals, 'UniformOutput', false); % convert to string  

%     ss = krusP <= alphaT;
%     ssi = find(ss);
    
% Update Kruskal-Wallis Statistical Test Plot
        handles.plotSS = imagesc( handles.krusWal, 'CData', krusP );
        handles.plotSS.Parent.Title.FontSize = 13.5;
        handles.plotSS.Parent.CLim = [0 0.01];
        handles.krusWal.CLim = [0 0.01];
        [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30, p31, p32, p33, p34, p35, p36, p37, p38, p39, p40, p41, p42, p43, p44, p45, p46, p47, p48, p49, p50, p51, p52, p53, p54, p55, p56, p57, p58, p59, p60, p61, p62] = deal(pSx{:});
        uistack([p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30, p31, p32, p33, p34, p35, p36, p37, p38, p39, p40, p41, p42, p43, p44, p45, p46, p47, p48, p49, p50, p51, p52, p53, p54, p55, p56, p57, p58, p59, p60, p61, p62], 'top');
        [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30, p31, p32, p33, p34, p35, p36, p37, p38, p39, p40, p41, p42, p43, p44, p45, p46, p47, p48, p49, p50, p51, p52, p53, p54, p55, p56, p57, p58, p59, p60, p61, p62] = deal(pSy{:});
        uistack([p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30, p31, p32, p33, p34, p35, p36, p37, p38, p39, p40, p41, p42, p43, p44, p45, p46, p47, p48, p49, p50, p51, p52, p53, p54, p55, p56, p57, p58, p59, p60, p61, p62], 'top');
    
        
        handles.g = g;
        
    for tt = 1:length(g)
        % Reset previous stat sig labels
        g(tt).String = ' ';
        alphaT = 0.05/size(filtdat,1);
        if str2double(chans_krus_wall_vals{tt})<= alphaT
            % Label images that are stat sig
             g(tt) = text(hori(tt), vert(tt), 'X', 'HorizontalAlignment', 'Center', 'Color','w', 'FontWeight', 'bold', 'FontSize', 11, 'Parent', handles.krusWal);
        else
            continue;
        end
    end

ind = tril(zeros(62,62));
for i = 1:length(g)
    if g(i).String == 'X'
        ind(i) = 1;
    end
end

handles.g = g;

[rowX,colX] = find(tril(ind,-1) == 1);
SS_chan_inter = [rowX colX];

% Save Statistically Significant channel interactions to workspace and take
% out the duplicates
assignin('base', 'SS_chan_inter', SS_chan_inter)
assignin('base', 'krusP', krusP)
assignin('base', 'chans_krus_wall_vals', chans_krus_wall_vals)
assignin('base', 'answers', answers)

% Choose default command line output for slidingPLV_Aug2021
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

'done'



% UIWAIT makes slidingPLV_Aug2021 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = slidingPLV_Aug2021_OutputFcn(hObject, eventdata, handles) 
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

global answers stud subj blocks RestOnset MoveOnset nchans krusP chans_krus_wall_vals alphaT RestOffset
global begF endF stepf meg
global pRx pRy pMx pMy pDx pDy pSx pSy g hori vert
global data srate megM3d megR3d roundnum

    % get the new parameter values
    centfreq = get(handles.slider_freq, 'Value');
    set(handles.text_centfreq, 'String', ['Center Frequency = ' num2str(handles.slider_freq.Value) ' Hz']);
    
    freqwidt = get(handles.slider_width, 'Value');
    set(handles.text_freqwidt, 'String', ['Frequency Width = ' num2str(handles.slider_width.Value)]);


    tws = get(handles.slider_timewin_size, 'Value');
    set(handles.text_timewinSize, 'String', ['Time Window Size = ' num2str(handles.slider_timewin_size.Value)]);

    timestart = get(handles.slider_time, 'Value');
    set(handles.text_time, 'String', ['Starting Time = ' num2str(handles.slider_time.Value) ' Sec']);

%%% CHANGE FOR LATER! +5 MUST REPRESENT CONTROL AS WELL
if RestOffset ~=0

        filtdatM = [];
        filtdatR = [];
        filtdat = [];
    for i = 1:size(megR3d,1)
        filtdatR =  [filtdatR,  filterFGx(megR3d{i}, srate, handles.slider_freq.Value,handles.slider_width.Value)];
        filtdatM =  [filtdatM,  filterFGx(megM3d{i}, srate, handles.slider_freq.Value, handles.slider_width.Value)];
    end
        % Filter data second to filter out edge artifacts
        % Make filtered move and rest data the same size, Round down from 50 on
        % the smallest matrix
        
            if size(filtdatR,2) <= size(filtdatM,2)
                    filtdatR = filtdatR(: , 1:roundnum*(floor(length(filtdatR)/roundnum)));
                    filtdatM = filtdatM(: , 1:roundnum*(floor(length(filtdatR)/roundnum)));
            else
                    filtdatM = filtdatM(: , 1:roundnum*(floor(length(filtdatM)/roundnum)));
                    filtdatR = filtdatR(: , 1:roundnum*(floor(length(filtdatM)/roundnum)));
            end
            
        if str2double(answers{9}) == 1
            % Round by time points
            filtdat = [reshape(filtdatR, nchans, roundnum, []) reshape(filtdatM, nchans, roundnum, [])];
        elseif str2double(answers{9}) == 2
            % Round by trials
            filtdat = [reshape(filtdatR, nchans, [] , roundnum) reshape(filtdatM, nchans, [], roundnum)];
        end
        

    
else
    
    % motor img - calib
    for ii = 1:size(RestOnset,1)
        data(:,:,ii) = [megR3d(:,:,ii) megM3d(:,:,ii)];
    end
    filtdat = filterFGx(data, srate, handles.slider_freq.Value, handles.slider_width.Value);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FOCUS ON HERE FOR CHECKING START TIME (timestart)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% time vector
timevec = 0+1/srate:(1/(srate)):size(filtdat,2)/srate;

%timewin{1} only checks the starting time (default 0 sec) to the Time Window Size (0.5 sec)
%For Rest PLV
timewin{1} = [timestart  timestart+tws];

%timewin{2} for Move PLV, at the end of the filtdat
timewin{2} = [timestart+(size(filtdat,2)/2)/srate    timestart+tws+(size(filtdat,2)/2)/srate];

% convert to time indices (sec)
% 
tidx1 = dsearchn(timevec',timewin{1}');
tidx2 = dsearchn(timevec',timewin{2}');

    synchmat1 = zeros( (size(filtdat,3)), (size(filtdat,1)), (size(filtdat,1)) );
    synchmat2 = zeros( (size(filtdat,3)), (size(filtdat,1)), (size(filtdat,1)) );

    nchans = size(data,1);
       for triali=1:size(data,3)
            % Applying hilbert transform to matrix inputs
            % This is where frequency-related data is put in, from filtdat
            angts(:,:,triali) = angle(hilbert(squeeze(filtdat(:,:,triali))').');
       end

         for chani=1:nchans
            for chanj=1:nchans
                for tria = 1:size(data,3)

                    %%% time window 1
                    % synchmat1 - Rest PLV Plot
                 synchmat1(tria,chani,chanj) = mean( abs(mean(exp(1i* (angts(chani,tidx1(1):tidx1(2),tria) - angts(chanj,tidx1(1):tidx1(2),tria))),2)) );
%                     % extract angles in channel i, all trials
%                     tmpAi = angts(chani,tidx1(1):tidx1(2),tria);
%                     % extract angles in channel j all trials
%                     tmpAj = angts(chanj,tidx1(1):tidx1(2),tria);
%                     % compute synchronization between the two channels
%                     % abs( mean( eulers formula( phase angle diff)), 2nd dim: time points)
%                     trialsynch = abs(mean(exp(1i*( tmpAi-tmpAj )),2));
%                     % average over trials - synchmat 1st Time Window
%                     synchmat1(tria,chani,chanj) = mean(trialsynch);

                    %%% time window 2
                    % synchmat2 - Move PLV Plot
                    % extract phase angles
                    tmpAi = angts(chani,tidx2(1):tidx2(2),tria);
                    tmpAj = angts(chanj,tidx2(1):tidx2(2),tria);
                    % Phase Locking Value on each trial
                    trialsynch = abs(mean(exp(1i*( tmpAi-tmpAj )),2));
                    % average over trials - synchmat 2nd Time Window

                    synchmat2(tria,chani,chanj) = mean(trialsynch); 
                end            
            end
         end
        % synchmat Difference - Time Window x channels x channels - These
        % must be the Rest Trials 
        synchmatD = [mean(synchmat1,1); mean(synchmat2,1)];
    %     synchmatDif(1,chani,chanj) = mean((synchmat2(:,chani,chanj) - synchmat1(:,chani,chanj)),1);
        % diffS = synchmatD(2,:,:) - synchmatD(1,:,:);
        % diffS == diff(synchmatD)
    %     subplot(1,4,1)
    %     title([ 'Synch: ' num2str(timewin{1}(1)) '-' num2str(timewin{1}(2)) ' Sec, ' num2str(centfreq) ' Hz \pm ' num2str(freqwidt)]);

    %% UPDATE THE IMAGES
%     % Update Rest Trials Plot
% %             set(handles.plotR, 'CData', squeeze(mean(synchmat1,1))' );
%             handles.restTW = imagesc( handles.restTW, 'CData', squeeze(mean(synchmat1,1))' ) 
%             handles.restTW.Parent.Title.String =  [ 'Rest Phase Synch: ' num2str(timewin{1}(1)) '-' num2str(timewin{1}(2)) ' Sec, ' num2str(centfreq) ' \pm ' num2str(freqwidt) ' Hz'];
%             handles.restTW.Parent.Title.FontSize = 11;
%         [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30, p31, p32, p33, p34, p35, p36, p37, p38, p39, p40, p41, p42, p43, p44, p45, p46, p47, p48, p49, p50, p51, p52, p53, p54, p55, p56, p57, p58, p59, p60, p61, p62, p63, p64, p65, p66, p67, p68, p69, p70, p71, p72, p73, p74, p75, p76, p77, p78, p79, p80, p81, p82, p83, p84, p85, p86, p87, p88, p89, p90, p91, p92, p93, p94, p95, p96] = deal(pRx{:});
%         uistack([p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30, p31, p32, p33, p34, p35, p36, p37, p38, p39, p40, p41, p42, p43, p44, p45, p46, p47, p48, p49, p50, p51, p52, p53, p54, p55, p56, p57, p58, p59, p60, p61, p62, p63, p64, p65, p66, p67, p68, p69, p70, p71, p72, p73, p74, p75, p76, p77, p78, p79, p80, p81, p82, p83, p84, p85, p86, p87, p88, p89, p90, p91, p92, p93, p94, p95, p96], 'top');
%         [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30, p31, p32, p33, p34, p35, p36, p37, p38, p39, p40, p41, p42, p43, p44, p45, p46, p47, p48, p49, p50, p51, p52, p53, p54, p55, p56, p57, p58, p59, p60, p61, p62, p63, p64, p65, p66, p67, p68, p69, p70, p71, p72, p73, p74, p75, p76, p77, p78, p79, p80, p81, p82, p83, p84, p85, p86, p87, p88, p89, p90, p91, p92, p93, p94, p95, p96] = deal(pRy{:});
%         uistack([p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30, p31, p32, p33, p34, p35, p36, p37, p38, p39, p40, p41, p42, p43, p44, p45, p46, p47, p48, p49, p50, p51, p52, p53, p54, p55, p56, p57, p58, p59, p60, p61, p62, p63, p64, p65, p66, p67, p68, p69, p70, p71, p72, p73, p74, p75, p76, p77, p78, p79, p80, p81, p82, p83, p84, p85, p86, p87, p88, p89, p90, p91, p92, p93, p94, p95, p96], 'top');
% 
% 
%     % Update Move Trials Plot
% %              set(handles.plotM, 'CData', squeeze(mean(synchmat2,1))' );
% 
%             handles.plotM = imagesc( handles.moveTW, 'CData', squeeze(mean(synchmat2,1))' )
%             handles.plotM.Parent.Title.String =  [ 'Move Phase Synch: ' num2str(timewin{2}(1)) '-' num2str(timewin{2}(2)) ' Sec, ' num2str(centfreq) ' Hz \pm ' num2str(freqwidt) ' Hz'];
%             handles.plotM.Parent.Title.FontSize = 11;
%         [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30, p31, p32, p33, p34, p35, p36, p37, p38, p39, p40, p41, p42, p43, p44, p45, p46, p47, p48, p49, p50, p51, p52, p53, p54, p55, p56, p57, p58, p59, p60, p61, p62, p63, p64, p65, p66, p67, p68, p69, p70, p71, p72, p73, p74, p75, p76, p77, p78, p79, p80, p81, p82, p83, p84, p85, p86, p87, p88, p89, p90, p91, p92, p93, p94, p95, p96] = deal(pMx{:});
%         uistack([p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30, p31, p32, p33, p34, p35, p36, p37, p38, p39, p40, p41, p42, p43, p44, p45, p46, p47, p48, p49, p50, p51, p52, p53, p54, p55, p56, p57, p58, p59, p60, p61, p62, p63, p64, p65, p66, p67, p68, p69, p70, p71, p72, p73, p74, p75, p76, p77, p78, p79, p80, p81, p82, p83, p84, p85, p86, p87, p88, p89, p90, p91, p92, p93, p94, p95, p96], 'top');
%         [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30, p31, p32, p33, p34, p35, p36, p37, p38, p39, p40, p41, p42, p43, p44, p45, p46, p47, p48, p49, p50, p51, p52, p53, p54, p55, p56, p57, p58, p59, p60, p61, p62, p63, p64, p65, p66, p67, p68, p69, p70, p71, p72, p73, p74, p75, p76, p77, p78, p79, p80, p81, p82, p83, p84, p85, p86, p87, p88, p89, p90, p91, p92, p93, p94, p95, p96] = deal(pMy{:});
%         uistack([p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30, p31, p32, p33, p34, p35, p36, p37, p38, p39, p40, p41, p42, p43, p44, p45, p46, p47, p48, p49, p50, p51, p52, p53, p54, p55, p56, p57, p58, p59, p60, p61, p62, p63, p64, p65, p66, p67, p68, p69, p70, p71, p72, p73, p74, p75, p76, p77, p78, p79, p80, p81, p82, p83, p84, p85, p86, p87, p88, p89, p90, p91, p92, p93, p94, p95, p96], 'top');
% 
% 
%     % Update Phase Difference Plot
% %             set(handles.plotD, 'CData', squeeze(diff(synchmatD)) );
%             handles.plotD = imagesc( handles.phaseDif, 'CData', squeeze(diff(synchmatD)) ) 
%     %         handles.plotD.Parent.Title.String = ['Phase Synchronization Difference: Move - Rest'];
%             handles.plotD.Parent.Title.FontSize = 13.5;
%         [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30, p31, p32, p33, p34, p35, p36, p37, p38, p39, p40, p41, p42, p43, p44, p45, p46, p47, p48, p49, p50, p51, p52, p53, p54, p55, p56, p57, p58, p59, p60, p61, p62, p63, p64, p65, p66, p67, p68, p69, p70, p71, p72, p73, p74, p75, p76, p77, p78, p79, p80, p81, p82, p83, p84, p85, p86, p87, p88, p89, p90, p91, p92, p93, p94, p95, p96] = deal(pDx{:});
%         uistack([p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30, p31, p32, p33, p34, p35, p36, p37, p38, p39, p40, p41, p42, p43, p44, p45, p46, p47, p48, p49, p50, p51, p52, p53, p54, p55, p56, p57, p58, p59, p60, p61, p62, p63, p64, p65, p66, p67, p68, p69, p70, p71, p72, p73, p74, p75, p76, p77, p78, p79, p80, p81, p82, p83, p84, p85, p86, p87, p88, p89, p90, p91, p92, p93, p94, p95, p96], 'top');
%         [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30, p31, p32, p33, p34, p35, p36, p37, p38, p39, p40, p41, p42, p43, p44, p45, p46, p47, p48, p49, p50, p51, p52, p53, p54, p55, p56, p57, p58, p59, p60, p61, p62, p63, p64, p65, p66, p67, p68, p69, p70, p71, p72, p73, p74, p75, p76, p77, p78, p79, p80, p81, p82, p83, p84, p85, p86, p87, p88, p89, p90, p91, p92, p93, p94, p95, p96] = deal(pDy{:});
%         uistack([p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30, p31, p32, p33, p34, p35, p36, p37, p38, p39, p40, p41, p42, p43, p44, p45, p46, p47, p48, p49, p50, p51, p52, p53, p54, p55, p56, p57, p58, p59, p60, p61, p62, p63, p64, p65, p66, p67, p68, p69, p70, p71, p72, p73, p74, p75, p76, p77, p78, p79, p80, p81, p82, p83, p84, p85, p86, p87, p88, p89, p90, p91, p92, p93, p94, p95, p96], 'top');
% 

        % Kruskal-Wallis test on Phase-Locking Value data
        krusP = zeros(nchans,nchans);
        for ii = 1:nchans
            for jj = 1:nchans
                %KRUSKALWALLIS Nonparametric one-way analysis of variance (ANOVA).
                [krusP(ii,jj), tbl{ii,jj}, stats{ii,jj}] = kruskalwallis([synchmat1(:,ii,jj) synchmat2(:,ii,jj)], [], 'off');
            end
        end
        % Arrange kruskal-wallis values for labeling in image
        chans_krus_wall_vals = num2cell(krusP);
        chans_krus_wall_vals = cellfun(@num2str, chans_krus_wall_vals, 'UniformOutput', false); % convert to string  

    %     ss = krusP <= alphaT;
    %     ssi = find(ss);
    
    % Update Kruskal-Wallis Statistical Test Plot
%         set(handles.plotSS, 'CData', krusP );
            handles.plotSS = imagesc( handles.krusWal, 'CData', krusP )
            handles.plotSS.Parent.Title.FontSize = 13.5;
            
        handles.krusWal.CLimMode = 'manual';
        colorbar(handles.krusWal)
        handles.krusWal.Colormap = lbmap(500, 'RedBlue');
        handles.krusWal.CLim = [0  0.01];
        [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30, p31, p32, p33, p34, p35, p36, p37, p38, p39, p40, p41, p42, p43, p44, p45, p46, p47, p48, p49, p50, p51, p52, p53, p54, p55, p56, p57, p58, p59, p60, p61, p62] = deal(pSx{:});
        uistack([p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30, p31, p32, p33, p34, p35, p36, p37, p38, p39, p40, p41, p42, p43, p44, p45, p46, p47, p48, p49, p50, p51, p52, p53, p54, p55, p56, p57, p58, p59, p60, p61, p62], 'top');
        [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30, p31, p32, p33, p34, p35, p36, p37, p38, p39, p40, p41, p42, p43, p44, p45, p46, p47, p48, p49, p50, p51, p52, p53, p54, p55, p56, p57, p58, p59, p60, p61, p62] = deal(pSy{:});
        uistack([p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30, p31, p32, p33, p34, p35, p36, p37, p38, p39, p40, p41, p42, p43, p44, p45, p46, p47, p48, p49, p50, p51, p52, p53, p54, p55, p56, p57, p58, p59, p60, p61, p62], 'top');

        
        for tt = 1:length(g)
            % Reset previous stat sig labels
            g(tt).String = ' ';
            alphaT = 0.05/size(data,1);
            if str2double(chans_krus_wall_vals{tt} ) <= alphaT
                % Label images that are stat sig
                 g(tt) = text(hori(tt), vert(tt), 'X', 'HorizontalAlignment', 'Center', 'Color','w', 'FontWeight', 'bold', 'FontSize', 12.5, 'Parent', handles.krusWal);
            else
                continue;
            end
            
        end

        ind = tril(zeros(62,62));
        for i = 1:length(g)
            if g(i).String == 'X'
                ind(i) = 1;
            end
        end
        
        handles.g = g;

        [rowX,colX] = find(tril(ind,-1) == 1);
        SS_chan_inter = [rowX colX];
        
        % Assign Statistically Significant channel interactions and take
        % out the duplicates
        assignin('base', 'SS_chan_inter', SS_chan_inter)
        assignin('base', 'krusP', krusP)
        assignin('base', 'chans_krus_wall_vals', chans_krus_wall_vals)
        
            
'Done'
