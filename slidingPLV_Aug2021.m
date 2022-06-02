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
global data megM3d megR3d

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
    % filtdat - 62 Channels x 16000 Data Points X 30 Trials (First 8000 points are
    % Rest, 2nd 8000 pnts are Move
    filtdat = filterFGx(data, srate, handles.slider_freq.Value, handles.slider_width.Value);

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
global pSx pSy g hori vert SS_chan_inter

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
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLV CALCULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    synchmat1 = zeros( (size(filtdat,3)), (size(filtdat,1)), (size(filtdat,1)) );
    synchmat2 = zeros( (size(filtdat,3)), (size(filtdat,1)), (size(filtdat,1)) );
     
     % Upgraded PLV test-code
     [nchans, ns_rest, nt] = size(filtdat(:, (tidx1(1):tidx1(2)), :) );
     [nchans, ns_move, nt] = size(filtdat(:, (tidx2(1):tidx2(2)), :) );
     
     ndat_rest = filtdat(:, (tidx1(1):tidx1(2)), :) ./ abs(filtdat(:, (tidx1(1):tidx1(2)), :) );
     ndat_move = filtdat(:, (tidx2(1):tidx2(2)), :) ./ abs(filtdat(:, (tidx2(1):tidx2(2)), :) );

     % PLV - channel x channel x time window
%      plv_fast_rest = zeros(nchans, nchans, nt);
%      plv_fast_move = zeros(nchans, nchans, nt);
     
    for t = 1: nt
 %         plv_fast_rest(:,:, t) = abs(ndat_rest(:, :, t) * ndat_rest(:, :, t)') / ns_rest;
%          plv_fast_move(:,:, t) = abs(ndat_move(:, :, t) * ndat_move(:, :, t)') / ns_move;
        synchmat1(t,:, :) = abs(ndat_rest(:, :, t) * ndat_rest(:, :, t)') / ns_rest;
        synchmat2(t,:, :) = abs(ndat_move(:, :, t) * ndat_move(:, :, t)') / ns_move;
    end
    
%     plv_fast_rest = permute(plv_fast_rest, [3,1,2]);
%     synchmat1 = plv_fast_rest;
%     
%     plv_fast_move = permute(plv_fast_move, [3,1,2]);
%     synchmat2 = plv_fast_move;
%     
%      
%     'test'
    
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

global answers RestOnset  nchans krusP chans_krus_wall_vals alphaT RestOffset
global  pSx pSy g hori vert
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

     [ns_rest] = size(filtdat(:, (tidx1(1):tidx1(2)), :), 2 );
     [ nchans, ns_move, nt] = size(filtdat(:, (tidx2(1):tidx2(2)), :) );
     
     ndat_rest = filtdat(:, (tidx1(1):tidx1(2)), :) ./ abs(filtdat(:, (tidx1(1):tidx1(2)), :) );
     ndat_move = filtdat(:, (tidx2(1):tidx2(2)), :) ./ abs(filtdat(:, (tidx2(1):tidx2(2)), :) );

     % PLV - channel x channel x time window
%      plv_fast_rest = zeros(nchans, nchans, nt);
%      plv_fast_move = zeros(nchans, nchans, nt);
     
    for t = 1: nt
%         plv_fast_rest(:,:, t) = abs(ndat_rest(:, :, t) * ndat_rest(:, :, t)') / ns_rest;
%          plv_fast_move(:,:, t) = abs(ndat_move(:, :, t) * ndat_move(:, :, t)') / ns_move;
        synchmat1(t,:, :) = abs(ndat_rest(:, :, t) * ndat_rest(:, :, t)') / ns_rest;
        synchmat2(t,:, :) = abs(ndat_move(:, :, t) * ndat_move(:, :, t)') / ns_move;
    end
    
    %%%%
    
        % Kruskal-Wallis test on Phase-Locking Value data
        krusP = zeros(nchans,nchans);
        for ii = 1:nchans
            for jj = 1:nchans
                %KRUSKALWALLIS Nonparametric one-way analysis of variance (ANOVA).
%                 [krusP(ii,jj), tbl{ii,jj}, stats{ii,jj}] = kruskalwallis([synchmat1(:,ii,jj) synchmat2(:,ii,jj)], [], 'off');
                    [krusP(ii,jj)] = kruskalwallis([synchmat1(:,ii,jj) synchmat2(:,ii,jj)], [], 'off');
            end
        end
        % Arrange kruskal-wallis values for labeling in image
        chans_krus_wall_vals = num2cell(krusP);
        chans_krus_wall_vals = cellfun(@num2str, chans_krus_wall_vals, 'UniformOutput', false); % convert to string  
    
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
                 % Label images that are statistiacally significiant with an X
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
