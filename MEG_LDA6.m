
% Load data and rearrange the values
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
            
            % Look for inputs just in the command line
questions = { '\fontsize{12} Select Task: [1] Motor Img Calib [2] Motor Img Control';...
                    '\fontsize{12} Which subject number would you like to analyze (1-9)?';...
                    '\fontsize{12} What frequency would you like to analyze first (Hz)?';...
                    '\fontsize{12} What frequency would like to analyze last (Hz)?';...
                    '\fontsize{12} How much should the frequency increment by?';...
                    '\fontsize{12} What should be the frequency width for the filter?';...
                    '\fontsize{12} What should be the time window size (sec)?';...
                    '\fontsize{12} What time do you want to start at (sec)?';...
                    '\fontsize{12} Do you want to arrange filtered data by time points [1] or trials [2]?';...
                    '\fontsize{12} What will be the size of the arrangement?'};

% Input Params -> Delete Unused Channels -> 
% Concatenate Both Blocks data -> 
% Check for Offset Data and make same number of indices as Onset Data ->
% Rearrange Channel Numbers to go from Top-Left to Bottom-Right on Scalp Map ->
% Cut and rearrange data for equal Rest and Move trial sizes
% Delete first 3000 MEG Data Points in both blocks



    % Input Parameters - use defaults for now - Motor Img Cal, Subj 9
    definput = cell(10,1);
    definput(:) = [{'1'};{'9'};{'23.5'};{'26'};{'0.5'};{'4'};{'1'};{'1'};{'2'};{'50'}];
    options.Interpreter = 'tex';        
    answers = inputdlg(questions, 'Parameters for PLV', [1 85], definput, options );
    srate = 1000;

     % stud{1} -> MIcal; stud{2} -> MIcont
    stud = {MIcal; MIcont};  
    stud = stud{str2double(answers{1})};
    subj = str2double(answers{2});
    blocks = [load(stud{subj}(1,:)); load(stud{subj}(2,:))]; 

    % Delete unused channels
    blocks(1).meg(:,[4 6 7 11 15 16 20 23]) = [];
    blocks(2).meg(:,[4 6 7 11 15 16 20 23]) = [];

    % stack meg data together, delete 1st 3000 points in each block
    meg = [blocks(1).meg(3001:end, :); blocks(2).meg(3001:end, :)]';
    
    % if input task is not Motor Image Control, get labels 
    % from the stimulusCode (take out 1st 3000 points
     if all(all( convertCharsToStrings( stud{subj}(1,:)) ~= convertCharsToStrings(MIcont{subj}(1,:)), 2))
       labels = [blocks(2).stimulusCode(3001:end); blocks(2).stimulusCode(3001:end)];
     else
        % else task is motor img control
        % Take out 1st 3000 points. Max trial size is 12401 points
        % in between trials 999 0s
        % I want the signals that are 1 or 2 but I want to keep the index in
        % line with the meg data
        labels = [blocks(1).stateSignals.targetCode(3001:end); blocks(2).stateSignals.targetCode(3001:end)];
        % Initialize clbl - class label
        clbl = cell(length(labels),1);
        % create temporary indices for all labels that are 1 and 2
        temp_ind = find(labels>0 & labels<3);
        for jj = 1:length(temp_ind)
            if labels(temp_ind(jj)) == 1
                clbl{temp_ind(jj)} = 'Move';
            elseif labels(temp_ind(jj)) == 2
                clbl{temp_ind(jj)} = 'Rest';
            end
        end
        % All I care about is meg & label data during the task trials (1 & 2)
        labels = labels(temp_ind);
        meg = meg(temp_ind);
    end
    

    % Concatenate RestOnset and MoveOnset times for the 2 blocks
    RestOnset = [blocks(1).RestOnset(:) - 3000; blocks(2).RestOnset(:)+length(blocks(1).meg(3001:end, :)) - 3000];
    MoveOnset = [blocks(1).MoveOnset(:) - 3000; blocks(2).MoveOnset(:)+length(blocks(1).meg(3001:end, :)) - 3000];
    %CHECK FOR OFFSET DATA. 
    % IF THERE'S OFFSET DATA  MAKE THE SAME SIZE OF OFFSET TIMES AS ONSET TIMES
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
                RestOffset = [blocks(1).RestOffset(:) - 3000; blocks(2).RestOffset(:)+length(blocks(1).meg(3001:end, :)) - 3000];
                MoveOffset = [blocks(1).MoveOffset(:) - 3000; blocks(2).MoveOffset(:)+length(blocks(1).meg(3001:end, :)) - 3000];
                trialStart = [blocks(1).trialStart(:) - 3000; blocks(2).trialStart(:)+length(blocks(1).meg(3001:end, :)) - 3000];
                trialEnd = [blocks(1).trialEnd(:) - 3000; blocks(2).trialEnd(:)+length(blocks(1).meg(3001:end, :)) - 3000];       
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
           % If there is no offset data, then catch the error and display message
            catch
                disp('There is no MoveOffset or RestOffset in these trial blocks');
             end


    % REARRANGE THE CHANNEL NUMBERS TO GO FROM TOP-LEFT TO BOTTOM-RIGHT ON
    % SCALP MAP
    megarr = [12 7 6 17 14 8 11 16 5 10 4 15 9 3 1 13 2]';
    megmat = meg([megarr], :);
    meg = megmat;

    % DEPENDING ON WHETHER THE DATA IS MOTOR CONTROL OR MOTOR IMAGERY,
    % CUT & REARRANGE DATA TO HAVE EQUAL REST AND MOVE EPOCH SIZES
        if all(all( convertCharsToStrings(stud{subj}) == convertCharsToStrings(MIcont{subj}), 2))
            % When the data is motor imagery - control
            megR3d = cell(size(MoveOnset, 1), 1);
            megM3d = cell(size(MoveOnset, 1), 1);

            filtdatM = [];
            filtdatR = [];
            
            % Cut out unrelated data in trials - 1 second before, 2 seconds after
            for i = 1:size(MoveOnset, 1)
                % Rest
                megR3d{i} = meg( :, (RestOnset(i) + 1000) : RestOffset(i) - 2000 );
                % Move
                megM3d{i} = meg( :, (MoveOnset(i) + 1000) : MoveOffset(i) - 2000 );
            end
            
            % Concatenate data from all trials based on Move or Rest Task
            megM = [];
            megM = [megM megM3d{:}];
            megR = [];
            megR = [megR megR3d{:}];
            % Filter data -> filtdat = filterFGx(data, srate, handles.slider_freq.Value, handles.slider_width.Value);

            % Filter the data - Morlet Wavelet
            filtdatM = filterFGx(megM, srate, str2num(answers{3}), str2num(answers{4}) );
            filtdatR = filterFGx(megR, srate, str2num(answers{3}), str2num(answers{4}) );
          
            % Move and Rest data become same size, Round down
            % on 2nd dim (Chans x spectData) from whichever matrix (Rest or
            % Move) is smaller
                    roundnum = str2double(answers{10});
                    if size(filtdatR,2) <= size(filtdatM,2)
                            filtdatR = filtdatR(: , 1:roundnum*(floor(length(filtdatR)/roundnum)));
                            filtdatM = filtdatM(: , 1:roundnum*(floor(length(filtdatR)/roundnum)));
                    else
                            filtdatM = filtdatM(: , 1:roundnum*(floor(length(filtdatM)/roundnum)));
                            filtdatR = filtdatR(: , 1:roundnum*(floor(length(filtdatM)/roundnum)));
                    end

            % Break up filtered Rest & Move data into (Chans x data points x Trials)
                if str2double(answers{9}) == 1
                    % Round by time points
                       filtdat = [reshape(filtdatR, size(filtdatR,1), roundnum, [])       reshape(filtdatM, size(filtdatM,1), roundnum, [])];
                elseif str2double(answers{9}) == 2
                    % Round by trials
                      filtdat = [reshape(filtdatR, size(filtdatR,1), [] , roundnum)      reshape(filtdatM, size(filtdatR,1), [], roundnum)];      
                end

    % When the data is motor imagery - calibration, create 3-d variables
    % channels x time_points x trials
        elseif all(all( convertCharsToStrings(stud{subj}(1,:) ) ~= convertCharsToStrings(MIcont{subj}(1,:) ), 2))
            megR3d = zeros(17,5000,size(RestOnset,1));
            megM3d = zeros(17,5000,size(MoveOnset,1));
            data = zeros(17,10000,size(RestOnset,1));

            % Extract Rest & Move data and filter
            % Morlet Wavelet, 
            for i = 1:size(RestOnset,1)
                megR3d(:,:,i) = meg(:,(RestOnset(i):MoveOnset(i)-1));
                megM3d(:,:,i) = meg(:,(MoveOnset(i):MoveOnset(i)+4999));
                data(:,:,i) = [megR3d(:,:,i) megM3d(:,:,i)];
            end
            
            filtdat = filterFGx(data, srate, str2num(answers{3}), str2num(answers{5}) );
        end

    % Sampling Frequency and filter coefficients
        fs = 1000/100;
    %     [hp2 hp1] = butter( 4, 8/(fs/2) , 'high');
        [lp2, lp1] = butter( 4, 2/(fs/2)  , 'low'); 
    %     
    %     [hp4 hp3] = butter( 4, 12.5/(fs/2) , 'high');
    %     [lp4 lp3] = butter( 4, 30/(fs/2)  , 'low'); 

    % Time window, shift, # of pnts
       tw = 500;
       shift = 100;
    % number of pnts is length of meg
       npnts = length(meg);
    % meg time window size
       meg_tw = cell((size(meg, 2) - tw)/shift + 1, 2 );
    % windowed meg frequency (Short-time fft)
       meg_f = cell((size(meg, 2) - tw)/shift + 1, 1 );
    % meg low pass
       x_lp = cell((size(meg, 2) - tw)/shift + 1, 1 );
       a = 'alpha';
       a_val = cell((size(meg, 2) - tw)/shift + 1, 2 );
       b = 'beta';
       b_val = cell((size(meg, 2) - tw)/shift + 1, 2 );
       d = 'data';
       d_val = cell((size(meg, 2) - tw)/shift + 1, 2 );
       
       % create a structure of meg-> time smoothed 0.5 sec -> normalized
       % all data; alpha band power features; beta band power features
       x_lp_norm = struct('alpha', {a_val}, 'beta', {b_val}, 'data', {d_val});

    %     x_lp_norm_medavg = cell((size(meg, 2) - tw)/shift + 1, 2 );
       feat = zeros((size(meg, 2) - tw)/shift + 1, 2);
       alpha = 8:13;
       beta = 25:30;
       bandpass = [alpha beta]
%        bandpass = length(alpha) + length(beta);
       % initialize group label cells - default val - size(group) = 4996 x 1
       
       group = cell( length(bandpass)/2 * ((size(meg, 2) - tw)/shift + 1), 1 );
% meg_tw{i,1} --> channels x time_window x trials
       for i = 1:(size(meg, 2) - tw)/shift + 1
           % If i == 1, save meg_tw{i,1} as the 1st 500 time points
           % else shift the 500 point time window every 100 points 
           % take fft of each time window
               if i == 1
                   % store meg data in meg_tw, needs to iterate with i
                    meg_tw{i,1} = meg(:, (1:tw + (i-1)*500), :);
                    meg_tw{i,2} = labels(tw);
               else
                    meg_tw{i,1} = meg( :, (1+(shift*(i-1) )) : (tw+(shift*(i-1) )), :);
                    meg_tw{i,2} = labels(tw+(shift*(i-1) ) );
               end

                if meg_tw{i,2} == 1
                     meg_tw{i,2} = 'Rest';
                else
                    meg_tw{i,2} = 'Move';
                end
               % take fft of the data
               meg_f{i, 1} = (fft(meg_tw{i})/npnts);
               % time smooth data at 2 hz, find fft power
               x_lp{i, 1} = 2*abs(filtfilt(lp2, lp1, meg_f{i, 1})).^2;
               % normalize data
               x_lp_norm.data{i, 1} = normalize(x_lp{i, 1}, 2, 'zscore');
               % separate alpha and beta data
               x_lp_norm.alpha{i, 1} = x_lp_norm.data{i, 1}(alpha);
               x_lp_norm.beta{i, 1} = x_lp_norm.data{i, 1}(beta);
               % add class labels to alpha and beta data
               x_lp_norm.alpha{i, 2} = repmat({meg_tw{i,2}}, length(alpha),1) ;  
               x_lp_norm.beta{i, 2} = repmat({meg_tw{i,2}}, length(beta),1) ;  
               x_lp_norm.data{i, 2} = repmat({meg_tw{i,2}}, length(bandpass),1) ;  


    %            x_lp_norm_medavg{i, 1} = median(mean(x_lp_norm{i}) );
    %            x_lp_norm_medavg{i, 2} = meg_tw{i,2};   
                for j = 1:length(bandpass)
                     group ([( (i-1) * length(bandpass))+(j)]) =  x_lp_norm.data{i, 2}(j);
                end
    %              group { [( (i-1) * length(bandpass))+(1:length(bandpass))] } =  x_lp_norm{i, 2};
       end

    dataA = []; dataB = [];
    dataA = [dataA x_lp_norm.alpha{:,1}];
    dataB = [dataB x_lp_norm.beta{:,1}];
% FIND CVPARTITION1 AND ADD TO PATH
    [ testIdx, trainIdx, nTest] = cvpartition1(size(x_lp_norm.data,1),5);
    g = kron(trainIdx,ones(1,12));
    gs = kron(trainIdx,ones(1,6));
    % Group must equal number of rows as training - Default vals -
    % size( find( g(i,:) ) ) = 1 x 47952
    % Initiate classify params
    C = cell(5,1); err = zeros(5,1); P = cell(5,1); logp = cell(5,1); coeff = cell(5,1);
    for i = 1:5 
        % Classify the data
        % classify( test data, training data, group labels );
        % size(group) -> 59952 x 1
        % size(g) -> 5 x 59952 (4996*12=59952);    size(find(g(1,:))) -> 1 x 47952
        
        [C{i, 1},err(i, 1),P{i, 1},logp{i, 1},coeff{i, 1}] = classify( [ cell2mat(x_lp_norm.alpha(testIdx(i,:)))' ;cell2mat(x_lp_norm.beta(testIdx(i,:)))' ], ...
                                                                    [cell2mat(x_lp_norm.alpha(trainIdx(i,:)))'; cell2mat(x_lp_norm.beta(trainIdx(i,:)))' ] , ...
                                                                    group(find(g(i,:))), 'linear' );

%         % Find min and max spectral power in Alpha and Beta freq bands
%         Extrema = [min( cell2mat(x_lp_norm.alpha(testIdx(i,:))) ); max(cell2mat( x_lp_norm.alpha(testIdx(i,:))) ); min( cell2mat(x_lp_norm.beta(testIdx(i,:))) ); max( cell2mat(x_lp_norm.beta(testIdx(i,:))) ) ];
%         % Create mesh grid for plots based on the Extrema of the alpha & beta bands
%         [X,Y] = meshgrid(linspace( floor(Extrema(1)), ceil(Extrema(2)), length( cell2mat(x_lp_norm.alpha(testIdx(i,:))) ) ), linspace( floor(Extrema(3)), ceil(Extrema(4)) ) );
%         X = X(:); Y = Y(:);
%         figure;
%         % use trainIdx data for gscatter
%         % error using gscatter -> X & Y must have the same length. 
%         % group must have the same # of values as rows of X
%         % Alpha is X and Beta is Y to distinguish the 2 features,
%         
%         % Need to modify groups. X is Alpha and Y is beta, their size is 23976 (4996 x 6)
% %         h1 = gscatter( [cell2mat(x_lp_norm.alpha(trainIdx(i,:)))'; cell2mat(x_lp_norm.beta(trainIdx(i,:)))' ] , [ cell2mat(x_lp_norm.alpha(testIdx(i,:)))' ;cell2mat(x_lp_norm.beta(testIdx(i,:)))' ], group(find(g(i,:))),'rb','v^',[],'off');
%         h1 = gscatter(  [cell2mat(x_lp_norm.alpha(trainIdx(i,:)))'], [cell2mat(x_lp_norm.beta(trainIdx(i,:)))' ] ,  group(find(gs(i,:))),'rb','v^',[],'off');
% 
%         set(h1,'LineWidth',2)
%         hold on;
%         gscatter(X,Y,C{i},'rb','.',1,'off');
%         K = coeff{i}(1,2).const;
%         L = coeff{i}(1,2).linear;
% 
%         f = @(x,y) K + [x y]*L;
% 
%         h2 = fimplicit(f, [ floor(Extrema(1)), ceil(Extrema(2)), floor(Extrema(3)), ceil(Extrema(4)) ] );
%         set(h2,'Color','m','LineWidth',2)
%         axis( [ floor(Extrema(1)) ceil(Extrema(2)) floor(Extrema(3)) ceil(Extrema(4))] )
%         xlabel('Normalized Alpha Power')
%         ylabel('Normalized Beta Power')
%         title('{\bf LDA of Rest vs Move with Alpha & Beta Band Power}')
%         legend('Rest Data','Move Data', 'Class Rest', 'Class Move', 'Line of Separation', 'Location','NW');
    end

    % Normalize the results


    % Run 5-fold cross-validation with 5 groups made up of 500ms time windows

    % Generate LDA Model and test it






