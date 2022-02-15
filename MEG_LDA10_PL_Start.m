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

% Input Parameters - use defaults for now - Motor Img Cal, Subj 9
definput = cell(10,1);
definput(:) = [{'1'};{'9'};{'8.5'};{'12'};{'0.5'};{'3'};{'0.5'};{'0'};{'2'};{'50'}];
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
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Sampling Frequency and filter coefficients
        fs = 1000;
        nfft = 1028;
        freq = 0:fs/nfft:fs/2;
        chans = size(meg,1);

        
    %     [hp2 hp1] = butter( 4, 8/(fs/2) , 'high');
        [lp2, lp1] = butter( 4, 2/(fs/2)  , 'low'); 
    %     
    %     [hp4 hp3] = butter( 4, 12.5/(fs/2) , 'high');
    %     [lp4 lp3] = butter( 4, 30/(fs/2)  , 'low'); 

    % Time window, shift, # of pnts
       tw = 500;
       shift = 100;fs
    % number of pnts is length of meg
       npnts = length(meg);
    % meg time window size
       meg_tw = cell((size(meg, 2) - tw)/shift + 1, 2 );
       
    % Filtered meg time window data
       Fmeg_tw = cell((size(meg, 2) - tw)/shift + 1, 2 );
    % windowed meg frequency (Short-time fft)
       meg_f = cell((size(meg, 2) - tw)/shift + 1, 1 );
       
       meg_fft_singleTW = zeros(chans, nfft);
       
       x_abs = zeros(chans, nfft);
       x_pow = zeros(chans, (nfft/2 + 1) );
       
       
    % meg low pass
%        x_abs = cell((size(meg, 2) - tw)/shift + 1, 1 );
%        x_pow = cell((size(meg, 2) - tw)/shift + 1, 1 );
       a = 'alpha';
       a_val = cell((size(meg, 2) - tw)/shift + 1, 1 );
       b = 'beta';
       b_val = cell((size(meg, 2) - tw)/shift + 1, 1 );
       d = 'data';
       d_val = cell( ((size(meg, 2) - tw)/shift + 1), 1 );
       
       % create a structure of meg-> time smoothed 0.5 sec -> normalized
       % all data; alpha band power features; beta band power features
       x_pow_str = struct('alpha', {a_val}, 'beta', {b_val}, 'data', {d_val});

    %     x_lp_str_medavg = cell((size(meg, 2) - tw)/shift + 1, 2 );
       feat = zeros((size(meg, 2) - tw)/shift + 1, 2);
       alpha = 8:13;
       beta = 25:30;
       % bandpass = length(alpha) + length(beta);
       % initialize group label cells - default val - size(group) = 4996 x 1
       
       w = hamming(tw);
       srate = 1000;
       freq_val = str2num(answers{3});
       freq_width = str2num(answers{6})
       
       PLV = cell( (size(meg, 2) - tw)/shift + 1, 1);
       d = cell( (size(meg, 2) - tw)/shift + 1, 1);
       B = cell( (size(meg, 2) - tw)/shift + 1, 1);
       I = cell( (size(meg, 2) - tw)/shift + 1, 1);
       d_sort = cell( (size(meg, 2) - tw)/shift + 1, 1);
       
       for i = 1:( (size(meg, 2) - tw)/shift + 1 )
           % If i == 1, save meg_tw{i,1} as the 1st 500 time points
           % else shift the 500 point time window every 100 points 
           % meg_tw{i,1} --> channels x time_window
               if i == 1
                   % store meg data in meg_tw, needs to iterate with i
                   
                    meg_tw{i,1} = meg(:, (1:tw) )' .* w;
                    meg_tw{i,2} = labels(tw);
               else
                    meg_tw{i,1} = meg( :, (1+(shift*(i-1) )) : (tw+(shift*(i-1) )), :)' .* w;
                    meg_tw{i,2} = labels(tw+(shift*(i-1) ) );
               end
               
               meg_tw{i,1} = meg_tw{i,1}';

               % Check labels for Rest or Move
                if meg_tw{i,2} == 1
                     meg_tw{i,2} = 'Rest';
                elseif meg_tw{i,2} == 2
                    meg_tw{i,2} = 'Move';
                end 
                
                % Angles from Hilbert Transform to measure phase difference
                Fmeg{i,1} = filterFGx(meg_tw{i,1}, srate, freq_val, freq_width);
                % filtdata - chans x data
                for j = 1: size(meg_tw{i,1},1)
                    for k = 1: size(meg_tw{i,1},2)
                        % angle is in radians - chans x data
                        % Check Hilbert Angle
                        angts(j,k) = angle(hilbert(Fmeg{i,1}(j,k)').');
                    end
                end
                
                % Phase Locking Value
                %Create variables
                PLV{i,1} = zeros(size(meg_tw{i,1},1), size(meg_tw{i,1},1));
%                 d{i,1} = zeros(size(meg_tw{i,1},1), size(meg_tw{i,1},1));
                d{i,1} = zeros(size(meg_tw{i,1},1), 1);
                B{i,1} = zeros(size(meg_tw{i,1},1), 1);
                I{i,1} = zeros(size(meg_tw{i,1},1), 1);
                d_sort{i,1} = zeros(size(10, 1));


                % chans
                for chani = 1: size(meg_tw{i,1},1)
                    % chans
                    for chanj = 1: size(meg_tw{i,1},1)
                        tmpAi = angts(chani, :);
                        tmpAj = angts(chanj, :);
                        % PLV - chans x data
                        PLV{i,1}(chani, chanj) = mean(abs(mean(exp(1i*(tmpAi-tmpAj )), 2)),1);
                    end
                             
                    % Mahal distance formula
                    % ?( PLV value - mean(PLV of whole channel) ) ---->
                    % 4996 time wins, 17x17 values
                    for chank = 1: size(meg_tw{i,1},1)
                        d{i,1}(chani, chank) = sqrt( (PLV{i,1}(chani, chank) - mean( PLV{i,1}(chani, :) ))' * (cov(PLV{i,1}(chani,:)))^-1 * ( PLV{i,1}(chani, chank) - mean( PLV{i,1}(chani, :) )) );
                    end                
%                     d{1,1} < mean(mean(d{1,1}))              
                end
                
%                 d{i,1} = abs(sqrt( (PLV{i,1}(chani, :) - mean( PLV{i,1},2 ))' * (cov(PLV{i,1}(chani,:)))^-1 * ( PLV{i,1}(chani) - mean( PLV{i,1},2 )) ));
                % Put the distances in order
                [ B{i,1},I{i,1} ] = sort(d{1,1}, 'ascend');
                
                % Keep only the top 10 values of 17
                d_sort{i,1} = d{i,1}(I{i,1}(1:10));
                
%                 [B,I] = sort(d{1,1}(1,:), 'ascend') 
%                 d_sort = d{1,1}(:,I(1:10))
                
%                 
%                % take fft of the data / 500 point time window
% %                meg_f{i, 1} = (fft(meg_tw{i}, nfft))/tw;
%                 for gg = 1 :  size(meg_tw{1},1)
%                     %Add in Hamming window w if used.
%                     meg_fft_singleTW(gg,:) = (fft(meg_tw{i}(gg,:), nfft)/(tw/sum(w))); 
%                    % take absolute value, find amplitude, Scale by 1 million or
%                    % the power will be too small for matlab (10^-18)
%                    x_abs(gg, :) = abs(meg_fft_singleTW(gg,:));
%                    x_pow(gg, :) = x_abs(gg, 1:nfft/2 + 1);
%                    x_pow(gg, 2:(end - 1)) = (2*1e09*x_pow(gg, 2:end-1)).^2;
%                    
% %                   % normalize data
% %                   x_pow_norm(gg,:) = normalize(x_pow{i, 1}(gg,:), 2, 'zscore');
%                    
% %                figure; plot(freq, x_pow{i,1});
% %                xlabel('Frequency (Hz)');
% %                ylabel('Power (something squared / Hz)');
%                
%                % separate alpha and beta data
%                    x_pow_str.alpha{i, 1} = [x_pow_str.alpha{i, 1}; sum(x_pow(gg, alpha))];
%                    x_pow_str.beta{i, 1} = [x_pow_str.beta{i, 1}; sum(x_pow(gg, beta))];
%                % 4996 cells, 2 freq bands x 17 channels
%                % I want 34 x 1 with [chan1alph; chan1beta; chan2alph;....chan17beta]
%                    x_pow_str.data{i, 1} = [x_pow_str.alpha{i, 1}; x_pow_str.beta{i,1}];
%                 end

       end
       
       
% Mahalanobis Distance


       
% fs_d = fs/shift;
%  
% maxdata = [];
% maxdata = [ maxdata    x_pow_str.data{:} ];
% maxdata = maxdata';
% 
% [lp2, lp1] = butter( 4, 2/(fs_d/2)  , 'low'); 
% maxdata_lp = filtfilt(lp2,lp1, maxdata);
% 
% maxdata_lp_z = zscore(maxdata_lp,2);
chek = zeros(17,1);

for z = 1:17
    top_d = z;
    for win = 1:( (size(meg, 2) - tw)/shift + 1 )

         d_sort{win,1} = zeros(size(top_d, 1));
         d_sort{win,1} = d{win,1}(I{win,1}(1:top_d));
    end

    feat = [];
    % feat = [feat   d_sort{:}];
    feat = [feat   d_sort{:}];
    feat = feat';

    CL = cell(size(meg_tw,1), 1);
    for zz = 1:size(meg_tw,1)
        CL{zz} = meg_tw{zz,2};
    end
    % size(CL) needs to be 4996 x 1
    % size(x_pow_str.data) = 4996 x 34
    % so we have 34 dim data

    % FIND CVPARTITION1 AND ADD TO PATH
    % cvpartition1 -> boolean value seperates train & test data 80/20 and
    % creates 5-fold data to test on
        [ testIdx, trainIdx, nTest] = cvpartition1(size(d_sort,1),5);

            C = cell(5,1); err = zeros(5,1); P = cell(5,1); logp = cell(5,1); coeff = cell(5,1);
        for i = 1:5 
            % Classify the data
            % classify( test data, training data, group labels );
            % size(group) -> 59952 x 1
            % size(g) -> 5 x 59952 (4996*12=59952);    size(find(g(1,:))) -> 1 x 47952

              [C{i, 1},err(i, 1),P{i, 1},logp{i, 1},coeff{i, 1}] = ...
                                                                        classify( feat(testIdx(i,:),:), ...
                                                                        feat(trainIdx(i,:),:) , ...
                                                                        CL(find(trainIdx(i,:))), 'linear' );


    %           [C{i, 1},err(i, 1),P{i, 1},logp{i, 1},coeff{i, 1}] = ...
    %                                                                             classify( maxdata_lp(testIdx(i,:),:), ...
    %                                                                             maxdata_lp(trainIdx(i,:),:) , ...
    %                                                                             CL(find(trainIdx(i,:))), 'linear' );

    %           [C{i, 1},err(i, 1),P{i, 1},logp{i, 1},coeff{i, 1}] = ...
    %                                                                             classify( maxdata(testIdx(i,:),:), ...
    %                                                                             maxdata(trainIdx(i,:),:) , ...
    %                                                                             CL(find(trainIdx(i,:))), 'linear' );



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

        chek(z) = mean(err);
end
chek