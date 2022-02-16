
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


definput = cell(10,1);
definput(:) = [{'1'};{'9'};{'23.5'};{'26'};{'0.5'};{'4'};{'1'};{'1'};{'2'};{'50'}];
options.Interpreter = 'tex';        
answers = inputdlg(questions, 'Parameters for PLV', [1 85], definput, options );

 % stud{1} -> MIcal; stud{2} -> MIcont
stud = {MIcal; MIcont};  
stud = stud{str2double(answers{1})};
subj = str2double(answers{2});
blocks = [load(stud{subj}(1,:)); load(stud{subj}(2,:))]; 

RestOnset = [blocks(1).RestOnset(:); blocks(2).RestOnset(:)+length(blocks(1).meg)];
MoveOnset = [blocks(1).MoveOnset(:); blocks(2).MoveOnset(:)+length(blocks(1).meg)];
% If data has Offset data, make the same number of trials as Onset data
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
         end

% Delete unused channels
blocks(1).meg(:,[4 6 7 11 15 16 20 23]) = [];
blocks(2).meg(:,[4 6 7 11 15 16 20 23]) = [];

meg = [blocks(1).meg; blocks(2).meg];
meg = meg';

% Rearrange the channel numbers to go sequentially from top-left to bottom-right
megarr = [12 7 6 17 14 8 11 16 5 10 4 15 9 3 1 13 2]';
megmat = meg([megarr], :);
meg = megmat;

% Depending on whether the data is Motor Control or Motor imagery,
% cut and rearrange data to have equal Move and Rest Trials
    if all(all( convertCharsToStrings(stud{subj}) == convertCharsToStrings(MIcont{subj}), 2))
        megR3d = cell(size(MoveOnset, 1), 1);
        megM3d = cell(size(MoveOnset, 1), 1);

        filtdatM = [];
        filtdatR = [];
        % Cut out unrelated data in trials
        % 1 second before, 2 seconds after
        for i = 1:size(RestOnset, 1)
            % Rest
            megR3d{i} = meg( :, (RestOnset(i) + 1000) : RestOffset(i)-2000 );
            % Move
            megM3d{i} = meg( :, (MoveOnset(i) + 1000) : MoveOffset(i)-2000 );
    %         
    %         filtdatR =  [filtdatR filterFGx(megR3d{i}, srate, handles.slider_freq.Value,handles.slider_width.Value)];
    %         filtdatM =  [filtdatM filterFGx(megM3d{i}, srate, handles.slider_freq.Value, handles.slider_width.Value)];
        end

        % Make filtered move and rest data the same size, Round down from roundnum on
        % the smallest matrix
                roundnum = str2double(answers{10});
                if size(filtdatR,2) <= size(filtdatM,2)
                        filtdatR = filtdatR(: , 1:roundnum*(floor(length(filtdatR)/roundnum)));
                        filtdatM = filtdatM(: , 1:roundnum*(floor(length(filtdatR)/roundnum)));
                else
                        filtdatM = filtdatM(: , 1:roundnum*(floor(length(filtdatM)/roundnum)));
                        filtdatR = filtdatR(: , 1:roundnum*(floor(length(filtdatM)/roundnum)));
                end

            if str2double(answers{9}) == 1
                 % Round by time points
                   filtdat = [reshape(filtdatR, size(filtdatR,1), roundnum, [])       reshape(filtdatM, size(filtdatM,1), roundnum, [])];
            elseif str2double(answers{9}) == 2
                % Round by trials
                  filtdat = [reshape(filtdatR, size(filtdatR,1), [] , roundnum)      reshape(filtdatM, size(filtdatR,1), [], roundnum)];      
            end

    elseif all(all( convertCharsToStrings(stud{subj}) ~= convertCharsToStrings(MIcont{subj}), 2))
        megR3d = zeros(17,5000,size(RestOnset,1));
        megM3d = zeros(17,5000,size(MoveOnset,1));
        data = zeros(17,10000,size(RestOnset,1));

        % motor img - calib
        for i = 1:size(RestOnset,1)
            megR3d(:,:,i) = meg(:,(RestOnset(i):MoveOnset(i)-1));
            megM3d(:,:,i) = meg(:,(MoveOnset(i):MoveOnset(i)+4999));
            data(:,:,i) = [megR3d(:,:,i) megM3d(:,:,i)];
        end
    %     filtdat = filterFGx(data, srate, handles.slider_freq.Value, handles.slider_width.Value);
    end
 
% BREAK UP DATA INTO 500ms TIME WINDOWS,
% Recheck Mike X Cohen vids.
% (DONT BREAK UP MOVE AND REST THEN SEPARATELY) DO abs(FFT, SQ, BP (feature), NORMALIZE (feature scaling), 
%   RANDPERM, PASS CELLS (Both move and
%  Rest) INTO 5 GROUPS for 5-FOLD CV. DO LDA

% Sampling Frequency and filter coefficients
    fs = 1000;
    [hp2 hp1] = butter( 4, 8/(fs/2) , 'high');
    [lp2 lp1] = butter( 4, 12.5/(fs/2)  , 'low'); 
    
    [hp4 hp3] = butter( 4, 12.5/(fs/2) , 'high');
    [lp4 lp3] = butter( 4, 30/(fs/2)  , 'low'); 
% Time window, shift, # of pnts
   tw = 500;
   shift = 100;
   npnts = length(meg);
   meg_tw = cell((size(meg, 2) - tw)/shift + 1, 1 );
   meg_f = cell((size(meg, 2) - tw)/shift + 1, 1 );
   x_lp{i} = cell((size(meg, 2) - tw)/shift + 1, 1 );
   x_lp_norm{i} = cell((size(meg, 2) - tw)/shift + 1, 1 );
   feat = zeros((size(meg, 2) - tw)/shift + 1, 2);
   
   for i = 1:(size(meg, 2) - tw)/overlap + 1
           if i == 1
                meg_tw{i} = meg(:, 1:tw, :);
           else
                meg_tw{i} = meg( :, (1+(overlap*(i-1) )) : (tw+(overlap*(i-1) )), :);
           end
           meg_f{i} = (2*abs(fft(meg_tw{i})/npnts)).^2;
              
%        x_hpA{i} = filtfilt(hp2, hp1, x);
           x_lp{i} = filtfilt(lp2, lp1, meg_f{i});
           x_lp_norm{i} = normalize(x_lp{i}, 2, 'zscore');
%        x_hp_lp_norm{i} = reshape(x_hp_lp_normA, size(x_hp_lpA) );

%        x_hpB = filtfilt(hp4, hp3, x);
%        x_hp_lpB = filtfilt(lp4, lp3, x_hpB);
%        x_hp_lp_normB = normalize(x_hp_lpB, 2, 'zscore');   
%        x_hp_lp_normB = reshape(x_hp_lp_normB, size(x_hp_lpB) );

       feat(i, 1) = median(mean(x_lp_norm{i}));
%        feat(i, 2) = median(mean(epochsB{i}));
        
   end
   feat(1:26, :) = [];

% Frequency X-Axis
%    hz = linspace(0, srate/2, floor(npnts/2)+1);
% MEG spectral power
%    x = (2*abs(fft(meg)/npnts)).^2;


% Extract Features such as alpha and beta band power
    %BP Filter alpha and beta band - butterworth, filtfilt
    % abs(x)^2
    % Average across trials mean(:,:,#)
    % Create time windows of 500ms
    % Average over time windows

   
  
%    % Convert normalized data back to epochs
%    % repmat 500, 5056, 1
%    epochsA = cell( (size(meg, 2) - tw)/shift + 1, 1);
%    epochsB = cell( (size(meg, 2) - tw)/shift + 1, 1);
%    feat = zeros((size(meg, 2) - tw)/shift + 1, 2);
%    for i = 1:(size(meg, 2) - tw)/shift + 1
%        if i == 1
%             epochsA{i} = x_hp_lp_normA(:, 1:tw);
%             epochsB{i} = x_hp_lp_normB(:, 1:tw);
% %             epochsA{i} = x_hp_lpA(:, 1:tw);
% %             epochsB{i} = x_hp_lpB(:, 1:tw);
% 
%        else
%             epochsA{i} =  x_hp_lp_normA( :, (1+(shift*(i-1) )) : (tw+(shift*(i-1) )));
%             epochsB{i} =  x_hp_lp_normB( :, (1+(shift*(i-1) )) : (tw+(shift*(i-1) )));
% %             epochsA{i} =  x_hp_lpA( :, (1+(shift*(i-1) )) : (tw+(shift*(i-1) )));
% %             epochsB{i} =  x_hp_lpB( :, (1+(shift*(i-1) )) : (tw+(shift*(i-1) )));
%        end
%        feat(i, 1) = median(mean(epochsA{i}));
%        feat(i, 2) = median(mean(epochsB{i}));
%    end
%    feat(1:26, :) = [];
   
   
   
% Randomly Rearrange epoch cells
%    epochcv = epochs(randperm(numel(epochs)));
 

% I dont want to add part of a window to any of these cells, it won't represent movement even if it has movement       
cv = cell(5,2);
for i = 1:5
cv{i,1} = [epochsA{((i-1)*floor(numel(epochsA)/5))+(1:floor(numel(epochsA)/5))}];
cv{i,2} = [epochsB{((i-1)*floor(numel(epochsB)/5))+(1:floor(numel(epochsB)/5))}];
end
   
   % Reorganize data in cells to 5 groups - 4 train & 1 test
   %classify(graph area, data, group labels, 'linear') 
   
labels = [blocks(1).stimulusCode; blocks(2).stimulusCode];
% delete zero labels
labels(1:3000) = [];
% graph area is min and max of values and indices
Extrema = [min(feat(:,1)); max(feat(:,1)); min(feat(:,2)); max(feat(:,2))];
[X,Y] = meshgrid(linspace( floor(Extrema(1)), ceil(Extrema(2)) ), linspace( floor(Extrema(3)), ceil(Extrema(4)) ) );
X = X(:); Y = Y(:);
labelsD = downsample(labels,100);
labelsN = cell( size(labelsD,1), 1);
for i = 1:size(labelsD,1)
    if labelsD(i) == 1
        labelsN{i} = 'Rest';
    else
        labelsN{i} = 'Move';
    end
end

figure;
h1 = gscatter(feat(:,1),feat(:,2),labelsN,'rb','v^',[],'off');
set(h1,'LineWidth',2)
   
[C,err,P,logp,coeff] = classify([X Y], [feat(:,1) feat(:,2)], labelsN, 'linear');

hold on;
gscatter(X,Y,C,'rb','.',1,'off');
K = coeff(1,2).const;
L = coeff(1,2).linear;
% Q = coeff(1,2).quadratic;
% Function to compute K + L*v + v'*Q*v for multiple vectors
% v=[x;y]. Accepts x and y as scalars or column vectors.
% f = @(x,y) K + [x y]*L + sum(([x y]*Q) .* [x y], 2);
f = @(x,y) K + [x y]*L, 2;

h2 = ezplot(f,[ floor(Extrema(1)) ceil(Extrema(2)) floor(Extrema(3)) ceil(Extrema(4))] );
set(h2,'Color','m','LineWidth',2)
axis( [ floor(Extrema(1)) ceil(Extrema(2)) floor(Extrema(3)) ceil(Extrema(4))] )
xlabel('Normalized Alpha Power')
ylabel('Normalized Beta Power')
title('{\bf LDA of Rest vs Move with Alpha & Beta Band Power}')
legend('Rest Data','Move Data', 'Class Rest', 'Class Move', 'Line of Separation', 'Location','NW');



% x as alpha pow and y as beta pow, groups as class labels
% Average over channels or grab either minimum or maximum of points from

% channel
% Average time points
% 2 Features

   
   % x_hp_lp_norm1 = (x_hp_lp(:) - min(x_hp_lp(:))) / (max(x_hp_lp(:)) - min(x_hp_lp(:) ));
   % x_hp_lp_norm1 = reshape(x_hp_lp_norm1, size(meg));
   % data (506000*5) - labels (503000*5) = 13000/500 = 26 time windows to
   % remove
   
   
%%

% Normalize the results


% Run 5-fold cross-validation with 5 groups made up of 500ms time windows

% Generate LDA Model and test it





