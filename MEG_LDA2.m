
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

% Add as input for later
 % stud{1} -> MIcal; stud{2} -> MIcont
stud = {MIcal; MIcont};  
stud = stud{str2double(answers{1})};
subj = str2double(answers{2});
blocks = [load(stud{subj}(1,:)); load(stud{subj}(2,:))]; 

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
        data_hp_lp_norm = zeros(17,10000,size(RestOnset,1));

        % motor img - calib
        for i = 1:size(RestOnset,1)
            megR3d(:,:,i) = meg(:,(RestOnset(i):MoveOnset(i)-1));
            megM3d(:,:,i) = meg(:,(MoveOnset(i):MoveOnset(i)+4999));
            data_hp_lp_norm(:,:,i) = [megR3d(:,:,i) megM3d(:,:,i)];
        end
    %     filtdat = filterFGx(data, srate, handles.slider_freq.Value, handles.slider_width.Value);
    end

% BREAK UP MOVE AND REST THEN SEPARATELY DO FFT (X), SQ(X), BP (feature)(X), NORMALIZE (feature scaling)(X), 
% GROUP FOR 5-FOLD CV... 2 Groups, Move vs Rest. BREAK UP DATA INTO TIME
% WINDOWS, PASS TIME WIN INTO 5 CELLS FOR EACH TASK
% DOC CLASSIFY 

RestFeat = abs(fft(meg).^2) ;
MoveFeat = abs(fft(meg).^2);

% Extract Features such as alpha and beta band power
    %BP Filter alpha and beta band - butterworth, filtfilt
    % abs(x)^2
    % Average across trials mean(:,:,#)
    % Create time windows of 500ms
    % Average over time windows
    
    % Alpha Band Power
    fs = 1000;
    [hp2 hp1] = butter( 3, 8/(fs/2) , 'high');
    [lp2 lp1] = butter( 3, 12.5/(fs/2)  , 'low');
    
   Arest_hp = filtfilt(hp2, hp1, RestFeat);
   Arest_hp_lp = filtfilt(lp2, lp1, Arest_hp);
   Arest_hp_lp_norm = normalize(Arest_hp_lp(:), 'range');
   Arest_hp_lp_norm = reshape(Arest_hp_lp_norm, size(meg) );
   
   Amove_hp = filtfilt(hp2, hp1, MoveFeat);
   Amove_hp_lp = filtfilt(lp2, lp1, Amove_hp);
   Amove_hp_lp_norm = normalize(Amove_hp_lp(:), 'range');
   Amove_hp_lp_norm = reshape(Amove_hp_lp_norm, size(meg) );
   
    % Beta Band Power
    fs = 1000;
    [hp2 hp1] = butter( 3, 12.5/(fs/2) , 'high');
    [lp2 lp1] = butter( 3, 30/(fs/2)  , 'low');
    
   Brest_hp = filtfilt(hp2, hp1, RestFeat);
   Brest_hp_lp = filtfilt(lp2, lp1, Brest_hp);
   Brest_hp_lp_norm = normalize(Brest_hp_lp(:), 'range');
   Brest_hp_lp_norm = reshape(Brest_hp_lp_norm, size(meg) );
   
   Bmove_hp = filtfilt(hp2, hp1, MoveFeat);
   Bmove_hp_lp = filtfilt(lp2, lp1, Bmove_hp);
   Bmove_hp_lp_norm = normalize(Bmove_hp_lp(:), 'range');
   Bmove_hp_lp_norm = reshape(Bmove_hp_lp_norm, size(meg) );
   
   
   
   % x_hp_lp_norm1 = (x_hp_lp(:) - min(x_hp_lp(:))) / (max(x_hp_lp(:)) - min(x_hp_lp(:) ));
   % x_hp_lp_norm1 = reshape(x_hp_lp_norm1, size(meg));
%%


% Normalize the results

  
   tw = 500;
   overlap = 100;
   megR3d_tw = cell((size(megR3d, 2) - tw)/overlap + 1, 1 );
   megM3d_tw = cell((size(megM3d, 2) - tw)/overlap + 1, 1 );
   
   for i = 1:(size(megR3d, 2) - tw)/overlap + 1
       if i == 1
            megR3d_tw{i} = megR3d(:, 1:tw, :);
            megM3d_tw{i} = megM3d(:, 1:tw, :);

       else
            megR3d_tw{i} = megR3d( :, (1+(overlap*(i-1) )) : (tw+(overlap*(i-1) )), :);
            megM3d_tw{i} = megM3d( :, (1+(overlap*(i-1) )) : (tw+(overlap*(i-1) )), :);
       end
   end

% Run 5-fold cross-validation with 5 groups made up of 500ms time windows


% Generate LDA Model and test it
