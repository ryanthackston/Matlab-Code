
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs = 1000;

t = [1/fs:1/fs:50];
% At = [1:1*fs   1+2*fs:3*fs   1+4*fs:5*fs   1+6*fs:7*fs   1+8*fs:9*fs];
% Bt = [1+1*fs:2*fs   1+3*fs:4*fs   1+5*fs:6*fs   1+7*fs:8*fs   1+9*fs:10*fs];
At = zeros(length(t)/2,1);
Bt = zeros(length(t)/2,1);

%Andrew version of task epochs, I'll make each cue last 5 seconds.
%define boundaries for the second class.
firstCut = 5:10:45;
lastCut = 10:10:50;
classVector = ones(length(t),1);
for i = 1:length(firstCut)
    classVector(t>firstCut(i) & t<lastCut(i)) =2;
end
At = find(classVector==1);
Bt = find(classVector==2);


% chan1 =  zeros(length(t),1)
chan1 = 0.1*randn(length(t),1);
for i = At
    % 10 Hz Sine Wave - 8 amplitude
    chan1(i) = chan1(i)+8*sin( 2*pi*10*t(i) )' ;
    
end
for j = Bt
    % 3 Hz Sine Wave - 3 amplitude
    chan1(j) = chan1(j) + 3*sin( 2*pi*27*t(j) )';
end

% chan1 =  zeros(length(t),1)
chan2 = 0.1*randn(length(t),1);
for i = At
    % 10 Hz Sine Wave - 3 amplitude
    chan2(i) = chan2(i)+3*sin(2*pi*10*t(i))';
end
for j = Bt
    % 3 Hz Sine Wave - 8 amplitude
    chan2(j) = chan2(j)+8*sin(2*pi*27*t(j))';
end

% figure; subplot(2,1,1); plot(t,chan1); 
% subplot(2,1,2); plot(t,chan2);
% ylim([-12 12]);
% legend('Chan1', 'Chan2');
% 
data = [chan1 chan2];
data = data';
%% looks parameter setting and initialization.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Time window, shift, # of pnts
       tw = 500;
       shift = 100;
    % number of pnts is length of data
       nfft = 1028;
       npnts = length(data);
    % meg time window size
       data_tw = cell((size(data, 1) - tw)/shift + 1, 2 );
    % windowed meg frequency (Short-time fft)
       data_fft_singleTW = cell((size(data, 1) - tw)/shift + 1, 1 );
    % meg low pass
       x_abs = cell((size(data, 2) - tw)/shift + 1, 1 );
       x_pow = cell((size(data, 2) - tw)/shift + 1, 1 );
       a = 'low';
       a_val = cell((size(data, 2) - tw)/shift + 1, 1 );
       b = 'high';
       b_val = cell((size(data, 2) - tw)/shift + 1, 1 );
       d = 'data';
       d_val = cell( ((size(data, 2) - tw)/shift + 1), 1 );
       
       group = cell( ( (size(data, 2) - tw)/shift + 1 ), 1 );
       
       % create a structure of meg-> time smoothed 0.5 sec -> normalized
       % all data; alpha band power features; beta band power features
       x_pow_str = struct('low', {a_val}, 'high', {b_val}, 'data', {d_val});

       feat = zeros((size(data, 1) - tw)/shift + 1, 2);
       low = 26:28;
       high = 9:11;      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
%% Actual processing
    for i = 1:( (size(data, 2) - tw)/shift + 1 )
           % If i == 1, save meg_tw{i,1} as the 1st 500 time points
           % else shift the 500 point time window every 100 points 
           % (Problem?) take fft of each time window ( STFFT )
           % normalize data
           % seperate low and high data
           % add class labels to low and high data
               if i == 1
                   % store meg data in meg_tw, needs to iterate with i
                    data_tw{i,1} =data(:, (1:tw) );
                    data_tw{i,2} = 0;
               else
                    data_tw{i,1} = data( :, (1+(shift*(i-1) )) : (tw+(shift*(i-1) )), :);
                    data_tw{i,2} = 0;
               end

               % CHECK IF ANY LABELS ARE 0 for MOVE
               % NEED TO FIX LABELS, MAKE IF/ELSE STATEMENT FOR CLASSIFYING
               % DATA
               
             % IF MAX(CHAN 1) > MAX(CHAN 2)
%                 if max(data_tw{i,1}(1,:) ) > max( data_tw{i,1}(2,:) )
%                      data_tw{i,2} = 'Rest';
%                 elseif max(data_tw{i,1}(1,:) ) < max( data_tw{i,1}(2,:) )
%                     data_tw{i,2} = 'Move';
%                 end
                %This is okay I guess.  You might run into
                %problems if the the time window is too small and you only
                %get the negative side of a sine curve, especially if the
                %frequency is slow.
                % I would recommend hard coding the class labels like you
                % did with the sine wave.
                % We do run into a curious problem of how to manage windows
                % with both class labels.  For now we can just pick
                % whichever class comes up most frequently.
               if i == 1
                   % store meg data in meg_tw, needs to iterate with i
                    data_tw{i,2} = 'Rest';
               else
                    class_singleTW = classVector((1+(shift*(i-1) )) : (tw+(shift*(i-1) )));
                    mostFreqClass = mode(class_singleTW);
                    if mostFreqClass ==1
                        data_tw{i,2} = 'Rest';
                    else
                        data_tw{i,2} = 'Move';
                    end
               end
                
                
               % take fft of the data / 500 point time window
                data_fft_singleTW1 = (fft(data_tw{i}(1,:), nfft))/tw; 
                data_fft_singleTW2 = (fft(data_tw{i}(2,:), nfft))/tw; 
                
               % take absolute value, find amplitude, Scale by 1 million or
               % the power will be too small for matlab (10^-18)
               x_abs1 = abs(data_fft_singleTW1);
               x_abs2 = abs(data_fft_singleTW2);
               
               x_pow1 = x_abs1(1:nfft/2 + 1);
               x_pow2 = x_abs2(1:nfft/2 + 1);

%                x_pow1(2:(end - 1)) = (2*1e09*x_pow1(2:end-1)).^2;
%                x_pow2(2:(end - 1)) = (2*1e09*x_pow2(2:end-1)).^2;
               %omitting the scaling part
              x_pow1(2:(end - 1)) = (2*x_pow1(2:end-1)).^2;
               x_pow2(2:(end - 1)) = (2*x_pow2(2:end-1)).^2;
               

%                % separate low and high data
%                % low - freq 1:5
%                % high - freq 8:12
%                x_pow_str.low{i, 1} =  [max(x_pow1(low)); max(x_pow2(low))] ;
%                x_pow_str.high{i, 1} = [max(x_pow1(high)); max(x_pow2(high))];
%                x_pow_str.data{i,1} = [x_pow_str.low{i, 1}; x_pow_str.high{i, 1} ] ;
% %                x_pow_str.data{2*i-1,1} =  x_pow_str.low{i, 1};
% %                x_pow_str.data{2*i,1} =  x_pow_str.high{i, 1};

               %Let's get the band power as usual.  We sum across the
               %2 frequency bands of interest: low being 1-5 Hz, 
               %    high being 8-12Hz
               
               numberPSDbins = nfft/2;
               freqVector = 0:fs/numberPSDbins:fs;
               lowfreq =  find(freqVector>1 & freqVector<5);
               highfreq  = find(freqVector>8 & freqVector<12);
               x_pow_str.low{i, 1} =  [sum(x_pow1(lowfreq)); sum(x_pow2(lowfreq))] ;
               x_pow_str.high{i, 1} = [sum(x_pow1(highfreq)); sum(x_pow2(highfreq))];
               x_pow_str.data{i,1} = [x_pow_str.low{i, 1}; x_pow_str.high{i, 1} ] ;


    % Group is the class labels for each low & high value in all time
    % windows
%               group{i} = data_tw{i,2};
    end
    %%
lowf=[]; highf = [];
% Store max values in frequencies of interest
% lowf

for i = 1:size(x_pow_str.low,1)
    % Low frequency is rest
    lowf = [lowf x_pow_str.low{i,1}];
    % High frequency is move
    highf = [highf x_pow_str.high{i,1}];
end


figure; subplot(3,1,1); plot(t,chan1); 
legend('Chan1');

subplot(3,1,2); plot(t,chan2);
ylim([-12 12]);
legend('Chan2');

tz = [ 1/(length(lowf)*2) : 50 / (length(lowf)*2) : 50];
subplot(3,1,3);
plot(tz, log(lowf(:))); 
hold on; 
plot(tz, log(highf(:)));
xlabel('Time Windows (100ms)')
ylabel('Spectral Power');
legend('Rest Freq', 'Move Freq');


% CL = cell(size(data,2), 1);

maxdata = [];
maxdata = [ maxdata x_pow_str.data{:} ];
% maxdata = maxdata(:); 
%I think this is a point you might be missing.  For each class prediction.,
%we're simultaneously looking at one sample from
% 4 features (2 channals, 2 power bands).
%Classifiers are design to reorganize the data into a hyperplane.  In this
%case, it's a 4-dimensional plane.
%When you stack them one on top another, you're framing the problem of
%trying to predict a class with only 1 feature at a time, and you're
%running that problem 4 times with a different feature.

[ testIdx, trainIdx, nTest] = cvpartition1(size(maxdata,2),5);

% Using maxdata for classifier
CL = cell(size(maxdata,1), 1);



% % FILL IN CL!!!!!!!!!!!!!
% % EACH DATA POINT CORRESPONDS TO GROUP LABEL
% % Based on it being
% for i = 1:size(data_tw,1)
%        CL{ 4*i-3}  = data_tw{ i,2 };
%        CL{ 4*i-2}  = data_tw{ i,2 };
%        CL{ 4*i-1}  = data_tw{ i,2 };
%        CL{ 4*i}  = data_tw{ i,2 };
%     
% end
%Again this part is incorrect. The classifier looks at 4 features in a
%hyperplane.


% So we don't need to replicate the class label 4 times (or x times
% depending on the number of features.)
CL = data_tw(:,2);



%%
[ CLtestIdx, CLtrainIdx, CLTest] = cvpartition1(size(CL,1),5);


%maxdata needs to be transposed since the classify function strictly
%treats each row as a sample
maxdata = maxdata';

% 2 Channels of data -> Time window - 500 window & 100 shift -> Assign Rest/Move based on Max
% amplitude in chan1 or chan2 -> find fft power -> extract high and low
% frequencies
% for i = 1:5
%     [C{i, 1},err(i, 1),P{i, 1},logp{i, 1},coeff{i, 1}] = classify( maxdata(testIdx(i,:)), ...
%                                                                   maxdata(trainIdx(i,:)) , ...
%                                                                   CL(find(trainIdx(i,:))), 'linear' );
% end


for i = 1:5
  [C{i, 1},err(i, 1),P{i, 1},logp{i, 1},coeff{i, 1}] = ...
      classify( maxdata(testIdx(i,:),:), ...
      maxdata(trainIdx(i,:),:) , ...
      CL(find(trainIdx(i,:))), 'linear' );
end






