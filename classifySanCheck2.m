% A = 10 + 2*randn(1000, 1);
% B = -7 + 2*randn(1000,1);
% 
% 
% [ testIdx, trainIdx, nTest] = cvpartition1(size(A,1),5);
% 
% CL = cell(size(A,1) + size(B,1),1);
% for i = 1:size(A,1)
%     CL{i} = 'A';
%     CL{ i + size(A,1) } = 'B';
% end
% 
% [ CLtestIdx, CLtrainIdx, CLTest] = cvpartition1(size(CL,1),5);
% 
% for i = 1:5
%     [C{i, 1},err(i, 1),P{i, 1},logp{i, 1},coeff{i, 1}] = classify( [ A(testIdx(i,:)) ; B(testIdx(i,:))], ...
%                                                                   [A(trainIdx(i,:)) ; B(trainIdx(i,:))] , ...
%                                                                   CL(find(CLtrainIdx(i,:))), 'linear' );
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs = 1000;

t = [1/fs:1/fs:50];
% At = [1:1*fs   1+2*fs:3*fs   1+4*fs:5*fs   1+6*fs:7*fs   1+8*fs:9*fs];
% Bt = [1+1*fs:2*fs   1+3*fs:4*fs   1+5*fs:6*fs   1+7*fs:8*fs   1+9*fs:10*fs];
At = zeros(length(t)/2,1);
Bt = zeros(length(t)/2,1);

% Create epochs for tasks
for i = 1 : length(t)/(2*fs)
        % At is even to odd
        % Bt is odd to even
        At([(i-1)*fs+1:i*fs]) = ((2*i-2)*fs)+1 : (2*i-1)*fs;
        Bt([(i-1)*fs+1:i*fs]) = ((2*i-1)*fs+1) : (2*i*fs);
end
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
chan2 = 0.1*randn(length(t),1)
for i = At
    % 10 Hz Sine Wave - 3 amplitude
    chan2(i) = chan2(i)+3*sin(2*pi*10*t(i))';
end
for j = Bt
    % 3 Hz Sine Wave - 8 amplitude
    chan2(j) = chan2(j)+8*sin(2*pi*27*t(j))';
end
figure; subplot(2,1,1); plot(t,chan1); 
subplot(2,1,2); plot(t,chan2);
ylim([-12 12]);
legend('Chan1', 'Chan2');

data = [chan1 chan2];
data = data';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Time window, shift, # of pnts
       tw = 500;
       shift = 100;fs
    % number of pnts is length of data
       npnts = length(data);
    % meg time window size
       data_tw = cell((size(data, 2) - tw)/shift + 1, 2 );
    % windowed meg frequency (Short-time fft)
       data_fft_singleTW = cell((size(data, 2) - tw)/shift + 1, 1 );
    % meg low pass
       x_abs = cell((size(data, 2) - tw)/shift + 1, 1 );
       x_pow = cell((size(data, 2) - tw)/shift + 1, 1 );
       a = 'alpha';
       a_val = cell((size(data, 2) - tw)/shift + 1, 1 );
       b = 'beta';
       b_val = cell((size(data, 2) - tw)/shift + 1, 1 );
       d = 'data';
       d_val = cell( ((size(data, 2) - tw)/shift + 1)*2, 1 );
       
       group = cell( size(data, 2), 1 );
       
       % create a structure of meg-> time smoothed 0.5 sec -> normalized
       % all data; alpha band power features; beta band power features
       x_pow_str = struct('alpha', {a_val}, 'beta', {b_val}, 'data', {d_val});

       feat = zeros((size(data, 1) - tw)/shift + 1, 2);
       alpha = 8:13;
       beta = 25:30;      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
       
    for i = 1:( (size(data, 1) - tw)/shift + 1 )
           % If i == 1, save meg_tw{i,1} as the 1st 500 time points
           % else shift the 500 point time window every 100 points 
           % (Problem?) take fft of each time window ( STFFT )
           % normalize data
           % seperate alpha and beta data
           % add class labels to alpha and beta data
               if i == 1
                   % store meg data in meg_tw, needs to iterate with i
                    data_tw{i,1} =data(:, (1:tw) );
                    data_tw{i,2} = labels(tw);
               else
                    data_tw{i,1} = data( :, (1+(shift*(i-1) )) : (tw+(shift*(i-1) )), :);
                    data_tw{i,2} = labels(tw+(shift*(i-1) ) );
               end

               % CHECK IF ANY LABELS ARE 0 for MOVE
                if data_tw{i,2} == 1
                     data_tw{i,2} = 'Rest';
                elseif data_tw{i,2} == 2
                    data_tw{i,2} = 'Move';
                end 
               % take fft of the data / 500 point time window
                data_fft_singleTW = (fft(data_tw{i}, nfft))/tw; 
               % take absolute value, find amplitude, Scale by 1 million or
               % the power will be too small for matlab (10^-18)
               x_abs = abs(data_fft_singleTW);
               x_pow = x_abs(1:nfft/2 + 1);
               x_pow(2:(end - 1)) = (2*1e09*x_pow(2:end-1)).^2;
               
               % separate alpha and beta data
               x_pow_str.alpha{i, 1} =  sum(x_pow(alpha)) ;
               x_pow_str.beta{i, 1} = sum(x_pow(beta));
               x_pow_str.data{2*i-1,1} =  x_pow_str.alpha{i, 1};
               x_pow_str.data{2*i,1} =  x_pow_str.beta{i, 1};


    % Group is the class labels for each alpha & beta value in all time
    % windows
              group{2*i-1, 1} = data_tw{i,2};
              group{2*i, 1} = data_tw{i,2};
    end

[ testIdx, trainIdx, nTest] = cvpartition1(size(data',1),5);

CL = cell(size(data,2),1);
% for i = 1:size(CL,1)
%     
% end

[ CLtestIdx, CLtrainIdx, CLTest] = cvpartition1(size(CL,1),5);

for i = 1:5
    [C{i, 1},err(i, 1),P{i, 1},logp{i, 1},coeff{i, 1}] = classify( data(testIdx(i,:))', ...
                                                                  data(trainIdx(i,:))' , ...
                                                                  CL( find(trainIdx(i,:) ) ), 'linear' );
end






