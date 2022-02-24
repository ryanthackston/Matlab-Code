function [meg_tw, Fmeg, meg_f, meg_fft_singleTW, x_abs, x_pow, x_pow_str, feat, alpha, beta, w, freq_val, freq_width, PLV, d, B, I, d_sort, d_s, PLV_Rest_I, PLV_Move_I] =meg_initvars(answers, meg,  fs,nfft, freq, chans, trials, Fmeg)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
  % Sampling Frequency and filter coefficients


        
        
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
       meg_tw = cell(trials, 2 );
       
    % Filtered meg time window data
       Fmeg = cell(trials, 2 );
    % windowed meg frequency (Short-time fft)
       meg_f = cell(trials, 1 );
       
       meg_fft_singleTW = zeros(chans, nfft);
       
       
       x_abs = zeros(chans, nfft);
       x_pow = zeros(chans, (nfft/2 + 1) );
       
       
    % meg low pass
%        x_abs = cell(trials, 1 );
%        x_pow = cell(trials, 1 );
       a = 'alpha';
       a_val = cell(trials, 1 );
       b = 'beta';
       b_val = cell(trials, 1 );
       d = 'data';
       d_val = cell( (trials), 1 );
       
       % create a structure of meg-> time smoothed 0.5 sec -> normalized
       % all data; alpha band power features; beta band power features
       x_pow_str = struct('alpha', {a_val}, 'beta', {b_val}, 'data', {d_val});

    %     x_lp_str_medavg = cell(trials, 2 );
       feat = zeros(trials, 2);
       alpha = 8:13;
       beta = 25:30;
       % bandpass = length(alpha) + length(beta);
       % initialize group label cells - default val - size(group) = 4996 x 1
       
       w = hamming(tw);
       srate = 1000;
       freq_val = str2num(answers{4});
       freq_width = str2num(answers{6});
       
       PLV = cell( trials, 2);
       d = cell( trials, 1);
       B = cell( trials, 1);
       I = cell( trials, 1);
       d_sort = cell( trials, 1);
       d_s = zeros(trials,10);
       
       PLV_Rest_I = zeros( trials,1);
       PLV_Move_I = zeros(trials,1);
       
end

