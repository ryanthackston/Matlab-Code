

% Make synthetic dataset

% Try feeding a 10Hz Sine Wave 1 channel and 2


% Rest if Channel 1 Amp is High & Channel 2 Amp is Low
% Move if Channel 1 Amp is Low & Channel 2 Amp is High

% Time Window -> FFT (TW data, nfft) / (TW Length) -> Abs Value --
% -> abs data ( (1:nfft/2) + 1) -> data(2:(end-1)) = (scale(data(2:end-1)))^2)

% You get power of certain frequencies in given time window 
% if Max Freq >threshold, then action 1
% elseif Max Freq < threshold then action 2






%%%%%%%%%%%%%%%%%%%%%%%%%%

% Feature Selection - Mahalanobis Distance
% Include only important features related to rest vs. movement. How do we
% figure out # of features?

% Mahalanobis - Calculates distance from 1 data point to cluster. Basically
% a z-score of one point to a cluster. A high z-score is easy to
% discriminate.

% Rank which features are best for each training set. 
% Do this with Mahalanobis & use 5 best features

% To accurately tell # of features to use, need to run another
% cross-validation set.




























