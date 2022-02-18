function [ featureMatrix ] = getWindowedSpectralPower(rawData, samplingRate, updateRate,windowLength,nfft,freqBand)
%getWindowSpectralPower - Windows the data and calculates the spectral
%power in the specific frequency band.  This uses the hamming window to
%help reduce spectral leakage.
%
% Input is framed as getWindowedSpectralPower(rawData, samplingRate,
%   updateRate,windowLength,nfft,freqBand)
%
% rawData - the data with columns representing channels or features, rows
% are samples.
%
% samplingRate - sampling rate of the data (samples/second)or Hz
%
% updateRate - The rate of the final feature matrix in estimates/second.
% The window shift is estimated as the (Sampling Rate) / (Update Rate).
%
% windowLength - The length of the time windows that are sliding along the
% data in samples.
%
% nfft - nfft for the FFT function.
%
% freqBand - a element vector that indicates the bounds of the frequency
% band you want to estimate (for alpha, put in [8 13]).
%%
[nSamp, nChan] = size(rawData);
    
%% Get the first and last indecies of the windows.
windowShift = samplingRate/updateRate;
window_firstIDX = 1:windowShift:(nSamp + 1 - windowLength);
window_lastIDX = windowLength:windowShift:nSamp;

numWindows = length(window_firstIDX);


%% Apply FFT on each window.

%get frequency vector.
freqVector = samplingRate*(0:(windowLength/2))/windowLength;

%Create the windowing function, and the corresponding matrix to
%multiply the data with.
hammingWindow = hamming(windowLength);
% NormFactor = (hammingWindow'*hammingWindow)*nfft;
NormFactor = nfft*2;

%Replicate the window for the number of channels.
hammingWindowMatrix = repmat(hammingWindow,1,nChan);
featureMatrix = zeros(numWindows,nChan);

for windowIDX = 1:numWindows
    X = rawData(window_firstIDX(windowIDX):...
        window_lastIDX(windowIDX),:);
    
    %Remove the linear trend
    X = detrend(X,'linear');
    
    %Multiply the data with the window function.
    X = X.*hammingWindowMatrix;

    %Calculate the FFT.
    Y = fft(X,nfft);

    %Calculate the power and normalize by power in the window function
    %and the nfft.
    Y = Y.*conj(Y)/NormFactor;
    Y = Y(1:nfft/2+1,:);
    Y(2:end-1,:) = 2*Y(2:end-1,:);
    
    %Find the indecies of frequency vector that corresponds to the
    %frequency band of interest.  
    getFreqBinIDX = (freqVector>=freqBand(1) &...
        freqVector<=freqBand(2));

    %Sum the spectral power values to get the frequeuncy band power. We'll
    %do this for all the channels.
    featureMatrix(windowIDX,:) = ...
        sum(Y(getFreqBinIDX,:));
    
end

% Pad the first values with zeros.
numSampleBlocks = floor(nSamp/samplingRate*updateRate);
numPadValues = numSampleBlocks - numWindows;
padMatrix = zeros(numPadValues,nChan);

featureMatrix = [padMatrix; featureMatrix];



end

