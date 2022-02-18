function [windowedClassVector] = getWindowedClass(classVector,...
    samplingRate, updateRate,windowLength,sampleType)
%getWindowedClass - To match the windowed estimates for other feature
%extraction methods, this function extracts the class label from the same
%sliding windows.
% Input is getWindowedClass(classVector,
% samplingRate,updateRate,windowLength,type)
%
% Output is a volumn vector of class labels that are undersampled to match the
% windowing update rate.
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

% 'sampleType' - Choose 'first','mode',or 'last'
%   'first' extracts the first sample
%   'mode' extracts the most common class label in the window.
%   'last' extracts the last sample in the window.


%% Get the first and last indecies of the windows.

nSamp = length(classVector);
windowShift = samplingRate/updateRate;
window_firstIDX = 1:windowShift:(nSamp + 1 - windowLength);
window_lastIDX = windowLength:windowShift:nSamp;

numWindows = length(window_firstIDX);

%% For each window, run the routine.

windowedClassVector = zeros(numWindows,1);

for windowIDX = 1:numWindows
    
    %Get the class label values for the window.
    windowClassData = classVector(window_firstIDX(windowIDX):...
        window_lastIDX(windowIDX));
    
    switch sampleType
        case 'first'
            windowedClassVector(windowIDX) = windowClassData(1);
        case 'mode'
            windowedClassVector(windowIDX) = mode(windowClassData);
        case 'last'
            windowedClassVector(windowIDX) = windowClassData(end);
    end
end


% Pad the first values with the first class value.
numSampleBlocks = floor(nSamp/samplingRate*updateRate);
numPadValues = numSampleBlocks - numWindows;
padMatrix = zeros(numPadValues,1);

windowedClassVector = [padMatrix; windowedClassVector];

end

