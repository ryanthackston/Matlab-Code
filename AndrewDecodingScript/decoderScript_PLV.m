%% Info on what calibration file you want to load.
subjectNo = 3;
conditionTag = 'imag'; %either 'real' or 'imag';
numBlocks = 1;

filepath = ['C:\Users\ryanr\Downloads\My Matlab Functions\Calibration-20210721T005419Z-001\Calibration\'];

%example of a filename 'S2_Block1_real_calib_NO.mat';

%% Parameters
samplingRate = 1000;    %sampling rate of the raw data.
updateRate = 10;        %"sampling" rate of the features.  Since they're based on windows that shift 100 ms, we have a rate of 10 Hz.
windowLength = 1000;    % length of the time window in samples.
nfft = 1024;            %nfft for the FFT used in he spectral power estimation.
alphaFreqBand = [8 13]; %Bounds of the frequency band - alpha and beta.
betaFreqBand = [13 30];
freqWidth = 2.5
windowedClassSampleType = 'last';   %parameter used for windowing the class labels.  We pick the last sample.

%% Load each block and extract the features.
%For each block, load the data and extract the spectral features and the
%class labels that correspond to each spectral power sample.

featureMatrixBlocks = cell(numBlocks,1);
classVectorBlocks = cell(numBlocks,1);
MoveOnsetBlocks = cell(numBlocks,1);
RestOnsetBlocks = cell(numBlocks,1);
for blockIDX = 1:numBlocks
    
    %Get the filename that I want to load.
    filename = ['S' int2str(subjectNo) ...
        '_Block' int2str(blockIDX) ...
        '_' conditionTag ...
        '_calib_NO.mat'];
    
    %Load the data. X
    load([filepath filename]);
    
    %PLV Matrix
    plvMatrix = GetPLV(meg,samplingRate,updateRate,windowLength,nfft,alphaFreqBand, freqWidth)
    
    %Get the windowed class labels.X
    windowedClassVector = getWindowedClass(stimulusCode,...
        samplingRate, updateRate,windowLength,windowedClassSampleType);
    
    %Stack the alpha and beta band features side by side.
    featureMatrixBlocks{blockIDX} = [plvMatrix];
    
    %Extract the class labels so that they match the same rate as
    %the spectral power features.
    classVectorBlocks{blockIDX} = windowedClassVector;
    
    %Convert the sample index for the event markers to correspond to the
    %feature index.  We're moving from raw data that has a sampling rate of
    %1000 Hz to features that vary at a sampling rate of 10 Hz.
    
%     %Save the onset times for the move and rest periods.
%     MoveOnsetBlocks{blockIDX} = floor(MoveOnset/100);
%     RestOnsetBlocks{blockIDX} = floor(RestOnset/100);
end


%% Organize the block data so they're  singular matrices.

%We can use cell2mat to append the two blocks together.
featureMatrix = cell2mat(featureMatrixBlocks);
classVector = cell2mat(classVectorBlocks);

%For the event markers, we need to add the first block's number of samples
%so the second block's event markers are aligned correctly.
if numBlocks == 2
    numSampFirstBlock = length(featureMatrixBlocks{numBlocks});
    MoveOnsetBlocks{numBlocks} = MoveOnsetBlocks{numBlocks} + numSampFirstBlock;
    RestOnsetBlocks{numBlocks} = RestOnsetBlocks{numBlocks} + numSampFirstBlock;
end

MoveOnset = cell2mat(MoveOnsetBlocks);
RestOnset = cell2mat(RestOnsetBlocks);

%We'll consolidate both as another vector that just indicates the start
%time for each cue.
trialStartTimes = sort([RestOnset; MoveOnset]);
numTrials = length(trialStartTimes);

%% Segment the data into trials.

%Orgnaize the data as cell arrays.
featureMatrixTrials = cell(numTrials,1);
classVectorTrials = cell(numTrials,1);

for trialIDX = 1:numTrials
    %We're interest in the getting data 1 to 8 seconds after the start of
    %the cue. Note the sampling rate of the features and classes are in 10
    %Hz.
    getSamples = trialStartTimes(trialIDX)+10:...
        trialStartTimes(trialIDX)+80;
    featureMatrixTrials{trialIDX} = ...
        featureMatrix(getSamples,:);
    classVectorTrials{trialIDX} = ...
        classVector(getSamples,:);
end

%% Run 10-fold cross validation.
numFeaturesToKeep = 10;
numFolds = 10;

%We get the trial indecies that correspond to training and testing sets for
%each fold.
[testIDX, trainIDX] = cvpartition1(numTrials,numFolds);

classPredictions_folds = cell(numFolds,1);
accuracies_folds = zeros(numFolds,1);

for foldIDX = 1:numFolds
    
    %Create the training set. We'll concatenate the trials that belong to
    %the training set.
    featureTraining = cell2mat(featureMatrixTrials(trainIDX(foldIDX,:)));
    classTraining = cell2mat(classVectorTrials(trainIDX(foldIDX,:)));
    
    %-------------The Mahal Ranking part -----------------------------
    %For the training matrix, rank the best ones based on the mahalonobis
    %distance.
    [featureRank] = mahalRank(featureTraining,classTraining);
    %featuresToUse is the list columns we'll use for classification.
    featuresToUse = featureRank(1:numFeaturesToKeep);
    
    %Create the testing set
    featureTesting = cell2mat(featureMatrixTrials(testIDX(foldIDX,:)));
    classTesting = cell2mat(classVectorTrials(testIDX(foldIDX,:)));

    %Only keep the top features that had were ranked.  The number we keep
    %is in 'featuresToUse', which depends on a manual parameter called
    %'numFeaturesToKeep');
    featureTraining_best = featureTraining(:,featuresToUse);
    featureTesting_best = featureTesting(:,featuresToUse);
    %-------------------------------------------------------------------
    
    %Use classify, which generates predicted classes in the testing set,
    %based on the training set of features and class labels.
    predictedClass = classify(featureTesting_best,featureTraining_best,...
        classTraining);
    
    %Now we'll see how the perdictedClass matches up with the actual class
    %labels.
    correctPredictionIDX = predictedClass == classTesting;
    
    %Calculate accuracy.
    acc = sum(correctPredictionIDX) ./ length(classTesting);
    
    %Save the class predictions and accuracy for each fold.
    classPredictions_folds{foldIDX} = predictedClass;
    accuracies_folds(foldIDX) = acc;
end

%%

disp(['median accuracies across all folds: ' num2str(median(accuracies_folds))]);
