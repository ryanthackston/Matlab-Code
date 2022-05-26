function [ testAcc, testAcc_rand] = meg_SVM(PLV_features, Class_Vec, folds, testIdx)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
X = PLV_features;
Y = Class_Vec;

% # of Support Vectors not the same as # of data points
SVMModel = fitcsvm(X,Y,'Standardize',true,'KernelFunction','RBF',...
'KernelScale','auto');

%classdef DisallowVectorOps - [varargout] = subsref(this,s) - finds
%NumObservations in SVMModel - "this" is ClassificationsSVM class - "s"
%is struct of NumObservations
%Returns 4496 in a cell as varargin

%classdef cvpartition -> cv = cvpartition(varargin) 
% varargin is now 1x3 cell {4496, 'k', 5)
% converts all string to chars
% After if-else goes to Regular Constructor

%cv.Impl = internal.stats.cvpartitionInMemoryImpl(varargin{:});
% Sets up the cvpartitionInMemoryImpl class. Gives type, NumTestSets,
% TrainSize, TestSize, NumObservations.
% Calculates amount of data points in each kFold for Test and Training
    cv = cvpartition(SVMModel.NumObservations, 'k', 5 );
    
 

    %Hidden functions, can't un-randomize
    %varargin is 1x2 cell - {'cvpartition', 1x1cvpartition class}
    % idxBaseArg should be empty
    % modelParams = this.ModelParams % Gives the model parameters of the
    % SVM model
    % modelParams.VerbosityLevel = 0;
    
    % Template fit from model parameters - 1x1 fitTemplate class
    %  temp = classreg.learning.FitTemplate.make(this.ModelParams.Method,...
    %            'type','classification','scoretransform',this.PrivScoreTransform,...
    %            'modelparams',modelParams,'CrossVal','on',varargin{:});
    % cross_validated_model = crossval(SVMModel,'cvpartition',cv);
    
    % Fitting SVM Model data with model parameters - 1x1 Classificantion class
    %  partModel = fit(temp,this.X,this.Y,'Weights',this.W,...
    %            'predictornames',this.DataSummary.PredictorNames,...
    %            'responsename',this.ResponseName,...
    %            'classnames',this.ClassNames,'cost',this.Cost,'prior',this.Prior);
    %        partModel.ScoreType = this.ScoreType;
    
    
    % cvpartition.m line 160            n = this.Impl.NumObservations;
        % classdef cvpartitionImpl.m line 121       methods
                                                                          % function n = get.NumObservations(this)
                                                                           %              n = this.N;
                                                                           %   end
 
       cross_validated_model = crossval(SVMModel,'cvpartition',cv);
    
  for i = 1:folds
    Predictions = predict(cross_validated_model.Trained{1}, X(testIdx(i,:),:) ) ;
%     Results = confusionmat(cross_validated_model.Y(X(testIdx(i,:),:), Predictions);

     PredictedVsActualClasses = strcmp(Predictions, Class_Vec(testIdx(i,:),:) );

 
     testAcc(i) = sum(PredictedVsActualClasses)/length(PredictedVsActualClasses);
     
     % - median Acc 84% with randomized training model and time-sequenced test data
          % - 76% randomized data - with Subj 13, Day 2, 14Hz, Freq Width 5, Start at 0.5 Sec
          
      Predictions_rand = predict(cross_validated_model.Trained{1}, X(test(cv, i), :) ) ;
      PredictedVsActualClasses_rand = strcmp(Predictions_rand, Class_Vec(test(cv, i)) );
      testAcc_rand(i) = sum(PredictedVsActualClasses_rand)/length(PredictedVsActualClasses_rand);
  end
    
end

