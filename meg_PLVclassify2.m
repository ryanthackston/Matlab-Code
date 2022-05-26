function [C err P logp coeff PredictedVsActualClasses testAcc] = meg_PLVclassify2(Class_Vec, PLV_features, featureTesting_best, featureTraining_best)
   

    folds = 5;
    [ testIdx, trainIdx, nTest] = cvpartition1(size(PLV_features,1), folds);

   C = cell(folds,1); err = zeros(folds,1); P = cell(folds,1); logp = cell(folds,1); coeff = cell(folds,1);
       
       
        for i = 1:folds
            
             [C{i, 1},err(i, 1),P{i, 1},logp{i, 1},coeff{i, 1}] = ...
                                                                    classify( PLV_features(testIdx(i,:),:), ...
                                                                    PLV_features(trainIdx(i,:),:) , ...
                                                                    Class_Vec(trainIdx(i,:),:), 'linear' );
                                                                
            PredictedVsActualClasses = strcmp(C{i}, Class_Vec(testIdx(i,:)));
            
            testAcc(i) = sum(PredictedVsActualClasses)/length(PredictedVsActualClasses);                                                    
                                                                
        end
end

