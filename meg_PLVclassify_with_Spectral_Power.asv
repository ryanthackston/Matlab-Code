function [C err P logp coeff PredictedVsActualClasses testAcc_LDA_PLV_PSD] = meg_PLVclassify_with_Spectral_Power(Class_Vec, PLV_features, featureTesting_best, featureTraining_best, featureMatrixTrials)


    folds = 5;
    
    PLV_features2 = PLV_features(  size(PLV_features,1)-size(cell2mat(featureMatrixTrials),1)+1  :  size(PLV_features, 1),:);
   Class_Vec2 = Class_Vec(1:size(cell2mat(featureMatrixTrials),1));
   
  [ testIdx, trainIdx, nTest] = cvpartition1(size(PLV_features2,1), folds);
  


    
    if (size(PLV_features(testIdx(1,:),:), 1)) <= (size( featureTesting_best, 1))
      featureTesting_best2=  featureTesting_best(1:size(PLV_features(testIdx(1,:),:), 1));       
       featureTraining_best2 =  featureTraining_best(1:size(PLV_features(trainIdx(1,:),:), 1));    
   else
       norm_test = normalize(featureTesting_best, 'range');
       norm_train = normalize(featureTraining_best, 'range');
       norm_PLV = normalize(PLV_features, 'range');
       
       
       
        
    end

    
   C = cell(folds,1); err = zeros(folds,1); P = cell(folds,1); logp = cell(folds,1); coeff = cell(folds,1);
       
%    Class_Vec3D = cat(3,Class_Vec2, Class_Vec2);
%    Class_Vec3D = permute(Class_Vec3D,[1,3,2]);
   
        for ii = 1:folds
            % windows X type_feature x trials
            train = cat(3, PLV_features2(trainIdx(ii,:),:), norm_train);
            train = permute(train,[1,3,2]);
            test = cat(3, PLV_features2(testIdx(ii,:),:), norm_test);
            test = permute(test,[1,3,2]);
            
            
             [C{ii, 1},err(ii, 1),P{ii, 1},logp{ii, 1},coeff{ii, 1}] = ...
                                                                    classify( test(:,:,ii), ...
                                                                    train(:,:,ii) , ...
                                                                    Class_Vec2(trainIdx(ii,:),:), 'linear' );
                                                                
            PredictedVsActualClasses = strcmp(C{ii}, Class_Vec2(testIdx(ii,:)));
            
           testAcc_LDA_PLV_PSD(ii) = sum(PredictedVsActualClasses)/length(PredictedVsActualClasses);       

%

%              [C{ii, 1},err(ii, 1),P{ii, 1},logp{ii, 1},coeff{ii, 1}] = ...
%                                                                     classify( PLV_features(testIdx(ii,:),:), ...
%                                                                     PLV_features(trainIdx(ii,:),:) , ...
%                                                                     Class_Vec(trainIdx(ii,:),:), 'linear' );
%                                                                 
%             PredictedVsActualClasses = strcmp(C{ii}, Class_Vec(testIdx(ii,:)));
%             
%             testAcc(ii) = sum(PredictedVsActualClasses)/length(PredictedVsActualClasses);          
                                                                
        end
end

