function [ testIdx, trainIdx, nTest] = cvpartition1(nTr,kFold)
%CVPARTITION1 creates boolian indecies of training and testing trials for
%k-fold cross validation.
%
%Inputs:
%nTr = number of trials.
%kFold = number of folds.
%
%
%Outputs:
%testIdx = Matrix of indecies for trials that belong to a testing trials for
%a particular fold. Folds correspond to each row.
%trainIdx - Matrix of indecies for trials that belong to a training folds for
%a particular fold. Folds correspond to each row.
%nTest - number of trials used for testing for a particular fold.
%

testIdx = false(kFold,nTr); 
nTest = floor(nTr/kFold)*ones(kFold,1);
if mod(nTr,kFold)
    for fold = 1:mod(nTr,kFold)
        nTest(fold) = nTest(fold)+1;
    end
end
nTest2 = [0; cumsum(nTest)];
for fold = 1:kFold
    testIdx(fold,nTest2(fold)+1:nTest2(fold+1)) = true;
end
trainIdx = ~testIdx;

end

