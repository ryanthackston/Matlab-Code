function [featureRankList] = mahalRank(featureMatrix, classVector)
%mahalRank - Returns rank of the features based on how separable the
%datapoints are between the two classes.
%   

%This should be 1 or 2.. but just in case we'll get the values.
classValues = unique(classVector);

if length(classValues)>2
    disp('WARNING: 3 classes found in the ranking');
end

featureValues_classOne = featureMatrix(classVector==classValues(1),:);
featureValues_classTwo = featureMatrix(classVector==classValues(2),:);


[~,numFeatures] = size(featureMatrix);

mahalDistances_median = zeros(numFeatures,1);
for featureIDX = 1:numFeatures
    mahalDistances = ...
        mahal(featureValues_classOne(:,featureIDX),...
        featureValues_classTwo(:,featureIDX));
    mahalDistances_median(featureIDX) = median(mahalDistances);
end


[rankedMahalDistances,featureRankList] = sort(mahalDistances_median,'descend');

%take off the sensors with "NaNs" since they likely had 0's.
nanFeatures = isnan(rankedMahalDistances);
featureRankList(nanFeatures) = [];


disp(featureRankList);



end

