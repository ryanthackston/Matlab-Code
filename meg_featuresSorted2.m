function [Class_Vec] = meg_featuresSorted2(PLV, trials)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    Class_Vec = cell(trials,1)
    for i = 1:trials
        Class_Vec{i} = PLV{i,2};
    end

end

