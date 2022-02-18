function [Mahal_Vec, Class_Vec] = meg_featuresSorted(trials_compared, top_features, PLV_Mahal_Sort, PLV_Mahal_Coord)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    Mahal_Vec = zeros(trials_compared*10, 3);
    Class_Vec = [];
    for i = 1:trials_compared
        %Sorted Mahal Distances - All PLV_Mahal_Sort cell values in one column
        Mahal_Vec( (i*10-10+1): i*10, 1 ) = PLV_Mahal_Sort{i};
        %Column Coord Channels
        Mahal_Vec( (i*10-10+1): i*10, 2 ) = PLV_Mahal_Coord{i}(1:10, 1);
        %Row Coord Channels
        Mahal_Vec( (i*10-10+1): i*10, 3 ) = PLV_Mahal_Coord{i}(1:10, 2);

        % Cell vector of Channel Interaction Classes
        for j = 1:top_features
            Class_Vec{j+(i-1)*10} = convertStringsToChars(PLV_Mahal_Coord{i,2}{j,1});
        end

    end

end

