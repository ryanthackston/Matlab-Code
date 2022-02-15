function [PLV_Mahal_Sort, PLV_Mahal_Coord, PLV_Diff_Coord, PLV_Diff_Sort, trials_compared] = meg_PLVfeatures(top_features, PLV_Move_I, PLV_Rest_I, PLV_cut, row, col)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% If there are less Rest PLV indices than Move,

    PLV_Rest_I = PLV_Rest_I';
    PLV_Move_I = PLV_Move_I';

    if size(PLV_Rest_I,1) <= size(PLV_Move_I,1)
        trials_compared = size(PLV_Rest_I,1);  
        % then cut out extra PLV Move Indices and create Mahal and Difference
        % variables.
        PLV_Move_I = PLV_Move_I( 1:trials_compared ) ;

        PLV_Rest = cell(trials_compared, 1);
        PLV_Move= cell(trials_compared, 1);
        PLV_Diff = cell(trials_compared, 1);

        PLV_Diff_Sort = cell(trials_compared, 1);
        PLV_Diff_Coord = cell(trials_compared, 2);

        PLV_Mahal_Sort = cell(trials_compared, 1);
        PLV_Mahal_Coord = cell(trials_compared, 2);

        PLV_Mahal = cell(trials_compared, 1);

        for c = 1:size(PLV_Rest_I,1)

                PLV_Rest{c,1} = PLV_cut{PLV_Rest_I(c), 1};
                PLV_Move{c,1} = PLV_cut{PLV_Move_I(c), 1};
                PLV_Diff{c} = abs(PLV_Move{c} - PLV_Rest{c});
                
                %This is the problem - need to fix copied values
                PLV_Mahal{c} = mahal(PLV_Rest{c}, PLV_Move{c} );

                % Sort out the furthest values between Rest and Move for a specific
                % trial and get the indices as well
                [tmp, PLV_Diff_Ind] = sort(PLV_Diff{c}, 'descend');
                [tmp2, PLV_Mahal_Ind] = sort(PLV_Mahal{c}, 'descend');
                
                % PLV_Diff_Sort saves the largest PLV differences between
                % rest and move trials among each of the 14 channels
                PLV_Diff_Sort{c} = tmp(1:10);
                PLV_Diff_Ind_Sort = PLV_Diff_Ind(1:10);
                % PLV_Diff_Coord 
                PLV_Diff_Coord{c,1} = [col(PLV_Diff_Ind_Sort), row(PLV_Diff_Ind_Sort)];

                %PLV_Mahal_Sort saves the 10 largest Mahal distances between equivalent Rest and Move trials 
                PLV_Mahal_Sort{c} = tmp2(1:10);
                
      % PROBLEM with PLV_Mahal_Ind - ALWAYS THE SAME VALUES!!! Max Values
      % were always these interactions!!!
                PLV_Mahal_Ind_Sort = PLV_Mahal_Ind(1:10);
                % Mahal Coord are low to high
                PLV_Mahal_Coord{c,1} = [col(PLV_Mahal_Ind_Sort), row(PLV_Mahal_Ind_Sort)];
%                 PLV_Mahal_Coord{c,2} = cell(10,1);     
                
                % PLV_Mahal_Coord - Shows the PLV_Mahal_Sort Channels Interactions
%                 for j = 1:10
%                     PLV_Mahal_Coord{c,2}{j,1} = strcat("Channels ", string(PLV_Mahal_Coord{c}(j,1)), "-", string(PLV_Mahal_Coord{c}(j,2)));
%                 end
                tmp2 = [];

        end
    

    % Else if there are more Rest indices than Move
    elseif size(PLV_Rest_I,1) > size(PLV_Move_I,1)
        trials_compared = size(PLV_Move_I,1);  
        % cut out extra PLV Rest Indices and create Mahal and Difference
        % variables.
        PLV_Rest_I = PLV_Rest_I( 1:trials_compared )

        PLV_Rest = cell(trials_compared, 1);
        PLV_Move= cell(trials_compared, 1);
        PLV_Diff = cell(trials_compared, 1);

        PLV_Diff_Sort = cell(trials_compared, 1);
        PLV_Diff_Coord = cell(trials_compared, 2);

        PLV_Mahal_Sort = cell(trials_compared, 1);
        PLV_Mahal_Coord = cell(trials_compared, 2);
        PLV_Mahal_Coord{:,2} = cell(top_features,1);

        for c = 1:size(PLV_Move_I,1)

                PLV_Rest{c} = PLV_cut{PLV_Rest_I(c), 1};
                PLV_Move{c} = PLV_cut{PLV_Move_I(c), 1};
                PLV_Diff{c} = abs(PLV_Move{c} - PLV_Rest{c});
                PLV_Mahal{c} = mahal(PLV_Rest{c}, PLV_Move{c});

                % Sort out the furthest values between Rest and Move for a specific
                % trial and get the indices as well
                [tmp, PLV_Diff_Ind] = sort(PLV_Diff{c}, 'descend');
                [tmp2, PLV_Mahal_Ind] = sort(PLV_Mahal{c}, 'descend');
                % Keep the top 10 indices
                % tril_ind is used
                PLV_Diff_Sort{c} = tmp(1:top_features);
                PLV_Diff_Ind_Sort = PLV_Diff_Ind(1:top_features);
                PLV_Diff_Coord{c,1} = [col(PLV_Diff_Ind_Sort), row(PLV_Diff_Ind_Sort)];

                PLV_Mahal_Sort{c} = tmp2(1:top_features);
                PLV_Mahal_Ind_Sort = PLV_Mahal_Ind(1:top_features);
                % Mahal Coord are low to high
                PLV_Mahal_Coord{c,1} = [col(PLV_Mahal_Ind_Sort), row(PLV_Mahal_Ind_Sort)];
                PLV_Mahal_Coord{c,2} = cell(top_features,1);            
                for j = 1:top_features
                    PLV_Mahal_Coord{c,2}{j,1} = strcat("Channels ", string(PLV_Mahal_Coord{c}(j,1)), "-", string(PLV_Mahal_Coord{c}(j,2)));
                end           
        end
        
    end

end

