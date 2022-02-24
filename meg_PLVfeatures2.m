function [PLV_features, PLV_Mahal, Median_Mahal, Median_Mahal_Sort, Median_Mahal_Ind, Top_Median_Mahal_Sort, Top_Median_Mahal_Ind, PLV_Diff_Coord, PLV_Diff_Sort, trials_compared, PLV_Move_I, PLV_Rest_I] = meg_PLVfeatures2(top_features, trials, PLV, PLV_Move_I, PLV_Rest_I, PLV_cut, row, col, tril_ind)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% If there are less Rest PLV indices than Move,

%     PLV_Rest_I = PLV_Rest_I';
%     PLV_Move_I = PLV_Move_I';

    if size(PLV_Rest_I,1) <= size(PLV_Move_I,1)
        trials_compared = size(PLV_Rest_I,1);  
        % then cut out extra PLV Move Indices and create Mahal and Difference
        % variables.
        PLV_Move_I = PLV_Move_I( 1:trials_compared ) ;

        PLV_Rest = zeros(trials_compared, size(PLV_cut,2) );
        PLV_Move= zeros(trials_compared, size(PLV_cut,2) );
        PLV_Diff = zeros(trials_compared, size(PLV_cut,2) );

        PLV_Diff_Sort = zeros(trials_compared, size(PLV_cut,2) );
        PLV_Diff_Coord = zeros(trials_compared, size(PLV_cut,2) );
        
        PLV_Mahal = zeros(trials_compared, size(PLV_cut,2) );
        Median_Mahal = zeros(size(PLV_cut,2), 1 );
        Avg_Mahal = zeros(size(PLV_cut,2), 1 );

%         PLV_Mahal_Sort = zeros(trials_compared, size(PLV_cut,2) );
%         PLV_Mahal_Coord = zeros(trials_compared, size(PLV_cut,2) );

        for chan_chan = 1:size(PLV_cut,2)

                PLV_Rest(:, chan_chan) = PLV_cut([PLV_Rest_I], chan_chan);
                PLV_Move(:, chan_chan) = PLV_cut([PLV_Move_I], chan_chan);
                PLV_Diff(:, chan_chan) = abs( PLV_Move(:, chan_chan)  -  PLV_Rest(:, chan_chan)) ;
                
                % PLV_Mahal - 
                PLV_Mahal(:, chan_chan) = mahal(PLV_Rest(:, chan_chan) ,  PLV_Move(:, chan_chan) );
                Median_Mahal(chan_chan) = median(PLV_Mahal(:, chan_chan));
                Avg_Mahal(chan_chan) = mean(PLV_Mahal(:, chan_chan));                
%                 PLV_Mahal_Sort{chan_chan}
        end

                %PLV_Mahal_Sort saves the 10 largest Mahal distances between equivalent Rest and Move trials 
                [Median_Mahal_Sort Median_Mahal_Ind] = sort(Median_Mahal, 'desc');
                [Avg_Mahal_Sort Avg_Mahal_Ind] = sort(Avg_Mahal, 'desc');
                
                Top_Median_Mahal_Sort = Median_Mahal_Sort(1:top_features)
                Top_Median_Mahal_Ind = Median_Mahal_Ind(1:top_features)
                
                % Mahal Coord are low to high
               Top_Median_Mahal_Chan = [col(Top_Median_Mahal_Ind), row(Top_Median_Mahal_Ind)];

               PLV_features = zeros(trials, top_features);
               
               for T = 1:trials
                   PLV_features(T,:) = PLV{T}(tril_ind(Top_Median_Mahal_Ind));
               end

    

    % Else if there are more Rest indices than Move
    elseif size(PLV_Rest_I,1) > size(PLV_Move_I,1)
        trials_compared = size(PLV_Move_I,1);  
        % then cut out extra PLV Move Indices and create Mahal and Difference
        % variables.
        PLV_Rest_I = PLV_Rest_I( 1:trials_compared ) ;

        PLV_Rest = zeros(trials_compared, size(PLV_cut,2) );
        PLV_Move= zeros(trials_compared, size(PLV_cut,2) );
        PLV_Diff = zeros(trials_compared, size(PLV_cut,2) );

        PLV_Diff_Sort = zeros(trials_compared, size(PLV_cut,2) );
        PLV_Diff_Coord = zeros(trials_compared, size(PLV_cut,2) );
        
        PLV_Mahal = zeros(trials_compared, size(PLV_cut,2) );
        Median_Mahal = zeros(size(PLV_cut,2), 1 );
        Avg_Mahal = zeros(size(PLV_cut,2), 1 );

%         PLV_Mahal_Sort = zeros(trials_compared, size(PLV_cut,2) );
%         PLV_Mahal_Coord = zeros(trials_compared, size(PLV_cut,2) );

        for chan_chan = 1:size(PLV_cut,2)

                PLV_Rest(:, chan_chan) = PLV_cut([PLV_Rest_I], chan_chan);
                PLV_Move(:, chan_chan) = PLV_cut([PLV_Move_I], chan_chan);
                PLV_Diff(:, chan_chan) = abs( PLV_Move(:, chan_chan)  -  PLV_Rest(:, chan_chan)) ;
                
                % PLV_Mahal - 
                PLV_Mahal(:, chan_chan) = mahal(PLV_Rest(:, chan_chan) ,  PLV_Move(:, chan_chan) );
                Median_Mahal(chan_chan) = median(PLV_Mahal(:, chan_chan));
                Avg_Mahal(chan_chan) = mean(PLV_Mahal(:, chan_chan));                
%                 PLV_Mahal_Sort{chan_chan}
        end

                %PLV_Mahal_Sort saves the 10 largest Mahal distances between equivalent Rest and Move trials 
                [Median_Mahal_Sort Median_Mahal_Ind] = sort(Median_Mahal, 'desc');
                [Avg_Mahal_Sort Avg_Mahal_Ind] = sort(Avg_Mahal, 'desc');
                
                Top_Median_Mahal_Sort = Median_Mahal_Sort(1:top_features)
                Top_Median_Mahal_Ind = Median_Mahal_Ind(1:top_features)
                
                % Mahal Coord are low to high
               Top_Median_Mahal_Chan = [col(Top_Median_Mahal_Ind), row(Top_Median_Mahal_Ind)];
               
               PLV_features = zeros(trials, top_features);
               
               for T = 1:trials
                   PLV_features(T,:) = PLV{T}(tril_ind(Top_Median_Mahal_Ind));
               end

    end

