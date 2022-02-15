function [meg_tw, Fmeg, PLV, PLV_Rest_I, PLV_Move_I, PLV_cut, row, col] = meg_PLV(meg, w, labels, trials, shift, tw, srate, freq_val, freq_width, PLV_Rest_I, PLV_Move_I)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    for i = 1:( trials )
           % If i == 1, save meg_tw{i,1} as the 1st 500 time points
           % else shift the 500 point time window every 100 points 
           % meg_tw{i,1} --> channels x time_window
               if i == 1
                   % store meg data in meg_tw, needs to iterate with i
                   
                    meg_tw{i,1} = meg(:, (1:tw) )' .* w;
                    meg_tw{i,2} = labels(tw);
               else
                    meg_tw{i,1} = meg( :, [(1+(shift*(i-1) )) : (tw+(shift*(i-1) ))], :)' .* w;
                    meg_tw{i,2} = labels(tw+(shift*(i-1) ) );
               end
               
               meg_tw{i,1} = meg_tw{i,1}';

               % Check labels for Rest or Move
                if meg_tw{i,2} == 1
                     meg_tw{i,2} = 'Rest';
                elseif meg_tw{i,2} == 2
                    meg_tw{i,2} = 'Move';
                end 
                
                % Angles from Hilbert Transform to measure phase difference
                Fmeg{i,1} = filterFGx(meg_tw{i,1}, srate, freq_val, freq_width);
                % filtdata - chans x data
                for j = 1: size(meg_tw{i,1}, 2)
                    for k = 1: size(meg_tw{i,1},2)
                        % angle is in radians - chans x data
                        % Check Hilbert Angle
                        angts(j,k) = angle(hilbert(Fmeg{i,1}(j,k)').');
                    end
                end
                
                % Phase Locking Value
                %Create variables
                PLV{i,1} = zeros(size(meg_tw{i,1},1), size(meg_tw{i,1},1));
%                 d{i,1} = zeros(size(meg_tw{i,1},1), size(meg_tw{i,1},1));
                d{i,1} = zeros(size(meg_tw{i,1},1), 1);
                B{i,1} = zeros(size(meg_tw{i,1},1), 1);
                I{i,1} = zeros(size(meg_tw{i,1},1), 1);
                d_sort{i,1} = zeros(size(10, 1));


                % chans
                for chani = 1: size(meg_tw{i,1},2)
                    % chans
                    for chanj = 1: size(meg_tw{i,1},2)
                        tmpAi = angts(chani, :);
                        tmpAj = angts(chanj, :);
                        % This is the Phase Lag Value formula
                        % PLV - chans x data
                        PLV{i,1}(chani, chanj) = mean(abs(mean(exp(1i*(tmpAi-tmpAj )), 2)),1);
                        PLV{i,2} = meg_tw{i,2};
                            
                    end
                end
                
                if PLV{i,2} == "Rest"
                    PLV_Rest_I(i) = i;
                elseif PLV{i,2} == "Move"
                    PLV_Move_I(i) = i;
                end

    end
  
    % Cleanup non zero indices
    PLV_Rest_I = find(PLV_Rest_I ~= 0);
    PLV_Move_I = find(PLV_Move_I ~= 0);

    
% PLV gives mirror values, use tril to get only the phase-locked value under the diagonals.

tril_I = tril(PLV{1},-1) ~= 0;
tril_ind = find(tril_I == 1);
[row,col] = find(tril_I);

%The best feature will have 2 clusters that are the farthest apart from each other.
% To do this, we look at each feature, and within each feature we take all the samples 
% from one class and compairing it with all the samples in the second class.
% After this we do this with all the features, you have a "score" for each feature 
% which is the median mahal distances.
 
PLV_cut = cell(trials,1);
for i = 1:trials
    PLV_cut{i} = PLV{i}(tril_I);
end

end

