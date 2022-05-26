function [meg_tw, Fmeg, PLV, PLV_Rest_I, PLV_Move_I, PLV_cut, row, col, tril_ind, tril_I] = meg_PLV3(meg, w, labels, trials, shift, tw, srate, freq_val, freq_width, PLV_Rest_I, PLV_Move_I)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    for T = 1:( trials )
           % If T == 1, save meg_tw{T,1} as the 1st 500 time points
           % else shift the 500 point time window every 100 points 
           % meg_tw{T,1} --> channels x time_window
               if T == 1
                   % store meg data in meg_tw, needs to iterate with i
                   
                    meg_tw{T,1} = meg(:, (1:tw) )' .* w;
                    meg_tw{T,2} = labels(tw);
               else
                    meg_tw{T,1} = meg( :, [(1+(shift*(T-1) )) : (tw+(shift*(T-1) ))], :)' .* w;
                    meg_tw{T,2} = labels(tw+(shift*(T-1) ) );
               end
               
               meg_tw{T,1} = meg_tw{T,1}';

               % Check labels for Rest or Move
                if meg_tw{T,2} == 1
                     meg_tw{T,2} = 'Rest';
                elseif meg_tw{T,2} == 2
                    meg_tw{T,2} = 'Move';
                end 
                
                % Angles from Hilbert Transform to measure phase difference
                % Fmeg is filtered MEG data
                % -10000 values are messing up the Fmeg
                Fmeg{T,1} = filterFGx(meg_tw{T,1}, srate, freq_val, freq_width);
                % filtdata - chans x data
%                 for j = 1: size(meg_tw{T,1}, 1)
%                     for k = 1: size(meg_tw{T,1},2)
%                         % angle is in radians - chans x data
%                         % Check Hilbert Angle
%                         h(j,k) = hilbert( (Fmeg{T,1}(j,k))' );
%                         angts(j,k) = angle(hilbert( (Fmeg{T,1}(j,k))' ) );
% %                         angts(j,k) = angle(hilbert(Fmeg{T,1}(j,k)').');
%                     end
%                 end
                angts = angle(hilbert(Fmeg{T,1}'))';
                
                % Phase Locking Value
                %Create variables
                PLV{T,1} = zeros(size(meg_tw{T,1},1), size(meg_tw{T,1},1));
%                 d{T,1} = zeros(size(meg_tw{T,1},1), size(meg_tw{T,1},1));
                d{T,1} = zeros(size(meg_tw{T,1},1), 1);
                B{T,1} = zeros(size(meg_tw{T,1},1), 1);
                I{T,1} = zeros(size(meg_tw{T,1},1), 1);
                d_sort{T,1} = zeros(size(10, 1));


                % chans 
%                 for chani = 1: size(meg_tw{T,1},2)
                    for chani = 1: size(angts,1)
                    % chans
%                     for chanj = 1: size(meg_tw{T,1},2)
                        for chanj = 1: size(angts,1)
                        tmpAi = angts(chani, :);
                        tmpAj = angts(chanj, :);
                        % This is the PLV formula
                        % PLV - chans x data
                        % 4996 trials, 17x17 PLV
                        PLV{T,1}(chani, chanj) = mean(abs(mean(exp(1i*(tmpAi-tmpAj )), 2)),1);
                        PLV{T,2} = meg_tw{T,2};
                            
                    end
                end
                
                if PLV{T,2} == "Rest"
                    PLV_Rest_I(T,1) = T;
                elseif PLV{T,2} == "Move"
                    PLV_Move_I(T,1) = T;
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

PLV_cut = zeros(trials, length(tril_ind) );
for T = 1:trials
    % Each row is a trial and each column is a channel interaction
    PLV_cut(T, :) = PLV{T}(tril_I);
end

end

