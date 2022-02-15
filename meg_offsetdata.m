function [RestOnset, MoveOnset, RestOffset, MoveOffset] = meg_offsetdata(blocks)
%Will give Rest Onset, Move Onset, Rest Offset, and Move Offset points
% Concatenate RestOnset and MoveOnset times for the 2 blocks

% Subtract 3000 points from No Action
    RestOnset = [blocks(1).RestOnset(:) - 3000; blocks(2).RestOnset(:)+length(blocks(1).meg(3001:end, :)) - 3000];
    MoveOnset = [blocks(1).MoveOnset(:) - 3000; blocks(2).MoveOnset(:)+length(blocks(1).meg(3001:end, :)) - 3000];
    %CHECK FOR OFFSET DATA. 
    % IF THERE'S OFFSET DATA  MAKE THE SAME SIZE OF OFFSET TIMES AS ONSET TIMES
            try
                blocks.MoveOffset;

                for i = 1:2
                    % Make Equal Onset, Offset, & trialEnd indexes
                    
                    if length(blocks(i).MoveOnset) ~= length(blocks(i).MoveOffset)
                        blocks(i).trialEnd(length(blocks(i).trialEnd) + 1) = length(blocks(i).meg);                        
                        blocks(i).MoveOffset(length(blocks(i).MoveOffset) + 1) = length(blocks(i).meg);
                    elseif length(blocks(i).RestOnset) ~= length(blocks(i).RestOffset)
                          blocks(i).trialEnd(length(blocks(i).trialEnd) + 1) = length(blocks(i).meg);   
                          blocks(i).RestOffset(length(blocks(i).RestOffset) + 1) = length(blocks(i).meg);
                    end
                end

                % Combine RestOffset, MoveOffset, trialStart, trialEnd,
                % trialLength index time points for 
                RestOffset = [blocks(1).RestOffset(:) - 3000; blocks(2).RestOffset(:)+length(blocks(1).meg(3001:end, :)) - 3000];
                MoveOffset = [blocks(1).MoveOffset(:) - 3000; blocks(2).MoveOffset(:)+length(blocks(1).meg(3001:end, :)) - 3000];
                trialStart = [blocks(1).trialStart(:) - 3000; blocks(2).trialStart(:)+length(blocks(1).meg(3001:end, :)) - 3000];
                trialEnd = [blocks(1).trialEnd(:) - 3000; blocks(2).trialEnd(:)+length(blocks(1).meg(3001:end, :)) - 3000];       
                trialLength = cell(size(trialStart,1),2);

                for i = 1:size(trialStart,1)
                        trialLength{i,1} = trialEnd(i) - trialStart(i);  

                        if any(trialStart(i) == MoveOnset)
                           trialLength{i, 2}= 'MoveOnset';
                        elseif any(trialStart(i) == RestOnset)
                            trialLength{i, 2} = 'RestOnset';
                        end
                end
                % For visuals
    %             trials = cat(2, num2cell(trialStart), num2cell(trialEnd), trialLength(:, 1), trialLength(:, 2));
           % If there is no offset data, then catch the error and display message
            catch
                disp('There is no MoveOffset or RestOffset in these trial blocks');
                RestOffset = [];
                MoveOffset = [];
             end
end

