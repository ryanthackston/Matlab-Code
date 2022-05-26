function  [RestOnset, MoveOnset, labels, meg, time_start_labels, time_start_meg] = meg_offsetdata_calib_aug2021(blocks, labels, answers,srate,meg)
%Will give Rest Onset, Move Onset, Rest Offset, and Move Offset points
% Concatenate RestOnset and MoveOnset times for the 2 blocks

    % Subtract 3100 points from No Action
        MoveOnset = [];
        RestOnset = [];
        MoveOffset = [];
        RestOffset = [];
    
        for b = 1:size(blocks,1)
            if b==1 
                % Subtract 3100 points for the beginning of every block and store the size of that into datsize
                datsize = size(blocks(b).meg,1)-3000;
                % Subtract 3100 points from every block of MoveOnset.
                MoveOnset = blocks(b).MoveOnset-3000;
                % Subtract 3100 points from every block of RestOnset.
                RestOnset = blocks(b).RestOnset-3000;
%                 % Subtract 3100 points from every block of MoveOffset.
%                 MoveOffset = blocks(b).MoveOffset-3100;
%                 % Subtract 3100 points from every block of RestOffset.
%                 RestOffset = blocks(b).RestOffset-3100;
            else
                MoveOnset = [MoveOnset; (blocks(b).MoveOnset-3000+datsize)];
                RestOnset = [RestOnset; (blocks(b).RestOnset-3000+datsize)];
%                 MoveOffset = [MoveOffset; (blocks(b).MoveOffset-3100+datsize)];
%                 RestOffset = [RestOffset; (blocks(b).RestOffset-3100+datsize)];
                datsize = datsize + size(blocks(b).meg,1)-3000;
            end
        end
        
        
        labels(RestOnset)=10;
%         labels(RestOffset)=22;
        labels(MoveOnset)=20;
%         labels(MoveOffset)=11;

        
% Off by 4000 points, 
% RestOnset(2) at 20001, in labels RestOnset(2) starts at 16001. 
% MoveOnset(1) at 12001, in labels MoveOnset(1) starts at 8001.

        for R = 1:length(RestOnset)
            
            if labels(RestOnset(R)) ~=0 && RestOnset(R) > 10 && labels(RestOnset(R)-1) == 0
                labels(RestOnset(R)-4000:RestOnset(R) -1) = [];
                 MoveOnset(find(MoveOnset > RestOnset(R) )) = MoveOnset(find(MoveOnset > RestOnset(R) )) - 4000;
                if R ~= length(RestOnset)
                    RestOnset( (R+1) : length(RestOnset)) = RestOnset( (R+1):length(RestOnset) ) - 4000;
                else
                end

            else
                continue;
            end
%              labels(RestOnset(R)+(str2num(answers{8})*srate):RestOnset(R) +7999)
%             meg(RestOnset(R)+(str2num(answers{8})*srate):RestOnset(R) +8000)
        end
        
        
        for M = 1:length(MoveOnset)
            if MoveOnset(M) ~= 0 && MoveOnset(M) > 10 && labels(MoveOnset(M)-1) == 0
                labels(MoveOnset(M) - 4000: MoveOnset(M) -1) = [];
                RestOnset(find(RestOnset > MoveOnset(M) )) = RestOnset(find(RestOnset > MoveOnset(M) )) - 4000;
                if M ~= length(MoveOnset)
                    MoveOnset( (M+1) : length(MoveOnset)) = (MoveOnset( (M+1) : length(MoveOnset))) - 4000;
                else
                end
    
            else
            end
%              labels(MoveOnset(M)+(str2num(answers{8})*srate):MoveOnset(R) -1)
%             meg(MoveOnset(M)+(str2num(answers{8})*srate):MoveOnset(R) -1)
        end
        
       labels( find(labels==0) ) = [];
       MoveOnset = MoveOnset-4000;
       
       TrialsOnset = sort([RestOnset; MoveOnset]);
       
       time_start_meg = [];
       time_start_labels = [];
       
       len_trial = size(TrialsOnset(1):TrialsOnset(2)-1,2);
       
       for T = 1:size(TrialsOnset)
           time_start_meg = [ time_start_meg; meg( [ (TrialsOnset(T) + (str2num(answers{8})*srate)):(TrialsOnset(T) + len_trial-1) ],:) ];
           time_start_labels = [ time_start_labels; labels( [ (TrialsOnset(T) + (str2num(answers{8})*srate)):(TrialsOnset(T) + len_trial-1) ]) ];
       end
       
       meg = time_start_meg;
       labels = time_start_labels;

end
            
%     %CHECK FOR OFFSET DATA. 
%     % IF THERE'S OFFSET DATA  MAKE THE SAME SIZE OF OFFSET TIMES AS ONSET TIMES
%             try
%                 blocks.MoveOffset;
% 
%                 for i = 1:2
%                     % Make Equal Onset, Offset, & trialEnd indexes
%                     
%                     if length(blocks(i).MoveOnset) ~= length(blocks(i).MoveOffset)
%                         blocks(i).trialEnd(length(blocks(i).trialEnd) + 1) = length(blocks(i).meg);                        
%                         blocks(i).MoveOffset(length(blocks(i).MoveOffset) + 1) = length(blocks(i).meg);
%                     elseif length(blocks(i).RestOnset) ~= length(blocks(i).RestOffset)
%                           blocks(i).trialEnd(length(blocks(i).trialEnd) + 1) = length(blocks(i).meg);   
%                           blocks(i).RestOffset(length(blocks(i).RestOffset) + 1) = length(blocks(i).meg);
%                     end
%                 end
% 
%                 % Combine RestOffset, MoveOffset, trialStart, trialEnd,
%                 % trialLength index time points for 
%                 RestOffset = [blocks(1).RestOffset(:) - 3000; blocks(2).RestOffset(:)+length(blocks(1).meg(3001:end, :)) - 3000];
%                 MoveOffset = [blocks(1).MoveOffset(:) - 3000; blocks(2).MoveOffset(:)+length(blocks(1).meg(3001:end, :)) - 3000];
%                 trialStart = [blocks(1).trialStart(:) - 3000; blocks(2).trialStart(:)+length(blocks(1).meg(3001:end, :)) - 3000];
%                 trialEnd = [blocks(1).trialEnd(:) - 3000; blocks(2).trialEnd(:)+length(blocks(1).meg(3001:end, :)) - 3000];       
%                 trialLength = cell(size(trialStart,1),2);
% 
%                 for i = 1:size(trialStart,1)
%                         trialLength{i,1} = trialEnd(i) - trialStart(i);  
% 
%                         if any(trialStart(i) == MoveOnset)
%                            trialLength{i, 2}= 'MoveOnset';
%                         elseif any(trialStart(i) == RestOnset)
%                             trialLength{i, 2} = 'RestOnset';
%                         end
%                 end
%                 % For visuals
%     %             trials = cat(2, num2cell(trialStart), num2cell(trialEnd), trialLength(:, 1), trialLength(:, 2));
%            % If there is no offset data, then catch the error and display message
%             catch
%                 disp('There is no MoveOffset or RestOffset in these trial blocks');
%                 RestOffset = [];
%                 MoveOffset = [];
%              end
% end

