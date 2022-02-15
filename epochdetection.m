

function [trigger, trigger_change, epoch_times, epoch_opm, standard_opm, epoch_trigger, standard_trigger] = epochdetection(xnow, thresh)
% Ex: [x.epoch.trial_start, x.trigger, x.cut.trigger x.cut.opm] = epochdetection(x, thresh)
%
% The input threshold for an ON/OFF stimulus will be the mean signal value
% thresh = mean(x.in(:,5));
%
% trigger - Creates a column vector of the trigger signals. If the trigger
% is 1 at that point, the stimulus is ON. Likewise, a value of 0 means the
% stimulus is OFF.
% 
% trigger_change - Creates a column vector of time points right before a
% trigger value change. 1 means the trigger value will change at the next
% point. 0 means no change in trigger value will occur
%
% epoch_times - Creates a double vector of time points when the trigger
% changed
%
% epoch_opm - Creates as many cells as there are epoch_times values. Stores
% measured signals from all channels. Each cell will store measured signals
% from the same indexed epoch_times to the time point right before the next
% epoch_times value.



% Create a column vector of zeros for the trigger signals
trigger = zeros(length(xnow.x.in(:,5)),1);

% Create a column vector of zeros for the trigger_change points
trigger_change = zeros(length(trigger),1);



% If the time point is below the threshold, keep the value as 0 and 
% restart the for loop. Otherwise, save save the time point as 1 and 
% restart the for loop.
for c = 1:length(xnow.x.in(:,5))
    
    if xnow.x.in(c,5) > thresh
        trigger(c) = 1;
    else
        trigger(c) = 0;
    end
    
    if c == 1
        trigger_change(c) = 0;
        continue;
    elseif trigger(c-1) == trigger(c)
        trigger_change(c-1) = 0;
    elseif trigger(c-1) ~= trigger(c)
        trigger_change(c-1) = 1;
    elseif c == length(xnow.x.in(:,5))
        trigger_change(c) = 0;
    end
end

% Find the time points right before the triggers change values and store in a
% double column vector
change_times = find(trigger_change == 1);

% Save where triggers change values and store every 2nd change
% in a double column vector. Each of these values will
% be the desired epoch where the stimuli will have 1 complete cycle of ON
% and OFF.
epoch_times = change_times(1:2:end);

% In each cell, store all OPM signals until the time point before the next
% epoch_times value. This will be the epoch.
epoch_opm = cell(length(epoch_times),1);
epoch_trigger = cell(length(epoch_times),1);
standard_opm = cell(length(epoch_times),1);
standard_trigger = cell(length(epoch_times),1);

for e = 1:length(epoch_times)
    % If you are storing the last epoch_times value, save the rest of the
    % OPM signals to the end
    if e == length(epoch_times)
        epoch_opm{e} = xnow.x.b(epoch_times(e)+1:1:end,:);
        epoch_trigger{e} = trigger(epoch_times(e)+1:1:end,:);
        standard_opm{e} = zscore(epoch_opm{e});
        standard_trigger{e} = zscore(epoch_trigger{e});
    else
        epoch_opm{e} = xnow.x.b(epoch_times(e)+1:1:epoch_times(e+1), :);
        epoch_trigger{e} = trigger(epoch_times(e)+1:1:epoch_times(e+1), :);
        standard_opm{e} = zscore(epoch_opm{e}); 
        standard_trigger{e} = zscore(epoch_trigger{e});
    end
    
    %Create average epoch_opm, make # of epochs to average and range
end
end





















