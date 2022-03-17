function [labels, meg, class_label, temp_ind] = meg_labels_calib_aug2021(answers, blocks, stud, subj)

% take out initial 3000 points due to No Action

meg = [];
for b = 1:size(blocks,1)
    meg = [meg; blocks(b).meg(3001:end, :)];
end

        % else task is motor img control
        % Take out 1st 3100 points
        % in between trials 999 0s
        % I want the signals that are 1 or 2 but I want to keep the index in
        % line with the meg data
        labels = [];
        for L = 1:size(blocks,1)
            labels = [ labels; blocks(L).stimulusCode(3001:end, :) ];
        end
        % Initialize class_label
        class_label = cell(length(labels),1);
        % create temporary indices for all labels that are 1 and 2
        temp_ind = find(labels>0 & labels<3);
        for jj = 1:length(temp_ind)
            if labels(temp_ind(jj)) == 1
                class_label{temp_ind(jj)} = 'Move';
            elseif labels(temp_ind(jj)) == 2
                class_label{temp_ind(jj)} = 'Rest';
            end
        end
        % All I care about is meg & label data during the task trials (1 & 2)
        
%         labels = labels(temp_ind);
        meg = meg(temp_ind,:);
        
        
        z = [];
        z = [z class_label(temp_ind)];
        class_label = z;
        
        
end

