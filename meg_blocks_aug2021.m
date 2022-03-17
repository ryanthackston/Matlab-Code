function [blocks,stud, subj] = meg_blocks_aug2021(answers)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
    
%     stud_type = {micon_imag_aug2021(); mrcon_real_aug2021() };  
    stud_type = {mical_imag_aug2021(); mrcal_real_aug2021() };  
    stud = stud_type{ str2double(answers{1} )   };
    subj = str2double(answers{2});
    
%   D = inputdlg('Which day of sessions do you want to study from the subject (number 1-5)');
  D = str2num(answers{11});
  
        if D > 0
            for Sess = 1: size(stud{subj-11}{D}{1},1)

                if (Sess==1)
                    blocks = load([stud{subj-11}{D}{1}(Sess,:)]);
                else
                    blocks = [blocks; load([stud{subj-11}{D}{1}(Sess,:)]) ];
                end

            end
    

        % if D==0, concatenate all blocks  of all days of the same subject
        elseif D == 0
            z = [];
            for D = 1: size(stud{subj-11}{1}{1},1)
                tmp = stud{subj-11}{D}{1};
                z = char(z,tmp);
            end
            z(1,:) = [];

            for Sess = 1:size(z,1)
                if Sess == 1
                    blocks = load([stud{subj-11}{D}{1}(S,:)]);
                else
                    blocks = [blocks; load([stud{subj-11}{D}{1}(S,:)]) ];
                end
            end
        else
            error('D cannot be negative')
        end
    
end