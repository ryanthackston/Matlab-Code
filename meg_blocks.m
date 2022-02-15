function [blocks,stud, subj] = meg_blocks(answers)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
stud_type = {mical(); micont() };  
stud = stud_type{str2double(answers{1})};
subj = str2double(answers{2});
blocks = [load(stud{subj}(1,:)); load(stud{subj}(2,:))]; 
end