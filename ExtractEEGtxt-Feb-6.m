
% Current method for extracting EEG data from clinical equipment
% 1. Open a text file containing trials recorded using EEG
% 2. Copy and paste each individual trial set into its own text file. 
% 3. Add an '/' on the last line of each text file to make the data easier
%    to import.
% 4. Name each file individually, distinguishing that trial from others.
% 5. Use 'File' variable to distinguish the different trials
% 6. 'StimFreq' & 'Trials' are for saving the variable names in a structure
%    to make the data easier to access and more organized.



File = {'2-Frame 7.5 Hz Epoch 1'; '2-Frame 7.5 Hz Epoch 2'; '2-Frame 7.5 Hz Epoch 3';...
        '2-Frame 7.5 Hz Epoch 4'; '2-Frame 7.5 Hz Epoch 5'; '2-Frame 7.5 Hz Epoch 6';...
        '1-Frame 15 Hz Epoch 1'; '1-Frame 15 Hz Epoch 2'; '1-Frame 15 Hz Epoch 3';...
        '1-Frame 15 Hz Epoch 4'; '1-Frame 15 Hz Epoch 5'; '1-Frame 15 Hz Epoch 6';...
        '2-Frame 1.9 Hz Avg 1'; '2-Frame 1.9 Hz Avg 2'; ...
        '2-Frame 10 Hz Avg 1'; '2-Frame 10 Hz Avg 2'; ...
        '1-Frame 10 Hz Avg 1'; '1-Frame 10 Hz Avg 2'};
    
StimFreq = {'F7o5', 'F7o5', 'F7o5', 'F7o5', 'F7o5', 'F7o5',...
            'F15', 'F15', 'F15', 'F15', 'F15', 'F15', ...
            'F1o9', 'F1o9',...
            'F10', 'F10', 'F10', 'F10',};
        
Trials = {'T1_Frame2_Epoch', 'T2_Frame2_Epoch', 'T3_Frame2_Epoch', 'T4_Frame2_Epoch', 'T5_Frame2_Epoch', 'T6_Frame2_Epoch',...
          'T1_Frame1_Epoch', 'T2_Frame1_Epoch', 'T3_Frame1_Epoch', 'T4_Frame1_Epoch', 'T5_Frame1_Epoch', 'T6_Frame1_Epoch',...
          'T1_Frame2_Avg', 'T2_Frame2_Avg',...
          'T1_Frame2_Avg', 'T2_Frame2_Avg',...
          'T1_Frame1_Avg', 'T2_Frame1_Avg'};


      
for jj = 1:length(File)
% Import the EEG values as a string 
B = importdata(['C:\Users\Ryan\Downloads\MATLAB Functions\AS_Trials\Feb-6\Feb-6 AS ' File{jj} '.txt'], ',');

% Add cells for every line including the data line. There should be no 'data'
% structure if you add '/' on the last line of each text file.
C = cell(length(B),1);

% Erase each '/'
for i = 1:length(B)
    C{i} = erase(B(i), '/');
end

comExpr = ',';
% Replace the comma ',' with a period '.'
% Note: Comma isn't an ordinary comma, you need to copy and paste the
% symbol from the text file.
% For 1 to length of lines
for lenc = 1:length(C)
    str = C{lenc};
    % Each line is now its own string
    comStartIndex = regexp(str,comExpr);
    
%     spaceVec=comStartIndex{1}(2:2:end);
    % Creates a vector of all indexed values with ',' on that line
    perVec = comStartIndex{1}(1:end);
    
    % Replace all values of the index with a period
    C{lenc}{1}([perVec]) = '.'; 
%     C{lenc}{1}([spaceVec]) = ' ';
end

% Converts strings into values, with each line inside its own cell
D = cell(length(C),1);
for n = 1:length(C)
    % For 1 to length of lines, convert strings 2 numbers and add to cell
    % in D
    D{n} = str2num(C{n}{1});
end

% convert data in D to one large vector of values p
p = [];
for z = 1:length(D)
    p = [p; D{z}];
end

%Scale p by 1000 to represent microvolts
p = p*1000;

% ASstruct.(EEGmatAS{jj,t}) = p;
% Name your session variable whatever you wish. 
% Ex AS750 for subject AS with epoch 750ms
% Use StimFreq and Trials to structure your data.
AS750.(StimFreq{jj}).(Trials{jj}) = p;

end

% Save data 'AS750' in location 'C:\...' with name '...\Feb6_AS_EEGtrials_Scaled1000.mat'
save(['C:\Users\Ryan\Downloads\MATLAB Functions\AS_Trials\Feb6_AS_EEGtrials_Scaled1000.mat'], 'AS750')




















