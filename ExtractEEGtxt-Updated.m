
STIM = [1.7; 3.4; 7];

StimFreq = {'F1o7', 'F3o4', 'F7'};

Trials = {'T1', 'T2', 'T3'};
% 
EEGmatAS = {'AS_1o7Hz_Trial1','AS_1o7Hz_Trial2', 'AS_1o7Hz_Trial3';...
            'AS_3o4Hz_Trial1','AS_3o4Hz_Trial2', 'AS_3o4Hz_Trial3';...
            'AS_7Hz_Trial1','AS_7Hz_Trial2','AS_7Hz_Trial3'};

        
for jj = 1:length(STIM)
for t = 1:3
% Import the EEG values as a string
B = importdata(['C:\Users\Ryan\Downloads\MATLAB Functions\AS_Trials\AS_' num2str(STIM(jj)) 'Hz_Trial' num2str(t) '.txt'], ',')

% Add cells for every line including the data line
C = cell(length(B),1);

for i = 1:length(B)
    C{i} = erase(B(i), '/');
end
% C{length(B.textdata)} = B.data;

comExpr = ',';
% Replace the comma with period
for lenc = 1:length(C)
    str = C{lenc};
    
    comStartIndex = regexp(str,comExpr);
    
%     spaceVec=comStartIndex{1}(2:2:end);
    perVec = comStartIndex{1}(1:end);
    
    C{lenc}{1}([perVec]) = '.'; 
%     C{lenc}{1}([spaceVec]) = ' ';
end

D = cell(length(C),1);
for n = 1:length(C)
    D{n} = str2num(C{n}{1});
end

p = [];
for z = 1:length(D)
    p = [p; D{z}];
end
%Scale p by 1000
p = p*1000;

% ASstruct.(EEGmatAS{jj,t}) = p;
AS.(StimFreq{jj}).(Trials{t}) = p;

end
end

save(['C:\Users\Ryan\Downloads\MATLAB Functions\AS_Trials\AS_EEGtrials_Scaled1000.mat'], 'AS')

%%
for i = 1
    load(['D:\Kiel\Kiel_MEG_Data\Dec-18_EEG_Measurements\Christin_EEG_Vector' num2str(1) '.mat']);
end

RyanEEGvaluesAll = [ReegMat1 ReegMat2 ReegMat3 ReegMat4 ReegMat5 ReegMat6 ReegMat7 ReegMat8];
ReegMatSumCheck = ReegMat1+ReegMat2+ReegMat3+ReegMat4;
ReegMatAvgCheck = ReegMatSumCheck/4;
ReegMatSumLED = ReegMat5+ReegMat6+ReegMat7+ReegMat8;
ReegMatAvgLED = ReegMatSumLED/4;
save(['D:\Kiel\Kiel_MEG_Data\Dec-18_EEG_Measurements\RyanEEGvaluesAll.mat'], 'RyanEEGvaluesAll', 'ReegMatSumCheck', 'ReegMatAvgCheck', 'ReegMatSumLED', 'ReegMatSumLED')



ChristinEEGvaluesAll = [CeegMat1 CeegMat2 CeegMat3 CeegMat4 CeegMat5 CeegMat6 CeegMat7 CeegMat8];
CeegMatSumCheck = CeegMat1+CeegMat2+CeegMat3+CeegMat4;
CeegMatAvgCheck = CeegMatSumCheck/4;
CeegMatSumLED = CeegMat5+CeegMat6+CeegMat7+CeegMat8;
CeegMatAvgLED = CeegMatSumLED/4;
save(['D:\Kiel\Kiel_MEG_Data\Dec-18_EEG_Measurements\ChristinEEGvaluesAll.mat'], 'ChristinEEGvaluesAll', 'CeegMatSumCheck', 'CeegMatAvgCheck', 'CeegMatSumLED', 'CeegMatSumLED')

