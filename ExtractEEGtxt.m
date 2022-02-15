
q = 8;

% Editing the EEG Values
B = importdata(['D:\Kiel\Kiel_MEG_Data\Christmas_2019_Assignment\VEP_Measurement\VEP_Measurement\data\Christin_EEG_Values' num2str(q) '.txt'], ',')
C = cell(length(B.textdata)+1,1);
for i = 1:length(B.textdata)
    C{i} = erase(B.textdata(i), '/');
end
C{length(B.textdata)+1} = B.data;

comExpr = ',';
for la = 1:length(C)-1
    str = C{la};
    
    comStartIndex = regexp(str,comExpr);
    
    spaceVec=comStartIndex{1}(2:2:end);
    perVec = comStartIndex{1}(1:2:end);
    
    C{la}{1}([perVec]) = '.'; 
    C{la}{1}([spaceVec]) = ' ';
end
D = cell(length(C{length(B.textdata)+1}),1);
for o = 1:length(C{length(B.textdata)+1})
  % C{length(B.textdata)+1}(o) = sprintf('%d', C{length(B.textdata)+1}(o));
   D{o} = num2str(C{length(B.textdata)+1}(o));
end

for oz = 1:length(C{length(B.textdata)+1})/2
    C{length(B.textdata)+1}(oz) = str2num(strcat(D{oz*2-1},'.',D{oz*2}))
end
C{length(B.textdata)+1}(length(C{length(B.textdata)+1})/2+1:end)=[];

for j = 1:length(C)-1
    C{j} = str2num(C{j}{1});
end

p = 0;
for k = 1:length(C)
    
    p = p + length(C{k});
end

ReegMat8 = zeros(p, 1);
ReegMat8 = cat(2, C{:})'

save(['D:\Kiel\Kiel_MEG_Data\Dec-18_EEG_Measurements\Christin_EEG_Vector' num2str(q) '.mat'], ['CeegMat' num2str(q)])

%%
for i = 1
    load(['D:\Kiel\Kiel_MEG_Data\Dec-18_EEG_Measurements\Christin_EEG_Vector' num2str(1) '.mat']);
end

RyanEEGvaluesAll = [ReegMat1 ReegMat2 ReegMat3 ReegMat4 ReegMat5 ReegMat6 ReegMat7 ReegMat8];
ReegMatSumCheck = ReegMat1+ReegMat2+ReegMat3+ReegMat4;
ReegMatAvgCheck = ReegMatSumCheck/4;
ReegMatSumLED = ReegMat5+ReegMat6+ReegMat7+ReegMat8;
ReegMatAvgLED = ReegMatSumLED/4
save(['D:\Kiel\Kiel_MEG_Data\Dec-18_EEG_Measurements\RyanEEGvaluesAll.mat'], 'RyanEEGvaluesAll', 'ReegMatSumCheck', 'ReegMatAvgCheck', 'ReegMatSumLED', 'ReegMatSumLED')



ChristinEEGvaluesAll = [CeegMat1 CeegMat2 CeegMat3 CeegMat4 CeegMat5 CeegMat6 CeegMat7 CeegMat8];
CeegMatSumCheck = CeegMat1+CeegMat2+CeegMat3+CeegMat4;
CeegMatAvgCheck = CeegMatSumCheck/4;
CeegMatSumLED = CeegMat5+CeegMat6+CeegMat7+CeegMat8;
CeegMatAvgLED = CeegMatSumLED/4;
save(['D:\Kiel\Kiel_MEG_Data\Dec-18_EEG_Measurements\ChristinEEGvaluesAll.mat'], 'ChristinEEGvaluesAll', 'CeegMatSumCheck', 'CeegMatAvgCheck', 'CeegMatSumLED', 'CeegMatSumLED')

