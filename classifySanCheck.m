% A = 10 + 2*randn(1000, 1);
% B = -7 + 2*randn(1000,1);
% 
% 
% [ testIdx, trainIdx, nTest] = cvpartition1(size(A,1),5);
% 
% CL = cell(size(A,1) + size(B,1),1);
% for i = 1:size(A,1)
%     CL{i} = 'A';
%     CL{ i + size(A,1) } = 'B';
% end
% 
% [ CLtestIdx, CLtrainIdx, CLTest] = cvpartition1(size(CL,1),5);
% 
% for i = 1:5
%     [C{i, 1},err(i, 1),P{i, 1},logp{i, 1},coeff{i, 1}] = classify( [ A(testIdx(i,:)) ; B(testIdx(i,:))], ...
%                                                                   [A(trainIdx(i,:)) ; B(trainIdx(i,:))] , ...
%                                                                   CL(find(CLtrainIdx(i,:))), 'linear' );
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = -7 + 2*randn(1000, 1);
B = 6 + 2*randn(1000,1);
data = [A;B];

[ testIdx, trainIdx, nTest] = cvpartition1(size(data,1),5);

CL = cell(size(A,1) + size(B,1),1);
for i = 1:size(A,1)
    CL{i} = 'A';
    CL{ i + size(A,1) } = 'B';
end

[ CLtestIdx, CLtrainIdx, CLTest] = cvpartition1(size(CL,1),5);

for i = 1:5
    [C{i, 1},err(i, 1),P{i, 1},logp{i, 1},coeff{i, 1}] = classify( data(testIdx(i,:)), ...
                                                                  data(trainIdx(i,:)) , ...
                                                                  CL(find(trainIdx(i,:))), 'linear' );
end