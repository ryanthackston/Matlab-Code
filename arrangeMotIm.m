function [data, megR3d, megM3d] = arrangeMotIm(blocks)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
RestOnset = [blocks(1).RestOnset(:); blocks(2).RestOnset(:)+253000];
MoveOnset = [blocks(1).MoveOnset(:); blocks(2).MoveOnset(:)+253000];

% Delete unused channels
blocks(1).meg(:,[4 6 7 11 15 16 20 23]) = [];
blocks(2).meg(:,[4 6 7 11 15 16 20 23]) = [];

meg = [blocks(1).meg; blocks(2).meg];
meg = meg';

% Rearrange the channel numbers to go sequentially from top-left to bottom-right
megMat = [];
megMat(1,:) = meg(12,:);
megMat(2,:) = meg(7,:);
megMat(3,:) = meg(6,:);
megMat(4,:) = meg(17,:);
megMat(5,:) = meg(14,:);
megMat(6,:) = meg(8,:);
megMat(7,:) = meg(11,:);
megMat(8,:) = meg(16,:);
megMat(9,:) = meg(5,:);
megMat(10,:) = meg(10,:);
megMat(11,:) = meg(4,:);
megMat(12,:) = meg(15,:);
megMat(13,:) = meg(9,:);
megMat(14,:) = meg(3,:);
megMat(15,:) = meg(1,:);
megMat(16,:) = meg(13,:);
megMat(17,:) = meg(2,:);
meg = megMat(:,:);

% Rearrange data to be channels x time points x trials 
% data(:, 1:5000, :) time points is rest task. 
% data(:, 5001:10000, :) time points is move task.
megR3d = zeros(17,5000,size(RestOnset,1));
megM3d = zeros(17,5000,size(MoveOnset,1));
data = zeros(17,10000,size(RestOnset,1));

% motor img - calib
for i = 1:size(RestOnset,1)
    megR3d(:,:,i) = meg(:,(RestOnset(i):MoveOnset(i)-1));
    megM3d(:,:,i) = meg(:,(MoveOnset(i):MoveOnset(i)+4999));
    data(:,:,i) = [megR3d(:,:,i) megM3d(:,:,i)];
end

end

