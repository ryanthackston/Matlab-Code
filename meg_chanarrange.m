function [meg] = meg_chanarrange(meg)
    % REARRANGE THE CHANNEL NUMBERS TO GO FROM TOP-LEFT TO BOTTOM-RIGHT ON
    % SCALP MAP
    megarr = [12 7 6 17 14 8 11 16 5 10 4 15 9 3 1 13 2]';
    megmat = meg([megarr], :);
    meg = megmat;
end

