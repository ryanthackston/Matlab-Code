
% Input Parameters - use defaults for now - Motor Img Cal, Subj 4
definput = cell(10,1);
definput(:) = [{'1'};{'4'};{'8.5'};{'12'};{'0.5'};{'3'};{'0.5'};{'0'};{'2'};{'50'}];
answers = meg_answer(definput);

srate = 1000;

 % stud{1} -> MIcal; stud{2} -> MIcont
[blocks, stud, subj] = meg_blocks_neu(answers);

% meg - channels x data
% Extract labels from MEG data for movement task
[labels, meg] = meg_labels_neu(answers, blocks, stud, subj);

% Find Onset and Offset Data of Tasks
[RestOnset, MoveOnset, RestOffset, MoveOffset] = meg_offsetdata(blocks);

% % REARRANGE THE CHANNEL NUMBERS TO GO FROM TOP-LEFT TO BOTTOM-RIGHT ON
% % SCALP MAP
% [meg] = meg_chanarrange(meg);

tw = 500;
shift = 100;
fs = 1000;
nfft = 1028;
freq = 0:fs/nfft:fs/2;
chans = size(meg,1);
trials = int64((size(meg, 2) - tw)/shift + 1);

% Initialize Variables
[meg_tw, Fmeg_tw, meg_f, meg_fft_singleTW, x_abs, x_pow, x_pow_str, feat, alpha, beta, w, freq_val, freq_width, PLV, d, B, I, d_sort, d_s, PLV_Rest_I, PLV_Move_I] =meg_initvars(answers, meg, fs,nfft, freq, chans, trials);

[meg_tw, Fmeg, PLV, PLV_Rest_I, PLV_Move_I, PLV_cut, row, col, tril_ind, tril_I] = meg_PLV3(meg, w, labels, trials, shift, tw, srate, freq_val, freq_width, PLV_Rest_I, PLV_Move_I);

top_features = 10;

[PLV_features, PLV_Mahal, Median_Mahal, Median_Mahal_Sort, Median_Mahal_Ind, Top_Median_Mahal_Sort, Top_Median_Mahal_Ind, PLV_Diff_Coord, PLV_Diff_Sort, trials_compared, PLV_Move_I, PLV_Rest_I] = meg_PLVfeatures2(top_features, trials, PLV, PLV_Move_I, PLV_Rest_I, PLV_cut, row, col, tril_ind);

[Class_Vec] = meg_featuresSorted2(PLV, trials);

[err] = meg_PLVclassify2(Class_Vec, PLV_features)
 
 
%  % Count frequency of channels being in top 10
%  z = [1: size(meg) ]; 
% count1 = countmember(z, Mahal_Vec(:,2));
% count2 = countmember(z, Mahal_Vec(:,3));
% chan_count = count1 + count2



