
% Input Parameters - use defaults for now - Motor Img Cal, Subj 4
definput = cell(11,1);
definput(:) = [{'1'};{'13'};{'14'};{'30'};{'0.5'};{'3'};{'0.5'};{'0'};{'2'};{'50'};{'2'}];
answers = meg_answer_aug2021(definput);

srate = 1000;

 % stud{1} -> MIcal; stud{2} -> MIcont
[blocks,stud, subj] = meg_blocks_aug2021(answers);

% meg - channels x data
% Extract labels from MEG data for movement task
% [labels, meg, class_label, temp_ind] = meg_labels_aug2021(answers, blocks, stud, subj);
[labels, meg, class_label, temp_ind] = meg_labels_calib_aug2021(answers, blocks, stud, subj);

% Find Onset and Offset Data of Tasks
% [RestOnset, MoveOnset, RestOffset, MoveOffset, labels] = meg_offsetdata_aug2021(blocks, labels);

 [RestOnset, MoveOnset, labels, meg, time_start_labels, time_start_meg] = meg_offsetdata_calib_aug2021(blocks, labels, answers,srate,meg);

% % REARRANGE THE CHANNEL NUMBERS TO GO FROM TOP-LEFT TO BOTTOM-RIGHT ON
% % SCALP MAP
% [meg] = meg_chanarrange(meg);

meg = meg';
tw = str2num(answers{5})*srate;
shift = 100;
nfft = 1028;
freq = 0:srate/nfft:srate/2;
chans = size(meg,1);
trials = int64((size(meg, 2) - tw)/shift + 1);

%Change onset labels back
labels(find(labels == 10)) = 1;
labels(find(labels == 20)) = 2;

% Initialize Variables
[meg_tw, Fmeg_tw, meg_f, meg_fft_singleTW, x_abs, x_pow, x_pow_str, feat, alpha, beta, w, freq_val, freq_width, PLV, d, B, I, d_sort, d_s, PLV_Rest_I, PLV_Move_I] =meg_initvars(answers, meg, srate,nfft, freq, chans, trials);

[meg_tw, Fmeg, PLV, PLV_Rest_I, PLV_Move_I, PLV_cut, row, col, tril_ind, tril_I] = meg_PLV3(meg, w, labels, trials, shift, tw, srate, freq_val, freq_width, PLV_Rest_I, PLV_Move_I);

top_features = 10;

[PLV_features, PLV_Mahal, Median_Mahal, Median_Mahal_Sort, Median_Mahal_Ind, Top_Median_Mahal_Sort, Top_Median_Mahal_Ind, PLV_Diff_Coord, PLV_Diff_Sort, trials_compared, PLV_Move_I, PLV_Rest_I] = meg_PLVfeatures2(top_features, trials, PLV, PLV_Move_I, PLV_Rest_I, PLV_cut, row, col, tril_ind);

[Class_Vec] = meg_featuresSorted2(PLV, trials);

folds = 5;

decoderScript_SpectralPower

[C err P logp coeff PredictedVsActualClasses testAcc] = meg_PLVclassify_with_Spectral_Power(Class_Vec, PLV_features, featureTesting_best, featureTraining_best, featureMatrixTrials)

[C err P logp coeff PredictedVsActualClasses testAcc]  = meg_PLVclassify2(Class_Vec, PLV_features)

Top_Channel_Pairs = [row([Top_Median_Mahal_Ind]) col([Top_Median_Mahal_Ind]) ];

[ testAcc, testAcc_rand] = meg_SVM(PLV_features, Class_Vec, folds, testIdx)

[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;