
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

folder = 'C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Calibration MAT files\';
chanfile = 'C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Aug2021_chanLocs.loc';

for j = 13

    for i = 4:5

        file1 = ['S', num2str(j),'_D', num2str(i), '_Block1_imag_calib'];
        file2 = ['S', num2str(j),'_D', num2str(i), '_Block2_imag_calib'];

        tmp1 = load([folder, file1]);
        tmp2 = load([folder, file2]);

        meg = [tmp1.meg; tmp2.meg]';

        save([folder, ['S', num2str(j),'_D', num2str(i),'_meg'] ], 'meg' );

        
        
    end
end
    
EEG.chanlocs=pop_chanedit(EEG.chanlocs, 'load',{ chanfile, 'filetype', 'autodetect'});


pop_loadset('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\S13-D2-imag-calib.set');

    
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

slidingPLV_Aug2021

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

pop_loadset('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\S13-D2-imag-calib.set');

topoGUI(krusP,SS_chan_inter, EEG, answers);