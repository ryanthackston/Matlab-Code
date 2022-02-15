% Create variables from all sessions to use for analysis
% Organize the variables cleanly

Trials = {'trial_1','trial_2'}; 
nov12_opm_noise.(Trials{1}) = load('C:\Users\Ryan\Desktop\Kiel MEG Data\Nov-12 Noise Measurement with LED\dat_meg_measurement_165544');
nov12_opm_noise.(Trials{2}) = load('C:\Users\Ryan\Desktop\Kiel MEG Data\Nov-12 Noise Measurement with LED\dat_meg_measurement_170159');

save('C:\Users\Ryan\Desktop\Kiel MEG Data\Data\Nov12_opm_noise.mat', 'nov12_opm_noise');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nov 21 Data      
Trials = {'FiberOptics_7Hz_OPMoff_Subject1Inside_3min',...
          'FiberOptics_7Hz_Subject1_3min',...
          'WhiteLED_7Hz_Subject1_3min',...
          'GreenLED_7Hz_Subject1_3min',...
          'GreenLED_10Hz_Subject1_3min',...
          'GreenLED_10Hz_Noise_3min',...
          'Subject2_Noise_3min',...
          'FiberOptic_10Hz_Subject2_3min',...
          'FiberOptic_10Hz_Noise_3min'};
      
      
      
      
      Folder = 'C:\Users\Ryan\Desktop\Kiel MEG Data\Nov-21-2019_Initial_MEG_Measurement\';
      
      Files = {'dat_dat_Human_inside_fiber_optics_on_noise_measurment144252',...
               'dat_dat_Human_inside_fiber_optics_on_person_2_10Hz_3min143408',...
               'dat_dat_Human_inside_fiber_optics_off_person_2_10Hz_3min143008',...
               'dat_noise_measurement_10_Hz_green140838',...
               'dat_Human_inside_green_led_on_person_1_3min_10Hz140043',...
               'dat_Human_inside_green_led_on_person_1_7_Hz_3min135530',...
               'dat_Human_inside_white_led_on_person_1_7_Hz_3min134811',...
               'dat_Human_inside_fiber_optics_on_person_1_7Hz_3min133833',...
               'dat_Human_inside_fiber_optics_off_person_1_7_Hz_3min133118'} 
           
 for i = 1:length(Files)
     nov21_opm.(Trials{i}) = load([Folder Files{i}])
 end
 
 save('C:\Users\Ryan\Desktop\Kiel MEG Data\Data\Nov21_opm.mat', 'nov21_opm');
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Nov 28 data
 Trials = {'WhiteLED_OPMoff_SubjectASInside',...
          'CurrentSource_WhiteLED_15Hz_SubjectAS',...
          'CurrentSource_WhiteLED_7Hz_SubjectAS',...
          'CurrentSource_WhiteLED_10Hz_SubjectAS',...
          'GreenLED_10Hz_SubjectAS',...
          'GreenLED_10Hz_OPMoff_SubjectASInside',...
          'WhiteLed_10Hz_OPMoff_SubjectJInside',...
          'WhiteLED_10Hz_SubjectJ',...
          'WhiteLED_10Hz_SubjectR'};
      
 Folder = 'C:\Users\Ryan\Desktop\Kiel MEG Data\Nov-28-2019_2nd_MEG_Measurements\';
 
 Files = {'dat_white_led_subject_a_off',...
          'dat_white_led_subject_a_15_Hz_on_current_source_on',...
          'dat_white_led_subject_a_7_Hz_on_current_source_on',...
          'dat_white_led_subject_a_10_Hz_on_current_source_on',...
          'dat_green_led_subject_a_10_Hz_on164908',...
          'dat_green_led_subject_a_10_Hz_off164510',...
          'dat_white_led_subject_j_10_Hz_off_160818',...
          'dat_white_led_subject_j_10_Hz_on_160414',...
          'dat_white_led_subject_r_10_Hz_on_154528'};
      
for i = 1:length(Files)
     nov28_opm.(Trials{i}) = load([Folder Files{i}])
end
 
 save('C:\Users\Ryan\Desktop\Kiel MEG Data\Data\Nov28_opm.mat', 'nov28_opm');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Dec 16 Data
load('C:\Users\Ryan\Desktop\Kiel MEG Data\Jan-30 Plots\matlabplotscodechristmasbreak2019\ChristinEEGvaluesAll')
load('C:\Users\Ryan\Desktop\Kiel MEG Data\Jan-30 Plots\matlabplotscodechristmasbreak2019\RyanEEGvaluesAll')

dec16_eeg.cb.checker = ChristinEEGvaluesAll(:,1:4);
dec16_eeg.cb.led = ChristinEEGvaluesAll(:,5:8);

dec16_eeg.rt.checker = RyanEEGvaluesAll(:,1:4);
dec16_eeg.rt.led = RyanEEGvaluesAll(:,5:8);

save('C:\Users\Ryan\Desktop\Kiel MEG Data\Data\Dec16_eeg.mat', 'dec16_eeg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jan 30 Data
jan30_eeg = load('C:\Users\Ryan\Desktop\Kiel MEG Data\Jan-30 Plots\Jan-30 Data\AS_EEGtrials_Scaled1000.mat');

save('C:\Users\Ryan\Desktop\Kiel MEG Data\Data\Jan30_eeg.mat', 'jan30_eeg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Feb-6 Data
feb6_eeg = load('C:\Users\Ryan\Desktop\Kiel MEG Data\Data\Feb6_AS_EEGtrials_Scaled1000.mat');
save('C:\Users\Ryan\Desktop\Kiel MEG Data\Data\Feb6_eeg.mat', 'feb6_eeg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alldata.nov12_opm_noise = nov12_opm_noise;
alldata.nov21_opm = nov21_opm;
alldata.nov28_opm = nov28_opm;
alldata.dec16_eeg = dec16_eeg;
alldata.jan30_eeg = jan30_eeg;
alldata.feb6_eeg = feb6_eeg;

save('C:\Users\Ryan\Desktop\Kiel MEG Data\Data\alldata.mat', 'alldata', '-v7.3');










