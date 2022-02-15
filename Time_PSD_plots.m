

i=1;
CHANNEL = 2;
Fs = 10000;
Matfile = fullfile(folder, [file{i}  '.mat']);
x{i} = load(Matfile);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FILE = 0;
CHANNEL = 2;
TRIALS = [1 2];
EPOCH = 1;
PSD = 0;
SPECT = 0;
FILTER = 1;
EPOCH_ANALYIS = 1;
EPOCH_PSD = 1;
PSD_NA = 1;

%can define the threshold here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if FILE == 0
 folder = 'C:\Users\Ryan\Desktop\Kiel MEG Data\Nov-12 Noise Measurement with LED';

 file = {'dat_meg_measurement_165544', 'dat_meg_measurement_170159'};
 
 titlename = {"MEG Noise Measurement at 16:55 Nov-12", ...
              "MEG Noise Measurement at 17:01 Nov-12"};

elseif FILE == 1
 folder = 'C:\Users\Ryan\Desktop\Kiel MEG Data\Nov-21-2019_Initial_MEG_Measurement';
    
 file = {'dat_Human_inside_fiber_optics_off_person_1_7_Hz_3min133118', ...
         'dat_Human_inside_fiber_optics_on_person_1_7Hz_3min133833', ...
         'dat_Human_inside_white_led_on_person_1_7_Hz_3min134811', ...
         'dat_Human_inside_green_led_on_person_1_7_Hz_3min135530', ...
         'dat_Human_inside_green_led_on_person_1_3min_10Hz140043', ...
         'dat_noise_measurement_10_Hz_green140838', ...
         'dat_dat_Human_inside_fiber_optics_off_person_2_10Hz_3min143008', ...
         'dat_dat_Human_inside_fiber_optics_on_person_2_10Hz_3min143408', ...
         'dat_dat_Human_inside_fiber_optics_on_noise_measurment144252'};
     
 titlename = {"Fiber Optic Blinking 7Hz & OPM Off Person 1 [13:31] Nov-21", ...
              "Fiber Optic Blinking 7Hz & OPMs On Person 1 [13:38] Nov-21", ...
              "White LED Blinking 7Hz & OPMs On Person 1 [13:48] Nov-21", ...
              "Green LED Blinking 7Hz & OPMs On Person 1 [13:55] Nov-21", ...
              "Green LED Blinking 10Hz & OPMs On Person 1 [14:04] Nov-21", ...
              "Green LED Blinking 10Hz & OPMs Noise Measurement [14:08] Nov-21", ...
              "Fiber Optic Blinking 10Hz & OPMs Off Person 2 [14:30] Nov-21", ...
              "Fiber Optic Blinking 10Hz & OPMs On Person 2 [14:34] Nov-21", ...
              "Fiber Optic Blinking 10Hz & OPMs Off Person 2 [14:42] Nov-21"};
     
elseif FILE == 2
folder =  'C:\Users\Ryan\Desktop\Kiel MEG Data\Nov-28-2019_2nd_MEG_Measurements';
    
 file = {'dat_white_led_subject_r_10_Hz_on_154528', ...
         'dat_white_led_subject_j_10_Hz_on_160414', ...
         'dat_white_led_subject_j_10_Hz_off_160818', ...
         'dat_green_led_subject_a_10_Hz_off164510', ...
         'dat_green_led_subject_a_10_Hz_on164908', ...
         'dat_white_led_subject_a_10_Hz_on_current_source_on', ...
         'dat_white_led_subject_a_7_Hz_on_current_source_on', ...
         'dat_white_led_subject_a_15_Hz_on_current_source_on', ...
         'dat_white_led_subject_a_off'};
     
 titlename = {"White LED Blinking 10Hz & OPMs on Subject R [15:45] Nov-28", ...
              "White LED Blinking 10Hz & OPMs on Subject J [16:04] Nov-28", ...
              "White LED Off & OPMs On Subject J [16:08] Nov-28", ...
              "Green LED Off & OPMs on Subject A [16:45] Nov-28", ...
              "Green LED Blinking 10Hz & OPMs on Subject A [16:49] Nov-28", ...
              "Current Source - White LED Blinking 10Hz & OPMs on Subject A Nov-28", ...
              "Current Source - White LED Blinking 7Hz & OPMs on Subject A Nov-28", ...
              "Current Source - White LED Blinking 15Hz & OPMs on Subject A Nov-28", ...
              "White LED Off & OPMs on Subject A Nov-28"};
          
 legNames = {'Arduino - W LED 10 Hz - Subject R', 'Arduino - W LED 10 Hz - Subject J', ...
             'Arduino - W No Blinking - Subject J', 'Arduino - G No Blinking - Subject A', ...
             'Arduino - G LED 10 Hz', 'Current Source - 10 Hz', 'Current Source - 7 Hz', ...
             'Current Source - 15 Hz', 'Current Source - No Blinking'};
end

ChannelName = {'Y ', 'V ', 'X ', 'W '};

pre_pxx_noise = zeros(length(file),1);
pre_fxx_noise = zeros(length(file),1);

    if FILE == 0
        for d = 1:2
            Matfile = fullfile(folder, [file{d}  '.mat']);
            x{d} = load(Matfile);
            
            [hp2, hp1] = butter(3, 1/(Fs/2), 'high');
            [lp2, lp1] = butter(3, 40/(Fs/2), 'low');
            x{d}.x.filter.hp = zeros(length(x{d}.x.b), length(x{d}.x.b(1,:)));
            x{d}.x.filter.hp_lp = zeros(length(x{d}.x.b), length(x{d}.x.b(1,:)));
                for j = 1:4
                    x{d}.x.filter.hp(:, j) = filtfilt(hp2, hp1, x{d}.x.b(:, j));
                    x{d}.x.filter.hp_lp(:, j) = filtfilt(hp2, hp1, x{d}.x.filter.hp(:, j));
                 end
            [pre_pxx_noise{d} pre_fxx_noise{d}] = pwelch(x{d}.x.filter.hp_lp(:,CHANNEL), 10000, 9500, 512*512, Fs);
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i = 1;
Fs = 10000;
FILE = 1;
CHANNEL = 2;
TRIALS = [1 2];
EPOCH = 1;
PSD = 0;
SPECT = 0;
FILTER = 1;
EPOCH_ANALYIS = 1;
EPOCH_PSD = 1;
PSD_NA = 1;

%can define the threshold here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if FILE == 0
 folder = 'C:\Users\Ryan\Desktop\Kiel MEG Data\Nov-12 Noise Measurement with LED';

 file = {'dat_meg_measurement_165544', 'dat_meg_measurement_170159'};
 
 titlename = {"MEG Noise Measurement at 16:55 Nov-12", ...
              "MEG Noise Measurement at 17:01 Nov-12"};

elseif FILE == 1
 folder = 'C:\Users\Ryan\Desktop\Kiel MEG Data\Nov-21-2019_Initial_MEG_Measurement';
    
 file = {'dat_Human_inside_fiber_optics_off_person_1_7_Hz_3min133118', ...
         'dat_Human_inside_fiber_optics_on_person_1_7Hz_3min133833', ...
         'dat_Human_inside_white_led_on_person_1_7_Hz_3min134811', ...
         'dat_Human_inside_green_led_on_person_1_7_Hz_3min135530', ...
         'dat_Human_inside_green_led_on_person_1_3min_10Hz140043', ...
         'dat_noise_measurement_10_Hz_green140838', ...
         'dat_dat_Human_inside_fiber_optics_off_person_2_10Hz_3min143008', ...
         'dat_dat_Human_inside_fiber_optics_on_person_2_10Hz_3min143408', ...
         'dat_dat_Human_inside_fiber_optics_on_noise_measurment144252'};
     
 titlename = {"Fiber Optic Off & OPM On Person 1 [13:31] Nov-21", ...
              "Fiber Optic Blinking 7Hz & OPMs On Person 1 [13:38] Nov-21", ...
              "White LED Blinking 7Hz & OPMs On Person 1 [13:48] Nov-21", ...
              "Green LED Blinking 7Hz & OPMs On Person 1 [13:55] Nov-21", ...
              "Green LED Blinking 10Hz & OPMs On Person 1 [14:04] Nov-21", ...
              "Green LED Blinking 10Hz & OPMs Noise Measurement [14:08] Nov-21", ...
              "Fiber Optic Blinking 10Hz & OPMs Off Person 2 [14:30] Nov-21", ...
              "Fiber Optic Blinking 10Hz & OPMs On Person 2 [14:34] Nov-21", ...
              "Fiber Optic Off & OPMs On Person 2 [14:42] Nov-21"};
          
legNames = {'Fiber Optic Off Subject 1', 'Fiber Optic W LED 7 Hz - Subject 1', ...
            'Arduino - W LED 7 Hz - Subject 1', 'Arduino - G LED 7 Hz - Subject 1', ...
            'Arduino - G LED 10 Hz Subject 1', 'Arduino - G LED 10 Hz Noise Measurement', 'Fiber Optic W LED 10 Hz - Subject A', ...
            'Fiber Optic W LED 10 Hz - Subject A', 'Fiber Optic W LED 10 Hz - Noise Measurement'};
     
elseif FILE == 2
folder =  'C:\Users\Ryan\Desktop\Kiel MEG Data\Nov-28-2019_2nd_MEG_Measurements';
    
 file = {'dat_white_led_subject_r_10_Hz_on_154528', ...
         'dat_white_led_subject_j_10_Hz_on_160414', ...
         'dat_white_led_subject_j_10_Hz_off_160818', ...
         'dat_green_led_subject_a_10_Hz_off164510', ...
         'dat_green_led_subject_a_10_Hz_on164908', ...
         'dat_white_led_subject_a_10_Hz_on_current_source_on', ...
         'dat_white_led_subject_a_7_Hz_on_current_source_on', ...
         'dat_white_led_subject_a_15_Hz_on_current_source_on', ...
         'dat_white_led_subject_a_off'};
     
 titlename = {"White LED Blinking 10Hz & OPMs on Subject R [15:45] Nov-28", ...
              "White LED Blinking 10Hz & OPMs on Subject J [16:04] Nov-28", ...
              "White LED Off & OPMs On Subject J [16:08] Nov-28", ...
              "Green LED Off & OPMs on Subject A [16:45] Nov-28", ...
              "Green LED Blinking 10Hz & OPMs on Subject A [16:49] Nov-28", ...
              "Current Source - White LED Blinking 10Hz & OPMs on Subject A Nov-28", ...
              "Current Source - White LED Blinking 7Hz & OPMs on Subject A Nov-28", ...
              "Current Source - White LED Blinking 15Hz & OPMs on Subject A Nov-28", ...
              "White LED Off & OPMs on Subject A Nov-28"};
          
 legNames = {'Arduino - W LED 10 Hz - Subject R', 'Arduino - W LED 10 Hz - Subject J', ...
             'Arduino - W No Blinking - Subject J', 'Arduino - G No Blinking - Subject A', ...
             'Arduino - G LED 10 Hz', 'Current Source - 10 Hz', 'Current Source - 7 Hz', ...
             'Current Source - 15 Hz', 'Current Source - No Blinking'};
end

Matfile = fullfile(folder, [file{i}  '.mat']);
x{i} = load(Matfile);

figure; subplot(2, 2, 1); 
t = (0:length(x{i}.x.b)-1)/Fs;
plot(t, x{i}.x.b(:, CHANNEL));
title(sprintf('Before Preprocessing - Time Domain Signal of Trial %s', num2str(i)));

[pxx fxx] = pwelch(x{i}.x.b(:,CHANNEL), 10000, 9500, 512*512, Fs);
subplot(2, 2, 2); 
plot(fxx, db(pxx/2));
title(sprintf('Before Preprocessing - Welch PSD of Trial %s', num2str(i)));
xlim([0 100]);

subplot(2, 2, 3);
% Parameters for zero-phase butterworth notch, highpass, and lowpass filters
    [hp2, hp1] = butter(3, 1/(Fs/2), 'high');
    [lp2, lp1] = butter(3, 40/(Fs/2), 'low');

     x{i}.x.filter.hp = zeros(length(x{i}.x.b), length(x{i}.x.b(1,:)));
     x{i}.x.filter.hp_lp = zeros(length(x{i}.x.b), length(x{i}.x.b(1,:)));
    for j = 1:4
        x{i}.x.filter.hp(:, j) = filtfilt(hp2, hp1, x{i}.x.b(:, j));
        x{i}.x.filter.hp_lp(:, j) = filtfilt(hp2, hp1, x{i}.x.filter.hp(:, j));
    end
    
    plot(t, x{i}.x.filter.hp_lp(:, CHANNEL));
    title(sprintf('After Preprocessing - Time Domain Signal of Trial %s', num2str(i)));

    
    [pre_pxx pre_fxx] = pwelch(x{i}.x.filter.hp_lp(:,CHANNEL), 10000, 9500, 512*512, Fs);
    
    subplot(2, 2, 4);
    p1 = plot(pre_fxx, db(pre_pxx)/2);
    p1.LineWidth = 0.75; p1.Color = 'blue';
    hold on; 
    p2 = plot(pre_fxx_noise{1}, db(pre_pxx_noise{1})/2);
    p2.LineWidth = 0.5; p2.Color = 'red'; p2.LineStyle = '--';
    p3 = plot(pre_fxx_noise{2}, db(pre_pxx_noise{2})/2);
    p3.LineWidth = 0.5; p3.Color = 'magenta'; p3.LineStyle = '-.';;
%     p4 = plot([10;10], [ylim]'); 
%     p4.LineWidth = 1.0; p4.Color = 'green';
%     h = text([7.3], [-280], {'Blinking', 'Frequency'}, 'VerticalAlignment', 'bottom',...
%         'HorizontalAlignment', 'left');
%     set(h, 'Rotation', 270);
    hold off;
    xlim([0 100]); 
    ylim([-320 -200]);
    legend([
    title(sprintf('After Preprocessing - Welch PSD of Trial %s', num2str(i)));
    
    suptitle(sprintf('Time Domain & PSD Plots of %s', titlename{i}));
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    