data = {AS750.F7o5.T1_Frame2_Epoch; AS750.F7o5.T2_Frame2_Epoch; AS750.F7o5.T3_Frame2_Epoch; ...
        AS750.F7o5.T4_Frame2_Epoch; AS750.F7o5.T5_Frame2_Epoch; AS750.F7o5.T6_Frame2_Epoch; ...
        AS750.F15.T1_Frame1_Epoch; AS750.F15.T2_Frame1_Epoch; AS750.F15.T3_Frame1_Epoch; ...
        AS750.F15.T4_Frame1_Epoch; AS750.F15.T5_Frame1_Epoch; AS750.F15.T6_Frame1_Epoch; 
        AS750.F1o9.T1_Frame2_Avg; AS750.F1o9.T2_Frame2_Avg; ...
        AS750.F10.T1_Frame2_Avg; AS750.F10.T1_Frame2_Avg; ...
        AS750.F10.T1_Frame1_Avg; AS750.F10.T1_Frame1_Avg}


[pxx fxx] = pwelch(repx, hann(100), 95, 512*512, 3200);
[pxxn fxxn] = pwelch(noise_x, hann(100), 95, 512*512, 3200);

figure; subplot(2,1,1);
plot(fxx, db(pxx)); hold on; plot(fxx, db(pxxn));
leg = legend('Signal', 'Signal + Noise');
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');
% Frequency @ 1.7 Hz, SNR between Signal & Signal+noise is -2.113 dB

[pxx fxx] = pwelch(repx, hann(3200), [], 2*3200, 3200);
[pxxn fxxn] = pwelch(noise_x, hann(3200), [], 2*3200, 3200);
subplot(2,1,2)
plot(fxx, db(pxx)); hold on; plot(fxx, db(pxxn));
leg = legend('Signal', 'Signal + Noise');
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency @ 1.7 Hz, SNR between Signal & Signal+noise is -2.147 dB 
%  
% load('C:\Users\Ryan\Downloads\morematlabvariables\Christin_EEG_Vector1')
x = ASmatrix.AS_1o7Hz_Trial1 - min(ASmatrix.AS_1o7Hz_Trial1);
repx = repmat(x, 48000/960, 1);
Fs = length(CeegMat1)/0.3;

repx_noise = repx + (50/70)*randn(length(repx),1);
%FFT Power
j = abs(fft(repx)).^2;
jc = abs(fft(x));
jc_noise = abs(fft(x + (50/70)*randn(length(x),1)));
% f = 0 : (Fs/(NFFT/2 - 1)) : 1600
f = 0:(Fs/(length(j) - 1)):1600;
df = Fs/(length(j) - 1);
% Round up f values
Ndec = 4;
rou = 10.^Ndec;
f = round(rou*f)/rou;
f = f';
% Take out every other point that is 0
j([find(j==0)]) = [];
jc([find(jc==0)]) = []
jdb = db(j);
v = var(j);
jnoise = db(j + v*(randn(length(j),1)));
%Try normalizing j
jrange = max(j)-min(j);
jnorm = (j - min(j))/jrange;
vnorm = var(jnorm);
jnorm_noise = jnorm + vnorm*(randn(length(jcnorm),1));
%Normalize jc
jcrange = max(jc)-min(jc);
jcnorm = (jc - min(jc))/jcrange;
jcrange_n = max(jc_noise)-min(jc_noise);
jcnorm_n = (jc_noise - min(jc_noise))/jcrange_n;

snrat = zeros(length(fDown),1);
bin = 0.1;

% finding the SNR for a specific frequency
for i = 2:length(fDown)-1
    % frequency measured is the index of vector f
    freq = fDown(i);
    % Make a vector of index values in the frequency bin.
    % Bin*2 gives bandwidth to measure the signal frequency against.
    vec = find(fDown <= freq + bin & fDown >= freq - bin);
    if i == 1
%         snrat(1) = jdb(1)/ (sum(jdb(vec(find(f==freq) > vec))));
%     snrat(1) = jnorm(1)/ (sum(jnorm(vec(find(f==freq) > vec))));
    snrat(1) = jcnorm(1)/ (sum(jcnorm(vec(find(fDown==freq) > vec))));
    elseif i == length(fDown)
%         snrat(length(f)) = jnorm(length(f))/ (sum(jnorm(vec(find(f==freq) < vec))));
     snrat(length(fDown)) = jcnorm(length(fDown))/ (sum(jcnorm(vec(find(fDown==freq) < vec))));
    else
    % SNR equation = (FFT value of specific frequency) / (FFT Values of Vec)
%     snrat(i) = jdb(find(f==freq))/( (sum(jnoise(vec(find(f==freq) > vec)))) + (sum(jnoise(vec(find(f==freq) < vec)))) );
%      snrat(i) = jnorm(find(f==freq))/( (sum(jnorm_noise(vec(find(f==freq) > vec)))) + (sum(jnorm_noise(vec(find(f==freq) < vec)))) );
    snrat(i) = jcnorm(find(fDown==freq))/( (sum(jcnorm_n(vec(find(fDown==freq) > vec)))) + (sum(jcnorm_n(vec(find(fDown==freq) < vec)))) );
    end
end
snrat(length(fDown)) = 0;

freqbase = 1;
for i = 2:length(f)-1
    % frequency measured is the index of vector f
    freq = f(i);
    % snrat(i) = jdb(find(f==freq))/( (sum(jnoise(vec(find(f==freq) > vec)))) + (sum(jnoise(vec(find(f==freq) < vec)))) );
%      snrat(i) = jnorm(find(f==freq))/( (jnorm_noise(find(f==freqbase))) );
     snrat(i) = jcnorm(find(f==freq))/( (jcnorm_n(find(f==freqbase))) );
end

% 
figure; plot(fDown, db(snrat));
xlim([0 20]);
xlabel('Frequency (Hz)');
ylabel('SNR (dB)');

% TD SNR - Not a good method
figure; hold on;
    x = CeegMat1 - min(CeegMat1);
    t = 0:1/(Fs-1):0.3;
    OPM = x*(7e-14);
    OPM_Noise = OPM + (50e-15)*randn(length(OPM), 1);

    p1 = plot(t, OPM, 'k');
    p2 = plot(t, OPM_Noise, 'r');
    ylim([0 2.1e-12]);

    ylabel('Magnetic Field (pT)');
    xlabel('Time (Seconds)');
    title('Time Domain - Christin Trial 1');
    leg = legend('Signal','Signal + 50 fT Noise'); 
    


% freq = 1.7;
% [minValue,closestIndex] = min(abs(bsxfun(@minus,f, freq)));


figure; plot(t, (OPM ./ OPM_Noise)); hold on;
% plot(t, OPM);
% plot(t, OPM_Noise);

TimeSNR = OPM ./ OPM_Noise;
repTimeSNR = repmat(TimeSNR, 48000/960, 1);

% FreqSNR = abs(fft(repTimeSNR));
% FreqSNR(find(FreqSNR == 0)) = [];
% figure; plot(f, db(FreqSNR));

fDown = downsample(f, length(f)/length(OPM));
figure; plot(fDown, db(abs(fft(TimeSNR))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency Domain SNR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    load('C:\Users\Ryan\Downloads\MATLAB Functions\AS_Trials\AS_EEGtrials_Scaled1000')
    
    x = [data{1:6}];
    x = x(:);
    x = (x-min(x));
    x_n = x + (50/70)*randn(length(x),1);
    
    Fs = length(data{1})/0.75;

    %FFT Power
    jc = abs(fft(x));
    jc_n = abs(fft(x_n));
    % Take out every other point that is 0
    jc_n(find(jc == 0)) = [];
    jc(find(jc == 0)) = [];

    f = 0:(Fs/(length(jc)*2 - 1)):(Fs/2);
    df = Fs/(length(jc)*2 - 1);
    % Round up f values
    Ndec = 4;
    rou = 10.^Ndec;
    f = round(rou*f)/rou;
    f = f';

    % snrat is signal-to-noise at specific frequency
    snrat = zeros(length(f),1);
    % bin is how wide you want to go with the surrounding bins. 
    % bin = 0.5 means you add frequency bins up to 0.5 Hz in front and 0.5 Hz
    % behind in the denominator
    bin = 0.5;
    
    %     snrat(1) = 0;

    % Edge artifacts happen when you copy the data several times.
    % Frequency resolution is low with only 960 points in each trial,
    % concatenating trials under the same conditions will help.
    for i = 2:length(f)
        % frequency measured is the index of vector f
        freq = f(i);
        % Make a vector of index values in the frequency bin.
        vec = find(f <= (freq + bin) & f >= (freq - bin));
        if i == 1
            % If i == 1 only add frequency bins in front
            snrat(i) = jc(i) / (sum(jc(vec(find(vec > i)))));
            % If i is the last index, only add frequency bins before
        elseif i == length(f)
            snrat(i) = jc(length(f))/ (sum(jc(vec(find(vec < i)))));
        else
        % SNR equation = (FFT value of specific frequency) / (FFT Values of Vec)
            snrat(i) = jc(find(f==freq))/( (sum(jc_n(vec(find(f==freq) > vec)))) + (sum(jc_n(vec(find(f==freq) < vec)))) );
        end
    end

    snrat(length(f)) = 0;

    figure; plot(f, db(snrat));
    xlim([0 200]);
    xlabel('Frequency (Hz)');
    ylabel('SNR (dB)');
    ylim([-inf inf]);
    title(['SSVEP SNR - AS 1.9 Hz 2-Frame Epoch Trial 13-14']);
    set(gca,'fontsize', 15)
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = AS.F7.T1 - mean(AS.F7.T1);
xnor = x ./ std(x);

repx = repmat(xnor, 100, 1);
jc = abs(fft(repx));

figure; plot(repx);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = [AS.F7.T1; AS.F7.T2; AS.F7.T3];
x = (x - mean(x))/std(x);
% figure; plot(db(abs(fft(x))));
% xlim([0 40])

f = 3200*(0:(length(x)/2))/length(x)


P2 = abs(fft(x)/length(x));
P1 = P2(1:(length(x)/2+1));
P1(2:end-1) = 2*P1(2:end-1);

z = [AS.F1o7.T1; AS.F1o7.T2; AS.F1o7.T3]
z = (z - mean(z))/std(z);
fz = 3200*(0:(length(z)/2))/length(z)
Pz2 = abs(fft(z)/length(z));
Pz1 = Pz2(1:(length(z)/2+1));
Pz1(2:end-1) = 2*Pz1(2:end-1);


% figure; subplot(3,1,1); plot(f, db(P1));
% xlim([0 60]); ylim([-inf inf]);
% subplot(3,1,2); plot(fz, db(Pz1));
% xlim([0 60]); ylim([-inf inf]);

figure; plot(f, db(P1)); hold on;
plot(fz, db(Pz1));
leg = legend ('7Hz', '1.7 Hz');
xlim([0 30]);

xlabel('Frequency (Hz)');
ylabel('Spectral Amplitude');
title('Spectral Amplitude of AS');
































