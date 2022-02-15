x = CeegMat1;
v = var(CeegMat1);

repx = repmat(x, 48000/960, 1);
Fs = length(CeegMat1)/0.3;

spectx= fft(CeegMat1);
dBspectx = db(spectx);

m = mean(spectx);
N = length(CeegMat1);

noise_x = repx' + v*randn(1,length(repx));
% noise_x = repx(1:960)' + v*randn(1,length(x)); noise_x = [noise_x repx(961:end)'];

noise_epochs = zeros(960,(48000/960));
for i = 1:(48000/960)
   noise_epochs(:,i) = repx((1:960)+((i-1)*960))+v*randn(length(repx(1:960)),1);
end
zz = mean(noise_epochs, 2)

% 6.505 dB @ 167 Hz, SNR 0.42 dB
figure; snr(zz, Fs);
title('SNR of Christin Trial 1 - Copied 50 trials, added noise each trial, averaged - SNR 0.42 dB');

% 6.664 dB @ 200 Hz, SNR -3.18 dB
figure(1); snr(repx, Fs);
title('SNR of Christin Trial 1 - Copied 50 trials, no noise - SNR -3.18 dB');

% 6.314 dB @ 200 Hz, SNR -17.03 dB
snr(noise_x', length(repx));
title('SNR of Christin trial 1 - Copied 50 trials, added noise - SNR -17.03 dB');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result = (repx.^2) ./ ((noise_x').^2);
figure; plot(db(result))
xlim([0 5]);

frepx = (fft(repx)).^2;
fnoise_x = (fft(noise_x')).^2;
fresult = (frepx) ./ (fnoise_x);
figure; plot(db(fresult));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sampleRate = 3200;
f0 = sampleRate/20;

% 0 to 300 milliseconds
time = 0:1/(sampleRate-1):0.3;

% noise = noise / norm(noise) * norm(signal) / 10^(SNR/20);

signal = CeegMat1;

SNR = 10;
noise = v*randn(size(CeegMat1)) + CeegMat1;
% noise = noise / norm(noise) * norm(signal) / 10^(SNR/20);
x = signal + noise;

actualSNR = 20*log10(norm(signal)/norm(x - signal));
disp(['SNR = ',num2str(actualSNR),'  dB']);

actualSNR = 20*log10(norm(signal)/norm(x - signal));
disp(['SNR = ',num2str(actualSNR),'  dB']) 

figure; subplot(2,1,1)
% plot the signal in time
plot(time,signal)
hold on; grid on
plot(time,x,'r')

% plot the noisy signal in frequency
NFFT = 512*512;
% fftshift - Shift zero-frequency component to center of spectrum
X = fftshift(fft(x,NFFT));
%Normalize X between 1 and 0, never reaches 0
X = X/max(abs(X));
SIG = fftshift(fft(signal,NFFT));
SIG = X/max(abs(SIG));
% 
f = sampleRate/2*linspace(-1,1,NFFT);
subplot(2,1,2)

plot(f,20*log10(abs(X)))
grid on;
hold on; 
plot(f,20*log10(abs(SIG)))
leg = legend('Signal + Noise', 'Signal');

ylim([ -120 5]);
xlim([0 40]);

% Signal + Noise @ 1.7 Hz is -7.11 dB
% Signal @ 1.7 Hz 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; 
noise_x = repx' + v*randn(1,length(repx));
noise_x = noise_x';
[pxx fxx] = pmtm(repx, 4, length(repx), 3200); 
[pxxn fxxn] = pmtm(noise_x, 4, length(repx), 3200);
plot(fxx, db(pxx)); hold on; plot(fxx, db(pxxn));
leg = legend('Signal', 'Signal + Noise');
ylabel('Power (dB/\surdHz)');
xlabel('Frequency (Hz)');

figure; subplot(2,1,1);
% SNR -3.1788, @ 3.4 Hz, -3.299 dB
snr(pxx, fxx, 'psd');

subplot(2,1,2);
% SNR -17.18, @ 3.4 Hz, -3.196 dB
snr(pxxn, fxx, 'psd');

% comparing here dB/Hz at 13.4 Hz, Signal/(Signal+Noise), I get an SNR of
% -0.103 dB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% Frequency @ 1.7 Hz, SNR between Signal & Signal+noise is -2.147 dB 
 
load('C:\Users\Ryan\Downloads\morematlabvariables\Christin_EEG_Vector1')
x = CeegMat1 - min(CeegMat1);
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

for k = 1
    load(['C:\Users\Ryan\Downloads\morematlabvariables\Christin_EEG_Vector' num2str(k)])
    x = CeegMat1 - min(CeegMat1);
    x_n = x + (50/70)*randn(length(repx),1);
    Fs = length(CeegMat1)/0.3;
    repx = repmat(x, 100, 1);
    repx_n = repx + (50/70)*randn(length(repx),1);
    %FFT Power
    jc = abs(fft(x));
    jc_n = abs(fft(x_n));
    % f = 0 : (Fs/(NFFT/2 - 1)) : 1600
    f = 0:(Fs/(length(jc)*2 - 1)):1600;
    df = Fs/(length(jc)*2 - 1);
    % Round up f values
    Ndec = 4;
    rou = 10.^Ndec;
    f = round(rou*f)/rou;
    f = f';
    % Take out every other point that is 0
    jc([find(jc==0)]) = [];

    snrat = zeros(length(f),1);
    bin = 0.7;

    % finding the SNR for a specific frequency - length(x) is too low for
    % frequency resolution
    for i = 1:length(f)
        % frequency measured is the index of vector f
        freq = f(i);
        % Make a vector of index values in the frequency bin.
        % Bin*2 gives bandwidth to measure the signal frequency against.
        vec = find(f <= (freq + bin) & f >= (freq - bin));
        if i == 1
            snrat(i) = jc(i) / (sum(jc(vec(find(vec > i)))));
        elseif i == length(f)
            snrat(i) = jc(length(f))/ (sum(jc(vec(find(vec < i)))));
        else
        % SNR equation = (FFT value of specific frequency) / (FFT Values of Vec)
            snrat(i) = jc(find(f==freq))/( (sum(jc_n(vec(find(f==freq) > vec)))) + (sum(jc_n(vec(find(f==freq) < vec)))) );
        end
    end

    snrat(length(f)) = 0;

    figure; plot(f, db(snrat));
    xlim([0 20]);
    xlabel('Frequency (Hz)');
    ylabel('SNR (dB)');
    ylim([-inf inf]);
    title(['SSVEP SNR - Christin Trial ' num2str(k)]);
end
















