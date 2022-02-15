
%% We'll make a signal with sinewaves at 2 frequencies.
%Signal is 10 seconds long, first 5 seconds has 10 Hz while the last 5
%seconds has 25 Hz.
sampRate = 1000;
t = 0:1/sampRate:10;
t = t';

signal = sin(2*pi*10*t);
tLast = t(5001:end);
signal(5001:end) = sin(2*pi*25*tLast);

%I'll add some noise.
signal = signal + 1*randn(length(t),1);

%Let's plot the raw signal,  first 5 seconds has slow oscillations (10 Hz)
%while the last 5 seconds has faster oscillations (25 Hz).
figure; plot(t,signal);
xlabel('Time (sec)');
ylabel('amplitude (something units)');

%% We'll make a class vector that corresponds to the frequency changes.
%Lets say we have a class vector of 1's and 2's that corrresponds to the
%changes in the sine wave frequency.
classVector = ones(length(t),1);
classVector(5001:end) = 2;
figure; plot(t,classVector);

%% Parse the signal and class vector into time windows.

windowLength = 500;
windowShift = 100;
signalLength = length(signal);

%This estimates the number of time windows we can get from the entire
%signal.
% Should be 1 added
numWindows = floor((signalLength-windowLength)/windowShift) + 1;

%Initialize the cell array for the signal and the class vector.  Each cell
%will contain a time window of data.
dataWindows = cell(numWindows,1);
classWindows = cell(numWindows,1);

% Get the signal into windows first...
firstIDX = 1;
for windowIDX = 1:numWindows
    getIDX = firstIDX:firstIDX+windowLength-1;
    dataWindows{windowIDX} = signal(getIDX);
    firstIDX = firstIDX+windowShift;
end


%Do the same with the class windows.
firstIDX = 1;
for windowIDX = 1:numWindows
    getIDX = firstIDX:firstIDX+windowLength-1;
    classWindows{windowIDX} = classVector(getIDX);
    firstIDX = firstIDX+windowShift;
end

%We can also do the same with the time vector.  This could help us keep
%track of "when" each time window is located in seconds.
timeWindowsSeconds = cell(numWindows,1);
firstIDX = 1;
for windowIDX = 1:numWindows
    getIDX = firstIDX:firstIDX+windowLength-1;
    timeWindowsSeconds{windowIDX} = t(getIDX);
    firstIDX = firstIDX+windowShift;
end

%So now at this point, we have cell arrays with the time windows of the
%signal, the class vector, and the time stamps.

%% Generating a PSD from the fast fourier transform. [SINGLE TIME WINDOW].

%We'll just look at data from a single window for now, this would need to be
%repeated for all windows.  
singleWindowSignal = dataWindows{1};

tWindow = 0:1/sampRate:windowLength/sampRate - 1/sampRate;
figure; plot(tWindow,singleWindowSignal);
xlabel('Time (sec)');
ylabel('amplitude (arbitrary)');


%Actual FFT portion.
nfft = 1028;
%Evaluate the FFT
fftSolution = fft(singleWindowSignal,nfft);
%Normalize it by the signal's power
fftSolution = fftSolution/windowLength;
%Convert the complex solution to a real one (We square the real and
%imaginary parts).
fftAbs = abs(fftSolution);

%Convert to the one sided solution.
PSD = fftAbs(1:nfft/2+1);
PSD(2:end-1) = 2*PSD(2:end-1).^2;

%Frequency vector.
freq = 0:sampRate/nfft:sampRate/2;

%Plotted here is the final PSD.
figure; plot(freq,PSD);
xlabel('Frequency (Hz)');
ylabel('Power (something squared / Hz)');

%% Alpha and Beta power extraction. [SINGLE TIME WINDOW].
%In this single time window of data, we have a PSD.  We need to extract
%alpha and beta band power. To do this, we need to integrate the PSD across the
%frequency band of interest (alpha and beta).

freq = 0:sampRate/nfft:sampRate/2;
figure; plot(freq,PSD,'color','k');
xlabel('Frequency (Hz)');
ylabel('Power (something squared / Hz)');
%The PSD does extend to 500 Hz, but I'm going to focus on <50 Hz for now.)
xlim([0 50]);
ylim([-0.5 1.5]);

rectObjectAlpha = rectangle('Position',[8 -0.5 6 6],'FaceColor',[0 1 1 .5],...
    'EdgeColor','none');

rectObjectBeta = rectangle('Position',[20 -0.5 10 6],'FaceColor',[1 0 0 .5],...
    'EdgeColor','none');

alphaText = text(10,1.4,{'Sum PSD values',' here for alpha'},...
    'HorizontalAlignment','center');

betaText = text(25,1.4,{'Sum PSD values',' here for beta'},...
    'HorizontalAlignment','center');

alphaIDX = find(freq>8 & freq<13);
betaIDX = find(freq>20 & freq<30);

%Here I'm integrating across the frequnce bands (or taking the sum)). Keep
%in mind after this we get a single value of alpha and beta from a single
%time window.
singleWindow_alpha = sum(PSD(alphaIDX));
singleWindow_beta = sum(PSD(betaIDX));

%% Now perform the FFT and power extrction for all time windows.

%Now we're making a new time series.  Keep in mind that with the time
%window size and the window shift, we're effectively making a new signals
%that represent alpha and beta band power, but at a rate of 10 Hz. Note
%that the number of elements is the number of time windows we used.
alphaPowerTimeSeries = zeros(numWindows,1);
betaPowerTimesSeries = zeros(numWindows,1);

for windowIDX = 1:numWindows
    
    %Get the signal from a single time window.  Note the window index.
    singleWindowSig = dataWindows{windowIDX};
    
    %------FFT Stuff -------------------
    %Evaluate the FFT
    fftSolution = fft(singleWindowSig,nfft);
    %Normalize it by the signal's power
    fftSolution = fftSolution/windowLength;
    %Convert the complex solution to a real one (We square the real and
    %imaginary parts).
    fftAbs = abs(fftSolution);
    %Convert to the one sided solution.
    PSD = fftAbs(1:nfft/2+1);
    PSD(2:end-1) = 2*PSD(2:end-1).^2;
    
    %-----Band power extraction ------------- 
    %Frequency vector.
    freq = 0:sampRate/nfft:sampRate/2;
    
    %Indecies corresponding to alpha and beta frequencies.
    alphaIDX = find(freq>8 & freq<13);
    betaIDX = find(freq>20 & freq<30);

    %Here I'm integrating across the frequnce bands (or taking the sum)). Keep
    %in mind after this we get a single value of alpha and beta from a single
    %time window.
    singleWindow_alpha = sum(PSD(alphaIDX));
    singleWindow_beta = sum(PSD(betaIDX));
    
    %We store the single spectral power value to the intialized times
    %series above this loop.  
    alphaPowerTimeSeries(windowIDX) = singleWindow_alpha;
    betaPowerTimesSeries(windowIDX) = singleWindow_beta;
end


%% We'll extract the last class sample in each time window.
%They will match that of the spectral power time series.

classVectorSync = zeros(numWindows,1);
for windowIDX = 1:numWindows
    classVectorSync(windowIDX) = classWindows{windowIDX}(end);
end

%% We'll do the same for the time vector.
%This will give us a sense of how the new power signal changes at a slower
%rate.  I'll also use this later to help plot the power time series with
%the raw signal we had from earlier.

timeWindowTimeStamps = zeros(numWindows,1);
for windowIDX = 1:numWindows
    timeWindowTimeStamps(windowIDX) = timeWindowsSeconds{windowIDX}(end);
end

%% Now let's plot it all together to see what we have in the end.

figure; hold on;
plot(t,signal,'k');
plot(timeWindowTimeStamps,alphaPowerTimeSeries,'b',...
    'LineWidth',2);
plot(timeWindowTimeStamps,betaPowerTimesSeries,'r',...
    'LineWidth',2);
legend('Raw Signal','Alpha Power','Beta Power',...
    'Location','southeast');

%So we got the raw signal, the time series for the alpha and beta band
%power.  As expected the first phase is purely 10 Hz, so we get high alpha
%in the beginning.  The last phase is purely 25 Hz so we get high beta
%in the end.


%Now I'll plot them as points rather than lines to emphasize this point.
%The way we used time windows makes the power time series vary at a slower
%rate than the raw signal we started with.  So this is the "new sampling
%rate" I mean where initially it was 1000 Hz, but with a shift in 100
%samples, the effective rate is now 10 Hz.
figure; hold on;
plot(t,signal,'k.');
plot(timeWindowTimeStamps,alphaPowerTimeSeries,'bo',...
    'LineWidth',2);
plot(timeWindowTimeStamps,betaPowerTimesSeries,'ro',...
    'LineWidth',2);
legend('Raw Signal','Alpha Power','Beta Power',...
    'Location','southeast');
xlim([4 6]);


%To help emphasize this point, I'll show you the raw class vector we had
%before, and the "windowed" class vector.  Keep in mind that the windowed
%class vector is synchronized with each point in the Power time series.
figure; hold on;
plot(t,classVector,'k.');
plot(timeWindowTimeStamps,classVectorSync,'bo',...
    'LineWidth',2);
legend('Raw class','extracted windowed classes',...
    'Location','southeast');
xlim([4 6]);

%The "powerTimeSeries" variables are what I intended to be fed into the
%input for the classifier. Note that here, the "classVectorSync" variable
%has the class labels synchronized with each power estimate. Using these
%features here should give you some decent accuracies, but there's an extra
%smoothing step you can take, which is shown next.
%% Now here is where you can apply the "smoothing" filters.
%I intended it to used on the spectral power time series.
%This is why I use the "new sampling rate" (10 Hz) rather than the sampling
%rate of the raw signal (1000 Hz).  I'm trying to filter the power time series, not
%the raw signal.

%Note that this is done wayyy after the band power extraction and quite a
%long ways away from the FFT step.

%Design a butterworth low pass filter at 2 Hz with an order of 1.
[b,a] = butter(1,2/(10/2),'low');

alphaFiltered = filtfilt(b,a,alphaPowerTimeSeries);
betaFiltered = filtfilt(b,a,betaPowerTimesSeries);

%I'll now plot the power series before and after low pass filtering.
figure; hold on;
plot(timeWindowTimeStamps,alphaPowerTimeSeries,'b',...
    'LineWidth',2);
plot(timeWindowTimeStamps,betaPowerTimesSeries,'r',...
    'LineWidth',2);
plot(timeWindowTimeStamps,alphaFiltered,'g',...
    'LineWidth',2);
plot(timeWindowTimeStamps,betaFiltered,'m',...
    'LineWidth',2);
legend('alpha','beta','alpha filtered','beta filtered',...
    'Location','east');

%We can see that smoothing it out just helps suppress some of the faster
%spikes.  We can probably go down to 1 Hz too, but it doesn't make it as
%"snappy" or responsive.  In any case, you can see how it just helps a
%little bit and isn't super necessary.

%% Done!

%These "alphaPowerTimeSeries" is probably the main thing that should be put
%into the classifier.  Note that the Power Time series Varabile, the
%classVectorSync and timeWindowTimeStamps all have the same number of
%samples and that they're properly synchronized to relate to the same time
%windows.

%Also note that the values in timeWindowTimeStamps are indeed the same as
%the raw time stamps we got before (in seconds). You could use them to help
%crop data from portions of trials.  I can find the window indecies
%corresponding to the first three seconds based on this variable.

