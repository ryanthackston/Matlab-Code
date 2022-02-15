%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all, clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tic Toc method
frequency = 10;
fs = 10e3;
duration = 20;
samples = 1:1:duration*fs;
t = (0:length(samples)-1)*1/fs;
signal =sin(2*pi*frequency*t);
% Create checkerboard image
m = 100; n = 160;
[C, K] = checkerb(m, n);
h = imshow(C, 'Border', 'tight', 'InitialMagnification', 'fit');
f = 1/(frequency*2);
trigger = zeros(8000000, 1);
% k is the number of periods/cycles for fully blinking on then off
k=1;
% n is the index of cycles in the tic toc function
n = 1;

while k < (frequency * duration)
        tic
        % d - pause time
        h.CData = K;
        drawnow('expose')
        while(toc < f)
            trigger(n) = 1;
            n = n+1;
        end
        
        tic
        h.CData = C;
        drawnow('expose')
        while(toc < f)
            trigger(n) = 2;
            n = n+1;
        end
        k = k+1;
end

del_trig = find(trigger == 0);
trigger(del_trig(1:end)) = [];

z = trigger;
zs = z(1:(length(z)/(1000*duration)):length(z));
thresh = mean(zs);
tictoc_change = zeros(length(zs),1);
 for c = 1:length(zs)
     if c==1
        tictoc_change(c) = 0;
        continue;
     elseif zs(c-1) == zs(c)
        tictoc_change(c-1) = 0;
     elseif zs(c-1) ~= zs(c)
        tictoc_change(c-1) = 1;
     elseif c == length(zs)
        tictoc_change(c) = 0;
     end
 end

tictoc_times = find(tictoc_change == 1);
for y = 2:length(tictoc_times)
    if y == 1
        tictoc_times(y,2) = 0;
    else
        tictoc_times(y,2) = tictoc_times(y,1) - tictoc_times(y-1,1);
    end
end

clearvars trigger z signal samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% frequency = 10;
% fs = 10e3;
% duration = 20;
% samples = 1:1:duration*fs;
% t = (0:length(samples)-1)*1/fs;
% signal =sin(2*pi*frequency*t);
% 
% 
% 
% %rectang is the trigger function
% rectang = zeros(length(signal), 1);
% for i = 1:length(rectang)
%     if signal(i) >= 0
%         rectang(i)=1;
%     else
%         rectang(i)=-1;
%     end
% end
% 
% 
% m = 100; n = 160;
% 
% 
% % i = 1;
% % while (i < length(rectang))    
% %     if rectang(i)==1
% %         h.CData = K;
% %         drawnow('expose')
% %     else
% %         h.CData = C;
% %         drawnow('expose')
% %     end    
% %     i= i+1;
% % end
% [C, K] = checkerb(m, n);
% h = imshow(C, 'Border', 'tight', 'InitialMagnification', 'fit');
% f = 1/(frequency*2);
% % k is the number of periods/cycles
% trigger = zeros(8000000, 1);
% k=1;
% n = 1;
% while k < (frequency * duration)
%         B = cputime;
%         % d - pause time
%         d = B + f;
%         h.CData = K;
%         drawnow('expose')
%         while(cputime < d)
%             trigger(n) = 1;
%             n = n+1;
%         end
%         
%         B = cputime;
%         d = B + f;
%         h.CData = C;
%         drawnow('expose')
%         while(cputime < d)
%             trigger(n) = 2;
%             n = n+1;
%         end
%         k = k+1;
% end
% 
% del = find(trigger == 0);
% trigger(del(1:end)) = [];
% 
% z = trigger;
% zs = z(1:(length(z)/(1000*duration)):length(z));
% thresh = mean(zs);
% cputime_change = zeros(length(zs),1);
%  for c = 1:length(zs)
%      if c==1
%         cputime_change(c) = 0;
%         continue;
%      elseif zs(c-1) == zs(c)
%         cputime_change(c-1) = 0;
%      elseif zs(c-1) ~= zs(c)
%         cputime_change(c-1) = 1;
%      elseif c == length(zs)
%         cputime_change(c) = 0;
%      end
%  end
% 
% cputime_times = find(cputime_change == 1);
% for y = 2:length(cputime_times)
%     if y == 1
%         cputime_times(y,2) = 0;
%     else
%         cputime_times(y,2) = cputime_times(y,1) - cputime_times(y-1,1);
%     end
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
















