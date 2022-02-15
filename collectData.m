function collectData(src,event)
     time = event.TimeStamps;
     data = event.Data;
     plot(time,data)
     save log.mat data time
     disp(‘Background Acquisition complete’);
end