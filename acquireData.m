function acquireData
global data
s = daq.createSession('ni');
s.addAnalogInputChannel('cDAQ1Mod2',0,'Voltage')
s.Rate = 10000
s.DurationInSeconds = 10
lh = s.addlistener('DataAvailable',@plotData);
s.startBackground();
% Do something
while(~s.IsDone)
    
    
end

close(gcf);
plot(data); %  plot global data

function plotData(src,event)
    persistent tempData;
    global data 
    if(isempty(tempData))
         tempData = [];
    end
    plot(event.TimeStamps, event.Data)
    tempData = [tempData;event.Data];
    data = tempData;
end
end