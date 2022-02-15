function plotData(src,event)
    global tempData;
    global data
    global time
    if(isempty(tempData))
         tempData = [];
     end
     plot(event.TimeStamps, event.Data)
     tempData = [tempData;event.Data];
     data = tempData;
     time = [time; event.TimeStamps];
end