function [ ] = topoGUI(krusP,SS_chan_inter, EEG, answers)

 
rowX = SS_chan_inter(:,1);
colX = SS_chan_inter(:,2);
% 1891 Channel pairs for 62 Channel dataset
% 1891 x 1
chan_pairs = (find(tril(krusP, -1) ~=0));

% Row # & Column # for the 1891 unique values
[row,col] = find(tril(krusP, -1));

tril_krusP = tril(krusP,-1);

%Creates vector of Kruskal-Wallis values with no duplicate kruskal-wallis
%values or diagonals

ds.chanPairs = [rowX, colX];

for i=1:length(rowX)
connectStrength(i) = krusP( rowX(i), colX(i) );
end

ds.connectStrength = connectStrength';
ds.connectStrengthLimits = [min(ds.connectStrength), max(ds.connectStrength)]

loc_file = 'C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\chanLocs_OPM_August2021.ced';

figure(2);
% topoplot(ds,EEG.chanlocs, 'electrodes', 'labels')
topoplot_connect(ds,EEG.chanlocs)

colormap(lbmap(11,'RedBlue'))
colorbar

fig2 = gcf;

for i = length(SS_chan_inter):length(fig2.Children(2).Children)
    if length(fig2.Children(2).Children(i).XData) == 62
        e = fig2.Children(2).Children(i);
    else
        continue;
    end
end

e.Color = [1 1 1]

X = e.XData;
Y = e.YData;

last = length(fig2.Children(2).Children);

fig2.Children(2).Children(last).FaceColor = [1 1 1]
fig2.Children(2).Children(last-1).Color = [1 1 1]
fig2.Children(2).Children(last-2).Color = [1 1 1]
fig2.Children(2).Children(last-3).Color = [1 1 1]

labels = cell(62,1);
for i = 1:size(EEG.chanlocs,2)
% %     labels{i} = EEG.chanlocs(i).labels;
   text(X(i), Y(i)-0.02, EEG.chanlocs(i).labels);
   fig2.Children(2).Children(1).FontWeight='bold';
   fig2.Children(2).Children(1).FontSize= 10;
   fig2.Children(2).Children(1).Color = [47/256 249/256 36/256];
   fig2.Children(2).Children(1).Position = [ X(i) Y(i)-0.02 2.2]
   
% annotation('textbox', [X(i), Y(i)-.02, .03 .01], 'String',EEG.chanlocs(i).labels,'FitBoxToText','on')
end


fig2.Children(2).Title.String = ['Statistically Significant Channel Pair PLV Differences Between Rest & Move at ', answers{3}, ' Hz' ] ;
fig2.Children(2).Title.FontSize = 22;
fig2.Children(2).Title.Color = [1 1 1]

fig2.Children(1).FontSize = 12
fig2.Children(1).Label.String = 'Normalized PLV Differences'
fig2.Children(1).Label.FontSize = 18;
fig2.Children(1).Label.Color = [0.9500 0.9500 0.9500]
fig2.Children(1).Color = [0.9500 0.9500 0.9500]


fig2.Children(2).XLim = [-.5 .5]
fig2.Children(2).YLim = [-.5 .5]'
fig2.Children(2).Color = [0 0 0]
fig2.Color = [0 0 0]
'done'


end

