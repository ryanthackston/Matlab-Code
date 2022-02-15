function gr = drawGrid(ax, linewidth)

    hold(ax, 'on');
    gr = hggroup();
    arrayfun(@(x) plot(gr, xlim(), [x x], 'k--', 'LineWidth', linewidth), get(gca,'YTick'));
    arrayfun(@(x) plot(gr, [x,x], ylim(), 'k--', 'LineWidth', linewidth), get(gca,'XTick'));
    hold(ax,'off');
end