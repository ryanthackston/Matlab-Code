fig = figure(1);

handle_patch = patch([0 1 1 0], [0 0 1 1],'b');

period = 0.1;
color  = 1.0;
seconds = 10;

tic;
loop_counter = 0;
loop_max     = (1/period)*seconds;

while (loop_counter < loop_max)
    t = toc;
    if (t > period/2)
        tic
        color = color * -1;
        loop_counter = loop_counter + 0.5;
        if (color == 1)
           set(handle_patch,'FaceColor','r');
        else
           set(handle_patch,'FaceColor','b');
        end
        drawnow;
    end
end
