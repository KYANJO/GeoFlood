%  create a movie of the simulation


% make a movie
mov = VideoWriter('movie_.avi');
open(mov);
% call a global routine plotclaw2.m via the terminal which will prompt you to hit enter
plotclaw2 
% capture current figure/frame
F = getframe(gcf);
% write the current frame to the movie
writeVideo(mov,F);

close(mov)