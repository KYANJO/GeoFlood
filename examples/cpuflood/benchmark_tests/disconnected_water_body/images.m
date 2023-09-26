% Set the folder containing the images
folder = './';

% Set the output video file name
outputVideo = VideoWriter('output.avi');

% Set the frame size (in pixels)
frameSize = [1050 300];

% Set the frame rate (in frames per second)
outputVideo.FrameRate = 10;

% Open the video file for writing
open(outputVideo);

% Get the list of image files in the folder
fileList = dir(fullfile(folder, '*.png'));

% Loop through each image file and add it to the video file
for i = 1:length(fileList)
    % Read the image file
    fileName = fullfile(folder, fileList(i).name);
    img = imread(fileName);

    % Resize the image to match the frame size of the video
    img = imresize(img, frameSize);

    % Add the image to the video file
    writeVideo(outputVideo, img);
end

% Close the video file
close(outputVideo);
