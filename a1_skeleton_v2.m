% Read the 3D tif images
clear;
imgname = "g2_large_s.tif";
tiffInfo = imfinfo(imgname);
numSlices = numel(tiffInfo);

imageData = zeros(tiffInfo(1).Height, tiffInfo(1).Width, numSlices);

for k = 1:numSlices
    imageData(:,:,k) = imread(imgname, k);
end
a = size(imageData);
fcrop = 0.95;
lim = floor([a' * 0.05, a' * 0.95]);
imageData = imageData(lim(1, 1): lim(1, 2), lim(2, 1): lim(2, 2), lim(3, 1): lim(3, 2));

% Set resolution
resolutionX = 5.01;
resolutionY = 5.01;
resolutionZ = 5.02;

[blackX, blackY, blackZ] = ind2sub(size(imageData), find(imageData == 0));

realX = blackX * resolutionX;
realY = blackY * resolutionY;
realZ = blackZ * resolutionZ;

% Save the coordinates of vertices as skeleton_coordinates.mat
coordinates = [realY, realX, realZ];
save('skeleton_coordinates.mat', "coordinates");

% Visualization
figure;
scatter3(realY, realX, realZ, 1);
axis equal;

