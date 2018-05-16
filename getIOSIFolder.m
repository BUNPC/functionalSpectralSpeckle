% getIOSIFolder - Reads raw tif files recorded by HCImage (Hamamatsu
% software) from the selected folder and calculates optical density changes
% and oxy, deoxy, and total hemoglobin changes.
% 
% Other m-files required: pathlengths.m, GetExtinctions.m
%
% Data is saved as a mat file within the selected raw images folder. 
%
% Smrithi Sunil
% email: ssunil@bu.edu
% BOAS Lab, Boston University


% Read the Data (images)
fs = input('What is the frame rate of the camera? '); % put in actual frequency (fps) here
trialTime = input('What is the length of each trial in seconds? ');
numTrials = input('How many trials were in the session? ');
wavelength = input('What wavelength was was imaged first? ');
if wavelength == 470
    lambda = [470 530 625];
elseif wavelength == 530
    lambda = [530 625 470];
elseif wavelength == 625
    lambda = [625 470 530];
else
    error('Enter a valid wavelength')
end
dataDir = uigetdir('Please select the Data folder');
dataFiles = dir([dataDir '/*.tif']);

% crop images
I=imread([dataDir '/' dataFiles(1).name]);
figure(1)
imagesc(I);
disp('Select crop region')
h = imrect(gca, [20 20 150 150]);
addNewPositionCallback(h,@(p) title(mat2str(p,3)));
fcn = makeConstrainToRectFcn('imrect',get(gca,'XLim'),get(gca,'YLim'));
setPositionConstraintFcn(h,fcn);
position = wait(h);

x(1) = int32(position(1));
x(2) = x(1)+int32(position(3));
y(1) = int32(position(2));
y(2) = y(1)+int32(position(4));
I = I(y(1):y(2),x(1):x(2));

sizeX = size(I,1);
sizeY = size(I,2);
sizeT = floor(length(dataFiles)/3);

data = uint16(zeros(3,sizeX,sizeY,sizeT));

for i = 1:1:sizeT
   waitbar(single(i)/single(sizeT));
   frame = imread([dataDir '/' dataFiles((i-1)*3+1).name]);
   data(1,:,:,i) = frame(y(1):y(2),x(1):x(2));
   frame = imread([dataDir '/' dataFiles((i-1)*3+2).name]);
   data(2,:,:,i) = frame(y(1):y(2),x(1):x(2));
   frame = imread([dataDir '/' dataFiles((i-1)*3+3).name]);
   data(3,:,:,i) = frame(y(1):y(2),x(1):x(2));
end
fclose('all');

% find optical density
numWave = size(data,1);
numFrame = size(data,4);
Ly = size(data,2);
Lx = size(data,3);

H = fspecial('gaussian',7 ,3); %For spacial filter of 5x5 pixels with 1.3 pixel standard deviation
opticalDensity = zeros([numWave numFrame Ly Lx]);
PL  = pathlengths(lambda);
for i = 1:1:numWave %iterate through each wavelength
    waitbar(single(i)/single(numWave))
    intensity0 = zeros(Ly,Lx,numTrials);
    for t = 1:numTrials
        intensity0(:,:,t) = mean(squeeze(data(i,:,:,1+(fs*trialTime/3*(t-1)):50+(fs*trialTime/3*(t-1)))),3);
    end
    intensity0 = mean(intensity0,3);
    for f = 1:numFrame
        opticalDensity_unfilt = -log(double(squeeze(data(i,:,:,f)))./intensity0)./PL(i);
        % spatially smooth optical density
        foo = conv2(opticalDensity_unfilt, H, 'same');
        opticalDensity(i,f,:,:) = reshape(foo, [1,1,Ly,Lx]);
    end
end

% Calculate changes in oxy and deoxy hemoglobin 
e = GetExtinctions(lambda);
opticalDensity = reshape(opticalDensity, [numWave,numFrame,Ly*Lx]);
Hb = zeros(numFrame,Ly*Lx,2);
for f = 1:numFrame
    waitbar(single(f)/single(numFrame))
    for i = 1:Ly*Lx
        Hb(f,i,:) = (e\[squeeze(opticalDensity(1,f,i)); squeeze(opticalDensity(2,f,i)); squeeze(opticalDensity(3,f,i))])';
    end
end

Hb = reshape(Hb, [numFrame,Ly,Lx,2]);
opticalDensity = reshape(opticalDensity, [numWave,numFrame,Ly,Lx]);
HbO = Hb(:,:,:,1);
HbR = Hb(:,:,:,2);
HbT = HbO + HbR;

opticalDensity = permute(opticalDensity, [1,3,4,2]);
HbO = permute(HbO, [2,3,1]);
HbR = permute(HbR, [2,3,1]);
HbT = permute(HbT, [2,3,1]);

save([dataDir,'\IOSI','.mat'],'fs','trialTime','numTrials','numFrame','opticalDensity','HbO','HbR','HbT','-v7.3');
