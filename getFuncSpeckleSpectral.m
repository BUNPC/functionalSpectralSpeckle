% getFuncActIOSILSCI - Functional activation analysis for laser speckle and
% intrinsic imaging
%
% Reads the mat files created individually for spectral and speckle raw
% data. mat files are created using getLSCIFolderIMG.m for speckle data
% and getIOSIFolder.m for spectral data.
%
% Smrithi Sunil
% email: ssunil@bu.edu
% BOAS Lab, Boston University

% Updated 09/19/2018
% Updated to register the LSCI and IOSI images

%% Load corresponding speckle and spectral files

LSCI = load('Z:\users\ssunil\Stroke\CS10\Functional Activation\Baseline\LSCI_forepaw');
IOSI = load('Z:\users\ssunil\Stroke\CS10\Functional Activation\Baseline\IOSI_forepaw');

%% Average all the trials and calculate baseline 

% LSCI average and response
% assuming that flow is 1/K^2
numTrials = input('Enter the number of trials in the session ');
trialTime = input('Enter the length of each trial in seconds ');
bin = 0.1;
LSCI.time = LSCI.time - LSCI.time(1);
LSCI.time_LSC = -5:bin:trialTime-bin-5;

% trial average mean response 
binLSCI = zeros(numTrials, size(LSCI.LSCI,1), size(LSCI.LSCI,2), trialTime/bin);
for j = 1:numTrials
    idx = find(LSCI.time>=trialTime*(j-1) & LSCI.time<=trialTime*j);
    LSCI.LSCIs.reshape(:,:,:) = LSCI.LSCI(:,:,idx);
    LSCI.LSCIs.respTime = LSCI.time(idx)-trialTime*(j-1);
    for k = 1:trialTime/bin
        idx2 = find(LSCI.LSCIs.respTime>bin*(k-1) & LSCI.LSCIs.respTime<=bin*k);
        binLSCI(j,:,:,k) = nanmean(squeeze(LSCI.LSCIs.reshape(:,:,idx2)),3);
    end
    LSCI.LSCIs = [];
end
LSCI.LSCIavg = nanmean(binLSCI(:,:,:,:),1);
LSCI.LSCIavg = permute(LSCI.LSCIavg, [2,3,4,1]);
LSCI.flowIndex = 1./(LSCI.LSCIavg.*LSCI.LSCIavg);

bStart = 1;
bEnd = find(LSCI.time_LSC < 0,1,'last');
LSCI.baseflowIndex = LSCI.flowIndex(:,:,bStart:bEnd); 
LSCI.mBaseflowIndex = squeeze(nanmean(LSCI.baseflowIndex,3));


% IOSI average and response 
IOSI.time_IOSfull = 0.1:0.1:IOSI.trialTime*IOSI.numTrials;
IOSI.HbTreshape = zeros([IOSI.numTrials,size(IOSI.HbT,1),size(IOSI.HbT,2),IOSI.numFrame/IOSI.numTrials]);
for t = 1:IOSI.numTrials
    IOSI.HbTreshape(t,:,:,:) = IOSI.HbT(:,:,(IOSI.numFrame/IOSI.numTrials*t)-(IOSI.numFrame/IOSI.numTrials-1):IOSI.numFrame/IOSI.numTrials*t);
end
IOSI.HbTavg = zeros([size(IOSI.HbT,1),size(IOSI.HbT,2),IOSI.numFrame/IOSI.numTrials]);
for i = 1:IOSI.numFrame/IOSI.numTrials 
    IOSI.HbTavg(:,:,i) = mean(squeeze(IOSI.HbTreshape(:,:,:,i)),1);  
end
IOSI.time_IOS = -4.9:0.1:IOSI.trialTime-5;
bStart = 1;
bEnd = find(IOSI.time_IOS < 0,1,'last');
IOSI.HbTavg = imrotate(IOSI.HbTavg,90);
IOSI.basevolIndex = IOSI.HbTavg(:,:,bStart:bEnd); 
IOSI.mBasevolIndex = squeeze(mean(IOSI.basevolIndex,3));

% Take the average of HbR
IOSI.HbRreshape = zeros([IOSI.numTrials,size(IOSI.HbR,1),size(IOSI.HbR,2),IOSI.numFrame/IOSI.numTrials]);
for t = 1:IOSI.numTrials
    IOSI.HbRreshape(t,:,:,:) = IOSI.HbR(:,:,(IOSI.numFrame/IOSI.numTrials*t)-(IOSI.numFrame/IOSI.numTrials-1):IOSI.numFrame/IOSI.numTrials*t);
end
IOSI.HbRavg = zeros([size(IOSI.HbR,1),size(IOSI.HbR,2),IOSI.numFrame/IOSI.numTrials]);
for i = 1:IOSI.numFrame/IOSI.numTrials 
    IOSI.HbRavg(:,:,i) = mean(squeeze(IOSI.HbRreshape(:,:,:,i)),1);  
end
IOSI.HbRavg = imrotate(IOSI.HbRavg,90);

%% Register LSCI and IOSI baseline images

moving=double(LSCI.mBaseflowIndex);
fixed=double(IOSI.mBasevolIndex);
figure(99)
colormap jet
subplot(1,2,1) 
imagesc(moving)
caxis([prctile(moving(:),5), prctile(moving(:),95)]); 
axis image
subplot(1,2,2)
imagesc(fixed)
caxis([prctile(fixed(:),5), prctile(fixed(:),95)]); 
axis image
[x2, y2] = ginput(10);
[D,Z,TRANSFORM] = procrustes([y2(2:2:10) x2(2:2:10)], [y2(1:2:10) x2(1:2:10)]);

IOSI.regHbTavg = zeros(size(LSCI.flowIndex));
IOSI.regHbRavg = zeros(size(LSCI.flowIndex));
for z = 1:size(LSCI.flowIndex,3)
    for u = 1:size(LSCI.mBaseflowIndex,1)
        for v = 1:size(LSCI.mBaseflowIndex,2)
            idx = round(TRANSFORM.b * [u v] * TRANSFORM.T + TRANSFORM.c(1,:));
            if idx(1) < 1 || idx(1) > size(IOSI.HbTavg,1) || idx(2) < 1 || idx(2) > size(IOSI.HbTavg,2)
                IOSI.regHbTavg(u,v,z) = 0;
                IOSI.regHbRavg(u,v,z) = 0;
            else
                IOSI.regHbTavg(u,v,z) = IOSI.HbTavg(idx(1),idx(2),z);
                IOSI.regHbRavg(u,v,z) = IOSI.HbRavg(idx(1),idx(2),z);
            end
        end
    end
end


img = zeros(size(LSCI.mBaseflowIndex));
for u = 1:size(LSCI.mBaseflowIndex,1)
    for v = 1:size(LSCI.mBaseflowIndex,2)
        idx = round(TRANSFORM.b * [u v] * TRANSFORM.T + TRANSFORM.c(1,:));
        if idx(1) < 1 || idx(1) > size(IOSI.HbTavg,1) || idx(2) < 1 || idx(2) > size(IOSI.HbTavg,2)
            img(u,v) = 0;
        else
            img(u,v) = IOSI.mBasevolIndex(idx(1),idx(2));
        end
    end
end

% crop images
figure(100)
colormap jet
subplot(1,2,1)
imagesc(moving) 
caxis([prctile(moving(:),5), prctile(moving(:),95)])
axis image
subplot(1,2,2)
imagesc(img)
caxis([prctile(img(:),5), prctile(img(:),95)])
axis image
title('Select crop region')
h = imrect(gca, [20 20 150 150]);
addNewPositionCallback(h,@(p) title(mat2str(p,3)));
fcn = makeConstrainToRectFcn('imrect',get(gca,'XLim'),get(gca,'YLim'));
setPositionConstraintFcn(h,fcn);
position = wait(h);
x(1) = int32(position(1));
x(2) = x(1)+int32(position(3));
y(1) = int32(position(2));
y(2) = y(1)+int32(position(4));
LSCI.flowIndex = LSCI.flowIndex(y(1):y(2),x(1):x(2),:);
IOSI.regHbTavg = IOSI.regHbTavg(y(1):y(2),x(1):x(2),:);
IOSI.regHbRavg = IOSI.regHbRavg(y(1):y(2),x(1):x(2),:);

clear LSCI.mBaseflowIndex LSCI.baseflowIndex
clear IOSI.mBasevolIndex IOSI.basevolIndex

%% Plot baseline and response data

bStart = 1;
bEnd = find(LSCI.time_LSC < 0,1,'last');
rStart = find(LSCI.time_LSC > 0,1,'first');
rEnd = find(LSCI.time_LSC > 5,1,'first');
LSCI.baseflowIndex = LSCI.flowIndex(:,:,bStart:bEnd); 
LSCI.mBaseflowIndex = squeeze(nanmean(LSCI.baseflowIndex,3));
LSCI.respflowIndex = LSCI.flowIndex(:,:,rStart:rEnd)./LSCI.mBaseflowIndex;
LSCI.mRespflowIndex = squeeze(nanmean(LSCI.respflowIndex,3));

HbT_base = 100*10^-6;  % assumed baseline from literature
bStart = 1;
bEnd = find(IOSI.time_IOS < 0,1,'last');
rStart = find(IOSI.time_IOS > 0,1,'first');
rEnd = find(IOSI.time_IOS > 5,1,'first');
IOSI.basevolIndex = IOSI.regHbTavg(:,:,bStart:bEnd); 
IOSI.mBasevolIndex = squeeze(mean(IOSI.basevolIndex,3));
IOSI.respvolIndex = 1 + IOSI.regHbTavg(:,:,rStart:rEnd)./HbT_base;
IOSI.mRespvolIndex = squeeze(mean(IOSI.respvolIndex,3));

% show flow map and select rectangular ROI
figure(1);
subplot(2,2,1)
imagesc(LSCI.mBaseflowIndex)
caxis([prctile(LSCI.mBaseflowIndex(:),5),prctile(LSCI.mBaseflowIndex(:),95)]);
title('Baseline flow','FontSize',16);
colormap jet
colorbar
axis image
% show response map
figure(1)
subplot(2,2,2)
imagesc(LSCI.mRespflowIndex)
caxis([prctile(LSCI.mRespflowIndex(:),5),prctile(LSCI.mRespflowIndex(:),95)]);
% caxis([90 120])
title('Response map, select ROI');
colorbar
axis image

% show vol map and select rectangular ROI
figure(2)
subplot(2,2,1)
imagesc(IOSI.mBasevolIndex)
caxis([prctile(IOSI.mBasevolIndex(:),5),prctile(IOSI.mBasevolIndex(:),95)]);
title('Baseline HbT','FontSize',16);
colormap jet
colorbar
axis image
% show response map
figure(2)
subplot(2,2,2)
imagesc(IOSI.mRespvolIndex)
caxis([prctile(IOSI.mRespvolIndex(:),5),prctile(IOSI.mRespvolIndex(:),95)]);
% caxis([-1 2.5])
title('Response map, select ROI');
colorbar
axis image

% select and plot ROIs
for i = 1:3
    figure(1)
    subplot(2,2,2)
    h = imrect;
    pos = wait(h);
    pos = round(pos);
    delete(h);
    hold on
    plot([pos(1),pos(1),pos(1)+pos(3),pos(1)+pos(3),pos(1)],[pos(2),pos(2)+pos(4),pos(2)+pos(4),pos(2),pos(2)],'LineWidth',3);
    hold off
    LSCI.flowIndexROI = LSCI.flowIndex(pos(2):1:pos(2)+pos(4),pos(1):1:pos(1)+pos(3),:);
    LSCI.respROI = LSCI.flowIndexROI./LSCI.baseflowIndex(pos(2):1:pos(2)+pos(4),pos(1):1:pos(1)+pos(3));
    LSCI.mRespRoi_LSC = squeeze(mean(mean(LSCI.respROI,1),2));
    figure(1)
    subplot(2,2,[3 4])
    hold on
    plot(LSCI.time_LSC,LSCI.mRespRoi_LSC,'LineWidth',3)
    hold off
    xlabel('Time, s','FontSize',16)
    ylabel('Fractional change','FontSize',16)
    title('ROI flow','FontSize',16)
    legend({'ROI 1','ROI 2','ROI 3'},'FontSize',12)
    
    figure(2)
    subplot(2,2,2)
    hold on
    plot([pos(1),pos(1),pos(1)+pos(3),pos(1)+pos(3),pos(1)],[pos(2),pos(2)+pos(4),pos(2)+pos(4),pos(2),pos(2)],'LineWidth',3);
    hold off
    IOSI.volIndexROI = IOSI.regHbTavg(pos(2):1:pos(2)+pos(4),pos(1):1:pos(1)+pos(3),:);
    IOSI.respROI = 1 + IOSI.volIndexROI./HbT_base;
    IOSI.mRespRoi_IOS(i).roi = squeeze(mean(mean(IOSI.respROI,1),2));
    figure(2)
    subplot(2,2,[3 4])
    hold on
    plot(IOSI.time_IOS',IOSI.mRespRoi_IOS(i).roi,'LineWidth',3)
    hold off
    xlabel('Time, s','FontSize',16)
    ylabel('Fractional change','FontSize',16)
    title('ROI HbT','FontSize',16)
    legend({'ROI 1','ROI 2','ROI 3'},'FontSize',12)
end
figure(1)
subplot(2,2,2)
title('Response map','FontSize',16);
figure(2)
subplot(2,2,2)
title('Response map','FontSize',16);


