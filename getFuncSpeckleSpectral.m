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

%% Load corresponding speckle and spectral files

LSCI = load('D:\FlowVolume_04252018\04252018_Charlie_LSCI_2sec\LSCI');
IOSI = load('D:\FlowVolume_04252018\04252018_Charlie_IOSI_2sec\IOSI');

%% Convert contrast to flow for laser speckle data

% assuming that flow is 1/K^2
numTrials = input('Enter the number of trials in the session ');
trialTime = input('Enter the length of each trial in seconds ');
bin = 0.25;

LSCI.time = LSCI.time - LSCI.time(1);
figure
subplot(3,2,[1 2])
plot(LSCI.time',squeeze(mean(mean(LSCI.LSCI,1),2)),'LineWidth',2);
xlabel('Time,s','FontSize',16)
ylabel('Contrast','FontSize',16)
title('Speckle contrast - full time series','FontSize',16)
LSCI.time_LSC = 0:bin:trialTime-bin;

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

subplot(3,2,4)
plot(LSCI.time_LSC',squeeze(mean(mean(LSCI.LSCIavg,1),2)),'LineWidth',3);
xlabel('Time,s','FontSize',16)
ylabel('Contrast','FontSize',16)
title('Select baseline time')
h = imrect;
pos = wait(h);
[~,bStart] = min(abs(LSCI.time_LSC-pos(1)));
[~,bEnd] = min(abs(LSCI.time_LSC-pos(1)-pos(3)));
if bStart < 1
    bStart = 1;
end
if bEnd > size(LSCI.flowIndex,3)
    bEnd = size(LSCI.flowIndex,3);
end

title('Select response time')
h = imrect;
pos = wait(h);
[~,rStart] = min(abs(LSCI.time_LSC-pos(1)));
[~,rEnd] = min(abs(LSCI.time_LSC-pos(1)-pos(3)));
if rStart < 1
    rStart = 1;
end
if rEnd > size(LSCI.flowIndex,3)
    rEnd = size(LSCI.flowIndex,3);
end
title('Selected baseline and response times','FontSize',16)

LSCI.baseflowIndex = LSCI.flowIndex(:,:,bStart:bEnd); 
LSCI.mBaseflowIndex = squeeze(mean(LSCI.baseflowIndex,3));

LSCI.respflowIndex = LSCI.flowIndex(:,:,rStart:rEnd)./LSCI.mBaseflowIndex;
LSCI.mRespflowIndex = squeeze(mean(LSCI.respflowIndex,3));

% show flow map and select rectangular ROI
subplot(3,2,3)
imagesc(LSCI.mBaseflowIndex)
caxis([prctile(LSCI.mBaseflowIndex(:),5),prctile(LSCI.mBaseflowIndex(:),95)]);
title('Baseline flow','FontSize',16);
colorbar
axis image

% show response map
subplot(3,2,5)
imagesc(LSCI.mRespflowIndex)
caxis([prctile(LSCI.mRespflowIndex(:),5),prctile(LSCI.mRespflowIndex(:),95)]);
title('Response map, select ROI');
colorbar
axis image

%select and plot ROIs
for i = 1:3
    subplot(3,2,5)
    h = imrect;
    pos = wait(h);
    delete(h);
    hold on
    plot([pos(1),pos(1),pos(1)+pos(3),pos(1)+pos(3),pos(1)],[pos(2),pos(2)+pos(4),pos(2)+pos(4),pos(2),pos(2)],'LineWidth',3);
    hold off

    pos = round(pos);
    LSCI.flowIndexROI = LSCI.flowIndex(pos(2):1:pos(2)+pos(4),pos(1):1:pos(1)+pos(3),:);
    LSCI.respROI = LSCI.flowIndexROI./LSCI.baseflowIndex(pos(2):1:pos(2)+pos(4),pos(1):1:pos(1)+pos(3));
    LSCI.mRespRoi_LSC = squeeze(mean(mean(LSCI.respROI,1),2));

    subplot(3,2,6)
    hold on
    plot(LSCI.time_LSC,LSCI.mRespRoi_LSC.*100,'LineWidth',3)
    hold off
    xlabel('Time, s','FontSize',16)
    ylabel('Response, %','FontSize',16)
    title('ROI normalized flow (trial averaged)','FontSize',16)
    legend({'ROI 1','ROI 2','ROI 3'},'FontSize',12)
end
subplot(3,2,5)
title('Response map','FontSize',16);

%% Convert HbT into activation map for spectral data

% IOSI post process
% Average all trials in the session

figure
subplot(3,2,[1 2])
IOSI.time_IOSfull = 0.1:0.1:IOSI.trialTime*20;
plot(IOSI.time_IOSfull',squeeze(mean(mean(IOSI.HbT,1),2)),'LineWidth',2);
xlabel('Time,s','FontSize',16)
ylabel('HbT change','FontSize',16)
title('HbT - full time series','FontSize',16)

IOSI.HbTreshape = zeros([IOSI.numTrials,size(IOSI.HbT,1),size(IOSI.HbT,2),IOSI.numFrame/IOSI.numTrials]);
for t = 1:IOSI.numTrials
    IOSI.HbTreshape(t,:,:,:) = IOSI.HbT(:,:,(IOSI.numFrame/IOSI.numTrials*t)-(IOSI.numFrame/IOSI.numTrials-1):IOSI.numFrame/IOSI.numTrials*t);
end
IOSI.HbTavg = zeros([size(IOSI.HbT,1),size(IOSI.HbT,2),IOSI.numFrame/IOSI.numTrials]);
for i = 1:IOSI.numFrame/IOSI.numTrials 
    IOSI.HbTavg(:,:,i) = mean(squeeze(IOSI.HbTreshape(:,:,:,i)),1);  
end

IOSI.time_IOS = 0.1:0.1:IOSI.trialTime;

subplot(3,2,4)
plot(IOSI.time_IOS',squeeze(mean(mean(IOSI.HbTavg,1),2)),'LineWidth',3);
xlabel('Time,s','FontSize',16)
ylabel('HbT change','FontSize',16)
title('Select baseline time')
h = imrect;
pos = wait(h);
[~,bStart] = min(abs(IOSI.time_IOS-pos(1)));
[~,bEnd] = min(abs(IOSI.time_IOS-pos(1)-pos(3)));
if bStart < 1
    bStart = 1;
end
if bEnd > size(IOSI.HbTavg,3)
    bEnd = size(IOSI.HbTavg,3);
end

title('Select response time')
h = imrect;
pos = wait(h);
[~,rStart] = min(abs(IOSI.time_IOS-pos(1)));
[~,rEnd] = min(abs(IOSI.time_IOS-pos(1)-pos(3)));
if rStart < 1
    rStart = 1;
end
if rEnd > size(IOSI.HbTavg,3)
    rEnd = size(IOSI.HbTavg,3);
end
title('Selected baseline and response times','FontSize',16)

IOSI.basevolIndex = IOSI.HbTavg(:,:,bStart:bEnd); 
IOSI.mBasevolIndex = squeeze(mean(IOSI.basevolIndex,3));
IOSI.respvolIndex = IOSI.HbTavg(:,:,rStart:rEnd)-IOSI.mBasevolIndex;
IOSI.mRespvolIndex = squeeze(mean(IOSI.respvolIndex,3));

% show vol map and select rectangular ROI
subplot(3,2,3)
imagesc(IOSI.mBasevolIndex)
caxis([prctile(IOSI.mBasevolIndex(:),5),prctile(IOSI.mBasevolIndex(:),95)]);
title('Baseline HbT','FontSize',16);
colorbar
axis image

% show response map
subplot(3,2,5)
imagesc(IOSI.mRespvolIndex)
caxis([prctile(IOSI.mRespvolIndex(:),5),prctile(IOSI.mRespvolIndex(:),95)]);
title('Response map, select ROI');
colorbar
axis image

%select and plot ROIs
for i = 1:3
    subplot(3,2,5)
    h = imrect;
    pos = wait(h);
    delete(h);
    hold on
    plot([pos(1),pos(1),pos(1)+pos(3),pos(1)+pos(3),pos(1)],[pos(2),pos(2)+pos(4),pos(2)+pos(4),pos(2),pos(2)],'LineWidth',3);
    hold off

    pos = round(pos);
    IOSI.volIndexROI = IOSI.HbTavg(pos(2):1:pos(2)+pos(4),pos(1):1:pos(1)+pos(3),:);
    IOSI.respROI = IOSI.volIndexROI-IOSI.basevolIndex(pos(2):1:pos(2)+pos(4),pos(1):1:pos(1)+pos(3));
    IOSI.mRespRoi_IOS = squeeze(mean(mean(IOSI.respROI,1),2));

    subplot(3,2,6)
    hold on
    plot(IOSI.time_IOS',IOSI.mRespRoi_IOS,'LineWidth',3)
    hold off
    xlabel('Time, s','FontSize',16)
    ylabel('Response','FontSize',16)
    title('ROI HbT change (trial averaged)','FontSize',16)
    legend({'ROI 1','ROI 2','ROI 3'},'FontSize',12)
end
subplot(3,2,5)
title('Response map','FontSize',16);
