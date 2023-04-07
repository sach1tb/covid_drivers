clearvars
close all
addpath(['boundedline', filesep, 'boundedline'])
addpath(['boundedline', filesep, 'Inpaint_nans'])
addpath(['cteUpdatedFiles', filesep])

load ukfOutput.mat  %size is 24, 24*2+1
sigmaPoints = sigmaPointAccumulutor;
datasetLength = size(sigmaPointAccumulutor,3);
P = covarianceMatrix;
meanValues = squeeze(sigmaPoints(:,1,:));
dxk = diff(sigmaPoints,1,3);


% windowSizeDays size for TE
windowSizeDays = [5*7,8*7,  12*7, 15*7, 20*7];
% # of samples from the UKF
nSigmaPoints = 21;
focussedSigmaPoints = ceil(linspace(1,49,nSigmaPoints));



for i =1:numel(windowSizeDays)
    j = 5;
    %str = sprintf('Sigma Point # %d of %d',i,nSigmaPoints);
    %disp(str);

    I = normalize(squeeze(dxk(7,j,:)))';
    Im = normalize(squeeze(dxk(8,j,:)))';
    Ih = normalize(squeeze(dxk(9,j,:)))';

    Itotdot = normalize(squeeze(dxk(7,j,:))+squeeze(dxk(8,j,:))+squeeze(dxk(9,j,:)))';

    %     xi2 = detrend(normalize(squeeze(sigmaPoints(16,j,2:end)))');
    %     xi1 = detrend(normalize(squeeze(sigmaPoints(15,j,2:end)))');
    %     phi1 = detrend(normalize(squeeze(sigmaPoints(18,j,2:end)))');
    %     phi2 = detrend(normalize(squeeze(sigmaPoints(19,j,2:end)))');
    %     sigma = detrend(normalize(squeeze(sigmaPoints(20,j,2:end)))');
    %     alpha = detrend(normalize(squeeze(sigmaPoints(17,j,2:end)))');
    %     kappa = detrend(normalize(squeeze(sigmaPoints(21,j,2:end)))');

    xi2 = detrendInAWindow(normalize(squeeze(sigmaPoints(16,j,2:end)))',windowSizeDays(i));
    xi1 = detrendInAWindow(normalize(squeeze(sigmaPoints(15,j,2:end)))',windowSizeDays(i));
    phi1 = detrendInAWindow(normalize(squeeze(sigmaPoints(18,j,2:end)))',windowSizeDays(i));
    phi2 = detrendInAWindow(normalize(squeeze(sigmaPoints(19,j,2:end)))',windowSizeDays(i));
    sigma = detrendInAWindow(normalize(squeeze(sigmaPoints(20,j,2:end)))',windowSizeDays(i));
    alpha = detrendInAWindow(normalize(squeeze(sigmaPoints(17,j,2:end)))',windowSizeDays(i));
    kappa = detrendInAWindow(normalize(squeeze(sigmaPoints(21,j,2:end)))',windowSizeDays(i));
    subplot(2,3,i)
    plot(xi2)
    hold on
    plot( xi1)
    plot(phi1)
    plot(phi2)
    plot(sigma)
    plot(alpha)
    plot(kappa)
    legend("\xi_2", "\xi_1", "\phi_1", "\phi_2", "\sigma", "\alpha", "\kappa")
    title(sprintf("Window %d",windowSizeDays(i) ))
    
end

function y = detrendInAWindow(y, windowSizeDays)
% Detrending along a window whose size is same as the TE window
% Alternatively we can decide the window size based on statioinarity
% criteria
y_smoothed = sma(y,windowSizeDays);
y = y-y_smoothed;

end