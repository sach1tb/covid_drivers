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
windowSizeDays = 12*7;
% # of samples from the UKF
nSigmaPoints = 24;

%CE_Stot_Itot_condE = zeros(nSigmaPoints,datasetLength-windowSizeDays); % (CE_X_Y_Z) conditional TE X->Y conditioned on Z
%CE_Stot_Itot_condEm = zeros(nSigmaPoints,datasetLength-windowSizeDays);
%CE_Stot_Itot_condEh = zeros(nSigmaPoints,datasetLength-windowSizeDays);
%TE_Stot_Itot= zeros(nSigmaPoints,datasetLength-windowSizeDays); % (TE_X_Y) TE X->Y


TE_phi1_Itot = zeros(nSigmaPoints,datasetLength-windowSizeDays);
TE_Itot_phi1 = zeros(nSigmaPoints,datasetLength-windowSizeDays);
NetTE_phi1_Itot = zeros(nSigmaPoints,datasetLength-windowSizeDays);

TE_xi2_Itot = zeros(nSigmaPoints,datasetLength-windowSizeDays);
TE_Itot_xi2 = zeros(nSigmaPoints,datasetLength-windowSizeDays);
NetTE_xi2_Itot = zeros(nSigmaPoints,datasetLength-windowSizeDays);


TE_phi2_Itot = zeros(nSigmaPoints,datasetLength-windowSizeDays);
TE_Itot_phi2 = zeros(nSigmaPoints,datasetLength-windowSizeDays);
NetTE_phi2_Itot = zeros(nSigmaPoints,datasetLength-windowSizeDays);

TE_xi1_Itot = zeros(nSigmaPoints,datasetLength-windowSizeDays);
TE_Itot_xi1 = zeros(nSigmaPoints,datasetLength-windowSizeDays);
NetTE_xi1_Itot = zeros(nSigmaPoints,datasetLength-windowSizeDays);

TE_sigma_Itot = zeros(nSigmaPoints,datasetLength-windowSizeDays);
TE_Itot_sigma = zeros(nSigmaPoints,datasetLength-windowSizeDays);
NetTE_sigma_Itot = zeros(nSigmaPoints,datasetLength-windowSizeDays);

TE_alpha_Itot = zeros(nSigmaPoints,datasetLength-windowSizeDays);
TE_Itot_alpha = zeros(nSigmaPoints,datasetLength-windowSizeDays);
NetTE_alpha_Itot = zeros(nSigmaPoints,datasetLength-windowSizeDays);

TE_kappa_Itot = zeros(nSigmaPoints,datasetLength-windowSizeDays);
TE_Itot_kappa = zeros(nSigmaPoints,datasetLength-windowSizeDays);
NetTE_kappa_Itot = zeros(nSigmaPoints,datasetLength-windowSizeDays);

%%
S = normalize(squeeze(dxk(1,:,:)))';
Sm = normalize(squeeze(dxk(2,:,:)))';
Sh= normalize(squeeze(dxk(3,:,:)))';

E = normalize(squeeze(dxk(4,:,:)))';
Em = normalize(squeeze(dxk(5,:,:)))';
Eh = normalize(squeeze(dxk(6,:,:)))';

I = normalize(squeeze(dxk(7,:,:)))';
Im = normalize(squeeze(dxk(8,:,:)))';
Ih = normalize(squeeze(dxk(9,:,:)))';

R = normalize(squeeze(dxk(10,:,:)))';
D = normalize(squeeze(dxk(11,:,:)))';
U = normalize(squeeze(dxk(12,:,:)))';
V = normalize(squeeze(dxk(13,:,:)))';
Stot = normalize((squeeze(dxk(1,:,:))+squeeze(dxk(2,:,:))+squeeze(dxk(3,:,:)))');
Itot = normalize((squeeze(dxk(7,:,:))+squeeze(dxk(8,:,:))+squeeze(dxk(9,:,:)))');
xi2 = normalize(squeeze(sigmaPoints(16,:,2:end))'.* squeeze(sigmaPoints(1,:,2:end))');
xi1 = normalize(squeeze(sigmaPoints(15,:,2:end)))'.* normalize(squeeze(sigmaPoints(3,:,2:end))');
phi1 = normalize(squeeze(sigmaPoints(18,:,2:end))'.*squeeze(sigmaPoints(1,:,2:end))');
phi2 = normalize(squeeze(sigmaPoints(19,:,2:end)))'.*normalize(squeeze(sigmaPoints(2,:,2:end))');
sigma = normalize(squeeze(sigmaPoints(20,:,2:end))'.*squeeze(sigmaPoints(12,:,2:end))');
alpha = normalize(squeeze(sigmaPoints(17,:,2:end))'.*(((squeeze(sigmaPoints(1,:,2:end)))+squeeze(sigmaPoints(2,:,2:end)))+normalize(squeeze(sigmaPoints(3,:,2:end))))');
kappa = (normalize(squeeze(sigmaPoints(21,:,2:end)))');

NetTE_phi1_Itot = zeros(nSigmaPoints, datasetLength-windowSizeDays);
for i =1:nSigmaPoints
    parfor k =1:datasetLength-windowSizeDays
        b = Itot(k:k+windowSizeDays-1,i)';
        NetTE_phi1_Itot(i,k) =   calcNetTE_XtoY(phi1(k:k+windowSizeDays-1,i)',b,windowSizeDays) ;
        NetTE_xi2_Itot(i,k) =   calcNetTE_XtoY(xi2(k:k+windowSizeDays-1,i)',b,windowSizeDays) ;
        NetTE_phi2_Itot(i,k) =   calcNetTE_XtoY(phi2(k:k+windowSizeDays-1,i)',b,windowSizeDays) ;
        NetTE_xi1_Itot(i,k) =   calcNetTE_XtoY(xi1(k:k+windowSizeDays-1,i)',b,windowSizeDays) ;
        NetTE_sigma_Itot(i,k) =   calcNetTE_XtoY(sigma(k:k+windowSizeDays-1,i)',b,windowSizeDays) ;
        NetTE_alpha_Itot(i,k) =   calcNetTE_XtoY(alpha(k:k+windowSizeDays-1,i)',b,windowSizeDays) ;
        NetTE_kappa_Itot(i,k) =   calcNetTE_XtoY(kappa(k:k+windowSizeDays-1,i)',b,windowSizeDays) ;
        
        fprintf('.');
        if rem(k,windowSizeDays) == 0
            fprintf('\n');
        end
    end
end

function y = calcNetTE_XtoY(X,Y,windowSizeDays) %net TE X to Y
addpath(['boundedline', filesep, 'boundedline'])
addpath(['boundedline', filesep, 'Inpaint_nans'])
addpath(['cteUpdatedFiles', filesep])

Hy = ent(Y, ceil(sqrt(windowSizeDays)), [-1 1], 'x');
Hx = ent(X, ceil(sqrt(windowSizeDays)), [-1 1], 'x');
TE_X_Y = ete_hist(Y,X,1,ceil(sqrt(windowSizeDays)),[-1 1]);
TE_X_Y = TE_X_Y/sqrt(Hx*Hy);
TE_Y_X = ete_hist(X,Y,1,ceil(sqrt(windowSizeDays)),[-1 1]);
TE_Y_X = TE_Y_X/sqrt(Hx*Hy);
NetTE_X_Y =  TE_X_Y - TE_Y_X;

y = NetTE_X_Y;

end


% function y = doTE()
% addpath(['boundedline', filesep, 'boundedline'])
% addpath(['boundedline', filesep, 'Inpaint_nans'])
% addpath(['cteUpdatedFiles', filesep])
% load ukfOutput.mat  %size is 24, 24*2+1
% sigmaPoints = sigmaPointAccumulutor;
% datasetLength = size(sigmaPointAccumulutor,3);
% P = covarianceMatrix;
% meanValues = squeeze(sigmaPoints(:,1,:));
% dxk = diff(sigmaPoints,1,3);
% 
% 
% % windowSizeDays size for TE
% windowSizeDays = 8*7;
% str = sprintf('Sigma Point # %d',i);
% disp(str);
% %normalizing over the whole dataset length
% % S = normalize(squeeze(dxk(1,i,:)))';
% % Sm = normalize(squeeze(dxk(2,i,:)))';
% % Sh= normalize(squeeze(dxk(3,i,:)))';
% %
% % E = normalize(squeeze(dxk(4,i,:)))';
% % Em = normalize(squeeze(dxk(5,i,:)))';
% % Eh = normalize(squeeze(dxk(6,i,:)))';
% %
% % I = normalize(squeeze(dxk(7,i,:)))';
% % Im = normalize(squeeze(dxk(8,i,:)))';
% % Ih = normalize(squeeze(dxk(9,i,:)))';
% %
% % R = normalize(squeeze(dxk(10,i,:)))';
% % D = normalize(squeeze(dxk(11,i,:)))';
% % U = normalize(squeeze(dxk(12,i,:)))';
% % V = normalize(squeeze(dxk(13,i,:)))';
% % Stot = normalize((squeeze(dxk(1,i,:))+squeeze(dxk(2,i,:))+squeeze(dxk(3,i,:)))');
% % Itot = normalize((squeeze(dxk(7,i,:))+squeeze(dxk(8,i,:))+squeeze(dxk(9,i,:)))');
% % xi2 = normalize(squeeze(sigmaPoints(16,i,2:end))'.* squeeze(sigmaPoints(1,i,2:end))');
% % xi1 = normalize(squeeze(sigmaPoints(15,i,2:end)))'.* normalize(squeeze(sigmaPoints(3,i,2:end))');
% % phi1 = normalize(squeeze(sigmaPoints(18,i,2:end))'.*squeeze(sigmaPoints(1,i,2:end))');
% % phi2 = normalize(squeeze(sigmaPoints(19,i,2:end)))'.*normalize(squeeze(sigmaPoints(2,i,2:end))');
% % sigma = normalize(squeeze(sigmaPoints(20,i,2:end))'.*squeeze(sigmaPoints(12,i,2:end))');
% % alpha = normalize(squeeze(sigmaPoints(17,i,2:end))'.*(((squeeze(sigmaPoints(1,i,2:end)))+squeeze(sigmaPoints(2,i,2:end)))+normalize(squeeze(sigmaPoints(3,i,2:end))))');
% % kappa = (normalize(squeeze(sigmaPoints(21,i,2:end)))');
% 
% %raw values, figure 2 add vaccination rate, figure 3 add loss of immunity
% %after vaccination.
% 
% for k = 1:1:datasetLength-windowSizeDays
% 
%     %TE_Stot_Itot(i,k) = ete_hist(Itot(k:k+windowSizeDays-1),Stot(k:k+windowSizeDays-1),1,ceil(sqrt(windowSizeDays)),[-1 1]);
% 
%     Hy = ent(Itot(k:k+windowSizeDays-1), ceil(sqrt(windowSizeDays)), [-1 1], 'x');
% 
%     Hx = ent(phi1(k:k+windowSizeDays-1), ceil(sqrt(windowSizeDays)), [-1 1], 'x');
% 
% 
% 
%     TE_phi1_Itot(i,k) = ete_hist(Itot(k:k+windowSizeDays-1),phi1(k:k+windowSizeDays-1),1,ceil(sqrt(windowSizeDays)),[-1 1]);
%     TE_phi1_Itot(i,k) = TE_phi1_Itot(i,k)/sqrt(Hx*Hy);
%     TE_Itot_phi1(i,k) = ete_hist(phi1(k:k+windowSizeDays-1),Itot(k:k+windowSizeDays-1),1,ceil(sqrt(windowSizeDays)),[-1 1]);
%     TE_Itot_phi1(i,k) = TE_phi1_Itot(i,k)/sqrt(Hx*Hy);
%     NetTE_phi1_Itot(i,k) =  TE_phi1_Itot(i,k) - TE_Itot_phi1(i,k);
% 
% 
%     Hx = ent(xi2(k:k+windowSizeDays-1), ceil(sqrt(windowSizeDays)), [-1 1], 'x');
%     TE_xi2_Itot(i,k) = ete_hist(Itot(k:k+windowSizeDays-1),xi2(k:k+windowSizeDays-1),1,ceil(sqrt(windowSizeDays)),[-1 1]);
%     TE_xi2_Itot(i,k) = TE_xi2_Itot(i,k)/sqrt(Hx*Hy);
%     TE_Itot_xi2(i,k) = ete_hist(xi2(k:k+windowSizeDays-1),Itot(k:k+windowSizeDays-1),1,ceil(sqrt(windowSizeDays)),[-1 1]);
%     TE_Itot_xi2(i,k) = TE_Itot_xi2(i,k)/sqrt(Hx*Hy);
%     NetTE_xi2_Itot(i,k) = TE_xi2_Itot(i,k)-TE_Itot_xi2(i,k) ;
% 
% 
%     Hx = ent(phi2(k:k+windowSizeDays-1), ceil(sqrt(windowSizeDays)), [-1 1], 'x');
%     TE_phi2_Itot(i,k) = ete_hist(Itot(k:k+windowSizeDays-1),phi2(k:k+windowSizeDays-1),1,ceil(sqrt(windowSizeDays)),[-1 1]);
%     TE_phi2_Itot(i,k) = TE_phi2_Itot(i,k)/sqrt(Hx*Hy);
%     TE_Itot_phi2(i,k) = ete_hist(phi2(k:k+windowSizeDays-1),Itot(k:k+windowSizeDays-1),1,ceil(sqrt(windowSizeDays)),[-1 1]);
%     TE_Itot_phi2(i,k) = TE_Itot_phi2(i,k)/sqrt(Hx*Hy);
%     NetTE_phi2_Itot(i,k) = TE_phi2_Itot(i,k) - TE_Itot_phi2(i,k);
% 
% 
%     Hx = ent(xi1(k:k+windowSizeDays-1), ceil(sqrt(windowSizeDays)), [-1 1], 'x');
%     TE_xi1_Itot(i,k) = ete_hist(Itot(k:k+windowSizeDays-1),xi1(k:k+windowSizeDays-1),1,ceil(sqrt(windowSizeDays)),[-1 1]);
%     TE_xi1_Itot(i,k) = TE_xi1_Itot(i,k)/sqrt(Hx*Hy);
%     TE_Itot_xi1(i,k) = ete_hist(xi1(k:k+windowSizeDays-1),Itot(k:k+windowSizeDays-1),1,ceil(sqrt(windowSizeDays)),[-1 1]);
%     TE_Itot_xi1(i,k) = TE_Itot_xi1(i,k)/sqrt(Hx*Hy);
%     NetTE_xi1_Itot(i,k) = TE_xi1_Itot(i,k) - TE_Itot_xi1(i,k);
% 
% 
%     Hx = ent(sigma(k:k+windowSizeDays-1), ceil(sqrt(windowSizeDays)), [-1 1], 'x');
%     TE_sigma_Itot(i,k) = ete_hist(Itot(k:k+windowSizeDays-1),sigma(k:k+windowSizeDays-1),1,ceil(sqrt(windowSizeDays)),[-1 1]);
%     TE_sigma_Itot(i,k) = TE_sigma_Itot(i,k)/sqrt(Hx*Hy);
%     TE_Itot_sigma(i,k) = ete_hist(sigma(k:k+windowSizeDays-1),Itot(k:k+windowSizeDays-1),1,ceil(sqrt(windowSizeDays)),[-1 1]);
%     TE_Itot_sigma(i,k) = TE_Itot_sigma(i,k)/sqrt(Hx*Hy);
%     NetTE_sigma_Itot(i,k) = TE_sigma_Itot(i,k) - TE_Itot_sigma(i,k);
% 
% 
%     Hx = ent(alpha(k:k+windowSizeDays-1), ceil(sqrt(windowSizeDays)), [-1 1], 'x');
%     TE_alpha_Itot(i,k) = ete_hist(Itot(k:k+windowSizeDays-1),alpha(k:k+windowSizeDays-1),1,ceil(sqrt(windowSizeDays)),[-1 1]);
%     TE_alpha_Itot(i,k) = TE_alpha_Itot(i,k)/sqrt(Hx*Hy);
%     TE_Itot_alpha(i,k) = ete_hist(alpha(k:k+windowSizeDays-1),Itot(k:k+windowSizeDays-1),1,ceil(sqrt(windowSizeDays)),[-1 1]);
%     TE_Itot_alpha(i,k) = TE_Itot_alpha(i,k)/sqrt(Hx*Hy);
%     NetTE_alpha_Itot(i,k) = TE_alpha_Itot(i,k) - TE_Itot_alpha(i,k);
% 
% 
%     Hx = ent(kappa(k:k+windowSizeDays-1), ceil(sqrt(windowSizeDays)), [-1 1], 'x');
%     TE_kappa_Itot(i,k) = ete_hist(Itot(k:k+windowSizeDays-1),kappa(k:k+windowSizeDays-1),1,ceil(sqrt(windowSizeDays)),[-1 1]);
%     TE_kappa_Itot(i,k) = TE_kappa_Itot(i,k)/sqrt(Hx*Hy);
%     TE_Itot_kappa(i,k) = ete_hist(kappa(k:k+windowSizeDays-1),Itot(k:k+windowSizeDays-1),1,ceil(sqrt(windowSizeDays)),[-1 1]);
%     TE_Itot_kappa(i,k) = TE_Itot_kappa(i,k)/sqrt(Hx*Hy);
%     NetTE_kappa_Itot(i,k) = TE_kappa_Itot(i,k) - TE_Itot_kappa(i,k);
% 
% 
% 
% end
% 
% end