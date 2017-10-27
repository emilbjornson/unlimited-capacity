%This Matlab script can be used to generate all the figures in the article:
%
%Emil Bjornson, Jakob Hoydis, Luca Sanguinetti, ?Massive MIMO has Unlimited
%Capacity,? IEEE Transactions on Wireless Communications, to appear.
%
%Download article: https://arxiv.org/pdf/1705.00538
%
%This is version 1.0 (Last edited: 2017-10-27)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original article listed above.


%Initialization
close all;
clear;


%% Define simulation scenario

%Select which figure to generate:
%
%Set simulation = 1 to generate Figure 4
%Set simulation = 2 to generate Figure 5
%Set simulation = 3 to generate Figure 6
%Set simulation = 4 to generate Figure 7
simulation = 4;


%Number of BSs in the area
L = 4;


%Define the variables in each of the simulations
if simulation == 1
    
    %Number of UEs per BS
    K = 2;
    
    %Select range of number of BS antennas
    Mrange = round(logspace(1,3,10));
    
    %Correlation factor in exponential correlation model
    correlationFactor = 0.5;
    
    %Number of points on the horizontal axis in the final figure
    nbrOfPoints = length(Mrange);
    
    %Number of setups with randomly generated channel statistics (only one
    %is needed in this case since the statistics are deterministic)
    nbrOfSetups = 1;
    
    %Number of different covariance matrix realizations per setup
    statisticsVariations = 1;
    
    %Determine maximum number of BS antennas
    Mmax = max(Mrange);
    
    %Select number of channel realizations per setup
    nbrOfChannelRealizations = 5000;
    
    %Standard deviation of large-scale fading variations over the array
    stdLSF = 0;
    
elseif simulation == 2
    
    %Number of UEs per BS
    K = 2;
    
    %Number of BS antennas
    M = 200;
    
    %Determine maximum number of BS antennas
    Mmax = M;
    
    %Number of setups with randomly generated statistics
    nbrOfSetups = 50;
    
    %Standard deviation of large-scale fading variations over the array
    stdLSF = linspace(0,5,11);
    
    %Number of measurement points along horizontal axis in final figure
    nbrOfPoints = length(stdLSF);
    
    %Number of different covariance matrix realizations per setup
    statisticsVariations = length(stdLSF);
    
    %Select number of channel realizations per setup
    nbrOfChannelRealizations = 1000;
    
elseif simulation == 3
    
    %Number of UEs per BS
    K = 2;
    
    %Select range of BS antennas
    Mrange = round(logspace(1,3,10));
    
    %Correlation factor in exponential correlation model
    correlationFactor = 0.5;
    
    %Number of measurement points along horizontal axis in final figure
    nbrOfPoints = length(Mrange);
    
    %Number of setups with randomly generated statistics
    nbrOfSetups = 50;
    
    %Number of different covariance matrix realizations per setup
    statisticsVariations = 1;
    
    %Determine maximum number of BS antennas
    Mmax = max(Mrange);
    
    %Select number of channel realizations per setup
    nbrOfChannelRealizations = 500;
    
    %Standard deviation of large-scale fading variations over the array
    stdLSF = 4;
    
elseif simulation == 4
    
    %Number of UEs per BS
    K = 10;
    
    %Select range of BS antennas
    Mrange = round(logspace(1,3,10));
    
    %Correlation factor in exponential correlation model
    correlationFactor = 0.5;
    
    %Number of measurement points along horizontal axis in final figure
    nbrOfPoints = length(Mrange);
    
    %Number of setups with randomly generated statistics
    nbrOfSetups = 1;
    
    %Number of different covariance matrix realizations per setup
    statisticsVariations = 1;
    
    %Determine maximum number of BS antennas
    Mmax = max(Mrange);
    
    %Select number of channel realizations per setup
    nbrOfChannelRealizations = 500;
    
    %Standard deviation of large-scale fading variations over the array
    stdLSF = 4;
    
end



%% Scenario setup (see Figure 3)

%Distance between BSs
interBSdistance = 300;

%Define BS positions using complex coordinates
BSlocations = interBSdistance*exp(1i*(pi/4+pi/2*(1:4)))/sqrt(2);

%Deploy UEs at the fixed locations in Figure 3
celledgeDistance = 0.33*interBSdistance/2;

if simulation == 1 || simulation == 2 || simulation == 3
    
    UElocations = celledgeDistance*[exp(1i*(pi/2-pi/80+pi/2*(1:4))); exp(1i*(pi/80+pi/2*(1:4)))];
    
end



%% Propagation parameters

%Pathloss exponent
alpha = 3.76;

%Constant term in pathloss model (assuming distances in meters)
constantTerm = -35.3;

%Communication bandwidth
bandwidth = 10e6;

%Define total uplink transmit power per UE (mW)
p = 100;

%Define total downlink transmit power per UE (mW)
rho = 100;

%Define noise figure at BS (in dB)
noiseFigure = 10;

%Compute total noise power
noiseVariancedBm = -174 + 10*log10(bandwidth) + noiseFigure;

%Select length of coherence block
tau_c = 200;



%Prepare to save simulation results
meanSE_MR = zeros(nbrOfPoints,nbrOfSetups);
meanSE_SMMSE = zeros(nbrOfPoints,nbrOfSetups);
meanSE_MMMSE = zeros(nbrOfPoints,nbrOfSetups);
meanSE_MZF = zeros(nbrOfPoints,nbrOfSetups);

if simulation == 1 || simulation == 4
    
    meanSE_singlecell = zeros(nbrOfPoints,nbrOfSetups);
    
    powerlevelsTotal = zeros(4,3,nbrOfPoints,nbrOfSetups);
    
elseif simulation == 3 || simulation == 4
    
    meanSE_MR_EW = zeros(nbrOfPoints,nbrOfSetups);
    meanSE_SMMSE_EW = zeros(nbrOfPoints,nbrOfSetups);
    meanSE_MMMSE_EW = zeros(nbrOfPoints,nbrOfSetups);
    meanSE_MZF_EW = zeros(nbrOfPoints,nbrOfSetups);
    
end


%% Go through all setups
for s = 1:nbrOfSetups
    
    %Output simulation progress
    disp([num2str(s) ' setups out of ' num2str(nbrOfSetups)]);
    
    %Prepare to store normalized channel covariance matrices
    R = zeros(Mmax,Mmax,K,L,L,statisticsVariations);
    
    %Prepare to store pathloss numbers (in dB)
    pathgaindB = zeros(K,L,L);
    
    if simulation == 4
        
        %Distribute UEs at random in the cell edge area
        UElocations = celledgeDistance*(rand(K,L)+1i*rand(K,L))*diag(exp(1i*[pi/2 pi -pi/2 0]));
        
    end
    
    
    %Go through all the cells
    for l = 1:L
        
        for j = 1:L
            
            %Compute distances between UEs in cell l and BS j
            distancesBSj = abs(UElocations(:,l)-BSlocations(j));
            
            %Compute angles between UEs in cell l and BS j
            angleBSj = angle(UElocations(:,l)-BSlocations(j));
            
            %Compute distant-dependent path gains (in dB)
            pathgaindB(:,l,j) = -constantTerm - alpha*10*log10(distancesBSj);
            
            
            %Compute angles between UEs in cell l and BS j, and generate
            %covariance matrices for all channels
            for k = 1:K
                
                %Generate realizations for large-scale fading variations
                fadingOverArray = randn(Mmax,1);
                
                %Go through all variations of the statistics
                for r = 1:statisticsVariations
                    
                    if simulation == 1 || simulation == 3 || simulation == 4
                        
                        %Generate large-scale fading variations over the array
                        largeScaleFadingD = diag(10.^(stdLSF*fadingOverArray/20));
                        
                        %Generate covariance matrix using the exponential
                        %correlation model, including large-scale
                        %variations over the array
                        R(:,:,k,l,j,r) = largeScaleFadingD*toeplitz((correlationFactor(r)*exp(1i*angleBSj(k))).^(0:Mmax-1))*largeScaleFadingD;
                        
                    elseif simulation == 2
                        
                        %Generate diagonal covariance matrix with
                        %large-scale variations over the array
                        R(:,:,k,l,j,r) = diag(10.^(stdLSF(r)*fadingOverArray/10));
                        
                    end
                    
                end
                
            end
            
        end
        
    end
    
    
    %Compute the normalized channel gains, where the normalization is by
    %noise power
    channelGaindB = pathgaindB -noiseVariancedBm;
    

    %Go through all points on the horizontal axis in the figure to be
    %plotted
    for m = 1:nbrOfPoints
        
        %Output simulation progress
        disp([num2str(m) ' points out of ' num2str(nbrOfPoints)]);
        
        if simulation == 1 || simulation == 3 || simulation == 4
            
            %Extract current number of antennas
            M = Mrange(m);
            
            %Extract covariance matrices of current size 
            Rmatrix = R(1:M,1:M,:,:,:,1);
            
        elseif simulation == 2
            
            %Extract current covariance matrices
            Rmatrix = R(:,:,:,:,:,m);
            
        end
        
        
        if simulation == 1
            
            %Generate MMSE channel estimates
            [Hhat,C,taupu,Rscaled] = functionChannelEstimates(Rmatrix,channelGaindB,nbrOfChannelRealizations,M,K,L,p);
            
            %Compute uplink SEs with different combining schemes
            [SE_MR,SE_MMMSE,SE_SMMSE,SE_MZF] = functionComputeSE_UL(Hhat,C,Rscaled,tau_c,taupu,nbrOfChannelRealizations,M,K,L,p);
            
            
            %Compute MMSE channel estimates in single-cell operation
            [Hhat,C,taupu,Rscaled] = functionChannelEstimates(Rmatrix(:,:,:,1,1),channelGaindB(:,1,1),nbrOfChannelRealizations,M,K,1,p);
            
            %Compute uplink SEs with MMSE combining
            [~,SE_singlecell] = functionComputeSE_UL(Hhat,C,Rscaled,tau_c,taupu,nbrOfChannelRealizations,M,K,1,p);
            
            %Save simulation result in single-cell operation, where each
            %cell is only active in 1/L of the coherence blocks
            meanSE_singlecell(m,s) = mean(mean(SE_singlecell,1))/L;
            
            
        elseif simulation == 2
            
            %Generate MMSE channel estimates
            [Hhat,C,taupu,Rscaled] = functionChannelEstimates(Rmatrix,channelGaindB,nbrOfChannelRealizations,M,K,L,p);
            
            %Compute uplink SEs and power levels with different combining
            %schemes
            [SE_MR,SE_MMMSE,SE_SMMSE,SE_MZF,powerlevels] = functionComputeSE_UL(Hhat,C,Rscaled,tau_c,taupu,nbrOfChannelRealizations,M,K,L,p);
            
            %Save simulation results on the power levels
            powerlevelsTotal(:,:,m,s) = powerlevels;
            
            
        elseif simulation == 3 || simulation == 4
            
            %Generate channel realizations and channel estimates, using
            %MMSE and EW-MMSE
            [Hhat_MMSE,C_MMSE,taupu,Rmatrix,H,Hhat_EW_MMSE,C_EW_MMSE] = functionChannelEstimates(Rmatrix,channelGaindB,nbrOfChannelRealizations,M,K,L,p);
            
            %Compute downlink SEs with different precoding schemes, using
            %the MMSE estimates
            [SE_MR,SE_MMMSE,SE_SMMSE,SE_MZF] = functionComputeSE_DL(H,Hhat_MMSE,C_MMSE,Rmatrix,tau_c,taupu,nbrOfChannelRealizations,M,K,L,p,rho);
            
            %Compute downlink SEs with different precoding schemes, using
            %the EW-MMSE estimates
            [SE_MR_EW,SE_MMMSE_EW,SE_SMMSE_EW,SE_MZF_EW] = functionComputeSE_DL(H,Hhat_EW_MMSE,C_EW_MMSE,Rmatrix,tau_c,taupu,nbrOfChannelRealizations,M,K,L,p,rho);
            
            %Save simulation results when using EW-MMSE estimation
            meanSE_MR_EW(m,s) = mean(mean(SE_MR_EW,1));
            meanSE_SMMSE_EW(m,s) = mean(mean(SE_SMMSE_EW,1));
            meanSE_MMMSE_EW(m,s) = mean(mean(SE_MMMSE_EW,1));
            meanSE_MZF_EW(m,s) = mean(mean(SE_MZF_EW,1));
            
        end
        
        
        %Save simulation results when using MMSE estimation
        meanSE_MR(m,s) = mean(mean(SE_MR,1));
        meanSE_SMMSE(m,s) = mean(mean(SE_SMMSE,1));
        meanSE_MMMSE(m,s) = mean(mean(SE_MMMSE,1));
        meanSE_MZF(m,s) = mean(mean(SE_MZF,1));

    end
    
end


%% Plot simulation results

if simulation == 1
    
    
    %Plot Figure 4
    figure;
    hold on; box on;
    
    plot(Mrange,meanSE_MMMSE,'rd-','LineWidth',1);
    plot(Mrange,meanSE_SMMSE,'k*--','LineWidth',1);
    plot(Mrange,meanSE_MR,'bs-','LineWidth',1);
    plot(Mrange,meanSE_MZF,'ko:','LineWidth',1);
    plot(Mrange,meanSE_singlecell,'k:','LineWidth',1);
    
    xlabel('Number of antennas (M)');
    ylabel('Spectral efficiency [bit/s/Hz/user]');
    
    legend('M-MMSE','S-MMSE','MR','M-ZF','Time splitting','Location','NorthWest');
    set(gca,'XScale','log');
    
    
elseif simulation == 2
    
    
    %Plot Figure 5a
    figure;
    hold on; box on;
    
    meanSE_MZF(1,:) = 0;
    
    plot(stdLSF,mean(meanSE_MMMSE,2),'rd-','LineWidth',1);
    plot(stdLSF,mean(meanSE_MZF,2),'ko:','LineWidth',1);
    plot(stdLSF,mean(meanSE_SMMSE,2),'k*--','LineWidth',1);
    plot(stdLSF,mean(meanSE_MR,2),'bs-','LineWidth',1);
    
    xlabel('Standard deviation of fading variations over the array');
    ylabel('Spectral efficiency [bit/s/Hz/user]');
    
    legend('M-MMSE','M-ZF','S-MMSE','MR','Location','NorthWest');
    ylim([0 4]);
    
    
    %Plot Figure 5b
    figure;
    bar(10*log10(mean(powerlevelsTotal(:,:,9,:),4)'));
    legend('MR','S-MMSE','M-MMSE','M-ZF');
    ylabel('Received power over noise power [dB]');
    colormap(hot);
    
    labels = {'   Desired signal  ', 'Interf: Same pilot ', 'Interf: Diff. pilot'};
    set(gca, 'XTick', 1:3, 'XTickLabel', labels);
    
    
elseif simulation == 3
    
    
    %Plot Figure 6a
    figure;
    hold on; box on;
    
    plot(Mrange,meanSE_MMMSE,'rd-','LineWidth',1);
    plot(Mrange,meanSE_MZF,'ko:','LineWidth',1);
    plot(Mrange,meanSE_SMMSE,'k*--','LineWidth',1);
    plot(Mrange,meanSE_MR,'bs-','LineWidth',1);
    
    xlabel('Number of antennas (M)');
    ylabel('Spectral efficiency [bit/s/Hz/user]');
    
    legend('M-MMSE','M-ZF','S-MMSE','MR','Location','NorthWest');
    set(gca,'XScale','log');
    title('MMSE estimation');
    
    
    %Plot Figure 5b
    figure;
    hold on; box on;
    
    plot(Mrange,meanSE_MMMSE_EW,'rd-','LineWidth',1);
    plot(Mrange,meanSE_MZF_EW,'ko:','LineWidth',1);
    plot(Mrange,meanSE_SMMSE_EW,'k*--','LineWidth',1);
    plot(Mrange,meanSE_MR_EW,'bs-','LineWidth',1);
    
    xlabel('Number of antennas (M)');
    ylabel('Spectral efficiency [bit/s/Hz/user]');
    
    legend('Approximate M-MMSE','M-ZF','Approximate S-MMSE','MRT','Location','NorthWest');
    set(gca,'XScale','log');
    title('EW-MMSE estimation');
    
    
elseif simulation == 4
    
    
    %Plot Figure 7a
    figure;
    hold on; box on;
    
    plot(Mrange,meanSE_MMMSE,'rd-','LineWidth',1);
    plot(Mrange,meanSE_MZF,'ko:','LineWidth',1);
    plot(Mrange,meanSE_SMMSE,'k*--','LineWidth',1);
    plot(Mrange,meanSE_MR,'bs-','LineWidth',1);
   
    xlabel('Number of antennas (M)');
    ylabel('Spectral efficiency [bit/s/Hz/user]');
    
    legend('M-MMSE','M-ZF','S-MMSE','MRT','Location','NorthWest');
    set(gca,'XScale','log');
    title('MMSE estimation');
    
    
    %Plot Figure 7b
    figure;
    hold on; box on;
    
    plot(Mrange,meanSE_MMMSE_EW,'rd-','LineWidth',1);
    plot(Mrange,meanSE_MZF_EW,'ko:','LineWidth',1);
    plot(Mrange,meanSE_SMMSE_EW,'k*--','LineWidth',1);
    plot(Mrange,meanSE_MR_EW,'bs-','LineWidth',1);
    
    xlabel('Number of antennas (M)');
    ylabel('Spectral efficiency [bit/s/Hz/user]');
    
    legend('Approximate M-MMSE','M-ZF','Approximate S-MMSE','MRT','Location','NorthWest');
    set(gca,'XScale','log');
    title('EW-MMSE estimation');
    

end
