function [Hhat_MMSE,C_MMSE,tau_p,R,H,Hhat_EW_MMSE,C_EW_MMSE] = functionChannelEstimates(R,channelGaindB,nbrOfRealizations,M,K,L,p)
%Generate the channel realizations and estimates of these channels for all
%UEs in the entire network. The channels are assumed to be correlated
%Rayleigh fading. The MMSE estimator and EW-MMSE estimator are used.
%
%This Matlab function is used in the article:
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
%
%
%INPUT:
%R                 = M x M x K x L x L matrix with spatial correlation
%                    matrices for all UEs in the network. R(:,:,k,j,l) is
%                    the correlation matrix for the channel between UE k
%                    in cell j and the BS in cell l. This such matrix can
%                    either include the average channel gain or can be
%                    normalized arbitrarily.
%channelGaindB     = K x L x L matrix containing the average channel gains
%                    in dB of all the channels, if these are not already
%                    included in the spatial channel correlation matrices.
%                    The product R(:,:,k,j,l)*10^(channelGaindB(k,j,l)/10)
%                    is the full spatial channel correlation matrix.
%nbrOfRealizations = Number of channel realizations
%M                 = Number of antennas per BS
%K                 = Number of UEs per cell
%L                 = Number of BSs and cells
%p                 = Uplink transmit power per UE (same for everyone)
%
%OUTPUT:
%Hhat_MMSE    = M x nbrOfRealizations x K x L x L matrix with the MMSE
%               channel estimates. The matrix Hhat_MMSE(:,n,k,j,l) is the
%               n:th channel estimate of the channel between UE k in cell j
%               and the BS in cell l.
%C_MMSE       = M x M x K x L x L matrix with estimation error correlation
%               matrices when using MMSE estimation. The matrix is
%               organized in the same way as R.
%tau_p        = Length of pilot sequences
%R            = Scaled version of the input spatial correlation matrices R,
%               where the channel gains from channelGaindB are included
%H            = M x nbrOfRealizations x K x L x L matrix with the true
%               channel realizations. The matrix is organized in the same
%               way as Hhat_MMSE.
%Hhat_EW_MMSE = Same as Hhat_MMSE, but using the EW-MMSE estimator
%C_EW_MMSE    = Same as C_MMSE, but using the EW-MMSE estimator


%% Generate channel realizations

%Generate uncorrelated Rayleigh fading channel realizations
H = (randn(M,nbrOfRealizations,K,L,L)+1i*randn(M,nbrOfRealizations,K,L,L));

%Store identity matrix of size M x M
eyeM = eye(M);

%Length of pilot sequences
tau_p = K;


%Go through all channels and apply path gain and channel covariance
betas = zeros(K,L,L);
Rdiag = zeros(size(R));

for j = 1:L
    
    for l = 1:L
        
        for k = 1:K
            
            %Extract path gain
            betas(k,j,l) = 10^(channelGaindB(k,j,l)/10);
            
            %Apply channel gain to covariance matrix
            R(:,:,k,j,l) = betas(k,j,l) * R(:,:,k,j,l);
            
            %Apply covariance structure to the uncorrelated realizations
            Rsqrt = sqrtm(R(:,:,k,j,l));
            H(:,:,k,j,l) = sqrt(0.5)*Rsqrt*H(:,:,k,j,l);
            
            %Compute a diagonalized version of the covariance matrix
            Rdiag(:,:,k,j,l) = diag(diag(R(:,:,k,j,l)));
            
        end
        
    end
    
end



%% MMSE estimation

%Generate realizations of normalized noise
Np = sqrt(0.5)*(randn(M,nbrOfRealizations,K,L) + 1i*randn(M,nbrOfRealizations,K,L));

%Prepare to store MMSE channel estimates
Hhat_MMSE = zeros(M,nbrOfRealizations,K,L,L);

%Prepare to store estimation error covariance matrices
C_MMSE = zeros(M,M,K,L,L);

%Go through all cells
for j = 1:L
    
    %Compute processed pilot signal for all UEs
    Yp = sqrt(p)*tau_p*sum(H(:,:,:,:,j),4) + sqrt(tau_p)*Np(:,:,:,j);
    
    %Go through all UEs in cell j
    for k = 1:K
        
        %Compute the matrix that is inverted in the MMSE channel estimation
        denomMatrix = (p*tau_p*sum(R(:,:,k,:,j),4) + eyeM);
        
        %Go through all cells
        for l = 1:L
            
            %Estimate channel between BS l and UE k in cell j
            numdenomMatrix = R(:,:,k,l,j) / denomMatrix;
            Hhat_MMSE(:,:,k,l,j) = sqrt(p)*numdenomMatrix*Yp(:,:,k);
            
            %Compute corresponding error covariance matrix
            C_MMSE(:,:,k,l,j) = R(:,:,k,l,j) - p*tau_p*numdenomMatrix*R(:,:,k,l,j);
            
        end
        
    end
    
end


%% EW-MMSE estimation
if nargout >= 5

    %Prepare to store EW-MMSE channel estimates
    Hhat_EW_MMSE = zeros(M,nbrOfRealizations,K,L,L);
    
    %Prepare to store estimation error covariance matrices
    C_EW_MMSE = zeros(M,M,K,L,L);
    
    %Go through all cells
    for j = 1:L
        
        %Compute processed pilot signal for all UEs
        Yp = sqrt(p)*tau_p*sum(H(:,:,:,:,j),4) + sqrt(tau_p)*Np(:,:,:,j);
        
        for k = 1:K
            
            %Compute the matrix that is inverted in EW-MMSE estimation
            denomDiagMatrix = (p*tau_p*sum(Rdiag(:,:,k,:,j),4) + eyeM);
            
            %Compute the matrix that is inverted in MMSE estimation
            denomMatrix = (p*tau_p*sum(R(:,:,k,:,j),4) + eyeM);
            
            %Go through cells that use same pilot as cell l
            for l = 1:L
                
                %Estimate channel between BS l and UE k in cell j
                A_EW_MMSE = sqrt(p)*Rdiag(:,:,k,l,j) / denomDiagMatrix;
                Hhat_EW_MMSE(:,:,k,l,j) = A_EW_MMSE*Yp(:,:,k);
                
                %Compute corresponding error covariance matrix
                productAR = A_EW_MMSE * R(:,:,k,l,j);
                Cmatrix = R(:,:,k,l,j) - (productAR + productAR') * sqrt(p)*tau_p + tau_p*A_EW_MMSE*denomMatrix*A_EW_MMSE';
                
                %Make the error covariance matrix diagonal
                C_EW_MMSE(:,:,k,l,j) = diag(diag(Cmatrix));
                
            end
            
        end
        
    end
    
end
