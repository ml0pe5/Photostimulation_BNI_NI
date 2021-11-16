function [network_zero,network_lag]=plv_method(data,sr,fmin)
% This function calculates the connectivity matrix based on the methods in 
% [1,2], except for the correction of indirect connections (spurious 
% connections). It is suitable for scalp EEG data, artefact-free, 
% background activity. 
% Connections are inferred using the PLV (also called PLF). Two functional 
% networks are computed: with and without zero lag PLF.
%
% Statistically significance of the connections is assessed by comparing to
% 99 surrogates. 
%
% Call function: [network_zero,network_lag]=plv_both_method(data,sr,fmin)
%
% Input: 'data'   : matrix containing scalp EEG data
%         Typically 20 seconds, at a sampling frequency of 256 Hz.
%        'sr'    : sampling rate
%        'fmin'  : smallest frequency in the data (e.g. if the data was
%        filtered between 6 and 9 Hz, then fmin = 6).
%
% Output: 'network_zero': connectivity matrix (undirected) with zero lag        
%         'network_lag' : connectivity matrix (directed) without zero lag
%
% Note: output elements ( i , j ) represent a connection from i to j. 
%
% This function uses the following function developed by others:
%    1. generate_iAAFT_it : to generate univariate surrogates
%    2. PLF_lag           : to compute the PLV
% -----------------------------------------
%
% References:
% [1] Schmidt, Helmut, et al. "Dynamics on networks: the role of local 
% dynamics and global networks on the emergence of hypersynchronous neural 
% activity." PLoS Comput Biol 10.11 (2014): e1003947. 
%
% [2] Tait, Luke, et al., "A Large-Scale Brain Network Mechanism for 
% Increased Seizure Propensity in Alzheimer's Disease", Accepted in 
% PLoS Comput. Biol., bioRxiv:10.1101/2021.01.19.427236 (2021).
%
% ------------------------------------------
% % M.A. Lopes, m.lopes@exeter.ac.uk, 2017
% ------------------------------------------
% Current version: v4_2021-07-07
%

n_ch = min(size(data)); % number of channels
n_samp=max(size(data)); % number of sample points

res = 2*pi*fmin/sr; % see update v1

if size(data,1)>size(data,2)
    data=data'; % Rearrange data as channels x time
end
n_surr=99;  % number of surrogates
alpha=0.05; % threshold of statistical significance 

surr=zeros(n_ch,n_samp,n_surr);     % surrogates
PLV=zeros(n_ch,n_ch);               % PLV 
PLV_surr=zeros(n_ch,n_ch,n_surr);   % PLV of surrogates
PLV_Z=zeros(n_ch,n_ch);             % PLV with zero lag 
PLV_surr_Z=zeros(n_ch,n_ch,n_surr); % PLV with zero lag of surrogates
network_lag=zeros(n_ch,n_ch);       % connectivity matrix
network_zero=zeros(n_ch,n_ch);      % connectivity matrix (with zero lag)

% 1. Make univariate surrogates
rng('shuffle');
for n=1:n_surr
    for ch=1:n_ch
        [surr(ch,:,n),~] = generate_iAAFT_it(data(ch,:),1,10);% 10 iterations!
    end
end

% 2. Compute PLF
%   2.1 PLF of original data
data=data'; 
DataMean = mean(data); 
DataStd = std(data);
data = (data-repmat(DataMean,[n_samp,1]))./repmat(DataStd,[n_samp,1]);
%data=(data-mean(data))./std(data); % normalization
for ch1=1:n_ch-1
    for ch2=ch1+1:n_ch
        [plf,lag]=PLF_lag(data(:,ch1),data(:,ch2));
        if lag>res %update v1 (in v0, it was 0 instead of res)
            PLV(ch1,ch2)=plf;
        elseif lag<-res %update v1 (in v0, it was 0 instead of res)
            PLV(ch2,ch1)=plf;
        end
        PLV_Z(ch1,ch2)=plf;
        PLV_Z(ch2,ch1)=plf;
    end
end

%   2.2 PLF of surrogates
for n=1:n_surr
    surrogate=surr(:,:,n);
    surrogate=surrogate';    
    SurrMean = mean(surrogate);
    SurrStd = std(surrogate);
    surrogate = (surrogate-repmat(SurrMean,[n_samp,1]))./repmat(SurrStd,[n_samp,1]);   
    %surrogate=(surrogate-mean(surrogate))./std(surrogate);
    for ch1=1:n_ch-1
        for ch2=ch1+1:n_ch
            [plf,lag]=PLF_lag(surrogate(:,ch1),surrogate(:,ch2));
            if lag>res %update v1 (in v0, it was 0 instead of res)
                PLV_surr(ch1,ch2,n)=plf;
            elseif lag<-res %update v1 (in v0, it was 0 instead of res)
                PLV_surr(ch2,ch1,n)=plf;
            end
            PLV_surr_Z(ch1,ch2,n)=plf;
            PLV_surr_Z(ch2,ch1,n)=plf;
        end
    end
end
    
% 3. Statistical correction using the surrogates
for ch1=1:n_ch
    for ch2=1:n_ch
        p=sum(PLV(ch1,ch2)>[squeeze(PLV_surr(ch1,ch2,:));PLV(ch1,ch2)])/(n_surr+1); % see update v2
        if p>1-alpha
           network_lag(ch1,ch2)=PLV(ch1,ch2);
        end
        p=sum(PLV_Z(ch1,ch2)>[squeeze(PLV_surr_Z(ch1,ch2,:));PLV_Z(ch1,ch2)])/(n_surr+1); % see update v2
        if p>1-alpha
           network_zero(ch1,ch2)=PLV_Z(ch1,ch2);
        end
    end
end