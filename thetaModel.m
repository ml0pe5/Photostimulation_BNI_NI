function BNI=thetaModel(network,T,w,I_0,I_sig,flag)
%#codegen
%          Dynamics on a network of theta-'neural masses'
%
%  Call function: BNI=thetaModel(network,T,w,I_0,I_sig,flag)
%
%  Inputs: 
%     - network = connectivity matrix NxN (N<=150)
%     - T       = number of time steps (scalar)
%     - w       = global coupling (scalar)
%     - I_0     = distance to the SNIC bifurcation* (vector Nx1, N<=150)
%     - I_sig   = intensity of noise (scalar)
%     - flag    = choose 0 for 'BNI' or 1 for 'NI' 
%                 Necessary for correct coupling normalisation 
% 
%  Outputs:
%     - BNI     = BNI of each node 
%
%  This function does not contain rng('shuffle')!!
%
%  In the case of directed networks, network(i,j)>0 means that 
%    there is a connection from i to j. 
%
%  *The SNIC bifurcation occurs at I_0 = 0
%   At I_0 < 0 there is a stable and an unstable point
%   At I_0 > 0 there is a limit cycle
%   If I_0 is a scalar, then all nodes are at the same distance
%
%  Typical values:
%    T     = 4*10^6;
%    w     = 10;
%    I_0   = -1.2*ones(N,1);
%    I_sig = 5*1.2*0.1;
%
% M.A.Lopes, 2017
if T < 40000*10^2
    disp('Warning: T should be large enough to get a reliable BNI.')
    disp('Suggestion: Use T > = 4*10^6.')
end

dt=10^-2;  % time step
threshold=0.9; % threshold for BNI
window_epochs=6*4/dt; % window for BNI

I_sig=I_sig/sqrt(dt); 
N=length(network); % number of nodes

% normalisation of coupling
if flag==0 %'BNI'
    wnet=w*network/N;
elseif flag==1 % 'NI'
    wnet=w*network/(N+1);
else
    disp('The flag argument should be either 0 (BNI) or 1 (NI).');  
    BNI=[];
    return;
end
wnet=wnet';

BNI=zeros(N,1);
signal=zeros(1,N);
x=false(T,N);
theta_s=-real(acos((1+I_0)./(1-I_0))); % stable point if I_0 < 0
theta_old=theta_s; % initial condition  

% Compute time series
for time=1:T-1
    I=I_0+I_sig*randn(N,1)+wnet*(1-cos(theta_old-theta_s));
    theta_new=theta_old+dt*(1-cos(theta_old)+(1+cos(theta_old)).*I);
    signal(1,:)=0.5*(1-cos(theta_old-theta_s));
    x(time+1,:)=signal>threshold;
    theta_old=theta_new;
end

% Compute BNI
for node=1:N
    aux=find(x(:,node));
    if numel(aux)==0
        BNI(node,1)=0;
    else
        seizure_index=zeros(length(aux),2);
        seizure_index(1,1)=aux(1);
        k=1;
        for i=2:length(aux)
            if aux(i)-aux(i-1)>window_epochs
                seizure_index(k,2)=aux(i-1);
                k=k+1;
                seizure_index(k,1)=aux(i);
            end
        end
        seizure_index(k,2)=aux(end);
        seizure_index(k+1:end,:)=[];
        time_seizure=0;
        for i=1:size(seizure_index,1)
            time_seizure=time_seizure+seizure_index(i,2)-seizure_index(i,1)+1;
        end
        BNI(node,1)=time_seizure/T;
    end
end
