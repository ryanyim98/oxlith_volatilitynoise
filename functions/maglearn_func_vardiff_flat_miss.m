function out=maglearn_func_vardiff_flat_miss(input,params)

% Optimal Bayesian Magnitude Learner with 5 variables
%
%
% Adapted from volatility learner of Meltem Sevgi.
% Also see the original .cpp code by Tim Behrens for the volatility learner.
% This learner estimates 5 variables from the following generative model
%
% Data (input) are continuous numbers between 0 and 1 (if the maximum is above 1 
% the model will assume the data are 0 to 100 and scale it accordingly).
% Data are transformed using the
% logit transform so that they are on the real line.

% This model assumes the (transformed) data are generated from a gaussian
% process with mean mu and SD exp(s).  The mu can change over time with probabilty
% of mu at time i+1 being given by a gaussian function centred on mu at
% time i and with a SD of exp(vmu). vmu itself changes over time with probability of value at
% time i+1 being given by a gaussian with mean of vmu at time i and a SD of
% exp(kmu). NB the part above is similar to Tim's volatility
% learner with the expception that the process is a gaussian rather than
% beta.

% In addition to the above s (log of the SD of the generative gaussian) can
% change over time. s at i+1 is given by s at i with sd exp(vs). 


% The params structure contains variables that influence the performance of
% the model. The model maintains a 5d pdf of its belief. The main descisions 
% that need to be made when running the model are the range of values of each
% parameter the model will consider and the number of values of each paramter
% it will use in the pdf. Generally wider ranges and more values for each parameter are
% best, but the higher the number of points, the more memory will be required
% and the longer the model will take to run. Similarly, there is little point in having
% values of each dimension outside those that are possible given the data/task.
% Our approach is to initially run the model with a coarse grid and, when 
% we are happy that it is working, run it on a linux server with higher
% numbers. For the data we have tested on, the estimates of the learner are stable
% from pdf sizes quite a bit below those we use (but this will probably depend
% on the nature of the data).Anyway:
 
%the dim_size parameter controls the number of points on each dimension (in
%the order mu, vmu, kmu, s, vs).

% the "range" parameters control the minimum and maximum of each dimension.
% The values chosen are well above/ below the values estimated for data we
% have tested on.

%flattenpoints causes the mu dimension to be flattened (keep information in
%the other dimensions) at specific time points. It isn't relevant
%for this study.

% The learner deals with missing data by performing the information leak (transitional
% matrices dealing with structure of world) in the absence of likelihood
% updating. Missing data should be represented as a NaN

%NB this version was written with a view to translating it to c++ which I thought
% might be more efficient (I did this but it didn't make the code more efficient).
% Basically this means that things like normal PDFs are written out
% explicitly rather than using matlab functions like normpdf.

%It might be possible to code this model using a particle filter, which
%could be more efficient in terms of memory/time. I haven't done this (as of 
%June 2020).

%default values for params. This will run quite a slow, high memory version
%of the model
if nargin<2; params=struct; end
if ~isfield(params,'murange'); params.murange=[inv_logit(-5,1) inv_logit(5,1)]; end
if ~isfield(params,'vmurange'); params.vmurange=[e-5 1]; end
if ~isfield(params,'kmurange'); params.kmurange=[5e-5 100]; end
if ~isfield(params,'srange'); params.srange=[0.01 10]; end
if ~isfield(params,'vsrange'); params.vsrange=[0.001 100]; end
if ~isfield(params,'dimsize'); params.dimsize=[36 27 21 27 27]; end
if ~isfield(params,'flattenpoints'); params.flattenpoints=[]; end

% define ranges for each dimension. mu is in logit space, other dimensions
% are in log space
out.muvec = inv_logit(params.murange(1)):((inv_logit(params.murange(2))-inv_logit(params.murange(1)))/(params.dimsize(1)-1)):inv_logit(params.murange(2));  
out.vmulog = log(params.vmurange(1)):((log(params.vmurange(2))-log(params.vmurange(1)))/(params.dimsize(2)-1)):log(params.vmurange(2));
out.kmulog = log(params.kmurange(1)):((log(params.kmurange(2))-log(params.kmurange(1)))/(params.dimsize(3)-1)):log(params.kmurange(2));
out.slog=log(params.srange(1)):((log(params.srange(2))-log(params.srange(1)))/(params.dimsize(4)-1)):log(params.srange(2));
out.vslog=log(params.vsrange(1)):((log(params.vsrange(2))-log(params.vsrange(1)))/(params.dimsize(5)-1)):log(params.vsrange(2));

% number of points on each dimension
out.musize=params.dimsize(1);
out.vmusize = params.dimsize(2);
out.kmusize = params.dimsize(3);
out.ssize=params.dimsize(4);
out.vssize=params.dimsize(5);

%number of timepoints in data
li=length(input);

out.ntrials=li;

%initialise output data
% for each dimension we output the marignal distribution (dist) and the expected
% value(est)
out.muDist=zeros(li,params.dimsize(1)); 
out.muEst=zeros(li,1);
out.vmuDist=zeros(li,params.dimsize(2));
out.vmuEst=zeros(li,1);
out.kmuDist=zeros(li,params.dimsize(3));
out.kmuEst=zeros(li,1);
out.sDist=zeros(li,params.dimsize(4));
out.sEst=zeros(li,1);
out.vsDist=zeros(li,params.dimsize(5));
out.vsEst=zeros(li,1);
out.volnoise=zeros(li,params.dimsize(2),params.dimsize(4)); % a 2 dimensional marginal pdf which allows visulisation of both volatility and noise
out.KLdiv=zeros(li,1); % the KL divergence between time points (an estimate of the information content of events)
out.entropy=zeros(li,1); % the entropy of the pdf



%%%%%%%%%%%%%%%%%%%%Transform Input Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data are assumed to lie either between 0 and 1 or 0 and 100. In the later
% case they are first transformed to 0 to 1 scale. In both cases the data
% are then transformed using inv_logit to the infinite real line (NB the actual
% data need to be >0 and <1 with this transform).

reward=input;
% if the max value fed to the model is greater than 1, it assumes it is a
% 0-100 scale (otherwise assume 0-1)
if max(reward)>1
    reward=reward./100;
end
% transform data to real line
trans_reward=inv_logit(reward);



%%%%%%%%%%%%%%%%%% Initialise Joint Distribution %%%%%%%%%%%%%%%%%%%%%%%%

% joint dist is the full joint distribution of all parameters. Initialised
% to a flat prior
jointdist = (ones(out.musize, out.vmusize, out.kmusize,out.ssize,out.vssize))./(out.musize*out.vmusize*out.kmusize*out.ssize*out.vssize);

jd_old=jointdist; % use this to hold the previous joint dist to calculare the KL divergence



%%%%%%%%%%%%%%%%%%%Calcuate Transistional Matrices for Updating Joint Dist
%%%%%%%%%%%%%%%%%%%

% these are used to account for the changes in the parameter values between
% time points.


% precompute p(vmui+1|vmui,k) for every vmu_{i+1}, vmu_i, kmu

vmup1gvmukmu = zeros(out.vmusize, out.vmusize, out.kmusize); % 3d array for vmu transition.
tmpvmu = zeros(out.vmusize, out.vmusize);  % to keep N(I_{i+1},k2

for k = 1 : out.kmusize
    for I = 1:out.vmusize
        for Ip1 = 1:out.vmusize
            var = exp(out.kmulog(k)*2); % k is stdev
            % below is the probability density function of a gaussian with
            % mean vmu and SD kmu at point vmu(i+1). It essentially provides the probability of
            % vmup1 given vmu and kmu
            tmpvmu(Ip1, I) = (exp(-power((out.vmulog(I) - out.vmulog(Ip1)),2)/(2*var))) / (sqrt(2*pi*var));
        end
        % normalise across range of vmu;
        tmpvmu(:,I) = tmpvmu(:,I)./sum(tmpvmu(:,I)); % this is a right stochastic (transitional) matrix.
    end
    %Probability of vmu on next trial given vmu and kmu on current trial.
    vmup1gvmukmu(:,:,k) = tmpvmu;
end


% precompute p(mui+1|mui,vmui+1) for every mu_{i+1}, mu_i, vmu_{i+1}
rmup1muvmup1 = zeros(out.musize, out.musize, out.vmusize);
tmpp = zeros(out.musize, out.musize);
for Ip1 = 1:out.vmusize
    for r = 1:out.musize
        for rp1 = 1:out.musize
            %Normal distribution of mu given previous mu and vmu
            var=exp(out.vmulog(Ip1)*2); % variance of normal pdf
            tmpp(rp1,r)=(exp(-power((out.muvec(r) - out.muvec(rp1)),2)/(2*var))) / (sqrt(2*pi*var));
            
        end
        tmpp(:,r) = tmpp(:,r)./sum(tmpp(:,r)); % normalise across each range of mu
    end
    rmup1muvmup1(:,:,Ip1) = tmpp; % place tmpp in it. Probability of r on next trial given r on this trial and I on next trial
end


% precompute update on s as above
rsp1svsp1 = zeros(out.ssize, out.ssize, out.vssize);
tmpp = zeros(out.ssize, out.ssize);
for vsp1 = 1:out.vssize
    for ss = 1:out.ssize
        for ss1 = 1:out.ssize
            %Normal distribution of mu given previous mu and vmu
            var=exp(out.vslog(vsp1)*2); % variance of normal pdf
            
            tmpp(ss1,ss)=(exp(-power((out.slog(ss) - out.slog(ss1)),2)/(2*var))) / (sqrt(2*pi*var));
            
        end
        tmpp(:,ss) = tmpp(:,ss)./sum(tmpp(:,ss)); % normalise across each range of s
    end
    rsp1svsp1(:,:,vsp1) = tmpp; % place tmpp in it. Probability of s on next trial given s on this trial and vs on next trial
end




% initialise arrays
pp1Ip1k = zeros(out.musize, out.vmusize, out.kmusize,out.ssize,out.vssize);
pvmup1kmu = zeros(out.musize, out.vmusize, out.kmusize,out.ssize,out.vssize);

%%%%%%%%%%%%Iterate through Trials, Updating Joint Distribution %%%%%%%%


for trial = 1:length(reward) 
    
    if sum(trial==params.flattenpoints)>0 % not relevant unless you want to flatten the mu dimension
        jointdist=repmat(mean(jointdist),out.musize,1,1,1,1);
        jd_old=jointdist;
    end
    
    
    %calculate entropy of pdf
    out.entropy(trial,1)=jointdist(jointdist>eps)'*(-log2(jointdist(jointdist>eps)));  % numbers less than 1e-15 are treated as 0
    
    % GET MARGINALS of each dimension
    %
    
    % mu
    
    out.muDist(trial,:) = sum(sum(sum(sum(jointdist,2),3),4),5); % sum over each row.
    out.muEst(trial,:)  = inv_logit(sum(out.muDist(trial,:).*out.muvec),1); % transform back to reward space after estimating expected value
    
    % vmu(volatility)
    out.vmuDist(trial,:) = sum(sum(sum(sum(jointdist,1),3),4),5);
    out.vmuEst(trial,:)  = sum(out.vmuDist(trial,:).*out.vmulog);
    
    %kmu (k for mu)
    out.kmuDist(trial,:) = sum(sum(sum(sum(jointdist,1),2),4),5);
    out.kmuEst(trial,:)  = sum(out.kmuDist(trial,:).*out.kmulog);
    
    % sd (expected uncertainty)
    out.sDist(trial,:)=sum(sum(sum(sum(jointdist,1),2),3),5);
    out.sEst(trial,:)=sum(out.sDist(trial,:).*out.slog);
    
    %sd vol
    out.vsDist(trial,:)=sum(sum(sum(sum(jointdist,1),2),3),4);
    out.vsEst(trial,:)=sum(out.vsDist(trial,:).*out.vslog);
    
    % surface of volatility vs. noise (i.e. the two sources of variabiity)
    out.volnoise(trial,:,:)=sum(sum(sum(jointdist,1),3),5);
    
    
    
    
    %%%%%%%%%%% Perform BAYESIAN UPDATE %%%%%%%%%%%%%%%%%%%%%%%
    %
    % this is liklihood of the observation given the model, multiplied by the
    % prior (i.e. joint distribution from previous trial). No correction for
    % other parameters yet.
    %
    if ~isnan(trans_reward(trial))
        for ss=1:out.ssize
            var=exp(out.slog(ss)*2); % var of normal pdf
            for r=1:out.musize
                jointdist(r,:,:,ss,:)=jointdist(r,:,:,ss,:).*((exp(-power((out.muvec(r) - trans_reward(trial)),2)/(2*var))) / (sqrt(2*pi*var)));
                               
            end
        end
    end
    
    % now do normalization (note likelihood does not sum to 1)
    jointdist = jointdist ./ sum(sum(sum(sum(sum(jointdist)))));
    
    %
    % Now account for other parameters. First deal with the change in vmu
    %
    
    % I) multiply jointdist (after bayes update) by vmup1gvmukmu (probability of vmu on
    %next trial given vmu and kmu on this trial), and integrate out vmu on this trial. This
    %will give pvmup1kmu (probability of vmu on next trial given kmu).

    
    for Ip1 = 1:out.vmusize
        tmp3=repmat(vmup1gvmukmu(Ip1,:,:),[1 1 1 out.ssize,out.vssize]);
        for r = 1:out.musize
            pvmup1kmu(r,Ip1,:,:,:) = sum(tmp3.*jointdist(r,:,:,:,:)); % for a given k calcuate the probability of r given the updated I
        end
    end
    
    % II) multiply by
    %rmup1muvmup1 (probability of mu on next trial given mu on this trial and
    %I on next trial), and integrate out mu on this trial. This will give pp1Ip1k
    %(probability of mu on next trial given I on next trial and kmu).
    
    for Ip1 = 1:out.vmusize
        tmp=pvmup1kmu(:,Ip1,:,:,:);
        for rp1 = 1:out.musize
            pp1Ip1k(rp1,Ip1,:,:,:) = sum(tmp.*repmat(rmup1muvmup1(rp1,:,Ip1)',[1 1 out.kmusize out.ssize,out.vssize]),1);
        end
    end
   
    % III) as above for the change in sd (rsp1svsp1= prob of s on next
    % trial given s on this trial and vs).
    
    for Ip1 = 1:out.vssize
        tmp2=pp1Ip1k(:,:,:,:,Ip1);
        for rp1 = 1:out.ssize
            jointdist(:,:,:,rp1,Ip1) = sum(tmp2.*permute(repmat(rsp1svsp1(rp1,:,Ip1)',[1 1 out.musize,out.vmusize,out.kmusize]),[3 4 5 1 2]),4);
        end
    end
    
    % calculate KL divergence
    idx=jointdist>eps & jd_old>eps;  % definition of kldiv when prob is 0 is 0 (faster doing it this way than across whole pdf)
   out.KLdiv(trial,1)=sum(sum(sum(sum(sum(jointdist(idx).*(log2(jointdist(idx))-log2(jd_old(idx))))))));  
 
    jd_old=jointdist;
    
end


%get last estimates (i.e. after the last update)
trial=trial+1;


    
    %entropy
    out.entropy(trial,1)=sum(sum(sum(sum(sum(jointdist(jointdist>eps).*(log2(1./(jointdist(jointdist>eps)))))))));  % numbers less than 1e-15 are 0
    
    % GET MARGINALS
    %
    
    % mu
    
    out.muDist(trial,:) = sum(sum(sum(sum(jointdist,2),3),4),5); % sum over each row.
    out.muEst(trial,:)  = inv_logit(sum(out.muDist(trial,:).*out.muvec),1); % transform back to reward space after estimating expected value
    
    % vmu(volatility)
    out.vmuDist(trial,:) = sum(sum(sum(sum(jointdist,1),3),4),5);
    out.vmuEst(trial,:)  = sum(out.vmuDist(trial,:).*out.vmulog);
    
    %kmu (k for mu)
    out.kmuDist(trial,:) = sum(sum(sum(sum(jointdist,1),2),4),5);
    out.kmuEst(trial,:)  = sum(out.kmuDist(trial,:).*out.kmulog);
    
    % sd (expected uncertainty)
    out.sDist(trial,:)=sum(sum(sum(sum(jointdist,1),2),3),5);
    out.sEst(trial,:)=sum(out.sDist(trial,:).*out.slog);
    
    %sd vol
    out.vsDist(trial,:)=sum(sum(sum(sum(jointdist,1),2),3),4);
    out.vsEst(trial,:)=sum(out.vsDist(trial,:).*out.vslog);
    
    % surface of volatility vs. noise
    out.volnoise(trial,:,:)=sum(sum(sum(jointdist,1),3),5);




