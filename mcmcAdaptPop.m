% Q1e - function to perform adaptive MCMC for all of the parameters
function [output] = mcmcAdaptPop(iters,startingValue,sigma0,n0,data)
% define acceptance rate monitors
accept=0;
reject=0;
reject1 = 0;
reject2 = 0;
% current value for the MCMC is theta
theta=startingValue;
% initial proposal covariance matrix is sigma0
% current estimates of posterior covariance matrix and mean are sigma and
% mu
sigma=sigma0;
mu=0;
epsilon=1e-6; % a small value -- zero will probably be ok.
% prepare storage of output
stored=zeros(iters,length(theta));
for i=1:iters
    % Generate a new proposal for beta - use optimal scaling result
    proposal=mvnrnd(theta,((2.38^(2))/3)*sigma+epsilon*eye(3));
    % check proposal is in range
    if max(proposal)<=10 && proposal(1)>=0 && proposal(2)>=0 && proposal(3)<=0 && proposal(3)>=-1 
        % calculate log acceptance ratio (may be bigger than 1)
        
        lar = LogLikelihoodBeta(data,proposal) - LogLikelihoodBeta(data,theta);       
        % since all priors are uniform, they cancel out and hence log
        % priors don't need to be added

        
        lar = min(lar,0);
        % generate a random number between 0 and 1;
        u = unifrnd(0,1);
        % accept if lar>log(u) iff ap>u 
        if lar>log(u)
            % the proposal becomes the new value of theta
            theta=proposal;
            accept=accept+1;
        else
            reject1=reject1+1;
        end
    else
        % automatically reject outside the range for theta (as it has prior
        % density zero).
        reject2=reject2+1;
    end
    reject = reject1+reject2;
    % store parameters for output every iteration (at the moment)
    stored(i,:)=theta;
    % recalculate and store trajectory for output
  


    % Update our estimate of sigma using the AM algorithm
    if i==n0
        mu=mean(stored(1:i,:),1); % mean over first component
        sigma=cov(stored(1:i,:))+epsilon*eye(3);
        
    elseif i>n0
        muPrevious=mu;
        mu=(i*mu+stored(i,:))/(i+1);

        sigma=((i-1)/i)*sigma+(stored(i,:)'*stored(i,:)+i*(muPrevious'*muPrevious)-(i+1)*(mu'*mu)+epsilon*eye(3))/i;
        % Matlab is sensitive to even tiny amounts of non-symmetry from numerical errors - this trick enforces symmetry 
        sigma=(sigma+sigma')/2;
        
    end      
end

% return the whole output matrix plus accept and reject counters
output=struct('theta',stored,'accept',accept,'reject',reject,'acceptanceRate',accept/(accept+reject),"sigma",sigma);


