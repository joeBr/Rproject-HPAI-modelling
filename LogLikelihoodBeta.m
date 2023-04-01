function [ll] = LogLikelihoodBeta(pop_obs,theta)

%Extract the parameter values from the list
alpha = theta(1);
beta = theta(2);
loc = theta(3); 
scale = 6.7;


    

%Initialise the log likelihood
ll = 0;

M = max(pop_obs);
m = min(pop_obs);

test1 = 1-(M-loc)/scale;



if test1>0
    %Loop through the data, adding the log density of that data points
    %value to the log likelihoods value
    for i = 1:length(pop_obs)
        k = pop_obs(i);
        log_pmf = (alpha-1)*log((k-loc)/scale) + (beta-1)*log(1-(k-loc)/scale) - log(gamma(alpha)) - log(gamma(beta)) + log(gamma(beta+alpha));
         
        ll = ll + log_pmf;
    end


else
    ll = -Inf;
end


end