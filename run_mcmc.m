t = readtable("nonzero_poultry_logcounts.txt");
data = table2array(t(:,3));





%Calculate the posterior mode 
theta0 = [5,5,0];
fun2 = @(theta) -LogLikelihoodBeta(data,theta);
posterior_mode = fminsearch(fun2,theta0);


%Calculate the variances of each prior distribution to go into the initial
%value of sigma
var_alpha = 100/12;
var_beta = 100/12;
var_loc = 1/12;


%Create initial covariance matric for MCMC
sigma0 = diag([var_alpha,var_beta,var_loc]);



%First run starting from the posterior mode using sigma0 above
output = mcmcAdaptPop(2000,posterior_mode,sigma0,100,data);
%Calculate the intial theta and sigma used in the next run from the values
%from the first
theta = output.theta(end,:);
sigma=2.38^2*cov(output.theta(101:end))/3;

%Second run uses output from first run as inputs 
output = mcmcAdaptPop(2000,theta,sigma,100,data);
theta = output.theta(end,:);
sigma=2.38^2*cov(output.theta(101:end))/3;

%Third and final longer run to give final estimates
iters = 10000;
output = mcmcAdaptPop(iters,theta,sigma,100,data);
sigma = output.sigma;

%Calculate the median posterior from the data after the burn in period
theta_hat = median(output.theta(101:end,:),1);
U = quantile(output.theta(101:end,:),[0.975],1);
L = quantile(output.theta(101:end,:),[0.025],1);



% Make some trace plots so that we can examine the output
figure(2)
clf
%Plot beta
subplot(3,1,1);
plot(1:iters,output.theta(:,1));
ylabel('\alpha');
title(strcat('\sigma = ',string(sigma(1,1)),', acceptance rate = ',string(output.acceptanceRate)));
%Plot gamma
subplot(3,1,2);
plot(1:iters,output.theta(:,2));
ylabel('\beta');
title(strcat('\sigma = ',string(sigma(2,2)),', acceptance rate = ',string(output.acceptanceRate)));
%Plot iota
subplot(3,1,3);
plot(1:iters,output.theta(:,3));
xlabel("iteration");
ylabel("loc");
title(strcat('\sigma = ',string(sigma(3,3)),', acceptance rate = ',string(output.acceptanceRate)));



