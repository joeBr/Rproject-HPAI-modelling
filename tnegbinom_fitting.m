data = importdata("data\nonzero_poultry_counts.txt");
%mle(data,'Distribution','burr')
pf_truncnbin = @(x1,r,p) nbinpdf(x1,r,p)./(1-nbincdf(0 ,r,p));
cf_truncnbin = @(x1,r,p) (nbincdf(x1,r,p)-nbincdf(0,r,p))./(1-nbincdf(0 ,r,p));
[start,ci] = mle(data,'Distribution','nbin');
[phat,pci] = mle(data,'pdf',pf_truncnbin,"cdf", cf_truncnbin,"Start",start,"LowerBound", [0,0], "UpperBound",[100,0.2],'TruncationBounds',[1,inf]);


phat(1);
phat(2);

bins = [];
for i = [0:6]
    for n = [1:9]
        bins(length(bins)+1) = n*10^i;
    end
end



figure
histogram(data,'Normalization','pdf', "BinEdges", bins)
set(gca, 'xscale','log')
xlabel("Poultry Premise Population")
    ylabel("Normalised Frequency/Fitted PMF")
hold on

plot(bins,pf_truncnbin(bins,phat(1),phat(2)))
legend('Observed Samples','Fitted Distribution')
hold off





figure
cdfplot(data)
set(gca, 'xscale','log')
hold on

plot(bins,cf_truncnbin(bins,phat(1),phat(2)))
legend('Observed Samples','Fitted Distribution')
hold off