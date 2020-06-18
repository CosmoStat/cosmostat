clear all

n  = 256;
MC = 1E3;
B  = 1E2;
sigma = 2.^linspace(-8,8,20);
alpha = linspace(1/B,1-1/B,100);
sigalpha = sqrt(alpha.*sqrt(1-alpha)/MC);

% Chi-square test of variance homogeneity (H0).
for rep=1:MC
 x = fft(randn(n,1))/sqrt(n);
 %[tmp,p(rep)] = vartest(x,1,0.05,'both');
 crithd(rep) = sum(abs(x - mean(x)).^2);
 tmp = gammainc(crithd(rep)/2,(n-1)/2);
 p(rep)=2*min(tmp,1-tmp);
 out = x(ceil(n*rand(n,B)));
 size(sum(abs(out - repmat(mean(out),n,1))))
 crithdb(rep,:) = sum(abs(out - repmat(mean(out),n,1)).^2)/var(x);
end

size(crithdb)

crithdb=sort(crithdb')';
for ia=1:length(alpha)
 q=ceil((B+1)*(1-alpha(ia)));
 %crtl=chi2inv(1-alpha(ia),n-1);
 %alphaH0(ia)  = length(find(p<=alpha(ia)))/MC;
 alphaest(ia) = length(find(crithd'>=crithdb(:,q)))/MC;
end

pause

subplot(121)
%errorbar(alpha,alphaH0,3*sigalpha);hold on
plot(alpha,alpha,'--k',alpha,alphaest,'b',alpha,alphaH0,'-.b');axis tight;hold on
p=patch([alpha fliplr(alpha)],[alphaest+3*sigalpha fliplr(alphaest-3*sigalpha)],[0.5 0 0]);set(p,'EdgeColor','none');
p=patch([alpha fliplr(alpha)],[alphaH0+3*sigalpha fliplr(alphaH0-3*sigalpha)],[0.5 0.5 0.5]);set(p,'EdgeColor','none','FaceAlpha',0.3);
%area(alpha,[alphaest'+3*sigalpha'],'FaceColor',[0.5 0 0],'EdgeColor','None');
%area(alpha,[alphaest'-3*sigalpha'],'FaceColor',[0.99 0.99 0.99],'EdgeColor','None');
%area(alpha,[alphaH0'+3*sigalpha'],'FaceColor',[0.5 0.5 0.5],'EdgeColor','None');
%area(alpha,[alphaH0'-3*sigalpha'],'FaceColor',[1 1 1],'EdgeColor','None');
set(gca,'Layer','Top');
plot(alpha,alpha,'--k',alpha,alphaest,'b',alpha,alphaH0,'-.b');hold off
legend('Expected','Bootstrap','Gaussian','+- 3 Error bars Bootstrap','+- 3 Error bars Gaussian','Location','Best');
xlabel('Type I error \alpha');ylabel('Observed FPF');
title('Specificity calibration under H_0: homogeneous variance');drawnow

% Chi-square test of variance homogeneity (H1).
for isig=1:length(sigma)
 for rep=1:MC
  x = fft(randn(n,1))/sqrt(n);
  pos = ceil(n*rand(1));
  x(pos) = sigma(isig)*x(pos);
  crithd(rep) = sum(abs(x - mean(x)).^2);
  out = x(ceil(n*rand(n,B)));
  crithdb(rep,:) = sum(abs(out - repmat(mean(out),n,1)).^2)/var(x);
  %tmp = gammainc(sum(abs(x-mean(x)).^2)/2,(n-1)/2);
  %p(rep)=2*min(tmp,1-tmp);
 end
 crithdb=sort(crithdb')';
 q=ceil((B+1)*(1-0.05));
 crtl=chi2inv(1-0.05,n-1);
 %alphaH1(isig)  = length(find(p <= 0.05))/MC;
 alphatheH1(isig)  = length(find(crithd'>=crtl))/MC;
 alphaestH1(isig)  = length(find(crithd'>=crithdb(:,q)))/MC;
end

subplot(122)
plot(log2(sigma),alphaestH1,'b',log2(sigma),alphatheH1,':b');grid on;axis tight;
legend('Bootstrap','Gaussian');
xlabel('log_2(\sigma)');ylabel('Observed TPF');
title('Sensitivity calibration under H_1 (\alpha=0.05): heteroscedastic profile at one sample');
