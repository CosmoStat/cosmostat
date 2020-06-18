clear all

n  = 256;
MC = 1E3;
B  = 1E2;
sigma = 2.^linspace(0,12,10);
alpha = linspace(1/B,1-1/B,100);
sigalpha = sqrt(alpha.*sqrt(1-alpha)/MC);

% UMPI Mean constancy test of mean after the log on the periodogram.
crithdb = zeros(MC,B);
for rep=1:MC
 x   = fft(randn(n,1))/sqrt(n);
 per = abs(x).^2;
 logper = log(per(2:end/2-1));
 crithd(rep) = sum((logper - mean(logper)).^2)/psi(1,1);
 xb   = std(x)*complex(randn(n,B),randn(n,B))/sqrt(2);
 perb = abs(xb).^2;
 outper = log(perb(2:end/2-1,:));
 size(outper)
 pause
 crithdb(rep,:) = sum((outper - repmat(mean(outper),length(logper),1)).^2)/psi(1,1);
%  for b=1:B
%    xb   = std(x)*complex(randn(n,B),randn(n,B))/sqrt(2);
%    perb = abs(xb).^2;
%    outper = log(perb(2:end/2-1));
%    crithdb(rep,b) = sum((outper - mean(outper)).^2)/psi(1,1);
%  end
  %outper = logper(ceil(length(logper)*rand(length(logper),B)));
  %crithdb(rep,:) = sum((outper - repmat(mean(outper),length(logper),1)).^2)/var(logper);
end

crithdb=sort(crithdb')';
for ia=1:length(alpha)
 q=ceil((B+1)*(1-alpha(ia)));
 crtl=chi2inv(1-alpha(ia),length(logper)-1);
 alphathe(ia) = length(find(crithd'>=crtl))/MC;
 alphaest(ia) = length(find(crithd'>=crithdb(:,q)))/MC;
end

subplot(121)
plot(alpha,alpha,'--k',alpha,alphaest,'b',alpha,alphathe,'-.b');axis tight;hold on
p=patch([alpha fliplr(alpha)],[alphaest+3*sigalpha fliplr(alphaest-3*sigalpha)],[0.5 0 0]);set(p,'EdgeColor','none');
p=patch([alpha fliplr(alpha)],[alphathe+3*sigalpha fliplr(alphathe-3*sigalpha)],[0.5 0.5 0.5]);set(p,'EdgeColor','none','FaceAlpha',0.3);
%area(alpha,[alphaest'+3*sigalpha'],'FaceColor',[0.5 0 0],'EdgeColor','None');
%area(alpha,[alphaest'-3*sigalpha'],'FaceColor',[0.99 0.99 0.99],'EdgeColor','None');
%area(alpha,[alphathe'+3*sigalpha'],'FaceColor',[0.5 0.5 0.5],'EdgeColor','None');
%area(alpha,[alphathe'-3*sigalpha'],'FaceColor',[1 1 1],'EdgeColor','None');
set(gca,'Layer','Top');
plot(alpha,alpha,'--k',alpha,alphaest,'b',alpha,alphathe,'-.b');hold off
legend('Expected','Bootstrap','Gaussian','+- 3 Error bars Bootstrap','+- 3 Error bars Gaussian','Location','Best');
xlabel('Type I error \alpha');ylabel('Observed FPF');
title('Specificity calibration under H_0: homogeneous variance');drawnow

% Chi-square test of variance homogeneity (H1).
for isig=1:length(sigma)
 for rep=1:MC
  x = fft(randn(n,1))/sqrt(n);
  pos = ceil((n/2-3)*rand(1)+2);
  x(pos) = sigma(isig)*x(pos);
  per = abs(x).^2;
  logper = log(per(2:end/2-1));
  crithd(rep) = sum((logper - mean(logper)).^2)/psi(1,1); 
  xb   = std(x)*complex(randn(n,B),randn(n,B))/sqrt(2);
  perb = abs(xb).^2;
  outper = log(perb(2:end/2-1,:));
  crithdb(rep,:) = sum((outper - repmat(mean(outper),length(logper),1)).^2)/psi(1,1);
%   for b=1:B
%     xb   = std(x)*complex(randn(n,1),randn(n,1))/sqrt(2);
%     perb = abs(xb).^2;
%     outper = log(perb(2:end/2-1));
%     crithdb(rep,b) = sum((outper - mean(outper)).^2)/psi(1,1);
%   end
  %outper = logper(ceil(length(logper)*rand(length(logper),B)));
  %crithdb(rep,:) = sum((outper - repmat(mean(outper),length(logper),1)).^2)/var(logper);
 end
 crithdb=sort(crithdb')';
 q=ceil((B+1)*(1-0.05));
 crtl=chi2inv(1-0.05,length(logper)-1);
 alphatheH1(isig)  = length(find(crithd'>=crtl))/MC;
 alphaestH1(isig)  = length(find(crithd'>=crithdb(:,q)))/MC;
end

subplot(122) 
V=log(sigma.^2).^2/psi(1,1)*(length(logper)-1)/length(logper);
plot(log2(sigma),1-ncx2cdf(crtl,length(logper)-1,V),'b',log2(sigma),alphaestH1,'b*',log2(sigma),alphatheH1,'b^');grid on;axis tight;
legend('Theoretical','Bootstrap','Gaussian','Location','Best');
xlabel('log_2(\sigma)');ylabel('Observed TPF');
title('Sensitivity calibration under H_1 (\alpha=0.05): heteroscedastic profile at one sample');

