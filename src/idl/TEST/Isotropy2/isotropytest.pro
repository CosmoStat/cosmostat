;
; Isotropy Test number 1
;
;  Written by Jalal Fadili for Matlab, translated into IDL by Francois-Xavier Dupe
;

; Approximation of the sum function of Matlab
function sum, a

  sa = size(a,/dimensions)
  suma = dblarr(sa[0])
  for i=0,sa[0]-1 do suma[i] = total(a[i,*])
  return,suma

end

; Approximation of the mean function of Matlab
function meanvec, a

  sa = size(a,/dimensions)
  suma = dblarr(sa[0])
  for i=0,sa[0]-1 do suma[i] = mean(a[i,*])
  return,suma

end

; Approximation of the sort function of matlab
function sortvec, a

  sa = size(a,/dimensions)
  suma = dblarr(sa[0],sa[1])
  for i=0,sa[0]-1 do suma[i,*] = a[i,sort(a[i,*])]
  return,suma

end

pro isotropytest

; Define some constants and vectors
  n  = 256
  MC = 1e3
  B  = 1e2
  sigma = 2.^linspace(-8.0,8.0,20)
  alpha = linspace(1.0/B,1.0-1.0/B,100)
  sigalpha = sqrt(alpha*sqrt(1.0-alpha)/MC)

; Initialize randomn
  tmp=randomn(seed,/normal)
 
;Chi-square test of variance homogeneity (H0)
  crithd = dblarr(MC)
  crithdb = dblarr(MC,B)
  p = dblarr(MC)
  for rep=0,MC-1 do begin

     x = n*fft( randomn(seed,n,/normal) )/sqrt(n)
     crithd[rep] = total(abs(x - mean(x))^2.0)
     tmp = igamma( (n-1.0)/2.0, crithd[rep]/2.0, /double)
     p[rep] = 2.0*min([tmp, 1.0-tmp])
     out = x[ceil(n*randomn(seed,n,B,/uniform))-1]
     crithdb[rep,*] = sum(transpose(abs(out - repmat(meanvec(out),n,1))^2.0)) $
                      /variance(real_part(x)+imaginary(x),/double)

  endfor

  size_a = size(alpha,/dimensions)
  crithdbt = sortvec(crithdb)
  alphaH0 = dblarr(size_a[0])
  alphaest = dblarr(size_a[0])
  for ia=0,size_a[0]-1 do begin

     q = ceil((B+1)*(1.0-alpha[ia]))-1
     ;crtl = chi2inv(1.0-alpha[ia],n-1)
     tmp1 = 0
     aix = where( p le alpha[ia], tmp1)
     alphaH0[ia] = tmp1/double(MC)
     aix = where( crithd ge crithdbt[*,q], tmp1 )
     alphaest[ia] = tmp1/double(MC)

  endfor

;plot(alpha,alpha,'--k',alpha,alphaest,'b',alpha,alphaH0,'-.b');hold off
;legend('Expected','Bootstrap','Gaussian','+- 3 Error bars Bootstrap','+- 3 Error bars Gaussian','Location','Best');
;xlabel('Type I error \alpha');ylabel('Observed FPF');
;title('Specificity calibration under H_0: homogeneous variance');drawnow

  window,0
  device, decomposed = 1
  plot,alpha,alpha
  oplot,alpha,alphaest,color=255 ; in red
  oplot,alpha,alphaH0,color=255*256 ; in blue

; Chi-square test of variance homogeneity (H1).
  size_sig = size(sigma,/dimensions)
  alphatheH1 = dblarr(size_sig[0])
  alphaestH1 = dblarr(size_sig[0])
  for isig=0,size_sig[0]-1 do begin
     for rep=0,MC-1 do begin

        x = n*fft( randomn(seed,n,/normal) )/sqrt(n)
        pos = ceil(n*randomn(seed,/uniform))-1
        x[pos] = sigma[isig]*x[pos]
        crithd[rep] = total(abs(x - mean(x))^2)
        out = x[ceil(n*randomn(seed,n,B,/uniform))-1]
        crithdb[rep,*] = sum(transpose(abs(out - repmat(meanvec(out),n,1))^2))$
                         /variance(real_part(x)+imaginary(x),/double)

     endfor

     crithdbt = sortvec(crithdb)
     q = ceil((B+1)*(1-0.05))-1
     crtl = chi2inv(1.0-0.05,n-1)
     tmp1 = 0
     aix = where( crithd ge crtl, tmp1 )
     alphatheH1[isig] = tmp1/double(MC)
     aix = where( crithd ge crithdbt[*,q], tmp1 )
     alphaestH1[isig] = tmp1/double(MC)
     
  endfor

;subplot(122)
;plot(log2(sigma),alphaestH1,'b',log2(sigma),alphatheH1,':b');grid on;axis tight;
;legend('Bootstrap','Gaussian');
;xlabel('log_2(\sigma)');ylabel('Observed TPF');
;title('Sensitivity calibration under H_1 (\alpha=0.05): heteroscedastic profile at one sample');

  window,1
  device, decomposed = 1
  plot,alog(sigma),alphaestH1
  oplot,alog(sigma),alphatheH1,color=255

end
