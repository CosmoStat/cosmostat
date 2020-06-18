;
; Isotropy Test number 2
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

pro isotropytest2

; Define some constants and vectors
 n  = 256
 MC = 1e3
 B  = 1e2
 sigma = 2.0^linspace(0,12,10)
 alpha = linspace(1.0/B,1.0-1./B,100)
 sigalpha = sqrt(alpha*sqrt(1.0-alpha)/MC)

; UMPI Mean constancy test of mean after the log on the periodogram
 crithd = dblarr(MC)
 crithdb = dblarr(MC,B)
 size_logp = 0
 for rep=0,MC-1 do begin

    x = fft( n*randomn(seed,n,/normal) )/sqrt(n)
    per = abs(x)^2
    size_per = size(per,/dimensions)
    logper = dblarr(size_per[0]/2-2)
    for i=1,(size_per[0]/2)-2 do logper[i-1] = alog(per[i])
    crithd[rep] = total((logper - mean(logper))^2)/psi(1.0,1.0)
    xb = stddev(real_part(x)+imaginary(x))* $
         complex(randomn(seed,n,B,/normal),randomn(seed,n,B,/normal))/sqrt(2.0)
    perb = abs(xb)^2
    outper = dblarr(size_per[0]/2-2,B)
    for i=1,size_per[0]/2-2 do outper[i-1,*] = alog(perb[i,*])
    size_logp = size(logper,/dimensions)
    crithdb[rep,*] = sum(transpose((outper - repmat(meanvec(outper),size_logp[0],1))^2)) $
                     /psi(1.0,1.0)

 endfor

 size_a = size(alpha,/dimensions)
 crithdbt=sortvec(crithdb)
 alphathe = dblarr(size_a[0])
 alphaest = dblarr(size_a[0])
 for ia=0,size_a[0]-1 do begin

    q = ceil((B+1)*(1-alpha[ia]))-1
    crtl = chi2inv(1.0-alpha[ia],size_logp[0]-1.0)
    tmp1 = 0
    aix = where( crithd ge crtl, tmp1 )
    alphathe[ia] = tmp1/double(MC)
    aix = where( crithd ge crithdbt[*,q], tmp1 )
    alphaest[ia] = tmp1/double(MC)

 endfor

 window,0
 device, decomposed = 1
 plot,alpha,alpha
 oplot,alpha,alphaest,color=255     ; in red
 oplot,alpha,alphathe,color=255*256 ; in blue

; Chi-square test of variance homogeneity (H1)
 size_sig = size(sigma,/dimensions)
 alphatheH1 = dblarr(size_sig[0])
 alphaestH1 = dblarr(size_sig[0])
 for isig=0,size_sig[0]-1 do begin
    for rep=0,MC-1 do begin

       x = fft( n*randomn(seed,n,/normal) )/sqrt(n)
       pos = ceil((n/2.0-3.0)*randomn(seed,/uniform)+2.0)-1
       x[pos] = sigma[isig]*x[pos]
       per = abs(x)^2
       size_per = size(per,/dimensions)
       logper = dblarr((size_per[0]-1)/2-2)
       for i=1,(size_per[0]-1)/2-2 do logper[i-1] = alog(per[i])
       crithd[rep] = total((logper - mean(logper))^2)/psi(1,1)
       xb = stddev(x)*complex(randomn(seed,n,B,/normal),randomn(seed,n,B,/normal))/sqrt(2)
       perb = abs(xb)^2
       outper = dblarr((size_per[0]-1)/2-2,B)
       for i=1,(size_per[0]-1)/2-2 do outper[i-1,*] = alog(perb[i,*])
       size_logp = size(logper,/dimensions)
       crithdb[rep,*] = sum(transpose((outper - repmat(meanvec(outper),size_logp[0],1))^2)) $
                        /psi(1,1)
       
    endfor

    q = ceil((B+1)*(1-0.05))-1
    crithdbt = sortvec(crithdb)
    crtl = chi2inv(1-0.05,size_logp[0]-1)
    tmp1 = 0
    aix = where( crithd ge crtl, tmp1 )
    alphatheH1[isig] = tmp1/MC
    aix = where( crithd ge crithdbt[*,q], tmp1 )
    alphaestH1[isig] = tmp1/MC 

 endfor

 window,1
 device, decomposed = 1
 plot,alog(sigma),alphaestH1
 oplot,alog(sigma),alphatheH1,color=255 ; in red

end
