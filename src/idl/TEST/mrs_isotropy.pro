;+
; NAME:
;        mrs_isotropy
;
; PURPOSE:
;        This function tests if the input image is isotropic. It return 1 if it is isotrop and 0 otherwise.
;        If the keyword Alm is set, it takes as an input the Alm coefficients instead of the healpix map
;
; CALLING:
;     Isotrop = mrs_isotropy(Imag, Cl=Cl)
;     Isotrop = mrs_isotropy(Alm=Alm, Cl=Cl)
;    
; INPUT:
;     Imag -- IDL array of healpix map: Input image be tested 
;
; INPUT/OUTPUT:
;     Alm -- IDL structure: Alm coefficients (see mrs_almtrans.pro)  
;
; INPUT KEYWORDS:
;     
; INPUT/OUTPUT KEYWORDS:
;       
;
; EXTERNAL CALLS:
;       
;
; EXAMPLE:
;        
;         
; HISTORY:
;       Written: Jalal Fadili & Jean-Luc Starck, April 15, 2010.
;-

 
;==================================================================

function chi2_variance_homogeneity_test, Signal, Alpha, Chi2Score=Chi2Score, Pval=Pval
; Test of homogeneity of variance
; Suppose que le spectre de puissance est connue, et de variance 1.
; 
; IN Alpha: risque de 1er espece = prob d'avoir une fausse detection sour H0
; OUT  Chi2Score = somme des carres centres = score chi2
; PVal = si PVal < Alpha, ce n'est pas homogene ==> rejet H0 et retourne 0
; sinon retourne 1
N = N_ELEMENTS(Signal)
Mean  = total(Signal) / float(N)
Chi2Score = total(abs(Signal - Mean)^2.0)
tmp = igamma( (n-1.0)/2.0,Chi2Score/2.0, /double)
Pval = 2.0*min([tmp, 1.0-tmp])
if Pval LT Alpha then Detect = 0 else Detect = 1
return, Detect
end

;==================================================================
function sum, a

  sa = size(a,/dimensions)
  suma = dblarr(sa[0])
  for i=0,sa[0]-1 do suma[i] = total(a[i,*])
  return,suma

end

;==================================================================

function meanvec, a

  sa = size(a,/dimensions)
  suma = dblarr(sa[1])
  for i=0,sa[1]-1 do suma[i] = mean(a[*,i])
  return,suma

end
;==================================================================
; Approximation of the sort function of matlab
function sortvec, a
 help, a
  sa = size(a,/dimensions)
  suma = dblarr(sa[0],sa[1])
  for i=0,sa[0]-1 do suma[i,*] = a[i,sort(a[i,*])]
  return,suma

end

;==================================================================

function chi2_mean_cst_log_periodo_test, Signal, Alpha, Chi2Score=Chi2Score 
;  Mean constancy test of mean after the log on the periodogram using bootstap

x = Signal
N = N_ELEMENTS(Signal)
B  = 1e3

per = abs(x)^2
size_per = size(Signal,/dimensions)
logper = dblarr(size_per[0]/2-2)
for i=1L,(size_per[0]/2)-2 do logper[i-1] = alog(per[i])
crithd = total((logper - mean(logper))^2)/psi(1.0,1.0)
    
xb = stddev(real_part(x)+imaginary(x))* $
         complex(randomn(seed,n,B,/normal),randomn(seed,n,B,/normal))/sqrt(2.0)
perb = abs(xb)^2
;print, size_per
 outper = dblarr(size_per[0]/2-2,B)
 
  for i=1L,size_per[0]/2-2 do outper[i-1,*] = alog(perb[i,*])
 size_logp = size(outper,/dimensions)
 ;help, outper
Temp = meanvec(outper)
;help, Temp

 ; crithdb = sum(transpose((outper - repmat(meanvec(outper),size_logp[0],1))^2))   / psi(1.0,1.0)
crithdb = sum(transpose((outper)^2))   / psi(1.0,1.0)

i = sort(crithdb)
crithdbt = crithdb(i)

; help, crithdb
   
    q = ceil((B)*(1-alpha))-1
    Chi2Score = crithdbt[q]
    if crithd  GT Chi2Score  then Detect = 0 else Detect = 1
return, Detect
end
  
 ;==================================================================
 
function cfsigma, taba
return, total( abs(taba)^2) / float( N_ELEMENTS(TabA))
 end
 
 ;==================================================================

function mrs_isotropy, Imag,  Alm=Alm, Cl=Cl, Alpha=Alpha, Lmin=Lmin, Lmax=Lmax, Chi2Score=Chi2Score, Pval=Pval, NSigma=NSigma, Verb=Verb, TabChi2Score=TabChi2Score, TabPVal=TabPVal, TabDet=TabDet
; Test NULL Hypethosis that the process is isotropic
; smaller is alpha, the higher is the confidence level, i.e. more we can trust the result 

COMMON C_PLANCK
 
if N_PARAMS() LT 1 and not keyword_set(Alm) then begin 
        print, 'CALLING SEQUENCE: Isotrop =   mrs_isotropy(Imag,  Alm=Alm, Cl=Cl, Alpha=Alpha) '
        goto, DONE
        end
        
if not keyword_set (Alm) then  mrs_almtrans, Imag, Alm, /complex, /tab, /norm
if not keyword_set (Cl) then  Cl = mrs_alm2spec(Alm)
if not keyword_set (Lmin) then  Lmin = 10
if Lmin LT 5 then Lmin=5
if not keyword_set (LMAX) then  LMAX = ALM.LMAX
if not keyword_set (Alpha) then  Alpha = 0.05
if keyword_set (NSigma) then  Alpha = 1. - errorf(double(Nsigma) / sqrt(double(2.)))

A = ALM.Alm
Coef = sqrt(Cl[2:*])   
ind = where(Coef EQ 0, c)
if c GT 0 then Coef[ind] = 1.
Coef = 1. / Coef
B = A
for l=2, ALM.LMAX do A[l,*] = A[l,*] * Coef[l-2]
TabA = reform(A[LMin,0:Lmin])
for l= LMIN+1, LMAX do  TabA = [TabA, reform(A[l,0:l]) ]
 
 ; Global Test of homogeneity of variance, assumung power spectrum Cl is known and equal to 1 (because of previous normalization)
Isotrop  = chi2_variance_homogeneity_test(TabA, Alpha, Chi2Score=Chi2Score, Pval=Pval)

Chi2ScoreB = 0

; Global test of mean constancy on the log periodogram using bootstrap (Cl unknown)
IsotropB  = chi2_mean_cst_log_periodo_test(TabA, Alpha, Chi2Score=Chi2ScoreB)

TabChi2Score = dblarr(LMAX-LMIN+1)
TabPVal = dblarr(LMAX-LMIN+1)
TabDet = intarr(LMAX-LMIN+1)
TabSig = dblarr(LMAX-LMIN+1)

TabDetB = intarr(LMAX-LMIN+1)
TabChi2ScoreB = dblarr(LMAX-LMIN+1)

; Local test per l
for l= LMIN, LMAX do begin
   TabM = reform(A[l,0:l]) 
   TabSig[l-LMIN] =  total(abs(TabM)^2) / float(l+1)
   ImZ  = chi2_variance_homogeneity_test(TabM, Alpha, Chi2Score=C2, Pval=P2)
  TabChi2Score[l-LMIN] = C2
  TabPVal[l-LMIN] = P2
  TabDet[l-LMIN] = ImZ
  
 TabM = reform(B[l,0:l]) 
 IsotropBL  = chi2_mean_cst_log_periodo_test(TabM, Alpha, Chi2Score=Chi2ScoreBL)
 TabChi2ScoreB[l-LMIN] = Chi2ScoreBL
  TabDetB[l-LMIN] = IsotropBL

end
 
if keyword_set(Verb) then begin
	 if Isotrop  EQ 1 then print, " Global test 1: Isotropic map: Prob = ", Pval * 100.  $
	 else  print, " Global test 1: Anisotropic map: Prob = ", (1. -Pval) * 100.  
	 if IsotropB  EQ 1 then print, " Global test 2: Isotropic map." $
	 else  print, " Global test 2: Anisotropic map. "  
	 ind = where (TabDet EQ 0, c )
	 print, 'Local test 1 (variance homogeneity): ', c, ' of positive anisotropy detections though the l'
	 ind = where (TabDetB EQ 0, c )
	 print, 'Local test 2 (mean constancy on  log PowSpectrum Density): ', c, ' of positive anisotropy detections though the l'
    xr = lindgen(LMAX+1)
    plot, xr[LMIN:*], tabdet
    oplot, xr[LMIN:*], TabDetB, line=2
 end
 
  Ret = {Alpha: Alpha, $
         A:A, $
         TabA: TabA, $
         Cl: Cl, $
         LMin: LMin, $
         LMax: LMax, $
 	     Isotrop: Isotrop, $
        TabSig: TabSig, $
         IsotropB: IsotropB, $
         Chi2Score: Chi2Score, $
         Chi2ScoreB: Chi2ScoreB, $
         Pval: Pval, $
         TabChi2Score: TabChi2Score, $
         TabPVal: TabPVal, $
         TabDet: TabDet, $
          TabChi2ScoreB: TabChi2ScoreB, $
         TabDetB: TabDetB }
         
DONE:
return, Ret
END

;==================================================================

pro testiso
nside=128L
C = getcmb(Cl=Cl, nside=nside)
C = randomn(seed, 128L^2*12)
mrs_almtrans, C, Alm, /complex, /tab, /norm
A = Alm.Alm
; A(100, 10) =   25.
Alm.Alm = A
Cl = fltarr(385)+1.
s = mrs_isotropy(Alm=Alm, nsigma=5, TabChi2Score=TabChi2Score, TabPVal=TabPVal, TabDet=TabDet, lmax=128, /verb, cl=cl)
end

;==================================================================

pro testaniso
; Test du modele de Bianchi
;   Notre test 1 est comparable a http://fr.arxiv.org/PS_cache/arxiv/pdf/1002/1002.3173v2.pdf
;    ==> sensibilite au bon l (a 5sigma), mais a partir d'une amplitude de 5 fois le modele de Bianchi.
;           Sensibilite a 2sigma a 1xBianchi mais avec un FPF comparable a l'hypothese H0
;   ==> test de homogeneite de la variance plus sensible que le test de constance de la moyenne sur le log PSD
;   
;  ==> est que la carte inpaintee introduit des anisotropies detectable a 5sigma par ce test  ?
;  ==> Faire des tests plus sensibles, HC-like ?

nside=256L
B = rims("$ISAP/data/bianchi_ns256_bmap_scaled_microK_CMB.fits")
Cl=0
C = getcmb(Cl=Cl, nside=nside)
mrs_almtrans, C, AlmC, /complex, /tab, /norm
Cl = Cl * alm.normval^2
Pc = mrs_alm2spec(AlmC)
s = mrs_isotropy(Alm=AlmC, nsigma=5., TabChi2Score=CTabChi2Score, TabPVal=CTabPVal, TabDet=CTabDet, lmax=128, /verb, Cl=Cl)
s = mrs_isotropy(Alm=AlmC, nsigma=5., TabChi2Score=CTabChi2Score, TabPVal=CTabPVal, TabDet=CTabDet, lmax=128, /verb, Cl=Pc)


D = 2.*B  +C 
mrs_almtrans, D, Alm, /complex, /tab , /norm
s = mrs_isotropy(Alm=Alm, nsigma=5., TabChi2Score=TabChi2Score, TabPVal=TabPVal, TabDet=TabDet,  lmax=128, /verb, Cl=Cl)
s = mrs_isotropy(Alm=Alm, nsigma=5., TabChi2Score=TabChi2Score, TabPVal=TabPVal, TabDet=TabDet,  lmax=128, /verb, Cl=Pc)

Pd = mrs_alm2spec(Alm)

end

;==================================================================

pro test_planck_simu_iso

for i=2, 20 do begin
   FN = 'inp_masked_cmb_' + STRC(i)
   map = rims(FN+'.fits')
   print, i, ' ', FN
   mrs_almtrans, Map, Alm, /complex, /tab, /norm
   s = mrs_isotropy(Alm=Alm, nsigma=5, TabChi2Score=TabChi2Score, TabPVal=TabPVal, TabDet=TabDet, /verb)
   save, filename='iso_' + FN +'.xdr', s, TabPVal, TabDet, TabChi2Score
end

end

;==================================================================

pro test_planck_dx5_iso

Nf=11
TabM = ['CCA', 'FastMEM', 'N-GenILC', 'SMICA', 'AltICA', 'LGMCA', 'N-ILC', 'Commander', 'Sevem']
TabFN = ['inp_masked_CCA_I_CMB', $
'inp_masked_CMB_FastMEM_2048_5arcmin_v20110304_uK_common_mask', $
'inp_masked_CMB_needlet_GenILC_2048_5arcmin_mK_masked_dx5full', $
'inp_masked_cmb_smica_5arcmin_microK', $
'inp_masked_cmb_altica', $
'inp_Temp_DX5_lGMCA_mKAntenna', $
'inp_masked_nilc5arcminv2', $
'inp_masked_dx5_mean_CMB_ns2048', $
'inp_masked_sevem_fulldata', $
'inp_Temp_DX5_lGMCA_LAg_mKAntenna', $
'inp_masked_CMB_FastMEM_2048_5arcmin_v20110304_uK_common_mask']

for f=8,Nf-1 do begin
   Map = rims(TabFN[f]+'.fits')
   print, f, ' ', TabFN[f]
   mrs_almtrans, Map, Alm, /complex, /tab, /norm
   s = mrs_isotropy(Alm=Alm, nsigma=5, TabChi2Score=TabChi2Score, TabPVal=TabPVal, TabDet=TabDet, /verb)
   save, filename='iso_' + TabFN[f] +'.xdr', s, TabPVal, TabDet, TabChi2Score
end
 
end
;==================================================================
