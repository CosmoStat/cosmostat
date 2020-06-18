;+
; NAME:
;        mrs_allstat
;
; PURPOSE:
;        Return statistical information relative to several multiscales transforms 
;        of a given data set. The six used multiscale transforms are: the pyramidal 
;        wavelet transform, the isotropic undecimated wavelet transform, the ridgelet 
;		 transform with a block size equals to 8, the ridgelet transform with a block
;		 size equals to 16, the ridgelet transform with a block size equals to 32 and
;		 the pyramidal curvelet transform.
;
;		 All statistical informations are computed with the survival option.
;
;        The return value is a IDL structure with several fields for six transforms statistics: 
;
;        Each of the statistics fields is a 2D IDL array of 9 elements times Number of scales and 
;        For each scale j, we have:
;			Tab[0,j] = standard deviation
;			Tab[1,j] = skewness
;			Tab[2,j] = Kurtosis
;			Tab[3,j] = Min
;			Tab[4,j] = Max 
;			Tab[5,j] = HC
;			Tab[6,j] = HC^+
;			Tab[7,j] = Cumulant of order 5
;			Tab[8,j] = Cumulant of order 6
;
;			If TabFile is set, then the statistic is computed on a set of images. Tab[*,*,f] will be the statistic related to the file TabFile[f]
;	CALLING:
;
;		StatData = mrs_allstat( Imag, TabFile=TabFile, NbrScale2D=NbrScale2D, TabStatName=TabStatName, verb=verb, 
;								normMad=normMad, iwt=iwt, owt=owt, rid8=rid8, rid16=rid16, rid32=rid32, cur=cur, 
;								all=all, save=save, TabTransformName=TabTransformName )
;
;	RETURN VALUE
;		StatData -- IDL structure with the following fields:
;				OWT : IDL float array [ 9, NbrScale ], statistics of the Orthogonal Wavelet Transform.
;				OWTSurv : TabSurvStat parameter for Orthogonal Wavelet Transform.
;				OWTSurvNu : TabSurvNu parameter for Orthogonal Wavelet Transform.
;				IWT : IDL float array [ 9, NbrScale ], statistics of the Isotropic Undecimated Wavelet Transform.
;				IWTSurv : TabSurvStat parameter for Isotropic Undecimated Wavelet Transform.
;				IWTSurvNu : TabSurvNu parameter for Isotropic Undecimated Wavelet Transform.
;				Rid8 : IDL float array [ 9, NbrScale ], statistics of the Ridgelet Transform (Length=8).
;				Rid8Surv : TabSurvStat parameter for Ridgelet Transform (Length=8).
;				Rid8SurvNu : TabSurvNu parameter for Ridgelet Transform (Length=8).
;				Rid16 : IDL float array [ 9, NbrScale ], statistics of the Ridgelet Transform (Length=16).
;				Rid16Surv : TabSurvStat parameter for Ridgelet Transform (Length=16).
;				Rid16SurvNu : TabSurvNu parameter for Ridgelet Transform (Length=16).
;				Rid32 : IDL float array [ 9, NbrScale ], statistics of the Ridgelet Transform (Length=32).
;				Rid32Surv : TabSurvStat parameter for Ridgelet Transform (Length=32).
;				Rid32SurvNu : TabSurvNu parameter for Ridgelet Transform (Length=32).
;				Cur : IDL float array [ 9, NbrScale ], statistics of the Curvelet Transform.
;				CurSurv : TabSurvStat parameter for Curvelet Transform.
;				CurSurvNu : TabSurvNu parameter for Curvelet Transform.
;				TabStatName : IDL table of string: TabStatName = [ "Sigma", "Skewness", "Kurtosis", "Min", "Max", "HC1", "HC2", "CUMULANT ORDER 5", "CUMULANT ORDER 6" ]
;				TabTransformName : IDL table of string: TabTransformName=[ 'Orthogonal Wavelet', 'Isotropic Undecimated Wavelet', 'Ridgelet Transform (Length=8)', 
;																			'Ridgelet Transform (Length=16)', 'Ridgelet Transform (Length=32)', 'Curvelet' ]
;
;				If a transform is not chosen, the three corresponding parameters fields are set to 0.
;       
; INPUTS:
;     Data -- IDL array of healpix map: Input data to analyze
;
; INPUT KEYWORDS:
;		verb: scalar -- if set, the calculated statistics are printed on the screen
;		NbrScale2D  -- int: number of scales. Default is 4
;		NormMad -- scalar, if set, a normalization is applied to the ridgelet and curvelet coefficients.
;		iwt -- if set, the statistics of the Isotropic Undecimated Wavelet Transform are computed
;		owt -- if set, the statistics of the Orthogonal Wavelet Transform are computed
;		rid8 -- if set, the statistics of the Ridgelet Transform (Length=8) are computed
;		rid16 -- if set, the statistics of the Ridgelet Transform (Length=16) are computed
;		rid32 -- if set, the statistics of the Ridgelet Transform (Length=32) are computed
;		cur -- if set, the statistics of the Curvelet Transform are computed
;		all -- if set, the statistics of all the 6 transforms are computed
;		save -- if set, the results are saved in separate files
;		TabFile -- IDL table of string, list of file where the function read the maps to be analized, in that case, on output, 
;					Imag is the last map that had been proceed and the fields in the return value structure are 3D arrays
;
; OUTPUT KEYWORDS: 
;		TabStatName -- IDL table of string: TabStatName = [ "Sigma", "Skewness", "Kurtosis", "Min", "Max", "HC1", "HC2", "CUMULANT ORDER 5", "CUMULANT ORDER 6" ]
;		TabTransformName -- IDL table of string: TabTransformName=[ 'Orthogonal Wavelet', 'Isotropic Undecimated Wavelet', 'Ridgelet Transform (Length=8)', 
;																	'Ridgelet Transform (Length=16)', 'Ridgelet Transform (Length=32)', 'Curvelet' ]
;
; EXAMPLE:
;
;       TabStat = mrs_allstat( Data, NbrScale2D=4, /verb, /all)
;		Compute the six transforms with 4 scales and compute statistical information relative to each scale and transform.
;
; EXTERNAL CALLS:
;       mrs_ridstat, mrs_curstat, mrs_wtstat
;
; HISTORY:
;	Written: Jean-Luc Starck, 2005
;	September, 2005 File creation
;-

function mrs_allstat, Imag, TabFile=TabFile, NbrScale2D=NbrScale2D, TabStatName=TabStatName,verb=verb, Ntrans=Ntrans, normMad=normMad, iwt=iwt, owt=owt, rid8=rid8, rid16=rid16, rid32=rid32, cur=cur, all=all, save=save, TabTransformName= TabTransformName, Mask=Mask

if N_PARAMS() LT 1  and not keyword_set(TabFile) then begin 
print, 'CALLING SEQUENCE: TabStat = mrs_allstat(Imag, TabFile=TabFile, NbrScale2D=NbrScale2D, TabStatName=TabStatName, verb=verb, Ntrans=Ntrans, normMad=normMad, iwt=iwt, owt=owt, rid8=rid8, rid16=rid16, rid32=rid32, cur=cur, all=all, save=save, TabTransformName=TabTransformName)'
        Stat=-1
	goto, DONE
        end

Ntrans=5
if not keyword_set(TabFile) then TabFile=0
TabTransformName=['Orthogonal Wavelet', 'Isotropic Undecimated Wavelet', 'Ridgelet Transform (Length=8)', $
                   'Ridgelet Transform (Length=16)', 'Ridgelet Transform (Length=32)', 'Curvelet']
IWTStat = 0
IWTTabSurvStat =0
 IWTTabSurvNu = 0
 RidStat8 = 0
 Rid8TabSurvStat = 0
 Rid8TabSurvNu = 0
 RidStat16 = 0
 Rid16TabSurvStat = 0
 Rid16TabSurvNu = 0
 RidStat32 = 0
 Rid32TabSurvStat = 0
 Rid32TabSurvNu = 0
 CurStat = 0
 CurTabSurvStat = 0
 CurTabSurvNu = 0
 OWTTabSurvNu=0
 OWTStat=0  
 OWTTabSurvStat=0
 StatLatWT=0

if not keyword_set(NbrScale2D) then NbrScale2D=4

if keyword_set(verb) then print, 'STEP 1: Isotropic Wavelet: NbrScale=', NbrScale2D
if keyword_set(iwt) or keyword_set(all) then begin
   TabSurvNu=0
   TabSurvStat=0
   TabMaskSurvStat=0
   TabMaskAllSurvStat=0
   TabMaskAllStat=0
   TabMaskStat =0 
   IWTStat = mrs_wtstat(Imag, tabfile=TabFile, NbrScale=NbrScale2D, wt=wt, verb=verb, /survival, TabSurvStat=IWTTabSurvStat, TabSurvNu=IWTTabSurvNu, TabStatName=TabStatName,  $
                        Mask=Mask, TabMaskStat=TabMaskStat, TabMaskSurvStat=TabMaskSurvStat, TabMaskAllSurvStat=TabMaskAllSurvStat, TabMaskAllStat=TabMaskAllStat, StatWT=StatLatWT)
end
if keyword_set(TabFile) then Imag = 0

if keyword_set(verb) then print, '`STEP 2: Orthogonal Wavelet: NbrScale=', NbrScale2D
if keyword_set(owt) or keyword_set(all) then begin
   OWTStat = mrs_owtstat(Imag, tabfile=TabFile, NbrScale=NbrScale2D, verb=verb, /survival, TabSurvStat=OWTTabSurvStat, TabSurvNu=OWTTabSurvNu, TabStatName=TabStatName)
   if keyword_set(save) then save, filename='xx_OWTStat.xdr', OWTStat
end
if keyword_set(TabFile) then Imag = 0


if keyword_set(verb) then print, 'STEP 3: Ridgelet: blocksize=8'
NbrScaleRid8=3
if keyword_set(rid8) or keyword_set(all) then begin
   RidStat8 = mrs_ridstat(Imag, tabfile=TabFile,  normMad=normMad, NbrScale=NbrScaleRid8,blocksize=8,verb=verb,/survival, TabSurvStat=Rid8TabSurvStat, TabSurvNu=Rid8TabSurvNu, TabStatName=TabStatName)
   if keyword_set(save) then save, filename='xx_Rid8Stat.xdr',  RidStat8
end
if keyword_set(TabFile) then Imag = 0

if keyword_set(verb) then print, 'STEP 4: Ridgelet: blocksize=16'
NbrScaleRid16 = 3
if keyword_set(rid16) or keyword_set(all) then begin
RidStat16 = mrs_ridstat(Imag, tabfile=TabFile,  normMad=normMad, NbrScale=NbrScaleRid16,blocksize=16,verb=verb,/survival, TabSurvStat=Rid16TabSurvStat, TabSurvNu=Rid16TabSurvNu, TabStatName=TabStatName)
   if keyword_set(save) then save, filename='xx_Rid16Stat.xdr',  RidStat16
end
if keyword_set(TabFile) then Imag = 0

if keyword_set(verb) then print, 'STEP 5: Ridgelet: blocksize=32'
NbrScaleRid32=4
if keyword_set(rid32) or keyword_set(all) then begin
   RidStat32 = mrs_ridstat(Imag,  tabfile=TabFile,  normMad=normMad, NbrScale=NbrScaleRid32,blocksize=32, /survival,verb=verb, TabSurvStat=Rid32TabSurvStat, TabSurvNu=Rid32TabSurvNu, TabStatName=TabStatName)
   if keyword_set(save) then save, filename='xx_Rid32Stat.xdr',  RidStat32
end
if keyword_set(TabFile) then Imag = 0

if keyword_set(verb) then print, 'STEP 6: Curvelet: NbrScale=',NbrScale2D
if keyword_set(cur) or keyword_set(all) then begin
  CurStat = mrs_curstat(Imag, tabfile=TabFile, NbrScale=NbrScale2D,verb=verb, normMad=normMad, /survival, TabSurvStat=CurTabSurvStat, TabSurvNu=CurTabSurvNu, TabStatName=TabStatName)
   if keyword_set(save) then save, filename='xx_CurStat.xdr',  CurStat
end

if keyword_set(TabFile) then Imag = 0

StatData = {OWT:OWTStat, OWTSurv:OWTTabSurvStat , OWTSurvNu: OWTTabSurvNu,   $
	        IWT:IWTStat, IWTSurv:IWTTabSurvStat , IWTSurvNu: IWTTabSurvNu,   $
            Rid8:RidStat8, Rid8Surv:Rid8TabSurvStat, Rid8SurvNu: Rid8TabSurvNu, $
            Rid16:RidStat16, Rid16Surv:Rid16TabSurvStat, Rid16SurvNu: Rid16TabSurvNu, $
            Rid32:RidStat32, Rid32Surv:Rid32TabSurvStat, Rid32SurvNu: Rid32TabSurvNu, $
            Cur:CurStat, CurSurv:CurTabSurvStat, CurSurvNu: CurTabSurvNu, $
            TabStatName:TabStatName, TabTransformName:TabTransformName, TabFile: TabFile, IWTLat:StatLatWT} 

DONE:
   return, StatData

END

 
;=======================================================================
;==============================================================================================
; Assuming a random data set, with sorted values,
; the function get_pvalue return the probability that
; that a given value is the realization of a variable
; which follows the distribution described by the 
; random data set

function get_pvalue, TabnormSort, Val
vs = size(TabnormSort)
n = vs[1]
Tind = where (TabnormSort GT Val, c)
if c eq n then p = 1. $
else if c eq 0 then p = 0. $
else p = 1. - double(Tind[0]) / double(n) 
return, p
end

;=======================================================

function mrs_cmp_onestat_allband, TabStatSimu, StatData, UseStat=UseStat, Band=Band, verb=verb, Nsigma=NSigma, Pval=PVal, HistoPVal=HistoPVal, TransName=TransName
TabStatName = ["Sigma", "Skewness", "Kurtosis", "Min", "Max", "HC1", "HC2", "CUMULANT Order 5", "CUMULANT Order 6"]
 
 if not keyword_set(Band) then Band=0
 B=Band

if not keyword_set(UseStat) then UseStat=2
s=UseStat

vs = size(StatData)
NbrBand = vs[2]
TabVal = dblarr(NbrBand)
TabNSigma = dblarr(NbrBand)
TabPVal = dblarr(NbrBand)
TabHistoNSigma = dblarr(NbrBand)
TabHistoPVal = dblarr(NbrBand)
TabMeanSimu = dblarr(NbrBand)
TabSigmaSimu = dblarr(NbrBand)

for b=0,NbrBand-1 do BEGIN
    Simu = reform( TabStatSimu[UseStat, b, *] )
    SigmaSimu = sigma(Simu)
    MeanSimu = mean(Simu)
    TabMeanSimu[b] = MeanSimu
    TabSigmaSimu[b] = SigmaSimu
        
    Data = StatData[UseStat, b]
    NSigma = ABS(Data - MeanSimu) / SigmaSimu
    Pval = double(1.) - (1. + errorf(double(NSigma/sqrt(2.)))) / 2.
    inds = sort(Simu)
    DataH0Sort = Simu(inds)
    vs = size(DataH0Sort)
    n = vs[1]
    p = DataH0Sort
    Pval2 = get_pvalue(DataH0Sort, Data)
    HistoPVal = Pval2
    HistoNSig = getnsig(Pval2)
    
    TabVal [b] = Data
    TabNSigma[b] = NSigma
    TabPVal[b] = Pval
    TabHistoNSigma[b] = HistoNSig
    TabHistoPVal[b] = HistoPval

    if keyword_set(verb) then print, "Band ", b+1, ' ', TabStatName[s], ' NSigma = ', NSigma, ' PVal = ', Pval, ' HistoPVal = ', Pval2, ' HistoNsigma = ', HistoNSig
end

CMP = {TabSigmaSimu:TabSigmaSimu, TabMeanSimu:TabmeanSimu, TabCoefData:TabVal, TabNSigma:TabNSigma, $
       TabPVal:TabPVal, TabHistoNSigma:TabHistoNSigma,  TabHistoPVal:TabHistoPVal, StatName:TabStatName[s], StatIndex:s, BandNumber: b} 
return, CMP
end

;==============================================================================================

pro testcmp, TabF, c, sc, s
Nsim=10
Nside=64
Nbrscale=3
TabF = strarr(Nsim)
for s=0, Nsim-1 do begin
 ; c = mrs_getcmb(nside=nside)
 FN = 'xx_test_' + strc(s+1) + '.fits'
 TabF[s] = FN
 ; mrs_write, FN, c
end

C = getcmb(nside=nside)
sc = mrs_owtstat(c, Nbrscale=Nbrscale, /survival, TabSurvStat=OWTSurvStat, TabSurvNu=TabSurvNu1)
s = mrs_owtstat(tabfile=TabF, Nbrscale=Nbrscale, /survival, TabSurvNu=OWTTabSurvNu, TabAllSurvStat=OWTTabAllSurvStat)
; CMB = mrs_cmp_one_band(s, sc, b=0, /verb)

sc = mrs_curstat(c, Nbrscale=Nbrscale, /survival, TabSurvStat=CurTabSurvStat, TabSurvNu=CurTabSurvNu)
 s = mrs_curstat(tabfile=TabF, Nbrscale=Nbrscale, /survival,  TabSurvNu=TabSurvNu, TabAllSurvStat=CurTabAllSurvStat)
 
 sc = mrs_wtstat(c, Nbrscale=Nbrscale, /survival, TabSurvStat=CurTabSurvStat, TabSurvNu=CurTabSurvNu)
 s = mrs_wtstat(tabfile=TabF, Nbrscale=Nbrscale, /survival,  TabSurvNu=TabSurvNu, TabAllSurvStat=CurTabAllSurvStat)
 
; CMB = mrs_cmp_one_band(s, sc, b=0, /verb)


 sc = mrs_ridstat(c, Nbrscale=Nbrscale, /survival, TabSurvStat=CurTabSurvStat, TabSurvNu=CurTabSurvNu)
 s = mrs_ridstat(tabfile=TabF, Nbrscale=Nbrscale, /survival,  TabSurvNu=TabSurvNu, TabAllSurvStat=CurTabAllSurvStat)

; CMB = mrs_cmp_one_band(s, sc, b=0, /verb)

; plotstat, sc,   UseStat=2, Title=Title, oplot=oplot, line=line, thick=thick, psym=psym
end


;==============================================================================================

pro mrs_statmeansigma, Tab, TabMean, TabSigma
vs = size(tab)
Nstat = vs[1]
Nb = vs[2]
Nsim = vs[3]
TabStatMean=fltarr(Nb,2)
TabStatSigma=fltarr(Nb,2)

if not keyword_set(usestat) then usestat=2

end

;==============================================================================================

pro plotstat, TabStatSimu, PS=PS, UseStat=UseStat, Title=Title, oplot=oplot, line=line, thick=thick, psym=psym, log=log
TabStatName = ["Sigma", "Skewness", "Kurtosis", "Min", "Max", "HC1", "HC2", "CUMULANT Order 5", "CUMULANT Order 6"]

; CMP = {SigmaSimu: sigmaSimu, MeanSimu:  meanSimu, StatVal: Data, StatName: TabStatName[s], StatIndex: s, BandNumber: b} 
vs = size(TabStatSimu)
NStat = vs[1]
NbrBand = vs[2]
Nj = NbrBand
NbrSimu = vs[2]
if not keyword_set(useStat) then useStat = 0
if not keyword_set(Title) then Title = TabStatName[useStat]

TabSigmaSimu = dblarr(Nj)
TabMeanSimu = dblarr(Nj)

for b=0,NbrBand-1 do BEGIN
    if keyword_set(DataStat) then ValDat = DataStat[UseStat, b]
    Simu = reform( TabStatSimu[UseStat, b, *] )
    if keyword_set(log) then Simu=alog(simu)
    ; Data = StatData[UseStat, b]
    TabSigmaSimu[b] = sigma(Simu)
    TabMeanSimu[b] = mean(Simu)
    if keyword_set(DataStat) then TabMeanSimu[b] = ValDat - TabMeanSimu[b]
END
    if not keyword_set(oplot) then ploterror, title=title, ytitle=TabStatName[useStat], xtitle='Scale', indgen(Nj) +1,  TabmeanSimu, TabSigmaSimu, line=line, thick=thick, psym=psym $
    else oploterror, indgen(Nj) +1,  TabmeanSimu, TabSigmaSimu, line=line, thick=thick, psym=psym
end

;==============================================================================================

pro plotallstat, AllStat, iwt=iwt, owt=owt, rid8=rid8, rid16=rid16, rid32=rid32, cur=cur, all=all, Title=Title, oplot=oplot, ps=ps
TabTransformName=['Orthogonal Wavelet', 'Isotropic Undecimated Wavelet', 'Ridgelet Transform (Length=8)', $
                   'Ridgelet Transform (Length=16)', 'Ridgelet Transform (Length=32)', 'Curvelet']
; TabStatName = ["Sigma", "Skewness", "Kurtosis", "Min", "Max", "HC1", "HC2", "CUMULANT Order 5", "CUMULANT Order 6"]

; TabTransformName = AllStat.TabTransformName
vs = size(TabTransformName)
Ntrans = vs[1]
TabTransOK = intarr(Ntrans)
if  keyword_set(all) then TabTransOK[*] = 1
if  keyword_set(owt) then TabTransOK[0] = 1
if  keyword_set(iwt) then TabTransOK[1] = 1
if  keyword_set(rid8) then TabTransOK[2] = 1
if  keyword_set(rid16) then TabTransOK[3] = 1
if  keyword_set(rid32) then TabTransOK[4] = 1
if  keyword_set(cur) then TabTransOK[5] = 1

for t=0, Ntrans-1 do begin
 Stat = AllStat.(t*3)
 vs = size(Stat)
 if vs[0] NE 0 then begin
    NameTrans = TabTransformName[t]
    NStat = (size(Stat))[1]
    for s=0, NStat-1 do begin
      NameStat =  AllStat.TabStatName[s]
      Tit = NameTrans + ': ' + NameStat
      if keyword_set(Title) then Tit = Title + ': ' + Tit
      plotstat, Stat, UseStat=s, Title=Tit, oplot=oplot, line=line, thick=thick, psym=psym
      wait, 3
    end
 end
end
end

;==============================================================================================

pro tt, T1, T2, usestat=usestat, title=title, log=log
plotstat, T1, us=usestat, title=title, thick=2, log=log
plotstat, T2, us=usestat, line=2, thick=2, /op, log=log
end

 
 
 

;==============================================================================================
