;+
; NAME:
;        mrs_ridstat
;
; PURPOSE:
;	 Return statistical information relative to the ridgelet transform of a given data set.
;        The return value is a 2D IDL array of 9 elements x Number of scales.
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
;			If the keyword normmad is set, the ridgelet coefficients are first normalized.
;
; CALLING:
;
;	TabStat = mrs_ridstat( Data, TabStatName=TabStatName, NbrScale=NbrScale, BlockSize=BlockSize, NormMad=NormMad, verb=verb, 
;							Ridtrans=Ridtrans, survival=survival, TabFile=TabFile, TabSurvStat=TabSurvStat, TabAllSurvStat=TabAllSurvStat, TabSurvNu=TabSurvNu ) 
;       
; INPUTS:
;     Data -- IDL array of healpix map: Input data to analyze
;
; INPUT KEYWORDS:
;		verb: scalar -- if set, the calculated statists are printed on the screen.
;		BlockSize     --  LONG: Block size used in the ridgelet transform.
;		NbrScale  -- int: number of scales. By default, it is automatically estimated.
;		NormMad -- int: if set, a normalization is applied to the ridgelet coefficient.
;		TabFile -- IDL table of string, list of file where the function read the maps to be analized, in that case, on output, Imag is the last map that had been proceed and the return value is a 3D array:
;			Tab[i,j,f]	statistic i for scale j and map TabFile[f]
;		survival: if set, use the survival function and activate the keywords parameters TabSurvStat, TabAllSurvStat and TabSurvNu for the results
;
; OUTPUT KEYWORDS: 
;		TabStatName -- IDL table of string: TabStatName = ["Sigma", "Skewness", "Kurtosis", "Min", "Max", "HC1", "HC2", "CUMULANT ORDER 5", "CUMULANT ORDER 6"]
;		Ridtrans: -- IDL structure containing the ridgelet transform of Data.
;		TabSurvStat -- 2D float array [*,j] survival value at scale j
;		TabAllSurvStat -- 3D double array, use together with TabFile keyword, TabAllSurvStat[*,*,f] is TabSurvStat parameter for map TabFile[f]
;		TabSurvNu -- 2D float array [*,j] nu survival value at scale j
;
; EXAMPLE:
;     Compute the pyramidal wavelet transform with 4 scales
;       TabStat = mrs_ridstat(Data, NbrScale=4, /verb)
;
; EXTERNAL CALLS:
;       mrs_ridtrans, mrs_ridget, get_stat
;
; HISTORY:
;	Written: Jean-Luc Starck, 2005
;	September, 2005 File creation
;--------------------------------------------------------------------------------------------------------------------------------------------------------

function mrs_ridstat, Imag, NbrScale=NbrScale, BlockSize=BlockSize, TabStatName=TabStatName, verb=verb, Ridtrans=Ridtrans, NormMad=NormMad, $
      survival=survival, TabFile=TabFile, TabSurvStat=TabSurvStat, TabAllSurvStat=TabAllSurvStat, TabSurvNu=TabSurvNu


if N_PARAMS() LT 1  and not keyword_set(TabFile) then begin 
        print, 'CALLING SEQUENCE:  TabStat = mrs_ridstat( Data, TabStatName=TabStatName, NbrScale=NbrScale, BlockSize=BlockSize, NormMad=NormMad, verb=verb, Ridtrans=Ridtrans, survival=survival, TabFile=TabFile, TabSurvStat=TabSurvStat, TabAllSurvStat=TabAllSurvStat, TabSurvNu=TabSurvNu )'
        goto, DONE
        end

if keyword_set(TabFile) then NbrFiles = (size(TabFile))[1]

if keyword_set(TabFile) then Imag = mrs_read(TabFile[0]) $
else NbrFiles = 1

mrs_ridtrans, Imag, RidTrans, NbrScale=NbrScale,  blocksize=blocksize
NbrBand=RidTrans.NBRSCALE

n = randomn(seed, 32,32)
Tn = get_stat(n, TabStatName=TabStatName, verb=verb)
NStat= (size(tn))[1]

TabStat = fltarr(NStat, NbrBand-1)
if keyword_set(survival) then begin
     NSurvStat=1001
     TabSurvStat = fltarr(NSurvStat, NbrBand-1)
     TabSurvNu = fltarr(NSurvStat, NbrBand-1)
end


for f=0,NbrFiles-1 do begin
  if f NE 0 then Imag = mrs_read(TabFile[f])
  if f NE 0 then mrs_ridtrans, Imag, RidTrans, NbrScale=NbrScale,  blocksize=blocksize 
  for j =0, RidTrans.NBRSCALE-2 do begin
    Scale =  mrs_ridget(RidTrans, j, NormMad=NormMad, ImaMean=ImaMean, ImaMad=ImaMad)
    if j EQ 0 then TabStat = fltarr(NStat, RidTrans.NBRSCALE-1)
    if keyword_set(verb) then print, '   SCALE ', j+1
    TabStat[*,j] = get_stat(Scale, verb=verb)
    ; help, tabstat, RidTrans.NBRSCALE
     if keyword_set(survival) then begin 
	      a = survival(Scale, /norm, np=Np)
          TabSurvStat[*,j] = a[*,1]
	      TabSurvNu[*,j] = a[*,0]
      end
    end
     if keyword_set(TabFile) then BEGIN
      if f EQ 0 then BEGIN 
         vs  = size(TabStat)
         Nx = vs[1]
         Ny = vs[2]
	 TabAllStat = dblarr(nx, ny, NbrFiles)
	 if  keyword_set(survival) then   begin
	    vs  = size(TabSurvStat)
            Nx = vs[1]
            Ny = vs[2]
	    TabAllSurvStat = dblarr(nx, ny, NbrFiles)
        end
        END
       TabAllStat[*,*,f] = TabStat
       if  keyword_set(survival) then  TabAllSurvStat[*,*,f] = TabSurvStat
    END   

end

DONE:

 if keyword_set(TabFile) then return, TabAllStat $
  else return, TabStat

END    
 


