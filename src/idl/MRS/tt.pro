;+
; NAME:
;        mrs_wtstat
;
; PURPOSE:
;	 Return statistical information relative to the isotropic  wavelet transform of a given data set.
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
;		If TabFile is set, then the statistic is computed on a set of images. Tab[*,*,f] will be the statistic related to the file TabFile[f]
;       If a Mask is given, statisics from wavelet coefficients inside the mask is also calculated and store in the IDL structure StatWT.
;       If the keyword survival is set, then the survival function is also calculated and store in the structure StatWT.
;
; CALLING:
;
;      TabStat = mrs_wtstat( Imag, TabStatName=TabStatName, NbrScale=NbrScale,  verb=verb, Mask=Mask
;							wt=wt, survival=survival, TabFile=TabFile,  StatWT=StatWT) 
;       
; INPUTS:
;     Imag -- IDL array of healpix map: Input data to analyze
;
; INPUT KEYWORDS:
;		NbrScale -- int: Number of scales. Default is   log(nside)) / log(2) 
 ;		verb -- scalar, if set, the calculated statistics are printed on the screen
;		TabFile -- IDL table of string, list of file where the function read the maps to be analized, in that case, on output, Imag is the last map that had been proceed and the return value is a 3D array:
;	    Tab[i,j,f]	statistic i for scale j and map TabFile[f]
;		survival: if set, use the survival function and activate the keywords parameters TabSurvStat, TabAllSurvStat and TabSurvNu for the results
;
; OUTPUT KEYWORDS: 
;     TabStatName -- IDL table of string: TabStatName = ["Sigma", "Skewness", "Kurtosis", "Min", "Max", "HC1", "HC2", "CUMULANT ORDER 5", "CUMULANT ORDER 6"]
;      wt -- IDL structure, wavelet transform of the data (see mrs_wttrans), if TabFile keyword is used, the last map proceed
;      StatWT -- IDL structure, = { TabStat -- IDL array [0:8, 0:NbrScale-2, 0:NbrFiles-1]  :  computed statistics.
;                                                 TabSurvStat  -- float array [*,j, 0:NbrFiles-1] survival value at scale j
;	                                              TabSurvNu --   float array [*,j, 0:NbrFiles-1] nu survival value at scale j
;                                                 TabStatName -- IDL table of string: TabStatName = ["Sigma", "Skewness", "Kurtosis", "Min", "Max", "HC1", "HC2", "CUMULANT ORDER 5", "CUMULANT ORDER 6"]
;									  		   	TabFile -- IDL table of string, list of file where the function read the maps to be analized, in that case, on output, Imag is the last map that had been proceed and the return value is a 3D array:
;												NbrScale -- int: Number of scales. 
;												NbrFiles : --int : number of files.
;                                               TabMaskSurvStat -- same as TabSurvStat, but statistics computed only in the mask.
;                                               TabMaskStat  -- same TabStat, but statistics computed only in the mask.
;                                               TabLatStat -- IDL array [0:8, 0:NbrScale-2, 0:NbrLat-1, 0:NbrFiles-1],  statistics per scale  and per latitude band
;												TabMaskLatStat:  IDL array [0:8, 0:NbrScale-2, 0:NbrLat-1, 0:NbrFiles-1],  statistics per scale  and per latitude band, but only in the mask.
;												TabBLat: -- 1D float array [0: NbrLat-1]: B latitude for each latitude band} 
;    
; EXAMPLE:
;     Compute the statistics of  wavelet transform with 5 scales
;       TabStat = mrs_wtstat(Data, NbrScale=5, /verb,  StatWT=StatWT)
;
; EXTERNAL CALLS:
;
; HISTORY:
;	Written: Jean-Luc Starck, 2005
;	September, 2005 File creation
;-

;==================================================================================
 
function mrs_wtstat, Imag, NbrScale= NbrScale,   TabStatName=TabStatName,verb=verb, wt=wt, $
      survival=survival, TabFile=TabFile, TabSurvStat=TabSurvStat, TabAllSurvStat=TabAllSurvStat, TabSurvNu=TabSurvNu, Mask=Mask, TabMaskStat=TabMaskStat,    TabMaskSurvStat=TabMaskSurvStat, $
      TabMaskAllSurvStat=TabMaskAllSurvStat, TabMaskAllStat=TabMaskAllStat, ApplyMask=ApplyMask, StatWT=StatWT, Pyr=Pyr


if N_PARAMS() LT 1  and not keyword_set(TabFile) then begin 
        print, 'CALLING SEQUENCE:  TabStat = mrs_wtstat(Data, NbrScale= NbrScale,  survival=survival, Pyr=Pyr, TabStatName=TabStatName, Mask=Mask, StatWT=StatWT, verb=verb)'
        trans=0
        TabStat=0
        goto, DONE
        end

 Undec = 1 
if keyword_set(Pyr) then Undec =0 else Pyr = 0
if keyword_set(Mask) then begin
  		Undec = 1 
  		Pyr = 0 
  end
 
 TabSurvStat=0
 TabSurvNu=0
 TabMaskSurvStat=0
 TabMaskStat=0
 TabBLat=0
 TabLatStat=0
 TabMaskLatStat=0
 TabMaskLatAllStat=0
 TabLatAllStat =0
if keyword_set(WT) then Trans = WT

if keyword_set(TabFile) then NbrFiles = (size(TabFile))[1] else TabFile = 0
if keyword_set(TabFile) then Imag = mrs_read(TabFile[0]) $
else NbrFiles = 1
Nside = gnside(Imag)

if not keyword_set(NbrScale) then NbrScale=  alog(float(gnside(imag))) / alog(2.)  ; nside=256 ==> 8, nside=2048 ==> 11
NbrScale = fix(NbrScale)
NbrBand=NbrScale

; We call get_set to known how many statistics are computed, in order to fix array sizes.
n = randomn(seed, 32,32)
Tn = get_stat(n, TabStatName=TabStatName, verb=verb)
NStat= (size(tn))[1]
TabStat = fltarr(NStat, NbrScale-1)
NbrLat = 25
TabLatStat = fltarr(NStat, NbrScale-1, NbrLat)

if keyword_set(survival) then begin
     NSurvStat=1001
     TabSurvStat = fltarr(NSurvStat, NbrBand-1)
     TabSurvNu = fltarr(NSurvStat, NbrBand-1)
end

if keyword_set(mask) then begin
      indMask = where( Mask EQ 1, c)
      TabMaskStat = TabStat
      TabMaskLatStat = TabLatStat
      if keyword_set(survival) then TabMaskSurvStat = TabSurvStat
end

for f=0,NbrFiles-1 do begin
       if keyword_set(verb) then print, ' WT Imag  ', f+1
       if f NE 0 then Imag = mrs_read(TabFile[f])
      if keyword_set(Mask) and keyword_set(ApplyMask) then Imag = Imag * Mask
      if f GT 0 or not keyword_set(wt) then begin
        
         if Pyr EQ 0 then mrs_wttrans, Imag, Trans, NbrScale=NbrScale $
         else begin print, "Pyramidal WT" & mrs_pwttrans, Imag, Trans, NbrScale=NbrScale & end
     end
      for j =0,NbrScale-2 do begin
         Scale = mrs_wtget(Trans,j,NormVal=NormVal)
         if keyword_set(Mask) then ScaleMask = Scale[indMask]
       
        if keyword_set(verb) then print, '      SCALE ', j+1
        TabStat[*,j] = get_stat(Scale, verb=verb)
        if keyword_set(Mask) then begin  print, '      MASK SCALE '  & TabMaskStat[*,j] = get_stat(ScaleMask, verb=verb) & end
        if keyword_set(survival) then begin 
	            a = survival(Scale, /norm, np=Np)
                TabSurvStat[*,j] = a[*,1]
	           TabSurvNu[*,j] = a[*,0]
	        	if keyword_set(Mask) then  begin
		              Mask_a = survival(ScaleMask, /norm, np=Np)
		              TabMaskSurvStat[*,j] = Mask_a[*,1]
	           endif
          endif ; survival
          for l=0, NbrLat-1 do begin
             IndLat =  mrs_lat_index(Nside, l, NbrLat=NbrLat, Mask=LatMask, MinLat=MinLat, MaxLat=MaxLat, MeanLat=MeanLat, TabBLat=TabBLat)
             ScaleLat =  Scale[IndLat]
            ; help, ScaleLat
            ;info,   ScaleLat
            ;zz  = get_stat(ScaleLat, verb=verb)
            ;help, zz, j, l
            ;print, zz[2]
	        ;  if keyword_set(verb) then print, "     LAT ", l,  TabBLat[l]
             TabLatStat[*,j,l] = get_stat(ScaleLat, verb=0)
             if keyword_set(mask) then begin
                   ind = where(   Mask EQ 1 and LatMask EQ 1, c)
	            ;  if keyword_set(verb) then print, "     MASKLAT ", l, TabBLat[l]
	              TabMaskLatStat[*,j,l] = 0
                 if c GT 20 then TabMaskLatStat[*,j,l]  = get_stat(Scale[ind], verb=0)  
             endif 
         endfor ; for l
     endfor ;  j =0,NbrScale-2
     if keyword_set(TabFile) then BEGIN
         if f EQ 0 then BEGIN 
       		  vs  = size(TabStat)
      		   Nx = vs[1]
      		   Ny = vs[2]
      		   TabAllStat = dblarr(nx, ny, NbrFiles)
     		   TabLatAllStat = dblarr(nx, ny, NbrLat, NbrFiles)
			 if keyword_set(mask) then begin  TabMaskAllStat = TabAllStat & TabMaskLatAllStat = TabLatAllStat & end 
			 if  keyword_set(survival) then   begin
			        vs  = size(TabSurvStat)
  		              Nx = vs[1]
   		             Ny = vs[2]
			        TabAllSurvStat = dblarr(nx, ny, NbrFiles)
			        if keyword_set(mask) then   TabMaskAllSurvStat = TabAllSurvStat
            end
      END
      TabAllStat[*,*,f] = TabStat
      TabLatAllStat [*,*,*,f] = TabLatStat
     if keyword_set(mask) then   TabMaskLatAllStat [*,*,*,f] = TabMaskLatStat
      if keyword_set(mask) then  TabMaskAllStat[*,*,f] = TabMaskStat
      if keyword_set(survival) then begin
         TabAllSurvStat[*,*,f] = TabSurvStat
	 if keyword_set(mask) then TabMaskAllSurvStat[*,*,f] = TabMaskSurvStat
      end
    END   
endfor ; for f

if keyword_set(TabFile) then begin
TabStat= TabAllStat
TabSurvStat  =TabAllSurvStat
TabMaskSurvStat = TabMaskAllSurvStat
TabMaskStat = TabMaskAllStat
TabLatStat = TabLatAllStat
TabMaskLatStat =TabMaskLatAllStat
end

StatWT = {TabStat: TabStat, TabSurvStat: TabSurvStat , TabSurvNu: TabSurvNu,   TabStatName:TabStatName,  TabFile: TabFile, NbrScale : NbrScale,  NbrFiles : NbrFiles, $
TabMaskSurvStat: TabMaskSurvStat, TabMaskStat: TabMaskStat, $
TabLatStat: TabLatStat, TabMaskLatStat: TabMaskLatStat, TabBLat: TabBLat} 

DONE:
    wt=trans
    return, TabStat
END    
 


