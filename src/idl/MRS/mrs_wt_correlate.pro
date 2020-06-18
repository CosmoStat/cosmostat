;+
; NAME:
;        mrs_wt_correlate
;
; PURPOSE:
;	 Return the correlation per wavelet scale of of two images.  If a mask is given, it calculates also the correlation in the mask only.
;
; CALLING:
;        WT_Correl= mrs_wt_correlate(Imag_1, Imag_2, WT1=WT1, WT2=WT2,  NbrScale=NbrScale, Mask=Mask,  NbrLat= NbrLat) 
;       
; OPTIONAL INPUTS:
;     Imag_1 -- IDL array of healpix map. If WT1 is set, Imag_1 is not used.
;     Imag_2 -- IDL array of healpix map.  If WT2 is set, Imag_2 is not used.
;
; INPUT KEYWORDS:
;		NbrScale -- int: Number of scales. Default is   log(nside)) / log(2) 
;
; INPUT/OUTPUT  KEYWORDS:
;      WT1 -- IDL structure: Wavelet transform of Imag_1
;      WT2 -- IDL structure: Wavelet transform of Imag_2
;	   NbrLat -- int: Number  of latitude bands.   
;
; RETURN: 
;   CorrelWT: IDL Stucture = {TabCorrel -- IDL array [0: NbrScale-1] :   TabCorrel[j] is the correlation between the two wavelets scale WT1.coef[*,j] and  WT2.coef[*,j] 
;                                            NbrScale -- int: Number of scales. 
;                                            TabMaskCorrel -- IDL array [0: NbrScale-1] :   TabMaskCorrel[j] is the correlation between the two wavelets scale WT1.coef[IndMask,j] and  WT2.coef[IndMask,j] 
;                                                                        where IndMask = where(Mask EQ 1)
;                                            TabLatCorrel -- IDL array [0: NbrScale-1, 0:NbrLat-1].   TabCorrel[j, l] is the correlation between the two wavelets scale WT1.coef[*,j] and  WT2.coef[*,j]  at a given 
; 								                                        latitute band
;											TabMaskLatCorrel -- IDL array [0: NbrScale-1, 0:NbrLat-1].    TabCorrel[j, l] is the correlation between the two wavelets scale WT1.coef[IndMask,j] and 
;																	    WT2.coef[IndMask,j]  at a given  latitute band
;																		TabBLat-- 1D float array [0: NbrLat-1]: B latitude for each latitude band }
; EXAMPLE:
;     Compute the multiscale correlation function with 5 scales
;       WT_Correl= mrs_wt_correlate(Imag_1, Imag_2, NbrScale=5, NbrLat= NbrLat) 
;
; EXTERNAL CALLS:
;            mrs_wttrans.pro,  mrs_lat_index.pro
;
; HISTORY:
;	Written: Jean-Luc Starck, April 2011
;-
;==================================================================================

function MCORRELATE, Ima1, Ima2
Ret=0.
if sigma(Ima1)  GT 0 and  sigma(Ima2) GT 0  then Ret = CORRELATE(Ima1, Ima2)
return, Ret
end

;==================================================================================
  
function mrs_wt_correlate, Imag_1, Imag_2, WT1=WT1, WT2=WT2,  NbrScale=NbrScale, Mask=Mask,  NbrLat= NbrLat

ind=-1
if N_PARAMS() LT 1   then begin 
       if not keyword_set(WT1) or not keyword_set(WT2) then begin
           print, 'CALLING SEQUENCE:  WT_Correl = mrs_wt_correlate(Imag_1, Imag_2, WT1=WT1, WT2=WT2,  NbrScale=NbrScale, Mask=Mask,  NbrLat= NbrLat) '
           goto, DONE
        end
end else  if N_PARAMS() LT 2   then begin 
        if not keyword_set(WT1) and not keyword_set(WT2) then begin
           print, 'CALLING SEQUENCE:  WT_Correl = mrs_wt_correlate(Imag_1, Imag_2, WT1=WT1, WT2=WT2,  NbrScale=NbrScale, Mask=Mask,  NbrLat= NbrLat) '
           goto, DONE
        end
end    

if keyword_set(WT1) then NbrScale =   WT1.NBRSCALE
if keyword_set(WT2) then NbrScale =   WT2.NBRSCALE

if not keyword_set(NbrScale) then begin
	  NbrScale=  alog(float(gnside(Imag_1))) / alog(2.)   ; nside=256 ==> 8, nside=2048 ==> 11
end

NbrScale = fix(NbrScale)
NbrBand=NbrScale
TabCorrel = dblarr(NbrScale)
TabMaskCorrel = dblarr(NbrScale)
 if keyword_set(Mask) then   indMask = where( Mask EQ 1, c)
if not keyword_set(NbrLat) then NbrLat = 25
TabLatCorrel = fltarr(NbrScale, NbrLat)
TabMaskLatCorrel = fltarr(NbrScale, NbrLat)

if not keyword_set(WT1) then mrs_wttrans, Imag_1, WT1, NbrScale=NbrScale 
if not keyword_set(WT2) then mrs_wttrans, Imag_2, WT2, NbrScale=NbrScale 

for j =0,NbrScale-1 do begin
          Scale1 = mrs_wtget(WT1,j,NormVal=NormVal)
          Nside = gnside(Scale1)
          Scale2 = mrs_wtget(WT2,j,NormVal=NormVal)
          TabCorrel[j] = MCORRELATE(Scale1, Scale2)
	  print, '    Scale ', j+1, ', Correl = ', TabCorrel[j]
          if keyword_set(Mask) then begin
              ScaleMask1 = Scale1[indMask]
              ScaleMask2 =  Scale2[indMask]
              TabMaskCorrel[j] = MCORRELATE(ScaleMask1, ScaleMask2)
          end
          for l=0, NbrLat-1 do begin
             IndLat =  mrs_lat_index(Nside, l, NbrLat=NbrLat, Mask=LatMask, MinLat=MinLat, MaxLat=MaxLat, MeanLat=MeanLat, TabBLat=TabBLat)
             ScaleLat1 =  Scale1[IndLat]
             ScaleLat2 =  Scale2[IndLat]
             TabLatCorrel[j,l] = MCORRELATE(ScaleLat1, ScaleLat2)
             if keyword_set(mask) then begin
                 ind = where(   Mask EQ 1 and LatMask EQ 1, c)
                 if c GT 0 then begin
                    ScaleLat1 =  Scale1[ind]
                    ScaleLat2 =  Scale2[ind]
                    TabMaskLatCorrel[j,l] = MCORRELATE(ScaleLat1, ScaleLat2)
                 end else TabMaskLatCorrel[j,l]=0
                 end 
             print, "         Scale ", j+1, ", Lat ", l+1, ", Correl = ", TabLatCorrel[j,l] , ',  MaskCorrel =  ', TabMaskLatCorrel[j,l]
         endfor
end

CorrelWT = {TabCorrel: TabCorrel,  NbrScale : NbrScale,  TabMaskCorrel: TabMaskCorrel,  TabLatCorrel: TabLatCorrel, TabMaskLatCorrel: TabMaskLatCorrel, TabBLat: TabBLat} 

DONE:

return, CorrelWT

end

;==================================================================================
 

