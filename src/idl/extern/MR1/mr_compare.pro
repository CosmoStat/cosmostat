PRO MR_COMPARE, ImaRef, Ima_Or_Cube, TabCorrel, TabRMS, TabSNR, Nscale=Nscale, plot=plot, title=title
;+ 
; NAME: 
;     MR_COMPARE
;
; PURPOSE: 
;     Comparison between a reference image and an image or a sequence of image.
;     The comparison is done in the multiresolution scales.
;     if Ima_Or_Cube is a cube,  the treament is repeated on each image of the
;     cube (*,*,i).
;
;     For each image ima_i of the cube Ima_Or_Cube
;  -     We compute the wavelet transform WaveRef of ImaRef
;  -     We compute the wavelet transform WaveIma of ima_i
;
;  -     we calculate the  correlation at each scale s:
;        if TabCorrel(s,i) = 1 then WaveRef(*,*,s) = WaveIma(*,*,s)
;                                    are identical
;        else TabCorrel(s,i) is less than 1, and some differences exit.
;
;  -     we calculate the  RMS at each scale:    
;            TabRMS(s,i) = sigma( WaveRef(*,*,s) - WaveIma(*,*,s))
; 
;  -     we calculate the normalized SNR at each scale (in dB): 
;            TabSNR(s,i) = 10 alog10 ( mean(WaveRef(*,*,s)^2) / TabRMS(s,i)^2)
;
; CALLING SEQUENCE: 
;   MR_COMPARE, ImaRef, Ima_Or_Cube, TabCorrel, TabRMS, TabSNR, Nscale=Nscale
;
; INPUTS: 
;   ImaRef -- IDl 2D array : reference image
;   Ima_Or_Cube -- IDl 2D or 3D  array: image or sequence of image 
;                                       to be compared
;
; KEYED INPUTS: 
;   Nscale -- scalar: number of scales for the comparison
;   PLOT -- scalar: if set, then plot the results
;   title -- string: is set, add the title to the plots
;
; OUTPUTS: 
;    TabCorrel -- 1D or 2D IDL array: correlation table
;    TabRMS  -- 1D or 2D IDL array: RMS table
;    TabSNR -- 1D or 2D IDL array: SNR table (in dB)
;
; ALGORITHM: 
;    see J.L. Starck, and A. Bijaoui, "Filtering and Deconvolution by
;          the Wavelet Transform", Signal Processing, 35, 195--211, 1994.
;
; EXAMPLE:
;    ; we want compare two images (Ima1 to ImaRef).
;    mr_compare, imaref, ima1, tc, tr, ts, /plot, title='Comparison ImaRef-Ima1'
;
; MODIFICATION HISTORY: 
;    22-Jan-1996 JL Starck 
;-
 
;------------------------------------------------------------
; parameters check
;------------------------------------------------------------
 
 IF N_PARAMS() LT 4 THEN BEGIN
   PRINT, 'CALLING SEQUENCE: ', $ 
    'MR_COMPARE, ImaRef, Ima_Or_Cube, TabCorrel, TabRMS, Nscale=Nscale, plot=plot, title=title'
   GOTO, CLOSING
 ENDIF
 
if not keyword_set(Nscale) then Nscale = 4

vsize =size(ImaRef)
if vsize[0] NE 2 then BEGIN
   PRINT, 'Error: first parameter must be an image' 
   GOTO, CLOSING
ENDIF
vsize =size(Ima_Or_Cube)
if vsize[0] NE 2 and vsize[0] NE 3 then BEGIN
   PRINT, 'Error: second parameter must be an image or a cube' 
   GOTO, CLOSING
ENDIF
iscube=0
if vsize[0] EQ 3 then iscube=1

;------------------------------------------------------------
; function body
;------------------------------------------------------------
 
if iscube EQ 0 then BEGIN
   dim = 1
   TabCorrel = fltarr(Nscale) 
   TabRMS = fltarr(Nscale)
   TabSNR = fltarr(Nscale)
END ELSE BEGIN
   dim = vsize[3]
   TabCorrel = fltarr(Nscale, dim) 
   TabRMS = fltarr(Nscale, dim)
   TabSNR = fltarr(Nscale, dim)

END

OPT = '-n '+strcompress(string(Nscale),/REMOVE_ALL)
mr_transform, ImaRef, WaveRef, OPT=Opt

for n=0, dim-1 do $
BEGIN
  if iscube EQ 0 then mr_transform, Ima_Or_Cube, WaveIma, OPT=Opt $
  else BEGIN
     Ima = Ima_Or_Cube[*,*,n]
     mr_transform, Ima, WaveIma, OPT=Opt
  END
 
  for s=0, Nscale-1 do BEGIN
     ScaleRef = WaveRef[*,*,s]
     ScaleIma = WaveIma[*,*,s]

     RMS = sigma(ScaleRef - ScaleIma)
     Correl = total(ScaleRef*ScaleIma) / $
                              sqrt( total(ScaleRef^2)*total(ScaleIma^2))
     SNR = 10. * alog10 (  mean( ScaleRef^2) / RMS^2 )
     if iscube EQ 1 then BEGIN
               TabCorrel[s,n] = Correl
               TabRMS[s,n] = RMS
               TabSNR[s,n] = SNR
     END ELSE  BEGIN
               TabCorrel[s] = Correl
               TabRMS[s] = RMS
               TabSNR[s] = SNR
     END
  END
  if keyword_set(plot) then BEGIN
     p = !p.multi
     !p.multi = [0,1,3]
     if keyword_set(title) then TIT = title+":"+" Correlation versus scale" $
     else TIT = " Correlation versus scale"
     plot,TabCorrel(*,n),xtitle='Wavelet scale',ytitle='Correlation',title=TIT

     if keyword_set(title) then TIT = title+":"+" RMS versus scale" $
     else TIT = " RMS versus scale"
     plot, TabRMS(*,n), xtitle='Wavelet scale', title=TIT, ytitle='RMS'

     if keyword_set(title) then TIT = title+":"+" SNR versus scale" $
     else TIT = " SNR versus scale"
     plot, TabSNR(*,n), xtitle='Wavelet scale', title=TIT, ytitle='SNR (dB)'
     !p.multi = p
  END
END
 
;------------------------------------------------------------
; closing
;------------------------------------------------------------
 
 CLOSING:
 
  RETURN
 
 END
