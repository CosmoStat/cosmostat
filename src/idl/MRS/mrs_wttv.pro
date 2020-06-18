;+
; NAME:
;        mrs_wttv
;
; PURPOSE:
;		Visualization of the wavelet transform obtained by the command mrs_wttrans or
;       mrs_pwttrans. If the keyword WRITE is set to a string, then all the scales
;       are written on the disk as PNG files, and the string is used as
;       a prefix for the file name of the different scales.
;          FILENAME = WRITE + ' scale_ ' + ScaleNumber + '.png'
;
; CALLING:
;
;     mrs_wttv,  Trans, Tit=Tit, write=write, graticule=graticule, min=min, max=max, big=big 
;       
; INPUTS:
;     Trans -- Trans -- IDL structures containing the wavelet transform 
;
; KEYWORDS:
;      Tit  : string -- Title of the plot
;      Write : string -- Prefix filename. If set, write to disk each scale of the wavelet transform in PNG format
;      graticule: int -- Mollview Healpix command graticule keyword
;	   min max : float Mollview Healpix command min or max keyword, new min or max value to be used for the display
;	   big : bool if set, Mollview Healpix keyword pxsize is set to 1500
;      mad: bool. If set, visualize with a dynamic range between [Mean - Kmad*MAD, Mean+ Kmad*MAD],  where 
;                      mean is the mean of the data, MAD is the median absolute deviation and Kmad a parameter with default value=3
;      kmad: float -- default value is 3. If set, visualize with a dynamic range between [Mean - Kmad*MAD, Mean+ Kmad*MAD]
;
; EXTERNAL CALLS:
;       mrs_wtget
;       mollview (healpix software)
;
; HISTORY:
;	Written: Jean-Luc Starck, 2005
;	February, 2005 File creation
;-

pro mrs_wttv, Data, Tit=Tit, write=write, graticule=graticule, min=min, max=max, big=big, mad=mad, kmad=kmad

for j=0,Data.NbrScale-1 do begin
   if not keyword_set(tit) then Tit=' '
   if keyword_set(write) then png=write+'_scale_' + strcompress(string(j+1),/remove_all) + '.png'
   Title= Tit + ' Scale ' + strcompress(string(j+1),/remove_all)
   if keyword_set(write) then print, 'Scale ', j+1, ' : Filename = ', png
   Scale = mrs_wtget(Data,j)
   if keyword_set(Mad) or keyword_set(kmad) then begin
      MadScale = mad(Scale)
      if not keyword_set(kmad) then kmad=3.
      min = mean(scale) - kmad * MadScale
      max = mean(scale) + kmad * MadScale
   end
   
   mrs_tv, Scale,  TITLE=Title, png=png,  graticule=graticule, /healpix, min=min, max=max, big=big
   end
end
