;+
;  NAME: 
;        MRS_PLOT_LATPOWSPEC
;
;  PURPOSE: 
;        Plot the power spectrum per latitude        
;
; CALLING: 
;        mrs_plot_latpowspec, cl, title=title, xr=xr, yr=yr, south=south, left_legend=left_legend, 
; right_legend=right_legend, top_legend=top_legend, bottom_legend=bottom_legend, xlog=xlog
;
; INPUTS: 
;        IDL structure cl_g: 
;             l --- modes
;             mcl --- mean power spectrum
;             vcl --- variance on the mean power spectrum
;             tcl --- power spectrum per patch
;             mlcl --- power spectrum per latitude 
;             vlcl --- variance on the power spectrum per latitude
;             nlat --- number of latitudes
;             lat --- array of latitudes
;             nlati --- index of patches per latitude
;	     
;	
;  OUTPUTS: 
;             plot of power spectrum per latitude
;
;  KEYWORDS: 
;       plot keywords
;
; HISTORY:
;	Written: Sandrine Pires 2010.
;-
;-------------------------------------------------------------------------------

pro mrs_plot_latpowspec, cl, title=title, xr=xr, yr=yr, south=south, left_legend=left_legend, right_legend=right_legend, top_legend=top_legend, bottom_legend=bottom_legend, xlog=xlog
p = cl.mlcl
l = cl.l
ll1 = cl.l*(cl.l+1)
lat =  cl.LATITUDE
if not keyword_set(yr) then begin
   y =  ll1 *  p(*,0)  /2. /!pi
   yr=[0, max(y)]
end 
if not keyword_set(xr) then xr=[500, 3000]

if not keyword_set(title) then tit='Power Spectrum' $
else tit = title + ' ' + ' Power Spectrum' 
tek_color

if not keyword_set(south) then begin
  col= [1,2,3,4,5,6,7,8,9]
  Tabi=[3,5,7,9,11,13,15,17,19]
  i=0
  items = ['Latitude ' + string( lat[ Tabi[i]], format='(f5.1)')]
  for i=1,8 do items=[items, ['Latitude ' + string( lat[Tabi[i]], format='(f5.1)')]]
   i=0
   plot, l, ll1*  p(*,Tabi[i])  /2. /!pi, color=col[i],  xtitle='l',  charsize=1.6,  xrange=xr, yrange=yr, title=tit, xlog=xlog
  for i=1,8 do oplot, l, ll1*  p(*,Tabi[i])  /2. /!pi, color=col[i]
  legend, textcol=col, items,   colors=col, charsize=2, spacing=3, left_legend=left_legend, right_legend=right_legend, top_legend=top_legend, bottom_legend=bottom_legend
endif else   begin
  col= [1,2,3,4,5,6,7,8,9]
  Tabi=[2,4,6,8,10,12,14,18,20]
  i=0
  items = ['Latitude ' + string( lat[ Tabi[i]], format='(f5.1)')]
  for i=1,8 do items=[items, ['Latitude ' + string( lat[Tabi[i]], format='(f5.1)')]]

  i=0
  plot, l, ll1*  p(*,Tabi[i])  /2. /!pi, color=col[i],  xtitle='l',  charsize=1.6,  xrange=xr, yrange=yr, title=tit, xlog=xlog
  for i=1,8 do oplot, l, ll1*  p(*,Tabi[i])  /2. /!pi, color=col[i]
  legend, textcol=col, items,   colors=col, charsize=2, spacing=3,left_legend=left_legend, right_legend=right_legend, top_legend=top_legend, bottom_legend=bottom_legend
endelse
end
