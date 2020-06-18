;+
; NAME:
;        mrs_tv
;
; PURPOSE:
;	Visualization  of a HEALPix image (NESTED data representation) or a GLESP image.
;
; CALLING:
;
;     mrs_tv, Data, graticule=graticule, gif=gif, png=png, TITLE=TITLE, COLT=COLT, NOBAR=NOBAR, Healpix=Healpix, PS=PS, log=log, min=min, max=max, pxsize=pxsize, big=big, x=x, pol=pol, UNITS=UNITS
;       
; INPUTS:
;     Data -- IDL array of healpix map or GLESP structure: Input image to be visualized 

; KEYWORDS:
;		min: float -- Data image is visualized with the new min set.
;		max: float -- Data image is visualized with the new max set.
;		log: if set, plot the image in log scale, Data must be positive.
;		graticule: int -- Mollview Healpix command graticule keyword.
;		png: string -- If set, write to the disk a PNG file with the filename given by the png keyword. 
;		gif: string -- If set, write to the disk a GIF file with the filename given by the gif keyword.
;		PS: string -- If set, write to the disk a Postcript file with the filename given by the ps keyword.
;		TITLE: IDL string -- Title of the plot in the Healpix  representation.
;		UNIT: string -- name of the image's unit to be plotted with the LUT in the Healpix representation
;		COLT: int --  IDL Color table.
;		NOBAR: int -- Do not plot the LUT in the Healpix representation
;		Healpix: scalar -- if set, then convert the GLESP image into a Healpix one, and use
;                          the Healpix representation for visualization.
;                          This keyword is not active for maps already in Healpix representation
;		pxsize: int set the number of horizontal pixel on the plot (it will be the same in vertical), default value is 800
;		big: set pxsize to 1500
;		pol: int -- for plotting polarized map T,Q,U default is 0 (no polarized map), see mollview help for more details.
;
;		x: if set, start interactive plot with the mapview prog, nside max=1024, data will be resized if greater
;
; EXTERNAL CALLS:
;      mollview - Healpix command
;      f2fig  -- GLESP C binary
;      
; HISTORY:
;	Written: Jean-Luc Starck, 2006
;	January, 2006 File creation
;-

;======================================================

pro tv_glesp_rot, g_in, theta=theta , phi=phi, COLT=COLT, maxc=maxc,minc=minc,title=title

 if N_params() ne 1 then begin print,'Error: no GLESP map parameter in tv_glesp_rot'
 goto,fin
  endif
  
x_sky = g_in.x_sky
y_sky = g_in.y_sky
t_sky = g_in.t_sky

write_glesp, 'tmp_rot.fits', g_in
if keyword_set(theta) or keyword_set(phi) then begin
   command = "difmap tmp_rot.fits"
   if keyword_set(theta) then command = command + " -dt "+string(theta)
   if keyword_set(phi) then command = command + " -dp "+string(phi)
   command = command + ' -o tmp_rot.fits'
   spawn,command
endif 

command = 'f2fig tmp_rot.fits -o tmp_rot.gif'
if keyword_set(minc) or keyword_Set(maxc) then command = command + ' -Cs '+strcompress(string(minc)+','+string(maxc),/remove_all)
spawn,command

read_gif,'tmp_rot.gif',imag,r,g,b

delete, 'tmp_rot.gif'
delete, 'tmp_rot.fits'

if defined(COLT) then loadct, COLT
tvlct,r,g,b
window, /free,xsize=800,ysize=440
tv,imag
fin:
end

;======================================================


pro tv_glesp_rot3, g_in, theta=theta , phi=phi, COLT=COLT,xs=xs,ys=ys,earth=earth,minc=minc,maxc=maxc,title=title

 if N_params() ne 1 then begin print,'Error: no GLESP map parameter in tv_glesp_rot'
 goto,fin
  endif
  
x_sky = g_in.x_sky
y_sky = g_in.y_sky
t_sky = g_in.t_sky
for i=0,2 do begin
write_glesp, 'tmp_rot.fits', g_in(i)
if keyword_set(theta) or keyword_set(phi) then begin
   command = "difmap tmp_rot.fits"
   if keyword_set(theta) then command = command + " -dt "+string(theta)
   if keyword_set(phi) then command = command + " -dp "+string(phi)
   command = command + ' -o tmp_rot.fits'
   spawn,command
endif 
command = 'f2fig tmp_rot.fits -o tmp_rot.gif'
if keyword_set(xs) and keyword_Set(ys) then command = command + ' -xs '+string(xs)+' -ys '+string(ys)
if keyword_set(minc) or keyword_Set(maxc) then command = command + ' -Cs '+strcompress(string(minc)+','+string(maxc),/remove_all)
;print,minc,maxc
commmand = command + ' >/dev/null'
spawn,command
print,command
if i eq 0 then begin
read_gif,'tmp_rot.gif',imag2,r,g,b
imag = fltarr((size(imag2))(1),(size(imag2))(2),3)
imag(*,*,0) = imag2
endif else begin
read_gif,'tmp_rot.gif',imag2,r,g,b
imag(*,*,i) = imag2
endelse

delete, 'tmp_rot.gif'
delete, 'tmp_rot.fits'
endfor
help,imag
info,imag
imag = 255*(imag - min(imag))/(max(imag)-min(imag))
if defined(COLT) then loadct, COLT
if keyword_set(earth) then imag = reverse(imag,1)
;tvlct,r,g,b

;tv,imag,true=3
temp_var=  findgen(256)
tvlct,temp_var,temp_var,temp_var
if not(keyword_set(xs) and keyword_Set(ys)) then begin

window, /free,xsize=800,ysize=440 
;tv,imag,true=3 
endif else begin 
    window, /free,xsize=xs,ysize=ys*1.1,title=title
    ;tv,imag,true=3
    v_band = ceil(ys*0.005)
    v_deb = floor(ys*0.02)
    h_min = floor(xs*0.2)
    h_max = floor(xs*0.8)
    h_band = h_max-h_min
    band = (fltarr(v_band)+1)##(findgen(h_band)*255./h_band)
    imag[0:xs-1,0:ys*0.1-1,*]=0
    imag[h_min:h_max-1,v_deb:v_deb+v_band-1,0]  = band
    
    imag[h_min:h_max-1,v_deb+(v_band):v_deb+(v_band*2)-1,1]  = band
    imag[h_min:h_max-1,v_deb+(v_band):v_deb+(v_band*2)-1,0]  = band
    
    imag[h_min:h_max-1,v_deb+(v_band*2):v_deb+(v_band*3)-1,1]  = band
    
    imag[h_min:h_max-1,v_deb+(v_band*3):v_deb+(v_band*4)-1,2]  = band
     imag[h_min:h_max-1,v_deb+(v_band*3):v_deb+(v_band*4)-1,1]  = band
  
    imag[h_min:h_max-1,v_deb+(v_band*4):v_deb+(v_band*5)-1,2]  = band
  
     imag[h_min:h_max-1,v_deb+(v_band*5):v_deb+(v_band*6)-1,2]  = band
     imag[h_min:h_max-1,v_deb+(v_band*5):v_deb+(v_band*6)-1,0]  = band
  
   
    imag[h_min:h_max-1,v_deb+(v_band*6):v_deb+(v_band*7)-1,0]  = band
    imag[h_min:h_max-1,v_deb+(v_band*6):v_deb+(v_band*7)-1,1]  = band
    imag[h_min:h_max-1,v_deb+(v_band*6):v_deb+(v_band*7)-1,2]  = band
    
    imag[h_min-1:h_max,v_deb-1,*] = 255
    imag[h_min-1:h_max,v_deb+(v_band*7),*] = 255
    imag[h_min-1,v_deb-1:v_deb+(v_band*7),*] = 255
    imag[h_max,v_deb-1:v_deb+(v_band*7),*] = 255    
endelse

tv,imag,true=3
  xyouts,h_min-150,v_deb+(v_band*2),string(minc),/device,color=255, CHARSIZE = 2
    
    xyouts,h_max,v_deb+(v_band*2),string(maxc),/device,color=255, CHARSIZE = 2
  ;   xyouts,h_min-20,v_deb+(v_band*5),string(minc),color=255
  ;  
 ;   xyouts,h_max+20,v_deb+(v_band*5),string(maxc),color=255
  ;;  xyouts,500,500,"toto toto",color=128
  ;  XYOUTS, [0, 200, 250], [200, 50, 100], $  
  ; ['This', 'is', 'text'], CHARSIZE = 3, /DEVICE  
;tv,imag,true=3


fin:
end

;======================================================

pro winb, win=win
if not keyword_set(win) then window,  xsize=1024, ysize=1024 $
else window, win, xsize=1024, ysize=1024
end

pro tvbf, map, Face=Face, NumF=NumF, x=x, win=win
if keyword_set(win) then wset, win
if not keyword_set(Face) then Face=H2F(map)
if not keyword_set(NumF) then NumF=1
Ima = Face[*,*,NumF]
tvscl, congrid(Ima, 1024, 1024)
if keyword_set(x) then atv, Ima
end

pro winbs, win=win
if not keyword_set(win) then win=0
window, win, xsize=1600, ysize=500
end

;======================================================

pro mrs_xtv, Data
ns = gnside(Data)
if ns GT 1024 then mapview, array=mrs_resize(Data, nside=1024) $
else mapview, array=Data
end

;======================================================
pro mrs_tv, Data, graticule=graticule, gif=gif, png=png, TITLE=TITLE, COLT=COLT, NOBAR=NOBAR, Healpix=Healpix, PS=PS, log=log, min=min, max=max, pxsize=pxsize, big=big, x=x, pol=pol, window=window, UNITS=UNITS, charsize=charsize
  if keyword_set(big) then pxsize=1500
  if type_code(data) EQ 8 and not keyword_set(Healpix) then tv_glesp_rot, data, minc=min, maxc=max $
  else begin 
        if type_code(data) eq 8 then H = glesp2healpix(Data) else H = Data
 	if keyword_set(gif) then png='tt.png'
	mollview, H, /nested, NOBAR=NOBAR, graticule=graticule, png=png, TITLEPLOT=TITLE, COLT=COLT, log=log, min=min, max=max, pxsize=pxsize, pol=pol, window=window, /silent, charsize=charsize
	if keyword_set(gif) then begin
	   x = read_png(png)
	   write_gif, gif, x
	   delete, png
       end
  end
  if keyword_set(PS) then xdump, filename=PS
  if keyword_set(x) then begin
    ns = gnside(Data)
    if ns GT 1024 then mapview, array=mrs_resize(Data, nside=1024) $
    else mapview, array=Data
  end
end

;======================================================
 
pro tvs, Data, graticule=graticule, gif=gif, png=png, TITLE=TITLE, COLT=COLT, NOBAR=NOBAR, Healpix=Healpix, PS=PS, log=log, min=min, max=max, pxsize=pxsize, big=big, pol=pol, UNITS=UNITS, window=window, charsize=charsize
mrs_tv, Data, graticule=graticule, gif=gif, png=png, TITLE=TITLE, COLT=COLT, NOBAR=NOBAR, Healpix=Healpix, PS=PS, log=log, min=min, max=max, pxsize=pxsize, big=big,  pol=pol, UNITS=UNITS, window=window, charsize=charsize
end

;======================================================

pro tvsf, Data, f, win=win, hd=hd

if not keyword_set(hd) then h = h2f(data)
if keyword_set(win) then begin
   vs = size(h)
   nx = vs[1]
   ny = vs[2]
   window, win, xsize=nx, ysize=ny
   end
   
loadct, 15
tvi, h[*,*,f], size=1024
end

;======================================================
 
pro tvs3, Data, earth=earth,minc=minc,maxc=maxc,title=title
if not keyword_set(minc) then minc = min([min(data[0].t_sky),min(data[1].t_sky),min(data[2].t_sky)])
if not keyword_set(maxc) then maxc = max([max(data[0].t_sky),max(data[1].t_sky),max(data[2].t_sky)])

tv_glesp_rot3,data,xs=1400,ys=700,earth=earth,minc=minc,maxc=maxc,title=title
;tv_glesp_rot3,data,phi=-180,xs=1400,ys=700,earth=earth,minc=minc,maxc=maxc

;tv_glesp_rot3,data,xs=1200,ys=600

end

pro diff3D,recons,ini,res,sort_sky,res2,loc,integre
res = ini
sort_sky = ini
for i =0,2 do begin
help,res[i].t_sky
;res[i].t_sky = abs((recons[i].t_sky-ini[i].t_sky)/ini[i].t_sky )
res[i].t_sky = abs((recons[i].t_sky-ini[i].t_sky));/ini[i].t_sky )

ind = sort(res[i].t_sky)
sort_sky[i].t_sky  = (res[i].t_sky)(ind)
endfor
tvs3,res,/earth,minc=0,maxc=10
loadct,5
window,/free
plot,sort_sky[0].t_sky,color=255,yrange=[0,1]
oplot,sort_sky[1].t_sky,color=128
oplot,sort_sky[2].t_sky,color=100
taille = n_elements(sort_sky[0].t_sky)
tab = fltarr(taille,3)
tab[*,0] = sort_sky[0].t_sky
tab[*,1] = sort_sky[1].t_sky
tab[*,2] = sort_sky[2].t_sky

tab = reform(tab)
help,tab
res2 = histogram(tab,nbins=10000,locations=loc,max=50)
window,/free
info,loc
help,loc
plot,loc,res2,xrange=[0,10],color=255
print,'erreur @ ' ,loc[100],total(res2[100:9999])/taille*3
print,'erreur @ ' ,loc[1000],total(res2[1000:9999])/taille*3
print,'erreur @ ' ,loc[9999],total(res2[9999:9999])/taille*3
print,'erreur @ ' ,loc[00],total(res2[0:9999])/taille*3

integre = total(res2,/cumulative)
window,/free,title='integrated'
;plot,loc,integre
plot,loc,100*(integre(9999)-integre)/integre(9999),color=255,xtitle='value',ytitle='pourcent of pixel error'
end

pro make_film,nb_imag,nom_gif


for i = 0 ,nb_imag-1 do begin

mask_g = mrs_read("mask_g.fits")
simu = mrs_read('simu_0.fits')
name = strcompress('film_'+string(i)+'.fits',/remove_all)
map = mrs_read(name)
map.t_sky = map.t_sky*(1-mask_g.t_sky) + simu.t_sky*mask_g.t_sky
mrs_write,'tmp_g.fits',map
name = 'tmp_g.fits'
command = 'f2fig '+name+' -o tmp_rot.gif'
spawn,command
read_gif,'tmp_rot.gif',imag,r,g,b
write_gif,nom_gif,imag,r,g,b,/multiple
write_gif,nom_gif,imag,r,g,b,/multiple


endfor
write_gif,nom_gif,imag,r,g,b,/close

end

