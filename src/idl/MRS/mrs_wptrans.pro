;+
; NAME:
;        mrs_wptrans
;
; PURPOSE:
;   Computes the undecimated isotropic wavelet PACKET transform on the sphere with meyer wavelets, 
;   using the HEALPix representation (NESTED data representation) or the GLESP representation
;
; CALLING:
;     mrs_wptrans, Imag, Trans, NbrScale=NbrScale
;
; INPUTS:
;     Imag -- IDL array of healpix map or GLESP IDL structure: Input image to be transformed 
; 
; OUTPUTS:
;     Trans -- IDL structures with the following fields: 
;                     UseGLESP : int = 1 if we use the GLESP pixelisation 
;                  NbrScale : int = number of wavelet packet scales 
;                     nside : int = Healpix nside parameter
;                        nx : int: number of rings in GLESP representation
;                        np : int: number of pixels on the central ring
;                       npix : int: total number of pixels in the original image
;                       Coef : fltarr[npix,NbrScale] = wavelet transform of the data 
;                      lmax: int= lmax parameter at the finest scale
;                       x_sky: float array: ring position in cos(theta) ( in [-1,1]
;                       y_sky: number of pixels per ring
;
;
; OUTPUT KEYWORDS:
;		NbrScale : int = number of scales used in the decomposition
;
; EXTERNAL CALLS:
;       mrs_almtrans  
;   	mrs_almrec
;
; EXAMPLE:
;
;       Compute the wavelet packet of an image I with default options
;        The result is stored in Trans
;              
;       mrs_wptrans, Imag, Trans
;
;
;	IMAGE RECONSTRUCTION:
;
;		The image reconstruction is done by using the procedure mrs_wtrec, see mrs_wtrec.pro
;         
; HISTORY:
;	Written: Pierrick Abrial & Jean-Luc Starck, 2007
;	December, 2006 File creation
;--------------------------------------------------------------------------------------------------------



pro get_mirror_basis, lmax, nbr_ech, res, final, lmax_ech
; Create the mirror basis filters
; lmax: input = MAX A_lm
; nbr_ech: number of wavelet scales: the number of WP will be Np = 2*(nbr_ech+1)+1
; Res = fltarray [0:lmax+1,  Np]     H filters related to each scale
; final = fltarray [0:lmax+1,  Np]   G filters related to each scale
;                                    Final can be used as the Wave input in mrs_wptrans
; lmax_ech = fltarray [Np]          lmax_ech[j] = lmax of scale j

res = fltarr(lmax+1,2*(nbr_ech+1)+1)
lmax_ech = fltarr(2*(nbr_ech+1)+1)
mil = (nbr_ech)
ech = 0.5
nlmax = lmax
;lmax =lmax/2

for i= 1, nbr_ech  do begin
lmax =lmax/2

compute_h, lmax/2, ech, h
res(lmax-lmax/2:lmax,mil-i+1) = h
res(0:lmax-lmax/2,mil-i+1)=1
lmax_ech(mil-i+1) = lmax

inv = lmax/2+1-findgen(lmax/2+1)
;res(0:lmax/2,mil+i) = 1
print,nlmax-(lmax-lmax/2),nlmax-(lmax),mil+i
res(0:nlmax-lmax,mil+i) = 1
res(nlmax-(lmax):nlmax-(lmax-lmax/2),mil+i) =1-h(inv)
lmax_ech(mil+i) = nlmax-(lmax-lmax/2)
;scale = scale +1
;ech =ech*2
;lmax =lmax/2
endfor
res(0:lmax/2,0) = h
res(0: nlmax-lmax/2  ,mil+i) = 1
res(nlmax-lmax/2:nlmax,mil+i) =1-h(inv)
res(*,2*(nbr_ech+1)) =1
final =res


lmax_ech(mil+i) = nlmax 
lmax(0) = lmax/2 

for i= 0, 2*nbr_ech+1  do begin

final(*,i+1) = res(*,i+1)-res(*,i)

endfor
end


;==============================================

pro plot_mirror, final
; PLot all filters
n_element= (size(final))(2)
window,/free
plot,final(*,0),color=10,yrange=[-1,1]
for i = 0 , n_element-1 do begin
oplot,final(*,i),color=10 + i *250 /n_element
endfor
end

;==============================================

; ondelette spline
pro get_wp_spline_filter, interval=interval, final, win, lmax=lmax, nside=nside, hfilter=hfilter
; interval = input keyword intarray [NScale]. interval[j] = lmax of scale j
;            LMAX = MAX( interval )
;            NbrScale = (size(interval)[1])
;            By default, NbrScale=22 
;                        LMAX = 3000
; final = output array [0:LMAX, NbrScale]
; win = intarray [NbrScale,2] :  The WP scale j is defined between  win[j,0] and  win[j,1]

cl_fine =[   64,128,192,256,320,384,512,1024,1280,1536,1664,1792,1920,2048,2176,2304,2432,2560, 2688, 2816,2944,3000]
other_cl =[   64,128,256,320,384, 512,1024, 1280, 1536,1792,  2048, 2304, 2560,2688, 2816,3000]
cl_1024 = [32,64,80,96, 112,128,144,160, 176, 192,256, 272, 288,304, 320,336,352, 384,448,480,512, 576, 640, 704,736,768, 896,960,1024]
cl_512  =  [32,64,80,96, 112,128,144,160, 176, 192,256, 272, 288,304, 320,336,352, 384,448,480,509]
cl_128 = [32,64,80,96, 112,128]
cl = cl_fine;128

if not keyword_set(interval) then cl = cl_fine else cl= interval

n = n_elements(cl)
lmax= cl(n-1)+1

res = fltarr(lmax,n+1)
final = fltarr(lmax,n+1)
win = fltarr(n+1,2)

compute_h,cl(0),0.5,h
res[0:cl[0],0] = h[0:cl[0]]
win(0,0) = 0
win(0,1) = cl[0]
win(1,0) = cl[0]
for i = 1,n-1 do begin
compute_h,cl(i)-cl(i-1), 0.5, h
;help,h,cl[i-1],cl[i]
;help,res[cl[i-1]:cl[i],i] , h[0:cl[i]-cl[i-1]]
res[cl[i-1]:cl[i],i] = h[0:cl[i]-cl[i-1]]
res[0:cl[i-1],i] = 1

endfor
res[*,n]=1

final[*,0] = res[*,0]
for i= 0, n-1  do begin

final(*,i+1) = res(*,i+1)-res(*,i)

endfor
hfilter=res[*,0:n-1]

for i= 2, n-1  do begin

win(i,0) = cl[i-2]
win(i,1) = cl[i]

endfor
win(0,0) = 0
win(0,1) = cl[0]
win(1,0) = cl[0]
win(1,0) = 0
win(1,1) = cl[1]
win(n-1,0) = cl[n-3]
win(n-1,1) = cl[n-1]
win(n,0) = cl[n-2]
win(n,1) = cl[n-1]

end

;===================================================

function get_cl_interval, lmax=lmax, nside=nside
; return the lmax value for each band of the wavelet packet decomposition
COMMON C_PLANCK

if not keyword_set(nside)  then nside = 2048

cl_1024 = [32,64,80,96, 112,128,144,160, 176, 192,256, 272, 288,304, 320,336,352, 384,448,480,512, 576, 640, 704,736,768, 896,960,1024]
; JALAL decomposition for nside=1024, many small bands

if not keyword_set(lmax)  then begin
  if keyword_set(nside) then lmax = 3*nside else lmax = P_LMAX
  if lmax GT P_LMAX then lmax = P_LMAX
end
if nside EQ 2048 then cl_fine =[64,128,192,256,320,384,512,768,1024,1280,1536,1664,1792,1920,2048,2176,2304,2432,2560, 2688, 2816,3000] $
else if nside EQ 1024 then cl_fine = [64,128,192,256,320,384,512,768,1024,1280,1536,1664,1792,1920,2048,2176,2304,2432,2560, 2688, 2816,3000] $
else if nside EQ 512 then cl_fine = [64,128,192,256,320,384,512,768,1024,1280,1536] $
else if nside EQ 256 then cl_fine =[64,128,192,256,320,384,512,768] $
else if nside EQ 128 then cl_fine =[64,80,96, 112,128,144,160, 176, 192,256, 272, 288,304, 320,336,352, 384] $
else if nside EQ 64 then cl_fine =[16,32,64,80,96, 112,128,144,160, 176, 192] $
else if nside EQ 32 then cl_fine =[8,16,32,64,92] $
else if nside EQ 16 then cl_fine =[8,16,48] $
else cl_fine =[4,8,16,24] 

ind = where(cl_fine LE lmax)
cl = cl_fine[ind]
vs = size(cl)
Nb = vs[1]

if  cl[Nb-1] eq  lmax then NewCl = Cl $
else begin
   NewCl = [Cl, lmax]
end

;print, nside, lmax
;print, "NCL = ", newcl
return, NewCl
end

;===================================================

;ondelette de meyer
pro get_wp_meyer_filter, interval=interval, final, win, lmax=lmax, nside=nside, hfilter=hfilter
; interval = input keyword intarray [NScale]. interval[j] = lmax of scale j
;            LMAX = MAX( interval )
;            NbrScale = (size(interval)[1])
;            By default, NbrScale=22 
;                        LMAX = 3000
; lmax = input keyword, if set and interval not set, then the interval is calculated from get_cl_interval
; final = output array [0:LMAX, NbrScale]
; win = output  intarray [NbrScale,2] :  The WP scale j is defined between  win[j,0] and  win[j,1]

if not keyword_set(interval) then cl = get_cl_interval(lmax=lmax, nside=nside) else cl = interval
interval=cl

;print, cl
;print, nside
n = n_elements(cl)
lmax= cl(n-1)+1
; print, cl

res = dblarr(lmax,n+1)
final = dblarr(lmax,n+1)
win = dblarr(n+1,2)

;compute_h,cl(0),0.5,h  ;--- Spline
; print, cl(0)
hgmey,cl(0),0.5, h, g  ;--- MeyerWave

res[0:cl[0],0] = h[0:cl[0]]
win(0,0) = 0
win(0,1) = cl[0]
win(1,0) = cl[0]

for i = 1,n-1 do begin
	;compute_h,cl(i)-cl(i-1), 0.5, h  ;--- Spline
	hgmey,cl(i)-cl(i-1),0.5, h, g  ;--- MeyerWave
	; help, h
	res[cl[i-1]:cl[i],i] = h[0:cl[i]-cl[i-1]]
	res[0:cl[i-1],i] = 1
endfor

res[*,n]=1
final[*,0] = res[*,0]
for i= 0, n-1  do begin
	final(*,i+1) = res(*,i+1)-res(*,i)
endfor
;print, n
;help, final
hfilter=res[*,0:n-1]

for i= 2, n-1  do begin
	win(i,0) = cl[i-2]
	win(i,1) = cl[i]
endfor

win(0,0) = 0
win(0,1) = cl[0]
win(1,0) = cl[0]
win(1,0) = 0
win(1,1) = cl[1]
if n GT 2 then win(n-1,0) = cl[n-3]
win(n-1,1) = cl[n-1]
win(n,0) = cl[n-2]
win(n,1) = cl[n-1]
; help, hfilter
end

;=================================================================


;=================================================================

;ondelette "haar"
pro get_wp_box_filter, final, interval=interval, win
; interval = input keyword intarray [NScale]. interval[j] = lmax of scale j
;            LMAX = MAX( interval )
;            NbrScale = (size(interval)[1])
;            By default, NbrScale=22 
;                        LMAX = 3000
; final = output array [0:LMAX, NbrScale]
 ; win = intarray [NbrScale,2] :  The WP scale j is defined between  win[j,0] and  win[j,1]

cl_fine =[   64,128,192,256,320,384,512,1024,1280,1536,1664,1792,1920,2048,2176,2304,2432,2560, 2688, 2816,2944,3000]
if not keyword_set(interval) then cl = cl_fine else cl= interval
;cl = [0,200,400,600,800,990]

cl0 = cl[0]
n =5
n = (size(cl))(1)
lmax= cl(n-1)+1

res = fltarr(lmax,n+1)
final = fltarr(lmax,n)
final_tot = fltarr(lmax)
win = fltarr(n,2)

compute_h,cl0,0.5,h
n_h = n_elements(h)
bande = dblarr(n_h*2)
bande[0:n_h-1] = reverse(h)
bande[n_h:2*n_h-1] =h
;help,bande
final[0:cl[0],0] = h[0:cl[0]]
win(0,0) = 0
win(0,1) = cl[0]
win(1,0) = cl[0]

for i = 0,n-1 do begin
mini = max([0,cl(i)-n_h])
maxi = min([max(cl),cl(i)+n_h-1])
;print,mini,maxi
mini2 = max([0,n_h-cl(i)])
maxi2 = min([2*n_h-1,max(cl)-cl(i)+n_h])
;print,mini2,maxi2
;help,final[mini:maxi,i] 
;help,bande[mini2:maxi2]

final[mini:maxi,i] = bande[mini2:maxi2]
win[i,0] = mini
win[i,1] = maxi+2
final_tot = final_tot + final(*,i)

;oplot,final(*,i)
endfor

;plot,final_tot

for i = 0,n-1 do begin
final(*,i) = final(*,i) / final_tot
;oplot,final(*,i)
endfor
;window,/free
;plot,final(*,0)

finalt_tot = final(*,0)
for i = 1,n-1 do begin
final_tot = final_tot + final(*,i)

;oplot,final(*,i)

endfor
;oplot,final_tot
;info,final_tot

end

;=================================================================
; ondelette overlapping using spline wavelet

pro get_overlapping_spine_filter, final, interval=interval, win
;Figure partitionCl.jpg:
cl_fine =[64,128,192,256,320,384,512,1024,1280,1536,1664,1792,1920,2048,2176,2304,2432,2560, 2688, 2816,2944,3000]
 
cl = cl_fine;128
cl_fine =[  [ 0,64,128,192,256,320,384,512,768],$
            [ 64,128,192,256,320,384,512,768,995] ]
 if not keyword_set(interval) then cl = cl_fine else cl= interval

cl0 = 300
n =5
n = (size(cl))(1)
lmax= cl(n-1)+1

res = fltarr(lmax,n+1)
final = fltarr(lmax,n-1)
final_tot = fltarr(lmax)
win = fltarr(n,2)

compute_h,cl0,0.5,h
n_h = n_elements(h)
bande = dblarr(n_h*2)
bande[0:n_h-1] = reverse(h)
bande[n_h:2*n_h-1] =h
win(0,0) = 0
win(0,1) = cl[0]
win(1,0) = cl[0]

for i = 0,n-2 do begin
mini = cl[i]
maxi = cl[i+1]-1
final[mini:maxi,i] = 1
win[i,0] = mini
win[i,1] = maxi+2
final_tot = final_tot + final(*,i)
endfor

plot,final_tot

for i = 0,n-1 do begin
;final(*,i) = final(*,i) 
oplot,final(*,i)
endfor
end


;=================================================================

pro mrs_wptrans, Imag, out, NbrScale=NbrScale
COMMON C_PLANCK

if N_PARAMS() LT 2  then begin 
        print, 'CALLING SEQUENCE: mrs_wptrans, Imag, Trans, NbrScale=NbrScale'
        goto, DONE
        end

if type_code(Imag) EQ 8 then begin
	GLESP = 1  
	Healpix_with_Glesp = 0
end else begin 
	GLESP = 0
	Healpix_with_Glesp = 0
endelse

if GLESP EQ 0 then begin
	npix = (size(imag))[1]
	nside = npix2nside(npix)
	;if not keyword_set(lmax) then 
	lmax = nside *3
	if lmax GT P_LMAX  then lmax = P_LMAX
	nx = 0
	np = 0
	x_sky = 0
	y_sky = 0
end else begin
	npix = (size(imag.t_sky))[1]
	nside = 0
	nx = imag.nx
	np = imag.np
	x_sky = imag.x_sky
	y_sky = imag.y_sky
	;if not keyword_set(lmax) then 
	lmax = min([(nx-1)/2,np/4])
end

get_wp_meyer_filter, Wave, lwin, lmax=lmax, nside=nside, hfilter=hfilter

n_wave = (size(Wave))(2)

Hscale = imag
 
mrs_almtrans, imag, ALM, lmax=lmax
ALM_HighResolImag = ALM.alm
 
TabWavelet = dblarr(12L*nside*nside, N_wave)
if GLESP EQ 0 then TabWavelet[*,0] = Imag else TabWavelet[*,0] = Imag.t_sky

for j=0,n_wave-2 do begin
   h = hfilter(*,n_wave-2-j)
   if j NE 0 then ALM.alm = ALM_HighResolImag
   alm_product2, ALM_HighResolImag, h, alm_h  
   ALM.alm = alm_h
   ALM.lmin = 0
   ALM.lmax = lwin[n_wave-2-j,1]
   mrs_almrec, ALM, LScale

	if GLESP EQ 0 then begin
		TabWavelet[*,j] = TabWavelet[*,j] - LScale
		TabWavelet[*,j+1] = LScale
	end else begin
		TabWavelet[*,j] = TabWavelet[*,j] - LScale.t_sky
		TabWavelet[*,j+1] = LScale.t_sky
	end
   
   if (j EQ 0) then begin
    TabFilterH = fltarr( N_ELEMENTS(h), n_wave-1)
    TabFilterG = fltarr( N_ELEMENTS(h), n_wave-1)
    TabPhi = fltarr( N_ELEMENTS(h), n_wave-1)
    TabPsi = fltarr( N_ELEMENTS(h), n_wave-1)
    TabFilterH[*,j] = h
    TabFilterG[*,j] = 1. - h
    TabPhi[*,j] = h
    TabPsi[*,j] =  1. - h
  end else begin
    TabFilterH[*,j] = h
    TabFilterG[*,j] =  TabFilterH[*,j-1] - h
    TabPhi[*,j] =  h
    TabPsi[*,j] = TabFilterH[*,j-1] - h 
  end

endfor
meyerwave =0
difinsh = 0
pyrtrans=0
Healpix_with_Glesp =0
TabNorm=[0.85,0.12,0.046,0.0224606,0.011,0.006]
out = {UseGLESP: GLESP, NbrScale : N_wave, nside : nside, nx: nx, np:np, npix:npix, Coef : TabWavelet, lmax:lmax, MeyerWave:MeyerWave, $
       TabFilterH:TabFilterH, TabFilterG:TabFilterG, TabPhi:TabPhi, TabPsi:TabPsi, $
       DifInSH:DifInSH, pyrtrans:pyrtrans, x_sky:x_sky,  y_sky :y_sky, TabNorm:TabNorm, Healpix_with_Glesp: Healpix_with_Glesp, Wave:Wave, NeedletWave:0,  B_NeedletParam:2}
NbrScale = N_wave
DONE:
end


;=================================================================

 
