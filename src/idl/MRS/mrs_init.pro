;+
; NAME:
;        mrs_init
;
; PURPOSE:
;	The MultiResolution on the Sphere (MRS) IDL code can work only if the file mrs_init
;       has been compile (.r mrs_init). All the MRS routine work with NESTED online HEALPIX
;       maps.
;
;       This file contains several routines
;
;			function getbeam, Fwhm=Fwhm, lmax=lmax -- Return a gaussian beam of size lmax+1 (lmax=400 by default) and Fwhm arc min (Fwhm=10 arc min by default).
;
;			function getidealbeam, beamdata, lmin=lmin, lmax=lmax, tozero=tozero -- From a beam beamdata, creates and retrun an ideal beam of same size, defined by lmin, lmax and option tozero.
;
;			function pixel_size, nside -- Return the pixel size in arc minutes of a healpix map with parameter nside.
;
;			function l2amin, l -- Convert a l parameter into an arc minute value.
;			function amin2l, a -- Convert an arc minute value into a l parameter.
;
;			function gettmpfilename -- Generates a random fits file name.
;
;			function get_dircontent, Dir, suffix=suffix -- Get the contents of directory Dir and return the results in a string table (Unix command ls -l).
;
;			function wts, Imag, keywords... -- Call mrs_wttrans, Imag, trans or mrs_pwttrans, Imag, trans and return the output structure trans.
;			function iwts, Trans, NbrScale=NbrScale, filter=filter -- Call inverse transform and return the reconstructed image.
;			function wtsf, Imag, keywords... -- Call mrs_wtfilter, Imag, Filter and return the filtered image Filter.
;
;			pro plotcl, Cl, tit=tit, indl=indl, xrange=xrange, yrange=yrange, line=line, thick=thick, nonorm=nonorm, lnorm=lnorm, color=color, background=background, err=err, xlog=xlog, nodata=nodata
;					-- Plot in a multipole l*(l+1) scheme
;			pro oplotcl, Cl,  indl=indl, line=line, thick=thick, nonorm=nonorm, lnorm=lnorm, color=color -- Overplot function
;
;			pro plotAm, a, l, log=log -- Plot the alm coefficients for a given l value. a is a structure, see mrs_almtrans.pro
;
;			pro tab2nest, nside, in, out
;
;			pro put_all_faces, CubeFace, HealpixIma -- Put the 12 faces CubeFace[*,*,12] into a healpix NESTED map
;			function f2h, CubeFace
;			pro get_all_faces, Imag, CubeFace -- Extract the 12 faces from an heapix NESTED map and put them in CubeFace[*,*,12]
;			function h2f, Imag
;			function get_one_face, Imag, NumFace -- Extract one of the 12 faces from an healpix NESTED map
;			pro put_one_face, Imag, Face, NumFace -- Put one of the 12 faces into a healpix NESTED map
;
;			pro softthreshold, Data, Lambda -- Aply a soft thesholding to the floatarray Data with threshold level Lambda  
;
;			pro tvso, Data -- Run the visualization orthview command for a healpix NESTED map
;
;			pro spline, size, tab
;			pro spline2, size, l, lc, tab
;
;			pro sym, tab, tab_sym -- From an array tab of size n, create an array tab_sym of size 2*n+1 and symmetric structure.
;
;			pro compute_h, size, lc, h
;			pro compute_htilde, nlmax, ech, htilde
;			pro compute_g, size, lc, g
;			pro compute_gtilde, nlmax, ech, gtilde
;
;			pro alm_product2, alm1, al0, result -- Compute product of alm floatarray with al0
;
;			function fctlow, npix
;
;			pro hgmey, npix1, scale, h, g, dif=dif
;
;			pro rotate_map_nest, in, a1, a2, a3, out
;
;			pro text2map, in, out
;
;			pro mak_map, n, out2, t_interpol = t_interpol
;
;			pro index2lm2, index, l, m -- From an index in a ALM list array, get the corresponding l and m values.
;			pro lm2index2, l, m, index -- From l and m value, get the index in ALM list array
;
;			pro tab2alm, res, alm, complex=complex -- Convert ALM coefficients sorted in table form into coefficients sorted in list array form with real and immaginary part as components. 
;			pro tab2alm_pola, res, alm, complex=complex -- Same but for polarized ALM T, E and B
;			pro alm2tab, lm, res, complex=complex, TabNbrM=TabNbrM, NbrL=NbrL -- Convert ALM coefficients sorted in list array form with real and immaginary part as components into coefficients sorted in table form.
;			pro alm2tab_complex_in, lm, res, complex=complex, TabNbrM=TabNbrM, NbrL=NbrL -- Convert ALM coefficients sorted in list array form with complex values into coefficients sorted in table form.
;			pro alm_pola2tab, lm, res, complex=complex, TabNbrM=TabNbrM, NbrL=NbrL -- Same but for polarized ALM T, E and B
;			pro alm_pola2tab_complex_in, lm, res, complex=complex, TabNbrM=TabNbrM, NbrL=NbrL -- Same but for polarized ALM T, E and B
;
;			pro ismr1 -- Check if MR1 package is installed.
;
;			function mrs_variance_stabilization, PowSpec, mu=mu, psi1=psi1, FirstL=FirstL -- Variance stabilisation of a power spectrum
;			function mrs_variance_l1_stabilization, PowSpec, mu, psi1
;			function mrs_inv_variance_stabilization, StabPowSpec, mu, psi1 -- Inverse variance stabilisation
;			function mrs_inv_variance_l1_stabilization, StabPowSpec, mu, psi1
;
;			pro reim2mp, re, im, m, p -- Real part, immaginary part -> modulus, phasis
;			pro reim2pp, re, im, m, p -- Real part, immaginary part -> squared modulus, phasis
;			pro mp2reim, m, p, re, im -- Modulus, phasis -> real part, immaginary part
;			pro pp2reim, m, p, re, im -- Squared modulus, phasis -> real part, immaginary part
;
;			pro pnside, data -- Print the nside value of the healpix map data.
;			function gnside, data -- Return the nside value of the healpix map data.
;
;			pro mrs_info, data, mes=mes -- Print statistics of the map data.
;			function mrs_sigma, data -- Return the standard deviation of the map data.
;			function mrs_max, data -- Return the maximum of the map data.
;			function mrs_absmax, data -- Return the maximum of the absolute value of the map data.
;			function mrs_min, data -- Return the minimum of the map data.
;			function mrs_mean, data -- Return the mean of the map data.
;			pro mrs_set, data, Value -- Set the whole map data to the constant value.
;			function mrs_diff, d1, d2 -- Return the difference between the two maps d1 and d2, pixel by pixel.
;			function mrs_mult, d1, d2 -- Return the product between the two maps d1 and d2, pixel by pixel.
;			function mrs_add, d1, d2 -- Return the sum between the two maps d1 and d2, pixel by pixel.
;			function mrs_getpix, Ima -- Return the full data array of Ima.
;			pro mrs_putpix, Ima, pix -- Copy the full data array pix into Ima.
;			function mrs_mad, Ima -- Return the median of the absolute value of the map.
;
;			function mrs_absthreshold, Ima, T, soft=soft, l2=l2 -- Return the result of a hard or soft thresholding of the image.
;
;			pro mrs_pos, Ima, T=T -- In the image Ima, set pixel lower than T to zero, T >= 0 default is T = 0
;
;			pro pixf2pix, nside, x, y, face, ipix -- Convert pixel position in a face into pixel index in Healpix list. 
;			pro pix2pixf, nside, ipix, x, y, face -- Convert pixel index in Healpix list into pixel position in a face and face number.
;			pro pixf2ang, nside, x, y, face, theta, phi -- Convert pixel position in a face into spherical coordinates in radian.
;			pro ang2pix, nside, theta, phi, x, y, face -- Convert spherical coordinates in radian into pixel position in a face and face number.
;			pro ang2lb, theta, phi, l, b -- Convert spherical coordinates in radian into latitude and longitude coordinates in degree.
;			pro lb2ang, l, b, theta, phi -- Convert latitude and longitude coordinates in degree into spherical coordinates in radian.
;			pro pix2lb, nside, ipix, l, b -- Convert pixel index in Healpix list into latitude and longitude coordinates in degree.
;			pro lb2pix, nside, l, b, ipix -- Convert latitude and longitude coordinates in degree into pixel index in Healpix list.
;			pro ang2radec, theta, phi, ra, dec, year=year, degree=degree -- Convert spherical coordinates in radian into galactical coordinates.
;			pro radec2ang, ra, dec, theta, phi, year=year, degree=degree -- Convert galactical coordinates into spherical coordinates in radian.
;
;
; HISTORY:
;	Written: Jean-Luc Starck and Pierrick Abrial, 2005
;	February, 2005 File creation
;-
;===============================================================

pro wadel
while !D.window NE -1 do wdelete
end
;===============================================================


function tol, Map, Lmax, lmin=lmin
mrs_almtrans, Map, a, lmax=lmax, /tab
if keyword_set(lmin) then a.alm[0:lmin-1, *,*] = 0.
mrs_almrec, a, Rec
return, rec
end

;=============================================
 
function toquad, Map
Lmin=2
Lmax=2
return, tol(Map, Lmax, lmin=lmin)
end

;=============================================
 
function tooct, Map
Lmin=3
Lmax=3
return, tol(Map, Lmax, lmin=lmin)
end

;=============================================

function getbeam, Fwhm=Fwhm, lmax=lmax
if not keyword_set(lmax) then lmax=4000
if not keyword_set(Fwhm) then Fwhm=10.  ; 10 arc minutes

F = Fwhm / 60. * !dtor
l = findgen(lmax+1)
ell = l*(l+1)
bl = exp(-ell*F^2. /16./alog(2.))
  
return, bl
end


;================================================================

function getdirac, nside=nside
if not keyword_set(nside) then nside=256
f = fltarr(nside, nside,12)
f[nside/2,nside/2, 4] = 1.
d = f2h(f)
return, d
end

;=============================================
; Create an ideal beam 
function getidealbeam, beamdata, lmin=lmin, lmax=lmax, tozero=tozero
bl = float(beamdata)
bl[0:lmin] = 1.
Np = lmax-lmin
x = findgen(Np) / float(Np-1)*!PI/2
; bl[lmin:lmax] = beamdata[lmin:lmax] + cos(x) * (1. - beamdata[lmin:lmax])
spline2, Np, 1, 1, t 
if keyword_set(tozero) then begin
   bl[lmin:lmax] = t 
   bl[lmax:*]=0.
end else bl[lmin:lmax] = beamdata[lmin:lmax] + t * (1. - beamdata[lmin:lmax])
return, bl
end

;===============================================================

function pixel_size, nside
; Return the pixel size of a healpix map in arc minutes
; SKI_SURFACE IN SQUARE DEGREES =  4. * !PI * (360. / (2*!PI))^2 = 41253
 psize = 41253. / (float(nside)^2.*12.) * 60.^2.
 return, sqrt(psize)
 end
 
;===============================================================

function l2amin, l
a = 1. / l
a  =  a * 180.* 60. / !PI
return, a
end

function amin2l, a
ar =  a / (180.* 60.) * !PI
l = 1. / ar
return, l
end

;===============================================================

function get_dircontent, Dir, suffix=suffix
 t = gettmpfilename()
 if not keyword_set(suffix) then cmd = 'ls -1 ' + Dir + '*  > ' + t $
 else cmd = 'ls -1 ' + Dir + '*' + suffix + '  > ' + t
 spawn, cmd
 readcol, t, TabFn, format='a'
 return, TabFN
 delete, t
end

;===============================================================
 
function wts, Imag, Pyr=Pyr, NbrScale=NbrScale, lmax=lmax, DifInSH=DifInSH, MeyerWave=MeyerWave 
if keyword_set(pyr) then mrs_pwttrans, Imag, out, NbrScale=NbrScale, lmax=lmax, DifInSH=DifInSH, MeyerWave=MeyerWave $
else  mrs_wttrans, Imag, out, NbrScale=NbrScale, lmax=lmax, DifInSH=DifInSH, MeyerWave=MeyerWave
return, out
end

;===============================================================

function iwts, Trans, NbrScale=NbrScale, filter=filter
if Trans.pyrtrans eq 1 then mrs_pwtrec, Trans, Imag, filter=filter $
else mrs_wtrec, Trans, Imag,  filter=filter
return, Imag
end

;===============================================================

function wtsf, Imag, NbrScale= NbrScale, NSigma=NSigma, SigmaNoise=SigmaNoise, $
                  mad=mad, KillLastScale=KillLastScale,Trans=Trans, Pyr=Pyr
mrs_wtfilter, Imag, Filter, NbrScale= NbrScale, NSigma=NSigma, SigmaNoise=SigmaNoise, $
                  mad=mad, KillLastScale=KillLastScale,Trans=Trans, Pyr=Pyr
return, Filter
end

;===============================================================

pro plotcl, Cl, tit=tit, indl=indl, xrange=xrange, yrange=yrange, line=line, thick=thick, nonorm=nonorm, lnorm=lnorm, color=color, background=background, err=err, xlog=xlog, nodata=nodata, ylog=ylog
if not keyword_set(tit) then tit='TEMPERATURE'
vs = size(cl)
ncols = vs[1]
;screen_size = get_screen_size()
;window, /free, title=fitsfile, xs=600<screen_size[0], ys=((ncols+1)/2)*300 < screen_size[1]
 nextra = 0
 l = findgen(ncols)
 xtitle = '!6 Multipole !12 l !6'
 fl2 = l*(l+1.)/(2.*!pi)
 ytitle = '!12  l (l+1) !6C!12!dl!n /!6 2!7p!n!6 '
 
 lnorm = fl2   
 if keyword_set(nonorm) then fl2 = 1.
 if keyword_set(nonorm) then print, "nonorm"
 if keyword_set(xlog) and not keyword_set(xrange) then xrange=[1,(size(cl))[1]]
 if not keyword_set(err) then plot,l,fl2*Cl, thick=thick, $
      xtitle=xtitle,ytitle=ytitle, title='!6'+tit, charsize=1.3, xrange=xrange, yrange=yrange, line=line, color=color, background=background, xlog=xlog, nodata=nodata, ylog=ylog $
 else ploterror,l,fl2*Cl, err*fl2, thick=thick, $
      xtitle=xtitle,ytitle=ytitle, title='!6'+tit, charsize=1.3, xrange=xrange, yrange=yrange, line=line, color=color, background=background, xlog=xlog, nodata=nodata, ylog=ylog
 indl=l
 oplotcl, cl*0
end

;===============================================================

pro oplotcl, Cl,  indl=indl, line=line, thick=thick, nonorm=nonorm, lnorm=lnorm, color=color
vs = size(cl)
ncols = vs[1]
 nextra = 0
 l = findgen(ncols)
 indl=l
 fl2 = l*(l+1.)/(2.*!pi)
 lnorm=fl2
 if keyword_set(nonorm) then fl2=1.    
 oplot,l,fl2*Cl, thick=thick, line=line, color=color
 l = fl2
end



;===============================================================

pro plotAm, a, l, log=log
Nm = a.tabnbrm[l]-1
if keyword_set(log) then plot, alog(a.alm[l, 0:Nm, 0]) $
else plot,  a.alm[l, 0:a.tabnbrm[l]-1,0]
end


;===============================================================

pro tab2nest, nside, in, out
i = long(in mod (nside/2))
j = long(in / (nside/2))
pas = long(1)
out = long(0)
nside2 = long(nside)
while (nside2 gt 2) do begin
out = long(out) +(i mod 2)*pas + 2 *pas * ( j mod 2)
;print ,out
i = i/2
j= j/2
nside2 = nside2/2
pas = pas*4
endwhile
;print,in
;print,'out ',out
;print,in/2 
end

;===============================================================

pro put_all_faces, CubeFace, HealpixIma
vs = size(CubeFace)
; N = vs[1]
;CubeFace(*)=0
;for i=0,11 do   CubeFace[*,N/2,i] = i+1

nside = vs[1]
taille = long(12.*float(nside)^2)
nside = npix2nside(taille)
HealpixIma = fltarr(taille)
taille2 = taille/12l
index = lindgen(taille2)
cote = sqrt(taille2)
tab2nest,nside*2,index,index2
;window,0
;plot,index
;window,1
;plot,index2
for i=0,11 do begin
F = CubeFace[*,*,i]
; put_face, HealpixIma, F, nside, i+1
HealpixIma[index2+taille2*i] = F[*]
end
; tvs, HealpixIma 
end
;===============================================================

; FACE to Healpix image
function f2h, CubeFace 
 put_all_faces, CubeFace, HealpixIma
return, HealpixIma
end

;===============================================================

pro get_all_faces, Imag, CubeFace

taille = (size(imag))[1]
nside = npix2nside(taille)
; print,nside
taille2 = taille/12l
index = lindgen(taille2)
cote = sqrt(taille2)
tab2nest,nside*2,index,index2
out = fltarr(taille)
out = reform(out,12,taille2)

CubeFace = fltarr(cote,cote,12)
N = fix(sqrt(taille2))

for i=0,11 do begin
    text = imag[index2+taille2*i]
    text = reform(text,long(cote),long(cote) )
    ; text = extract_face(Imag,nside,i+1)  
    ;load, text
    ;pause, 1
    CubeFace[*,*,i] = text
end

end
;===============================================================

function get_one_face, Imag, NumFace
taille = (size(imag))[1]
nside = npix2nside(taille)
; print,nside
taille2 = taille/12l
index = lindgen(taille2)
cote = sqrt(taille2)
tab2nest,nside*2,index,index2
out = fltarr(taille)
out = reform(out,12,taille2)

Face = fltarr(cote,cote)
N = fix(sqrt(taille2))
i = NumFace
Face[*,*]  =  reform( imag[index2+taille2*i], long(cote),long(cote) )
return, Face
end

;===============================================================

pro put_one_face, Imag, Face, NumFace
nside = gnside(imag)
taille = long(12.*float(nside)^2)
taille2 = taille/12l
index = lindgen(taille2)
cote = sqrt(taille2)
tab2nest,nside*2,index,index2
Imag[index2+taille2*NumFace] = Face[*]
end


;===============================================================
; Healpix to Faces cube
function h2f, Imag
get_all_faces, Imag, CubeFace
return, CubeFace
end

;===============================================================

pro softthreshold, Data, Lambda
if N_ELEMENTS(Lambda) NE N_ELEMENTS(Data) then begin
SoftThresold = dblarr(N_ELEMENTS(Data))
SoftThresold[*] = Lambda[0]
end else SoftThresold = Lambda

ind = where ( Data GT 0, c)
if c GT 0 then begin
   Pos = Data[ind] - SoftThresold[ind]
   indN = where (Pos LT 0, c)
   if c GT 0 then Pos[indN] = 0
   Data[ind] = Pos
end 
ind = where ( Data LT 0, c)
if c GT 0 then begin
   Neg = Data[ind] + SoftThresold[ind]
   indP = where (Neg GT 0, c)
   if c GT 0 then Neg[indP] = 0
   Data[ind] = Neg
end 
end



;===============================================================

;pro tvs, Data, graticule=graticule, png=png, TITLE=TITLE, COLT=COLT, NOBAR=NOBAR
;mrs_tv, Data, graticule=graticule, png=png, TITLE=TITLE, COLT=COLT, NOBAR=NOBAR
;end

;===============================================================

pro tvso, Data
orthview, Data, /online, /nested
end
;===============================================================

pro spline, size, tab
res  =findgen((2*size)+1)

res =res -size
res = 2.0 * res/ size
;res = 4.0 * res / size
tab = 1.0 /12.0 * (( abs(res-2))^3 - 4.0* (abs(res-1))^3 + 6 *(abs(res))^3 - 4.0 *( abs(res+1))^3+(abs(res+2))^3)

; window,/free
; plot,res,tab
;tab2 = 1/sqrt(2)*exp(-res^2)

;res2= res*2.0
;tab3 = 2*1.0 /12.0 * (( abs(res2-2))^3 - 4.0* (abs(res2-1))^3 + 6 *(abs(res2))^3 - 4.0 *( abs(res2+1))^3+(abs(res2+2))^3)
;tab4 = tab3  - tab

;plot,res,tab,yrange= [-0.5,2]
;oplot,res,tab2,color=3
;oplot,res,tab3,color=4
;oplot,res,tab4,color=5

;window,1
;wset,1
;plot,res,tab3,title='range/2'

;window,2
;wset,2
;plot,res,tab4,title='phi  - phi/2'


;print,total(tab4)
;Splot ,res
end

;==================================
pro sym, tab, tab_sym
size=size(tab)
size=size[1]
tab=reform(tab)
tab_sym=dblarr(2*size-1)
for j=0,size-1 do begin
	tab_sym(j)=tab(size-1-j)
	tab_sym(j+size-1)=tab(j)
endfor
end

;==================================

pro spline2, size, l, lc, tab
res =dindgen(size+1)

res = 2.0 * l * res / (lc *size)
tab = (3.0/2.0)*1.0 /12.0 * (( abs(res-2))^3 - 4.0* (abs(res-1))^3 + 6 *(abs(res))^3 - 4.0 *( abs(res+1))^3+(abs(res+2))^3)

;plot, tab

;window,1
;wset,1
;plot,res,tab,title='spline'


end

;==================================

pro compute_h, size, lc, h

spline2,size,2.*lc,1,tab1
spline2,size,lc,1,tab2

;plot, alog(abs(tab1))

h = tab1/(tab2+0.000001)
h[size/(2.*lc):size]=0.
;plot,h
end

;==================================

pro compute_htilde, nlmax, ech, htilde

compute_g,nlmax,ech,g
compute_h,nlmax,ech,h

htilde = h / (g^2+h^2)
htilde[nlmax/(2.*ech):nlmax]=0
end
;===
pro compute_filter_htilde,size,lc,filter_htilde
compute_filter_g,size,lc,filter_g
compute_filter_h,size,lc,filter_h
filter_htilde = filter_h / (filter_g^2+filter_h^2)
end
;===
pro compute_filter_gtilde,size,lc,filter_gtilde
compute_filter_g,size,lc,filter_g
compute_filter_h,size,lc,filter_h
filter_gtilde = filter_g / (filter_g^2+filter_h^2)
end
;===
PRO compute_filter_h,taille,lc,filter_h
filter_h=dblarr(taille)
parite=0
if ((taille)/2 eq float(taille)/2.0) then parite=1    ;test de parite
if (parite eq 1) then begin

	compute_h,taille/2,lc,h
	sym,h,hsym
	filter_h(0:taille/2)=h
	filter_h(taille/2+1:taille-1)=hsym(1:(taille/2)-1)
endif else begin
	compute_h,taille/2,lc,h
	sym,h,hsym
	filter_h(0:taille/2)=h
	filter_h(taille/2+1:taille-1)=hsym(0:taille/2-1)
endelse
END
;===
PRO compute_filter_g,taille,lc,filter_g
filter_g=dblarr(taille)
parite=0
if ((taille)/2 eq float(taille)/2.0) then parite=1    ;test de parite
if (parite eq 1) then begin
	compute_g,taille/2,lc,g
	sym,g,gsym
	filter_g(0:taille/2)=g
	filter_g(taille/2+1:taille-1)=gsym(1:(taille/2)-1)
endif else begin
	compute_g,taille/2,lc,g
	sym,g,gsym
	filter_g(0:taille/2)=g
	filter_g(taille/2+1:taille-1)=gsym(0:taille/2-1)
endelse

END

;==================================

pro compute_g, size, lc, g

spline2,size,2.*lc,1,tab1
spline2,size,lc,1,tab2


g = (tab2-tab1)/(tab2+0.000001)
g[size/(2.*lc):size]=1
;plot,g
end

;==================================

pro compute_gtilde, nlmax, ech, gtilde

compute_g,nlmax,ech,g
compute_h,nlmax,ech,h

gtilde = g / (g^2+h^2+0.0001)
gtilde[nlmax/(2.*ech):nlmax]=1
end

;==================================


; compute product of alm with al0.

pro alm_product2, alm1, al0, result

; print,'produit de alm. al0 doit etre etre a symetrie azimutal, ie alm = 0 pour m!=0
taille= size(al0)
taille = taille[1]
pi = !dpi
taille2 = (size(alm1))(1)
result = fltarr(taille2,2)

    j=0l
    for i=0l,taille-1 do begin 
      tmpval = al0[i,0]
      ;for j=((i)*(i+1))/2,((i+1)*(i+2))/2-1 do begin
      for j=((i)*(i+1))/2,(((i+1)*(i+2))/2)-1 do begin   
	 ; Calcul 1/l ou l/1 ??
;         result [j,0] = 2.0*pi*sqrt((4.0*pi)/((2.0*i)+1.0)) * tmpval * alm1[j,0]
;         result [j,1] = 2.0*pi*sqrt((4.0*pi)/((2.0*i)+1.0)) * tmpval * alm1[j,1]
         result [j,0] =  tmpval * alm1[j,0]
         result [j,1] =  tmpval * alm1[j,1]
        
	
	
	endfor
    endfor
;plot,result[*,0]
end

;==================================

function fctlow, npix 
 h = dblarr(npix)
 H(*) = 1
 n = double(npix)/2.
 i = indgen(Npix)
 x = double(i) ; double(i) - npix/2
 r = (ABS(x) - N) / N
 ind = where (r le 0,c)
 if c gt 0 then H(ind) = 1
 ind = where (r ge 1, c)
 if c gt 0 then H(ind) = 0
 ind = where ( (r GT 0) and (r LT 1), c)
 if c GT 0 then begin
         xw = 1. - r(ind)
	 lw = exp(1. -1./(1-exp(1.-1./(1-xw))))
         rw = exp(1. -1./(1-exp(1.-1./xw)))
         norm = double(sqrt(lw*lw+rw*rw))
         lw = lw/  norm
         H(ind) = lw
	 ;info, h
	end 
;plot, h
;info, h
return, h
end

;==================================
 
pro hgmey, npix1, scale, h, g, dif=dif
  npix = 2 * npix1/2 + 1
  h = dblarr(npix)
  n1 = npix / 2.^(scale+1)
  hl = fctlow(n1)
  h(0:n1-1) = hl
  if keyword_set(dif) then g = 1d - h $
  else g = sqrt( 1d - h^2)
  ;plot, h
  ;oplot, g
end


;==================================

pro rotate_map_nest, in, a1, a2, a3, out

nb_pix = (size(in))[1]

nside = npix2nside(nb_pix)
eul = euler_matrix_new(a1,a2,a3)

t = lindgen(nb_pix)
pix2vec_nest,nside,t,vect
vect2 = rotate_coord(vect,Euler_Matrix=eul)
vec2pix_nest,nside,vect2,out
out = in[out]

end

;==================================

pro text2map, in, out
taille = (size(in))[2]
nside = long((sqrt(taille)))
index = lindgen(taille)
tab2nest,nside*2,index,index2
out = fltarr(nside2npix(nside))
t= indgen(12)
for i=0,11 do begin
out[(long(taille)*i)+index2] = reform(in[i,*],taille) ;+index2

 endfor
end

;==================================
 
pro mak_map, n, out2, t_interpol = t_interpol
; t_interpol = 0 : out[i,j] = sqrt((i-mil)^2+(j-mil)^2)
; t_interpol = 1 : out[i,j] = min([abs(i-mil)/2,abs(i-mil)/2])abs(i-mil)/2+abs(j-mil)/2

if not keyword_set (t_interpol) then t_interpol = 0

n= long(n)
out = fltarr(n,n)

if t_interpol eq 0 then begin 
mil = (n-1.)/2.
for i=0l,n-1 do begin
  for j=0l,n-1 do begin
     out[i,j] = sqrt((i-mil)^2+(j-mil)^2)
  endfor
 endfor
out = 1 - (out/sqrt(2.*mil*mil))

endif else begin 
 mil = n/2.  
 for i=0l,n-1 do begin
    for j=0l,n-1 do begin
      out[i,j] =  max([abs(i-mil),abs(j-mil)])
 endfor
endfor
out = 1 - (out/mil)
endelse 

out2 = fltarr(12,n,n)
for i=0,11 do begin  
  out2[i,*,*] = out
  endfor
  out2 = reform(out2,12,n*n)

text2map,out2,maphealpixnested
out2 = maphealpixnested
end

;===========================================================
;===========================================================

pro index2lm2, index, l, m
; PB index=10619136  ==> m= -1  !
; 
l = long(sqrt(double(1l+8l*index))-1)/2
m = index - ((l*(l+1L))/2L)
end

pro lm2index2, l, m, index
index = ((long(l)*(l+1l)) / 2l) + m
end

;==========================================

function nbr_alm, Lmax
l = Lmax
N =  ((long(l)*(l+1l)) / 2l) 
return, N
end

;==========================================

pro tab2alm, res, alm, complex=complex 

vs = size(res)
taille = vs[1]
taille2 = nbr_alm(taille)
 
alm = dblarr(taille2,2)

if keyword_set(complex) then begin

for i=0l,taille2-1 do begin
  index2lm2,i,l,m
  alm[i,0] = real_part(res[l,m])
  alm[i,1] = imaginary(res[l,m])
endfor
endif else begin

for i=0l,taille2-1 do begin
  index2lm2,i,l,m
  alm[i,0] = res[l,m,0] 
  alm[i,1] = res[l,m,1] 
endfor
endelse
end

;==========================================

pro tab2alm_pola, res, alm, complex=complex

taille = size(res)
taille = taille[1]

lm2index2,taille,0,taille2

;index2lm2,taille,l,m
;print,l,m   
;l = l+1

alm = dblarr( taille2, 2, 3 )

if keyword_set(complex) then begin

	for i=0l,taille2-1 do begin

		index2lm2,i,l,m
		
		alm[i,0,0] = real_part( res[l,m,0] );	ALM T
		alm[i,1,0] = imaginary( res[l,m,0] )
		
		alm[i,0,1] = real_part( res[l,m,1] );	ALM E
		alm[i,1,1] = imaginary( res[l,m,1] )
		
		alm[i,0,2] = real_part( res[l,m,2] );	ALM B
		alm[i,1,2] = imaginary( res[l,m,2] )

	endfor
	
endif else begin

	for i=0l,taille2-1 do begin

		index2lm2,i,l,m
		
		alm[i,0,0] = res[l,m,0,0];	ALM T
		alm[i,1,0] = res[l,m,1,0]
		
		alm[i,0,1] = res[l,m,0,1];	ALM E
		alm[i,1,1] = res[l,m,1,1]
		
		alm[i,0,2] = res[l,m,0,2];	ALM B
		alm[i,1,2] = res[l,m,1,2]

	endfor

endelse
end

;==========================================

pro alm2tab, lm, res, complex=complex, TabNbrM=TabNbrM, NbrL=NbrL
taille = long(size(lm))
taille = taille[1]
; print,taille
index2lm2,taille,l,m
; print,l,m   
; l = l+1
TabNbrM = lonarr(l)
NbrL=long(l)-1

if keyword_set(complex) then begin
res = complexarr(l,l)
for i=0l,taille-1 do begin
  index2lm2,i,l,m
  if  l LT 0 or l GT NbrL then print, 'ERROR in alm2tab : NbrMaxL = ', NbrL, ' Current l = ', l
  res[l,m] = complex(lm[i,0],lm[i,1])
  TabNbrM[l] = TabNbrM[l] + 1
endfor
endif else begin
res = dblarr(l,l,2)

; help, TabNbrM
for i=0l,taille-1 do begin
  index2lm2,i,l,m
  if l LT 0 or l GT NbrL then print, 'ERROR in alm2tab : NbrMaxL = ', NbrL, ' Current l = ', l
  if m LT 0 or m GT NbrL then print, 'ERROR in alm2tab : NbrMaxL = ', NbrL, ' Current l = ', l, ' Current m = ', m
  res[l,m,0] = lm[i,0]
  res[l,m,1] = lm[i,1]
  TabNbrM[l] = TabNbrM[l] + 1
endfor
endelse
end

;==========================================

pro alm2tab_complex_in, lm, res, complex=complex, TabNbrM=TabNbrM, NbrL=NbrL

taille = long(size(lm))
taille = taille[1]
; print,taille
index2lm2,taille,l,m
; print,l,m   
; l = l+1
TabNbrM = lonarr(l)
NbrL=long(l)-1

if keyword_set(complex) then begin
res = complexarr(l,l)
for i=0l,taille-1 do begin
  index2lm2,i,l,m
  if  l LT 0 or l GT NbrL then print, 'ERROR in alm2tab : NbrMaxL = ', NbrL, ' Current l = ', l
  res[l,m] = lm[i]
  TabNbrM[l] = TabNbrM[l] + 1
endfor
endif else begin
res = dblarr(l,l,2)

; help, TabNbrM
for i=0l,taille-1 do begin
  index2lm2,i,l,m
  if l LT 0 or l GT NbrL then print, 'ERROR in alm2tab : NbrMaxL = ', NbrL, ' Current l = ', l
  if m LT 0 or m GT NbrL then print, 'ERROR in alm2tab : NbrMaxL = ', NbrL, ' Current l = ', l, ' Current m = ', m
  res[l,m,0] = real_part(lm[i])
  res[l,m,1] = imaginary(lm[i])
  TabNbrM[l] = TabNbrM[l] + 1
endfor
endelse
end

;==========================================

pro alm_pola2tab, lm, res, complex=complex, TabNbrM=TabNbrM, NbrL=NbrL

taille = long(size(lm))
taille = taille[1]

; print,taille
index2lm2,taille,l,m
; print,l,m   
; l = l+1
TabNbrM = lonarr(l)
NbrL=long(l)-1

if keyword_set(complex) then begin

	res = complexarr(l,l,3)
	
	for i=0l,taille-1 do begin

		index2lm2,i,l,m

		if l LT 0 or l GT NbrL then print, 'ERROR in alm_pola2tab : NbrMaxL = ', NbrL, ' Current l = ', l

		res[l,m,0] = complex( lm[i,0,0], lm[i,1,0] );ALM T
		
		res[l,m,1] = complex( lm[i,0,1], lm[i,1,1] );ALM E
		
		res[l,m,2] = complex( lm[i,0,2], lm[i,1,2] );ALM B
		
  		TabNbrM[l] = TabNbrM[l] + 1
  		
	endfor
	
endif else begin

	res = dblarr(l,l,2,3)

	; help, TabNbrM

	for i=0l,taille-1 do begin

		index2lm2,i,l,m
		
		if l LT 0 or l GT NbrL then print, 'ERROR in alm_pola2tab : NbrMaxL = ', NbrL, ' Current l = ', l
  		if m LT 0 or m GT NbrL then print, 'ERROR in alm_pola2tab : NbrMaxL = ', NbrL, ' Current l = ', l, ' Current m = ', m
  
  		res[l,m,0,0] = lm[i,0,0];ALM T
  		res[l,m,1,0] = lm[i,1,0]
  		
  		res[l,m,0,1] = lm[i,0,1];ALM E
  		res[l,m,1,1] = lm[i,1,1]
  		
  		res[l,m,0,2] = lm[i,0,2];ALM B
  		res[l,m,1,2] = lm[i,1,2]
  		
  		TabNbrM[l] = TabNbrM[l] + 1
	
	endfor

endelse

end

;==========================================

pro alm_pola2tab_complex_in, lm, res, complex=complex, TabNbrM=TabNbrM, NbrL=NbrL

taille = long(size(lm))
taille = taille[1]

; print,taille
index2lm2,taille,l,m
; print,l,m   
; l = l+1
TabNbrM = lonarr(l)
NbrL=long(l)-1

if keyword_set(complex) then begin

	res = complexarr(l,l,3)
	
	for i=0l,taille-1 do begin

		index2lm2,i,l,m

		if l LT 0 or l GT NbrL then print, 'ERROR in alm_pola2tab : NbrMaxL = ', NbrL, ' Current l = ', l

		res[l,m,0] = lm[i,0];ALM T
		
		res[l,m,1] = lm[i,1];ALM E
		
		res[l,m,2] = lm[i,2];ALM B
		
  		TabNbrM[l] = TabNbrM[l] + 1
  		
	endfor
	
endif else begin

	res = dblarr(l,l,2,3)

	; help, TabNbrM

	for i=0l,taille-1 do begin

		index2lm2,i,l,m
		
		if l LT 0 or l GT NbrL then print, 'ERROR in alm_pola2tab : NbrMaxL = ', NbrL, ' Current l = ', l
  		if m LT 0 or m GT NbrL then print, 'ERROR in alm_pola2tab : NbrMaxL = ', NbrL, ' Current l = ', l, ' Current m = ', m
  
  		res[l,m,0,0] = real_part(lm[i,0]);ALM T
  		res[l,m,1,0] = imaginary(lm[i,0])
  		
  		res[l,m,0,1] = real_part(lm[i,1]);ALM E
  		res[l,m,1,1] = imaginary(lm[i,1])
  		
  		res[l,m,0,2] = real_part(lm[i,2]);ALM B
  		res[l,m,1,2] = imaginary(lm[i,2])
  		
  		TabNbrM[l] = TabNbrM[l] + 1
	
	endfor

endelse

end

;==========================================

pro ismr1
COMMON MR1ENV
if keyword_set(mr1ok) then print, "MR1 OK for MRS" else print, "NO MR1 for MRS"
end

;==========================================

function mrs_variance_stabilization, PowSpec, mu=mu, psi1=psi1, FirstL=FirstL
if not keyword_set(FirstL) then FirstL=0

vs = size(PowSpec)
NbrAlm=vs[1]
Eps=double(1.e-15)
ind = where(PowSpec LT Eps, c)
if c GT 0 then PowSpec[ind] = 0.
spec1d = alog(PowSpec + Eps)
;PowSpec[0] = 1
;spec1d = alog(PowSpec)

 if not keyword_set(mu) or not keyword_set(psi1) then begin
   l = findgen(NbrAlm) + FirstL
   ; l = 2*l+1
   l = l + 0.5
   mr_prog, 'im1d_der_loggamma', l, f0, opt='-o 0'
   mr_prog, 'im1d_der_loggamma', l, f1, opt='-o 1'
   mu = f0 - alog(l)
   slog = spec1d
   psi1 = sqrt(f1)  ;  * sqrt(2.)
end
; mu[*]=0
spec1d = (spec1d  - mu) / psi1
return, spec1d
end

;==========================================

function mrs_variance_l1_stabilization, PowSpec, mu, psi1
Eps=double(1.e-15)
ind = where(PowSpec LT Eps, c)
if c GT 0 then PowSpec[ind] = 0.
spec1d = alog(PowSpec + Eps)
spec1d = (spec1d  - mu) / psi1
return, spec1d
end

;==========================================

function mrs_inv_variance_l1_stabilization, StabPowSpec, mu, psi1
Eps=double(1.e-15)
PowSpec = StabPowSpec * psi1  + mu
PowSpec = exp(PowSpec) - Eps
return, PowSpec 
end

;==========================================

function mrs_inv_variance_stabilization, StabPowSpec, mu, psi1
Eps=double(1.e-15)
PowSpec = StabPowSpec * psi1 + mu
PowSpec = exp(PowSpec) - Eps

return, PowSpec 
end

;==========================================

pro reim2mp, re, im, m, p
   m = sqrt(re^2 + im^2)
   P = atan(im, re)
end

pro reim2pp, re, im, m, p
   m = re^2 + im^2
   P = atan(im, re)
end

pro mp2reim, m, p, re, im 
  re = m * cos (p)
  im = m * sin(p)
end

pro pp2reim, m, p, re, im 
   re = sqrt(m) * cos (P)
   im = sqrt(m) * sin(P)
end

;==========================================

pro pnside, data
vs = size(data)
npix = vs[1]
nside = npix2nside(npix)
print, nside
end

;==========================================

function gnside, data
vs = size(data)
npix = vs[1]
nside = npix2nside(npix)
return, nside
end

;==========================================

pro mrs_info, data, mes=mes
if  type_code(Data) EQ 8 then d = data.T_Sky else d=data
m1 = min(d)
m2 = max(d)
m = mean(d)
if not keyword_set(mes) then mes= 'Stat: '
if m1 NE m2 then print, mes, ' Min = ', m1, ' Max = ', m2, ' Mean = ', m, ' Sigma = ', sigma(d) $
else  print, mes, ' Min = ', m1, ' Max = ', m2, ' Mean = ', m
end

;==========================================

function mrs_sigma, data
if  type_code(Data) EQ 8 then m = sigma(data.T_Sky) else m= sigma(data)
return,m 
end

;==========================================

function mrs_max, data
if  type_code(Data) EQ 8 then m = max(data.T_Sky) else m= max(data)
return,m 
end

;==========================================

function mrs_absmax, data
if  type_code(Data) EQ 8 then m = max(ABS(data.T_Sky)) else m= max(ABS(data))
return,m 
end

;==========================================

function mrs_min, data
if  type_code(Data) EQ 8 then m = min(data.T_Sky) else m= min(data)
return,m 
end

;==========================================

function mrs_mean, data
if  type_code(Data) EQ 8 then m = mean(data.T_Sky) else m= mean(data)
return,m 
end

;==========================================

pro mrs_set, data, Value
if  type_code(Data) EQ 8 then  data.T_Sky[*] = Value else data[*]= Value
end

;==========================================

function mrs_diff, d1, d2
m = d1
if  type_code(d1) EQ 8 then m.T_Sky = d1.T_Sky - d2.T_Sky  else m = d1 - d2
return,m 
end

;==========================================

function mrs_mult, d1, d2
m = d1
if  type_code(d1) EQ 8 then m.T_Sky = d1.T_Sky * d2.T_Sky  else m = d1 * d2
return,m 
end

;==========================================

function mrs_add, d1, d2
m = d1
if  type_code(d1) EQ 8 then m.T_Sky = d1.T_Sky + d2.T_Sky  else m = d1 + d2
return,m 
end

;==========================================

function mrs_getpix, Ima
if  type_code(Ima) EQ 8 then pix = Ima.t_sky else pix = Ima
return, pix
end

;==========================================

pro mrs_putpix, Ima, pix
if  type_code(Ima) EQ 8 then  Ima.t_sky = pix else Ima = pix
end

;==========================================

function mrs_mad, Ima
 if  type_code(Ima) EQ 8 then  Mad = median( abs(Ima.t_sky)) / 0.6745 $
 else Mad = median( abs(Ima)) / 0.6745
return, Mad
end

;==========================================

function mrs_absthreshold, Ima, T, soft=soft, l2=l2
Tima = Ima
Glesp=0
if type_code(Ima) EQ 8 then Glesp=1

if not keyword_set(l2) then begin
	if not keyword_set(soft) then begin
		if Glesp EQ 1 then begin
			index = where ( ABS(Ima.t_sky) LT T, count )
			if count GT 0 then TIma.t_sky[index] = 0
		end else begin
			index = where ( ABS(Ima) LT T, count )
			if count GT 0 then TIma[index] = 0
		end
	end else begin
		if type_code(Ima) EQ 8 then begin 
			d = Ima.t_sky
			softthreshold, d, T  
			TIma.t_sky = d
		end else begin
			softthreshold, Tima, T
		end
	end
end else begin
	Tima = Ima / (1. + T)
end

return, Tima
end

;==========================================

pro mrs_pos, Ima, T=T
if not defined(T) then T=0.  

if type_code(Ima) EQ 8 then begin
      index = where (Ima.t_sky LT T, count) 
      if count GT 0 then Ima.t_sky[index] = 0
end else begin
      index = where (Ima LT T, count) 
      if count GT 0 then Ima[index] = 0
end  
end

;==========================================

pro pixf2pix, nside, x, y, face, ipix
    nside = float(nside)
    ipix = face*nside*nside
    itmp =   x + y *nside
    tab2nest,2*nside,itmp,itmp2
    ipix = ipix + itmp2
end

pro pix2pixf, nside, ipix, x, y, face
    x = 0
    y = 0
    nside2 =long(nside)
    ipix2 = long(ipix)
    face = floor(ipix/(nside2*nside2))
    ;print,'face ',face
    index = face*nside2*nside2
    ipix = ipix - index ; ou un modulo...
    while nside2 ge 1 do begin
         quot  = ipix / (nside2*nside2)
	 ;print, quot
	 if quot eq 3 then begin 
	                  x = x +nside2
			  y = y +nside2
	 endif
	 if quot eq 2 then  y = y +nside2
	 if quot eq 1 then  x = x +nside2
	 ipix = ipix mod (nside2*nside2)
	 nside2 = nside2 /2
    
    endwhile
	  
 ; print,x,y,face  
 
end


pro pixf2ang, nside, x, y, face, theta, phi
; theta,phi in Radian
    pixf2pix,nside,long(x), long(y), long(face),ipix
    pix2ang_nest, nside, ipix, theta, phi 
end

pro ang2pix, nside, theta, phi, x, y, face
; theta,phi in Radian
     ang2pix_nest, nside, theta, phi, ipix
     pix2pixf,nside,ipix,x,y,face
end

pro ang2lb, theta, phi, l, b
; Theta(z-axis) 
; Phi (xy plane)  (Theta,Phi) radian -> (l,b) degrees
    l = phi * 180. / !PI
    b = 90. - theta*180./!PI
end

pro lb2ang, l, b, theta, phi
 ; (l,b) degrees ->  (Theta,Phi) radian

    phi = l * !PI / 180.
    theta = (90. - b)*!PI / 180.
end

pro pix2lb, nside, ipix, l, b
; healpix pixels --> (l,b) degrees
pix2ang_nest, nside, ipix, theta, phi 
ang2lb, theta,phi,l ,b
end

pro lb2pix, nside, l, b, ipix
; (l,b) degrees --> healpix pixels
  lb2ang, l ,b, theta,phi
  ang2pix_nest, nside, theta, phi, ipix
end


pro ang2radec, theta, phi, ra, dec, year=year, degree=degree
   if not keyword_set(year) then year=2000
   ang2lb, theta,phi,l ,b
   ; print, 'lb = ', l, b
   glactc, ra, dec, year, l, b, 2, degree=degree
end
pro radec2ang, ra, dec, theta, phi, year=year, degree=degree
   if not keyword_set(year) then year=2000
   glactc, ra, dec, year, l, b, 1, degree=degree
   ; print, 'lb = ', l, b
   lb2ang, l ,b, theta,phi
end

pro testradec
t1 = 10. / 180. * !PI
t2 = 40. / 180. * !PI
print, 'Theta-Phi = ', t1, t2
ang2radec, t1, t2, ra, dec, /degree
print, 'RADEC = ', ra, dec
radec2ang, ra, dec, theta, phi, /degree
print, 'Theta-Phi = ', theta / !PI * 180. , phi / !PI * 180.
end
pro testtethaphi, x, y, f
nside=32
print, x,y,f
pixf2ang,nside,x,y,f,theta,phi
print, theta,phi
ang2pix_ring, nside, Theta, Phi, poslist
RING2NEST, Nside, poslist, in
pix2pixf,nside,in,x,y,f
print, x,y,f
end

pro testangle_jl, a, f, f1
nside=64L
f = fltarr(nside,nside, 12)
f1 = f
theta = float(a) / 180. *  !PI
for i=0,359 do begin
   phi = float(i) / 180. *  !PI
   ; print,'phi ', phi
   ang2pix,nside,theta,phi,x,y,face
   F(x,y,face) = 1
   
   pixf2ang,nside,x,y,face,theta1,phi1
   ang2pix,nside,theta1,phi1,x1,y1,face1
   F1(x1,y1,face1) = 1
end
map = f2h(f)
map1 = f2h(f1)

tvs, map
tvs, map1

end

pro testangle, a, f
nside=64L
f = fltarr(nside,nside, 12)

theta = float(a) / 180. *  !PI
for i=0,359 do begin
   phi = float(i) / 180. *  !PI
   ;print,'phi ', phi
   ang2pix,nside,theta,phi,x,y,face
   F(x,y,face) = 1
end
map = f2h(f)
tvs, map

f2 = f(*,*,1)   
ind = where(f2 eq 1,count) & print,count
      y = ind / nside
      x = ind mod nside

       pixf2ang,64L,x,y,1,theta,phi
      print,'theta ',theta
      print,'phi ',phi
     
      
      help,ind
      
;ind = where f eq 1
;pix2ang,nside

f = fltarr(nside,nside, 12)

phi = float(a) / 180. *  !PI
for i=0,179 do begin
   theta = float(i) / 180. *  !PI
   ;print,'phi ', phi
   ang2pix,nside,theta,phi,x,y,face
   F(x,y,face) = 1
end
map = f2h(f)
tvs, map
for face=0,11 do begin
f2 = f(*,*,face)   
print,'face ',face
ind = where(f2 eq 1,count) & print,count
      if not (count eq 0) then  begin
      y = ind / nside
      x = ind mod nside

       pixf2ang,64L,x,y,1,theta,phi
      print,'theta ',theta
      print,'phi ',phi
     endif
      
      help,ind
endfor
end

pro test2
nside = 64.
npix = nside2npix(nside)
map  = fltarr(npix)
nbtest = 1
test1=0
if test1 eq 1 then begin
 for i=1,nbtest do begin
   map = map*0
   theta = randomu(seed)*!dpi
   phi  = randomu(seed)*2.*!dpi
   print,'theta :',theta,' phi : ', phi
  
   ang2pix_nest, nside,theta,phi,ipix
   map(ipix) =1
   get_all_faces,map,faces
   ind = where(faces eq 1)
   face  = ind/(nside*nside)
   print,'face :',face
   y = (ind mod (nside*nside))/nside
   x = ind mod nside
   print,'x y :',x,y
   ang2pix,nside,theta,phi,x1,y1,face1
   print,'x y estimé',x1,y1,face1
 endfor
end

 for i=1,nbtest do begin
   map = map*0
   x = long(randomu(seed)*nside)
   y = long(randomu(seed)*nside)
   face  = 0 ; long(randomu(seed)*12)
   get_all_faces,map,faces
   help,faces
   faces(x,y,face) = 1
   put_all_faces,faces,map
   ind = where (map eq 1)
   pix2ang_nest,nside,ind,theta,phi
   print,'theta reel ',theta
   print,'phi reel :', phi
   
   pixf2ang,nside,x,y,face,theta,phi
   print,'theta calculé ',theta
   print,'phi calculé :', phi
   
   ang2pix_ring, nside, theta, phi, ipring_r
   ; pix2vec_ring, nside,[ipring_r],vector

   ang2pix_ring, nside, Theta, Phi, poslist
   pix2vec_ring, nside, poslist, vector

 endfor


end


;======================================================================
  
function mrs_get_block_rms, Map, Cut=Cut, Norm=Norm, TabSig=TabSig, Mad=Mad, HighFreq=HighFreq, bin_size=bin_size
 
if not keyword_set(cut) then cut=2500
if not keyword_set(bin_size) then bin_size =16

if keyword_set(HighFreq) then begin
  Map1 = TOL(Map, Cut) 
  D = Map - Map1
end

if not keyword_set(Norm) then Norm=0.935  

nside = gnside(Map)
nb_bin = double(nside)/bin_size
TotalNbin = float(nb_bin) * nb_bin * 12.

if keyword_set(HighFreq) then ScaleData =  H2F(D) else  ScaleData =  H2F(Map)
ScaleRMRS = ScaleData
ScaleRMRS[*] = 0
TabSig = fltarr(TotalNbin)
			
var_scale = 0.
Ind = 0L
for f=0, 11 do begin
  for i=0,nb_bin-1 do begin
  for j=0,nb_bin-1 do begin    
        if not keyword_set(Mad) then sig_im =  sigma(ScaleData[i*bin_size:(i+1)*bin_size-1, j*bin_size:(j+1)*bin_size-1, f])   $
        else sig_im =  MAD(ScaleData[i*bin_size:(i+1)*bin_size-1, j*bin_size:(j+1)*bin_size-1, f])
        TabSig[Ind] = sig_im
        ScaleRMRS[i*bin_size:(i+1)*bin_size-1, j*bin_size:(j+1)*bin_size-1, f] = TabSig[Ind]
        Ind = Ind + 1
 endfor
 endfor
 endfor 
 
if keyword_set(Mad) then  Rms  = F2H(ScaleRMRS) / 0.935
return, ScaleRMRS
end

;======================================================================

function mrs_get_local_rms, Imag, halfsizeArcMin=halfsizeArcMin, WindowSize=WindowSize, mad=mad
rmap=-1
if N_PARAMS() NE 1  then begin 
        print, 'CALLING SEQUENCE: RmsMap = mrs_get_local_rms(imag, halfsizeArcMin=halfsizeArcMin, WindowSize=WindowSize)'
        goto, DONE
        end

nested=1
npix  = n_elements(imag)
nside = npix2nside(npix)
imag2 = imag

imag2 = reorder(imag,out='ring',in='nest')
rmap = imag2

if not keyword_set( WindowSize) then WindowSize = 5.

if not keyword_set( halfsize) then halfsize=pixel_size(nside) * WindowSize / 2.

; Apply erosion operator
for i = 0ul, npix-1 do begin
     ; Get the xyz coordinate
     pix2vec_ring, nside, i, vec

     ; Get the disk of radius halfsize
     query_disc, nside, vec, float(halfsize)/60.0, lpix, /deg

     if not keyword_set(mad) then rmap[i]  = sigma( imag2[lpix] ) $
     else rmap[i]  = MAD( imag2[lpix] ) 
end     

rmap = reorder(rmap,in='ring',out='nest')

DONE:

return,rmap
end


;======================================================================
 

function  mr1d_get_local_rms, Imag,  WindowSize=WindowSize, mad=mad
rmap=-1
if N_PARAMS() NE 1  then begin 
        print, 'CALLING SEQUENCE: RmsMap = mr1d_get_local_rms(Signal,  WindowSize=WindowSize)'
        goto, DONE
        end

 npix  = n_elements(imag)
 rmap = imag

if not keyword_set( WindowSize) then WindowSize = 5.
W2 = long(WindowSize/2)
 for i = 0L, npix-1 do begin
     minp = max( [0, i-W2])
     maxp = min( [Npix-1, i+W2])
     if not keyword_set(mad) then rmap[i]  = sigma( imag[minp:maxp] ) $
     else rmap[i]  = MAD( imag2[minp:maxp] ) 
end     
 
DONE:

return,rmap
end


;======================================================================

function mrs_bandpass_sigma, Cl, Filter,  lmax=lmax, npix=npix, NormVal= NormVal

if not keyword_set(NormVal) then  NormVal = 1.
vs = size(filter)
lm = vs[1]
vs = size(Cl)
lm1 = vs[1]
if not keyword_set(lmax) then lmax = min([lm, lm1]) -1
 
 CoefN = NormVal* NormVal
 IndL = lindgen(LMAX+1)
IndL = 2. * IndL + 1.
Sig = total( IndL * Cl[0:Lmax]  * Filter[0:Lmax]^2 / CoefN)
Sig = Sig / (4. * !DPI)
return, sqrt(Sig)
end

;======================================================================
  
