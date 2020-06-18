;+
; NAME:
;        mrs_glesp
;
; PURPOSE:
;	The MultiResolution on the Sphere (MRS) IDL code can work only if the file mrs_init
;       has been compile (.r mrs_glesp). All the MRS routine work with nested online HEALPIX
;       maps.
;
;       This file contains several routines
;
;    
;    function healpix2glesp, HMap,alm=alm,optNx= optNx,optNp = optnp
;     function h2g, h
;    function glesp2healpix, GMap,direct_map = direct_map,alm=alm
; 
;    read_glesp, file, data
;    
;    tvs, Data, graticule=graticule, png=png, TITLEPLOT=TITLEPLOT, COLT=COLT, NOBAR=NOBAR
;    view_glesp,g_in
;    view_glesp_rot,g_in,theta = theta , phi = phi
;    glesp_interpol,t_sky,out
;    glesp_desinterpol,g_in,out
;
;    write_glesp, filename, data
;    write_glesp2, filename, imag
;    write_glesp3, filename, imag
;    function diff( g1, g2)
;
;; HISTORY:
;	Written: Jean-Luc Starck and Pierrick Abrial, 2005
;	February, 2005 File creation
;-
;===============================================================


function healpix2glesp, HMap, alm=alm, optNx=optNx, optNp=optnp, lmax=lmax
alm=1  ; cmap routine is not anymore supported by GLESP. We need to use the Alm decomposition
if keyword_set(alm) then begin
   mrs_almtrans,Hmap, alm2, lmax=lmax
   mrs_almrec, alm2, Gmap, /to_glesp, nx= optNx, np = optnp
end else begin
    HMap2 = HMap ; copy a the original image preventing changment

   npix = (size(HMap2))[1]
   nside = npix2nside(npix)
   FileName = gettmpfilename()
   GlespFitsFile = gettmpfilename()
   GlespFitsErrFile = gettmpfilename()

;OptNx = " -nx " +  STRCOMPRESS(string(nside-1), /REMOVE_ALL) 
;OptNp = " -np " +   STRCOMPRESS(string(nside*2), /REMOVE_ALL) 
   command = "cmap -sp " + FileName + " -o " + GlespFitsFile 
   if keyword_set(OptNx) then  begin 
                         if optnx ge 4*nside then begin
                                print,'ud_grade'
                                nside = 2.^(1+floor(alog(optnx/4)/alog(2)))
						       print,nside
			                   ud_grade,HMap2,HMap2,nside_out=nside,order_in='nested',order_out='nested'
                         endif
                         command = command + " -nx "+string(optnx)
  endif
  if keyword_set(OptNp) then  command = command + " -np "+string(optnp)

;command = "cmap " + FileName + " -o " + GlespFitsFile + OptNx + OptNp

   write_fits_map, FileName, HMap2, /nested
; print, command
  OutFileStdOut=gettmpfilename()

; print, command
   spawn, command + ' >& ' + OutFileStdOut  

;print, command
;print, GlespFitsFile
   read_glesp, GlespFitsFile, GMap


   delete, FileName
  delete, GlespFitsFile
  delete, OutFileStdOut
endelse
 
return, GMap
end

;===============================================================

function h2g, HMap, alm=alm, Nx=Nx, Np=np, lmax=lmax
if keyword_set(alm) then alm=1
return,  healpix2glesp(HMap, alm=alm,optNx=Nx, optNp=np, lmax=lmax)
end

;===============================================================


function glesp2healpix, GMap,direct_map = direct_map, alm=alm, nside=nside
;print,'toto'
alm=1  ; cmap routine is not anymore supported by GLESP. We need to use the Alm decomposition

if keyword_set(nside) then alm =1 
if keyword_set(alm) then begin
    mrs_almtrans,Gmap,alm2
    if keyword_set(nside) then mrs_almrec,alm2,Hmap,to_healpix=nside $ 
    else mrs_almrec,alm2,Hmap,/to_healpix
   return,HMap
endif

;npix = (size(HMap))[1]
;nside = npix2nside(npix)
FileName = gettmpfilename()
GlespFitsFile = gettmpfilename()
GlespFitsErrFile = gettmpfilename()

write_glesp,GlespFitsFile,gmap

command = "f2map "+GlespFitsFile + " -o " +FileName
if keyword_set(direct_map) then command = command + " -rp"

OutFileStdOut=gettmpfilename()
spawn, command + ' >& ' + OutFileStdOut  

read_fits_map,FileName,hmap

delete, FileName
delete, GlespFitsFile
delete, GlespFitsErrFile
delete, OutFileStdOut

return, reorder(hmap,in='ring',out='nested')
end


function g2h, GMap,direct_map = direct_map, alm=alm, nside=nside
return, glesp2healpix(GMap,direct_map = direct_map, alm=alm, nside=nside)
end

;===============================================================

pro read_glesp, file, data
  ; print,file
  g_read_fits_map,file,t_sky,x_sky,y_sky
  nx = (size(x_sky))[1]
  np = max(y_sky)

data = {t_sky : t_sky , x_sky : x_sky , y_sky : y_sky , nx : nx , np : np }
end
 
;===============================================================
;pro tvs, Data, graticule=graticule, png=png, TITLEPLOT=TITLEPLOT, COLT=COLT, NOBAR=NOBAR
;
;  if type_code(data) eq 8 then  view_glesp_rot,data    else mollview, Data, /online, /nested, NOBAR=NOBAR, graticule=graticule, png=png, TITLEPLOT=TITLEPLOT, COLT=COLT 
;end

;===============================================================


function g2ima, g_in

x_sky = g_in.x_sky
y_sky = g_in.y_sky
t_sky2 = g_in.t_sky
nb_lat = (size(x_sky))[1]
imag = fltarr(max(y_sky)+10,nb_lat)
count = 0
mil = max(y_sky)/2
for lat=0,nb_lat-1 do begin
    	deb = mil-(y_sky[lat]-1)/2
    	imag[deb:deb+y_sky[lat]-1,lat] = t_sky2[count:count+y_sky[lat]-1]
    	count = count + y_sky[lat]
endfor
return, imag
end

;===============================================================

pro view_glesp,g_in

x_sky = g_in.x_sky
y_sky = g_in.y_sky
t_sky2 = g_in.t_sky
if N_params() ne 1 then begin print,'view_glesp,glesp_struc'
   goto,done
endif


nb_lat = (size(x_sky))[1]

imag = fltarr(max(y_sky)+10,nb_lat)

count = 0
mil = max(y_sky)/2
for lat=0,nb_lat-1 do begin
deb = mil-(y_sky[lat]-1)/2
imag[deb:deb+y_sky[lat]-1,lat] = t_sky2[count:count+y_sky[lat]-1]
count = count + y_sky[lat]
endfor
window,/free,xsize=800,ysize=400
tvscl,congrid(imag,800,400)
done:
end
;===============================================================
;===============================================================

pro view_glesp_rot,g_in,theta = theta , phi = phi

 if N_params() ne 1 then begin print,'view_glesp,glesp_struc '
 goto,fin
endif
x_sky = g_in.x_sky
y_sky = g_in.y_sky
t_sky = g_in.t_sky
write_glesp,'tmp_rot.fits',g_in
if keyword_set(theta) or keyword_set(phi) then begin
   command = "difmap tmp_rot.fits"
   if keyword_set(theta) then command = command + " -dt "+string(theta)
   if keyword_set(phi) then command = command + " -dp "+string(phi)
   command = command + ' -o tmp_rot.fits'
   spawn,command

endif 

spawn,'f2fig tmp_rot.fits -o tmp_rot.gif'

read_gif,'tmp_rot.gif',imag,r,g,b


delete, 'tmp_rot.gif'
delete, 'tmp_rot.fits'
tvlct,r,g,b

window,/free,xsize=800,ysize=440
tv,imag
fin:
end
;===============================================================



pro glesp_interpol,t_sky,out

x_sky = t_sky.x_sky
y_sky = t_sky.y_sky
t_sky2 = t_sky.t_sky
 
 if N_params() ne 2 then print,'glesp_interpol,glesp_struc,out'

nb_lat = (size(x_sky))[1]
nb_lon  = max(y_sky)+1
imag = fltarr(nb_lon,nb_lat)

count = 0
mil = max(y_sky)/2
for lat=0,nb_lat-1 do begin
deb = 0; mil-(y_sky[lat]-1)/2
imag[deb:deb+y_sky[lat]-1,lat] = t_sky2[count:count+y_sky[lat]-1]
inter = (findgen(nb_lon)*y_sky[lat])/nb_lon

;imag[0:nb_lon-1,lat] = congrid (imag[deb:deb+y_sky[lat]-1,lat],nb_lon,cubic=1)
imag[0:nb_lon-1,lat] = imag[inter,lat]
count = count + y_sky[lat]
endfor
window,/free,xsize=800,ysize=400
tvscl,congrid(imag,800,400)
out = {t_sky:imag,x_sky:x_sky,y_sky:lonarr(nb_lat)+nb_lon,nx:t_sky.nx,np:t_sky.np,old_ysky:y_sky}

 
end
;===============================================================
pro glesp_desinterpol,g_in,out

x_sky = g_in.x_sky
y_sky = g_in.y_sky
t_sky = g_in.t_sky

if N_params() ne 2 then print,'glesp_desinterpol,glesp_struc,out '

nb_lat = (size(x_sky))[1]
nb_lon  = max(y_sky)+1
imag = fltarr(total(g_in.old_ysky))

nb_lat = g_in.nx

count = 0
deb = long(0)
for lat=0,nb_lat-1 do begin

imag[deb:deb+g_in.old_ysky[lat]-1] = congrid(t_sky[*,lat],g_in.old_ysky[lat])
deb = deb + g_in.old_ysky[lat]
; print,deb

;inter = (findgen(nb_lon)*y_sky[lat])/nb_lon
;imag[0:nb_lon-1,lat] = congrid (imag[deb:deb+y_sky[lat]-1,lat],nb_lon,cubic=1)
;imag[0:nb_lon-1,lat] = imag[inter,lat]

endfor
;window,/free,xsize=800,ysize=400
;tvscl,congrid(imag,800,400)


out = {t_sky:imag,x_sky:x_sky,y_sky:g_in.old_ysky,nx:g_in.nx,np:g_in.np}
end


;===============================================================
; -----------------------------------------------------------------------------
;
;
; -----------------------------------------------------------------------------
pro write_glesp, filename, data
; write a glesp structure in fits glesp format
; data structure is :
;{t_sky : t_sky , x_sky : x_sky , y_sky : y_sky , nx : nx , np : np }
; t_sky : data map
; number of pix per row
; lat of each row
; nx number of row (optionnal)
; np number of pixel on equator


 t_sky = data.t_sky
 x_sky = data.x_sky
 y_sky = data.y_sky

defsysv, '!DEBUG', EXISTS = i  ; check if astrolib variables have been set-up
if (i ne 1) then astrolib       ; if not, run astrolib to do so
;
;t_sky = data.t_sky
;x_sky = data.x_sky
;y_sky = data.y_sky


if N_params() lt 2 or N_params() gt 5 then begin
    print, 'syxtax : write_gles,filename,data_struc'
    return
endif


; ------- primary unit ----------
; opens the file, write a minimal header and close it
WRITEFITS,filename,0

FXHMODIFY, FILENAME, 'BITPIX', 8
FXHMODIFY, FILENAME, 'COMMENT', "'GLESP MAP'"
; create the minimal extension header
nrows = 1

FXBHMAKE,xthdr,nrows


FXBADDCOL,index,xthdr,x_sky
;print,index
FXBADDCOL,index,xthdr,y_sky
;print,index
FXBADDCOL,index,xthdr,t_sky
;print,index
;print,xthdr
FXBCREATE, UNIT, FILENAME, xthdr
; add data
FXBWRITE, UNIT, x_sky, 1, 1
FXBWRITE, UNIT, y_sky, 2, 1
FXBWRITE, UNIT, t_sky, 3, 1
; modify header
FXHMODIFY, FILENAME, 'TTYPE1', 'COS(THETA)',EXTENSION=1
FXHMODIFY, FILENAME, 'TUNIT1', 'COS(RAD)',EXTENSION=1
FXHMODIFY, FILENAME, 'TTYPE2', 'DIM',EXTENSION=1
FXHMODIFY, FILENAME, 'TUNIT2', 'UNITLESS',EXTENSION=1
FXHMODIFY, FILENAME, 'TTYPE3', 'TEMPERATURE',EXTENSION=1
FXHMODIFY, FILENAME, 'TUNIT3', 'K',EXTENSION=1
FXHMODIFY, FILENAME, 'OBJECT', 'CMB map',EXTENSION=1
FXHMODIFY, FILENAME, 'ORIGIN', 'GLESP',EXTENSION=1
FXHMODIFY, FILENAME, 'PIXTYPE', 'GLESP',EXTENSION=1
; close file
FXBFINISH,unit
end 


;=======================================================================


pro write_glesp2, filename, imag
; write an image in glesp fits format
; data structure is :
;{t_sky : t_sky , x_sky : x_sky , y_sky : y_sky , nx : nx , np : np }
; t_sky : data map
; number of pix per row
; lat of each row
; nx number of row (optionnal)
; np number of pixel on equator

taille_x = (size(imag))[1]
taille_y = (size(imag))[2]

t_sky = float(reform(imag,taille_x*taille_y))

y_sky = lonarr(taille_y)+taille_x
x_sky = findgen(taille_y)/(taille_y/2) -1
;x_sky = findgen(taille_y)/(taille_y/4) -0.5

 ;t_sky = data.t_sky
 ;x_sky = data.x_sky
 ;y_sky = data.y_sky

defsysv, '!DEBUG', EXISTS = i  ; check if astrolib variables have been set-up
if (i ne 1) then astrolib       ; if not, run astrolib to do so
;
;t_sky = data.t_sky
;x_sky = data.x_sky
;y_sky = data.y_sky


if N_params() lt 2 or N_params() gt 5 then begin
    print, 'syxtax : write_gles,filename,data_struc'
    return
endif

; ------- primary unit ----------
; opens the file, write a minimal header and close it
WRITEFITS,filename,0

FXHMODIFY, FILENAME, 'BITPIX', 8
FXHMODIFY, FILENAME, 'COMMENT', "'GLESP MAP'"
; create the minimal extension header
nrows = 1

FXBHMAKE,xthdr,nrows


FXBADDCOL,index,xthdr,x_sky
;print,index
FXBADDCOL,index,xthdr,y_sky
;print,index
FXBADDCOL,index,xthdr,t_sky
;print,index
;print,xthdr
FXBCREATE, UNIT, FILENAME, xthdr
; add data
FXBWRITE, UNIT, x_sky, 1, 1
FXBWRITE, UNIT, y_sky, 2, 1
FXBWRITE, UNIT, t_sky, 3, 1
; modify header
FXHMODIFY, FILENAME, 'TTYPE1', 'COS(THETA)',EXTENSION=1
FXHMODIFY, FILENAME, 'TUNIT1', 'COS(RAD)',EXTENSION=1
FXHMODIFY, FILENAME, 'TTYPE2', 'DIM',EXTENSION=1
FXHMODIFY, FILENAME, 'TUNIT2', 'UNITLESS',EXTENSION=1
FXHMODIFY, FILENAME, 'TTYPE3', 'TEMPERATURE',EXTENSION=1
FXHMODIFY, FILENAME, 'TUNIT3', 'K',EXTENSION=1
FXHMODIFY, FILENAME, 'OBJECT', 'CMB map',EXTENSION=1
FXHMODIFY, FILENAME, 'ORIGIN', 'GLESP',EXTENSION=1
FXHMODIFY, FILENAME, 'PIXTYPE', 'GLESP',EXTENSION=1
; close file
FXBFINISH,unit
end 

;=======================================================================


pro write_glesp3, filename,data,nb_pix_per_ring
npix = (size(data))[1]
;print,npix
nside = npix2nside(npix)
;print,nside
nb_pix_per_ring = lonarr(4*nside)
nb_pix_per_ring[indgen(nside)] = 4 * indgen(nside)
nb_pix_per_ring[nside+indgen(2*nside)] = fltarr(2*nside)+ 4 * nside 
nb_pix_per_ring[3*nside+indgen(nside)] = 4 * (nside-indgen(nside))

x_sky = findgen(4*nside)/(2*nside) -1

data2 = { t_sky : data,x_sky : x_sky, y_sky : nb_pix_per_ring}
write_glesp, filename, data2
;plot ,nb_pix_per_ring
end


;===============================================================

;===============================================================
 
;===============================================================

;pro glesp_alm_trans, imag, TabAlm, complex=complex, TabNbrM=TabNbrM, NbrL=NbrL
pro glesp_alm_trans, data, TabAlm, complex=complex, TabNbrM=TabNbrM, NbrL=NbrL, index=index 
if not keyword_set(NbrL)  then NbrL = min([data.nx/2,data.np/4])
print, "MY ALM"
  filename_tmp = "xxxtmp.fits"
  write_glesp, filename_tmp, data
  command = "cl2map -map "+  filename_tmp
  if keyword_set(NbrL) then command = command + " -lmax "+ string(NbrL)
  command = command + " -ao alm.fits"
  ;$cl2map -map map.glesp.fits -lmax 255 -ao alm.fits -o cl2.txt
  print,command
  spawn,command
  fits2alm,index,tabAlm,'alm.fits'
 ; alm2tab,lm,tabalm,complex=complex;, TabNbrM=TabNbrM;, NbrL=NbrL
end

;===============================================================
function diff, Data1, data2
   toto = {t_sky : data1.t_sky-data2.t_sky,x_sky: data1.x_sky,y_sky : data1.y_sky,nx:data1.nx,np:data1.np} 
   return,toto
end

function glesp_convol, Data1, Data2
    glesp_alm_trans, Data1, TabAlmData1, /complex
    glesp_alm_trans, Data2, TabAlmData2, /complex
    TabAlmData2 = TabAlmData1 * TabAlmData2
    glesp_alm_itrans, TabAlmData2, DataConvol
    return, DataConvol
end

pro mod_phase,alm,module,phase

   module = sqrt(alm[*,0]^2+alm[*,1]^2)
   phase = atan (alm[*,1],alm[*,0])

end

;===============================================================


pro world2glesp, world, map_glesp
; project image carre en glesp

;help,world
; world must be [*,*,3]
;help,world

write_glesp2,'map_g3D.fits',world[*,*,0]
map_tmp = (mrs_read('map_g3D.fits'))
mrs_almtrans,map_tmp,alm
mrs_almrec,alm,map_0
map_glesp = REPLICATE(map_0,3)

for i = 1 ,2 do begin
write_glesp2,'map_g3D.fits',world[*,*,i]

map_tmp = (mrs_read('map_g3D.fits'))
mrs_almtrans,map_tmp,alm
mrs_almrec,alm,map_tmp
map_glesp(i) = map_tmp
endfor

delete,'map_g3D.fits'
end


pro make_mask,gl,mask

mask = gl(0)
taille = mask.nx* mask.np
pourcent = .5
random = floor(randomu(seed,taille*pourcent)*(taille-1))

mask.t_sky(*)= 0
mask.t_sky(random) = 1


end

;===============================================================

pro inpaint_world,world,mask=mask,world2,map_masked,wave,alm
world2 = world
wave = world
alm =world
map_masked = world
if not keyword_set(mask) then make_mask,world,mask

for i=0,2 do begin
map_masked[i].t_sky = world[i].t_sky * mask.t_sky

mrs_mca,map_masked[i],mask = mask,tabc,selecttrans=[4],niter=10,/expo,/cstsigma
world2[i].t_sky = tabc[0].t_sky;+tabc[1].t_sky

;wave[i].t_sky = tabc[0].t_sky
;alm[i].t_sky = tabc[1].t_sky
endfor

end


