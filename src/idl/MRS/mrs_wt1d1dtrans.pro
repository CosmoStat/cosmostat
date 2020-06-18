; pro mrs_wt1d1dtrans, Data, outWT,  ecliptic=ecliptic, NbrScaleX=NbrScaleX, NbrScaleY=NbrScaleY, Dirac=Dirac, healpix=healpix, HWT=HWT, gen2=gen2
;+
; NAME:
;       MRS_WT1D1DTRANS
;
; PURPOSE:
;       Apply a 1D-1D wavelet transform on a healpix image. An undecimated 1D wavelet transform is first applied along the longitude for each 
;       wavelet band, and then another undecimated 1D wavelet transform is applied along the latitude for each longitude.
;        The output is therefore a four dimension array (x,y,sx,sy), i.e. position in the image (x,y). scale sx on the longitude, and scale sy on the latitude.
;        If the keyword ecliptic is set, the coordinate are first transform from galactic to ecliptic before the transform is applied.
;        The GLESP package MUST be installed in order to have this routine working.
;
; CALLING SEQUENCE:
;        mrs_wt1d1dtrans, map_in, outWT,  ecliptic=ecliptic, NbrScaleX=NbrScaleX, NbrScaleY=NbrScaleY, healpix=healpix, HWT=HWT, gen2=gen2
;
; INPUT:
;       map_in = input HEALPix map (nested format)
;
; OUTPUT:
;      outWT: IDL structure  with the following fields:  
;                  coef:  4D IDL array -- [*,*, 0:NbrScaleX-1,  0: NbrScaleY-1] is the wavelet image at scale sx along longitude and scale sy along latitude.
;                  ecliptic: scalar: if set, an ecliptic to galactic coordinate conversion  is applied before transform.
;                  NbrScaleX: scalar -- number of scales along the longitude
;                  NbrScaleY: scalar -- number of scales along the latiude

; INPUT KEYWORDS:
;       gen2:  scalar -- if set, the second starlet generation algorithm (i.e. with positive filters for the reconstruction) is applied.
;       healpix: scalar -- if set, convert also each wavelet scale in healpix format. 
;
; OUTPUT KEYWORDS:
;       HWT -- IDL 3D array -- Conversion of the coef array into healpix maps:  [*,0:NbrScaleX-1,  0: NbrScaleY-1] 
;
; EXAMPLES:
;    Compute the wavelet transform 1D-1D, and visualisation of the wavelet band 2-3  (scale 2 along longitude and 3 along latitude) 
;         ISAP> mrs_wt1d1dtrans, map_in, outWT
;         ISAP> tvscl, outWT[*,*,2,3]  
;
;    Same transformation, but convert also all wavelet bands into healpix format.
;         ISAP> mrs_wt1d1dtrans, map_in, outWT, /healpix, HWT=HWT
;         ISAP> mrs_tv, HWT[*,2,3]  
;
; EXTERNAL CALLS:
;         Call glesp routine healpix2glesp.
;
; HISTORY:
;    File create date: April 2013,  Jean-Luc Starck
;-
;======================================================

pro spherical_IUWT1D, X, NbrScale, W, out_N=out_N, gen2=gen2
; X = input signal to transform
; W = output wavelet scale [0: out_N-1,  NbrScale], wavelet coeff are interpolated if out_N is larger than the size of X
;        if NbrScale is too large compare to the size of input data, then a smaller number NS of wavelet decomposition is done, 
;       but the output array still have NbrScale for the second dimension, and the largest wavelet scale is always stored in W[*, NbrScale-1]
; 
D = NbrScale-1
vs=size(X)
n=vs(1)
MaxScale = long(alog(n)/alog(2))-1
if D GT MaxScale-1 then D = MaxScale-1
if not keyword_set(out_N) then out_N = N
; print, 'MaxScale=', MaxScale, N, NbrScale, out_n

;filters definition
h1=double([1, 4, 6, 4, 1])       ;low-band filter
h1=h1/double(16)
n1=5L

cj = double(X)        ;high-band filter
W = dblarr(out_N, NbrScale)
 
for i=0,D-1 do begin
  ck=convol(cj,h1, /CENTER, /EDGE_WRAP)
  if keyword_set(Gen2) then begin
  		   Im_Aux =convol(ck, h1, /CENTER, /EDGE_WRAP)
 		   W(*,i)= congrid(cj-Im_Aux, out_N, /cubic)
   end else  W(*,i)= congrid(cj-ck, out_N, /cubic)
  cj=ck ;update the convolution input signal 
  ;update the kernel filter bank
  h2= fltarr(n1*2-1) 
  for j=0,n1-2 do begin
    h2[2*j]=h1[j];
    h2[2*j+1]=0;
 endfor
  h2[2*n1-2]=h1[n1-1]
  h1=h2
  n1=n1*2-1  
endfor
W(*,NbrScale-1)= congrid(Cj, out_N, /cubic) 
end

;======================================================

pro mrs_wt1d1dtrans, Data, outWT,  ecliptic=ecliptic, NbrScaleX=NbrScaleX, NbrScaleY=NbrScaleY, Dirac=Dirac, healpix=healpix, HWT=HWT, gen2=gen2

Np = N_ELEMENTS(Data)
Nside = gnside(Data)
if not keyword_set(ecliptic) then ecliptic=0
if not keyword_set(gen2) then gen2 =0

if not keyword_set(Dirac) then begin
if keyword_set(ecliptic) then Map = coord_conversion(Data, in_coordinate='G', out_coordinate='E') $
else Map = Data
end

if keyword_set(Dirac) then begin
   Map = getdirac(nside=64)
   Map = mrs_resize(Map, nside=128)
   Data = Map
   Np = N_ELEMENTS(Data)
   Nside = gnside(Data)
end

g = h2g(Map)
x_sky = G.x_sky
y_sky = G.y_sky
nb_lon  = max(y_sky)+1
nb_lat = G.nx
out_N = max(y_sky) + 1 ; +10
if not keyword_set(NbrScaleX) then NbrScaleX = long(alog(out_N)/alog(2))-2
if not keyword_set(NbrScaleY) then NbrScaleY = long(alog(nb_lat)/alog(2))-2
; print, nb_lat, NbrScaleY

WTLat = dblarr(out_N, nb_lat, NbrScaleX)
WT = dblarr(out_N, nb_lat, NbrScaleX, NbrScaleY)

mil = max(y_sky)/2
count = 0L
for lat=0L,nb_lat-1 do begin

   deb = mil-(y_sky[lat]-1)/2
   X = (G.t_sky)[count:count+y_sky[lat]-1]
   ; vs = size(x)
   ; print, lat, ' ', vs[1]

   spherical_IUWT1D, X, NbrScaleX, WtX, out_N=out_N, gen2=gen2
   for jx=0,NbrScaleX-1 do  WTLat[*,lat,jx] = Wtx[*,jx]
   count = count + y_sky[lat]
   ;  if lat EQ 0 then    help, WtX
   ;  if lat EQ 0 then    help, X, out_N

endfor

for lon=0L,nb_lon-1 do begin
for jx=0L,NbrScaleX-1 do begin
   X = reform(WTLat[lon, *, jx])
   spherical_IUWT1D, X, NbrScaleY, WtY, gen2=gen2
   for jy=0,NbrScaleY-1 do  WT[lon,*,jx, jy] = Wty[*,jy]
endfor
endfor

; Concert Glesp scales to Healpix scale in Galactic coordinates
HWT = 0
if keyword_set(healpix) then begin
HWT = dblarr(Np, NbrScaleX, NbrScaleY)
mil = max(y_sky)/2

for jx=0L,NbrScaleX-1 do  begin
for jy=0L,NbrScaleY-1 do begin
   Scale = WT[*,*,jx, jy]
   ; print, 'Gjx = ', jx, ', jy = ', ', Min = ', Min(G.t_sky),  ', Max = ', Max(G.t_sky),  ', Sigma = ', sigma(G.t_sky)

   ; if keyword_set(Dirac) then 
   count = 0L
   G.T_Sky[*]=0
   for lat=0L,nb_lat-1 do begin
       deb = mil-(y_sky[lat]-1)/2
       G.T_Sky[count:count+y_sky[lat]-1] = congrid( reform(Scale[*,lat]), y_sky[lat])
       count = count + y_sky[lat]
   endfor
   ; view_glesp,g
   HealpixScale = G2H(G, nside=nside, /dir)
   ; print, 'jx = ', jx, ', jy = ', ', Min = ', Min(G.t_sky),  ', Max = ', Max(G.t_sky),  ', Sigma = ', sigma(G.t_sky)
   ; tvscl, g2ima(G)
   
  ; if keyword_set(ecliptic) then HealpixScale = coord_conversion(HealpixScale, in_coordinate='E', out_coordinate='G')
   HWT[*, jx, jy] = HealpixScale
endfor
endfor
endif

; Shift the map so they look similar to healpix mollview maps
for jx=0L,NbrScaleX-1 do  begin
for jy=0L,NbrScaleY-1 do begin
   WT[*,*,jx, jy] = shift( WT[*,*,jx, jy], out_N/2,0)
  for lat=0L,nb_lat-1 do  WT[*,lat,jx, jy] = reverse( WT[*,lat,jx, jy])
end
end

vs = size( g.t_sky)
NpixGlesp = vs[1]
outWT = {coef:wt, HWT:HWT, ecliptic:ecliptic, NbrScaleX:NbrScaleX, NbrScaleY:NbrScaleY, Nside:Nside, x_sky: x_sky, y_sky: y_sky, nb_lon: nb_lon, nb_lat: nb_lat, out_N: out_N, Nx: G.Nx, Np: G.Np, NpixGlesp:NpixGlesp, gen2:gen2}
end


;======================================================


pro mrs_tvwt1d, w
ind=1
vs = size(w.coef)
Nx= vs[1]
Ny = vs[2]
Sx = vs[3]
Sy = vs[4]
for i=3,sx-1 do begin
for j=3,sy-1 do begin
  window, ind, xsize=nx, ysize=ny, title='Scale ' + strc(i+1)+'-'+strc(j+1)
  tvscl, w.coef[*,*,i,j]
  ind = ind + 1
end
end

end

