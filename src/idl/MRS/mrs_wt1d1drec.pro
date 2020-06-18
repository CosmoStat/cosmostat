; pro mrs_wt1d1drec,  WT1D1D,  RecData 
;+
; NAME:
;       MRS_WT1D1DREC
;
; PURPOSE:
;       Apply a 1D-1D wavelet recpnstruction. GLESP package must be installed.
;
; CALLING SEQUENCE:
;        mrs_wt1d1drec,  WT1D1D,  RecData
;
; INPUT:
;       WT1D1D = input IDL structure obtained from the routine mrs_wt1d1dtrans.pro
;
; OUTPUT:
;     RecData -- IDL 1D array: Healpix reconstructed map (nested)  
;
; EXAMPLES:
;    Compute the wavelet transform 1D-1D, and its reconstruction
;         ISAP> mrs_wt1d1dtrans, map_in, outWT
;         ISAP> mrs_wt1d1drec,  outWT, RecIma
;
; EXTERNAL CALLS:
;         Call glesp routine glesp2healpix.
;
; HISTORY:
;    File create date: April 2013,  Jean-Luc Starck
;-
;======================================================

pro inv_spherical_IUWT1D, W,  Xrec,  in_N=in_N, wt1d=wt1d, gen2=gen2

vs = size(W)
out_N = vs[1]
N = out_N
if  keyword_set(in_N) then N= in_N

NbrScale = vs[2]
MaxScale = long(alog(N)/alog(2))-1
Xrec = dblarr(N)

D = NbrScale-1
if D GT MaxScale-1 then D = MaxScale-1

; Recopy the wavelet transform in the buffer with the correct size
if D LE 0 then Xrec = congrid (W[*, NbrScale-1], N) $
else begin

   ;filters definition
   if keyword_set(Gen2) then  begin
      h0=double([1, 4, 6, 4, 1])       ;low-band filter
      h0=h0 / double(16)
   end

   WT1D = dblarr(N, D+1)
   WT1D[*,D] = congrid (W[*, NbrScale-1], N)
   for j=0,D-1 do WT1D[*,j] = congrid (W[*, j], N)

   ;filters definition
  
  Xrec = WT1D[*, D] 
   for j=D-1,0,-1 do begin
      if keyword_set(Gen2) then  begin
          h1  = dblarr(2^(j+2)+1)
          for k=0,4  do h1[k*2^j] = h0[k]
          Im_Out =convol(Xrec, h1, /CENTER, /EDGE_WRAP)
          Xrec = Im_Out + WT1D[*,j]
      end else  Xrec = Xrec + WT1D[*,j]
    ; print, j+1, ': ',  sigma(Xrec - total( WT1D[*, j:*], 2))
endfor
endelse
end

;======================================================

pro test1d
n = randomn(seed, 1000)*100.
spherical_IUWT1D, n, 5, w, /gen2
print, " rec "
inv_spherical_IUWT1D, W,  Xrec,  /gen2
plot, n-xrec
info, n-xrec
end


;======================================================

pro old_spherical_lat_IUWT1D, Data, NbrScale=NbrScale, WTLat
x_sky = Data.x_sky
y_sky = Data.y_sky
nb_lon  = max(y_sky)+1
nb_lat = Data.nx
out_N = max(y_sky)+10

WTLat = dblarr(out_N, nb_lat, NbrScale)
mil = max(y_sky)/2
count = 0L
for lat=0L,nb_lat-1 do begin
   deb = mil-(y_sky[lat]-1)/2
   X = (Data.t_sky)[count:count+y_sky[lat]-1]
   spherical_IUWT1D, X, NbrScale, WtX, out_N=out_N
   for j=0,NbrScale-1 do  WTLat[*,lat,j] = Wtx[*,j]
   count = count + y_sky[lat]
endfor
end


;======================================================

pro test2d, x=x, y=y, w=w, sx=sx, sy=sy, NBRSCALEX= NBRSCALEX, NBRSCALEY= NBRSCALEY, reset=reset
if keyword_set(reset) then begin
NBRSCALEX=0
NBRSCALEY=0
x=0
y=0
sx=0
sy=0
y=0
w=0
end

if not keyword_set(w) then mrs_wt1d1d, Data, w,  /ecl, /dirac, /gen, NBRSCALEX=NBRSCALEX, NBRSCALEY=NBRSCALEY
help, w.coef
vs =size(w.coef)
Npx = vs[1]
Npy = vs[2]
if not keyword_set(sx) then sx = NBRSCALEX-1
if not keyword_set(sy) then sy = NBRSCALEY-1
if not keyword_set(x) then x = Npx/2
if not keyword_set(y) then y = Npy/2
w.coef [*]=0
w.coef[x,y, sx,sy] = 100.
print, x,y,sx,sy
mrs_invwt1d1d,  W,  R, wtlat=wtlat, g=g
tvs, r
end

;======================================================

pro mrs_wt1d1drec,  WT1D1D,  RecData, WTLat=WTLat, g=g

NbrScaleX = WT1D1D.NbrScaleX
NbrScaleY = WT1D1D.NbrScaleY
Np = WT1D1D.Np
Nx = WT1D1D.Nx
Nside = WT1D1D.Nside
nb_lon = WT1D1D.nb_lon
nb_lat= WT1D1D.nb_lat
y_sky = WT1D1D.y_sky
x_sky = WT1D1D.x_sky
out_N = WT1D1D.out_N
WTLat = dblarr(out_N, nb_lat, NbrScaleX)
WT = WT1D1D.coef
gen2=WT1D1D.gen2
ecliptic=WT1D1D. ecliptic

; inverse shift
for jx=0L,NbrScaleX-1 do  begin
for jy=0L,NbrScaleY-1 do begin
   for lat=0L,nb_lat-1 do  WT[*,lat,jx, jy] = reverse( WT[*,lat,jx, jy])
   WT[*,*,jx, jy] = shift( WT[*,*,jx, jy], -out_N/2,0)
end
end

Wty = dblarr(nb_lat, NbrScaleY)
in_N= nb_lat
for lon=0L,nb_lon-1 do begin
   for jx=0L,NbrScaleX-1 do begin
   for jy=0,NbrScaleY-1 do  Wty[*,jy] = WT[lon,*,jx, jy]
   inv_spherical_IUWT1D, WtY,  X,  in_N=in_N, gen2=gen2
   WTLat[lon, *, jx] =  X
endfor
endfor

mil = max(y_sky)/2
count = 0L
Wtx = dblarr(out_N, NbrScaleX)
g = {t_sky: dblarr(WT1D1D.NpixGlesp), x_sky: x_sky, y_sky: y_sky, Nx:Nx, Np:Np}
for lat=0L,nb_lat-1 do begin
    in_N= y_sky[lat]
  ;  print, lat, ' ', in_N
  ;  help, Wtx
    X = dblarr(in_N)
    for jx=0,NbrScaleX-1 do Wtx[*,jx] =  WTLat[*,lat,jx]
    inv_spherical_IUWT1D, WtX,  X,  in_N=in_N, gen2=gen2
    deb = mil-(y_sky[lat]-1)/2
    endcount = count+y_sky[lat]-1
   G.t_sky[count: endcount] = X
   count = count + y_sky[lat]
endfor

RecData = g2h(g, /direct, nside=Nside)
if keyword_set(ecliptic) then RecData = coord_conversion(RecData, in_coordinate='E', out_coordinate='G') 
end

