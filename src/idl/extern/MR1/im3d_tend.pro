;-------------------------------------------------------------------------------
;+
; NAME:
;        IM3D_TEND
;
; PURPOSE:
;     For each spatial position (x,y) of the cube, we extract the 1D signal along   
;     the z-axis and remove the baseline (tendancy).
;
; CALLING:
;      RES = IM3D_TEND(Data,WinSize=WinSize)
;
; INPUTS:
;     Data : 3D fltarr = input data    
;
; OUTPUTS:
;     RES: 3D fltarr = output data. 
;
; KEYWORDS:
;      WinSize: integer: window size for the polynomial fit.
;                        Default is 100.
;
; EXTERNAL CALLS:
;       im1d_tend (IDL program)
;
; EXAMPLE:
;       Remove the baseline using a window size of 60 pixels
;              RES = IM3D_TEND(Data, WinSize=60)
;
; HISTORY:
;	Written: Jean-Luc Starck 2004.
;	July, 2004 File creation
;-

function im3d_tend, d, opt=opt
dout = d
vs = size(d)
nx = vs[1]
ny = vs[2]
nz = vs[3]

for i=0, nx-1 do begin
for j=0, ny-1 do begin
  v = reform(d[i,j,*])
  im1d_tend, v, s1, tend=tend, WinSize=WinSize
  dout[i,j,*] = s1
end
end
return, dout
end

 
