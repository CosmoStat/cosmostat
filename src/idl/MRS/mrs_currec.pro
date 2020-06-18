;+
; NAME: 
;       MRS_CURREC
;
; PURPOSE: 
;        Reconstruct an image on the Sphere from its curvelet transform (see mrs_curtrans).   
;
; CALLING:
;       MRS_CURREC, Cur_Struct, result
;
; INPUT:
;       Cur_Struct -- IDL structure: Curvelet transform structure (see MRS_CURTRANS) 
;          
; OUTPUTS:
;      Result --  1D IDL array: Healpix image (nested format)
;             
; EXTERNAL CALLS
;           mrs_ridrec (IDL prog)
;           mrs_wtput (IDL prog)
;           mrs_curput (IDL prog)
;           mrs_curget (IDL prog)
;
; EXAMPLE:
;   mrs_curtrans, Imag, Cur 
;   mrs_currec, Cur, RecIma
;       Curvelet transform and reconstruction
;
; HISTORY:
;       Written: Jean-Luc Starck 2005.
;       February, 2005 File creation
;-
;-----------------------------------------------------------------

pro mrs_currec, Cur, HealpixIma 

if N_PARAMS() LT 2 then begin 
         print, 'CALL SEQUENCE: mrs_currec, Cur_Struct, result'
        goto, DONE
        end

NS = Cur.NbrScale
for j=0,NS-2 do begin
     my_command = 'mrs_ridrec, Cur.RidScale'+strcompress(string(j+1), /remove_all) + ', Scale'
     my_command = strcompress( my_command, /remove_all)
     ; print, 'cmd = ',  my_command
     ACK = EXECUTE( my_command) 
     ; info, scale
     ; mrs_wtput, Cur.WT, Scale, j
     if Cur.WT.pyrtrans EQ 0 then Cur.WT.coef[*,j] = Scale $
     else begin
        my_command = 'Cur.WT.scale'+strcompress(string(j+1), /remove_all)
        my_command =  my_command + '= Scale '  
        my_command = strcompress( my_command, /remove_all)
        ; print, 'cmd = ',  my_command
        ACK = EXECUTE( my_command) 
    end 
end

j=NS-1
if Cur.WT.pyrtrans EQ 0 then Cur.WT.coef[*,j] = Cur.LastScale $
else begin
        my_command = 'Cur.WT.scale'+strcompress(string(j+1), /remove_all)
        my_command =  my_command + '= Cur.LastScale '  
        my_command = strcompress( my_command, /remove_all)
        ; print, 'cmd = ',  my_command
        ACK = EXECUTE( my_command) 
end 
; mrs_wtput, Cur.WT, Cur.LastScale, NS-1

if Cur.pyrtrans EQ 0 then mrs_wtrec, Cur.WT,  HealpixIma $
else mrs_pwtrec, Cur.WT, HealpixIma 

DONE:
end



pro fig_curbackproj, Cur, c
; Example of backproject of Dirac
i = fltarr(196608L)
; mrs_curtrans, i, Cur, nbrscale=5, opt=' -O -t5 ', /SSR
; mrs_curtrans, i, Cur, nbrscale=5, opt=' -O', /undec

; mrs_curtrans, i, Cur, nbrscale=5, opt=' -O', /SSR
mrs_curtrans, i, Cur, nbrscale=5, opt=' -O -t5 ', /SSR, /undec

vs = size(Cur.TabNorm)
NxTNorm = vs[1]
NyTNorm = vs[2]

Cur.RIDSCALE1.coef[*] = 0
Cur.RIDSCALE2.coef[*] = 0
Cur.RIDSCALE3.coef[*] = 0
Cur.RIDSCALE4.coef[*] = 0
Cur.LASTSCALE[*] = 0
NS=Cur.NBRSCALE
 for j=0,Cur.NBRSCALE-2 do begin
  if j eq 0 then RidTrans = Cur.RIDSCALE1 $
  else if j eq 1 then RidTrans = Cur.RIDSCALE2 $
  else if j eq 2 then RidTrans = Cur.RIDSCALE3 $
  else if j eq 3 then RidTrans = Cur.RIDSCALE4 $
  else  RidTrans = Cur.RIDSCALE5
  s2d = j
  s1d = j  
  s2d1 = s2d
  s1d1 = s1d
  b = Cur.TabBlockSize[j]
    s1d1 = s1d
    if s2d1 GE NyTNorm then s2d1 = NyTNorm - 1
    if s1d1 GE NxTNorm then s1d1 = NxTNorm - 1
    if s1d1 EQ Cur.TabNbrScaleRid[s2d]-1 and Cur.TabNbrScaleRid[s2d] lt NxTNorm then s1d1 = s1d1-1
    Norm = Cur.TabNorm[s1d1,s2d1] * sqrt(Cur.TabBlockSize[s2d])
    
    
  CubeBand = mrs_curget(Cur, j,0) 
  CubeBand1 = mrs_curget(Cur, j,1) 
  vs = size(CubeBand)
  Nx = vs[1] / RidTrans.nxb
  Ny = vs[2] / RidTrans.nyb
  CubeBand[Nx/2, b/4, j] = 100.  / Norm
  mrs_curput, Cur, CubeBand, j, 0
 
  vs = size(CubeBand1)
  Nx = vs[1] / RidTrans.nxb
  Ny = vs[2] / RidTrans.nyb
  CubeBand1[Nx/2, b/2, NS+j] = 100.  / Norm
  CubeBand1[Nx/2, 3*b/2, 11-j] = 100./ Norm
  mrs_curput, Cur, CubeBand1, j, 1
  end
 
mrs_currec,  Cur, c 
tvs, c
end
