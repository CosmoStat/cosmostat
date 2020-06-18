;+
; NAME: 
;       MRS_RIDGET
;
; PURPOSE: 
;        Extract a ridgelet band from the ridgelet transform (see mrs_ridtrans).
;        If the keyword NormMad is set, a normalization is applied to all coefficients.   
;
; CALLING:
;       result = MRS_RIDGET( Rid_Struct, ScaleRid, NormMad=NormMad, ImaMean=ImaMean, ImaMad=ImaMad )
;
; INPUT:
;       Rid_Struct -- IDL structure: Ridgelet transform structure (see MRS_RIDTRANS) 
;       ScaleRid  -- int:  Ridgelet band number (must be between 0 and Rid_Struct.NBRSCALE-1)
;       NormMad -- scalar: if set, normalize the coefficients by the Median Absolution Deviation
;                          of all coefficients at a given position in the block.
;       ImaMad  -- 2D fltarr: Image containing the normalization parameter
;       ImaMean -- 2D fltarr: Image containing  the mean value for all coefficients at a given position.
;           
;
; OUTPUTS:
;      Result -- 4D IDL float array[*,*,*,12]: Ridgelet band
;             
; EXTERNAL CALLS
;           mrs_ridget
;
; EXAMPLE:
;   mrs_ridtrans, Imag, Rid 
;   Band = mrs_ridget(Rid, 1)
;       Extract the  second ridgelet scale  
;
; HISTORY:
;       Written: Jean-Luc Starck 2005.
;       September, 2005 File creation
;-
;-----------------------------------------------------------------


pro getblock, ima, tab,  blocsizex, blocsizey
vs = size(ima)
Nx = vs[1]
Ny = vs[2]
Nbx = Nx/blocsizex 
Nby = Ny/blocsizey
Nb = Nbx * Nby
tab = fltarr(blocsizex,blocsizey,Nb,12)
; print, Nx, Ny, Nb, Nx*Ny, Nb*blocsizex*blocsizey
for f=0, 11  do begin
for j=0,Ny-1 do begin
for i=0,Nx-1 do begin
  ii = i / blocsizex
  jj = j / blocsizey
  b = jj * Nbx + ii
  ii = i mod blocsizex
  jj = j mod blocsizey
  tab(ii, jj, b, f) = ima(i,j,f)
  end
  end
  end
end

pro putblock, ima, tab,  blocsizex, blocsizey
vs = size(ima)
Nx = vs[1]
Ny = vs[2]
Nbx = Nx/blocsizex 
Nby = Ny/blocsizey
Nb = Nbx * Nby
for f=0, 11  do begin
for j=0,Ny-1 do begin
for i=0,Nx-1 do begin
  ii = i / blocsizex
  jj = j / blocsizey
  b = jj * Nbx + ii
  ii = i mod blocsizex
  jj = j mod blocsizey
  ima(i,j,f) = tab(ii, jj, b,f)  
  end
  end
end
end

function mrs_ridget, Rid, j, NormMad=NormMad, ImaMean=ImaMean, ImaMad=ImaMad
   if N_PARAMS() NE 2  then begin 
        print, 'CALLING SEQUENCE: Band = mrs_ridget( Rid, j, NormMad=NormMad, ImaMean=ImaMean, ImaMad=ImaMad )'
        goto, DONE
        end
 if j LT 0 or j GE Rid.NBRSCALE then begin
     print, 'Bad second parameter: ', j, ', Min Val = 0', ' Max Val = ', Rid.NBRSCALE-1
     goto, DONE
     end

  if not keyword_set(NormMad) then BEGIN 
    
     if Rid.bin EQ 1 then return, Rid.coef[ Rid.TabDepX[j]:Rid.TabDepX[j]+Rid.TabBandNx[j]-1,*,*] $
     else return, rid.coef.(11+j)
  END else begin
      ; if Rid.NXB * Rid.NYB LE 4 then begin
      ; print, 'ERROR: The mad normalization needs more than 4 blocks'
      ; return, -1
      ; end
       
       
      if Rid.bin EQ 1 then BEGIN
          Scale =  Rid.coef[ Rid.TabDepX[j]:Rid.TabDepX[j]+Rid.TabBandNx[j]-1,*,*]
	 vs = size(Scale)
         Nx = vs[1]
         Ny = vs[2]
	 BS = Rid.Bsize
	 if Rid.overlap GT 0 then BS = BS * 2
	 blocsizex = Ny / BS * 2
         blocsizey = Ny / BS
 	 getblock, Scale, Band,  blocsizex, blocsizey
     end else begin   
            Band = rid.coef.(11+j)
	    end
     vs = size(Band)
     Nx = vs[1]
     Ny = vs[2]
     ImaMean=fltarr(Nx,Ny)
     ImaMad=fltarr(Nx,Ny)
     for x=0,Nx-1 do BEGIN
     for y=0,Ny-1 do BEGIN
         Vect = reform(Band(x,y,*,*))
	 ImaMean[x,y] = mean(Vect)
	 ImaMad[x,y] = mad(Vect)
	 if  ImaMad[x,y]  GT 0 then Band(x,y,*,*) = (Band(x,y,*,*) - ImaMean[x,y]) / ImaMad[x,y] $
	 else Band(x,y,*,*) = (Band(x,y,*,*) - ImaMean[x,y])
     end
     end  
     return, Band   
  end
  
DONE:
  
end
