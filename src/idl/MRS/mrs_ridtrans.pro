;+
; NAME:
;        mrs_ridtrans
;
; PURPOSE:
;		Computes the local ridgelet transform on the sphere, using the healPix pixel representation (NESTED data
;       representation). The standard ridgelet transform is applied on the 12 faces of the Healpix image. Each
;		face is an nside * nside quadrilateral image where nside is a power of two.
;
;       The output is an IDL structure.
;       A band at scale j (j=0..NBRSCALE-1) can be extracted using the function
;       function mrs_ridget(Rid, j) (ex: Scale2 = mrs_ridget(RidTrans, 2))
;       and a band can be inserted in the transformation using the routine  mrs_ridput
;       (ex:  mrs_ridput, RidTrans, Scale2, 2).
;
; CALLING:
;
;     mrs_ridtrans, Imag, RidTrans, NbrScale=NbrScale, overlap=overlap, blocksize=blocksize, Opt=Opt
;
; INPUTS:
;     Imag -- IDL array of healpix map: Input image be transformed 
;    
; OUTPUTS:
;     RidTrans -- IDL structures with the following fields:  
;         NBRSCALE  -- LONG: Number of scales of the ridgelet transform
;         COEF      --  Table of Ridgelet coefficients
;         BSIZE     --  LONG: Block size used in the ridgelet transform
;         NXB       -- LONG: Number of blocks in the x-axis direction
;         NYB       -- LONG: Number of blocks in the y-axis direction
;         OVERLAP   -- LONG: is equal to 1 if blocks are overlapping
;         TABNORM   -- FLOAT Array[0:NBRSCALE-1]: Normalization value for each scale
;
; KEYWORDS and optional inputs:
;     NbrScale : number of scales in the ridgelet decomposition. By default it is automatically estimated
;		
;	overlap : if this keyword is set, the blocks in the local ridgelet transform overlap by half there size.
;			  if not set, then the blocks do not overlap.
;
;   blocksize : long, this is the size of the square blocks on which the local ridgelet transform is computed.
;				if not set, then the blocksize is taken to be half the size of the individaul faces ie blocksize = nside/2
;				N.B. Blocksize is required to be a power of two smaller than nside which is also a power of two. 
;				There is no testing of this requirement.
;					
;   opt : string, sets options to used if the mre package is available. If mre is not available, then opt is useless, 
;		 and only procedures in mrs are used.
;
;
; EXTERNAL CALLS:
;        rid2d,  lino_grid
;
; EXAMPLE:
;       Compute the ridgelet transform of an image I with default options
;        The result is stored in Output
;               mrs_ridtrans, Imag, Output 
;         
; HISTORY:
;	Written:  Jean-Luc Starck & Yassir Moudden & Ludovic Poupard
;	September, 2005 File creation
;----------------------------------------------------------------------------------------------------------------

pro mrs_ridtrans, Imag, RidTrans, NbrScale=NbrScale, overlap=overlap, blocksize=blocksize, Opt=Opt

COMMON MR1ENV, mr1ok


;checking for the correct calling sequence
if N_PARAMS() LT 2  then begin 
        print, 'CALLING SEQUENCE: mrs_ridtrans, Imag, RidTrans, NbrScale=NbrScale, overlap=overlap, blocksize=blocksize, Opt=Opt'
        goto, DONE
        end
	
;recovering a few constants and parameters 	 	   
npix = (size(imag))[1]
nside = npix2nside(npix)
if keyword_set(blocksize) then blocksize = long(blocksize) else blocksize=long(nside/2.)			;setting default block size
if not keyword_set(NbrScale) then NbrScale = fix( alog( float(blocksize) / 4. * 3. / alog(2.)))		;setting default number of scales


;building a data cube from the 12 faces of the Healpix map
get_all_faces, Imag, CubeFace

;----------------------------------------------------------------------------------------------------------------
;applying the local ridgelet transform on the faces of the Healpix map
;----------------------------------------------------------------------------------------------------------------

										;-----------------------------
if keyword_set(mr1ok) then begin		;checking if mre is available
										;-----------------------------
	bin=1
	if not keyword_set(Opt) then Opt = ' ' 
	;if keyword_set(blocksize) then Opt = Opt + ' -b ' + STRCOMPRESS(string(blocksize), /REMOVE_ALL) 
	;if keyword_set(NbrScale) then Opt = Opt + ' -n ' + STRCOMPRESS(string(NbrScale), /REMOVE_ALL) 
	Opt = Opt + ' -b ' + STRCOMPRESS(string(blocksize), /REMOVE_ALL)  + ' -n ' + STRCOMPRESS(string(NbrScale), /REMOVE_ALL)
	if keyword_set(overlap) then Opt = Opt + ' -O '													;default is with overlapping blocks
	for f=0,11 do begin
		rid_trans, CubeFace[*,*,f], Trans, Opt=Opt 
		if f eq 0 then  Coef = fltarr(Trans.NXRID, Trans.NYRID, 12)
		Coef[*,*,f] = Trans.Coef
	endfor	

	RidTrans = { NbrScale : Trans.NbrScale, Nxrid : Trans.Nxrid, Nyrid : Trans.Nyrid, Coef : Coef, Bsize : Trans.Bsize, Nxb : Trans.Nxb, Nyb : Trans.Nyb, $
            Overlap : Trans.Overlap, HeadTrans : Trans.HeadTrans, TabDepX : Trans.TabDepX, TabBandNx : Trans.TabBandNx, TabBandNy : Trans.TabBandNy, TabNorm : Trans.TabNorm, bin : bin } 

										;-----------------------------
endif else begin						;using only procedures in mrs
										;-----------------------------
	bin=0
	
	;computing an index table for this block size, to be used for the extraction of lines 
	;in the linogram which is one step of the ridgelet transform on each block. 
	lino_grid, blocksize+1, index_table		
	
	for f=0,11 do begin
		rid2d, CubeFace[*,*,f], Trans, NbrScale, blocksize=blocksize, index_table=index_table, overlap=overlap
		if f eq 0 then  Coef = replicate(Trans, 12)
		Coef[f] = Trans
	endfor	
	
	RidTrans = { NbrScale : Trans.NbrScale, Coef : Coef, Bsize : Trans.Bsize, Nxb : Trans.Nxb, Nyb : Trans.Nyb, Overlap : Trans.Overlap, TabNorm : Trans.TabNorm, bin : bin }
	
	;in the future, try to include index_table to be used for recontruction.

endelse


DONE:

END

