;+
; NAME: 
;       MRSP_MATMASK_DECONV
;
; PURPOSE:
;        Invert the power spectra coupling due to masks (potentially different for temperature/polarization) using a standard MASTER pseudo Cl approach.
;
; CALLING:
;     PSPECMASTER = mrsp_matmask_deconv(PolImag, Mask, Lmax=Lmax,Niter=Niter,MaskP=MaskP,NoiseRea=NoiseRea,EBMAP=EBMAP,PS_MSKMAP=PS_MSKMAP,PS_MSKNOISE=PS_MSKNOISE,PS_MASK=PS_MASK,MatMask= MatMask,SpecRad=SpecRad,InvMatMask= InvMatMask,tebflag=tebflag,Positivity=Positivity,ZeroMonopDip=ZeroMonopDip,NormAlm=NormAlm, PolaFastALM = PolaFastALM,WorkEB=WorkEB,SVD=SVD,uncorr=uncorr,FastDeconv=FastDeconv,verb=verb,Timer=Timer)
;
; INPUTS:
;     PolImag -- Healpix 1D/2D array: Input Healpix (polarized) image the spectra of which are looked after
;	  Mask --  Healpix 1D array: weighting coefficients for each region in the sky (used either only for T, for TQU or for TEB, see WorkEB flag and MaskP input)
;     
;          
; INPUT KEYWORDS:
;      Lmax      : Number of spherical harmonics computed in the decomposition
;                                       (default is 3*nside, should be between 2*nside and 3*nside)
;      niter: int: number of iterations used in the iterative reconstruction (default: C++ default (40)- check DEF_NITER_ITER_INV in the C++ code)
;      EBMAP: Compute and return the EBMAP if TQU maps is used as inputs
;	   MaskP: Polarization Mask for Q/U or E/B (see WorkEb flag) if different from temperature mask (default: Mask)
;	   NoiseRea: Realisation of Noise -Gaussian, same size as PolImag, spectra to be subtracted prior to MASTER deconvolution
;
; INPUT/OUTPUT KEYWORDS:
;      EBMAP: Compute and return the EBMAP if TQU maps is used as inputs
;  	   PS_MSKMAP: Power spectra of the masked input maps
;	   PS_MSKNOISE: Power spectra of the masked noise maps
;	   PS_MASK: Power spectra of the masks
;	   MatMask: Master coupling matrices
;	   SpecRad: Spectral radii for each coupling matrices (iterative inversion)
;	   InvMatMask: inverse of coupling matrices (SVD)
;
; INPUT FLAGS:
;      tebflag: Input map is in TEB
;      Positivity: Positivity of the solution enforced (iterative approach)
;      ZeroMonopDip: Monopole and Dipole are assumed equal to zero (temperature power spectrum, both iterative and SVD)
;      NormAlm: Alm l2 normalized
;      PolaFastALM: Fast Alm Xform allowed
;      WorkEB: Mask the TEB maps and not the TQU
;      SVD: Perform inversion of coupling matrix using SVD (unregularized)
;	   uncorr: Uncorrelated noise in pixel space (means no TE, TB, EB noise subtraction, only TT,EE and BB).
;      FastDeconv: fast computation of either spectral radius or inverse matrix for EEEE,EEBB matrices
;      Timer: Time information printed
;
;
; OUTPUTS:
;     PSPECMASTER -- IDL 1D/2D array: Power spectra of the masks after deconvolution of the mask effects using MASTER  
;	  
;       
; EXTERNAL CALLS:
;       mrsp_matmask (C++ program)
;
; EXAMPLE:
;      
; HISTORY:
;       Written : Florent Sureau   2013.
;-
;-----------------------------------------------------------------
function mrsp_matmask_deconv,PolIma,Mask,Lmax=Lmax,Niter=Niter,tebflag=tebflag,Positivity=Positivity,ZeroMonopDip=ZeroMonopDip,verb=verb,NormAlm=NormAlm, PolaFastALM = PolaFastALM,WorkEB=WorkEB,EBMAP=EBMAP,MaskP=MaskP,NoiseRea=NoiseRea,uncorr=uncorr,PS_MSKMAP=PS_MSKMAP,PS_MSKNOISE=PS_MSKNOISE,PS_MASK=PS_MASK,MatMask= MatMask,InvMatMask= InvMatMask,SVD=SVD,SpecRad=SpecRad,FastDeconv=FastDeconv,Timer=Timer,Thr=Thr

if keyword_set(Timer) then start_it_time =systime(1)

MasterPspec=0d
Nside=0S
if keyword_set(tebflag) then ftqu =0S else ftqu =1S 
if keyword_set(Positivity) then fPositivity =1S else fPositivity =0S 
if keyword_set(ZeroMonopDip) then fZeroMonopDip =1S else fZeroMonopDip =0S 
if keyword_set(Verb) then fverb=1S else fverb=0S 
if keyword_set(NormALM) then fNormALM =1S else fNormALM =0S 
if keyword_set(PolaFastALM) then fPolaFastALM =1S else fPolaFastALM =0S 
if keyword_set(WorkEB) then fWorkEB =1S else fWorkEB =0S 
if keyword_set(uncorr) then funcorr =1S else funcorr =0S 
if keyword_set(FastDeconv) then ffastdec =1S else ffastdec =0S 
if keyword_set(Timer) then fTimer =1S else fTimer =0S 
if keyword_set(SVD) then fSVD =1S else fSVD =0S 
if (N_ELEMENTS(Lmax) eq 0) then Lmax=0
if (N_ELEMENTS(Niter) eq 0)then Niter=0S
if (N_ELEMENTS(Thr) eq 1) then lThr =[fix(Thr),fix(Thr)]
if (N_ELEMENTS(Thr) eq 0) then lThr =[0S,0S]

if(Niter lt 0) then Niter=0
fSaveEB=0S
EBMAP=0d

;INPUT MAPS/SPECS
if((N_ELEMENTS(PolIma) gt 0) && (N_ELEMENTS(PS_MSKMAP) eq 0)) then begin
	NpixMap=(size(PolIma,/dim))[0]
	NsideMap=fix(NPIX2NSIDE(NpixMap))
	Nside=NsideMap
	if((size(PolIma))[0] eq 1) then fOnlyT=1S $
	else if(((size(PolIma,/dim))[1]) eq 3) then fOnlyT=0S $
	else fOnlyT=1S
	fMSKMAP_PS=0S
	if(Lmax eq 0) then Lmax=3* NsideMap
	if ARG_PRESENT(EBMAP) then begin
	 	fSaveEB=1S 
	 	EBMAP=dblarr(size(PolIma,/dim))
	endif
endif else begin
	NsideMap =0S
	NpixMap =0S
	if (N_ELEMENTS(PS_MSKMAP) eq 0) then begin
		print,"Should Specify either input map or masked input map power spectra"
		goto,END_MRSP_MATMASK_DECONV
	endif
	if((size(PS_MSKMAP))[0] eq 2) then NMapspec=(size(PS_MSKMAP,/dim))[1] else NMapspec=1
	if(NMapspec eq 6) then fOnlyT=0S else fOnlyT=1S
	fMSKMAP_PS=1S
	if(Lmax gt 0) then Lmax=Min([(size(PS_MSKMAP,/dim))[0]-1,Lmax]) else Lmax=(size(PS_MSKMAP,/dim))[0]-1
endelse

;INPUT MASKS
if((N_ELEMENTS(Mask) gt 0) && (N_ELEMENTS(PS_MASK) eq 0)) then begin
	 NpixMask=(size(Mask,/dim))[0]
	 NsideMask=fix(NPIX2NSIDE(NpixMask))
	 if((NsideMap gt 0) && (NsideMap ne NsideMask)) then begin
		print,"mask and input map dimensions should agree"
		goto,END_MRSP_MATMASK_DECONV	
	 endif else Nside= NsideMask
	 fMSK_PS =0S
	 if (N_ELEMENTS(MaskP) eq NpixMask) then fMSKP=1S else fMSKP=0S
endif else begin
	fMSKP=0S
	NpixMask =(size(Mask,/dim))[0]
	NsideMask =fix(NPIX2NSIDE(NpixMask))
	if (N_ELEMENTS(PS_MASK) eq 0) then begin
		print,"Should Specify either input mask or mask power spectra"
		goto,END_MRSP_MATMASK_DECONV
	endif
	fMSK_PS =1S
	if(fOnlyT eq 0) then begin
		if(((size(PS_MASK,/dim))[1]) ne 3) then begin
			PS_MASK_OLD= PS_MASK
			PS_MASK = [[PS_MASK[*,0]],[PS_MASK[*,0]],[PS_MASK[*,0]]]
		endif
	endif
	if(Lmax gt 0) then Lmax=Min([(size(PS_MASK,/dim))[0]-1,Lmax]) else Lmax=(size(PS_MASK,/dim))[0]-1
endelse

;INPUT NOISE
if ((N_ELEMENTs(NoiseRea) gt 0)&&(N_ELEMENTS(PS_MSKNOISE) eq 0)) then begin
	print,"NOISE SUBTRACTION REQUIRED"
	fNoiseSub=1S
	fMSKNOISE_PS=0S
	NpixNoiseMap=(size(NoiseRea,/dim))[0]
	NsideNoiseMap=fix(NPIX2NSIDE(NpixNoiseMap))
	if(NsideMask le 0) then begin
		print,"BEWARE: No Noise Subtraction before MASTER DECONV (No valid mask)" 
		fNoiseSub=0S
		fMSKNOISE_PS=0S
	endif else if(NsideNoiseMap ne NsideMask) then begin
		print,"BEWARE: No Noise Subtraction before MASTER DECONV; NSides=", NsideNoiseMap," ", NsideMask 
		fNoiseSub=0S
		fMSKNOISE_PS=0S
	endif
endif else begin
	NpixNoiseMap=0S
	NsideNoiseMap=0S
	if ((N_ELEMENTS(PS_MSKNOISE) eq 0)) then begin
		print,"BEWARE: No Noise Subtraction before MASTER DECONV" 
		fNoiseSub=0S
		fMSKNOISE_PS=0S
	endif else begin
		if((size(PS_MSKNOISE))[0] eq 1) then begin
			if(fOnlyT eq 0) then begin
				print,"BEWARE: No Noise Subtraction before MASTER DECONV (only T noise for polarization data)" 
				fNoiseSub=0S
				fMSKNOISE_PS=0S
			endif
			if(Lmax gt 0) then Lmax=Min([(size(PS_MSKNOISE,/dim))[0]-1,Lmax]) else Lmax=(size(PS_MSKNOISE,/dim))[0]-1
		endif else begin
			if((size(PS_MSKNOISE))[1] ne 6) then begin
				print,"BEWARE: No Noise Subtraction before MASTER DECONV (should have 1 or 6 spectra)" 
				fNoiseSub=0S
				fMSKNOISE_PS=0S
			endif else begin 	
				fNoiseSub=1S
				fMSKNOISE_PS=1S
				if(Lmax gt 0) then Lmax=Min([(size(PS_MSKNOISE,/dim))[0]-1,Lmax]) else Lmax=(size(PS_MSKNOISE,/dim))[0]-1
			endelse
		endelse
	endelse
endelse

;MATRICES AND RELATED PARAMETERS
if(N_ELEMENTS(InvMatMask) gt 0) then begin
	if(fOnlyT) then begin
		if(((size(InvMatMask))[0] eq 3) && ((size(InvMatMask,/dim))[0] eq (size(InvMatMask,/dim))[1]) && ((size(InvMatMask,/dim))[2] eq 5)) then begin 
			if(Lmax gt 0) then Lmax=Min([(size(InvMatMask,/dim))[0]-1,Lmax]) else Lmax=(size(InvMatMask,/dim))[0]-1
			fCoupling =0S
			fInvCoupling=1S
			fSpecRadii=0S
			fSVD=1S
		endif else begin ;Not compatible matrix 
			print,"INV MATRIX DO NOT HAVE THE SPECIFIED DIMENSIONS: [Nls,Nls]"
			goto,END_MRSP_MATMASK_DECONV
		endelse
	endif else begin
		if(((size(InvMatMask))[0] eq 3) && ((size(InvMatMask,/dim))[0] eq (size(InvMatMask,/dim))[1]) && ((size(InvMatMask,/dim))[2] eq 5)) then begin 
			if(Lmax gt 0) then Lmax=Min([(size(InvMatMask,/dim))[0]-1,Lmax]) else Lmax=(size(InvMatMask,/dim))[0]-1
			fCoupling =0S
			fInvCoupling=1S
			fSpecRadii=0S
			fSVD=1S
		endif else begin
			print,"INV MATRIX DO NOT HAVE THE SPECIFIED DIMENSIONS: [Nls,Nls,5]"
			goto,END_MRSP_MATMASK_DECONV
		endelse
	endelse
endif else begin
	fInvCoupling=0S
	if(N_ELEMENTS(MatMask) gt 0) then begin
		fCoupling=1S
		if(fOnlyT) then begin
			if((size(MatMask,/dim))[0] eq (size(MatMask,/dim))[1]) then begin 
				if(Lmax gt 0) then Lmax=Min([(size(MatMask,/dim))[0]-1,Lmax]) else Lmax=(size(MatMask,/dim))[0]-1
			endif else begin
				print,"MATRIX DO NOT HAVE THE SPECIFIED DIMENSIONS: [Nls,Nls]"
				goto,END_MRSP_MATMASK_DECONV
			endelse
		endif else begin
			if(((size(MatMask))[0] eq 3) && ((size(MatMask,/dim))[0] eq (size(MatMask,/dim))[1]) && ((size(MatMask,/dim))[2] eq 5)) then begin 
				if(Lmax gt 0) then Lmax=Min([(size(MatMask,/dim))[0]-1,Lmax]) else Lmax=(size(MatMask,/dim))[0]-1
			endif else begin
				print,"MATRIX DO NOT HAVE THE SPECIFIED DIMENSIONS: [Nls,Nls,5]"
				goto,END_MRSP_MATMASK_DECONV
			endelse
		endelse
	endif else fCoupling=0S
	if((N_ELEMENTS(SpecRad) gt 0)) then begin
		fSpecRadii=1S
		if((fOnlyT eq 0)&& (N_ELEMENTS(SpecRad) ne 6)) then begin
			print,"Should specify spec rad for TTTT,EEEE, EEBB, TETE, EBEB, EB4Blocks"
			goto,END_MRSP_MATMASK_DECONV
		endif
	endif else fspecRadii=0S
endelse

Lmax=fix(Lmax)
if(Lmax eq 0) then begin
	print,"Problem setting power spectra/mixing matrices: Lmax should not be zero at that stage"
	goto, END_MRSP_MATMASK_DECONV
end else print,"USE LMAX=",Lmax

;Create Structure with all flags
outert=lThr[0]
innert=lThr[1]
MATMASK_STRUC={IDL_MATMASK_STRUC, Nside: NsideMap, Nested:1S,Verbose:fverb,NormALM:fNormALM,PolaFastALM:fPolaFastALM, OnlyT:fOnlyT, SaveEB:fSaveEB,WorkEB:fWorkEB,MSKP: fMSKP, $
	NoiseSub: fNoiseSub, UncorrNoise:funcorr, MSKMAP_PS:fMSKMAP_PS, MSKNOISE_PS: fMSKNOISE_PS,MSK_PS:fMSK_PS, $
	Coupling:fCoupling,InvCoupling:fInvCoupling,FastEst:ffastdec, SpecRadii:fSpecRadii,SVD:fSVD, $
	Positivity: fPositivity, ZeroMonopDip: fZeroMonopDip,Timer:fTimer, Lmax:Lmax,NIT_INV: Niter,TQUFLAG: ftqu, $
	outthr: outert,innthr: innert}

;Now make sure every input is correctly allocated
Nls=Lmax+1

;Start With MAP
if(size(PolIma,/type) ne 5) then begin
	PolIma_OLD= PolIma
	PolIma=dblarr(size(PolIma,/dim))+PolIma
endif
if((fNoiseSub)&&(size(NoiseRea,/type) ne 5)) then begin
	NoiseRea_OLD= NoiseRea
	NoiseRea =dblarr(size(NoiseRea,/dim))+ NoiseRea
endif
if(size(Mask,/type) ne 5) then begin
	Mask_OLD= Mask
	Mask =dblarr(size(Mask,/dim))+ Mask
endif
if((fMSKP)&&(size(MaskP,/type) ne 5)) then begin
	MaskP_OLD= MaskP
	MaskP =dblarr(size(MaskP,/dim))+ MaskP
endif

;Then POWER SPECTRA
if(fOnlyT) then MasterPspec = dblarr(Nls) else MasterPspec = dblarr(Nls,6)
if(N_ELEMENTS(PS_MSKMAP) eq 0) then begin
	if(fOnlyT) then PS_MSKMAP = dblarr(Nls) else PS_MSKMAP = dblarr(Nls,6)
endif else if (size(PS_MSKMAP,/type) ne 5) then begin
	PS_MSKMAP_OLD= PS_MSKMAP
	PS_MSKMAP =dblarr(size(PS_MSKMAP,/dim))+ PS_MSKMAP
endif
if(N_ELEMENTS(PS_MASK) eq 0) then begin
	if(fOnlyT) then PS_MASK = dblarr(Nls) else PS_MASK = dblarr(Nls,3)
endif else if (size(PS_MASK,/type) ne 5) then begin
	PS_MASK_OLD= PS_MASK
	PS_MASK =dblarr(size(PS_MASK,/dim))+ PS_MASK
endif
if(N_ELEMENTS(PS_MSKNOISE) eq 0) then begin
	if(fOnlyT) then PS_MSKNOISE = dblarr(Nls) else PS_MSKNOISE = dblarr(Nls,6)
endif else if (size(PS_MSKNOISE,/type) ne 5) then begin
	PS_MSKNOISE_OLD= PS_MSKNOISE
	PS_MSKNOISE =dblarr(size(PS_MSKNOISE,/dim))+ PS_MSKNOISE
endif

;Then Matrices and related parameter
if(fInvCoupling eq 0) then begin
	if(fSVD eq 1S) then  begin
		if(fOnlyT) then InvMatMask =dblarr(Nls,Nls) else InvMatMask = dblarr(Nls,Nls,5) 
		if((N_ELEMENTS(SpecRad) eq 0)) then SpecRad=0d
	endif else begin
		if((N_ELEMENTS(SpecRad) eq 0)) then SpecRad=dblarr(6)
		if((N_ELEMENTS(InvMatMask) eq 0)) then InvMatMask=0d
	endelse
	if(fCoupling eq 0) then begin
		if(fOnlyT) then MatMask =dblarr(Nls,Nls) else MatMask = dblarr(Nls,Nls,5) 
	endif
endif else begin
	if((N_ELEMENTS(MatMask) eq 0)) then MatMask=0d
	if((N_ELEMENTS(SpecRad) eq 0)) then SpecRad=0d
endelse

;int idl_cpp_wrapper_mrsp_matmask(IDL_MATMASK_STRUC* StrucIn, double* InputMap, double* InputNoiseMap,double* MaskT,double* MaskP, double *MasterPspec, double *Mat,double *InvMat,double *SpecRad,double* PSPEC_MSKMAP,double* PSPEC_MSKNOISE, double* PSPEC_MASKS,double *MapEB) {
   tvs,PolIma[*,0]
   tvs,PolIma[*,1]
   tvs,PolIma[*,2]
	if keyword_set(Timer) then 	begin
		it_time=systime(1)
		print,"IT TIME to get to call: ",it_time-start_it_time
	endif
	LOBJ=FILE_SEARCH(getenv("ISAP")+'/cxx/sapdev/build/libmrsp_idl*')
	zero=0d
	res=CALL_EXTERNAL(LOBJ[0], 'mrsp_idl_fast_matmask_interface', MATMASK_STRUC, PolIma,(fNoiseSub)? NoiseRea :zero, (N_ELEMENTS(Mask) gt 0) ? Mask: zero , (fMSKP) ? MaskP: zero, MasterPspec, MatMask ,InvMatMask, SpecRad, PS_MSKMAP, PS_MSKNOISE, PS_MASK, EBMAP,/CDECL)

   fs_plot, PS_MASK
	if keyword_set(Timer) then 	begin
		it_time=systime(1)
		print,"IT TIME to get to end program: ",it_time-start_it_time
	endif

END_MRSP_MATMASK_DECONV:

 if(N_ELEMENTS(PS_MASK_OLD) gt 0) then PS_MASK = PS_MASK_OLD ;Case T mask spectra duplicated twice
 if(N_ELEMENTS(PolIma_OLD) gt 0) then PolIma = PolIma_OLD ;Case not dblarr
 if(N_ELEMENTS(NoiseRea_OLD) gt 0) then NoiseRea = NoiseRea_OLD ;Case not dblarr
 if(N_ELEMENTS(Mask_OLD) gt 0) then Mask = Mask_OLD ;Case not dblarr
 if(N_ELEMENTS(MaskP_OLD) gt 0) then MaskP = MaskP_OLD ;Case not dblarr
 if(N_ELEMENTS(PS_MSKMAP_OLD) gt 0) then PS_MSKMAP = PS_MSKMAP_OLD ;Case not dblarr
 if(N_ELEMENTS(PS_MASK_OLD) gt 0) then PS_MASK = PS_MASK_OLD ;Case not dblarr
 if(N_ELEMENTS(PS_MSKNOISE_OLD) gt 0) then PS_MSKNOISE = PS_MSKNOISE_OLD ;Case not dblarr


return, MasterPspec

 
end

;==============================================================================================
 
