;+
; NAME: 
;       inpainting_fft3d
;
; PURPOSE: 
;       Sparse Inpainting of 3D data using  the Fourier dictionary.
;       Two algorithms are implemented:
;                  type=0:    X^(n+1) =  A^(-1) Thresh(  A[Xn+(Y-Xn)M] )     (Default case)
;                  type=1:    X^(n+1) =  X^n +  A^(-1) Thresh(  A[(Y-Xn)M])
;        where A is the Fourier transform operator, and Thresh is the hard thresholding operator.
;       
; CALLING:
;       Rec = inpainting_fft3d(cube,  niter=niter,  type=type,  linear=linear)
;
; INPUTS:
;       Cube: IDL 3D array --  3D cube with missing data (marked as zeros).
;
; INPUT KEYWORDS:
;           niter:  long   -- Number of iterations. Default is 20.
;           type: int -- type of algoritm (0 or 1).  Default is 0. 
;           linear: int -- linear decreasing of the threshold parameter. By default, it is an exponential decreasing.  
;
; OUTPUTS:
;           REC: IDL 3D array --  3D inpainted cube.
;
;; EXAMPLE:
;      RecCube = inpainting_fft3d(DataCube)
;-

FUNCTION inpainting_fft3d, cube, niter=niter, type=type, linear=linear, plotim=plotim,postscript=postscript
;Take a 3D cube with missing data (marked as zeros) in input 
;Return a 3D image reconstructed
;Algorithm type
;    type=0: X^(n+1)=A^(-1)TA[Xn+(Y-Xn)M] (Default case)
;    type=1: X^(n+1)=X^n+A^(-1)TA[(Y-Xn)M]

IF NOT(KEYWORD_SET(niter)) THEN niter=20
IF NOT(KEYWORD_SET(type )) THEN type=0 ;Algorithm type


size_cube=SIZE(cube)
nx=size_cube(1) & ny=size_cube(2) & nz=size_cube(3)

mask=fltarr(nx,ny,nz) & mask(*,*,*)=0.
a=where(cube EQ 0) 
mask(a)=1.

fft_cube=FFT(cube)

lambda=max(abs(fft_cube))
a=where(abs(fft_cube) LT lambda)

fft_cube(a)=0.

dcubei=FFT(fft_cube,/inverse)
cubei=cube+mask*dcubei

FOR iter=2,niter DO BEGIN

    print,iter

    IF (type EQ 0) THEN BEGIN
        cubetmp=cube+dcubei*mask
    ENDIF ELSE BEGIN
        cubetmp=cubei-dcubei
    ENDELSE

    fft_cubetmp=FFT(cubetmp)

    IF KEYWORD_SET(linear) THEN BEGIN
        lambda_i=lambda/iter
    ENDIF ELSE BEGIN
        lambda_i=lambda*(1.-erf(2.8*float(iter)/(float(niter))))
    ENDELSE
    a=where(abs(fft_cubetmp) LT lambda_i)

    fft_cubetmp(a)=0.
    dcubetmp=FFT(fft_cubetmp,/inverse)
    IF (type EQ 0) THEN BEGIN
        dcubei=dcubetmp
    ENDIF ELSE BEGIN
        dcubei=dcubei+dcubetmp
    ENDELSE
    cubei=cube+mask*dcubei

ENDFOR

; help,cube,cubei

IF KEYWORD_SET(plotim) THEN BEGIN
    x=findgen(nx) & y=findgen(ny) & z=findgen(nz)
    visu_2d,cube (*,*,0),x,y          ,scale=[3,3],postscript=postscript,filename='cube_init.eps'
    visu_2d,cubei(*,*,0),x,y,iwindow=1,scale=[3,3],postscript=postscript,filename='cube.eps'
    image=fltarr(nx,nz) & image(*,*)=cube(*,0,*)
    visu_2d,image,x,z,iwindow=2,scale=[3,3]
    image(*,*)=cubei(*,0,*)
    visu_2d,image,x,z,iwindow=3,scale=[3,3]
ENDIF

RETURN, cubei
END
