;+
; NAME:
;        MKGAUSS
;
; PURPOSE:
;	Simulate a Gaussian Random Field
;
; CALLING:
;
;      MKGAUSS, Data, OPT=Opt
;       
;
; INPUTS:
;    
; OUTPUTS:
;     Data -- 3D IDL array: Input data cube 
;
; KEYWORDS:
;
;      Opt -- string: string which contains the differents options. Options are:
;
;        where options = 
;
;          [-c Std]
;              Convolve the simulated data with a Gaussian width sigma=Std.
; 
;          [-A SpectralAmplitude]
;              Spectral amplitude. Default is  1.
;
;          [-G IndexVal]
;            Power index in the simulated Gaussian field distribution.
;            Default is -1.
;
;          [-D Dimension]
;            Cube size dimension in the simulation.
;            Default is 32.
;
;          [-I InitRandomVal]
;            Value used for random value generator initialization.
;
;          [-v]
;             Verbose.
;
;
; EXTERNAL CALLS:
;       simgauss (C++ program)
;
; EXAMPLE:
;
;       Compute a Gaussian Random Field with dimension 64^3
;              MKGAUSS, Cube, OPT='-D 64'
;
; HISTORY:
;-

pro mkgauss, Data, OPT=Opt

if N_PARAMS() LT 1 then begin 
        spawn, 'simgauss'
        print, 'CALLING SEQUENCE: mkgauss, Data,  OPT=Opt'
        goto, DONE
        end

if not keyword_set(Opt) then Opt = ' '

NameCube='xx_output.fits'
 
com = 'simgauss -v ' + Opt + ' ' + NameCube
spawn, com
Data = readfits(NameCube)

delete, NameCube
DONE:

END
