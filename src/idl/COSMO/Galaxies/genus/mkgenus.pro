;+
; NAME:
;        MKGENUS
;
; PURPOSE:
;	Computes the Genus curve
;
; CALLING:
;
;      GENUS, Data, Nu, G, OPT=Opt
;       
;
; INPUTS:
;     Data -- 3D IDL array: Input data cube 
;    
; OUTPUTS:
;     Nu -- 1D array: X-axis coordinate of the genus
;      G -- 1D array: genus
;
; KEYWORDS:
;
;      Opt -- string: string which contains the differents options. Options are:
;
;        where options = 
;
;          [-n VfStep]
;              VfStep.Default is 0.02 
;
;          [-P]
;              Periodic data. Default is nonperiodic. 
;
;          [-c Std]
;              Convolve the input data with a Gaussian width sigma=Std.
; 
;          [-l Lambda]
;             Convolve the input data with a Gaussian width sigma=sqrt(2)Lambda.
;
;          [-p]
;             Poisson distribution simulation.
;
;          [-g]
;            Gaussian field distribution simulation.
;
;          [-G IndexVal]
;            Power index in the simulated Gaussian field distribution.
;            Default is -1.
;
;          [-D Dimension]
;            Cube size dimension in the simulation.
;            Default is 32.
;
;          [-s Step]
;             Bin size. For ASCII catalogue only.
;             default is 1.000000. 
;
;          [-v]
;             Verbose.
;
;
; EXTERNAL CALLS:
;       genus (C++ program)
;
; EXAMPLE:
;
;       Compute the genus curve, convolving the data by a Gaussian
;       with standard deviation equals to 2.
;              MKGENUS, Cube, Nu, G, OPT='-c 2'
;
; HISTORY:
;-

pro mkgenus, Data, nu, g, OPT=Opt, plot=plot, filename=filename, simu=simu

if N_PARAMS() LT 2 then begin 
        spawn, 'genus'
        print, 'CALLING SEQUENCE: mkgenus, Data,  nu, g, OPT=Opt'
        goto, DONE
        end

if not keyword_set(filename) and not keyword_set(simu) then begin
 vsize = size(Data)
 if vsize(0) NE 3 then begin
        print, 'ERROR: bad first parameter ...'
        print, 'CALLING SEQUENCE: mkgenus, Data,  nu, g, OPT=Opt'
	goto, DONE
        end
 NameCube='xx_input.fits'
 writefits, NameCube, Data
end else NameCube='xx_input.fits'

filename = 'xx_temp.fits'

if not keyword_set(Opt) then Opt = ' '

com = 'genus -v ' + Opt + ' ' + NameCube + ' ' +  filename
print, com
spawn, com
Minkov = readfits(filename)
nu = Minkov(*,6)
g = Minkov(*,3)

if keyword_set(plot) then plot, nu, g

delete, NameCube
delete, filename
DONE:

END
