

!quiet=1

Hp=getenv('HEALPIX')
HealpixOK=0
if Hp NE '' then HealpixOK=1
if HealpixOK EQ 1 then HEALPix_path = '+' + '$HEALPIX/src/idl'
Gp=getenv('GLESP')
GlespOK=0
if Gp NE '' then GlespOK =1

ISAP_path = '+' + '$ISAP/idl'
UserPath = '.'  
PLANCK_path = '+' + '$ISAP/planck_idl'
ICosmo_path = '$ISAP'+ '/idl/extern/icosmo_v1.1'

if GlespOK EQ 1 then !PATH =   EXPAND_PATH("+$GLESP/idl")  +  ':' + !PATH 
if HealpixOK EQ 1 then !PATH =  EXPAND_PATH(HEALPix_path) +  ':' + !PATH 
!PATH = EXPAND_PATH(UserPath)+ ':' + ':' +  !PATH  + ':' + EXPAND_PATH(ISAP_path) 

; !help_path= '$MRS/help:' + !help_path 
if GlespOK EQ 1 then  !help_path= '$ISAP/idl/HELP:' +  '$GLESP:' + !help_path $
else   !help_path= '$ISAP/idl/HELP:' +  !help_path

; resize IDL memory in order to have enough memory for compilation
; .size 65000 20000

; SET the MRS variable to $ISAP
MRS_PATH = getenv('ISAP')
CMD = 'MRS=' + MRS_PATH
SETENV, CMD
; prompt definition

print, " "
print, "**********************************************************"
print, "** Interactive Sparse Astronomical Data Analysis (V3.0) **"
print, "**********************************************************"
print, " "

print, "HELP: type isaph for the iSAP help ... "

; LINUX DELL LAPTOP
; device, TRUE_COLOR=64, RETAIN=4
; window,/pixmap & wdelete
; device,bypass_translation=0

; LINUX NEC  LAPTOP
device, TRUE_COLOR=24, RETAIN=4, DECOMPOSED=0
mapref = 'MAP_REF=' + getenv('ISAP')+'/idl/extern/mapview/ref'
setenv, mapref

!PROMPT='ISAP> '

astrolib
COMMON MR1ENV, MR1OK, HealpixCXX, DEF_ALM_FAST, DEF_ALM_NITER, DEF_NORM_POWSPEC, ISAPCXX, BIN_ISAPCXX
MR1OK=0
HealpixCXX=1
ISAPCXX=1
BIN_ISAPCXX = '$ISAP/bin'
DEF_ALM_FAST=1
DEF_ALM_NITER=10
DEF_NORM_POWSPEC=0

COMMON SMICAenv , Statistics, number_s, number_c, number_b
COMMON C_PLANCK, P_NbrCannels, P_nside, P_nuGhz, P_nuMicron, P_nuMeter, P_Tcmb, P_h, P_k, P_c, P_Fwhm, P_Lmax, P_EffBeam_Min, P_EffBeam_Max, P_SigmaNoise_mKAntenna, P_ALLNbrCannels, P_TEAM, P_WG2_OUTVERS, P_OUTRES, P_NoiseFileName, P_DataFileName, P_DIR_INPUT_MAP, P_DIR_INPUT_NOISEMAP, P_PowSpec_NbrStep, P_PowSpec_From, P_PowSpec_To, P_PowSpec_Step, P_WAIT

COMMON STAR2D_DECONV,  Huv, Htuv, Starlet_deconv_Nscale, Starlet_deconv_Gen2, Starlet_deconv_Nx, Starlet_deconv_Ny
 
;======  Initialization for PLANCK

;const= {c:299792458d0, h:6.6260693d-34,$
;hbar:1.05457168d-34,k:1.3806505d-23,g:6.6742d-11,msol:1.98844d30,tcmb:2.725}
;defsysv, '!const', const

P_WAIT = 0
P_nside=2048
P_NbrCannels = 9
P_ALLNbrCannels = 14
P_nuGHz=[30,44,70,100,143,217,353,545,857, $
         23,33,41,61,94]  ; 'K','Ka','Q','V','W'  WMAP Freq   
	 
P_c=2.99792458e8                 ; m / s
P_k = 1.38066e-23 ; J/K (constante de Boltzmann)
P_h = 6.62617e-34  ; J.s (constante de Planck
P_Tcmb=2.726                     ; K
P_nuMicron = P_c / (P_nuGHz*10.^9.) * 10.^6
P_nuMeter = P_c / (P_nuGHz*10.^9.)
P_Fwhm = [33,24,14,10,7.1,5,5,5,5,   $   ; PLANCK PSF in arc min
           0.88*60.,0.66*60.,0.51*60.,0.35*60.,0.22*60.]  ; 'K','Ka','Q','V','W'   WMAP PSF in arc min (52.8, 39.6, 30.6, 21., 13.2)
P_Lmax = 3200
P_EffBeam_Min = [500,500,1200,1200, 3000, 3000, 3000, 3000, 3000,  300, 400, 500, 800, 1200]
P_EffBeam_Max = [650,650,1400,1400, 3200, 3200, 3200, 3200, 3200,  450, 550, 650,1000,1400]
P_SigmaNoise_mKAntenna = [1027.,1434.,2383., 1245.,753.6,609.1,424.5, 154.9, 71.8, $
                          1419.6, 1423.3, 2097.1,  2841.1,  5208.5] / 1000. ; Noise RMS in milliK Antenna unit

P_PowSpec_NbrStep = 7
P_PowSpec_From  = [2,  11,  31, 151, 421,  1201, 2501]
P_PowSpec_To    = [10, 30, 150, 420, 1200, 2500, 3000]
P_PowSpec_Step  = [1,   2,   5,  10,  20,  50,  100]

P_TEAM = 'sap-saclay'
P_WG2_OUTVERS = 'v1'
P_OUTRES = '.'
  
;P_NoiseFileName = P_CH2_V1_NoiseFileName
;P_DataFileName = P_CH2_V1_DataFileName
;P_DIR_INPUT_MAP = P_CH2_V1_DIR_INPUT_MAP
;P_DIR_INPUT_NOISEMAP =P_CH2_V1_DIR_INPUT_NOISEMAP

;================== Compile files =========================

.r alias.pro
.r mrs_alias.pro
.r at1dwt.pro
.r at2dwt.pro
.r conjgrad_fit.pro
.r bwt01_direct.pro
.r bwt01_inverse.pro
.r fastica.pro
.r get_stat.pro
.r dct.pro
.r invrid2d.pro
.r jade.pro
.r rid2d.pro
.r test_data_sph.pro
.r survival.pro
.r at2dwt
 .r read_cambcl
 .r star1d
 
.r alias.pro
.r mrs_alias.pro
.r at1dwt.pro
.r at2dwt.pro
.r conjgrad_fit.pro
.r bwt01_direct.pro
.r bwt01_inverse.pro
.r fastica.pro
.r get_stat.pro
.r dct.pro
.r invrid2d.pro
.r jade.pro
.r rid2d.pro
.r test_data_sph.pro
.r survival.pro
.r at2dwt
.r read_cambcl
.r star1d

;=============== Initalization iComo info =================

generalpath = ICosmo_path
newpath = strcompress(generalpath+'/general/')
!PATH=!PATH+':'+expand_path(newpath)+':' 
 
newpath = strcompress(generalpath+'/cosmo/')
!PATH=!PATH+':'+expand_path(newpath)+':' 
 
newpath = strcompress(generalpath+'/bao/')
!PATH=!PATH+':'+expand_path(newpath)+':' 
 
newpath = strcompress(generalpath+'/lensing/')
!PATH=!PATH+':'+expand_path(newpath)+':' 
 
newpath = strcompress(generalpath+'/plotting/')
!PATH=!PATH+':'+expand_path(newpath)+':' 
 
newpath = strcompress(generalpath+'/sne/')
!PATH=!PATH+':'+expand_path(newpath)+':' 
 
newpath = strcompress(generalpath+'/fisher/')
!PATH=!PATH+':'+expand_path(newpath)+':' 
 
newpath = strcompress(generalpath+'/icosmo_tests/')
!PATH=!PATH+':'+expand_path(newpath)+':' 
 
newpath = strcompress(generalpath+'/expt/')
!PATH=!PATH+':'+expand_path(newpath)+':' 
 
newpath = strcompress(generalpath+'/docs_help_and_eg/')
!PATH=!PATH+':'+expand_path(newpath)+':' 

newpath = strcompress(generalpath+'/camb/')
!PATH=!PATH+':'+expand_path(newpath)+':' 

print, 'IDL LIBRARIES: Sparse2D V1.0 MSVST V1.0 MRS V3.2, MRSP V1.0  DarthFader V1.0 ISW V1.0 iCosmo V1.1, Astrolib Vers Dec2000, Mapview V2.0 Master V1.0'
print, 'HEALPix IDL Lib = ', getenv('HEALPIX')

;=============== 

