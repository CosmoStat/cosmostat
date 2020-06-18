
;+
; NAME:
;	 MR_UPRESOL
;
; PURPOSE:
;        From an image which have been compressed by mr_comp or mr_lcomp,
;        we can extract a low resolution image. The principle of MR_UPRESOL
;        is to improve the resolution by a factor two of the low resolution image,
;        using the compressed file (MRC format). Only wavelet coefficients
;        at the corresponding scale and position are extracted and decompressed.
;        If the keyword init is specified, there is no low resolution image,
;        and it must be created from the MRC file.
;        Two files are created:
;          'xx_resol.inf' : ascii file which contains information about Imag
;                 the file contains the three following lines:
;                   L_Nl  L_Nc  Nl  Nc
;                   Resol BlockSize KeepResi
;                   FirstBlockNl  LastBlockNl  FirstBlockNc  LastBlockNc
;                 where L_Nl x L_Nc is the size of the full resolution image
;                         Nl x Nc is the size of Imag
;                       Resol is the resolution level of Imag
;                       BlockSize is the block size used for the compression
;                       KeepResi indicates whether the noise map is stored
;                                in the compressed file
;                       FirstBlockNl  LastBlockNl  FirstBlockNc  LastBlockNc
;                                indicates which blocks are decompressed.
;                   
;          'xx_resol.fits': fits file which contains the image
;       Using the DELFILE keyword, these files are deleted.
;
;       The new image has a better resolution, and a size compatible with
;       the window size. It means that if the new image should be larger
;       than the window size, only a part (few blocks) of it will be decompressed.
;       The keywords Zx,Zy indicates the aera to zoom. This is possible only
;       if the compression was performed using the block option ("-C"). 
;
; CALLING SEQUENCE:
;	 MR_UPRESOL, FileName, Imag, init=init, Hd=Hd, WinSize=WinSize, $
;                    OPT=OPT, BI=BI, Zx=Zx, Zy=Zy, delfile=delfile
;	
; INPUT
;       FileName: string array: MRC file name
;
; INPUT-OUTPUT
;       Imag -- 2D float array: if init is set, it is only an ouput parameter
;                      if init is not, it is the image we want to improve 
;                      the resolution.
;                      At the end, Imag contains the image (or a part of the 
;                      image at a better resolution)
;
; INPUT KEYWORD:
;       init -- int: if set, Imag does not exist, and create the image at
;                    low resolution.
;       Zx -- int: X position in Imag where to zoom
;       Zy -- int: Y position in Imag where to zoom
;       WinSize -- int: window size
;       OPT -- string array: other option to give  to mr_upresol C++ program
;       delfile -- int: if set, the created files are deleted
;
; INPUT-OUTPUT KEYWORD:
;       Hd  --  string array: Fits header of the image Imag. At the end, it
;                      contains the fits header of the new image.
;       BI -- structure: Block information of Imag. It allows to know which
;                      blocks of the input image are decompressed, the 
;                      resolution level of Imag, the size of Imag and
;                      the size of the full resolution image.
;
;          BLOCKINFO = { 
;             L_Nl -- int: number of lines of the full resolution image  
;             L_Nc -- int: number of columns of the full resolution image
;             Nl -- int: number of lines of Imag
;             Nc -- int: number of columns of Imag
;             BlockSize -- int: Block size used in the compression
;             Resol -- int: Resolution level of Imag
;                            Resol = j ==> Pixel is Imag have a size of 2^j pixels
;                                          in the original image
;                            Resol = 0 ==> full resolution
;                            Resol = -1 ==> The noise map is added in Imag
;                                           only possible if a noise model
;                                           was used for the compression.
;            LastResol -- int: Lowest resolution level
;            KeepResi -- int : 0 if a noise map is stored in the compressed file
;                             > 0 otherwise
;            NbrBlock -- int: number of blocks which are decompressed
;            FirstBlockNl --  int: First Block number in y direction
;            FirstBlockNc --  int: First Block number in x direction
;                  FirstBlockNl and FirstBlockNc starts at zero (IDL convention)
;            LastBlockNl  --  int: Last Block number in y direction
;            LastBlockNc  --  int: Last Block number in x direction
;
; EXTERNAL CALLS:
;       mr_upresol (C++ program)
;
; EXAMPLE:
;            > mr_upresol, 'ngc2997.fits.MRC', imag, BI=BI, /init
;            > help, imag
;                   IMAG            FLOAT     = Array[32, 32]
;                decompress the low resolution image and put it Imag.
;               
;            > mr_upresol, 'ngc2997.fits.MRC', imag, BI=BI
;                   IMAG            FLOAT     = Array[64, 64]
;                improve the resolution by a factor of two.
;       
;
; MODIFICATION HISTORY:
; 	Written by:	JL Starck 13/10/98;
;
;-

      
;===================================================================

pro write_info,  NameInfo, BI
openw, Unit, NameInfo, /get_lun
cmd = strcompress(BI.L_Nl,/REMOVE_ALL)+strcompress(BI.L_Nc)+strcompress(BI.Nl)+strcompress(BI.Nc)
printf, Unit,  cmd
cmd = strcompress(BI.Resol,/REMOVE_ALL)+strcompress(BI.BlockSize)+strcompress(BI.KeepResi)
printf, Unit, cmd
cmd = strcompress(BI.FirstBlockNl,/REMOVE_ALL)+strcompress(BI.LastBlockNl)+strcompress(BI.FirstBlockNc)+strcompress(BI.LastBlockNc)
printf, Unit, cmd
close, Unit
free_lun, Unit
end
;===================================================================

pro read_info, NameInfo, BLOCKINFO
; read the information file
openr, Unit, NameInfo, /get_lun 
readf, Unit, L_Nl, L_Nc, Nl, Nc
readf, Unit, Resol, BlockSize, KeepResi
readf, Unit, FirstBlockNl, LastBlockNl, FirstBlockNc, LastBlockNc
close, Unit
free_lun, Unit


BLOCKINFO = { L_Nl: long(L_Nl), $
              L_Nc: long(L_Nc), $
	      Nl: long(Nl), $
              Nc: long(Nc), $
	      BlockSize: long(BlockSize), $
	      Resol: long(Resol), $
	      LastResol:  long(Resol), $
	      KeepResi: long(KeepResi), $
	      NbrBlock: long(LastBlockNl-FirstBlockNl+1)*long(LastBlockNc-FirstBlockNc+1), $
	      FirstBlockNl: long(FirstBlockNl), $
	      FirstBlockNc: long(FirstBlockNc), $
	      LastBlockNl: long(LastBlockNl), $
	      LastBlockNc: long(LastBlockNc)}
END

;===================================================================

pro mr_upresol,  FileName, Imag, init=init, Hd=Hd, WinSize=WinSize, OPT=OPT, $
      BI=BI, Zx=Zx, Zy=Zy, DELFILE=DELFILE


if N_PARAMS() NE 2 then BEGIN
 print, "Error: bad number of parameters ... "
 print, "CALLING SEQUENCE: mr_upresol, FileName, Imag, init=init, Hd=Hd, WinSize=WinSize, OPT=OPT, "
 print, "                              BI=BI, Zx=Zx, Zy=Zy, DELFILE=DELFILE"
 goto, done
 END

; File names
MRFIle=FileName
NameInfo='xx_resol.inf'
NameImag='xx_resol.fits'

; set the option to the C++ mr_upresol program
OptC = ' '
if keyword_set(OPT) then OptC = Opt
if keyword_set(WinSize) then begin
  OptC = OptC + ' -W ' + strcompress(string(WinSize),/REMOVE_ALL)
  end

; test if the image must be written to the disk
if not keyword_set(init) then begin
  read_info, NameInfo, OLD_BLOCKINFO
  if (OLD_BLOCKINFO.Resol NE BI.Resol) then BEGIN
     writefits, NameImag, Imag
     write_info, NameInfo, BI
  END
  OptC = OptC + ' -l ' + NameImag 
end

; Set zoom parameters
if keyword_set(Zx) then OptC = OptC + ' -x ' + strcompress(string(Zx),/REMOVE_ALL)
if keyword_set(Zy) then OptC = OptC + ' -y ' + strcompress(string(Zy),/REMOVE_ALL)

; run MR_UPRESOL
com = 'mr_upresol ' + OptC + ' ' + MRFIle  + ' ' +  NameImag
; print, com
spawn, com


if keyword_set(DELFILE) then BEGIN
 delete, NameImag
 if not keyword_set(init) then delete, NameInfo
END

; Read the decompressed image
Imag = readfits(NameImag,Hd)
read_info, NameInfo, BLOCKINFO


if keyword_set(BI) then BLOCKINFO.LastResol = BI.LastResol
BI = BLOCKINFO

; help, /struct, BI
print, "Resolution level: ", strcompress(string(BI.Resol)), $
       "    Nx = ", strcompress(string(BI.Nc)), "   Ny = ",   strcompress(string(BI.Nl))

DONE:
end
