;
;+
; NAME:
;	XLIVE
;
; PURPOSE:
;	 XLIVE is a widget program for large image analysis. The large image
;        have to be compressed  by mr_comp or mr_lcomp before being 
;        analysed. Then XLIVE data format is the MRC format.
;        When a MRC file is read, an image at very low resolution is displayed
;        and the user can improve the resolution of the image, or of a part
;        the image using the RESOLUP button. When RESOLUP is called,
;        XLIVE read in the MRC file the  wavelet cofficient needed for improving
;        the resolution. 
;        If the large image has been compressed by block (-C option), only  
;        needed blocks are decompressed, and the image in memory has always
;        a size compatible with the window size.
;        The Dat parameter (if given) contains the last displayed image.
;        If it is not given, the user can however get it using the global
;        variable XIMA.
;
;        The following operations can be done by pressing mouse buttons:
;           QUIT: quit the application
;           LOAD: load an image. Image formats can be either fits or midas
;                 format
;           LUT: modify the LUT
;           PROFILE: examine lines or columns 
;           CURSOR: examine pixel values. If image format is fits and if
;                   the header contains astrometric position, then pixel
;                   position in the sky is given (Right ascension and
;                   declinaison.
;           HISTO: plot the histogram
;           INFO: print the min, max, mean, sigma of the image
;           FFT: computes the Fourier transform and display either the power 
;                spectrum, the phase, the real part or the imaginary part.
;           PAN: make a zoom. Zoom factor are 2,4,8,1/2,1/4,1/8
;           RESOLUP: Improve the image resolution. If the new image size is
;                greater than the window size, the user must click in the aera
;                he want to see, and only this part of the image will be 
;                decompressed.
;           RESOLDOWN: decrease the resolution.
;           RESOL: goto a given resolution.
;
; CALLING SEQUENCE:
;	XLIVE, [Dat,], FILEName=FileName, WindowSize=WindowSize
;	
; KEYWORD:
;       MRCFileName -- string array :	MRC file name
;       WindowSize  -- int: window size. Default is 512.
;
; EXTERNAL CALLS:
;       mr_upresol
;
; EXAMPLE:
;               XLIVE
;                The MRC file will be read using the LOAD button
;
;		XLIVE, FileName='ngc2997.fits.MRC', WindowSize=256
;                run XLIVE. Analys the image ngc2997.fits with a window 
;                of size 256.
;
; MODIFICATION HISTORY:
; 	Written by:	JL Starck 13/10/98
;
;-

;===================================================================

pro common_xlive
COMMON  XLIVE_BLOCK, WIDGET_SET, MAIN1, DRAW13, TEXT11,$
         WDATA, Imag, FitsHeader, KFits, Zoom, ImagZoom, TabWin, MRCFile, $
	 BLOCKINFO, WinSize, FirstResol
	 
COMMON  XRESOL,	ln_ima, l0_ima, l1_ima, l2_ima, l3_ima, l4_ima, l5_ima, l6_ima, l7_ima, l8_ima, l9_ima, l10_ima, $
         BILN, BIL0, BIL1, BIL2, BIL3, BIL4, BIL5, BIL6, BIL7, BIL8, BIL9, BIL10, $
         HDn_hd, HD0_hd, HD1_hd, HD2_hd, HD3_hd, HD4_hd, HD5_hd, HD6_hd, HD7_hd, HD8_hd, HD9_hd, HD10_hd
 
COMMON XLIVE, XIma

COMMON XINFO, id_menu_resol, id_menu_down, id_menu_up, id_menu_pan, id_menu_fft, $
              id_menu_load, id_menu_quit, id_menu_info, id_menu_plot, id_menu_get, id_menu_print, id_menu_lut
end

;===================================================================

pro sensitive_init, keepone=keepone, all=all
COMMON  XLIVE_BLOCK
COMMON  XINFO
COMMON XRESOL

if keyword_set(keepone) then ko = keepone else ko = 0

if keyword_set(all) then WIDGET_CONTROL,id_menu_load[0],sensitive=0
if keyword_set(all) then WIDGET_CONTROL,id_menu_quit[0],sensitive=0

if ko NE 1 then WIDGET_CONTROL,id_menu_lut[0],sensitive=0  ; happend NO
if ko NE 2 then WIDGET_CONTROL,id_menu_print[0],sensitive=0  ; happend NO
if ko NE 3 then WIDGET_CONTROL,id_menu_get[0],sensitive=0  ; happend NO
if ko NE 4 then WIDGET_CONTROL,id_menu_plot[0],sensitive=0  ; happend NO
if ko NE 5 then WIDGET_CONTROL,id_menu_info[0],sensitive=0  ; happend NO
if ko NE 6 then WIDGET_CONTROL,id_menu_fft[0],sensitive=0  ; happend NO
if ko NE 7 then WIDGET_CONTROL,id_menu_pan[0],sensitive=0  ; happend NO

if ko NE 8 then WIDGET_CONTROL,id_menu_up[0],sensitive=0  ; happend NO
if ko NE 9 then WIDGET_CONTROL,id_menu_down[0],sensitive=0  ; happend NO
if ko NE 10 then WIDGET_CONTROL,id_menu_resol[0],sensitive=0  ; happend NO

end

;===================================================================

pro sensitive_resol, Resolfrom, Resolto
COMMON  XINFO
COMMON XRESOL
COMMON  XLIVE_BLOCK

; print, "Resolfrom = ", Resolfrom,  "KR = ", BLOCKINFO.KeepResi
WIDGET_CONTROL,id_menu_resol[0], /sensitive

from = Resolfrom + 2
to = Resolto + 2
for i=1,from-1 do WIDGET_CONTROL,id_menu_resol[i], sensitive=0
for i=from,to do WIDGET_CONTROL,id_menu_resol[i], /sensitive
for i=to+1,12 do WIDGET_CONTROL,id_menu_resol[i], sensitive=0

if  BLOCKINFO.Resol EQ BLOCKINFO.LastResol then $
                         WIDGET_CONTROL,id_menu_down[0],sensitive=0 $
else WIDGET_CONTROL,id_menu_down[0],/sensitive

if (((BLOCKINFO.KeepResi GE 1) AND  (BLOCKINFO.Resol   EQ -1)) OR $
   ((BLOCKINFO.KeepResi EQ 0) AND  (BLOCKINFO.Resol   EQ 0))) then $
                          WIDGET_CONTROL,id_menu_up[0],sensitive=0 $
else WIDGET_CONTROL,id_menu_up[0],/sensitive

end

;===================================================================

pro sensitive_ima
COMMON  XLIVE_BLOCK
COMMON  XINFO
COMMON XRESOL

WIDGET_CONTROL,id_menu_load[0],/sensitive  ; happend NO
WIDGET_CONTROL,id_menu_quit[0],/sensitive  ; happend NO
WIDGET_CONTROL,id_menu_get[0],/sensitive  ; happend NO
WIDGET_CONTROL,id_menu_print[0],/sensitive  ; happend NO
WIDGET_CONTROL,id_menu_plot[0],/sensitive  ; happend NO
WIDGET_CONTROL,id_menu_info[0],/sensitive  ; happend NO
WIDGET_CONTROL,id_menu_fft[0],/sensitive  ; happend NO
WIDGET_CONTROL,id_menu_pan[0],/sensitive  ; happend NO
WIDGET_CONTROL,id_menu_lut[0],/sensitive  ; happend NO
sensitive_resol, FirstResol, BLOCKINFO.LastResol 
end

;===================================================================

pro load_ima
COMMON  XLIVE_BLOCK
COMMON XRESOL

 Zoom = 1.
 erase
 vsize = size(Imag)
 Nl = vsize[2]
 Nc = vsize[1]
 N = Nl
 if N LT Nc then N = Nc
 while N*Zoom LT 512 do Zoom = Zoom * 2
 if Zoom GT 1 and N*Zoom GT 512 then Zoom = Zoom / 2
 
 Nl = Nl*Zoom
 Nc = Nc*Zoom
 Imagzoom = congrid(Imag, Nc, Nl)
 tvscl, ImagZoom
 WIDGET_CONTROL, DRAW13, SET_DRAW_VIEW=[0,0]
end

;===================================================================

PRO PDMENU22_Event, Event
COMMON XLIVE_BLOCK
COMMON XLIVE
COMMON XRESOL

  CASE Event.Value OF 


  'QUIT': BEGIN
     ; defsysv, "!ImaLive", Imag
     Xima = imag
     widget_control, MAIN1, /destroy
    END
  ENDCASE
END
 
;===================================================================

pro copy_ima
COMMON  XLIVE_BLOCK
COMMON XRESOL

if BLOCKINFO.Resol GE 0 then C=strcompress(string(BLOCKINFO.Resol),/REMOVE_ALL) $
else C='n'
NameVarIma = 'l'+C+'_ima'
cmd = NameVarIma + ' = Imag'
ACK = EXECUTE(cmd)
; print, cmd

NameVarBI = 'BIL'+C
cmd = NameVarBI + ' = BLOCKINFO'
ACK = EXECUTE(cmd)
;print, cmd
if KFits EQ 1 then BEGIN
       NameVarHD = 'HDL'+C+'_hd'
       cmd = NameVarHD  + ' = FitsHeader'
       ACK = EXECUTE(cmd)
END
end

;===================================================================

pro get_ima, resol
COMMON  XLIVE_BLOCK
COMMON XRESOL

if Resol GE 0 then C=strcompress(string(Resol),/REMOVE_ALL) $
else C='n'
NameVarIma = 'l'+C+'_ima'
cmd = 'Imag=' + NameVarIma
ACK = EXECUTE(cmd)
; print, cmd

NameVarBI = 'BIL'+C
cmd = 'BLOCKINFO = ' +  NameVarBI
ACK = EXECUTE(cmd)
;print, cmd
if KFits EQ 1 then BEGIN
       NameVarHD = 'HDL'+C+'_hd'
       cmd = 'FitsHeader = ' + NameVarHD
       ACK = EXECUTE(cmd)
END
end

;===================================================================

PRO PDMENU23_Event, Event
COMMON  XLIVE_BLOCK
COMMON XRESOL

WIDGET_CONTROL,/HOURGLASS

  CASE Event.Value OF 
  'LOAD': BEGIN
          case !version.os of
                  'vms': filt = '*.MRC*'
                  'windows': filt = '*.MRC'
;                  'MacOS': goto, notsup
                  else: filt = '*.MRC'
               endcase

          FileName = pickfile(/READ, FILTER=filt, GET_PATH=GET_PATH)
          if strlen(FileName) GT 0 then BEGIN
               MRCFile=FileName
	       mr_upresol, MRCFile, Imag, Hd=FitsHeader, WinSize=WinSize, /init, BI=BLOCKINFO
               KFits = 0
	       FirstResol = BLOCKINFO.LastResol
               if type_of(FitsHeader) EQ 'STRING' then BEGIN
                  extast,FitsHeader,astr,noparams
                  if (noparams ge 0) then KFits = 1
                  END
             load_ima
	     copy_ima
	     sensitive_ima
             END
     END
  ENDCASE
END

;===================================================================

PRO PDMENU24_Event, Event

  CASE Event.Value OF 
  'LUT': BEGIN
    xloadct
    END
  ENDCASE
END

;===================================================================


PRO PDMENU36_Event, Event
COMMON  XLIVE_BLOCK
COMMON XRESOL

WIDGET_CONTROL,/HOURGLASS

vsize = size(imag)
if vsize[0] EQ 2 then BEGIN

  CASE Event.Value OF 
  'PRINT': BEGIN
     set_plot, 'PS' 
     device, filename='idl.ps'
     tvscl, ImagZoom
     device, /close
     set_plot, 'X'
    END
  ENDCASE
  ENDIF
END

 ;===================================================================


PRO PDMENU25_Event, Event
COMMON  XLIVE_BLOCK
COMMON XRESOL

vsize = size(imag)
if vsize[0] EQ 2 then BEGIN
  sensitive_init, /all
  CASE Event.Value OF 

  'GET.PROFILE': BEGIN
    str = 'Left mouse button to toggle between rows and columns. Left to Display and Right mouse button to Exit.'
      WIDGET_CONTROL,  TEXT11, SET_VALUE=str
    plot_profile, Imag, zoom=zoom, /order
    str = ' '
    WIDGET_CONTROL,  TEXT11, SET_VALUE=str
    END
  'GET.CURSOR': BEGIN
           Nx = (size(imag))(1)
           Ny = (size(imag))(2)
          str = [' Press left or center mouse button for new output line', $
                 '... right mouse button to exit.']
          WIDGET_CONTROL,  TEXT11, SET_VALUE=str
          if KFits EQ 0  then get_cursor, Imag, zoom=zoom $
          else  BEGIN
                get_cursor, FitsHeader, Imag, zoom=zoom
                END
          str = '  '
          WIDGET_CONTROL,  TEXT11, SET_VALUE=str
    END
  ENDCASE
  sensitive_ima
ENDIF
END

 ;===================================================================

PRO PDMENU27_Event, Event
COMMON  XLIVE_BLOCK
COMMON XRESOL

WIDGET_CONTROL,/HOURGLASS

vsize = size(imag)
if vsize[0] EQ 2 then BEGIN
  CASE Event.Value OF 
  'PLOT.HISTO': BEGIN
              str = 'Right mouse button to Exit histogram.' 
          WIDGET_CONTROL,  TEXT11, SET_VALUE=str

           wsize = .75
           window,/free ,xs=wsize*640, ys=wsize*512,title='Histogram'  
           new_w = !d.window
           plothist, Imag
           print,'Right mouse button to Exit histogram.'
           Run = 1
	   while Run do begin
	       cursor,x,y,2,/dev	;Read position
               if !err eq 4 then  Run=0
           end
           wdelete, new_w
          str = ' '
          WIDGET_CONTROL,  TEXT11, SET_VALUE=str
    END
  'PLOT.CONTOURS': BEGIN
    a = !d.window
    xcont, Imag, window=9
    TabWin[9] = 1
    END
  'PLOT.3D VISU': BEGIN
    a = !d.window
    x3view, Imag
    wset, TabWin[11]
    END
  'PLOT.LOAD IN ANOTHER WINDOW.0': BEGIN
    window, 0, xsize=vsize[1]*Zoom, ysize=vsize[2]*Zoom
    tvscl, ImagZoom
    TabWin[0] = 1
    END
  'PLOT.LOAD IN ANOTHER WINDOW.1': BEGIN
    window, 1, xsize=vsize[1]*Zoom, ysize=vsize[2]*Zoom
    tvscl, ImagZoom
    TabWin[1] = 1
    END
  'PLOT.LOAD IN ANOTHER WINDOW.2': BEGIN
    window, 2, xsize=vsize[1]*Zoom, ysize=vsize[2]*Zoom
    tvscl, ImagZoom
    TabWin[2] = 1
    END
  'PLOT.LOAD IN ANOTHER WINDOW.3': BEGIN
    window, 3, xsize=vsize[1]*Zoom, ysize=vsize[2]*Zoom
    tvscl, ImagZoom  
    TabWin[3] = 1  
    END
  'PLOT.LOAD IN ANOTHER WINDOW.4': BEGIN
    window, 4, xsize=vsize[1]*Zoom, ysize=vsize[2]*Zoom
    tvscl, ImagZoom
    TabWin[4] = 1
    END
  'PLOT.LOAD IN ANOTHER WINDOW.5': BEGIN
    window, 5, xsize=vsize[1]*Zoom, ysize=vsize[2]*Zoom
    tvscl, ImagZoom
    TabWin[5] = 1
    END
  'PLOT.LOAD IN ANOTHER WINDOW.6': BEGIN
    window, 6, xsize=vsize[1]*Zoom, ysize=vsize[2]*Zoom
    tvscl, ImagZoom
    TabWin[6] = 1
    END
  'PLOT.LOAD IN ANOTHER WINDOW.7': BEGIN
    window, 7, xsize=vsize[1]*Zoom, ysize=vsize[2]*Zoom
    tvscl, ImagZoom
    TabWin[7] = 1
    END
  'PLOT.LOAD IN ANOTHER WINDOW.8': BEGIN
    window, 8, xsize=vsize[1]*Zoom, ysize=vsize[2]*Zoom
    tvscl, ImagZoom
    TabWin[8] = 1
    END
  'PLOT.LOAD IN ANOTHER WINDOW.9': BEGIN
    window, 9, xsize=vsize[1]*Zoom, ysize=vsize[2]*Zoom
    tvscl, ImagZoom
    TabWin[9] = 1
    END
  'PLOT.LOAD IN ANOTHER WINDOW.10': BEGIN
    window, 10, xsize=vsize[1]*Zoom, ysize=vsize[2]*Zoom
    tvscl, ImagZoom
    TabWin[10] = 1
    END
  ENDCASE
    wset, TabWin[11]
  ENDIF
   wset, TabWin[11]
END
 
;===================================================================

function i_to_s, i, all=all
if keyword_set(all) then Val =  strcompress(string(i),/REMOVE_ALL) $
else Val =  strcompress(string(i))
return, val
end

 ;===================================================================

PRO PDMENU28_Event, Event
COMMON  XLIVE_BLOCK
COMMON XRESOL

vsize = size(imag)
if vsize[0] EQ 2 then BEGIN
WIDGET_CONTROL,/HOURGLASS


  CASE Event.Value OF 


  'INFO': BEGIN
  
    FBC = BLOCKINFO.FirstBlockNc 
    FBL = BLOCKINFO.FirstBlockNl
    LBC = BLOCKINFO.LastBlockNc 
    LBL = BLOCKINFO.LastBlockNl
    BS = BLOCKINFO.BlockSize
    L_Nl = BLOCKINFO.L_Nl
    L_Nc = BLOCKINFO.L_Nc
    NbrBlockNl = long(L_Nl/BS)
    if (NbrBlockNl*BS NE L_Nl) then NbrBlockNl = NbrBlockNl + 1
    NBRBLOCKNC = long(L_Nc/BS)
    
    print, " Full resolution image size: ", i_to_s(BLOCKINFO.L_NC), i_to_s(BLOCKINFO.L_NL)
    print, " Image size: ",  i_to_s(BLOCKINFO.Nc), i_to_s(BLOCKINFO.NL)
    print, " BlockSize: ", i_to_s(BS), "  Number of blocks: ", i_to_s(NBRBLOCKNC), i_to_s(NBRBLOCKNL)
    print, " First decompressed block : ", i_to_s(FBC+1),  i_to_s(FBL+1) 
    print, " Last decompressed block : ", i_to_s(LBC+1),  i_to_s(LBL+1)
    print, " Position of pixel (0,0) in the large image: ",  i_to_s(FBC*BS),  i_to_s(FBL*BS)
    print, " Information about the decompressed image part: "
    min = min(Imag)
    max = max(Imag)
    Mean = total(Imag) / N_ELEMENTS(Imag)
    Sigma = sqrt(total( (Imag-mean)^2) /N_ELEMENTS(Imag))
    str = '        min='+strcompress(string(min),/REMOVE_ALL)+ $
      ' max='+strcompress(string(max),/REMOVE_ALL)+ $
      ' mean='+strcompress(string(mean),/REMOVE_ALL)+ $
      ' sigma='+strcompress(string(sigma),/REMOVE_ALL)
    print, str
    if Zoom NE 1 then str = str + ' zoom='+strcompress(string(zoom),/REMOVE_ALL)
      WIDGET_CONTROL,  TEXT11, SET_VALUE=str

    END
  ENDCASE
ENDIF
END

 ;===================================================================

PRO PDMENU7_Event, Event
COMMON  XLIVE_BLOCK
COMMON XRESOL

vsize = size(Imag)
if vsize[0] EQ 2 then BEGIN

vsize = size(ImagZoom)
Nl = vsize[2]
Nc = vsize[1]
WIDGET_CONTROL,/HOURGLASS

  CASE Event.Value OF 


  'FFT.POWER SPECTRUM': BEGIN
        fftImag = power_spectrum(Imag)
    END
  'FFT.PHASE': BEGIN
    F = dft(Imag)
    fftImag = atan(imaginary(f), float(f))
    END
  'FFT.REAL': BEGIN
    Imag = float(dft(Imag))
    END
  'FFT.IMAGINARY': BEGIN
    fftImag = imaginary(dft(Imag))
    END
  ENDCASE
    erase
    if Zoom NE 1 then fftImagzoom = congrid(fftImag, Nc, Nl)  $
    else fftImagzoom = fftImag
    tvscl, fftImagzoom
ENDIF
END

;===================================================================

PRO PDMENU9_Event, Event
COMMON  XLIVE_BLOCK
COMMON XRESOL

vsize = size(imag)
if vsize[0] EQ 2 then BEGIN
  CASE Event.Value OF 
  'PAN.PAN 2': BEGIN
    fact = 2.
    END
  'PAN.PAN 4': BEGIN
   fact= 4.
    END
  'PAN.PAN 8': BEGIN
    fact = 8
    END
  'PAN.PAN 1/2': BEGIN
    fact = 0.5
    END
  'PAN.PAN 1/4': BEGIN
    fact = 0.25
    END
  'PAN.PAN 1/8': BEGIN
    fact = 0.125
    END
  ENDCASE
    vsize = size(ImagZoom)
    Nl = fix(vsize[2]*fact)
    Nc = fix(vsize[1]*fact)
    erase
    widget_control,/hourglass
    Imagzoom = congrid(Imag, Nc, Nl)
    zoom = zoom * fact
    tvscl, ImagZoom
END
END

 ;===================================================================
 ;   MR UP RESOL
 ;===================================================================
 
PRO PDMENU29_Event, Event
COMMON  XLIVE_BLOCK
COMMON XINFO
COMMON XRESOL

if (BLOCKINFO.resol GT 0) AND (BLOCKINFO.NBRBLOCK GT 1) AND $
   (BLOCKINFO.Nl GT WinSize/2 OR BLOCKINFO.Nc GT WinSize/2) then $
BEGIN
   print, "Click the position to zoom in the image ... "
   cursor,X0,y0,2,/dev,/down
   X0 = X0 / zoom
   y0 = Y0 / zoom
   zx=long(x0)
   zy=long(y0)
   print, 'Zoom at position: ', zx, zy
END ELSE BEGIN
   zx=0
   zy=0
END

if (BLOCKINFO.resol GT 0) OR $
	   (BLOCKINFO.resol EQ 0 AND BLOCKINFO.KeepResi  GT 0) then $
BEGIN

    widget_control,/hourglass 
    mr_upresol, MRCFile, Imag, Hd=Hd, WinSize=WinSize, BI=BLOCKINFO,zx=zx,zy=zy
    copy_ima
    FirstResol = BLOCKINFO.resol
    if type_of(Hd) EQ 'STRING' then $
    BEGIN
           extast,hd,astr,noparams
           if (noparams ge 0) then $
	   BEGIN
                KFits = 1
                FitsHeader = Hd
           END
    END
    load_ima
    sensitive_resol, BLOCKINFO.Resol, BLOCKINFO.LastResol
END ELSE  print, "The image is already decompressed ... "
END

 ;===================================================================

pro goto_resol, Resol
COMMON  XLIVE_BLOCK
COMMON XRESOL

if (Resol LT FirstResol OR Resol GT BLOCKINFO.Lastresol) then print, "Error: cannot get this resolution ... " $
ELSE BEGIN
    get_ima, Resol
    load_ima
    sensitive_resol, FirstResol, BLOCKINFO.LastResol
    print, "Resolution level: ", Resol
END  
END

;===================================================================
; DOWN resol
;===================================================================

PRO PDMENU30_Event, Event
COMMON  XLIVE_BLOCK
COMMON XINFO
COMMON XRESOL

if (BLOCKINFO.resol LT  BLOCKINFO.Lastresol) then $
BEGIN
   goto_resol, BLOCKINFO.Resol+1  
END ELSE  print, "Last resolution is already displayed ... "
 
END

;===================================================================
; GOTO resol
;===================================================================

PRO PDMENU31_Event, Event
COMMON  XLIVE_BLOCK
COMMON  XINFO
COMMON XRESOL

Resol=-1
  CASE Event.Value OF 
  'RESOL.RESOL Ima': Resol = - 1
  'RESOL.RESOL 0': Resol = 0
  'RESOL.RESOL 1': Resol = 1
  'RESOL.RESOL 2': Resol = 2
  'RESOL.RESOL 3': Resol = 3
  'RESOL.RESOL 4': Resol = 4
  'RESOL.RESOL 5': Resol = 5
  'RESOL.RESOL 6': Resol = 6
  'RESOL.RESOL 7': Resol = 7
  'RESOL.RESOL 8': Resol = 8
  'RESOL.RESOL 9': Resol = 9
  'RESOL.RESOL 10': Resol = 10
  ENDCASE
 if Resol GE FirstResol AND Resol NE BLOCKINFO.Resol then BEGIN
      goto_resol, Resol

  END
END
 
 ;===================================================================
 ;   MAIN EVENT
 ;===================================================================

PRO MAIN1_Event, Event
COMMON  XLIVE_BLOCK
COMMON XRESOL

  WIDGET_CONTROL,Event.Id,GET_UVALUE=Ev
  wset, TabWin[11]

  CASE Ev OF 

  ; Event for PDMENU22
  'PDMENU22': PDMENU22_Event, Event
  ; Event for PDMENU23
  'PDMENU23': PDMENU23_Event, Event
  ; Event for PDMENU24
  'PDMENU24': PDMENU24_Event, Event
  ; Event for PDMENU36
  'PDMENU36': PDMENU36_Event, Event
  ; Event for PDMENU25
  'PDMENU25': PDMENU25_Event, Event
  ; Event for PDMENU27
  'PDMENU27': PDMENU27_Event, Event
  ; Event for PDMENU28
  'PDMENU28': PDMENU28_Event, Event
  ; Event for PDMENU7
  'PDMENU7': PDMENU7_Event, Event
  ; Event for PDMENU9
  'PDMENU9': PDMENU9_Event, Event
  'TEXT11': BEGIN
      Print, 'Event for TEXT11'
      END
  'PDMENU29': PDMENU29_Event, Event  
  'PDMENU30': PDMENU30_Event, Event
  'PDMENU31': PDMENU31_Event, Event
  'DRAW13': BEGIN
            vsize = size(imag)
      if vsize[0] EQ 2 then BEGIN
           x = fix(event.x / zoom)
           y = fix(event.y / zoom)
           Nl = vsize[2]
           Nc = vsize[1]
           if (x ge 0 and x lt Nc) AND (y ge 0 and y lt Nl) then $
               Inten = Imag(x,y) else Inten = 0

           str =  '('+strcompress(string(x),/REMOVE_ALL)+','+$
                 strcompress(string(y),/REMOVE_ALL)+') = '+$
                 strcompress(string(Inten),/REMOVE_ALL)
           WIDGET_CONTROL,  TEXT11, SET_VALUE=str
          END
      END
  'PDMENU31': Print, 'Event for TEXT11'
  ENDCASE
END
 
;===================================================================

PRO xlive, Dat, FileName=FileName, WindowSize=WindowSize, GROUP=Group 
COMMON  XLIVE_BLOCK
COMMON XINFO
COMMON XRESOL

  astrolib
  Zoom = 1
  Imag = 0
  TabWin = intarr(12)
  KFits  = 0

  if N_PARAMS() GT 1 then BEGIN
     print, "Error: too many parameters ... "
     print, "CALLING SEQUENCE: xlive, Dat, FileName=FileName, WindowSize=WindowSize"
     goto, DONE
  END
  
  if Keyword_set(WindowSize) then  WinSize = WindowSize else  WinSize = 512

  if Keyword_set(FileName) then BEGIN 
         MRCFile=FileName
         mr_upresol, FileName, Imag, Hd=Hd, WinSize=WinSize, /init, BI=BLOCKINFO
         if (size(Imag))[0] NE 2 then begin
	   print, "Error: cannot read ", FileName
	   goto, DONE
	 END
 	 FirstResol = BLOCKINFO.LastResol
         if type_of(Hd) EQ 'STRING' then BEGIN
	          FitsHeader = Hd
                  extast,FitsHeader,astr,noparams
                  if (noparams ge 0) then KFits = 1
                  END

  END

  IF N_ELEMENTS(Group) EQ 0 THEN GROUP=0

  junk   = { CW_PDMENU_S, flags:0, name:'' }

  MAIN1 = WIDGET_BASE(GROUP_LEADER=Group, $
      /COLUMN, $
      MAP=1, $
      TITLE='XLIVE', $
      UVALUE='MAIN1', $
      XSIZE=600, $
      YSIZE=700)

  BASE3 = WIDGET_BASE(MAIN1, $
      ROW=2, $
      MAP=1, $
      UVALUE='BASE3')

  BASE6 = WIDGET_BASE(BASE3, $
      ROW=1, $
      MAP=1, $
      UVALUE='BASE6')

  BASE7 = WIDGET_BASE(BASE3, $
      ROW=1, $
      MAP=1, $
      UVALUE='BASE7')
      
  MenuDesc1116 = [ $
      { CW_PDMENU_S,       2, 'QUIT' } $  ;      0

  ]


  PDMENU22 = CW_PDMENU( BASE6, MenuDesc1116, /RETURN_FULL_NAME, $
      UVALUE='PDMENU22', IDS=id_menu_quit)

  MenuDesc1117 = [ $
      { CW_PDMENU_S,       2, 'LOAD' } $  ;      0

  ]


  PDMENU23 = CW_PDMENU( BASE6, MenuDesc1117, /RETURN_FULL_NAME, $
      UVALUE='PDMENU23', IDS=id_menu_load)

  MenuDesc1118 = [ $
      { CW_PDMENU_S,       2, 'LUT' } $  ;      0

  ]


  PDMENU24 = CW_PDMENU( BASE6, MenuDesc1118, /RETURN_FULL_NAME, $
      UVALUE='PDMENU24', IDS=id_menu_lut)

  MenuDesc1119 = [ $
      { CW_PDMENU_S,       2, 'PRINT' } $  ;      0

  ]


  PDMENU36 = CW_PDMENU( BASE6, MenuDesc1119, /RETURN_FULL_NAME, $
      UVALUE='PDMENU36', IDS=id_menu_print)

  MenuDesc1120 = [ $
      { CW_PDMENU_S,       3, 'GET' }, $ ;        0
        { CW_PDMENU_S,       0, 'PROFILE' }, $ ;        1
        { CW_PDMENU_S,       2, 'CURSOR' } $  ;      2

  ]


  PDMENU25 = CW_PDMENU( BASE6, MenuDesc1120, /RETURN_FULL_NAME, $
      UVALUE='PDMENU25', IDS=id_menu_get)

  MenuDesc1122 = [ $
      { CW_PDMENU_S,       3, 'PLOT' }, $ ;        0
        { CW_PDMENU_S,       0, 'HISTO' }, $ ;        1
        { CW_PDMENU_S,       0, 'CONTOURS' }, $ ;        2
        { CW_PDMENU_S,       0, '3D VISU' }, $ ;        3
        { CW_PDMENU_S,       3, 'LOAD IN ANOTHER WINDOW' }, $ ;        4
          { CW_PDMENU_S,       0, '0' }, $ ;        5
          { CW_PDMENU_S,       0, '1' }, $ ;        6
          { CW_PDMENU_S,       0, '2' }, $ ;        7
          { CW_PDMENU_S,       0, '3' }, $ ;        8
          { CW_PDMENU_S,       0, '4' }, $ ;        9
          { CW_PDMENU_S,       0, '5' }, $ ;       10
          { CW_PDMENU_S,       0, '6' }, $ ;       11
          { CW_PDMENU_S,       0, '7' }, $ ;       12
          { CW_PDMENU_S,       0, '8' }, $ ;       13
          { CW_PDMENU_S,       0, '9' }, $ ;       14
          { CW_PDMENU_S,       2, '10' } $  ;     15

  ]


  PDMENU27 = CW_PDMENU( BASE6, MenuDesc1122, /RETURN_FULL_NAME, $
      UVALUE='PDMENU27', IDS=id_menu_plot)

  MenuDesc1125 = [ $
      { CW_PDMENU_S,       2, 'INFO' } $  ;      0

  ]


  PDMENU28 = CW_PDMENU( BASE6, MenuDesc1125, /RETURN_FULL_NAME, $
      UVALUE='PDMENU28', IDS=id_menu_info)

  MenuDesc1126 = [ $
      { CW_PDMENU_S,       3, 'FFT' }, $ ;        0
        { CW_PDMENU_S,       0, 'POWER SPECTRUM' }, $ ;        1
        { CW_PDMENU_S,       0, 'PHASE' }, $ ;        2
        { CW_PDMENU_S,       0, 'REAL' }, $ ;        3
        { CW_PDMENU_S,       2, 'IMAGINARY' } $  ;      4

  ]


  PDMENU7 = CW_PDMENU( BASE6, MenuDesc1126, /RETURN_FULL_NAME, $
      UVALUE='PDMENU7', IDS=id_menu_fft)

  MenuDesc1128 = [ $
      { CW_PDMENU_S,       3, 'PAN' }, $ ;        0
        { CW_PDMENU_S,       0, 'PAN 2' }, $ ;        1
        { CW_PDMENU_S,       0, 'PAN 4' }, $ ;        2
        { CW_PDMENU_S,       0, 'PAN 8' }, $ ;        3
        { CW_PDMENU_S,       0, 'PAN 1/2' }, $ ;        4
        { CW_PDMENU_S,       0, 'PAN 1/4' }, $ ;        5
        { CW_PDMENU_S,       2, 'PAN 1/8' } $  ;      6

  ]


  PDMENU9 = CW_PDMENU( BASE6, MenuDesc1128, /RETURN_FULL_NAME, $
      UVALUE='PDMENU9', IDS=id_menu_pan)
      
      
  TextVal1130 = [ '' ]
  TEXT11 = WIDGET_TEXT( BASE7,VALUE=TextVal1130, $
      UVALUE='TEXT11', $
      YSIZE=1)

  MenuDesc1129 = [ $
      { CW_PDMENU_S,       2, 'RESOL UP' } $  ;      0

  ]
  
  PDMENU29 = CW_PDMENU( BASE7,  MenuDesc1129, /RETURN_FULL_NAME, $
      UVALUE='PDMENU29', IDS=id_menu_up)
      
  MenuDesc1130 = [ $
      { CW_PDMENU_S,       2, 'RESOL DOWN' } $  ;      0

  ]      
  
  PDMENU30 = CW_PDMENU( BASE7, MenuDesc1130, /RETURN_FULL_NAME, $
      UVALUE='PDMENU30', IDS=id_menu_down)      


 MenuDesc1131 = [ $
        { CW_PDMENU_S,       3, 'RESOL' }, $ ;        0
	{ CW_PDMENU_S,       0, 'RESOL Ima' }, $ ;        1
	{ CW_PDMENU_S,       0, 'RESOL 0' }, $ ;        2
	{ CW_PDMENU_S,       0, 'RESOL 1' }, $ ;        3
        { CW_PDMENU_S,       0, 'RESOL 2' }, $ ;        4
        { CW_PDMENU_S,       0, 'RESOL 3' }, $ ;        5
        { CW_PDMENU_S,       0, 'RESOL 4' }, $ ;        6
        { CW_PDMENU_S,       0, 'RESOL 5' }, $ ;        7
	{ CW_PDMENU_S,       0, 'RESOL 6' }, $ ;        8
        { CW_PDMENU_S,       0, 'RESOL 7' }, $ ;        9
        { CW_PDMENU_S,       0, 'RESOL 8' }, $ ;        10
        { CW_PDMENU_S,       0, 'RESOL 9' }, $ ;        11
        { CW_PDMENU_S,       2, 'RESOL 10' } $  ;       12
  ]
  PDMENU31 = CW_PDMENU(BASE7, MenuDesc1131, /RETURN_FULL_NAME, $
      UVALUE='PDMENU31', IDS=id_menu_resol)    
  
     
  DRAW13 = WIDGET_DRAW( MAIN1, $
      BUTTON_EVENTS=1, $
      MOTION_EVENTS=1, $
      RETAIN=2, $
      UVALUE='DRAW13', $
      XSIZE=1024, $
      X_SCROLL_SIZE=512, $
      YSIZE=1024, $
      Y_SCROLL_SIZE=512, $
      /VIEWPORT_EVENTS)

  WIDGET_CONTROL, MAIN1, /REALIZE
  
  ; Get drawable window index

  COMMON DRAW13_Comm, DRAW13_Id
  WIDGET_CONTROL, DRAW13, GET_VALUE=DRAW13_Id
  sensitive_init

  vsize = size(Imag)
  if vsize[0] EQ 2 then BEGIN
             load_ima
	     copy_ima
	     sensitive_ima
              END
  TabWin[11] = !d.window
  WIDGET_CONTROL, DRAW13, SET_DRAW_VIEW=[0,0]
  XMANAGER, 'MAIN1', MAIN1
  
  DONE:
   if N_PARAMS() EQ 1 and (size(Imag))[0] EQ 2 then Dat = Imag
END
