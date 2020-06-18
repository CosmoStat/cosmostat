;+
; NAME:
;	XMR_TRANSFORM
;
; PURPOSE:
;	Computes a multiresolution transform of an image. Parameters are
;       choosen through an widget interface. Creates a file which contains
;       the multiresolution transform of the image. The transform is made
;       by calling the routine mr_transform. If there is an image parameter
;       the user don't need to load an image with the interface.  
;
; CALLING SEQUENCE:
;
;       XMR_TRANSFORM, TransfData, input=input
;
; OUTPUT:
;       TransfData -- 2D or 3D array: multiresolution of an image
;
; INPUT KEYWORDS:
;	input -- 2D array: image which will be transformed	 
;
; RESTRICTIONS:
;	only  fits, midas and raw matrice image can be read
;
; PROCEDURE:
;       mr_transform
;
; MODIFICATION HISTORY:
; 	Written by:	Jean-Luc Starck
;	December, 1995  File creation
;------------------------------------------------------------------

pro common_pro
COMMON SHARE, TR_SLIDER1, TR_FIELD_NSCA, MainTrans, TR_FIELD_IM,  TR_FIELD_MR, $
              TR_BUTTON_CAN, TR_BUTTON_APP, $
              TR_FIND_SEL, TR_BASE0, TR_BASE1, TR_BASE2, TR_BASE3, TR_BASE4, $
              TypeSelect, NumSelect, NbrScale, NameImage, NameMR_File, $
              KEYIM, Dat, MResol
END


PRO  TR_MENU1_Event, Event
COMMON SHARE

  CASE Event.Value OF 
  'Select.linear wavelet transform: a trous algorithm': BEGIN
      TypeSelect='linear wavelet transform: a trous algorithm'
      NumSelect = 1
      WIDGET_CONTROL,  TR_FIND_SEL, SET_VALUE=TypeSelect
;     WIDGET_CONTROL,  TR_FIND_SEL, SENSITIVE=0 
    END
  'Select.bspline wavelet transform: a trous algorithm': BEGIN
     TypeSelect='bspline wavelet transform: a trous algorithm'
      NumSelect = 2 
      WIDGET_CONTROL,  TR_FIND_SEL, SET_VALUE=TypeSelect
     END
  'Select.wavelet transform in Fourier space': BEGIN
      TypeSelect='wavelet transform in Fourier space' 
      WIDGET_CONTROL,  TR_FIND_SEL, SET_VALUE=TypeSelect
      NumSelect = 3
    END
  'Select.morphological median transform': BEGIN
      TypeSelect='morphological median transform' 
      WIDGET_CONTROL,  TR_FIND_SEL, SET_VALUE=TypeSelect
      NumSelect = 4
    END
  'Select.morphological minmax transform': BEGIN
      TypeSelect='morphological minmax transform' 
      WIDGET_CONTROL,  TR_FIND_SEL, SET_VALUE=TypeSelect
      NumSelect = 5
    END
  'Select.pyramidal linear wavelet transform': BEGIN
      TypeSelect='pyramidal linear wavelet transform' 
      WIDGET_CONTROL,  TR_FIND_SEL, SET_VALUE=TypeSelect
      NumSelect = 6
    END
  'Select.pyramidal bspline wavelet transform': BEGIN
      TypeSelect='pyramidal bspline wavelet transform' 
      WIDGET_CONTROL,  TR_FIND_SEL, SET_VALUE=TypeSelect 
      NumSelect = 7
    END
  'Select.pyramidal wavelet transform in Fourier space: algo 1 (diff. between two resolutions)': BEGIN
      TypeSelect='pyramidal wavelet transform in Fourier space: algo 1 (diff. between two resolutions)' 
      WIDGET_CONTROL,  TR_FIND_SEL, SET_VALUE=TypeSelect
      NumSelect = 8
    END
  'Select.pyramidal wavelet transform in Fourier space: algo 2 (diff. between the square of two resolutions)': BEGIN
      TypeSelect='pyramidal wavelet transform in Fourier space: algo 2 (diff. between the square of two resolutions)' 
      WIDGET_CONTROL,  TR_FIND_SEL, SET_VALUE=TypeSelect
      NumSelect = 9
    END
  'Select.pyramidal median transform': BEGIN
      TypeSelect='pyramidal median transform' 
      WIDGET_CONTROL,  TR_FIND_SEL, SET_VALUE=TypeSelect
      NumSelect = 10
    END
  'Select.pyramidal laplacian': BEGIN
      TypeSelect='pyramidal laplacian' 
      WIDGET_CONTROL,  TR_FIND_SEL, SET_VALUE=TypeSelect
      NumSelect = 11
    END
  'Select.morphological pyramidal minmax transform': BEGIN
      TypeSelect='morphological pyramidal minmax transform' 
      WIDGET_CONTROL,  TR_FIND_SEL, SET_VALUE=TypeSelect
      NumSelect = 12
    END
  'Select.decomposition on scaling function': BEGIN
      TypeSelect='decomposition on scaling function' 
      WIDGET_CONTROL,  TR_FIND_SEL, SET_VALUE=TypeSelect
      NumSelect = 13
    END
  'Select.Mallat''s wavelet transform': BEGIN
      TypeSelect= 'Mallat''s wavelet transform'
      WIDGET_CONTROL,  TR_FIND_SEL, SET_VALUE=TypeSelect
      NumSelect = 14
    END
  'Select.Feauveau''s wavelet transform': BEGIN
      TypeSelect= 'Feauveau''s wavelet transform'
      WIDGET_CONTROL,  TR_FIND_SEL, SET_VALUE=TypeSelect
      NumSelect = 15
    END
  'Select.Feauveau''s wavelet transform without undersampling': BEGIN
      TypeSelect=  'Feauveau''s wavelet transform without undersampling'
      WIDGET_CONTROL,  TR_FIND_SEL, SET_VALUE=TypeSelect
      NumSelect = 16
    END
  'Select.G transform (morphological min-max algorithm)': BEGIN
      TypeSelect=  'G transform (morphological Haar algorithm)'
      WIDGET_CONTROL,  TR_FIND_SEL, SET_VALUE=TypeSelect
      NumSelect = 17
    END
  'Select.Haar''s wavelet transform': BEGIN
      TypeSelect= 'Haar''s wavelet transform'
      WIDGET_CONTROL,  TR_FIND_SEL, SET_VALUE=TypeSelect
      NumSelect = 18
     END
  'Select.Half-pyramidal transform': BEGIN
      TypeSelect=  'Half-pyramidal transform'
      WIDGET_CONTROL,  TR_FIND_SEL, SET_VALUE=TypeSelect
      NumSelect = 19
    END
  'Select.Mixed Half-pyramidal WT and Median method': BEGIN
      TypeSelect=  'Mixed Half-pyramidal WT and Median method'
      WIDGET_CONTROL,  TR_FIND_SEL, SET_VALUE=TypeSelect
      NumSelect = 20
    END
   'Select.diadic wavelet transform': BEGIN
      TypeSelect=  'diadic wavelet transform'
      WIDGET_CONTROL,  TR_FIND_SEL, SET_VALUE=TypeSelect
      NumSelect = 21
    END 
  'Select.Mixed WT and PMT method (WT-PMT)': BEGIN
      TypeSelect=  'Mixed WT and PMT method (WT-PMT)'
      WIDGET_CONTROL,  TR_FIND_SEL, SET_VALUE=TypeSelect
      NumSelect = 22
    END
   'Select.undecimated Haar transform: a trous algorithm': BEGIN
      TypeSelect=  'undecimated Haar transform: a trous algorithm'
      WIDGET_CONTROL,  TR_FIND_SEL, SET_VALUE=TypeSelect
      NumSelect = 23
    END 
   'Select.undecimated Mallat transform': BEGIN
      TypeSelect=  'undecimated Mallat transform'
      WIDGET_CONTROL,  TR_FIND_SEL, SET_VALUE=TypeSelect
      NumSelect = 24
    END 
   'Select.Wavelet transform via lifting scheme': BEGIN
      TypeSelect=  'Wavelet transform via lifting scheme'
      WIDGET_CONTROL,  TR_FIND_SEL, SET_VALUE=TypeSelect
      NumSelect = 25
    END 
  ENDCASE
END

;-----------------------------------------------------------

PRO MainTrans_Event, Event
COMMON SHARE

WIDGET_CONTROL,Event.Id,GET_UVALUE=Ev
  
Debug = 0
if Debug then print, 'Event = ', Ev

; event manadgement
  CASE Ev OF 
  'TR_FIND_SEL': $
      BEGIN
       ; Print, 'Event for Transform'
      END
  'TR_MENU1':  TR_MENU1_Event, Event
  'TR_FIELD_NSCA': $
      BEGIN
        ; read the field NbrScale
       WIDGET_CONTROL,  TR_FIELD_NSCA, GET_VALUE=NbrScale
        ; set the slider to the good value
       if NbrScale GT 1 and NbrScale LT 11 then $
                              WIDGET_CONTROL,  TR_SLIDER1, SET_VALUE=NbrScale
       if NbrScale GT 10  then $
         BEGIN
           WIDGET_CONTROL,  TR_FIELD_NSCA, SET_VALUE=4
           NbrScale = 4
         ENDIF
      END
  'TR_SLIDER1': $
      BEGIN
        ; read the slider and print it in the field ScaleNumber
        WIDGET_CONTROL,  TR_SLIDER1, GET_VALUE=NbrScale
        WIDGET_CONTROL, /HOUR
        WIDGET_CONTROL,  TR_FIELD_NSCA, SET_VALUE=NbrScale
      END
  'TR_FIELD_IM': $
      BEGIN
      ; read the image name
        WIDGET_CONTROL,  TR_FIELD_IM, GET_VALUE=XData
        NameImage = XData(0)
      END
  'TR_BUTTON_LOAD': $
      BEGIN
         ; select a file by using the widget pickfile
         ; and print the name in the field NameImage
        FileName = pickfile(/READ, FILTER='*.fits')
        if strlen(FileName) GT 0 then $
         BEGIN
           NameImage = FileName
           WIDGET_CONTROL,  TR_FIELD_IM, SET_VALUE=FileName
           KEYIM = 0
         END
      END
  'TR_FIELD_MR': $
      BEGIN
         ; read the multiresolution file name
       WIDGET_CONTROL,  TR_FIELD_MR, GET_VALUE=XData
       NameMR_File = XData(0)
      END
  'TR_BUTTON_CAN': WIDGET_CONTROL, MainTrans, /destroy
  'TR_BUTTON_APP': $
       BEGIN
        if strlen(NameImage) GT 1 then Format=format_image(NameImage) $
        else Format = 'unknown'    
          ; verifies all parameters
        if Debug EQ 1 then $
        BEGIN
           print, 'NbrScale = ', NbrScale
           print, 'Numselect = ', NumSelect, ': ', TypeSelect
           print, 'NameImage = ', NameImage
           print, 'Format = ', Format
           print, 'NameMR_File = ', NameMR_File
           print, 'KEYIM = ', KEYIM
        END

      Error = 0
      if KEYIM EQ 0 AND strlen(NameImage) LT 1 then $
        BEGIN
           Error = 1
             ; WIDGET_CONTROL, MainTrans, SENSITIVE=0 
           werr, 'Image File Name must be set'
        ENDIF
      if KEYIM EQ 0 AND Error EQ 0 AND Format EQ "unknown" then $
        BEGIN
           Error = 1
           werr, 'Unknown image format'
        ENDIF
      if Error EQ 0 AND strlen(NameMR_File) LT 1 then $
        BEGIN
           Error = 1
           werr, 'Multiresolution File Name must be set'
        ENDIF
      if Debug EQ 1 AND Error EQ 0 then print, 'OK '
      if Error EQ 0 and KEYIM EQ 0 then $
        BEGIN
           Dat = rim(NameImage)
           if (size(Dat))(0) NE 2 then  $
             BEGIN
                Error = 1
                werr, 'Unable to read image '+NameImage
             END
        ENDIF

       ; run the program 
       if (size(Dat))(0) EQ 2 then $ 
         BEGIN
           WIDGET_CONTROL, /HOUR
           Option = '-t ' + strcompress(string(Numselect))+ ' ' + $
                    '-n '+ strcompress(string(NbrScale))
           mr_transform, Dat, MResol, MR_File_Name=NameMR_File, OPT=Option 
           WIDGET_CONTROL, MainTrans, /destroy
         END
      END
  ENDCASE
END

;------------------------------------------------------------------------

PRO XMR_TRANSFORM, TransfData, input=input, GROUP=Group
COMMON SHARE

common_pro

NbrScale = 4
NameImage= ''
NameMR_File = 'xx_temp.mr'
NumSelect = 2
TypeSelect = 'b spline wavelet transform: a trous algorithm '
DataLoaded = 0

if N_PARAMS() LT 1 then begin 
        print, 'CALLING SEQUENCE: xmr_transform, TransfData, input=input'
        goto, DONE
        end

if keyword_set(input) then BEGIN
         KEYIM = 1 
         Dat = input
         DataLoaded = 1
         ENDIF else KEYIM = 0
  
  IF N_ELEMENTS(Group) EQ 0 THEN GROUP=0

if DataLoaded EQ 1 then BEGIN
   vsize = size(Dat)
   if vsize(0) NE 2 then begin
        print, 'ERROR: bad input parameter ...'
        print, 'CALLING SEQUENCE: xmr_transform, DataTransf, input=input'
        goto, DONE
        end
   end

  junk   = {  TR_MENUS, flags:0, name:' ' }

  MainTrans = WIDGET_BASE(GROUP_LEADER=Group, $
      TITLE = 'MULTIRESOLUTION TRANSFORM', $
      COLUMN=1, $
      MAP=1, $
      UVALUE='MainTrans', $
      XSIZE=700, $
      YSIZE=512)

   TR_BASE0 = WIDGET_BASE(MainTrans, $
      COLUMN=2, $
      XPAD=5, $
      YPAD=20, $
      MAP=1, $
      UVALUE='TR_BASE0')

  FieldVal123 = [ 'b spline wavelet transform: a trous algorithm ' ]
   TR_FIND_SEL = CW_FIELD(TR_BASE0,VALUE=FieldVal123, $
      ROW=1, $
      STRING=1, $
      TITLE='Transform', $
      UVALUE='TR_FIND_SEL', $
      XSIZE=55)

   TRMenuDesc1 = [ $
     {TR_MENUS, 3, 'Select' }, $ ; 0
     {TR_MENUS, 0, 'linear wavelet transform: a trous algorithm' },  $ ; 1
     {TR_MENUS, 0, 'bspline wavelet transform: a trous algorithm' }, $ ; 2
     {TR_MENUS, 0, 'wavelet transform in Fourier space' }, $           ; 3
     {TR_MENUS, 0, 'morphological median transform' }, $               ; 4
     {TR_MENUS, 0, 'morphological minmax transform' }, $               ; 5
     {TR_MENUS, 0, 'pyramidal linear wavelet transform' }, $           ; 6
     {TR_MENUS, 0, 'pyramidal bspline wavelet transform' }, $          ; 7
     {TR_MENUS, 0, 'pyramidal wavelet transform in Fourier space: algo 1 (diff. between two resolutions)' }, $                                      ; 8
     {TR_MENUS, 0, 'pyramidal wavelet transform in Fourier space: algo 2 (diff. between the square of two resolutions)' }, $                        ; 9
     {TR_MENUS, 0, 'pyramidal median transform' }, $                   ; 10
     {TR_MENUS, 0, 'pyramidal laplacian' }, $                          ; 11
     {TR_MENUS, 0, 'morphological pyramidal minmax transform' }, $     ; 12
     {TR_MENUS, 0, 'decomposition on scaling function' }, $            ; 13
     {TR_MENUS, 0, 'Mallat''s wavelet transform' }, $                  ; 14
     {TR_MENUS, 0, 'Feauveau''s wavelet transform' }, $                ; 15
     {TR_MENUS, 0, 'Feauveau''s wavelet transform without undersampling'}, $ ;16
     {TR_MENUS, 0, 'G transform (morphological min-max algorithm)' }, $; 17
     {TR_MENUS, 0, 'Haar''s wavelet transform' } , $                     ; 18
     {TR_MENUS, 0, 'Half-pyramidal transform' }, $                  ; 19
     {TR_MENUS, 0, 'Mixed Half-pyramidal WT and Median method' }, $ ;20                 ; 14
     {TR_MENUS, 0, 'diadic wavelet transform' }, $                  ; 21
     {TR_MENUS, 0, 'Mixed WT and PMT method (WT-PMT)' }, $                  ; 22
     {TR_MENUS, 0, 'undecimated Haar transform: a trous algorithm' }, $                  ; 23
     {TR_MENUS, 0, 'undecimated Mallat transform' }, $                  ; 24
     {TR_MENUS, 2, 'Wavelet transform via lifting scheme' } $                  ; 25
  ]

   TR_MENU1=CW_PDMENU(TR_BASE0,TRMenuDesc1,/RETURN_FULL_NAME,UVALUE='TR_MENU1')

   TR_BASE1 = WIDGET_BASE(MainTrans,  ROW=1,  SPACE=50,  XPAD=5, $
                          YPAD=20,  MAP=1,  UVALUE='TR_BASE1')

   FieldVal595 = ['4']
   TR_FIELD_NSCA = CW_FIELD (TR_BASE1,VALUE=FieldVal595,  ROW=1, $
                           INTEGER=1,  TITLE='Number of Scales', $
                           UVALUE='TR_FIELD_NSCA', $
                           XSIZE=2, /ALL_EVENTS)
   TR_SLIDER1 = WIDGET_SLIDER(TR_BASE1,  MAXIMUM=10,  MINIMUM=2, $
                              UVALUE='TR_SLIDER1',  VALUE=4)


   TR_BASE2 = WIDGET_BASE(MainTrans,  ROW=1,  XPAD=5,  YPAD=20, $
                          MAP=1,  UVALUE='TR_BASE2')
   if DataLoaded EQ 1 then FieldVal5435 = ['Image already loaded'] $
   else FieldVal5435 =  [NameImage]
   TR_FIELD_IM = CW_FIELD (TR_BASE2,VALUE=FieldVal5435,  ROW=1, $
                           STRING=1,  TITLE='Image File Name', $
                           UVALUE='TR_FIELD_IM', XSIZE=50, /ALL_EVENTS)
   TR_BUTTON_LOAD = WIDGET_BUTTON(TR_BASE2,  UVALUE='TR_BUTTON_LOAD', $
                                   VALUE='LOAD')


   TR_BASE3 = WIDGET_BASE(MainTrans,  ROW=1,  XPAD=5,  YPAD=2, $
                         MAP=1,  UVALUE='TR_BASE3')
   FieldVal5757 = ['xx_temp.mr']
   TR_FIELD_MR = CW_FIELD (TR_BASE3,VALUE=FieldVal5757, $
                           ROW=1,  STRING=1, $
                           TITLE='Output Multiresolution File Name (.mr)', $
                           UVALUE='TR_FIELD_MR', $
                           XSIZE=40, /ALL_EVENTS)


   TR_BASE4 = WIDGET_BASE(MainTrans,  ROW=1,  SPACE=100,  XPAD=150, $
                          YPAD=90,  MAP=1,  UVALUE='TR_BASE4')
   TR_BUTTON_CAN = WIDGET_BUTTON(TR_BASE4, FRAME=5, $
                                 UVALUE='TR_BUTTON_CAN', VALUE='CANCEL')
   TR_BUTTON_APP = WIDGET_BUTTON(TR_BASE4,  FRAME=5, $
                                  UVALUE='TR_BUTTON_APP', VALUE='APPLY')

  WIDGET_CONTROL, MainTrans, /REALIZE

  XMANAGER, 'MainTrans', MainTrans
  if NameMR_File EQ 'xx_temp.mr' then delete, NameMR_File
  TransfData = MResol

;---------------------------------------------------

DONE:

return
END
