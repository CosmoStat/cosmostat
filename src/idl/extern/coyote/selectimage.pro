;+
; NAME:
;   SELECTIMAGE
;
; PURPOSE:
;
;   The purpose of this program is to allow the user to select
;   an image file for reading. The image data is returned as the
;   result of the function. The best feature of this program is
;   the opportunity to browse the image before reading it.
;
; AUTHOR:
;
;   FANNING SOFTWARE CONSULTING
;   David Fanning, Ph.D.
;   1645 Sheely Drive
;   Fort Collins, CO 80526 USA
;   Phone: 970-221-0438
;   E-mail: davidf@dfanning.com
;   Coyote's Guide to IDL Programming: http://www.dfanning.com/
;
; CATEGORY:
;
;   General programming.
;
; CALLING SEQUENCE:
;
;   image = SelectImage()
;
; INPUT PARAMETERS:
;
;   None. All input is via keywords.
;
; INPUT KEYWORDS:
;
;   BMP -- Set this keyword to select BMP files.
;
;   DICOM -- Set this keyword to select DICOM files.
;
;   DIRECTORY -- The initial input directory name. The current directory by default.
;
;   FILENAME -- The initial filename. If the initial directory has image files of the
;               correct type, the default is to display the first of these files. Otherwise, blank.
;
;   _EXTRA -- This keyword is used to collect and pass keywords on to the FSC_FILESELECT object. See
;             the code for FSC_FILESELECT for details.
;   GIF -- Set this keyword to select GIF files. This capability is not available in IDL 5.4 and higher.
;
;   GROUP_LEADER -- Set this keyword to a widget identifier group leader. This keyword MUST be
;                   set when calling this program from another widget program to guarantee modal operation.
;
;   JPEG -- Set this keyword to select JPEG files.
;
;   PICT -- Set this keyword to select PICT files.
;
;   PNG -- Set this keyword to select PNG files.
;
;   PREVIEWSIZE -- Set this keyword to the maximum size (in pixels) of the preview window. Default is 150.
;
;   TIFF -- Set this keyword to select TIFF files. (This is the default filter selection.)
;
; OUTPUT KEYWORDS:
;
;   CANCEL -- This keyword is set to 1 if the user exits the program in any way except hitting the ACCEPT button.
;             The ACCEPT button will set this keyword to 0.
;
;   FILEINFO -- This keyword returns information about the selected file. Obtained from the QUERY_**** functions.
;
;   OUTDIRECTORY -- The directory where the selected file is found.
;
;   OUTFILENAME -- The short filename of the selected file.
;
; COMMON BLOCKS:
;
;   None.
;
; RESTRICTIONS:
;
;   Probably doesn't work correctly on VMS systems :-( If you can help, please
;   contact me. I don't have a VMS system to test on.
;
; OTHER COYOTE LIBRARY FILES REQUIRED:
;
;  http://www.dfanning.com/programs/fsc_fileselect.pro
;  http://www.dfanning.com/programs/tvimage.pro
;
; EXAMPLE:
;
;   To read JPEG files from the directory:
;
;      IDL> image = SelectImage(/JPEG)
;
; MODIFICATION HISTORY:
;
;   Written by: David Fanning, 18 Jan 2001.
;   Added modification to read both 8-bit and 24-bit BMP files. 27 Jan 2001. DWF.
;   Added FlipImage keyword and changed the action of the FLIP IMAGE button
;     to actually reverse the Y order instead of just change the !Order setting. 16 March 2001. DWF.
;   This version modified slightly to run PETROGRAPH_CONFIG for SHELL. 16 March 2001. DWF.
;   Added OutDirectory and OutFilename output keywords. 11 June 2001. DWF.
;-
;
;###########################################################################
;
; LICENSE
;
; This software is OSI Certified Open Source Software.
; OSI Certified is a certification mark of the Open Source Initiative.
;
; Copyright © 2000 Fanning Software Consulting
;
; This software is provided "as-is", without any express or
; implied warranty. In no event will the authors be held liable
; for any damages arising from the use of this software.
;
; Permission is granted to anyone to use this software for any
; purpose, including commercial applications, and to alter it and
; redistribute it freely, subject to the following restrictions:
;
; 1. The origin of this software must not be misrepresented; you must
;    not claim you wrote the original software. If you use this software
;    in a product, an acknowledgment in the product documentation
;    would be appreciated, but is not required.
;
; 2. Altered source versions must be plainly marked as such, and must
;    not be misrepresented as being the original software.
;
; 3. This notice may not be removed or altered from any source distribution.
;
; For more information on Open Source Software, visit the Open Source
; web site: http://www.opensource.org.
;
;###########################################################################


FUNCTION SelectImage_Dimensions, image, $

; This function returns the dimensions of the image, and also
; extracts relevant information via output keywords. Works only
; with 2D and 3D (24-bit) images.

   XSize=xsize, $          ; Output keyword. The X size of the image.
   YSize=ysize, $          ; Output keyword. The Y size of the image.
   TrueIndex=trueindex, $  ; Output keyword. The position of the "true color" index. -1 for 2D images.
   XIndex=xindex, $        ; Output keyword. The position or index of the X image size.
   YIndex=yindex           ; Output keyword. The position or index of the Y image size.

   ; Get the number of dimensions and the size of those dimensions.

ndims = Size(image, /N_Dimensions)
dims =  Size(image, /Dimensions)

   ; Is this a 2D or 3D image?

IF ndims EQ 2 THEN BEGIN
   xsize = dims[0]
   ysize = dims[1]
   trueindex = -1
   xindex = 0
   yindex = 1
ENDIF ELSE BEGIN
   IF ndims NE 3 THEN Message, /NoName, 'Unknown image dimensions. Returning.'
   true = Where(dims EQ 3, count)
   trueindex = true[0]
   IF count EQ 0 THEN Message, /NoName, 'Unknown image type. Returning.'
   CASE true[0] OF
      0: BEGIN
         xsize = dims[1]
         ysize = dims[2]
         xindex = 1
         yindex = 2
         ENDCASE
      1: BEGIN
         xsize = dims[0]
         ysize = dims[2]
         xindex = 0
         yindex = 2
         ENDCASE
      2: BEGIN
         xsize = dims[0]
         ysize = dims[1]
         xindex = 0
         yindex = 1
         ENDCASE
   ENDCASE
ENDELSE
RETURN, dims
END; ----------------------------------------------------------------------------------------


PRO SelectImage_FlipImage, event

; This event handler reverses the Y dimension of the image and re-displays it.

Widget_Control, event.top, Get_UValue=info, /No_Copy

dims = SelectImage_Dimensions(*(*(info.storagePtr)).image, YIndex=yindex)
*(*(info.storagePtr)).image = Reverse(*(*(info.storagePtr)).image, yindex + 1)
WSet, info.previewWID
TVLCT, info.r, info.g, info.b
IF (Min(*(*(info.storagePtr)).image) LT 0) OR (Max(*(*(info.storagePtr)).image) GT (!D.Table_Size-1)) THEN $
   TVImage, BytScl(*(*(info.storagePtr)).image, Top=!D.Table_Size-1), /Keep_Aspect, /NoInterpolation, /Erase ELSE $
   TVImage, *(*(info.storagePtr)).image, /Keep_Aspect, /NoInterpolation, /Erase

Widget_Control, event.top, Set_UValue=info, /No_Copy
END; ----------------------------------------------------------------------------------------



PRO SelectImage_SetFilter, event

; This event handler sets the filter for image data files.

Widget_Control, event.top, Get_UValue=info, /No_Copy

   ; The filter is in the User Value of the button. Store it.

Widget_Control, event.id, Get_UValue=theFilter
info.filter = theFilter

   ; Get the current filename.

Widget_Control, info.filenameID, Get_Value=filename

   ; Set the new filter in the Filename compound widget.

info.filenameObj->SetProperty, Filter=theFilter

   ; Look in the data directory for the files.

CD, info.dataDirectory, Current=thisDirectory

   ; Did you find any files?

theFiles = Findfile(info.filter, Count=fileCount)
IF fileCount EQ 0 THEN BEGIN
   theFiles = ""
   filename = theFilter
ENDIF ELSE BEGIN
   IF (filename EQ "") THEN BEGIN
      filename = theFiles[0]
   ENDIF ELSE BEGIN
      filename = theFilter
   ENDELSE
ENDELSE

   ; Update the widget interface according to what you found.

Widget_Control, info.filenameID, Set_Value=filename
Widget_Control, info.fileListID, Set_Value=theFiles
*info.theFiles = theFiles

   ; Is this a valid image file name. If so, go get the image.

image = BytArr(info.previewsize, info.previewsize)
fileInfo = {channels:0, dimensions:[info.previewsize, info.previewsize]}

IF filename NE "" THEN BEGIN

   CASE info.filter OF

      "*.bmp": BEGIN
         ok = Query_BMP(filename, fileInfo)
         IF ok THEN IF fileInfo.channels EQ 3 THEN image = Read_BMP(filename, /RGB) ELSE $
                                                   image = Read_BMP(filename, r, g, b)
         ENDCASE

      "*.dcm": BEGIN
         ok = Query_DICOM(filename, fileInfo)
         IF ok THEN image = Read_Dicom(filename, r, g, b)
         ENDCASE

      "*.gif": BEGIN
         ok = Query_GIF(filename, fileInfo)
         IF ok THEN Read_GIF, filename, image, r, g, b
         ENDCASE

      "*.pict": BEGIN
         ok = Query_PICT(filename, fileInfo)
         IF ok THEN Read_PICT, filename, image, r, g, b
         ENDCASE

      "*.png": BEGIN
         ok = Query_PNG(filename, fileInfo)
         IF ok THEN image = Read_PNG(filename, r, g, b)
         ENDCASE

      "*.jpg": BEGIN
         ok = Query_JPEG(filename, fileInfo)
         IF ok THEN Read_JPEG, filename, image, True=1
         ENDCASE

      "*.jpeg": BEGIN
         ok = Query_JPEG(filename, fileInfo)
         IF ok THEN Read_JPEG, filename, image, True=1
         ENDCASE

      "*.tif": BEGIN
         ok = Query_TIFF(filename, fileInfo)
         IF ok THEN BEGIN
            CASE fileInfo.channels OF
               3: image = Read_TIFF(filename)
               ELSE: image = Read_TIFF(filename, r, g, b)
            ENDCASE
         ENDIF
         ENDCASE

      "*.tiff": BEGIN
         ok = Query_TIFF(filename, fileInfo)
         IF ok THEN BEGIN
            CASE fileInfo.channels OF
               3: image = Read_TIFF(filename)
               ELSE: image = Read_TIFF(filename, r, g, b)
            ENDCASE
         ENDIF
         ENDCASE

      ELSE: BEGIN
         Message, 'File type unrecognized', /Infomational
         ENDCASE

   ENDCASE

ENDIF

   ; Store RGB vectors if they got set.

IF N_Elements(r) NE 0 THEN info.r = r ELSE info.r = Bindgen(!D.Table_Size)
IF N_Elements(g) NE 0 THEN info.g = g ELSE info.g = Bindgen(!D.Table_Size)
IF N_Elements(b) NE 0 THEN info.b = b ELSE info.b = Bindgen(!D.Table_Size)

   ; What kind of image is this?

CASE fileinfo.channels OF
   3: imageType = "True-Color Image"
   0: imageType = "No Image"
   ELSE: imageType = '8-Bit Image'
ENDCASE

   ; Get the file sizes.

dimensions = SelectImage_Dimensions(image, XSize=xsize, YSize=ysize, YIndex=yindex)

   ; Flip the image if required.

IF info.flipimage THEN image = Reverse(image, yindex+1)

   ; Calculate a window size for the image preview.

aspect = Float(xsize) / ysize
IF aspect GT 1 THEN BEGIN
   wxsize = Fix(info.previewSize)
   wysize = Fix(info.previewSize / aspect)
ENDIF ELSE BEGIN
   wysize = Fix(info.previewSize)
   wxsize = Fix(info.previewSize / aspect)
ENDELSE

   ; If you don't have an image, then get sensible numbers for the labels.

IF imageType EQ 'No Image' THEN BEGIN
   xsize = 0
   ysize = 0
   minval = 0
   maxval = 0
ENDIF

   ; Update the display with what you have.

IF imageType EQ 'No Image' THEN imageDataType = 'NONE' ELSE imageDataType = Size(image, /TNAME)
Widget_Control, info.labelTypeID, Set_Value=imageType
Widget_Control, info.labelXSizeID, Set_Value="X Size: " + StrTrim(xsize, 2)
Widget_Control, info.labelYSizeID, Set_Value="Y Size: " + StrTrim(ysize, 2)
Widget_Control, info.labelDataTypeID, Set_Value="Type: " + imageDataType
Widget_Control, info.labelminvalID, Set_Value="Min Value: " + StrTrim(Fix(Min(image)), 2)
Widget_Control, info.labelmaxvalID, Set_Value="Max Value: " + StrTrim(Fix(Max(image)), 2)
Widget_Control, info.previewID, Draw_XSize=wxsize, Draw_YSize=wysize

   ; Draw the preview image.

WSet, info.previewWID
TVLCT, info.r, info.g, info.b
IF (Min(image) LT 0) OR (Max(image) GT (!D.Table_Size-1)) THEN $
   TVImage, BytScl(image, Top=!D.Table_Size-1), /Keep_Aspect, /NoInterpolation, /Erase ELSE $
   TVImage, image, /Keep_Aspect, /NoInterpolation, /Erase
IF imageDataType EQ 'NONE' THEN image = 0

   ; Save the image data for later retrieval.

*(*(info.storagePtr)).image = image
*(*(info.storagePtr)).fileInfo = fileInfo

   ; Clean up.

CD, thisDirectory
Widget_Control, event.top, Set_UValue=info, /No_Copy
END; ----------------------------------------------------------------------------------------



PRO SelectImage_FilenameEvents, event

Widget_Control, event.top, Get_UValue=info, /No_Copy

   ; Get the name of the file.

filename = event.basename

   ; Set the filter name, if required.

IF !Version.Release LT 5.3 THEN dot = RStrPos(filename, ".") ELSE $
   dot = StrPos(filename, ".", /Reverse_Search)
thisfilter = '*' + StrLowCase(StrMid(filename, dot))

CD, event.directory, Current=thisDirectory
theFiles = Findfile(thisfilter, Count=fileCount)
IF fileCount EQ 0 THEN theFiles = "" ELSE BEGIN
   IF (filename EQ "" OR filename EQ thisfilter) AND (thisDirectory NE info.dataDirectory) THEN filename = theFiles[0]
ENDELSE
info.dataDirectory = event.directory

Widget_Control, info.fileListID, Set_Value=theFiles
*info.theFiles = theFiles

   ; Is this a valid image file name. If so, go get the image.

image = BytArr(info.previewsize, info.previewsize)
fileInfo = {channels:0, dimensions:[info.previewsize, info.previewsize]}

IF filename NE "" THEN BEGIN

   CASE thisfilter OF

      "*.bmp": BEGIN
         ok = Query_BMP(filename, fileInfo)
         IF ok THEN IF fileInfo.channels EQ 3 THEN image = Read_BMP(filename, /RGB) ELSE $
                                                   image = Read_BMP(filename, r, g, b)
         IF thisfilter NE info.filter THEN BEGIN
            info.filter = thisfilter
            info.filenameObj->SetProperty, Filter=thisFilter
         ENDIF
         ENDCASE

      "*.dcm": BEGIN
         ok = Query_DICOM(filename, fileInfo)
         IF ok THEN image = Read_Dicom(filename, r, g, b)
         IF thisfilter NE info.filter THEN BEGIN
            info.filter = thisfilter
            info.filenameObj->SetProperty, Filter=thisFilter
         ENDIF
         ENDCASE

      "*.gif": BEGIN
         ok = Query_GIF(filename, fileInfo)
         IF ok THEN Read_GIF, filename, image, r, g, b
         IF thisfilter NE info.filter THEN BEGIN
            info.filter = thisfilter
            info.filenameObj->SetProperty, Filter=thisFilter
         ENDIF
         ENDCASE

      "*.pict": BEGIN
         ok = Query_PICT(filename, fileInfo)
         IF ok THEN Read_PICT, filename, image, r, g, b
         IF thisfilter NE info.filter THEN BEGIN
            info.filter = thisfilter
            info.filenameObj->SetProperty, Filter=thisFilter
         ENDIF
         ENDCASE

      "*.png": BEGIN
         ok = Query_PNG(filename, fileInfo)
         IF ok THEN image = Read_PNG(filename, r, g, b)
         IF thisfilter NE info.filter THEN BEGIN
            info.filter = thisfilter
            info.filenameObj->SetProperty, Filter=thisFilter
         ENDIF
         ENDCASE

      "*.jpg": BEGIN
         ok = Query_JPEG(filename, fileInfo)
         IF ok THEN Read_JPEG, filename, image, True=1
         IF thisfilter NE info.filter THEN BEGIN
            info.filter = thisfilter
            info.filenameObj->SetProperty, Filter=thisFilter
         ENDIF
         ENDCASE

      "*.jpeg": BEGIN
         ok = Query_JPEG(filename, fileInfo)
         IF ok THEN Read_JPEG, filename, image, True=1
         IF thisfilter NE info.filter THEN BEGIN
            info.filter = thisfilter
            info.filenameObj->SetProperty, Filter=thisFilter
         ENDIF
         ENDCASE

      "*.tif": BEGIN
         ok = Query_TIFF(filename, fileInfo)
         IF ok THEN BEGIN
            CASE fileInfo.channels OF
               3: image = Read_TIFF(filename)
               ELSE: image = Read_TIFF(filename, r, g, b)
            ENDCASE
         ENDIF
         IF thisfilter NE info.filter THEN BEGIN
            info.filter = thisfilter
            info.filenameObj->SetProperty, Filter=thisFilter
         ENDIF
         ENDCASE

      "*.tiff": BEGIN
         ok = Query_TIFF(filename, fileInfo)
         IF ok THEN BEGIN
            CASE fileInfo.channels OF
               3: image = Read_TIFF(filename)
               ELSE: image = Read_TIFF(filename, r, g, b)
            ENDCASE
         ENDIF
         IF thisfilter NE info.filter THEN BEGIN
            info.filter = thisfilter
            info.filenameObj->SetProperty, Filter=thisFilter
         ENDIF
         ENDCASE

      ELSE:

   ENDCASE

ENDIF

   ; Store RGB vectors if they got set.

IF N_Elements(r) NE 0 THEN info.r = r ELSE info.r = Bindgen(!D.Table_Size)
IF N_Elements(g) NE 0 THEN info.g = g ELSE info.g = Bindgen(!D.Table_Size)
IF N_Elements(b) NE 0 THEN info.b = b ELSE info.b = Bindgen(!D.Table_Size)

   ; What kind of image is this?

CASE fileinfo.channels OF
   3: imageType = "True-Color Image"
   0: imageType = "No Image"
   ELSE: imageType = '8-Bit Image'
ENDCASE

   ; Get the file sizes.

dimensions = SelectImage_Dimensions(image, XSize=xsize, YSize=ysize, YIndex=yindex)

   ; Flip the image if required.

IF info.flipimage THEN image = Reverse(image, yindex+1)

   ; Calculate a window size for the image preview.

aspect = Float(xsize) / ysize
IF aspect GT 1 THEN BEGIN
   wxsize = Fix(info.previewSize)
   wysize = Fix(info.previewSize / aspect)
ENDIF ELSE BEGIN
   wysize = Fix(info.previewSize)
   wxsize = Fix(info.previewSize / aspect)
ENDELSE

   ; If you don't have an image, then get sensible numbers for the labels.

IF imageType EQ 'No Image' THEN BEGIN
   xsize = 0
   ysize = 0
   minval = 0
   maxval = 0
ENDIF

   ; Update the display with what you have.

IF imageType EQ 'No Image' THEN imageDataType = 'NONE' ELSE imageDataType = Size(image, /TNAME)
Widget_Control, info.labelTypeID, Set_Value=imageType
Widget_Control, info.labelXSizeID, Set_Value="X Size: " + StrTrim(xsize, 2)
Widget_Control, info.labelYSizeID, Set_Value="Y Size: " + StrTrim(ysize, 2)
Widget_Control, info.labelDataTypeID, Set_Value="Type: " + imageDataType
Widget_Control, info.labelminvalID, Set_Value="Min Value: " + StrTrim(Fix(Min(image)), 2)
Widget_Control, info.labelmaxvalID, Set_Value="Max Value: " + StrTrim(Fix(Max(image)), 2)
Widget_Control, info.previewID, Draw_XSize=wxsize, Draw_YSize=wysize

   ; Draw the preview image.

WSet, info.previewWID
TVLCT, info.r, info.g, info.b
IF (Min(image) LT 0) OR (Max(image) GT (!D.Table_Size-1)) THEN $
   TVImage, BytScl(image, Top=!D.Table_Size-1), /Keep_Aspect, /NoInterpolation, /Erase ELSE $
   TVImage, image, /Keep_Aspect, /NoInterpolation, /Erase
IF imageDataType EQ 'NONE' THEN image = 0

   ; Store the image data for later retrieval.

*(*(info.storagePtr)).image = image
*(*(info.storagePtr)).fileInfo = fileInfo

   ; Clean up.

CD, thisDirectory
Widget_Control, event.top, Set_UValue=info, /No_Copy
END ;---------------------------------------------------------------------------------



PRO SelectImage_ListEvents, event

   ; Only handle single click events.

IF event.clicks NE 1 THEN RETURN

Widget_Control, event.top, Get_UValue=info, /No_Copy

   ; Get the name of the file.

filename = (*info.theFiles)[event.index]
CD, info.dataDirectory, Current=thisDirectory

   ; Set it in the Filename widget.

Widget_Control, info.filenameID, Set_Value=filename

   ; Is this a valid image file name. If so, go get the image.

image = BytArr(info.previewsize, info.previewsize)
fileInfo = {channels:2, dimensions:[info.previewsize, info.previewsize]}

IF filename NE "" THEN BEGIN

   CASE info.filter OF

      "*.bmp": BEGIN
         ok = Query_BMP(filename, fileInfo)
         IF ok THEN IF fileInfo.channels EQ 3 THEN image = Read_BMP(filename, /RGB) ELSE $
                                                   image = Read_BMP(filename, r, g, b)
         ENDCASE

      "*.dcm": BEGIN
         ok = Query_DICOM(filename, fileInfo)
         IF ok THEN image = Read_Dicom(filename, r, g, b)
         ENDCASE

      "*.gif": BEGIN
         ok = Query_GIF(filename, fileInfo)
         IF ok THEN Read_GIF, filename, image, r, g, b
         ENDCASE

      "*.pict": BEGIN
         ok = Query_PICT(filename, fileInfo)
         IF ok THEN Read_PICT, filename, image, r, g, b
         ENDCASE

      "*.png": BEGIN
         ok = Query_PNG(filename, fileInfo)
         IF ok THEN image = Read_PNG(filename, r, g, b)
         ENDCASE

      "*.jpg": BEGIN
         ok = Query_JPEG(filename, fileInfo)
         IF ok THEN Read_JPEG, filename, image, True=1
         ENDCASE

      "*.jpeg": BEGIN
         ok = Query_JPEG(filename, fileInfo)
         IF ok THEN Read_JPEG, filename, image, True=1
         ENDCASE

      "*.tif": BEGIN
         ok = Query_TIFF(filename, fileInfo)
         IF ok THEN BEGIN
            CASE fileInfo.channels OF
               3: image = Read_TIFF(filename)
               ELSE: image = Read_TIFF(filename, r, g, b)
            ENDCASE
         ENDIF
         ENDCASE

      "*.tiff": BEGIN
         ok = Query_TIFF(filename, fileInfo)
         IF ok THEN BEGIN
            CASE fileInfo.channels OF
               3: image = Read_TIFF(filename)
               ELSE: image = Read_TIFF(filename, r, g, b)
            ENDCASE
         ENDIF
         ENDCASE

      ELSE: BEGIN
         Message, 'File type unrecognized', /Infomational
         ENDCASE

   ENDCASE

ENDIF

   ; Store RGB vectors if they got set.

IF N_Elements(r) NE 0 THEN info.r = r ELSE info.r = Bindgen(!D.Table_Size)
IF N_Elements(g) NE 0 THEN info.g = g ELSE info.g = Bindgen(!D.Table_Size)
IF N_Elements(b) NE 0 THEN info.b = b ELSE info.b = Bindgen(!D.Table_Size)

   ; What kind of image is this?

CASE fileinfo.channels OF
   3: imageType = "True-Color Image"
   0: imageType = "No Image"
   ELSE: imageType = '8-Bit Image'
ENDCASE

   ; Get the file sizes.

dimensions = SelectImage_Dimensions(image, XSize=xsize, YSize=ysize, YIndex=yindex)

   ; Flip the image if required.

IF info.flipimage THEN image = Reverse(image, yindex+1)

   ; Calculate a window size for the image preview.

aspect = Float(xsize) / ysize
IF aspect GT 1 THEN BEGIN
   wxsize = Fix(info.previewSize)
   wysize = Fix(info.previewSize / aspect)
ENDIF ELSE BEGIN
   wysize = Fix(info.previewSize)
   wxsize = Fix(info.previewSize / aspect)
ENDELSE

   ; If you don't have an image, then get sensible numbers for the labels.

IF imageType EQ 'No Image' THEN BEGIN
   xsize = 0
   ysize = 0
   minval = 0
   maxval = 0
ENDIF

   ; Update the display with what you have.

IF imageType EQ 'No Image' THEN imageDataType = 'NONE' ELSE imageDataType = Size(image, /TNAME)
Widget_Control, info.labelTypeID, Set_Value=imageType
Widget_Control, info.labelXSizeID, Set_Value="X Size: " + StrTrim(xsize, 2)
Widget_Control, info.labelYSizeID, Set_Value="Y Size: " + StrTrim(ysize, 2)
Widget_Control, info.labelDataTypeID, Set_Value="Type: " + imageDataType
Widget_Control, info.labelminvalID, Set_Value="Min Value: " + StrTrim(Fix(Min(image)), 2)
Widget_Control, info.labelmaxvalID, Set_Value="Max Value: " + StrTrim(Fix(Max(image)), 2)
Widget_Control, info.previewID, Draw_XSize=wxsize, Draw_YSize=wysize

   ; Draw the preview image.

WSet, info.previewWID
TVLCT, info.r, info.g, info.b
IF (Min(image) LT 0) OR (Max(image) GT (!D.Table_Size-1)) THEN $
   TVImage, BytScl(image, Top=!D.Table_Size-1), /Keep_Aspect, /NoInterpolation, /Erase ELSE $
   TVImage, image, /Keep_Aspect, /NoInterpolation, /Erase
IF imageDataType EQ 'NONE' THEN image = 0

   ; Store the image data for later retrieval.

*(*(info.storagePtr)).image = image
*(*(info.storagePtr)).fileInfo = fileInfo

   ; Clean up.

CD, thisDirectory
Widget_Control, event.top, Set_UValue=info, /No_Copy
END ;---------------------------------------------------------------------------------



PRO SelectImage_Action, event

; This event handler responds to CANCEL and ACCEPT buttons.

Widget_Control, event.top, Get_UValue=info, /No_Copy
Widget_Control, event.id, Get_Value=buttonValue
IF buttonValue EQ 'Accept' THEN (*info.storagePtr).cancel = 0
info.filenameObj->GetProperty, Directory=outdirectory, Filename=outfilename
(*info.storagePtr).outdirectory = outdirectory
(*info.storagePtr).outfilename = outfilename
Widget_Control, event.top, Set_UValue=info, /No_Copy
Widget_Control, event.top, /Destroy
END ;---------------------------------------------------------------------------------



PRO SelectImage_Cleanup, tlb

; Program pointers are cleaned up here.

Widget_Control, tlb, Get_UValue=info, /No_Copy
IF N_Elements(info) EQ 0 THEN RETURN
Ptr_Free, info.theFiles
END ;---------------------------------------------------------------------------------



FUNCTION SelectImage, $
   BMP=bmp, $                      ; Set this keyword to select BMP files.
   Cancel=cancel, $                ; An output keyword. Returns 0 if the ACCEPT button is used, 1 otherwise.
   Dicom=dicom, $                  ; Set this keyword to select DICOM files
   Directory=directory, $          ; Initial directory to search for files.
   FileInfo=fileInfo, $            ; An output keyword containing file information from the Query_*** routine.
   Filename=filename, $            ; Initial file name of image file.
   Flipimage=flipimage, $          ; Set this keyword to flip the Y indices of the image. Set to 1 by default.
   _Extra=extra, $                 ; This is used to pass keywords on to FSC_FILESELECT. See that documentation for details.
   GIF=gif, $                      ; Set this keyword to select GIF files
   Group_Leader=group_leader, $    ; The group leader ID of this widget program.
   JPEG=jpeg, $                    ; Set this keyword to select JPEG files
   OutDirectory=outdirectory, $    ; The directory name of the selected image file.
   OutFilename=outfilename, $      ; The short filename (without directory) of the selected image file.
   PICT=pict, $                    ; Set this keyword to select PICT files
   PNG=png, $                      ; Set this keyword to select PNG files
   TIFF=tiff, $                    ; Set this keyword to select TIFF files
   PreviewSize=previewsize         ; The maximum size of the image preview window. 150 pixels by default.


   ; Error handling.

Catch, theError
IF theError NE 0 THEN BEGIN
   Catch, /Cancel
   Cancel = 1
   ok =Dialog_Message(!Error_State.Msg + ' Returning...')
   RETURN, 0
ENDIF

   ; Check for availability of GIF files.

thisVersion = Float(!Version.Release)
IF thisVersion LT 5.3 THEN haveGif = 1 ELSE haveGIF = 0

   ; Set up the filter.

filter = "*.tif"
IF Keyword_Set(bmp) THEN filter = "*.bmp"
IF Keyword_Set(dicom) THEN filter = "*.dcm"
IF Keyword_Set(gif) THEN BEGIN
   IF havegif THEN filter = "*.gif" ELSE $
      ok = Dialog_Message('GIF files not supported in this IDL version. Replacing with TIFF.')
ENDIF
IF Keyword_Set(pict) THEN filter = "*.pict"
IF Keyword_Set(png) THEN filter = "*.png"
IF Keyword_Set(jpeg) THEN filter = "*.jpg"
IF Keyword_Set(tiff) THEN filter = "*.tif"
IF N_Elements(flipimage) EQ 0 THEN flipimage = 1 ELSE flipimage = Keyword_Set(flipimage)

   ; Get the current directory. Some processing involved.

CD, Current=startDirectory

   ; Run configuration program to establish working directories.

Petrograph_Config

IF N_Elements(directory) EQ 0 THEN directory = !PETROGRAPH_DATA_DIRECTORY ELSE BEGIN
   IF StrMid(directory, 0, 2) EQ ".." THEN BEGIN
      CASE StrUpCase(!Version.OS_Family) OF
      'MACOS': BEGIN
         CD, '..'
         CD, Current=basename
         directory = basename + StrMid(directory, 3)
         END
      'VMS': BEGIN
         CD, '..'
         CD, Current=basename
         directory = basename + StrMid(directory, 3)
         END
      ELSE: BEGIN
         CD, '..'
         CD, Current=basename
         directory = basename + StrMid(directory, 2)
         END
      ENDCASE
   ENDIF
   IF StrMid(directory, 0, 1) EQ "." THEN BEGIN
      CASE StrUpCase(!Version.OS_Family) OF
      'MACOS': BEGIN
         CD, Current=basename
         directory = basename + StrMid(directory, 2)
         END
      'VMS': BEGIN
         CD, Current=basename
         directory = basename + StrMid(directory, 2)
         END
      ELSE: BEGIN
         CD, Current=basename
         directory = basename + StrMid(directory, 1)
      END
      ENDCASE
   ENDIF
ENDELSE
CD, directory

   ; Check other keyword values.

IF N_Elements(filename) EQ 0 THEN filename = ""
IF N_Elements(previewSize) EQ 0 THEN previewSize = 150

   ; Locate appropriate files.

theFiles = Findfile(filter, Count=fileCount)
IF fileCount EQ 0 THEN theFiles = "" ELSE BEGIN
   IF filename EQ "" THEN filename = theFiles[0]
ENDELSE

   ; Is this a valid image file name. If so, go get the image.

image = BytArr(previewsize, previewsize)
fileInfo = {channels:2, dimensions:[previewsize, previewsize]}

IF filename NE "" THEN BEGIN

   CASE filter OF

      "*.bmp": BEGIN
         ok = Query_BMP(filename, fileInfo)
         IF ok THEN IF fileInfo.channels EQ 3 THEN image = Read_BMP(filename, /RGB) ELSE $
                                                   image = Read_BMP(filename, r, g, b)
         ENDCASE

      "*.dcm": BEGIN
         ok = Query_DICOM(filename, fileInfo)
         IF ok THEN image = Read_Dicom(filename, r, g, b)
         ENDCASE

      "*.gif": BEGIN
         ok = Query_GIF(filename, fileInfo)
         IF ok THEN Read_GIF, filename, image, r, g, b
         ENDCASE

      "*.pict": BEGIN
         ok = Query_PICT(filename, fileInfo)
         IF ok THEN Read_PICT, filename, image, r, g, b
         ENDCASE

      "*.png": BEGIN
         ok = Query_PNG(filename, fileInfo)
         IF ok THEN image = Read_PNG(filename, r, g, b)
         ENDCASE

      "*.jpg": BEGIN
         ok = Query_JPEG(filename, fileInfo)
         IF ok THEN Read_JPEG, filename, image, True=1
         ENDCASE

      "*.jpeg": BEGIN
         ok = Query_JPEG(filename, fileInfo)
         IF ok THEN Read_JPEG, filename, image, True=1
         ENDCASE

      "*.tif": BEGIN
         ok = Query_TIFF(filename, fileInfo)
         IF ok THEN BEGIN
            CASE fileInfo.channels OF
               3: image = Read_TIFF(filename)
               ELSE: image = Read_TIFF(filename, r, g, b)
            ENDCASE
         ENDIF
         ENDCASE

      "*.tiff": BEGIN
         ok = Query_TIFF(filename, fileInfo)
         IF ok THEN BEGIN
            CASE fileInfo.channels OF
               3: image = Read_TIFF(filename)
               ELSE: image = Read_TIFF(filename, r, g, b)
            ENDCASE
         ENDIF
         ENDCASE

      ELSE: BEGIN
         Message, 'File type unrecognized', /Infomational
         ENDCASE

   ENDCASE
ENDIF

   ; Get the file sizes.

dimensions = SelectImage_Dimensions(image, XSize=xsize, YSize=ysize, YIndex=yindex)

   ; Flip the image if required.

IF flipimage THEN image = Reverse(image, yindex+1)

   ; Create the widgets.

IF N_Elements(group_leader) NE 0 THEN BEGIN
   tlb = Widget_Base(Title='Select Image File', Column=1, /Base_Align_Center, $
      /Modal, Group_Leader=group_leader)
ENDIF ELSE BEGIN
   tlb = Widget_Base(Title='Select Image File', Column=1, /Base_Align_Center)
ENDELSE

fileSelectBase = Widget_Base(tlb, Column=1, Frame=1)
buttonBase = Widget_Base(tlb, Row=1)

   ; Define file selection widgets.

filenameID = FSC_FileSelect(fileSelectBase, Filename=filename, ObjectRef=filenameObj,$
   Directory=directory, Event_Pro='SelectImage_FilenameEvents', Filter=filter, _Extra=extra)
fsrowbaseID = Widget_Base(fileSelectBase, Row=1, XPad=10)
filelistID = Widget_List(fsrowbaseID, Value=theFiles, YSize = 10, Event_Pro='SelectImage_ListEvents')
spacer = Widget_Label(fsrowbaseID, Value="     ")
previewID = Widget_Draw(fsrowbaseID, XSize=previewSize, YSize=previewSize)
spacer = Widget_Label(fsrowbaseID, Value="     ")
labelBaseID = Widget_Base(fsrowbaseID, Column=1, /Align_Left)
IF fileInfo.channels EQ 3 THEN imageType = "True-Color Image" ELSE imageType = '8-Bit Image'
imageDataType = Size(image, /TNAME)
labeltypeID = Widget_Label(labelBaseID, Value=imageType, /Dynamic_Resize)
labelxsizeID = Widget_Label(labelBaseID, Value="X Size: " + StrTrim(xsize, 2), /Dynamic_Resize)
labelysizeID = Widget_Label(labelBaseID, Value="Y Size: " + StrTrim(ysize, 2), /Dynamic_Resize)
labeldataTypeID = Widget_Label(labelBaseID, Value="Type: " + imageDataType, /Dynamic_Resize)
labelminvalID = Widget_Label(labelBaseID, Value="Min Value: " + StrTrim(Fix(Min(image)),2), /Dynamic_Resize)
labelmaxvalID = Widget_Label(labelBaseID, Value="Max Value: " + StrTrim(Fix(Max(image)),2), /Dynamic_Resize)

   ; Size the draw widget appropriately.

IF xsize NE ysize THEN BEGIN
   aspect = Float(xsize) / ysize
   IF aspect GT 1 THEN BEGIN
      wxsize = previewSize
      wysize = previewSize / aspect
   ENDIF ELSE BEGIN
      wysize = previewSize
      wxsize = previewSize / aspect
   ENDELSE
   Widget_Control, previewID, Draw_XSize=wxsize, Draw_YSize=wysize
ENDIF

   ; Define buttons widgets.

button = Widget_Button(buttonBase, Value='Cancel', Event_Pro='SelectImage_Action')
filterID = Widget_Button(buttonBase, Value='Image Type', /Menu, Event_Pro='SelectImage_SetFilter')
button = Widget_Button(filterID, Value='BMP Files', UValue='*.bmp')
button = Widget_Button(filterID, Value='DICOM Files', UValue='*.dcm')
IF havegif THEN button = Widget_Button(filterID, Value='GIF Files', UValue='*.gif')
button = Widget_Button(filterID, Value='PICT Files', UValue='*.pict')
button = Widget_Button(filterID, Value='PNG Files', UValue='*.png')
button = Widget_Button(filterID, Value='JPEG Files', UValue='*.jpg')
button = Widget_Button(filterID, Value='TIFF Files', UValue='*.tif')
button = Widget_Button(buttonBase, Value='Flip Image', Event_Pro='SelectImage_FlipImage')
acceptID = Widget_Button(buttonBase, Value='Accept', Event_Pro='SelectImage_Action')

Widget_Control, tlb, /Realize
Widget_Control, previewID, Get_Value=previewWID
Widget_Control, acceptID, Input_Focus=1

   ; Set up RGB color vectors.

IF N_Elements(r) EQ 0 THEN r = Bindgen(!D.Table_Size)
IF N_Elements(g) EQ 0 THEN g = Bindgen(!D.Table_Size)
IF N_Elements(b) EQ 0 THEN b = Bindgen(!D.Table_Size)
WSet, previewWID
TVLCT, r, g, b

   ; Display the image.

IF (Min(image) LT 0) OR (Max(image) GT (!D.Table_Size-1)) THEN $
   TVImage, BytScl(image, Top=!D.Table_Size-1), /Keep_Aspect, /NoInterpolation ELSE $
   TVImage, image, /Keep_Aspect, /NoInterpolation

   ; Set up information to run the program.

storagePtr = Ptr_New({cancel:1, image:Ptr_New(image), fileInfo:Ptr_New(fileInfo), outdirectory:"", outfilename:""})

info = { storagePtr: storagePtr, $           ; The "outside the program" storage pointer.
         previewID: previewID, $             ; The ID of the preview draw widget.
         previewWID: previewWID, $           ; The window index number of the preview draw widget.
         r:r, $                              ; The R color vector.
         g:g, $                              ; The G color vector.
         b:b, $                              ; The B color vector.
         theFiles: Ptr_New(theFiles), $      ; The current list of files in the directory.
         filenameID: filenameID, $           ; The identifier of the FileSelect compound widget.
         fileListID: fileListID, $           ; The identifier of the file list widget.
         flipimage:flipimage, $              ; A flag to flip the image Y order.
         previewSize: previewSize, $         ; The default size of the preview window.
         filter: filter, $                   ; The file filter.
         filenameObj: filenameObj, $         ; The FileSelect compound widget object reference.
         dataDirectory: directory, $         ; The current data directory.
         labelmaxvalID: labelmaxvalID, $     ; The ID of the Max Value label.
         labelminvalID: labelminvalID, $     ; The ID of the Max Value label.
         labelTypeID: labelTypeID, $         ; The ID of the Image Type label.
         labelXSizeID: labelXSizeID, $       ; The ID of the X Sizee label.
         labelYSizeID: labelYSizeID, $       ; The ID of the Y Size label.
         labelDataTypeID: labelDataTypeID $  ; The ID of the Data Type label.
       }
Widget_Control, tlb, Set_UValue=info, /No_Copy

   ; Blocking or modal widget mode, depending upon presence of GROUP_LEADER.

XManager, "selectimage", tlb, Cleanup='SelectImage_Cleanup'

   ; Return collected information.

cancel = (*storagePtr).cancel
fileInfo = *(*storagePtr).fileInfo
image = *((*storagePtr).image)
outDirectory = (*storagePtr).outDirectory
outFilename = (*storagePtr).outFilename
Ptr_Free, (*storagePtr).image
Ptr_Free, (*storagePtr).fileInfo
Ptr_Free, storagePtr

   ; Restore start directory.

CD, startDirectory
IF cancel EQ 1 THEN RETURN, 0 ELSE RETURN, image
END