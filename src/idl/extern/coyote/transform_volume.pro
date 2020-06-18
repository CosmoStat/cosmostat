;+
; NAME:
;       TRANSFORM_VOLUME
;
; PURPOSE:
;
;       The purpose of this program is to transform (e.g., rotate,
;       scale, and translate) a 3D array or volume.
;
; AUTHOR:
;
;       Martin Downing,
;       Clinical Research Physicist,
;       Grampian Orthopaedic RSA Research Centre,
;       Woodend Hospital, Aberdeen, AB15 6LS.
;       Pnone: 01224 556055 / 07903901612
;       Fa: 01224 556662
;       E-mail: m.downing@abdn.ac.uk
;
; CATEGORY:
;
;      Mathematics, graphics.
;
; CALLING SEQUENCE:
;
;      result = TRANSFORM_VOLUME( volume )
;
; INPUTS:
;
;       volume:    The 3D array or volume to be transformed.
;
; OPTIONAL KEYWORDS:
;
;      BUFFER_SIZE: To reduce memory overhead the routine processes the job in chunks, the number
;         of elements of which can be set using the BUFFER_SIZE keyword, set this keyword to
;         0 to force the whole array to be processed at one time. The default value is 128.
;
;      MISSING: The value to return for transformed values outside the bounds of
;         the volume. (Passed to the INTERPOLATE function.) Default is 0.
;
;      T3DMAT: The homogeneous transforamtion matrix. If this keyword is not present,
;         the following keywords can be used to create a homogeneous transformation matrix:
;
;         ROTATION - The rotation vector [rx,ry,rz]. The order of rotation is ZYX.
;         TRANSLATE - The translation vector [tx,ty,tz].
;         SCALE - The scale vector [sx,sy,sz].
;         CENTRE_ROTATION - The centre of rotation [cx,cy,cz].
;
; OUTPUTS:
;
;       result:    The transformed array or volume.
;
; COMMON BLOCKS:
;
;       None.
;
; DEPENDENCIES:
;
;       The program uses the library INTERPLOLATE routine, which currently (IDL 5.4)
;       uses linear interpolation. Note that the operation is performed in chunks,
;       each of which is independant of the result of the others, so the operation
;       could easiliy be parallelised.
;
; MODIFICATION HISTORY:
;
;       Written by: Martin Downing, 16 September 2001.
;       Added MISSING keyword. Removed INPLACE keyword. 25 Nov 2001. MD
;-
;###########################################################################
;
; LICENSE
;
; This software is OSI Certified Open Source Software.
; OSI Certified is a certification mark of the Open Source Initiative.
;
; Copyright © 2001 Martin Downing.
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


FUNCTION Transform_Volume, volume, Rotation=rotation, $
    Scale=scale, Translate=translate, Centre_Rotation=centre_rotation, $
    T3Dmat=t3dmat, Buffer_Size=buffer_size, Missing=missing

   ; Error handling.

Catch, theError
IF theError NE 0 THEN BEGIN
   Catch, /Cancel
   ok = Dialog_Message(!Error.State.Msg)
   RETURN, -1
ENDIF

   ; Find the dimensions of the volume.

 s = Size(volume)
 sx=s[1] & sy=s[2] & sz=s[3]
 st = sx*sy*sz

vol_t = volume
IF N_Elements(missing) THEN missing = 0

   ; Create a transform matrix, if one is not provided.

IF N_Elements(t3dmat) EQ 0 THEN begin

   IF N_Elements(rotation) EQ 0 THEN rotation =[0,0,0]
   IF N_Elements(centre_rotation) EQ 0  THEN centre_rotation=[(sx-1)/2.0,(sy-1)/2.0,(sz-1)/2.0]
   IF N_Elements(translate) EQ 0 THEN translate =[0,0,0]
   IF N_Elements(scale) EQ 0 THEN scale =[1,1,1]

   T3D, /Reset, Translate = -centre_rotation
   T3D, Rotate=rotation
   T3D, Translate= centre_rotation + translate, Scale=scale
   t3dmat = !P.T

ENDIF

   ; Check buffer size. The size 128 is optimim on my system, You may
   ; want to try other values.

 IF N_Elements(buffer_size) EQ 0 THEN buffer_size = 128
 IF buffer_size LE 0 THEN buffer_size = st

   ; Perform the transformations.

 FOR j=0L,(st-1),buffer_size DO BEGIN

      ; Account for possible odd last chunk.

   bufsize = buffer_size < (st-j)

      ; Generate volume coordinates by interpolating temporary array of volume indices.

   i = j + Lindgen(bufsize)
   coords = [ [(i MOD sx)],[((i / sx) MOD (sy))], [(i / (sx * sy))], [Replicate(1b, bufsize)]]
   coords = Temporary(coords) # t3dmat
   vol_t[j:j+bufsize-1] = Interpolate(volume, coords[*,0], coords[*,1], coords[*,2], Missing=missing)
ENDFOR

   ; Return the transformed volume.

RETURN, vol_t
END
