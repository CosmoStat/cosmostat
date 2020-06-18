;+
; NAME:
;     CIRCOPLOT
;
; PURPOSE:
;     Plots a circle of user-requested radius over an 
;     image projection. Circle is centered on the cursor position
;     specified either by the user clicking the left mouse button (screen) 
;     or providing a coordinate list (postscript).
;     Option also exists to indicate MAP-visible portions of circle.
;     
; CALLING SEQUENCE:
;     circoplot [,radius=radius] [,proj={'M','Z'} [,color=num] 
;               [ /visible] [,coord={'G','E','C'}, 
;               [/ps [,image=image] [,ctable=num] [,coorlist=coorlist]
;
; INPUTS:
;   Screen mode:
;     User positions cursor over image and clicks left mouse button to
;      specify circle center.
;     Routine remains active until the user clicks the right mouse button.
;   Postscript mode:
;     User provides a list of centers via coorlist keyword.
;
; OUTPUTS:  
;   Screen mode:
;     Output consists of a circle outline overlayed on an already existing
;     image. A small cross is drawn at the requested the circle center 
;     position. Multiple circles may be drawn during the one call to circoplot.
;     If /visible is specified, triangles are plotted over the MAP-visible 
;     segments of the circle.
;   Postscript mode:
;     Outputs similar to screen mode and written to file circoplot.ps.
;     Plot symbols have been tailored for hardcopy.
;     Postscript is written in portrait deliberately to allow user the option
;      of using eps and preview options within ghostview.
;     
;
; KEYWORDS:
;     radius - float     - Radius of the circle, in degrees.
;                          Default value = 141 deg.
;
;     proj - char string - Single character specifying the projection type
;                          of the underlying image.  Either 'M' (Mollweide) or
;                          'Z' (Zenithal Equal area) are allowed.  Case 
;                          insensitive.   Defaults to 'M'.
; 
;     color - byte -       The color (range 0-255) used to plot the 
;                          circle and plot symbols. Defaults to 0
;
;     /visible -   Set this keyword to have circoplot indicate those
;                          portions of the circle in which the spin axis
;                          is within 22.5 deg of the ecliptic.  
;                          The result is coordinate system
;                          dependent -- see the coord keyword.
;
;     coord - char string -THIS KEYWORD IS MEANINGFUL ONLY IF /VISIBLE 
;                           IS SPECIFIED. 
;                          Single character specifying the coordinate 
;                          system of the projection. Three systems are
;                          recognized: 'E' (Ecliptic J2000), 'G' (Galactic) or
;                          'C' (Celestial J2000).  Case insensitive.
;                          Defaults to 'G'.
;
;     ps - value=0 or 1 -  Set this keyword to direct the plot to a postscript
;                          file.  The output file is named 'circoplot.ps'.  
;                          If not set, plot will be directed to the screen 
;                          (windows device).
;
;     coorlist - fltarr  - A list of circle center coordinates, in the 
;                          coordinate system OF THE UNDERLYING IMAGE.
;                          Format is [2,N] where N is the number of centers.
;                          Must be present when /ps is requested, and is only
;                          implemented for ps at present.
;
;     image - bytarr -     A byte-scaled image over which the circles will
;                          be plotted.  At present, only implemented for
;                          /ps, since user can already plot over an existing
;                          image on the screen.
;
;     ctable - long  -     A number specifying which standard IDL color table
;                          to load.  Only implemented for /ps.
;
;
;
; COMMON BLOCKS:
;     None.
;
; ROUTINES CALLED:
;     projxy2coord, mollweide_xy, zea_xy
;     
; EXAMPLE:
;     circoplot,radius=10,proj='Z',color=20
;     circoplot,/vis,color=100,coord='e'
;     circoplot,image=myimage,coorlist=[[120,50],[130,80]],/ps
;     
; COMMENTS:
;     The routine assumes that the image over which the circle is being
;     plotted fills the active screen window. If the
;     user is trying to overlay an image on the screen which does not
;     fill the window, then the overlay will be scaled incorrectly.
;     
; MODIFICATION HISTORY:
;     initial version, J. Weiland, 30 Apr 1999
;     added /visible option, JW, 17 Sep 1999
;     added postscript options, JW, Aug 2000
;     Use !DTOR and !RADEG  WL  Dec 2002
;
;-
;======================================================================
pro circoplot,radius=radius,proj=proj,color=color,visible=visible,coord=coord,$
              ps=ps,image=image,ctable=ctable,coorlist=coorlist

;
; draw a circle of angular radius r (units= degrees) on the screen or
;  to postscript file.
; circle is centered on [x0,y0] in NORMALIZED coordinates (screen), or
;  on user-specified center coordinates (ps).
;
; In screen mode, image on screen MUST fill the plot window in order for 
;  this to work.
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;
; default Radius, in degrees 
;
if (keyword_set(radius) eq 0) then radius = 141.

;
; Projection Type
;
if keyword_set(proj) then begin
   proj = strupcase(strtrim(proj,2))
   if ((proj ne 'M') and (proj ne 'Z')) then begin
     print, 'Only Mollweide (M) or Zenithal Equal area (Z) projections allowed.'
     return
   endif
endif else begin
   proj = 'M'                      ;default = Mollweide
endelse
;
; Coordinate Type
;
if keyword_set(coord) then begin
   coord = strupcase(strtrim(coord,2))
   if ((coord ne 'E') and (coord ne 'G') and (coord ne 'C')) then begin
     print, 'Ecliptic(E), Celestial(C) or Galactic(G) coordinates only '
     return
   endif
endif else begin
   coord = 'G'                      ;default = Galactic
endelse
;
; Color (defaults to 0 if not set)
;
if keyword_set(color) eq 0 then color=0


IF keyword_set(ps) THEN BEGIN                   ;hardcopy
   ;
   ; check to make sure that the center coordinate list is provided in the
   ; correct format
   ;
   if keyword_set(coorlist) then begin
     diminfo = size(coorlist)     
     if (diminfo[0] eq 2) then begin
       n_center = diminfo[2]
     endif else if (diminfo[0] eq 1 AND diminfo[1] eq 2) then begin
       n_center = 1
     endif else begin
       print,'Coorlist needs to be array of format [2,N]'
       return
     endelse   
   endif else begin
     print, 'A coordinate list must be provided in /ps mode!'
     return
   endelse
   ;
   ;  Set the plotting environment
   ;
   user_device = !d.name
   set_plot,'ps'
   ;
   ;  the preview keyword is added for importation into powerpoint
   ; 
   device,/port,/color,bits=8,/inches,xs=8,ys=4,yoff=4,xoff=0.25
     ; /encapsulated,/preview     ;doesn't seem to work for powerpoint
   device,file='circoplot.ps'
   if keyword_set(ctable) then begin
      loadct,ctable
      tvlct,r,g,b,/g
      r[0]=255 & b[0]=255 &g[0]=255
      tvlct,r,g,b
   endif else begin
      loadct,27
   endelse
   if keyword_set(image) then tv,image,/inches,xs=8,ys=4
   ;
   ; loop over the provided center coordinates
   ;
   for i=0L,n_center-1L do begin
      
     ; make darn sure the centers are floating point numbers.  Integers will
     ; not necessarily convert correctly to xy coords.
     ;
     lon0 = float(coorlist[0,i])
     lat0 = float(coorlist[1,i])

     ;
     ; define the circle perimeter points 
     ;
     xi  = findgen(361) * !DTOR
     rr  = radius *!DTOR
     ln0 = lon0*!RADEG
     lt0 = lat0*!RADEG

     if (radius gt 5) then begin   ; use rigorous defn
  
      lat = asin((cos(xi)*sin(rr) +  $
            cos(rr)*tan(lt0))/(cos(lt0)+ tan(lt0)*sin(lt0)))
      lon = ln0 +  $
            atan ((sin(xi)*sin(rr)*cos(lt0)),(cos(rr) - sin(lt0)*sin(lat)))
      lat = lat *!RADEG
      lon = lon *!RADEG

     endif else begin         ; use sloppier defn, more robust for small radius
 
      lat = lat0 + radius*cos(xi)
      lon = lon0 + radius*sin(xi)/cos(lat*!dtor)

     endelse

     ; Ensure that lon is in the range [0,360) before feeding to conversion 
     ; routines
     ;
     neg = where(lon lt 0, Nneg)
     if (Nneg GT 0) then lon[neg] = lon[neg] + 360.

     ;
     ; If VISIBLE is set, then find the spin axis locations within 22.5 degrees
     ; of the ecliptic plane.  Assume the spin axis locations occur along a
     ; circle of radius= radius/2.
     ; Assumes the user will not do something unusual like set the radius to
     ; a small number when asking for this option.
     ;

     IF keyword_set(visible) then begin
     
        rad2 = (radius/2.)*!DTOR
        lat2 = asin((cos(xi)*sin(rad2) +   $
                      cos(rad2)*tan(lt0))/(cos(lt0)+ tan(lt0)*sin(lt0)))
        lon2 = ln0 + atan ((sin(xi)*sin(rad2)*cos(lt0)),   $
                           (cos(rad2) - sin(lt0)*sin(lat2)))
        lat2 = lat2 *!RADEG
        lon2 = lon2 *!RADEG
        ;
        ; Ensure that lon2 is in the range [0,360) before feeding to conversion
        ; routines
        ;
        test2 = where(lon2 lt 0)
        if (test2[0] ne -1) then lon2[test2] = lon2[test2] + 360.
        ;       
        ;        Convert the r/2 circle from the 'coord' system into ecliptic
        ;        if necessary, and to unit vectors.
        ;
        if (coord eq 'E') then begin
           trans_code = 'll2u'
           coortrans,[[lon2],[lat2]],uvec2,trans_code
        endif else begin
           trans_code = coord+'2E'
           coortrans,[[lon2],[lat2]],uvec2,trans_code
        endelse
        ;
        ;        Take dot product with normal to ecliptic plane
        ;
        normal = fltarr(361,3) & normal[*,2] = 1.0
        dotp = total(uvec2 * normal,2)
        limit = cos((90.-22.5)*!dtor)
        vis  = where(abs(dotp) le limit)
     ENDIF
     ;
     ; 
  
     ;
     ; Convert [lon,lat] array back into [x,y]
     ;

     case proj of
  
     'M': begin
       
         ; convert the center position (lon0,lat0) to (x0,y0) and plot
         ; a cross there.
         ;
         mollweide_xy,lon0,lat0,x0,y0
         x0 = (x0 +2.)/4.
         y0 = (y0 +1.)/2.
         plots,[x0-.01,x0+.01],[y0,y0],/norm,color=color
         plots,[x0,x0],[y0-.02,y0+.02],/norm,color=color
              ;try to make cross thicker
           plots,[x0-.01,x0+.01],[y0+0.001,y0+0.001],/norm,color=color
           plots,[x0-.01,x0+.01],[y0-0.001,y0-0.001],/norm,color=color
           plots,[x0+0.001,x0+0.001],[y0-.02,y0+.02],/norm,color=color
           plots,[x0-0.001,x0-0.001],[y0-.02,y0+.02],/norm,color=color
      
         ;
         ; Now plot the circle.
         ; call Gary's routine to get x in [-2,2] and y in [-1,1]
         ;
         mollweide_xy,lon,lat,x,y
         ;
         ; Scale to normalized coordinates
         ;
         x = (x +2.)/4.
         y = (y +1.)/2.
         ;
         ; Plot the grid point by point and suppress the line when crossing 
         ; the 180 degree meridian
         PLOTS, X[0], Y[0], /norm, Color=color
         FOR J = 1,n_elements(x)-1 DO BEGIN
             IF( ABS(lon[J]-lon[J-1]) LT 180. )THEN BEGIN
               PLOTS, X[J], Y[J], /norm, /Continue, Color=color,thick=2
             ENDIF ELSE BEGIN
               PLOTS, X[J], Y[J], /norm, Color=color,thick=2 
             ENDELSE
         ENDFOR
  
         IF (KEYWORD_SET(VISIBLE))  THEN BEGIN
   
            aa = findgen(16)*(!pi*2./16.)
            usersym,cos(aa),sin(aa),/fill

            IF (vis[0] ne -1) then begin
              mollweide_xy,lon[vis],lat[vis],x2,y2
              x2 = (x2 +2.)/4.
              y2 = (y2 +1.)/2.
              plots,x2,y2,/norm,color=color,psym=8,symsize=0.75
            ENDIF
  
         ENDIF

       end

  
     'Z': begin
         ; convert the center position (lon0,lat0) to (x0,y0) and plot
         ; a cross there.
         ;
         zea_xy,lon0,lat0,x0,y0
         x0 = (x0 +2.)/4.
         y0 = (y0 +1.)/2.
         plots,[x0-.005,x0+.005],[y0,y0],/norm,color=color
         plots,[x0,x0],[y0-.01,y0+.01],/norm,color=color
         ;
       ; Now plot the circle.
       ;
       ; given lon, lat in degrees, get x in [-2,2] and y in [-1,1]
       ;
       zea_xy,lon,lat,x,y
       ;
       ; Scale to normalized coordinates
       ;
       x = (x +2.)/4.
       y = (y +1.)/2.
       ;
       ; Plot the grid point by point and suppress the line when crossing the 
       ; equator
       PLOTS, X[0], Y[0], /norm, Color=color
       FOR J = 1,n_elements(x)-1 DO BEGIN
          IF(lat[J]*lat[J-1] gt 0) THEN BEGIN
            PLOTS, X[J], Y[J], /norm, /Continue, Color=color 
          ENDIF ELSE BEGIN
            PLOTS, X[J], Y[J], /norm, Color=color 
          ENDELSE
       ENDFOR
  
       IF (KEYWORD_SET(VISIBLE)) THEN BEGIN
  
          IF (vis[0] ne -1) THEN BEGIN
            zea_xy,lon[vis],lat[vis],x2,y2
            x2 = (x2 +2.)/4.
            y2 = (y2 +1.)/2.
            plots,x2,y2,/norm,color=color,psym=5
          ENDIF

       ENDIF

       end

      endcase

   endfor

   device,/close
   set_plot, user_device


ENDIF ELSE BEGIN                                ;screen
   ;
   print,'Place cursor on central position and click left mouse button'
   print,' in order for a circle of specified radius (degrees) to be drawn.'
   print,'Click on right mouse button to quit.'
   !mouse.button=0


   WHILE !mouse.button ne 4 DO BEGIN
 
     cursor,x0,y0,1,/norm               ; return x,y in normalized coordinates
     if (!mouse.button eq 4) then goto, cleanup
   ;
   ; Place a cross at the central position
   ; x:y = 2:1
   ;
   plots,[x0-.005,x0+.005],[y0,y0],/norm,color=color
   plots,[x0,x0],[y0-.01,y0+.01],/norm,color=color

   ;
   ; Convert (x0,y0) to (lon0,lat0) for specified projection
   ;

   projxy2coord,x0,y0,lon0,lat0,proj=proj

   ;
   ; define the circle perimeter points 
   ;
   xi  = findgen(361) * !dtor
   rr  = radius *!dtor
   ln0 = lon0*!dtor
   lt0 = lat0*!dtor

   if (radius gt 5) then begin   ; use rigorous defn

    lat = asin((cos(xi)*sin(rr) +  $
          cos(rr)*tan(lt0))/(cos(lt0)+ tan(lt0)*sin(lt0)))
    lon = ln0 + atan ((sin(xi)*sin(rr)*cos(lt0)),(cos(rr) - sin(lt0)*sin(lat)))
    lat = lat *!radeg
    lon = lon *!radeg

   endif else begin         ; use sloppier defn, more robust for small radius
 
    lat = lat0 + radius*cos(xi)
    lon = lon0 + radius*sin(xi)/cos(lat*!dtor)

   endelse

   ; Ensure that lon is in the range [0,360) before feeding to conversion 
   ; routines
   ;
   test = where(lon lt 0)
   if (test[0] ne -1) then lon[test] = lon[test] + 360.

   ;
   ; If VISIBLE is set, then find the spin axis locations within 22.5 degrees
   ; of the ecliptic plane.  Assume the spin axis locations occur along a
   ; circle of radius= radius/2.
   ; Assumes the user will not do something unusual like set the radius to
   ; a small number when asking for this option.
   ;

   IF keyword_set(visible) then begin
   
      rad2 = (radius/2.)*!dtor
      lat2 = asin((cos(xi)*sin(rad2) +   $
                    cos(rad2)*tan(lt0))/(cos(lt0)+ tan(lt0)*sin(lt0)))
      lon2 = ln0 + atan ((sin(xi)*sin(rad2)*cos(lt0)),   $
                         (cos(rad2) - sin(lt0)*sin(lat2)))
      lat2 = lat2 * !radeg
      lon2 = lon2 * !radeg
      ;
      ; Ensure that lon2 is in the range [0,360) before feeding to conversion
      ; routines
      ;
      neg2 = where(lon2 lt 0, Nneg2)
      if (Nneg2 GT 0) then lon2[neg2] = lon2[neg2] + 360.
      ;       
      ;        Convert the r/2 circle from the 'coord' system into ecliptic
      ;        if necessary, and to unit vectors.
      ;
      if (coord eq 'E') then begin
         trans_code = 'll2u'
         coortrans,[[lon2],[lat2]],uvec2,trans_code
      endif else begin
         trans_code = coord+'2E'
         coortrans,[[lon2],[lat2]],uvec2,trans_code
      endelse
      ;
      ;        Take dot product with normal to ecliptic plane
      ;
       normal = fltarr(361,3) & normal[*,2] = 1.0
       dotp = total(uvec2 * normal,2)
       limit = cos((90.-22.5)*!dtor)
       vis  = where(abs(dotp) le limit)
   ENDIF
   ;
   ; 

   ;
   ; Convert [lon,lat] array back into [x,y]
   ;

   case proj of

   'M': begin
     
       ;
       ; call Gary's routine to get x in [-2,2] and y in [-1,1]
       ;
       mollweide_xy,lon,lat,x,y
       ;
       ; Scale to normalized coordinates
       ;
       x = (x +2.)/4.
       y = (y +1.)/2.
       ;
       ; Plot the grid point by point and suppress the line when crossing the 
       ; 180 degree meridian
       PLOTS, X[0], Y[0], /norm, Color=color
       FOR J = 1,n_elements(x)-1 DO BEGIN
           IF( ABS(lon[J]-lon[J-1]) LT 180. )THEN BEGIN
             PLOTS, X[J], Y[J], /norm, /Continue, Color=color 
           ENDIF ELSE BEGIN
             PLOTS, X[J], Y[J], /norm, Color=color 
           ENDELSE
       ENDFOR

       IF (KEYWORD_SET(VISIBLE))  THEN BEGIN
 
          IF (vis[0] ne -1) then begin
            mollweide_xy,lon[vis],lat[vis],x2,y2
            x2 = (x2 +2.)/4.
            y2 = (y2 +1.)/2.
            plots,x2,y2,/norm,color=color,psym=5
          ENDIF

       ENDIF

     end


   'Z': begin
     ;
     ; given lon, lat in degrees, get x in [-2,2] and y in [-1,1]
     ;
     zea_xy,lon,lat,x,y
     ;
     ; Scale to normalized coordinates
     ;
     x = (x +2.)/4.
     y = (y +1.)/2.
     ;
     ; Plot the grid point by point and suppress the line when crossing the 
     ; equator
     PLOTS, X[0], Y[0], /norm, Color=color
     FOR J = 1,n_elements(x)-1 DO BEGIN
        IF(lat[J]*lat[J-1] gt 0) THEN BEGIN
          PLOTS, X[J], Y[J], /norm, /Continue, Color=color 
        ENDIF ELSE BEGIN
          PLOTS, X[J], Y[J], /norm, Color=color 
        ENDELSE
     ENDFOR

     IF (KEYWORD_SET(VISIBLE)) THEN BEGIN

        IF (vis[0] ne -1) THEN BEGIN
          zea_xy,lon[vis],lat[vis],x2,y2
          x2 = (x2 +2.)/4.
          y2 = (y2 +1.)/2.
          plots,x2,y2,/norm,color=color,psym=5
        ENDIF

     ENDIF

     end

    endcase
   ;

   ENDWHILE

   ;
   cleanup:
   
   !mouse.button=0

ENDELSE

;
end




