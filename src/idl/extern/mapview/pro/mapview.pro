;+
; NAME:
;	MAPVIEW
;
; PURPOSE:
;	Interactive HEALPix map viewer.
;
; CALLING SEQUENCE:
;	mapview,  [file=filename, array=array, /ring]
;
; INPUTS:
;	None
;
; OPTIONAL INPUT KEYWORD PARAMETER:
;	file  - filename of HealPix map to view, scalar string
;	array - data array containing map
;
;       If no keywords are supplied then the user must enter a HealPix FITS 
;       filename into the widget.
; OUTPUTS:
;	None
;
; EXAMPLES:
;	Display map(s) from a file:
;
;		mapview, file='map_w1_imap_yr1_v1.fits'
;
;	Display map(s) from an array:
;
;		mapview, array=array, [/ring]
;
;	If viewing a single map, array should have dimensions [n_pix].
;	For multiple maps, array should have dimensions [n_pix,n_maps].
;	The array is assumed to be in 'nested' format. If not, specify
;	the /ring keyword.
;
; COMMON BLOCKS:
;	mvCommon
;
; COMMENTS:
;       This widget displays a HealPix sky map, and lets the user zoom in 
;       on regions of the map on the fly.
;
;       The widget displays two images on the right: the upper image is an 
;       all-sky map (a Mollweide projection in Galactic coordinates), with a 
;       circle drawn over the zoomed region.    The lower 512 x 512 image 
;       displays the  zoomed region.    The images are manipulated using the 
;       mouse cursor and the control panels on the left.   
;       
;       The Control panel functions are (from top to bottom)
;       (1) "File"  -- Enter a new HealPix FITS file to display
;       (2) "Select Map" -- Choose whether to display temperature or number of
;                           observations
;       (3) "Cursor Info" -- Informational panel giving the Galactic 
;            coordinates, Healpix pixel number and pixel value of the current
;            cursor position in either the all-sky or zoomed image.   (The 
;            panel is disabled if the cursor is not over either image.)
;       (4) "Color Scaling" - lets the user choose linear, logarithmic or 
;             histogram equalization scaling, and a color lookup table.  Note
;             temperature maps may appear entirely black with a linear scaling. 
;       (5) "Zoom Region" - User can type in the Galactic coordinates of the
;             region to be zoomed, and then press the draw button.    
;            Alternatively, if the zoom region is selected on the all-sky map 
;            with the cursor, then the chosen coordinates are displayed.    A 
;            slider allows the user to choose the area of the zoomed region; 
;            the larger the region, the smaller the zoom factor.   The actual 
;            pixel scale of the zoomed region is displayed at the terminal.
;       (6)  "Print buttons" - Push selected button to write either the all-sky
;            map to a PNG file 'map.png' or the zoomed map to a file 'zoom.png' 
;       (7) "Quit" button - exits the widget
;
;       The mouse cursor functions are as follows:
;        Move the cursor over either the all-sky or zoomed image to have values
;        displayed in the "Cursor Info" panel.   Press any mouse button on 
;        the all sky image to define a new center for the zoomed image.  Drag
;        the outline of the circle defining the zoomed image to increase or
;        decrease the size of the circle.    Press any mouse button on the 
;        zoomed image to move that position to the center of the zoomed image.
;
; MODIFICATION HISTORY:
;       Written, Mike Nolta	May 2002      
;       Modify calls to EULER_MATRIX_NEW() and WMAP_DATA2GNOM for HEALPix 1.2 
;               compatibility    W. Landsman    14 Feb 2003 
;	Added array keyword	M. Nolta	07 Apr 2003
;	Added file output, fixed bug with ring maps
;				M. Nolta	15 Jul 2003
;       Added filename prompt, improve speed by bypassing read_fits_s
;                               W. Landsman     07 Jul 2004
;
;-
;

function int2str, n
	return, string( format='(I0)', n )
end

function keyword_def, x
	return, n_elements(x) gt 0
end

function num2str, num
	return, strcompress( string(num), /remove_all )
end
;
;+
; NAME:
;	Pixmap
;
; PURPOSE:
;	Pixmap object
;
; CATEGORY:
;	Object
;
; CALLING SEQUENCE:
;	pixmap = obj_new( 'Pixmap', xsize, ysize )
;
; INPUTS:
;	xsize:	width in pixels
;	ysize:	height in pixels
;
; OPTIONAL INPUT PARAMETERS:
;	None
;
; KEYWORD PARAMETERS:
;	None
;
; OUTPUTS:
;
; OPTIONAL OUTPUT PARAMETERS:
;
; COMMON BLOCKS:
;	None
;
; SIDE EFFECTS:
;	None
;
; RESTRICTIONS:
;	None
;
; PROCEDURE:
;
; MODIFICATION HISTORY:
;
;	$Log$
;
;-

pro sky2vec, l, b, vec, deg=deg, double=double

	n = n_elements(l)

	if keyword_set(double) then begin
		d2r = !dpi/180
		vec = dblarr( n, 3 )
	endif else begin
		d2r = !pi/180
		vec = fltarr( n, 3 )
	endelse

	if not keyword_set(deg) then begin
		theta = 90*d2r - b 
		phi = l
	endif else begin
		theta = (90. - b) * d2r
		phi = l * d2r
	endelse

	sin_theta = sin(theta)
	vec[*,0] = sin_theta * cos(phi)
	vec[*,1] = sin_theta * sin(phi)
	vec[*,2] = cos(theta)

	if size(l,/n_dimen) EQ 0 then vec = transpose(vec)
end

function Pixmap::init, xsize, ysize

	old_window = !d.window

	window, /free, /pixmap, xsize=xsize, ysize=ysize
	self.pixmap = !d.window
	self.xsize = xsize
	self.ysize = ysize

	wset, old_window
	return, 1
end

pro Pixmap::cleanup

	wdelete, self.pixmap
end

pro Pixmap::draw_begin

	self.old_window = !d.window
	wset, self.pixmap
end

pro Pixmap::draw_end

	wset, self.old_window
end

pro Pixmap::copy_to_device

	device, copy=[0,0,self.xsize,self.ysize,0,0,self.pixmap]
end

pro Pixmap__define

	struct = { Pixmap, $
		old_window	: -1L, $
		pixmap		: -1L, $
		xsize		: -1L, $
		ysize		: -1L $
	}
end

function MollewideProjection::init

	self.coordsys = 'G'
	self.flip = 0
	return, 1
end

pro MollewideProjection::hmap_to_image, hmap, image, $
	index=index, name=name, _extra=extra

	;map = hmap->get_map( index=index, name=name )
	hmap->read_map, map, index=index, name=name
	if hmap->ring() then map = reorder( map, in='RING', out='NESTED' )

	healpix_to_image, map, image, project=1, /silent, _extra=extra
end

pro MollewideProjection::uv_to_lb, u, v, lon_deg, lat_deg

;pro mv_moll2pix, nx, ny, hmap, x, y, id_pix, lon_deg, lat_deg
;+
; convert a uv position on a Mollweide map into a pixel number and
; (long, lat)
; only for scalar input
;
;-
	; longitude increase leftward by default (astro convention)
	if self.flip eq 1 then flipconv = +1 else flipconv = -1

	id_pix  = -1
	lon_deg = -1000.
	lat_deg = -1000.

	;xpos = 4.*float(x)/float(nx) - 2.
	;ypos = 2.*float(y)/float(ny) - 1.

	xpos = u
	ypos = v

	if (xpos^2 / 4. + ypos^2 le 1.) then begin

		nested = keyword_set(nest)
		v1 = ypos
		u1 = xpos
		s1 = SQRT((1.d0 - v1)*(1.d0 + v1))
		a1 = ASIN( v1 )

		z = 2./!PI * (a1 + v1*s1)

		; lon in [-pi,pi], minus sign is here to fit astro convention
		phi = (flipconv*!Pi/2.) * u1/s1
		sz = SQRT( (1. - z)*(1. + z) )
		vector = [[sz * COS(phi)], [sz * SIN(phi)], [z]]
		vec2ang, vector, lat_deg, lon_deg, /astro

	endif
end

pro MollewideProjection::uv_circle, lon0, lat0, r, u, v, lon=lon

	self->circle, lon0, lat0, r, u, v, lon=lon
end

pro MollewideProjection::circle, lon0, lat0, r, u, v, lon=lon, normal=normal

	;; define the circle perimeter points 
	xi  = findgen(361) * !pi/180.
	rr  = r *!pi/180.
	ln0 = lon0*!pi/180.
	lt0 = lat0*!pi/180.

	if ( r gt 5 ) then begin	;; use rigorous defn

		lat = asin((cos(xi)*sin(rr) +  $
			cos(rr)*tan(lt0))/(cos(lt0)+ tan(lt0)*sin(lt0)))
		lon = ln0 +  $
			atan ((sin(xi)*sin(rr)*cos(lt0)),(cos(rr) - sin(lt0)*sin(lat)))
		lat = lat *180./!pi
		lon = lon *180./!pi

	endif else begin

		;; use sloppier defn, more robust for small radius
		lat = lat0 + r*cos(xi)
		lon = lon0 + r*sin(xi)/cos(lat*!pi/180.)
	endelse

	; Ensure that lon is in the range [0,360) before feeding to conversion 
	; routines
	;
	test = where(lon lt 0)
	if (test[0] ne -1) then lon[test] = lon[test] + 360.

	mollweide_xy, lon, lat, u, v

	if keyword_set(normal) then begin
		u = (u + 2.)/4.
		v = (v + 1.)/2.
	endif
end

pro MollewideProjection__define

	struct = { MollewideProjection, $
		coordsys	: '', $
		flip		: 0 $
	}

end
;
; $Id$
;
;+
; NAME:
;    IMAGE_WIDGET
; PURPOSE:
;
; CATEGORY:
;
; CALLING SEQUENCE:
;
; INPUTS:
;	None
;
; KEYWORD PARAMETERS:
;	None
;
; OUTPUTS:
;	None
;
; COMMON BLOCKS:
;	None
;
; SIDE EFFECTS:
;	None
;
; RESTRICTIONS:
;	None
;
; MODIFICATION HISTORY:
;
;	$Log$
;-

function ImageWidget::init, parent, xsize, ysize, _extra=extra

	self.pixmap = obj_new( 'Pixmap', xsize, ysize )
	self.widget = ptr_new( widget_draw( parent, xsize=xsize, ysize=ysize, $
		_extra=extra ), /no_copy )
	self.xsize = xsize
	self.ysize = ysize

	return, 1
end

pro ImageWidget::cleanup

	ptr_free, self.image
	obj_destroy, self.pixmap
	ptr_free, self.widget
end

pro ImageWidget::erase

	old_window = !d.window
	wset, self->window_id()
	erase
	wset, old_window
end

pro ImageWidget::info, xsize=xsize, ysize=ysize

	xsize = self.xsize
	ysize = self.ysize
end

pro ImageWidget::load_image, image

	if ptr_valid( self.image ) then ptr_free, self.image
	self.image = ptr_new( image, /no_copy )
end

pro ImageWidget::draw_image

	if ptr_valid(self.image) then begin
		self.pixmap->draw_begin
		tv, *self.image
		plots, [100,100], [200,200], /device
		self.pixmap->draw_end
		self->refresh
	endif
end

pro ImageWidget::write_png, fn
	print, 'Writing', fn
	tvlct, r, g, b, /get
	write_png, fn, *self.image, r, g, b
end

pro ImageWidget::refresh

	old_window = !d.window
	wset, self->window_id()
	self.pixmap->copy_to_device
	wset, old_window
end

function ImageWidget::window_id

	widget_control, *self.widget, get_value=id
	return, id
end

pro ImageWidget::window_set

	wset, self->window_id()
end

pro ImageWidget__define

	struct = { ImageWidget, $
		image		: ptr_new(), $
		pixmap		: obj_new(), $
		widget		: ptr_new(), $
		xsize		: 0L, $
		ysize		: 0L $
	}
end
;
; $Id$
;
;+
; NAME:
;
;	HealpixMap
;
; PURPOSE:
;
;	Object representation of a HEALPix map.
;
; CATEGORY:
;
;	Object
;
; CALLING SEQUENCE:
;
;	hmap = obj_new( 'HealpixMap' )
;
; INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;
;	$Log$
;
;-

pro HealpixMap::from_fits_file, filename

	print, 'Reading ', filename

        fits_open,filename,fcb
        fits_read, fcb, dummy, phdr
	fits_read, fcb, tab, ehdr,exten=1,/no_pdu
        fits_close, fcb
	tbinfo, ehdr, tb_str
	
;	read_fits_s, filename, prim_st, exten_st

;	phdr = prim_st.hdr
;	ehdr = exten_st.hdr

	;self.coordsys = fxpar( ehdr, 'COORDSYS' )
	self.filename = filename
	self.ordering = strmid( fxpar(ehdr, 'ORDERING'), 0, 1 )
	self.nside = strtrim( fxpar(ehdr, 'NSIDE'),2 )
	print, self.nside, self.ordering

	self.n_maps = N_elements( tb_str.tbcol ) 

	;; get maps
	self.map = ptr_new( ptrarr(self.n_maps), /no_copy )
	for i = 0, self.n_maps-1 do begin
		(*self.map)[i] = ptr_new( tbget(tb_str, tab,i+1), /no_copy )
	endfor
	
	;; get map names
	tags = tb_str.ttype
	self.map_name = ptr_new( tags, /no_copy )

	;; get unit strings
;	unit = strarr( self.n_maps )
;	for i = 0, self.n_maps-1 do begin
;		unit_str = 'TUNIT' + int2str(i+1)
;		unit[i] = fxpar( ehdr, unit_str )
;	endfor
	unit = tb_str.tunit
	ptr_free,tb_str.tscal
	ptr_free,tb_str.tzero
	self.map_unit = ptr_new( unit )

	;; save headers
	self.phdr = ptr_new( phdr, /no_copy )
	self.ehdr = ptr_new( ehdr, /no_copy ) 
end

pro HealpixMap::get_disc_pix, pos, radius, pix, vec=vec, deg=deg

	pix = -1

	if not keyword_set(vec) then begin
		l = pos[0]
		b = pos[1]
		sky2vec, l, b, pos_vec, deg=deg
	endif else begin
		pos_vec = pos
	endelse

	getdisc_ring, self.nside, pos_vec, radius, pix_ring, deg=deg

	if self->nested() then $
		ring2nest, self.nside, pix_ring, pix $
	else $
		pix = pix_ring
end

function HealpixMap::get_header_val, name, _extra=extra

	return, strtrim( fxpar(*self.phdr, name, _extra=extra) )
end

function HealpixMap::get_map, index=index, name=name, ptr=ptr, unit=unit

	i = self->get_map_index( index=index, name=name )

	if keyword_set(unit) then $
		unit = self->get_unit( index=i ) 

	if keyword_set(ptr) then $
		return, (*self.map)[i] $
	else $
		return, *(*self.map)[i]
end

function HealpixMap::get_rotated_map, rot, deg=deg, index=index, name=name

	i = self->get_map_index( index=index, name=name )

	self->pix_to_ang, lindgen( self->npix() ), theta, phi, deg=deg
	phi = phi + rot
	self->ang_to_pix, theta, phi, pix, deg=deg

	return, (*(*self.map)[i])[pix]
end

pro HealpixMap::read_map, map, index=index, name=name

	i = self->get_map_index( index=index, name=name )
	map = *(*self.map)[i]
end

function HealpixMap::get_map_index, index=index, name=name, asd=asd

	i_map = 0

	if keyword_def(index) then begin
		i_map = index
	endif

	if keyword_def(name) then begin
		for i = 0, self.n_maps-1 do $
			if (*self.map_name)[i] eq name then i_map = i
	endif

	return, i_map
end

function HealpixMap::get_map_names

	return, *self.map_name
end

function HealpixMap::get_map_units

	return, *self.map_unit
end

function HealpixMap::get_ordering, full=full

	if keyword_set(full) then begin

		case self.ordering of
		'N'	: return, 'nest'
		'R'	: return, 'ring'
		else	: return, 'oops'
		endcase

	endif else return, self.ordering
end

function HealpixMap::get_unit, index=index, name=name

	i_map = self->get_map_index( index=index, name=name )
	
	return, (*self.map_unit)[i_map]
end

pro HealpixMap::info, nside=nside, ordering=ordering

	nside = self.nside
	ordering = self.ordering
end

function HealpixMap::map_exists, index=index, name=name

	i = self->get_map_index( index=index, name=name )
	if i eq -1 then return, 0 else return, 1
end

function HealpixMap::nested

	if self.ordering eq 'N' then return, 1 else return, 0
end

function HealpixMap::ring

	if self.ordering eq 'R' then return, 1 else return, 0
end

function HealpixMap::npix

	return, nside2npix( self.nside )
end

function HealpixMap::nside

	return, self.nside
end

function HealpixMap::pixval, pix, index=index, name=name

	map = self->get_map( index=index, name=name, /ptr )
	return, (*map)[pix]
end

pro HealpixMap::set_pixval, pix, val, index=index, name=name

	map = self->get_map( index=index, name=name, /ptr )
	(*map)[pix] = val
end

pro HealpixMap::add_pixval, pix, val, index=index, name=name

	map = self->get_map( index=index, name=name, /ptr )
	(*map)[pix] = (*map)[pix] + val
end

pro HealpixMap::set_map, index, map, name=name, unit=unit

	(*self.map)[index] = ptr_new( map )

	if keyword_set(name) then (*self.map_name)[index] = name
	if keyword_set(unit) then (*self.map_unit)[index] = unit
end

pro HealpixMap::set_map_name, index, name

	(*self.map_name)[index] = name
end

pro HealpixMap::add_to_primary_header, xhdr

	sxaddpar, xhdr, 'COMMENT', 'prim_st' ; doesn't like empty hdr
end

pro HealpixMap::write_fits_file, filename

	print, 'Writing ', filename

	map_name = self->get_map_names()
	map_unit = self->get_map_units()

	;; primary header
	info_xhdr = strarr(1)
;	sxaddpar, info_xhdr, 'COMMENT', 'prim_st' ; doesn't like empty hdr
	self->add_to_primary_header, info_xhdr
	prim_st = create_struct( 'HDR', info_xhdr )

	;; extension header
	info_xhdr = strarr(1)
	sxaddpar, info_xhdr, 'COMMENT', 'exten_st'
	for i = 0,self.n_maps-1 do begin
		istr = int2str[i+1]
		sxaddpar, info_xhdr, 'TTYPE' + istr, map_name[i]
		sxaddpar, info_xhdr, 'TUNIT' + istr, map_unit[i]
	endfor

	exten_st = create_struct( 'HDR', info_xhdr )
	map_names = self->get_map_names()
	for i = 0,self.n_maps-1 do begin
		map = self->get_map( index=i )
		exten_st = create_struct( exten_st, map_name[i], map )
	endfor

	ordering = self->get_ordering( /full )

	write_fits_sb, filename, prim_st, exten_st, ordering=ordering
end

pro HealpixMap::coadd_map, map_obj

	self->read_map, x0, index=0
	self->read_map, wt0, index=1

	map_obj->read_map, x1, index=0
	map_obj->read_map, wt1, index=1

	wt = wt0 + wt1
	iwh = where( wt gt 0 )

	x = x0*wt0 + x1*wt1
	x[iwh] = x[iwh]/wt[iwh]

	iwh = where( wt le 0, ict )
	if ict gt 0 then begin
		x[iwh] = 0.
		wt[iwh] = 0.
	endif

	iwh = lindgen(n_elements(x))
	self->set_pixval, iwh, x, index=0
	self->set_pixval, iwh, wt, index=1
end

;; coordinate conversion routines ---------------------------------------------

pro HealpixMap::pix_to_ang, pix, theta, phi, deg=deg

	if self->nested() then $
		pix2ang_nest, self.nside, pix, theta, phi $
	else $
		pix2ang_ring, self.nside, pix, theta, phi

	if keyword_set(deg) then begin
		theta = theta*!RADEG
		phi = phi*!RADEG
	endif
end

pro HealpixMap::ang_to_pix, theta, phi, pix, deg=deg

	if keyword_set(deg) then begin
		th = theta*!DTOR
		ph = phi*!DTOR
	endif else begin
		th = theta
		ph = phi
	endelse

	if self->nested() then $
		ang2pix_nest, self.nside, th, ph, pix $
	else $
		ang2pix_ring, self.nside, th, ph, pix
end

pro HealpixMap::pix_to_sky, pix, l, b, deg=deg

	self->pix_to_ang, pix, theta, phi, deg=deg

	if keyword_set(deg) then b0 = 90. else b0 = 0.5*!pi

	l = phi
	b = b0 - theta
end

pro HealpixMap::sky_to_pix, l, b, pix, deg=deg

	if keyword_set(deg) then begin
		theta = (90. - b) * (!pi/180)
		phi = l * (!pi/180)
	endif else begin
		theta = 0.5*!pi - b
		phi = l
	endelse

	iwh = where( theta lt 0. or theta gt !pi, ict )
	if ict gt 0 then theta[iwh] = 0.

	self->ang_to_pix, theta, phi, pix

	if ict gt 0 then pix[iwh] = -1
end

pro HealpixMap::pix_to_vec, pix, vec

	if self->nested() then $
		pix2vec_nest, self.nside, pix, vec $
	else $
		pix2vec_ring, self.nside, pix, vec
end

pro HealpixMap::vec_to_pix, vec, pix

	if self->nested() then $
		vec2pix_nest, self.nside, vec, pix $
	else $
		vec2pix_ring, self.nside, vec, pix
end

;; lifecycle routines ---------------------------------------------------------

function HealpixMap::init, file=file, nside=nside, ordering=ordering, n_maps=n_maps

	if not keyword_def(nside) then nside = 512
	if not keyword_def(ordering) then ordering = 'N'
	if not keyword_def(n_maps) then n_maps = 1

	if keyword_set(file) then begin
		self->from_fits_file, file
	endif else begin
		self.ordering = ordering
		self.nside = nside
		self.n_maps = n_maps
		self.map = ptr_new( ptrarr(self.n_maps), /no_copy )
		self.map_name = ptr_new( strarr(self.n_maps), /no_copy )
		self.map_unit = ptr_new( strarr(self.n_maps), /no_copy )
;		npix = nside2npix( self.nside )
;		for i = 0,n_maps-1 do begin
;			(*self.map)[i] = ptr_new( fltarr(npix), /no_copy )
;		endfor
	endelse

	return, 1
end

pro HealpixMap::cleanup

	for i = 0, self.n_maps-1 do ptr_free, (*self.map)[i]
	ptr_free, *self.map
	ptr_free, self.map

	ptr_free, self.phdr, self.ehdr
	ptr_free, self.map_name, self.map_unit
end

pro HealpixMap__define

	struct = { HealpixMap, $
		phdr		: ptr_new(), $
		ehdr		: ptr_new(), $
		filename	: '', $
		ordering	: '', $
		nside		: 0L, $
		n_maps		: 0, $
		map		: ptr_new(), $
		map_name	: ptr_new(), $
		map_unit	: ptr_new() $
	}
end


function GnomicProjection::init

	return, 1
end

pro GnomicProjection::cleanup

	ptr_free, self.euler_matrix
end

pro GnomicProjection::set, l0=l0, b0=b0, rot=rot

	if keyword_set(l0) then self.l0 = l0
	if keyword_set(b0) then self.b0 = b0
	if keyword_set(rot) then self.rot = rot

	rot_ang = [self.l0, self.b0, self.rot]
	self.do_rot = (total(abs(rot_ang)) gt 1.e-5)

	eul_mat = euler_matrix_new( self.l0, -self.b0, self.rot, /deg, /zyx )
	if ptr_valid(self.euler_matrix) then ptr_free, self.euler_matrix
	self.euler_matrix = ptr_new( eul_mat, /no_copy )
end

pro GnomicProjection::hmap_to_image, hmap, image, index=index, name=name, $
	diameter=diameter, xsize=xsize, ysize=ysize, _extra=extra
	;l0=l0, b0=b0, rot=rot, diameter=diameter, $

	;@gnomcom
	;loadsky

;	if not keyword_set(l0) then l0 = 0.
;	if not keyword_set(b0) then b0 = 0.
;	if not keyword_set(rot) then rot = 0.
;	eul_mat = euler_matrix( -self.l0, -self.b0, -self.rot, /deg, /zyx )
;	do_rot = (total(abs(rot_ang)) gt 1.e-5)

	rot_ang = [self.l0, self.b0, self.rot]

	hmap->info, nside=nside, ordering=ordering
	;map = hmap->get_map( index=index, name=name, unit=unit, /ptr )
	hmap->read_map, map, index=index, name=name
	unit = hmap->get_unit( index=index, name=name )

	; reso_arcmin is arcmin / pixel
	; diameter is in degrees
	reso = (60.*diameter)/sqrt(xsize*ysize)
	
	wmap_data2gnom, map, dummy, ordering, nside, $
		0, $			; do_conv
		self.do_rot, $		; do_rot
		'G', $			; coord_in
		'G', $			; coord_out
		*self.euler_matrix, $	; eul_mat
		image, tmax, tmin, color_bar, dx, $
		pxsize=xsize, pysize=ysize, $
		rot=rot_ang, reso_arcmin=reso, $
		_extra = extra
end

pro GnomicProjection::uv_to_lb, xpos, ypos, lon_deg, lat_deg
;+
; convert a uv position on a Gnomic map into (long, lat)
; only for scalar input
;-
	lon_deg = -1000.
	lat_deg = -1000.

;	if ( xpos ge !x.crange(0) and xpos le !x.crange(1) and $
;	     ypos ge !y.crange(0) and ypos le !y.crange(1) ) then begin

	vector = [1., -xpos, ypos] ; minus sign = astro convention
	if self.do_rot then vector = vector # (*self.euler_matrix)
	vector1 = vector
	vec2ang, vector1, lat_deg, lon_deg, /astro

		; we go from the final Gnomonic map (system coord_out) to
		; the original one (system coord_in)
		;if (do_conv) then vector = $
		;	SKYCONV(vector, inco=coord_out, outco=coord_in)
;	endif
end

pro GnomicProjection__define

	struct = { GnomicProjection, $
		l0		: 0., $
		b0		: 0., $
		rot		: 0., $
		do_rot		: 0, $
		euler_matrix	: ptr_new() $
	}
end

pro mv_cleanup, w

	widget_control, w, get_uvalue=state, /no_copy

	;; restore initial parameters
	tvlct, state.initial_ct
	!p.background = state.initial_background
end

pro mv_event_handler, event

	common mvCommon, $
		c_WidgetRoot, $
		c_MollewideProjection, $
		c_FullImageWidget, $
		c_GnomicProjection, $
		c_ZoomImageWidget, $
		c_Map_Select, $
		c_LabelL, c_LabelB, c_LabelPix, c_LabelVal, c_LabelUnit, $
		c_TextL, c_TextB, c_SliderSize, $
		c_TextFnFull, c_TextFnZoom, $
		c_zoom_select, $
		c_zoom_l, c_zoom_b, c_zoom_size_arcmin, $
		c_filename, $
		c_hmap, $
		c_hmap_index, $
		c_mouse_x0, c_mouse_y0, $
		c_mouse_x1, c_mouse_y1, $
		c_transfer, $
		c_fn_full, c_fn_zoom

	widget_control, event.id, get_uvalue=uvalue

	case uvalue of

		'button_draw' : begin

		;	mv_image_view, c_FullImageWidget, $
		;		c_MollewideProjection, c_hmap, $
		;		transfer=c_transfer
			mv_draw_full
		end

		'button_erase' : begin

			c_FullImageWidget->erase
			c_ZoomImageWidget->erase
		end

		'color_scaling' : begin

			print, event.select
		end
	endcase
end

pro mv_image_view, widget, projection, hmap, _extra=extra

	common mvCommon

	widget_control, /hourglass
	projection->hmap_to_image, hmap, image, index=c_hmap_index, _extra=extra
	widget->load_image, image
	widget->draw_image
end

pro mv_draw_full

	common mvCommon

	case c_transfer of
		'linear'	: transfer = 1
		'log'		: transfer = 2
		'histogram'	: transfer = 3
	endcase

	mv_image_view, c_FullImageWidget, c_MollewideProjection, c_hmap, $
		scale=1, transfer=transfer, size=1, coord=1
end

pro mv_draw_zoom

	common mvCommon

	c_ZoomImageWidget->info, xsize=xsize, ysize=ysize

	log = c_transfer eq 'log'
	hist_equal = c_transfer eq 'histogram'

	diameter = c_zoom_size_arcmin/60.
	mv_image_view, c_ZoomImageWidget, c_GnomicProjection, c_hmap, $
		l0=c_zoom_l, b0=c_zoom_b, rot=0., diameter=diameter, $
		xsize=xsize, ysize=ysize, log=log, hist_equal=hist_equal
end

pro mv_zoom_set_pos, l, b

	common mvCommon

	c_zoom_l = l
	widget_control, c_TextL, set_value=num2str(c_zoom_l)

	c_zoom_b = b
	widget_control, c_TextB, set_value=num2str(c_zoom_b)

	c_GnomicProjection->set, l0=l, b0=b
end

pro mv_zoom_set_size, size

	common mvCommon

	c_zoom_size_arcmin = size
	widget_control, c_SliderSize, set_value=c_zoom_size_arcmin
end

pro mv_zoom_change

	mv_draw_zoom
	mv_zoom_draw_circle
end

pro mv_zoom_update, event

	common mvCommon

	widget_control, c_TextL, get_value=l_str
	widget_control, c_TextB, get_value=b_str
	widget_control, c_SliderSize, get_value=size_arcmin

	mv_zoom_set_pos, float(l_str[0]), float(b_str[0])
	mv_zoom_set_size, size_arcmin

	mv_zoom_change
end

pro mv_draw_box, image, projection, x0, y0, x1, y1

	image->refresh
	image->window_set
	;print, x0, y0, x1, y1
	plots, [x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], /device
end

pro draw_circle, x, y, r, _extra=extra

	n = 100
	theta = 2.*!pi*findgen(n)/(n-1)
	xvec = r*cos(theta) + x
	yvec = r*sin(theta) + y

	plots, [x-2,x+2], [y,y], _extra=extra
	plots, [x,x], [y-2,y+2], _extra=extra
	plots, xvec, yvec, _extra=extra
end

pro draw_mollewide_circle, l, b, r, _extra=extra

	common mvCommon

	c_MollewideProjection->uv_circle, l, b, r, u, v, lon=lon

	;; crosshair
	;plots, [x-2,x+2], [y,y], _extra=extra
	;plots, [x,x], [y-2,y+2], _extra=extra

	x = (u + 2.)/4.
	y = (v + 1.)/2.

	PLOTS, X[0], Y[0], /normal, _extra=extra
	FOR J = 1,n_elements(x)-1 DO BEGIN
		IF ( ABS(lon[J]-lon[J-1]) LT 180. ) THEN BEGIN
			PLOTS, X[J], Y[J], /normal, /continue, _extra=extra
		ENDIF ELSE BEGIN
			PLOTS, X[J], Y[J], /normal, _extra=extra
		ENDELSE
	ENDFOR
end

pro mv_zoom_draw_circle

	common mvCommon

	c_FullImageWidget->refresh
	c_FullImageWidget->window_set

	;r0 = sqrt( (c_mouse_x0 - c_mouse_x1)^2 + (c_mouse_y0 - c_mouse_y1)^2 )
	;draw_circle, c_mouse_x0, c_mouse_y0, r0, color=255, /device

	draw_mollewide_circle, c_zoom_l, c_zoom_b, c_zoom_size_arcmin/60., color=255
end

pro mv_full_mouse_event, event

	common mvCommon

	x = event.x
	y = event.y

	c_FullImageWidget->info, xsize=xsize, ysize=ysize
	u = 4.*float(x)/float(xsize) - 2.
	v = 2.*float(y)/float(ysize) - 1.

	c_MollewideProjection->uv_to_lb, u, v, l, b
	c_hmap->sky_to_pix, l, b, pix, /deg

	pix = pix[0]
	l = l[0]
	b = b[0]

	badpix = pix lt 0

	case event.type of

		0 : begin	;; press

			if not badpix then begin 
				c_zoom_select = 1
				mv_zoom_set_pos, l, b
				c_mouse_x0 = x
				c_mouse_y0 = y
			endif
		end

		1 : begin	;; release

			c_zoom_select = 0
			mv_draw_zoom
			mv_zoom_draw_circle
		end

		2 : begin	;; motion

			if c_zoom_select then begin

				c_mouse_x1 = x
				c_mouse_y1 = y

				;mv_draw_box, c_FullImageWidget, $
				;	c_MollewideProjection, $
				;	c_mouse_x0, c_mouse_y0, $
				;	c_mouse_x1, c_mouse_y1

			;	r0 = sqrt( (c_mouse_x0 - c_mouse_x1)^2 + $
			;		(c_mouse_y0 - c_mouse_y1)^2 )
			;
			;	mv_draw_circle, c_FullImageWidget, $
			;		c_mouse_x0, c_mouse_y0, r0

				u0 = 4.*float(c_mouse_x0)/float(xsize) - 2.
				v0 = 2.*float(c_mouse_y0)/float(ysize) - 1.

				r = 90.*sqrt( (u - u0)^2 + (v - v0)^2 )

				mv_zoom_draw_circle
				mv_zoom_set_size, r*60
			endif

			if not badpix then begin
				val = c_hmap->pixval( pix, index=c_hmap_index )
				l_str = num2str(l)
				b_str = num2str(b)
				pix_str = int2str(pix)
				val_str = num2str(val)
			endif else begin
				l_str = 'n/a'
				b_str = 'n/a'
				pix_str = 'n/a'
				val_str = 'n/a'
			endelse

			widget_control, c_LabelL, set_value=l_str
			widget_control, c_LabelB, set_value=b_str
			widget_control, c_LabelPix, set_value=pix_str
			widget_control, c_LabelVal, set_value=val_str
		end
	endcase
end

pro mv_zoom_mouse_event, event

	common mvCommon

	x = event.x
	y = event.y

	c_ZoomImageWidget->info, xsize=xsize, ysize=ysize
	reso = (!pi/180.)*(c_zoom_size_arcmin/60.)/sqrt(xsize*ysize)
	dx = x - 0.5*xsize
	dy = y - 0.5*ysize

	u = reso*dx
	v = reso*dy

	c_GnomicProjection->uv_to_lb, u, v, l, b
	c_hmap->sky_to_pix, l, b, pix, /deg

	pix = pix[0]
	l = l[0]
	b = b[0]

	case event.type of

		0 : begin	;; press

			mv_zoom_set_pos, l, b
			mv_zoom_change
		end

		1 : begin	;; release

			;print, 'release'
		end

		2 : begin	;; motion

			val = c_hmap->pixval( pix, index=c_hmap_index )
			l_str = num2str(l)
			b_str = num2str(b)
			pix_str = int2str(pix)
			val_str = num2str(val)

			widget_control, c_LabelL, set_value=l_str
			widget_control, c_LabelB, set_value=b_str
			widget_control, c_LabelPix, set_value=pix_str
			widget_control, c_LabelVal, set_value=val_str
		end
	endcase
end

pro mv_full_write_file, event

	common mvCommon

	widget_control, c_TextFnFull, get_value=c_fn_full
	c_FullImageWidget->write_png, c_fn_full
end

pro mv_zoom_write_file, event

	common mvCommon

	widget_control, c_TextFnZoom, get_value=c_fn_zoom
	c_ZoomImageWidget->write_png, c_fn_zoom
end

pro mv_redraw

	common mvCommon

	mv_draw_full
	mv_draw_zoom
	mv_zoom_draw_circle
end

function mv_color_scale, event

	common mvCommon

	c_transfer = event.value
	mv_redraw

	return, 1
end

pro mv_color_table, event

	common mvCommon

	notify1 = {object:c_FullImageWidget,method:'draw_image'}
	notify2 = {object:c_ZoomImageWidget,method:'draw_image'}
	xcolors, notifyobj=[notify1,notify2], notifypro='mv_zoom_draw_circle', $
		group=event.top
end

pro mv_filename, event

       common mvCommon 
      
	obj_destroy, c_hmap
	c_filename = event.filename
        c_hmap = obj_new('HealpixMap', file=event.filename) 
        found_input = 1
	widget_control, c_map_select , $
		set_value=c_hmap->get_map_names()
	mv_redraw
end

pro mv_change_map, event

	common mvCommon

	c_hmap_index = event.index
	widget_control, c_LabelUnit, $
		set_value=c_hmap->get_unit( index=c_hmap_index )
	mv_redraw
end

pro mv_quit, event

	common mvCommon

;	print, 'Current heap...'
;	help, /heap

;	print, 'Destroying...'
	obj_destroy, c_MollewideProjection
	obj_destroy, c_GnomicProjection
	obj_destroy, c_hmap
	obj_destroy, c_FullImageWidget
	obj_destroy, c_ZoomImageWidget
	widget_control, c_WidgetRoot, /destroy

;	print, 'Current heap...'
;	help, /heap
;	print, 'Checking for leaks...'
;	heap_gc, /verbose
end

pro mapview, file=file, array=array, ring=ring

	common mvCommon

	program_name = 'MAPVIEW v1.3 (07 Jul 2004)'
	print, program_name + ' -- WMAP project HEALPix map viewer'

	;; only one instance can run at any one time due to common block
	if xregistered('mv') then begin
		print, 'Can only run one mapview per IDL session -- quitting'
		return
	endif

	;if not keyword_set(file) then begin
	;	print, 'No file specified -- quitting'
	;	return
	;endif

	found_input = 0

	if keyword_set(file) then begin
		if not file_test(file) then file = ''
		c_filename = file
                if file NE '' then begin
		input_hmap = obj_new( 'HealpixMap', file=file )
		found_input = 1
		fdecomp,file,disk,dir,fname
		cdir = disk + dir
		if cdir EQ '' then cd,current=cdir
		endif
	endif else  begin
	        fname = ''
		cd,current=cdir
		no_array = N_elements(array) EQ 0
	        if no_array then begin
		     c_filename = '' 
		     array = randomu(seed,786432,2)
		endif else c_filename = '(command line)'
		dims = size( array, /dim )
		npix = dims[0]
		nside = npix2nside( npix )
		ordering = keyword_set(ring) ? 'R' : 'N'
		n_maps = 1
		if n_elements(dims) eq 2 then n_maps = dims[1] 
		input_hmap = obj_new( 'HealpixMap', $
			nside=nside, ordering=ordering, n_maps=n_maps )

                if no_array then $
		      print,'Waiting for file name to be input' $
		else begin	
		print, 'Reading data from command line:'
		print, '  nside = ', nside
		print, '  n_maps = ', n_maps
                endelse
			
		for i = 0,n_maps-1 do begin
			name = 'array' + int2str(i)
			input_hmap->set_map, i, array[*,i], name=name
		endfor
		found_input = 1
      endelse


	;; save initial info --------------------------------------------------

	swin = !d.window

	initial_background = !p.background
	tvlct, ir, ig, ib, /get
	initial_ct = [[ir],[ig],[ib]]

	;; initialize common block --------------------------------------------

	device, decomp=0
	loadct, 13

  	c_zoom_b = 0.
	c_zoom_l = 0.
	c_zoom_size_arcmin = 10*60
	c_zoom_select = 0
	c_transfer = 'linear'
	if found_input EQ 1 then c_hmap = input_hmap
	c_hmap_index = 0
	c_fn_full = 'map.png'
	c_fn_zoom = 'zoom.png'

	;; layout -------------------------------------------------------------

	c_WidgetRoot = widget_base( /row, $
		title='mv: ' + c_filename, xoffset=10, yoffset=10 )
	widget_control, c_WidgetRoot, /managed

	;; --- controls ---

	control_base = widget_base( c_WidgetRoot, /col )

	tmp = widget_label( control_base, value=program_name )

	;; map selection
      
        file_base = widget_label(control_base, /frame, value= 'FILE')
          file1 = fsc_fileselect(control_base, event_pro='mv_filename', /nomax, $
	       /mustexist, filter = '*.fits', /read, filename=fname,dir=cdir )
	     
	map_base = widget_base( control_base, /col, /frame )

	tmp = widget_label( map_base, value='Select map:', /align_left )
	
	map_names = c_hmap->get_map_names()
	C_map_select = widget_droplist( map_base, value=map_names, $
	           event_pro='mv_change_map')

	;; --- cursor info ---

	cursor_base = widget_base( control_base, /col, /frame )

	tmp = widget_label( cursor_base, value='Cursor info:', /align_left )

	cursor_l_base = widget_base( cursor_base, /row )
	tmp = widget_label( cursor_l_base, value='l =', /align_right )
	c_LabelL = widget_label( cursor_l_base )
	widget_control, c_LabelL, set_value='000000000000'

	cursor_b_base = widget_base( cursor_base, /row )
	tmp = widget_label( cursor_b_base, value='b =', /align_right )
	c_LabelB = widget_label( cursor_b_base )
	widget_control, c_LabelB, set_value='000000000000'

	cursor_pix_base = widget_base( cursor_base, /row )
	tmp = widget_label( cursor_pix_base, value='pix =', /align_right )
	c_LabelPix = widget_label( cursor_pix_base )
	widget_control, c_LabelPix, set_value='00000000000'

	cursor_val_base = widget_base( cursor_base, /row )
	tmp = widget_label( cursor_val_base, value='val =', /align_right )
	c_LabelVal = widget_label( cursor_val_base )
	widget_control, c_LabelVal, set_value='000000000000'
	c_LabelUnit = widget_label( cursor_val_base, value='xxxxx' )
	widget_control, c_LabelUnit, $
		set_value=c_hmap->get_unit( index=c_hmap_index )

	;; --- color controls ---

	color_base = widget_base( control_base, /col, /frame )

	scale_names = ['linear','log','histogram']
	color_scale_cw = cw_bgroup( color_base, scale_names, $
		button_uvalue=scale_names, event_func='mv_color_scale', $
		/exclusive, label_top='Color scaling:', $
		/no_release, set_value=0 )

	color_button = widget_button( color_base, $
		value='Color Table', event_pro='mv_color_table' )

	;; --- zoom controls ---

	zoom_base = widget_base( control_base, /col, /frame )

	tmp = widget_label( zoom_base, value='Zoom region:', /align_left )

	text_l_base = widget_base( zoom_base, /row )
	tmp = widget_label( text_l_base, value='Center l [deg]:' )
	c_TextL = widget_text( text_l_base, /editable, xsize=8, $
		event_pro='mv_zoom_update', value=num2str(c_zoom_l) )

	text_b_base = widget_base( zoom_base, /row )
	tmp = widget_label( text_b_base, value='Center b [deg]:' )
	c_TextB = widget_text( text_b_base, /editable, xsize=8, $
		event_pro='mv_zoom_update', value=num2str(c_zoom_b) )

	c_SliderSize = widget_slider( zoom_base, $
		title='Size [arcmin]', minimum=1, maximum=4200, $
		event_pro='mv_zoom_update', value=c_zoom_size_arcmin )

	zoom_button = widget_button( zoom_base, value='Draw', $
		event_pro='mv_zoom_update' )

	;; --- write to file ---

	write_base = widget_base( control_base, /col, /frame )

	tmp = widget_label( write_base, value='Write img to file:', /align_left )

	text_fn_full_base = widget_base( write_base, /row )
	tmp = widget_label( text_fn_full_base, value=' Map:' )
	;text_fn_full_subbase = widget_base( text_fn_full_base, /row )
	c_TextFnFull = widget_text( text_fn_full_base, /editable, xsize=15, $
		event_pro='mv_full_write_file', value=c_fn_full )
	write_full_button = widget_button( text_fn_full_base, value='x', $
		event_pro='mv_full_write_file' )

	text_fn_zoom_base = widget_base( write_base, /row )
	tmp = widget_label( text_fn_zoom_base, value='Zoom:' )
	c_TextFnZoom = widget_text( text_fn_zoom_base, /editable, xsize=15, $
		event_pro='mv_zoom_write_file', value=c_fn_zoom )
	write_zoom_button = widget_button( text_fn_zoom_base, value='x', $
		event_pro='mv_zoom_write_file' )

	;; --- buttons ---

	button_base = widget_base( control_base, /row )

	quit_button = widget_button( button_base, $
		value='Quit', event_pro='mv_quit' )

	;; --- drawing regions ---

	draw_base = widget_base( c_WidgetRoot, /col )

	tmp = widget_label( draw_base, value=c_filename )

	c_MollewideProjection = obj_new( 'MollewideProjection' )
	c_FullImageWidget = obj_new( 'ImageWidget', draw_base, 512, 256, $
		/frame, /button_events, /motion_events, $
		event_pro='mv_full_mouse_event' )

	c_GnomicProjection = obj_new( 'GnomicProjection' )
	c_GnomicProjection->set, l0=c_zoom_l, b0=c_zoom_b

	c_ZoomImageWidget = obj_new( 'ImageWidget', draw_base, 512, 512, $
		/frame, /button_events, /motion_events, $
		event_pro='mv_zoom_mouse_event' )

	;; go -----------------------------------------------------------------

	state = { $
		initial_background	: initial_background, $
		initial_ct		: initial_ct }

	;; show
	widget_control, c_WidgetRoot, set_uvalue=state, /realize
	mv_draw_full

	;; restore initial window
	wset, swin

	xmanager, 'mv', c_WidgetRoot, $
		cleanup		= 'mv_cleanup', $
		event_handler	= 'mv_event_handler', $
		/no_block
end
