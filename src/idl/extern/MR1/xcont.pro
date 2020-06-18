;PRO XCONT, Imag_In, Window=Window
;+ 
; NAME: 
;      XCONT
;
; PURPOSE: 
;       widget for plotting contours of an image. The interface allows
;       to select the contours, to annotate the image, print the result
;       to a poscript file, ...
;
; CALLING SEQUENCE: 
;   XCONT, Imag_In, Window=Window
;
; INPUTS: 
;   Imag_In -- 2 d array : input array
;
; KEYED INPUTS: 
;   Window: -- integer : IDL window number 
;           define the window in which we plot the contours
;
; EXAMPLE: 
;       XCONT, Image, window=7
;
; MODIFICATION HISTORY: 
;    9-Feb-1995 JL Starck  
;-


FUNCTION cw_label_choix_event,ev

base=ev.handler
WIDGET_CONTROL, WIDGET_INFO(base, /CHILD), $
		GET_UVALUE=state, /NO_COPY

if state.value eq state.value1 then state.value = state.value2 $
else state.value = state.value1
WIDGET_CONTROL,state.id,set_value=state.value

ret={ ID:base, TOP:ev.top, HANDLER:0L, VALUE:state.value, IDBUT:state.id}

WIDGET_CONTROL, WIDGET_INFO(base, /CHILD), $
		SET_UVALUE=state,/NO_COPY
return,ret

end


FUNCTION cw_label_choix,parent, $
	FRAME = frame, $
	UVALUE = uval, $
	VALUE1 = val1, $
	VALUE2 = val2, $
	TITLE = title

; P.G., dernieres modifs le 19/05/94
; ce compound widget comprend un label (valeur de TITLE), 
; et un boutton qui prend alternativement les deux valeurs 
; definies dans VALUE1 et VALUE2 lorsqu'on clique dessus.
; la valeur est initialisee a VALUE1
; les seules options acceptees sont UVALUE et FRAME.
; l'evenement retourne par ce widget est du type :
; { ID:, TOP:, HANDLER:, VALUE:, IDBUT:}
; ID, TOP et HANDLER sont comme pour tous les widgets
; VALUE contient la valeur actuelle du boutton
; IDBUT contient l'IDentificateur du boutton


;parametres par defaut :
IF NOT KEYWORD_SET(frame) THEN frame = 0
IF NOT KEYWORD_SET(uval)  THEN uval = 0
IF NOT KEYWORD_SET(val1)   THEN val1= '1'
IF NOT KEYWORD_SET(val2)   THEN val2= '2'
IF NOT KEYWORD_SET(title) THEN title = ''


base=WIDGET_BASE(parent,FRAME=frame,/ROW,UVALUE=uval,$
		EVENT_FUNC='cw_label_choix_event' )
label=WIDGET_LABEL(base,VALUE=title)
butt=WIDGET_BUTTON(base,VALUE=val1)

state={ id:butt, value:val1, value1:val1, value2:val2}
WIDGET_CONTROL, WIDGET_INFO(base, /CHILD), $
		SET_UVALUE=state,/NO_COPY

return, base
end




pro common_xcontour
common share1,image,nbniv,label,taillechar,man_auto,niveaux,max,min
common share_id,auto,manu	,manniv,idniv
common share2,flag,couleurs,coul,ech
common share3, oldwin, contwin
common share4,FIELDMIN,FIELDMAX,MINREAL,MAXREAL,Imag_In
end

pro regniv

common share1
common share_id
common share2
common share4

manniv=WIDGET_BASE(/column,title='Manual Levels',/frame, $
	xsize=255,yoffset=380,xoffset=92)
idniv=lonarr(nbniv)

minv = max([0,min])
maxv = max([1,max])
if max LE minv then maxv = minv+1
for i=0,nbniv-1 do begin $
	Val = max([niveaux(i+nbniv),1])
	if Val LT minv then Val = minv
	if Val GT maxv then Val = maxv

	idniv(i)=WIDGET_SLIDER(manniv,min=minv, max=maxv,value=Val, $
		               uvalue=strtrim(i,1),/drag,/frame) & end
min = minv
max = maxv
if min LT MinReal then Min = MinReal
if max GT MaxReal then Max = MaxReal


;case 1 of
;
;(label eq 1):contour,image,c_charsize=taillechar,C_LINESTYLE=(niveaux lt 0.), $
;			      nlevels=2*nbniv,levels=niveaux,/follow,$
;			      c_labels=replicate(1,2*nbniv),$
;			      c_colors=couleurs
;(label eq 0) : contour,image,nlevels=2*nbniv,levels=niveaux,$
;			C_LINESTYLE=(niveaux lt 0.), c_colors=couleurs 
;endcase

WIDGET_CONTROL,manniv,/REALIZE
XMANAGER,'regniv',manniv

end


pro regniv_event,ev

common share1
common share_id
common share2
common share4

WIDGET_CONTROL,ev.id,get_uvalue=action
Value = ev.value
if Value LT MinReal then Value = MinReal
if Value GT MaxReal then Value = MaxReal
case action of

'0': begin niveaux(nbniv)=Value & niveaux(nbniv-1)=-Value & end
'1': begin niveaux(nbniv+1)=Value & niveaux(nbniv-2)=-Value & end
'2': begin niveaux(nbniv+2)=Value & niveaux(nbniv-3)=-Value & end
'3': begin niveaux(nbniv+3)=Value & niveaux(nbniv-4)=-Value & end
'4': begin niveaux(nbniv+4)=Value & niveaux(nbniv-5)=-Value & end
'5': begin niveaux(nbniv+5)=Value & niveaux(nbniv-6)=-Value & end
'6': begin niveaux(nbniv+6)=Value & niveaux(nbniv-7)=-Value & end
'7': begin niveaux(nbniv+7)=Value & niveaux(nbniv-8)=-Value & end
'8': begin niveaux(nbniv+8)=Value & niveaux(nbniv-9)=-Value & end
'9': begin niveaux(nbniv+9)=Value & niveaux(nbniv-10)=-Value & end
'10': begin niveaux(nbniv+10)=Value & niveaux(nbniv-11)=-Value & end
'11': begin niveaux(nbniv+11)=Value & niveaux(nbniv-12)=-Value & end
'12': begin niveaux(nbniv+12)=Value & niveaux(nbniv-13)=-Value & end
'13': begin niveaux(nbniv+13)=Value & niveaux(nbniv-14)=-Value & end
'14': begin niveaux(nbniv+14)=Value & niveaux(nbniv-15)=-Value & end


endcase

if ev.value lt fix(action+1) then $
BEGIN
   niveaux( fix(action)+nbniv)=fix(action+1) 
   niveaux(nbniv-fix(action)-1)=-fix(action+1)
   Val = niveaux(nbniv+fix(action))
   if Val LT 1 then Val = 1
print, Val
   WIDGET_CONTROL,ev.id,set_value=Val 
END
if ev.value gt 1-nbniv+fix(action)+max then $
BEGIN
   niveaux(fix(action)+nbniv)=1-nbniv+fix(action)+max 
   niveaux(-1-fix(action)+nbniv)=-(1-nbniv+fix(action)+max) 
   Val = niveaux(fix(action)+nbniv)
   if Val LT 1 then Val = 1
   WIDGET_CONTROL,ev.id,set_value=Val 
end

; cas ou on pousse d'autres sliders par le haut :
for i=fix(action),nbniv-2 do $
BEGIN
   if niveaux(nbniv+i) ge niveaux(nbniv+i+1) then $
   BEGIN  
     niveaux(nbniv+i+1)=niveaux(nbniv+i)+1 
     niveaux(nbniv-i-2)=-niveaux(nbniv+i)-1 
     Val = niveaux(nbniv+i+1)
     if Val LT 1 then Val = 1
     WIDGET_CONTROL,idniv(i+1),set_value=Val 
   END
END

; cas ou on pousse d'autres sliders par le bas :
for i=fix(action),1,-1 do $
BEGIN
   if niveaux(i+nbniv) le niveaux(i-1+nbniv) then $
   BEGIN
      niveaux(nbniv+i-1)=niveaux(nbniv+i)-1 
      niveaux(nbniv-i)=-niveaux(nbniv+i)+1 
      Val = niveaux(nbniv+i-1)
      if Val LT 1 then Val = 1
      WIDGET_CONTROL,idniv(i-1),set_value=Val  
   END
END

end

; =========================================================


pro cwcontour_event,ev

common share1
common share_id
common share2
common share3
common share4

WIDGET_CONTROL,ev.id,get_uvalue=action

oldwin = !d.window
wset, contwin

case action of
'Quit'	: begin WIDGET_CONTROL, ev.top, /DESTROY &$
	  if man_auto eq 1 then WIDGET_CONTROL,manniv,/DESTROY 
	  wdelete
	  oldwin=-1
	  end
'on_off': begin  label=(label xor 1) &$
	  if (label eq 0) then WIDGET_CONTROL,ev.id,set_value='OFF' $
	  else WIDGET_CONTROL,ev.id,set_value='ON'& end
	
'slidnb' : if ev.value ne nbniv then $
	   begin nbniv=ev.value &$
	   if man_auto eq 1 then WIDGET_CONTROL,manniv,/DESTROY &$
	   if ech eq 'Log ' then niveaux=[-max(image)*$
		replicate(0.5,nbniv)^findgen(nbniv),max(image)*reverse($
		replicate(0.5,nbniv)^findgen(nbniv))]$
	   else $
             niveaux=[-reverse((findgen(nbniv)/(nbniv+1)+1./(nbniv+1)) $
                     *(max-min)+min),$
		     (findgen(nbniv)/(nbniv+1)+1./(nbniv+1))*(max-min)+min] & $
	  if coul eq 0 then couleurs=replicate(255,nbniv) $
	  else couleurs=findgen(nbniv)/(nbniv-1)*200+55 &$
	  if man_auto eq 1 then begin regniv & WIDGET_CONTROL,/HOURGLASS $
          & flag=0 &$
	  endif & end $
	  else flag=0

'slidtaille': taillechar=ev.value
'Auto'	: if (man_auto eq 0) then WIDGET_CONTROL,ev.id,set_button=1 $
	  else begin man_auto=0 &$
	                  WIDGET_CONTROL,manu,set_button=0 &$
		  WIDGET_CONTROL,manniv,/DESTROY & end
'Manu'	:if (man_auto eq 1) then WIDGET_CONTROL,ev.id,set_button=1 $
	 else begin man_auto=1 &$
		WIDGET_CONTROL,auto,set_button=0 &$
		regniv & flag=0 & end
'Color'	: begin  coul=(coul xor 1) & $
	  if (coul eq 0) then begin WIDGET_CONTROL,ev.id, $
                                    set_value='Couleurs OFF' &$
			           couleurs=replicate(255,nbniv) & end $
	  else begin WIDGET_CONTROL,ev.id,set_value='Couleurs ON'&$
		  couleurs=findgen(nbniv)/(nbniv-1)*200+55 & end & end
'Scale' : begin 
	    ech=ev.value
	    case ech of
	     'Log ': niveaux=[-max(image)*$
		replicate(0.5,nbniv)^findgen(nbniv),max(image)*reverse($
		replicate(0.5,nbniv)^findgen(nbniv))]
             '  Linear   ':    levels=[-reverse((findgen(nbniv)/(nbniv+1)+1./ $ 
                        (nbniv+1))*(max-min)+min), $
	                (findgen(nbniv)/(nbniv+1)+1./(nbniv+1))*(max-min)+min] 
            endcase
	    if man_auto eq 1 then begin $
                    WIDGET_CONTROL,manniv,/DESTROY & regniv &  end
	    flag=0
	    end
'Apply' : begin
            case 1 of
             (man_auto eq 0 and label eq 1) : $
                       contour,image,c_charsize=taillechar, $
			       nlevels=nbniv,/follow,$
			      c_labels=replicate(1,nbniv),$
			      c_colors=couleurs
             (man_auto eq 0 and label eq 0) : $
                       contour,image,nlevels=nbniv,c_colors=couleurs
             (man_auto eq 1 and label eq 1) : $
                       contour,image,c_charsize=taillechar, $
			      nlevels=2*nbniv,levels=niveaux,/follow,$
			      c_labels=replicate(1,2*nbniv),$
			      c_colors=couleurs,C_LINESTYLE=(niveaux lt 0.)
             (man_auto eq 1 and label eq 0) : $
                        contour,image,nlevels=2*nbniv,levels=niveaux,$
				c_colors=couleurs,C_LINESTYLE=(niveaux lt 0.)
            endcase 
         end
'Print' : begin
                        oldwin = !d.window
            xdump 
;           set_plot, 'PS' 
;           device, filename='idl.ps'
;            case 1 of
;             (man_auto eq 0 and label eq 1) : $
;                       contour,image,c_charsize=taillechar, $
;			       nlevels=nbniv,/follow,$
;			      c_labels=replicate(1,nbniv),$
;			      c_colors=couleurs
;             (man_auto eq 0 and label eq 0) : $
;                       contour,image,nlevels=nbniv,c_colors=couleurs
;             (man_auto eq 1 and label eq 1) : $
;                       contour,image,c_charsize=taillechar, $
;			      nlevels=2*nbniv,levels=niveaux,/follow,$
;			      c_labels=replicate(1,2*nbniv),$
;			      c_colors=couleurs,C_LINESTYLE=(niveaux lt 0.)
;             (man_auto eq 1 and label eq 0) : $
;                        contour,image,nlevels=2*nbniv,levels=niveaux,$
;				c_colors=couleurs,C_LINESTYLE=(niveaux lt 0.)
;            endcase 
;            device, /close
;            set_plot, 'X'
         end
'Annotate': annotate
'FIELDMIN': BEGIN
            WIDGET_CONTROL, FIELDMIN, GET_VALUE=min
            if min LT MinReal then Min = MinReal
            WIDGET_CONTROL, FIELDMIN, SET_VALUE=Min
            niveaux=[-reverse((findgen(nbniv)/(nbniv+1)+1./ $
                    (nbniv+1))*(max-min)+min), $
	            (findgen(nbniv)/(nbniv+1)+1./(nbniv+1))*(max-min)+min] 
            image = Imag_In > min
            END
'FIELDMAX': BEGIN
            WIDGET_CONTROL, FIELDMAX, GET_VALUE=max
            if max GT MaxReal then Max = MaxReal
            WIDGET_CONTROL, FIELDMAX, SET_VALUE=Max
            niveaux=[-reverse((findgen(nbniv)/(nbniv+1)+1./ $
                    (nbniv+1))*(max-min)+min), $
	            (findgen(nbniv)/(nbniv+1)+1./(nbniv+1))*(max-min)+min]
            image = Imag_In < max
            END
endcase
   if oldwin GE 0 then wset, oldwin
end

pro xcont,Data, window=window

IF N_PARAMS() LT 1 THEN BEGIN
        PRINT, 'CALLING SEQUENCE: ', 'xcontour, Imag_In, Window=Window'
        return
        END
vsize = size(Data)
if vsize(0) NE 2 THEN BEGIN
        PRINT, 'Error in xcontour: First parameter must be an image'
        return
        END

common share1
common share_id
common share2
common share3
common share4

Imag_In = Data   
image=Imag_In 
MinReal = min(image)
MaxReal = max(image)

nbniv=6 & label=1 & taillechar=0.8 & min=min(image) & max=max(image)
man_auto=0 & flag=1
niveaux=fltarr(30)
niveaux=[-reverse((findgen(nbniv)/(nbniv+1)+1./(nbniv+1))*(max-min)+min), $
	(findgen(nbniv)/(nbniv+1)+1./(nbniv+1))*(max-min)+min] 
couleurs=findgen(nbniv)/(nbniv-1)*200+55
coul=1
ech='  Linear   '

base=WIDGET_BASE(/column,TITLE='Contours Parameters')
BASE_1_1 = WIDGET_BASE(base,ROW=1, MAP=1, UVALUE='BASE_1_1')
FieldMinVal = [ string(min) ]
FieldMaxVal = [ string(max) ]
FIELDMIN = CW_FIELD(BASE_1_1,VALUE=FieldMinVal, ROW=1, INTEGER=1, $
      TITLE='Min Value', UVALUE='FIELDMIN', XSIZE=4, /RETURN_EVENTS)
FIELDMAX = CW_FIELD(BASE_1_1,VALUE=FieldMaxVal, ROW=1, INTEGER=1, $
      TITLE='Max Value', UVALUE='FIELDMAX', XSIZE=4, /RETURN_EVENTS)

slidnbniv=WIDGET_SLIDER(base,/drag,xsize=150,maximum=15,minimum=2,value=6,$
	title='Number of contours ',/frame,uvalue='slidnb')

labels=WIDGET_base(base,/row,/frame)
lablab=WIDGET_LABEL(labels,value='Labels')
labtools=WIDGET_base(labels,/column)
togglebase = WIDGET_BASE(labtools, /COLUMN, /FRAME, /NONEXCLUSIVE)
onbutton = WIDGET_BUTTON(togglebase,VALUE='ON', UVALUE='on_off')
						
lab_taille=CW_fslider(labtools,minimum=0.5,maximum=2.5,value=0.8,$
		      title='Size',/frame,/drag,uvalue='slidtaille')

nivbase=WIDGET_base(base,/row,/frame)
nivlab=WIDGET_LABEL(nivbase,value='Levels')
toggniv=WIDGET_BASE(nivbase, /COLUMN, /FRAME)
autobase1=WIDGET_BASE(toggniv,/ROW)
autobase=WIDGET_BASE(autobase1,/NONEXCLUSIVE)
auto = WIDGET_BUTTON(autobase, $	
			VALUE='Automatic', $	
			UVALUE='Auto')
		manbase=WIDGET_BASE(toggniv,/ROW,/FRAME)
		manbasbutt=WIDGET_BASE(manbase,/NONEXCLUSIVE)
		manu= WIDGET_BUTTON(manbasbutt, VALUE='Manual', UVALUE='Manu')
		echmanu=cw_label_choix(manbase,title='Scale',VALUE2='Log ', VALUE1='  Linear   ',UVALUE='Scale')
				
coulbase=WIDGET_BASE(base,/row,/frame,/exclusive)
	coulbutt=WIDGET_BUTTON(coulbase, $	
			VALUE='Colors ON', $	
			UVALUE='Color')

apply = WIDGET_BUTTON(base,value='Apply',uvalue='Apply')
apply = WIDGET_BUTTON(base,value='Annotate',uvalue='Annotate')
apply = WIDGET_BUTTON(base,value='Print',uvalue='Print')
quit = WIDGET_BUTTON(base,value='Quit',uvalue='Quit')

if keyword_set(window) then win = window $
else win = !d.window+1
if win LT 32 then BEGIN
        t = 'ISOPHOT  '+ strcompress(string(win),/REMOVE_ALL)
        window, win, title=t
        END
contour,image,c_charsize=taillechar, nlevels=nbniv,$
	c_labels=replicate(1,nbniv),/follow,c_colors=couleurs
contwin = !d.window

WIDGET_CONTROL,base,/REALIZE
WIDGET_CONTROL,onbutton,set_button=1
WIDGET_CONTROL,auto,set_button=1
WIDGET_CONTROL,coulbutt,set_button=1

if MinReal LT 0 or (MaxReal - MinReal) LT 2 then $
     WIDGET_CONTROL,manu,sensitive=0 
XMANAGER ,'cwcontour',base
end

