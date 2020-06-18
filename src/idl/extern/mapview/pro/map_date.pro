;+
;  NAME:   
;      MAP_DATE
;  PURPOSE: 
;      Widget to convert the various MAP time and date formats, without
;      having to worry about calling sequences and names of routines.
; CALLING SEQUENCE:
;       MAP_DATE
; USAGE:    
;   Operation is pretty much self-explanatory.  Enter the values you
;   want to change into the appropriate text field(s) and hit the 
;   Enter key.  The value will propagate through in a reasonable way,
;   and a self-consistent final set of values is displayed.
;
;   The widget comes up showing the current time, which is updated
;   every second.  There are 3 buttons:  a QUIT button, a HELP button,
;   and the button to toggle between update and calculation modes.
;
;  PROCEDURES USED:
;       DAYOFYEAR, GMT2JUL(),GMT2TS, JUL2DATE, LOAD_MAP_PARAMS, TS2GMT, TS2JUL()
;  MODIFICATION HISTORY:
;       Written by Robert S. Hill, Sci Sys & Appl Inc, 25 Oct 2002
;       Added time zone droplist.  RSH, SSAI, 30 Oct 2002
;       Use systime(/JUlIAN) instead of now2jul()   W. Landsman 4 Feb 2003
;-


pro map_date_help_event, event
widget_control, event.id, get_value=val
v = strtrim(val[0],2)
if v eq 'Dismiss' then widget_control, event.top, /destroy
return
end

pro map_date_help, ev

widget_control, ev.top, get_uvalue=info

help_text = [ $
'MAP_DATE HELP', $
' ',  $
'   This widget computes and displays dates and times in the various', $
'   systems used internally by the MAP project.  ',  $
' ',  $
'   The widget has two modes:  continuous update, and calculator. ', $
' ',  $
'   The default mode is to update to the current time once per second.',  $
'   To enter your own data for computations, you need to stop the update,', $
'   which you can do using the button "Press for Calculator Mode." ',  $
' ',  $
'   Click on text fields to enter values in them.  When you hit the',  $
'   ENTER key on your keyboard, the changes will propagate through all',  $
'   the quantities in a reasonable way.',  $
' ',  $
'   Self-contradictory data will not mess up the computation, because ',  $
'   the program will just use the last type of date and time that ',  $
'   you altered.',  $
' ',  $
'   Because continuous update mode is based on system (local) time, the ',  $
'   time zone needs to be selected with the droplist at the upper right.']

if widget_info((*info).help, /valid) then begin
    
    widget_control, (*info).help, /show

endif else begin
    base = widget_base(/col, $
                        group_leader=(*info).tlb, title='MAP_DATE HELP')

    row1 = widget_base(base, /row)
    row2 = widget_base(base, /row)

    quitbutton = widget_button(row1, xsize=100, $
                    value='Dismiss')

    textdisplay = widget_text(row2, xsize=72, ysize=22, $
                               value=help_text, font='7x13')

    widget_control, base, /realize

    (*info).help = base

    xmanager, 'map_date_help', base, $
        event='map_date_help_event', /no_block
endelse

return
end

pro map_date_set_all, info, rjd, gmt, ts

widget_control, (*info).gmtyear_val, set_value=strmid(gmt,0,4)
widget_control, (*info).gmtday_val, set_value=strmid(gmt,4,3)
widget_control, (*info).gmthr_val, set_value=strmid(gmt,7,2)
widget_control, (*info).gmtmn_val, set_value=strmid(gmt,9,2)
scstr = strmid(gmt,11,2) + '.' + strmid(gmt,13,4)
widget_control, (*info).gmtsc_val, set_value=scstr

jul2date, rjd, yr, mon, day, hr, mn, sc
mon_text = (*info).months[mon-1]
widget_control, (*info).gmtmonth_val, set_value=mon_text
widget_control, (*info).gmtdom_val,  $
    set_value=string(round(day), format='(I3)')

widget_control, (*info).tsday_val,  $
    set_value=string(round(ts[0]), format='(I6)')
widget_control, (*info).tstms_val,  $
    set_value=string(round(ts[1]), format='(I10.9)')

widget_control, (*info).rjd_val,  $
    set_value=string(rjd, format='(F15.9)')

launch_rjd = gmt2jul((*info).params.true_launch_gmt)
met = rjd - launch_rjd
absmet = abs(met)
metd = floor(absmet)
meth = floor((absmet-metd)*24)
metm = floor((absmet-metd-meth/24.0d0)*24.0d0*60.0d0)
mets = (absmet-metd-meth/24.0d0-metm/(24.0d0*60.0d0))*86400.0d0
metdstr = string(metd, format='(I4.4)')
if met lt 0 then metdstr = '-' + metdstr
met_str =  metdstr + ':' $
    + string(meth, format='(I2.2)') + ':'  $
    + string(metm, format='(I2.2)') + ':'  $
    + string(mets, format='(F7.4)')
met_str = repchr(met_str,' ','0')
widget_control, (*info).met_val, set_value=met_str
return
end

pro map_date_quit, ev
tlb = ev.top
widget_control, tlb, get_uvalue=info
ptr_free, info
widget_control, tlb, /destroy
return
end

pro map_date_timer_event, ev
tlb = ev.top
widget_control, tlb, get_uvalue=info
sel = widget_info((*info).tz_droplist, /droplist_select)
now_jul = systime(/JULIAN) - 2450000.0d + (sel + 4)/24.0d0
widget_control, (*info).rjd_val,  $
    set_value=string(now_jul,format='(F15.9)')
phony_event = {id:(*info).rjd_val, top:(*info).tlb, handler:0L}
widget_control, (*info).rjd_val, send_event=phony_event
widget_control, (*info).mode_button, get_value=val
v = strtrim(val[0],2)
;  The button is showing the opposite of the current mode,
;  since you press it to get that.
if strpos(v, 'Calcul') ge 0 then begin
    widget_control, (*info).buttons, timer=1.0
endif
end

pro map_date_mode, ev
tlb = ev.top
widget_control, tlb, get_uvalue=info
widget_control, ev.id, get_value=val
v = strtrim(val[0],2)
if strpos(v,'Calcul') ge 0 then begin
    widget_control, ev.id, set_value='Press for Continuous Update'
    widget_control, (*info).buttons, /clear_events
endif
if strpos(v,'Continuous Update') ge 0 then begin
    widget_control, ev.id, set_value='      Press for Calculator Mode      '
    widget_control, (*info).buttons, timer=0.01
endif
return
end

pro map_date_event, ev
tlb = ev.top
widget_control, tlb, get_uvalue=info

if ev.id eq (*info).rjd_val then begin
    widget_control, (*info).rjd_val, get_value=rjd_text
    rjd = double(rjd_text[0])
    gmt = jul2gmt(rjd)
    ts = jul2ts(rjd)
    map_date_set_all, info, rjd, gmt, ts
endif

if ev.id eq (*info).tsday_val or ev.id eq (*info).tstms_val then begin
    widget_control, (*info).tsday_val, get_value=tsday_txt
    widget_control, (*info).tstms_val, get_value=tstms_txt
    ts = [ long(tsday_txt), long(tstms_txt) ]
    rjd = ts2jul(ts)
    rjd = rjd[0]
    ts2gmt, gmt, ts
    map_date_set_all, info, rjd, gmt, ts
endif

if ev.id eq (*info).gmtyear_val or ev.id eq (*info).gmtday_val  $
    or ev.id eq (*info).gmthr_val or ev.id eq (*info).gmtmn_val  $
    or ev.id eq (*info).gmtsc_val then begin
    
    widget_control, (*info).gmtyear_val, get_value=yr
    widget_control, (*info).gmtday_val, get_value=day
    widget_control, (*info).gmthr_val, get_value=hr
    widget_control, (*info).gmtmn_val, get_value=mn
    widget_control, (*info).gmtsc_val, get_value=sc

    yr = strmid('000000'+yr, 3, /rev)
    day = strmid('000000'+day, 2, /rev)
    hr = strmid('000000'+hr, 1, /rev)
    mn = strmid('000000'+mn, 1, /rev)
    sc2 = string(floor(double(sc)*1.0d4), format='(I6.6)')

    gmt = yr + day + hr + mn + sc2
    gmt2ts, ts, gmt
    rjd = gmt2jul(gmt)
    map_date_set_all, info, rjd, gmt, ts
endif
    
if ev.id eq (*info).gmtmonth_val or ev.id eq (*info).gmtdom_val then begin

    widget_control, (*info).gmtmonth_val, get_value=monstr
    widget_control, (*info).gmtdom_val, get_value=domstr
    widget_control, (*info).gmtyear_val, get_value=yr
    widget_control, (*info).gmthr_val, get_value=hr
    widget_control, (*info).gmtmn_val, get_value=mn
    widget_control, (*info).gmtsc_val, get_value=sc

    for i=0,11 do begin
        if strupcase(strmid((*info).months[i],0,3))  $
            eq strupcase(strmid(strtrim(monstr[0],2),0,3)) then mon_num = i+1
    endfor

    date = [ long(yr), mon_num, long(domstr) ]

    dayofyear, date, dd & day = string(floor(dd))

    yr = strmid('000000'+yr, 3, /rev)
    day = strmid('000000'+day, 2, /rev)
    hr = strmid('000000'+hr, 1, /rev)
    mn = strmid('000000'+mn, 1, /rev)
    sc2 = string(floor(double(sc)*1.0d4), format='(I6.6)')

    gmt = yr + day + hr + mn + sc2
    gmt2ts, ts, gmt
    rjd = gmt2jul(gmt)
    map_date_set_all, info, rjd, gmt, ts
endif

return
end

pro map_date

widget_control, default_font='7x13'
tlb = widget_base(col=1, title='MAP DATE AND TIME CONVERTER')

buttons = widget_base(tlb, row=1, event_pro='map_date_timer_event')
dates = widget_base(tlb, row=1)

gmt_col = widget_base(dates, col=1, frame=3)
map_col = widget_base(dates, row=3)

ts_col = widget_base(map_col, col=1, frame=3)
rjd_col = widget_base(map_col, col=1, frame=3)
met_col = widget_base(map_col, col=1, frame=3)

gmt_lab = widget_label(gmt_col, value='GMT')
ts_lab = widget_label(ts_col, value='MAP TIME STAMP')
rjd_lab = widget_label(rjd_col, value='MAP REDUCED JD')
met_lab = widget_label(met_col, value='TRUE MET')

gmtyear = widget_base(gmt_col, row=1)
gmtday = widget_base(gmt_col, row=1)
gmttime = widget_base(gmt_col, row=1)
gmtmonth = widget_base(gmt_col, row=1)
gmtdom = widget_base(gmt_col, row=1)

gmtyear_lab = widget_label(gmtyear, value='Year: ', xsize=50)
gmtyear_val = widget_text(gmtyear, value='2001', xsize=6, /editable)
gmtday_lab = widget_label(gmtday, value='Day:  ', xsize=50)
gmtday_val = widget_text(gmtday, value='200', xsize=6, /editable)
gmttime_lab = widget_label(gmttime, value='Time:  ', xsize=50)
gmthr_val = widget_text(gmttime, value='01', xsize=3, /editable)
gmtc1 = widget_label(gmttime, value=':', xsize=10)
gmtmn_val = widget_text(gmttime, value='01', xsize=3, /editable)
gmtc2 = widget_label(gmttime, value=':', xsize=10)
gmtsc_val = widget_text(gmttime, value='59.9999', xsize=8, /editable)
gmtmonth_lab = widget_label(gmtmonth, value='Month: ', xsize=50)
gmtmonth_val = widget_text(gmtmonth, value='July', xsize=10, /editable)
gmtdom_lab = widget_label(gmtdom, value='Day of Month:  ', xsize=105)
gmtdom_val = widget_text(gmtdom, value='30', xsize=4, /editable)

ts = widget_base(ts_col, row=1)

tsday_lab = widget_label(ts, value='Day:  ', xsize=40)
tsday_val = widget_text(ts, value='500', xsize=6, /editable)
tstms_lab = widget_label(ts, value='   1/10 ms:', xsize=80)
tstms_val = widget_text(ts, value='864000000', xsize=12, /editable)

rjd_val = widget_text(rjd_col, value='2100.8000001', xsize=14, /editable)
met_val = widget_label(met_col, value='1111:12:23:46.0001', xsize=150)

quit_button = widget_button(buttons, value='Quit',  $
    event_pro='map_date_quit')

help_button = widget_button(buttons, value='Help',  $
    event_pro='map_date_help')

mode_button = widget_button(buttons,  $
    value='      Press for Calculator Mode      ',  $
    event_pro='map_date_mode')

tz_droplist = widget_droplist(buttons,  $
    value=['4=EDT', '5=EST/CDT', '6=CST/MDT', '7=MST/PDT', '8=PST'])

widget_control, tz_droplist, set_droplist_select=1

load_map_params, params

months = [  $
    'January',  $
    'February',  $
    'March',  $
    'April',  $
    'May',  $
    'June',  $
    'July',  $
    'August',  $
    'September',  $
    'October',  $
    'November',  $
    'December' ]

info    = ptr_new({  $
        tlb: tlb,  $
        gmtyear_val: gmtyear_val,  $
        gmtday_val: gmtday_val,  $
        gmthr_val: gmthr_val,  $
        gmtmn_val: gmtmn_val,  $
        gmtsc_val: gmtsc_val,  $
        gmtmonth_val: gmtmonth_val,  $
        gmtdom_val: gmtdom_val,  $
        tsday_val: tsday_val,  $
        tstms_val: tstms_val,  $
        rjd_val: rjd_val,  $
        met_val: met_val,  $
        params: params,  $
        months: months,  $
        buttons: buttons,  $
        mode_button: mode_button,  $
        help_button: help_button,  $
        help: 0L,  $
        tz_droplist:tz_droplist,  $
        quit_button: quit_button })

widget_control, tlb, set_uvalue=info
widget_control, tlb, /realize

xmanager, 'map_date', tlb, event_handler='map_date_event', /no_block

widget_control, buttons, timer=0.01

return
end
