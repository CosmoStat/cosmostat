;
; Auto Save File For /Users/starck/Main/MRS/idl/./xmrs_tv.pro
;



; CODE MODIFICATIONS MADE ABOVE THIS COMMENT WILL BE LOST.
; DO NOT REMOVE THIS COMMENT: BEGIN HEADER




; DO NOT REMOVE THIS COMMENT: END HEADER
; CODE MODIFICATIONS MADE BELOW THIS COMMENT WILL BE LOST.


; CODE MODIFICATIONS MADE ABOVE THIS COMMENT WILL BE LOST.
; DO NOT REMOVE THIS COMMENT: BEGIN MAIN13




PRO MAIN13_Event, Event


  WIDGET_CONTROL,Event.Id,GET_UVALUE=Ev

  CASE Ev OF 

  'BUTTON3': BEGIN
      Print, 'Event for ReadFile'
      END
  'BUTTON4': BEGIN
      Print, 'Event for Reset'
      END
  'BUTTON5': BEGIN
      Print, 'Event for LUT'
      END
  'BUTTON11': BEGIN
      Print, 'Event for Quit'
      END
  'BUTTON14': BEGIN
      Print, 'Event for BUTTON14'
      END
  ENDCASE
END


; DO NOT REMOVE THIS COMMENT: END MAIN13
; CODE MODIFICATIONS MADE BELOW THIS COMMENT WILL BE LOST.



PRO xmrs_tv, GROUP=Group


  IF N_ELEMENTS(Group) EQ 0 THEN GROUP=0

  junk   = { CW_PDMENU_S, flags:0, name:'' }


  MAIN13 = WIDGET_BASE(GROUP_LEADER=Group, $
      ROW=3, $
      MAP=1, $
      UVALUE='MAIN13', $
      XSIZE=600, $
      YSIZE=600)

  BASE2 = WIDGET_BASE(MAIN13, $
      ROW=1, $
      SPACE=20, $
      FRAME=1, $
      FUNC_GET_VALUE='FirstLine_get', $
      KILL_NOTIFY='FirstLine_not', $
      MAP=1, $
      PRO_SET_VALUE='FirstLine_set', $
      UVALUE='BASE2', $
      XSIZE=585, $
      YSIZE=100)

  BUTTON3 = WIDGET_BUTTON( BASE2, $
      FRAME=10, $
      UVALUE='BUTTON3', $
      VALUE='Read', $
      XSIZE=100, $
      YSIZE=100)

  BUTTON4 = WIDGET_BUTTON( BASE2, $
      FRAME=10, $
      UVALUE='BUTTON4', $
      VALUE='Reset', $
      XSIZE=100, $
      YSIZE=100)

  BUTTON5 = WIDGET_BUTTON( BASE2, $
      FRAME=10, $
      UVALUE='BUTTON5', $
      VALUE='LUT', $
      XSIZE=100, $
      YSIZE=100)

  BUTTON11 = WIDGET_BUTTON( BASE2, $
      FRAME=10, $
      UVALUE='BUTTON11', $
      VALUE='Quit', $
      XSIZE=100, $
      YSIZE=100)


  BASE13 = WIDGET_BASE(MAIN13, $
      ROW=1, $
      FRAME=1, $
      MAP=1, $
      TITLE='tvline', $
      UVALUE='BASE13', $
      XSIZE=585, $
      YSIZE=100)

  BUTTON14 = WIDGET_BUTTON( BASE13, $
      UVALUE='BUTTON14')


  WIDGET_CONTROL, MAIN13, /REALIZE

  XMANAGER, 'MAIN13', MAIN13
END
