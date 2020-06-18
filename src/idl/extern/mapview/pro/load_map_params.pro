Pro Load_MAP_Params, params, File=file
;+
; NAME:
;       Load_MAP_Params
; PURPOSE:
;       Reads the contents of the general MAP program parameters file into
;       an IDL structure.
; CALLING SEQUENCE:
;       Load_MAP_Params, params
; OUTPUTS:
;       params - The structure.  See the body of this routine for details.
; KEYWORDS:
;       File   - The name of the parameters file.  Defaults to
;                $MAP_REF/programs.pars
; COMMENTS:
;       This is the IDL equivalent of the Fortran Load_Map_Params routine.
; PROCEDURES USED:
;       CONCAT_DIR(), GMT2JUL(), KEYFILE
; MODIFICATION HISTORY:
;       Written by Michael R. Greason, Raytheon ITSS, 02 July 1999.
;       Pass identifier added.  MRG, RITSS, 14 May 2001.
;       Updated for c1/flight_0.  RSH, RITSS, 26 June 2001.
;       Added question_event.  RSH, RITSS, 28 June 2001
;       Added true_launch_gmt. JP, 03 Jul 2001.
;       Added gal_mask_file. JP, 10 Aug 2001.
;       Raw and calibrated archive directories.  MRG, SSAI, 01 March 2002.
;       Synched with Fortran version.  RSH, SSAI, 10 June 2002.
;   	Use getenv to extract environment variables.  20 Aug 2002.
;       Use CONCAT_DIR() to concatenate directories  WL  May 2003
;-
on_error, 2
;
;                       Check arguments.
;
If (n_params() LE 0) Then message, 'Syntax:  Load_MAP_Params, params'
;
If (n_elements(file) LE 0) Then file = $
    concat_dir( strtrim(getenv('MAP_REF'),2),  'programs.pars')
;
;                       Define the output structure.
;
params = {Pass:               'N/A',                   $ ; Pass identifier.
          Mod_GMT:            '',                      $ ; Modification date
          CCSDS_GMT:          '',                      $ ; CCSDS ref. date
          True_Launch_GMT:    '',                      $ ; Actual Launch date
          Launch_GMT:         '',                      $ ; Launch date used as TS reference
          CCSDS2Launch:       0.0d0,                   $ ; CCSDS->MET Conversion.
          Sim_UTCF:           0.0d0,                   $ ; UTCF value for SC sims
          Sim_LeapSec:        0L,                      $ ; LeapSec val for SC sims
          UTCF_File:          '',                      $ ; dir+name of utcf table
          Archive_Raw:        '',                      $ ; Raw archive directory
          Archive_Cal:        '',                      $ ; Calibrated archive directory
          Archive_Index_File: '',                      $ ; Index file
          Telemetry_Index_File: '',                    $ ; Index file
          Questionable_Event_File: '',                 $ ; Questionable event file
          Gal_Mask_File:'',                            $ ; Galaxy Mask file
          DA_Name:            replicate('', 10),       $ ; DA mnemonic
          Chan_Name:          replicate('', 4, 10),    $ ; Channel mnemonic
          Frame_Duration:     0.0,                     $ ; Seconds
          TOD_Per_Frame:      replicate(0L, 10),       $ ; Science points per frame
          Spin_Period:        0.0d0,                   $ ; Seconds
          Precession_Period:  0.0d0,                   $ ; Seconds
          Precession_Angle:   0.0d0,                   $ ; Degrees
          Dir_A_LOS:          replicate(0.0d0, 3, 10), $ ; A side, each d/a
          Dir_B_LOS:          replicate(0.0d0, 3, 10), $ ; B side, each d/a
          Dir_A_Pol:          replicate(0.0d0, 3, 10), $ ; A side, each d/a
          Dir_B_Pol:          replicate(0.0d0, 3, 10), $ ; B side, each d/a
          Gain:               replicate(0.0, 4, 10),   $ ; mK per count
          Baseline:           replicate(0.0, 4, 10),   $ ; counts
          Noise_RMS:          replicate(0.0, 4, 10),   $ ; mK per observation
          F_Knee:             replicate(0.0, 4, 10),   $ ; power spec 1/f knee
          R34:                replicate(0.0, 2, 10),   $ ; linear corr, ch 3 & 4
          T_Monopole:         0.0,                     $ ; K
          T_Dipole:           replicate(0.0, 3),       $ ; mK, galactic J2000 [X,Y,Z] coordinates
          Ant_to_Th:          replicate(0.0, 10),      $ ; convert ant to th temp.
          Mars_Cut:           replicate(0.0, 10),      $ ; Degrees from beam center
          Jupiter_Cut:        replicate(0.0, 10),      $ ; Degrees from beam center
          Saturn_Cut:         replicate(0.0, 10),      $ ; Degrees from beam center
          Galaxy_Cut:         replicate(0.0, 10),      $ ; Degrees from beam center
          Uranus_Cut:         replicate(0.0, 10),      $ ; Degrees from beam center
          Neptune_Cut:        replicate(0.0, 10),      $ ; Degrees from beam center
          Sun_Shield_Cut:     0.0,                     $ ; Degrees from beam center
          Earth_Shield_Cut:   0.0,                     $ ; Degrees from beam center
          Moon_Shield_Cut:    0.0,                     $ ; Degrees from beam center
          DistFromEarth:      0.0d0                    $ ; MAP-Earth distance, km.
          }
;
;                       Read the file.
;
KeyFile, file, keyvals, Comment='#', Continue='&', Separator='='
;
;                       Fill the structure from the file data.
;
sc = -32768 & sc_s = -32768 & sc_e = -32768 & sc_m = -32768
n = n_elements(keyvals) - 1
For i = 0, n Do Begin
        Case (strupcase(keyvals[i].Key)) Of
          'PASS'              : params.Pass               = keyvals[i].Value
          'MOD_GMT'           : params.Mod_GMT            = keyvals[i].Value
          'CCSDS_GMT'         : params.CCSDS_GMT          = keyvals[i].Value
          'TRUE_LAUNCH_GMT'   : params.True_Launch_GMT    = keyvals[i].Value
          'LAUNCH_GMT'        : params.Launch_GMT         = keyvals[i].Value
          'SIM_UTCF'          : params.Sim_UTCF           = double(keyvals[i].Value)
          'SIM_LEAPSEC'       : params.Sim_LeapSec        = keyvals[i].Value
          'UTCF_FILE'         : params.UTCF_File          = keyvals[i].Value
          'ARCHIVE_RAW'       : params.Archive_Raw        = keyvals[i].Value
          'ARCHIVE_CAL'       : params.Archive_Cal        = keyvals[i].Value
          'ARCHIVE_INDFIL'    : params.Archive_Index_File = keyvals[i].Value
          'TELEM_INDFIL'      : params.Telemetry_Index_File = keyvals[i].Value
          'QUESTION_EVENTS'   : params.Questionable_Event_File = keyvals[i].Value
          'GAL_MASK_FILE'     : params.Gal_Mask_File      = keyvals[i].Value
          'SPIN_PERIOD'       : params.Spin_Period        = double(keyvals[i].Value)
          'PRECESSION_PERIOD' : params.Precession_Period  = double(keyvals[i].Value)
          'PRECESSION_ANGLE'  : params.Precession_Angle   = double(keyvals[i].Value)
          'DIST_FROM_EARTH'   : params.DistFromEarth      = double(keyvals[i].Value)
          'FRAME_DURATION'    : params.Frame_Duration     = float(keyvals[i].Value)
          'TOD_PER_FRAME'     : params.TOD_Per_Frame      = long(KeyArray(keyvals[i].Value))
          'T_MONOPOLE'        : params.T_Monopole         = float(keyvals[i].Value)
          'T_DIPOLE'          : params.T_Dipole           = float(KeyArray(keyvals[i].Value))
          'ANT_TO_TH'         : params.Ant_to_Th          = float(KeyArray(keyvals[i].Value))
          'GAIN'              : params.Gain               = float(KeyArray(keyvals[i].Value))
          'BASELINE'          : params.Baseline           = float(KeyArray(keyvals[i].Value))
          'NOISE_RMS'         : params.Noise_RMS          = float(KeyArray(keyvals[i].Value))
          'F_KNEE'            : params.F_Knee             = float(KeyArray(keyvals[i].Value))
          'R34'               : params.R34                = float(KeyArray(keyvals[i].Value))
          'DIR_A_LOS'         : params.Dir_A_LOS          = double(KeyArray(keyvals[i].Value))
          'DIR_B_LOS'         : params.Dir_B_LOS          = double(KeyArray(keyvals[i].Value))
          'DIR_A_POL'         : params.Dir_A_Pol          = double(KeyArray(keyvals[i].Value))
          'DIR_B_POL'         : params.Dir_B_Pol          = double(KeyArray(keyvals[i].Value))
          'MARS_CUT'          : params.Mars_Cut           = replicate(keyvals[i].Value, 10)
          'JUPITER_CUT'       : params.Jupiter_Cut        = replicate(keyvals[i].Value, 10)
          'SATURN_CUT'        : params.Saturn_Cut         = replicate(keyvals[i].Value, 10)
          'URANUS_CUT'        : params.Uranus_Cut         = replicate(keyvals[i].Value, 10)
          'NEPTUNE_CUT'       : params.Neptune_Cut        = replicate(keyvals[i].Value, 10)
          'GALAXY_CUT'        : params.Galaxy_Cut         = replicate(keyvals[i].Value, 10)
          'SHIELD_CUT'        : sc                        = keyvals[i].Value
          'SUN_SHIELD_CUT'    : sc_s                      = keyvals[i].Value
          'EARTH_SHIELD_CUT'  : sc_e                      = keyvals[i].Value
          'MOON_SHIELD_CUT'   : sc_m                      = keyvals[i].Value

          'ARCHIVE_DIR'       : a = keyvals[i].Value

          Else : print,' not recognized: ',strupcase(keyvals[i].Key)
        EndCase
EndFor
;
;                       Now fill the "computed" values.
;
params.DA_Name   = ['K', 'Ka', 'Q1', 'Q2', 'V1', 'V2', 'W1', 'W2', 'W3', 'W4']
params.Chan_Name = transpose([[params.DA_Name + 'A'], $
                              [params.DA_Name + 'B'], $
                              [params.DA_Name + 'C'], $
                              [params.DA_Name + 'D']])

If sc LT 0 Then Begin
    params.sun_shield_cut = 0.0
    params.moon_shield_cut = 0.0
    params.earth_shield_cut = 0.0
EndIf Else Begin
    params.sun_shield_cut = sc
    params.moon_shield_cut = sc
    params.earth_shield_cut = sc
EndElse

If sc_s GE 0 Then params.sun_shield_cut = sc_s
If sc_m GE 0 Then params.moon_shield_cut = sc_m
If sc_e GE 0 Then params.earth_shield_cut = sc_e
    
;
jc = GMT2Jul(params.CCSDS_GMT)
jm = GMT2Jul(params.Launch_GMT)
params.CCSDS2Launch = (jm - jc) * (24.0d0 * 3600.0d4)
;
Return
End
