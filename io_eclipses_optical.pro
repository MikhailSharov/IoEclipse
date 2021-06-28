
; Analysis of Io's optical spectra in eclipse from datasets with Apache Point and Large Binocular telescopes
; Written by C. Schmidt, M. Sharov. BU Center for Space Physics, 2019-2021


FUNCTION Match_Reference_Spectrum, X, P
  ; Multiply P[0], Add P[1] and Smooth P[2] a reference spectrum (X) until it best matches Y
  return, P[0]*gauss_smooth(X, P[2], /EDGE_TRUNCATE) + P[1]
end

FUNCTION Match_Telluric_Absorption, X, P
  ; Multiply P[0], Add P[1] and exponent P[2] a reference spectrum (X) until it best matches Y
  return, P[0]*(X^P[2]) + P[1]
end

FUNCTION Jup_Smooth_Shift, X, S
  ; Smooth S[0], then Shift S[1] a jupiter scatter spectrum (X) until it best matches Y
  return, interpolate(gauss_smooth(X, S[0], /EDGE_TRUNCATE), findgen(n_elements(X)) - S[1])
end

FUNCTION Match_Jup_Combo, X, P
  ; Smooth P[2], the Shift P[3] some transformed (X) with multiplicative P[0] and offset P[1]
  ; This is a combination of the previous two functions such that they occur at the same time
  return, interpolate(gauss_smooth((P[0]*X)+ P[1], P[2], /EDGE_TRUNCATE), findgen(n_elements(X)) - P[3])
end

function Gaussian_For_MPFIT, p, x=x, y=y, err=err, fit=fit
  z = (x - p[1])/p[2]
  fit = p[0]*exp(-z^2/2d)
  return, (y - fit)/err
end

FUNCTION MX_plus_B, X, P
  ; multipliy P[0] and offset P[1]
  return, P[0]*X+ P[1]
end

function shift_smooth, shift_and_smooth_by
  Common shift_smooth_common, spec, jup, correl_indicies, shift_only, debug

  ; shift_and_smooth_by[0] = # of pixels to smooth
  ; shift_and_smooth_by[1] = # of pixels to shift

  ;widen            = minmax(correl_indicies) + [-10, 10] ; broaden the potion of the spectrum we're manipulating to avoid edge effects, while preserving speed

  A                = spec[correl_indicies]
  match_to_scatter = interpolate( GAUSS_SMOOTH(jup, shift_and_smooth_by[0]), findgen(n_elements(jup)) - shift_and_smooth_by[1], Missing = !values.F_NaN)
  ;match_to_scatter = interpolate( GAUSS_SMOOTH(jup[widen[0]:widen[1]], shift_and_smooth_by[0]), findgen(n_elements(jup[widen[0]:widen[1]])) - shift_and_smooth_by[1], Missing = !values.F_NaN)
  B                = match_to_scatter[correl_indicies]
  correl           = correlate( A, B, /double )
  ;print, 'Smooth:', shift_and_smooth_by[0], ', Shift:', shift_and_smooth_by[1], ', Correlation:', correl
  ;if keyword_set(debug) then print, 'Smooth:', shift_and_smooth_by[0], ', Shift:', shift_and_smooth_by[1], ', Correlation:', correl
  ;if shift_and_smooth_by[0] lt 0. then correl = 0.1 ; GAUSS_SMOOTH by < 0. always returns the input anyhow.
  shift_only = shift_and_smooth_by[1]
  return, 1.D / correl
end

FUNCTION Shift_spec, X, P
  return, interpolate(X, findgen(n_elements(X)) - P[0])
end

FUNCTION med_filter, X, P ; Sigma filter in 1-dimension, if an array differs from it's median by P[1] standard deviations replace it with the median
  ; p[0] = width to evaluate the median
  ; p[1] = standard deviation about which to reject
  x_out = x
  med = median(x, P[0], /even)
  replace_pixels = where(abs(med - x) gt p[1]*stddev(x), /null)
  x_out[replace_pixels] = med[replace_pixels]
  return, x_out
end

; Notes from Relevant Papers:
;  cross-sections for ** 25 eV ** electrons on SO2:
;       1302 ---- 1.20e-19 at 30eV Vatti-Palle et al. 2004
;       1304 ---- 7.40e-20
;       1306 ---- 3.0 e-20
;1356 + 1358 ---- 1.88e-19 at 25eV Vatti-Palle et al. 2004 see their table 6.
;1356 + 1358 ---- 2.59e-19 at 30eV Vatti-Palle et al. 2004
;       5577 ---- 8.25e-19 Kedzierski et al. (2000) from Table 1 of Mcconkey & Kedzierski, Wladyslaw. (2014). Advances In Atomic, Molecular, and Optical Physics, 63, 1-46.
;       7774 ---- 3.90e-19 Ajello et al. (2008)
;       8446 ---- 2.15e-19 Ajello et al. (2008)
;       9225 ---- 6.63e-19 ( Individual lines of the triplet are 3.63e-19 + 1.76e-19 + 1.24e-19; Ajello et al. (2008) )
; And for ** 100 eV ** electrons on SO2
;       1302 ---- 1.28e-18 Vatti-Palle et al. 2004
;       1304 ---- 7.30e-19
;       1306 ---- 2.90e-19
;       5577 ---- 2.10e-18 Kedzierski et al. Can. J. Phys. 78: 617–624 (2000. Mcconkey & Kedzierski, Wladyslaw. (2014). Advances In Atomic, Molecular, and Optical Physics, 63, 1-46.
;       7774 ---- 4.20e-18
;       8446 ---- 2.60e-18
;       9225 ---- 3.40e-18 ( Individual lines of the triplet are 1.80e-18 + 9.21e-19 + 6.78e-19 )

; Geissler et al. (1999) Galileo *disk-integrated* brightness: (no narrow band disk-integrated brightess in 2001 or 2004 paper)
;       Green (510-605nm) --- 8.0 kR
;       Red   (615-710nm) --- 6.8 kR

Pro Io_Eclipses_optical, Part=Part, Date=Date, Sunlit=Sunlit
  Common shift_smooth_common, spec, jup, correl_indicies, shift_only, debug

  if date ne !null then begin
    Case 1 of
      date eq 'UT180320': begin
        ingress               = 1
        Penumbra_UTC          = '2018-Mar-20 10:46:46'
        Umbra_UTC             = '2018-Mar-20 10:50:22'
        dir                   = 'D:\DATA\Apache Point\Echelle\Io Eclipses\'+date+'\'
        reduced_dir           = 'D:\DATA\Apache Point\Echelle\Io Eclipses\'+date+'\Reduced\'
        calibration_dir       = 'D:\DATA\Apache Point\Echelle\Io Eclipses\UT190812\Reduced\'  ; Directory with Stellar Spectra for Telluric Corrections
        Jovian_Scatter_files  = [ 'Jovian_Scatter.0001', 'Jovian_Scatter.0002', 'Jovian_Scatter.0003' ]
        Eclipse_files         = [ 'Io_eclipsed.000'+strcompress(indgen(6)+1, /remove_all) ]
        Jupiter_Center_File   = 'Jupiter_Center_Spectrum.0001'
        Standard_Star_files   = ['HD_159975_BFR.0001','HD_155379_A0V.0001']
      end
      date eq 'UT190812': begin      ; Good pointing very high arimass. Nice for 6300A
        ingress               = 0
        Penumbra_UTC          = '2019-Aug-12 03:44:40'
        Umbra_UTC             = '2019-Aug-12 03:48:19'
        dir                   = 'D:\DATA\Apache Point\Echelle\Io Eclipses\'+date+'\'
        reduced_dir           = 'D:\DATA\Apache Point\Echelle\Io Eclipses\'+date+'\Reduced\'
        calibration_dir       = 'D:\DATA\Apache Point\Echelle\Io Eclipses\'+date+'\Reduced\'  ; Directory with Stellar Spectra for Telluric Corrections
        if sunlit eq 1 then Eclipse_files         = ['Io_Penumbra.0001','Io_Sunlit.0001','Io_Sunlit.0002','Io_Sunlit.0003', 'Io_Sunlit.0004'] ;use this for Na D in sunlight
        if sunlit eq 0 then Eclipse_files         = ['Io_eclipsed.000'+strcompress(indgen(7)+3, /remove_all), 'Io_eclipsed.00'+strcompress(indgen(1)+10, /remove_all) ]
        Jupiter_Center_File   = 'Jupiter_Disk_Center.0002'
        Jovian_Scatter_files  = ['Jupiter_Approaching_Limb.0001']
        Standard_Star_files   = ['HD_159975_BFR.0001','HD_155379_A0V.0001']
      end
      date eq 'UT200823': begin      ; **IFFY pointing** --> After checking pointing appears alright
        ingress               = 0
        Penumbra_UTC          = '2020-Aug-23 03:10:29'
        Umbra_UTC             = '2020-Aug-23 03:14:02'
        dir                   = 'D:\DATA\Apache Point\Echelle\Io Eclipses\'+date+'\'
        reduced_dir           = 'D:\DATA\Apache Point\Echelle\Io Eclipses\'+date+'\Reduced\'
        calibration_dir       = 'D:\DATA\Apache Point\Echelle\Io Eclipses\'+date+'\Reduced\'  ; Directory with Stellar Spectra for Telluric Corrections
        if sunlit eq 0 then Eclipse_files         = ['Io_eclipsed.000'+strcompress(indgen(5)+3, /remove_all)]
        if sunlit eq 1 then Eclipse_files         = ['Io_Sunlit.0001','Io_Sunlit.0002','Io_Sunlit.0003', 'Io_Sunlit.0004'] ;use this for sunlight
        Jupiter_Center_File   = 'Jupiter_Disk_Center.0003'
        Jovian_Scatter_files  = ['Jupiter_Scatter.0001', 'Jupiter_Scatter.0002', 'Jupiter_Scatter.0003', 'Jupiter_Scatter.0004']
        Standard_Star_files   = ['A0V_HD177213.0003', 'A0V_HD177213.0004']
      end
      date eq 'UT200908': begin      ; **IFFY pointing** 19, 20, 23 are good, rest iffy at best
        ingress               = 0
        Penumbra_UTC          = '2020-Sep-08 01:29:13'
        Umbra_UTC             = '2020-Sep-08 01:32:46'
        dir                   = 'D:\DATA\Apache Point\Echelle\Io Eclipses\'+date+'\'
        reduced_dir           = 'D:\DATA\Apache Point\Echelle\Io Eclipses\'+date+'\Reduced\'
        calibration_dir       = 'D:\DATA\Apache Point\Echelle\Io Eclipses\'+date+'\Reduced\'  ; Directory with Stellar Spectra for Telluric Corrections
        if sunlit eq 0 then Eclipse_files         = ['Io_eclipsed.0019', 'Io_eclipsed.0020', 'Io_eclipsed.0023']  ;  until frame 22 Frame 24 is clearly in penumbra despite the name
        if sunlit eq 1 then Eclipse_files         = ['Io_Free_and_Clear.0028', 'Io_Free_and_Clear.0029', 'Io_Free_and_Clear.0031', 'Io_Free_and_Clear.0033'] ;use this for sunlight
        Jupiter_Center_File   = 'Jupiter_Disk_Center.0016'
        Jovian_Scatter_files  = ['Jupiter_Scatter.0017', 'Jupiter_Scatter.0030', 'Jupiter_Scatter.0032']
        Standard_Star_files   = ['HD_163336.0013', 'HD_163336.0014']
      end
      date eq 'UT201001': begin       ; Good data
        ingress               = 0
        Penumbra_UTC          = '2020-Oct-01 01:43:52'
        Umbra_UTC             = '2020-Oct-01 01:47:23'
        dir                   = 'D:\DATA\Apache Point\Echelle\Io Eclipses\'+date+'\'
        reduced_dir           = 'D:\DATA\Apache Point\Echelle\Io Eclipses\'+date+'\Reduced\'
        calibration_dir       = 'D:\DATA\Apache Point\Echelle\Io Eclipses\'+date+'\Reduced\'  ; Directory with Stellar Spectra for Telluric Corrections
        if sunlit eq 1 then Eclipse_files         = ['Io_Penumbra.0001','Io_Free_and_Clear.0001','Io_Free_and_Clear.0002','Io_Free_and_Clear.0003', 'Io_Free_and_Clear.0004', 'Io_Free_and_Clear.0005', 'Io_Free_and_Clear.0006', 'Io_Free_and_Clear.0008']     ; Use this for sunlight
        if sunlit eq 0 then Eclipse_files         = ['Io_eclipsed.000'+strcompress(indgen(7)+1, /remove_all)]
        Jupiter_Center_File   = 'Jupiter_Disk_Center.0001'
        Jovian_Scatter_files  = ['Jupiter_Scatter.0001', 'Jupiter_Scatter.0002', 'Jupiter_Scatter.0003']
        Jovian_Scatter_files  = ['Jupiter_Scatter.0001', 'Jupiter_Scatter.0002', 'Jupiter_Scatter.0003']
        Standard_Star_files   = ['Chi_Cap_A0V.0001', 'Chi_Cap_A0V.0002']
      end
      date eq 'UT201017': begin      ; **Bad pointing** 22 & 23 are good
        ingress               = 0
        Penumbra_UTC          = '2020-Oct-17 00:03:33'
        Umbra_UTC             = '2020-Oct-17 00:07:04'
        dir                   = 'D:\DATA\Apache Point\Echelle\Io Eclipses\'+date+'\'
        reduced_dir           = 'D:\DATA\Apache Point\Echelle\Io Eclipses\'+date+'\Reduced\'
        calibration_dir       = 'D:\DATA\Apache Point\Echelle\Io Eclipses\UT201001\Reduced\'  ; Directory with Stellar Spectra for Telluric Corrections
        if sunlit eq 1 then Eclipse_files         = ['Io_Penumbra.0025','Io_Free_and_Clear.0026','Io_Free_and_Clear.0027','Io_Free_and_Clear.0028', 'Io_Free_and_Clear.0029', 'Callisto.0015']
        if sunlit eq 0 then Eclipse_files         = ['Io_eclipsed.0022', 'Io_eclipsed.0023']  ; until frame 22 Frame 24 is clearly in penumbra despite the name
        Jupiter_Center_File   = 'Jupiter_Disk_Center.0014'
        Jovian_Scatter_files  = ['Jupiter_Disk_Center.0014', 'Io_eclipsed.0020'] ; io eclipsed 20 has no sign of io emissions and may be a good background
        Standard_Star_files   = ['Chi_Cap_A0V.0001', 'Chi_Cap_A0V.0002']
      end
      date eq 'UT210601': begin      ; 
        ingress               = 1    ;
        Penumbra_UTC          = '2020-Oct-17 00:03:33' ;fix
        Umbra_UTC             = '2020-Oct-17 00:07:04' ;fix
        dir                   = 'D:\DATA\Apache Point\Echelle\Io Eclipses\'+date+'\'
        reduced_dir           = 'D:\DATA\Apache Point\Echelle\Io Eclipses\'+date+'\Reduced\'
        calibration_dir       = 'D:\DATA\Apache Point\Echelle\Io Eclipses\UT210601\Reduced\'  ; Directory with Stellar Spectra for Telluric Corrections
        Eclipse_files         = ['Ganymede_Full_Shadow.0008', 'Ganymede_Full_Shadow.0009', 'Ganymede_Full_Shadow.00'+strcompress(indgen(9)+10, /remove_all)] 
        Jupiter_Center_File   = 'Jupiter_Disk_Center.0024'
        Jovian_Scatter_files  = ['Jupiter_Scatter.0019', 'Jupiter_Scatter.0020', 'Jupiter_Scatter.0025', 'Jupiter_Scatter.0026'] ;
        Standard_Star_files   = ['SAO146044_A0V.0021', 'SAO146044_A0V.0022']
      end
      date eq 'UT210609': begin      ; Somehow pointing drifted after the first frame. Working theory is that a shift in wind direction moved the telescope. 
        ingress               = 1    ; bumping the pointing with a rough guess based on Europa clearly increase the 6300 signal. This data should not be used for a time series. Frustrating. 
        Penumbra_UTC          = '2020-Oct-17 00:03:33' ;fix
        Umbra_UTC             = '2020-Oct-17 00:07:04' ;fix
        dir                   = 'D:\DATA\Apache Point\Echelle\Io Eclipses\'+date+'\'
        reduced_dir           = 'D:\DATA\Apache Point\Echelle\Io Eclipses\'+date+'\Reduced\'
        calibration_dir       = 'D:\DATA\Apache Point\Echelle\Io Eclipses\UT210609\Reduced\'  ; Directory with Stellar Spectra for Telluric Corrections
        Eclipse_files         = ['Io_FullShadow.0014', 'Io_FullShadow.0015'];, 'Io_FullShadow.0016', 'Io_FullShadow.0017', 'Io_FullShadow_Bumped_pointing.0018', $
        ;                          'Io_FullShadow_Bumped_pointing.0019', 'Io_FullShadow_Bumped_pointing.0020', 'Io_FullShadow_Bumped_pointing.0021'] 
        Jupiter_Center_File   = 'Jupiter_Disk_Center.0027'
        Jovian_Scatter_files  = ['Jupiter_Scatter.0023', 'Jupiter_Scatter.0024'] ;
        Standard_Star_files   = ['gam_Aqr_A0V.0010', 'SAO146044_A0V.0022']
      end

      date eq 'UT180807': begin      ; Keck Data from Katherine de Kleer
        Directory             = 'D:\DATA\Keck\Io Eclipse HIRES\Katherine\'
        reduced_dir           = 'D:\DATA\Apache Point\Echelle\Io Eclipses\Reduced\'
        Penumbra_UTC          = '2018-Aug-07 06:26:18.600'
        Umbra_UTC             = '2018-Aug-07 06:29:48.000'                ; This is a rough number
        Io_frames             = [287,289,291,293,295,297,299,301,303,304] ; These are the ones not totally swamped by Jupiter. The first frame has massive scattered light, probably only useful for 6300A, maybe Na
        Jupiter_frames        = [283, 305]                                ; Jupiter disk center before and after eclipse. Each gives a few percent difference, not worth worrying about.
        Eclipse_files         = Io_frames
      end
      date eq 'UT190424': begin
        ingress               = 1
        dir                   = 'D:\DATA\LBT\'
        reduced_dir           = 'D:\DATA\LBT\Reduced\'
        calibration_dir       = 'D:\DATA\LBT\Reduced\'
        Penumbra_UTC          = '2019-Apr-24 10:15:12'
        Umbra_UTC             = '2019-Apr-24 10:18:52'
        Jupiter_Center_File   = ['pepsib.20190424.069.dxt.ffc.rec', 'pepsib.20190424.070.dxt.ffc.rec',  $   ; taken over 4 cross-disperser settings (2,3,4,6) in two exposures 3 minutes apart
          'pepsir.20190424.049.dxt.ffc.rec', 'pepsir.20190424.050.dxt.ffc.rec']      ; [4265-4800A, 4800-5441A, 5441-6278A, 7419-9067A ]

        Sunlit_files_D        = [['pepsib.20190424.055.dxt.ffc.rec', 'pepsir.20190424.035.dxt.ffc.rec'], $  ; [2,4]  09:41:44.7 Sunlit
          ['pepsib.20190424.056.dxt.ffc.rec', 'pepsir.20190424.036.dxt.ffc.rec'], $  ; [3,6]  09:46:57.1 Sunlit
          ['pepsib.20190424.059.dxt.ffc.rec', 'pepsir.20190424.039.dxt.ffc.rec'], $  ; [2,4]  10:10:09.0 Sunlit
          ['pepsib.20190424.060.dxt.ffc.rec', 'pepsir.20190424.040.dxt.ffc.rec']]    ; [3,6]  10:15:24.5 Penumbra

        Eclipse_files_D      = [['pepsib.20190424.061.dxt.ffc.rec', 'pepsir.20190424.041.dxt.ffc.rec'], $   ; [2,4]  10:20:40.1
          ['pepsib.20190424.062.dxt.ffc.rec', 'pepsir.20190424.042.dxt.ffc.rec'], $   ; [3,6]  10:25:55.6
          ['pepsib.20190424.063.dxt.ffc.rec', 'pepsir.20190424.043.dxt.ffc.rec'], $   ; [2,4]  10:32:53.1
          ['pepsib.20190424.064.dxt.ffc.rec', 'pepsir.20190424.044.dxt.ffc.rec'], $   ; [3,6]  10:40:08.3
          ['pepsib.20190424.065.dxt.ffc.rec', 'pepsir.20190424.045.dxt.ffc.rec'], $   ; [2,4]  10:47:36.4
          ['pepsib.20190424.066.dxt.ffc.rec', 'pepsir.20190424.046.dxt.ffc.rec'], $   ; [3,6]  10:54:52.7
          ['pepsib.20190424.067.dxt.ffc.rec', 'pepsir.20190424.047.dxt.ffc.rec'], $   ; [2,4]  11:02:49.6
          ['pepsib.20190424.068.dxt.ffc.rec', 'pepsir.20190424.048.dxt.ffc.rec']]     ; [3,6]  11:10:04.7
        Eclipse_files        = Eclipse_files_D
      end
    endcase
  endif else print, 'Plotting combined results from all eclipses...'

  ;-------------------------------------------Load SPICE Data-----------------------------------------------------
  ; Clean any lingering kernels out of memory here:
    cspice_ktotal, 'all', count
    Print, 'Deleting ', strtrim(string(count),2), ' old SPICE kernels from memory'
    i=0
    while i lt count do begin
      cspice_kdata, 0, 'all', file, type, source, handle, found
      cspice_unload, file
      i = i + 1
    endwhile

  ; Load New Kernels
    CSPICE_FURNSH, STRCOMPRESS('C:\SPICE\generic_kernels\lsk\naif0010.tls')         ; leap seconds kernel
    CSPICE_FURNSH, STRCOMPRESS('C:\SPICE\generic_kernels\pck\pck00010.tpc')         ; Planet rotational states
    CSPICE_FURNSH, STRCOMPRESS('C:\SPICE\Jupiter_System\jup310.bsp')
    CSPICE_FURNSH, STRCOMPRESS('C:\SPICE\generic_kernels\spk\planets\de421.bsp')    ; SPK (ephemeris kernel) for planets
    CSPICE_FURNSH, STRCOMPRESS('C:\SPICE\generic_kernels\spk\satellites\sat319.bsp'); SPK (ephemeris kernel) for satellites
    cspice_ktotal, 'all', count
    Print, 'Loaded ', strtrim(string(count),2), ' new Spice kernels'

  ; If the last call of this code crashed in the middle of writing a postscipt, close the postscript and set things to a plot window
    if !D.NAME eq 'PS' then begin
      device, /close
      set_plot, 'win'
    endif

  ;--------------------------------------

  ;-------------------------------------------Load Constants-----------------------------------------------------
  ; Define Rest wavelengths
    Na1 = 5889.95095
    Na2 = 5895.92424
    Na3 = 8183.256
    Na4 = 8194.824
    K1 = 7664.89913
    K2 = 7698.96456
    O1 = 5577.330
    O2 = 6300.304
    O3 = 6363.776  ; Should be fixed equal to 6300 / 2.997 if optically thin. Storey and Zeippen (2000) theory = Sharpee and Slanger (2006) measured nightglow
    O4 = 7771.944  ; Directly produced by e- smashing of SO2, See Ajello et al. (2008)
    O5 = 7774.166
    O6 = 7775.388
    O7 = 8446.25   ; Directly produced by e- smashing of SO2, See Ajello et al. (2008)
    O8 = 8446.36
    O9 = 8446.76
    S1 = 9212.865  ; Directly produced by e- smashing of SO2, See Ajello et al. (2008)
    S2 = 9228.092
    S3 = 9237.538
    Cl1 = 8084.51
    Cl2 = 8085.56
    Cl3 = 8086.67
    Cl4 = 8087.73
    Cl5 = 8428.254
    Cl6 = 8575.24
    Cl7 = 8585.97
    Cl8 = 8375.94
    Cl9 = 8333.307
    Cl10 = 8221.742
    Cl11 = 8212.038

  ; Define Rest wavelengths
    Na = [5889.95095, 5895.92424, 8183.256, 8194.824]
    K  = [7664.89913, 7698.96456]
    O  = [5577.330, 6300.304, 6363.776, 7771.944, 7774.166, 7775.388, 8446.25, 8446.36, 8446.76]
    S  = [9212.865, 9228.092, 9237.538]                         ; See Ajello et al. 2008
    SO = [9549.18, 9626.21]                                     ; 0-0 and 1-1 band heads. See Setzer et al. Journal of Molecular Spectroscopy 198, 163–174 (1999), converted to Air wavelength
    C  = [8335.15, 9405.73]                                     ; Worth a check but nothing here
    Cl = [8085.56, 8086.67, 8375.94, 8585.97, 9121.15, 9592.22] ; Worth a check but nothing here

  ; MUST Run part zero if you redefine this
    Io_Airglow         = [Na, K, O, S, SO]                      ; Get line list for individual emission components
    Undefine, Io_Airglow_params
    if date ne !null then begin
      Io_Airglow_params =  {line:         Io_Airglow, $
        brightness:     fltarr(N_elements(Io_Airglow), N_elements(Eclipse_files)), $
        err_brightness: fltarr(N_elements(Io_Airglow), N_elements(Eclipse_files)), $
        linewidth:      fltarr(N_elements(Io_Airglow), N_elements(Eclipse_files)), $
        linecenter:     fltarr(N_elements(Io_Airglow), N_elements(Eclipse_files)), $
        exptime:        fltarr(N_elements(Eclipse_files)), $
        torus_latitude: fltarr(N_elements(Eclipse_files)), $
        T_P_Shadow:     fltarr(N_elements(Eclipse_files)), $
        Umbra_ET:       0.d, $
        Penumbra_ET:    0.d }
    endif

  ; Get the telluric airglow line list
    Airglow_threshold  = 2.e2 ; we don't care much about faint telluric airglow, so threshold it here
    Files = file_search('D:\DATA\Solar and Telluric Spectra\Telluric_airglow\UVES_Sky\*.tfits', count = n_files) ; Load a telluric emission line list
    AG_WL = [] & AG_flux = [] & AG_FWHM = []
    for i = 0, N_files-1 do begin
      UVES    = MRDFITS(files[i], 1, header, /USE_COLNUM, /silent )
      AG_WL   = [AG_WL,UVES.c1]
      AG_flux = [AG_flux,UVES.c2]
      AG_FWHM = [AG_FWHM,UVES.c3]
    endfor
    keep   = where(AG_flux gt Airglow_threshold, /Null)
    Telluric_Airglow = AG_WL[Keep] & Telluric_Airglow_Flux = AG_Flux[Keep]

  ; Read in torus latitude data
    READCOL,'D:\DATA\___Calibration___\Carl_centrifugal_equator.txt', torus_lat_out, torus_deg, skipline = 1, /Silent
 
  ; plot the results of all nights if "part" is not specified   
    if part eq !null then goto, plot_combined

  ;================================Part 0 Generate Telluric Absorption Spectrum==================================================================================
  if part eq 0 then begin

    ; Define ingrees / egress times
    cspice_UTC2ET, PenUmbra_UTC, PenUmbra_ET
    cspice_UTC2ET, Umbra_UTC, Umbra_ET

    ;----------------------Compare 2 methods to get a decent telluric absorption spectrum---------------------------

    ; Method 1: Assume absorption *is* a unity-normalized blue fast rotator stellar spectrum
    ;           Does pretty well for Oxygen 6300
    ;           Beware that interstellar sodium absorption will screw this up
      Telluric_Absorption_BFR = MRDFITS(Calibration_Dir+'\fullspec' +Standard_Star_files[0]+'.ec.fits', 0, header, /fscale, /silent, /unsigned )
      Telluric_WL_BFR         = sxpar(header, 'CRVAL1')+findgen(N_elements(Telluric_Absorption_BFR))*sxpar(header, 'Cdelt1')
      Telluric_AIRMASS_BFR    = sxpar(header, 'AIRMASS')

    ; Method 2: Divide an A0V stellar type model emission spectrum with a measured A0V
      Telluric_Absorption_A0V = MRDFITS(Calibration_Dir+'\sumfitcont' +Standard_Star_files[1]+'.ec.fits', 0, header, /fscale, /silent, /unsigned )
      Telluric_WL_A0V         = sxpar(header, 'CRVAL1')+findgen(N_elements(Telluric_Absorption_A0V))*sxpar(header, 'Cdelt1')
      Telluric_AIRMASS_A0V    = sxpar(header, 'AIRMASS')

    READCOL,'D:\DATA\___Calibration___\Solar_and_Stellar_Calibration_Spectra\vegallpr25.30000', F='A,A', Model_A0V_WL, Model_A0V_Flux, Skipline = 700000, numline = 500000, /Silent
    Model_A0V_WL = Model_A0V_WL*10.
    Model_A0V_Flux = float(shift(Model_A0V_Flux, -100))
    Model_A0V_Flux = INTERPOL(Model_A0V_Flux, Model_A0V_WL, Telluric_WL_A0V)

    result = ROBUST_POLY_FIT(Telluric_WL_A0V , Telluric_Absorption_A0V, 9)
    Echelle_Order_Interweave = Telluric_Absorption_A0V / poly(Telluric_WL_A0V, result)
    result = ROBUST_POLY_FIT(Telluric_WL_A0V, Model_A0V_flux, 6)
    Model_A0V_norm = Model_A0V_flux / poly(Telluric_WL_A0V, result)
    A0V_norm = MRDFITS(Calibration_Dir+'\Fullspec' +Standard_Star_files[1]+'.ec.fits', 0, header, /fscale, /silent, /unsigned )
    Telluric_Absorption_A0V = A0V_norm / Model_A0V_norm

    window, 0, title = 'Unity Normalized A0V Spectral Model (Blue) vs and Measured A0V (Black: HD 155379)'
    cgplot, Telluric_WL_A0V, Model_A0V_norm, color = 'blue', xr = [5570,9300]
    cgplot, Telluric_WL_A0V, A0V_norm, color = 'green', /overplot

    window, 1, title = 'Comparison of Telluric Absorption Spectra: Blue fast rotator (black) vs. A0V Stellar Model (Red)'
    cgplot, Telluric_WL_BFR, Telluric_Absorption_BFR, xr = [5000, 6270]
    cgplot, Telluric_WL_A0V, Telluric_Absorption_A0V, /overplot, color = 'red'
    cgplot, Telluric_WL_A0V, Echelle_Order_Interweave - 0.5, /overplot, color = 'green'

    ;----------------------Flux Calibrate Using Distance and Absolute Spectral Reflectivity of Jupiter---------------------------

    ; Compare publications of Jupiter's spectral albedo at disk center

    ; absolute brightness: Digitized Plot from Woodman et al. 1979.
      READCOL,'C:\IDL\Io\Woodman_et_al_1979_plot1.txt', F='A,A', WL, Albedo, STRINGSKIP = '#', /Silent;wavelength in Angstroms I/F unitless
      READCOL,'C:\IDL\Io\Woodman_et_al_1979_plot2_new.txt', F='A,A', WL_2, Albedo_2, STRINGSKIP = '#', /Silent ;wavelength in Angstroms I/F unitless
      Woodman_WL = float([WL, WL_2])              ; STITCH THESE TOGETHER
      Woodman_Albedo = Float([albedo , albedo_2]) ; STITCH THESE TOGETHER
  
      ; absolute brightness: from Karkoschka (1998) Icarus on the PDS as ID # ESO-J/S/N/U-SPECTROPHOTOMETER-4-V1.0
      READCOL,'C:\IDL\Io\Karkoschka_1995low.tab', F='X,A,X,A', Karkoschka_WL, Karkoschka_Albedo, STRINGSKIP = '#', /Silent ;wavelength in nm I/F unitless
      READCOL,'C:\IDL\Io\Karkoschka_1995high.tab', F='X,A,X,A', Karkoschka_WL_2, Karkoschka_Albedo_2, STRINGSKIP = '#', /Silent ;wavelength in Angstroms I/F unitless
      Karkoschka_WL = float([Karkoschka_WL, Karkoschka_WL_2])             ; STITCH THESE TOGETHER
      Karkoschka_Albedo = Float([Karkoschka_albedo, Karkoschka_albedo_2]) ; STITCH THESE TOGETHER
      Karkoschka_Albedo = Karkoschka_Albedo[sort(Karkoschka_WL)]
      Karkoschka_WL = Karkoschka_WL[sort(Karkoschka_WL)]

    ; compare the two
      cgplot, Woodman_WL / 10., Woodman_Albedo, color = 'blue', xstyle = 1., psym = 3, Xtitle = 'Wavelength (nm)', ytitle = 'I/F Reflectivity'
      cgplot, Karkoschka_WL, Karkoschka_Albedo*1.35, color = 'red', /overplot
      cgtext, 340, .1, 'EQUATOR AT CENTRAL MERIDIAN (Woodman et al. 1979)', color = 'blue'
      cgtext, 340, .16, 'FULL DISK scaled by 1.35 (Karkoschka 1998)', color = 'red'

    ; make an informed choice via scaling
      Karkoschka_wl = Karkoschka_wl * 10. ;nm to Angstroms
      Karkoschka_Albedo = Karkoschka_Albedo*1.35 ;Scale the "Full disk albedo" to the "Central Meridian Equatorial Absolute Reflectivity"

    ; Get the Solar Spectral Irradiance at 1 AU
      READCOL,'C:\IDL\Io\Kurucz_2005_irradthuwl.dat', F='A,A', WL_nm, flux, STRINGSKIP = '#', /Silent ;flux is in W/m2/nm
      start  = where(WL_nm eq '299.100')
      WL_nm = float(WL_nm[start:*])
      flux = float(flux[start:*])

    ; change flux units from W/m^2/nm to photons/(cm^2 s A)
    ; multiply by ((lambda / hc) / 1 W)  * (1 m^2 / 1e4 cm^2) * (1 nm / 10 A)
      conversion = ((WL_nm*1.e-9)/(6.62606957e-34*299792458.D)) * (1./1.e4) * (1./10.)
      flux = flux * conversion      ; Cross-checked this result against Huebner et al. 1992.
      WL_A = temporary(WL_nm) * 10. ; Wavelength from nm into angstroms
      WL_A = shift(WL_A, 1)         ; Evidently this is needed to line up with e.g., Bass 2000
      WL_A = WL_A[1:-1]
      flux = flux[1:-1]             ; trim things to clean up any edge effect due to the wavelength shift
      WL_Vacuum = WL_A              ; Vacuum wavelength in A for the solar spectrum (need vaccum for g-value calculations)
      VACTOAIR, WL_A, WL_A_Air      ; Vacuum to air wavelength conversion      
      WL_A = temporary(WL_A_Air)      

    ; Scale the solar flux to Jupiter's instantaneous distance
      jup_center           = MRDFITS(reduced_dir+'fullspec'+Jupiter_Center_File+'.ec.fits', 0, Jupiter_center_header, /fscale, /silent, /unsigned )     ; fullspec, normalized to one
      jup_center_err       = MRDFITS(reduced_dir+'sig_fullspec'+Jupiter_Center_File+'.ec.fits', 0, Jupiter_center_err_header, /fscale, /silent, /unsigned ) ; 1 sig uncertainty in fullspec
      WL                   = sxpar(Jupiter_center_header, 'CRVAL1')+findgen(N_elements(jup_center))*sxpar(Jupiter_center_header, 'Cdelt1')
      cspice_UTC2ET, sxpar(Jupiter_center_header, 'DATE-OBS'), ET
      cspice_spkezr, 'Jupiter', ET, 'J2000', 'LT+S', 'Sun', Jupiter_Sun_State, ltime
      cspice_spkezr, 'Jupiter', ET, 'J2000', 'LT+S', 'Earth', Jupiter_Earth_State, ltime
      solar_distance = norm(Jupiter_Sun_State[0:2]) / 149597871.
      flux_at_jupiter = flux / solar_distance^2.

    ; Multiply incident solar irradiance x spectral albedo to determine the theoreitcal brightness of Jupiter at disk center
      Albedo = INTERPOL(Karkoschka_Albedo, Karkoschka_WL, WL_A)
      Rayleighs_per_angstrom = 4.*flux_at_jupiter*albedo / 1.e6
      if keyword_set(debug) then window, Title = 'Instantaneous Rayleighs per Angstrom: Center of Jupiter''s Disk'
      if keyword_set(debug) then plot, WL_A, Rayleighs_per_angstrom, xr = [5885, 5900], charsize = 2 ;compare to 5.5 MR per angstrom by Brown & Schneider, 1981

    ; Adjust the expected absolute flux for Jupiter's Instantaneous Doppler shift
      theta  = cspice_vsep(Jupiter_Earth_State[0:2], Jupiter_Earth_State[3:5])
      Dopplershift = cos(theta) * norm(Jupiter_Earth_State[3:5])   ; scalar projection of the relative velocity along the line of sight
      WL_A = WL_A + WL_A*Dopplershift/cspice_clight()

    ; ***********Correct the Unity-Normalized Jupiter spectrum for Telluric Absorption**********

    ; Make a master telluric correction that combines the various spectral ranges used and the two different stars A0V, and BFR
      Master_telluric = replicate(1., N_elements(jup_center))

    ; Try both the BFR and the A0V star methods.
      rough_absorption_A0V = interpol(Telluric_Absorption_A0V, Telluric_WL_A0V, WL) ; Using A0V star correction instead of BFR
      rough_absorption_BFR = interpol(Telluric_Absorption_BFR, Telluric_WL_BFR, WL) ; Using BFR star correction instead of AOV

      Io_header_First = headfits(reduced_dir+'fullspec'+Eclipse_files[0]+'.ec.fits', /silent)
      Io_header_Last  = headfits(reduced_dir+'fullspec'+Eclipse_files[-1]+'.ec.fits', /silent)
      airmass_range   = float( [sxpar(Io_header_First, 'AIRMASS'), sxpar(Io_header_last, 'AIRMASS') ] )
      print, 'Data Span an airmass range:', airmass_range

    ; Where telluric absorption corrections are needed, carefully identify wavelength indicies where we can isolate the absorption
      for i = 0, N_elements(Io_Airglow)-1 do begin
        telluric_absorption = 0            ; let's default to no correction for Earth's atmosphere, and handle corrections case by case with wavelength
        Case 1 of
          (Io_Airglow[i] eq  Na1) or (Io_Airglow[i] eq  Na2): begin                                                          ; correct telluric absorption only a high airmass?


            if Date eq '' then begin ; INSERT DATE HERE FOR WHICH DATES need to be corrected. Part 0 needsto be rerun ebfore Part 1 for telluric corrections to take effect


              telluric_absorption = 1
              spectral_range = [5885, 5901]  ; Wavelength region for inspection. 
              telluric_fitting_ind = where( (wl gt 5891.) and (wl lt 5895.), /NULL)                                          ; dominated by telluric lines if airmass is high enough
            endif
          end   
          (Io_Airglow[i] eq  Na3) or (Io_Airglow[i] eq  Na4): begin
            telluric_absorption = 0
            spectral_range = [8173, 8204]  ; Unfortunately the standard star lines here begin to saturate (e.g. 8189A), making W's in the corrected spectra
            telluric_fitting_ind = where( (wl gt 8176.) and (wl lt 8182.) or ((wl gt 8185.) and (wl lt 8194.)), /NULL)       ; dominated by telluric lines
            ; Can't fix this without better standard star data
          end
          (Io_Airglow[i] eq  K1): begin
            telluric_absorption = 1
            spectral_range = [7650, 7680]  ; Wavelength region for inspection
            telluric_fitting_ind = where( ((wl gt 7665.) and (wl lt 7667.)) or ((wl gt 7670.) and (wl lt 7673.)), /NULL)     ; dominated by telluric lines
          end
          (Io_Airglow[i] eq  K2): begin
            telluric_absorption = 1
            spectral_range = [7680, 7720]  ; Wavelength region for inspection
            telluric_fitting_ind = where( (wl gt 7695.5) and (wl lt 7697.5), /NULL)                                          ; dominated by telluric lines
            ;telluric_fitting_ind = where( ((wl gt 7659.) and (wl lt 7661.)) or ((wl gt 7670.) and (wl lt 7673.)), /NULL)    ; dominated by telluric lines
          end
          (Io_Airglow[i] eq O1):                                                                                             ; no telluric absorption
          (Io_Airglow[i] eq O2): begin
            telluric_absorption = 1
            spectral_range = [6292, 6308]  ; Wavelength region for inspection
            telluric_fitting_ind = where( ((wl gt 6295.0) and (wl lt 6296.5)) or ((wl gt 6305.5) and (wl lt 6307.0)), /NULL) ; dominated by telluric lines
          end
          (Io_Airglow[i] eq O3):                                                                                             ; no tellurics noticable near 6364A
          else:                                                                                                              ; Certainly some tellurics but just too messy to remove them!
        endcase
        
      if telluric_absorption eq 0 then continue ; otherwise fix the tellurics

      ; adjust absorption depths to match Jupiter
        Rough_Absorption_A0V = Rough_Absorption_A0V^( float(sxpar(Jupiter_center_header, 'AIRMASS')) / Telluric_AIRMASS_A0V )
        Rough_Absorption_BFR = Rough_Absorption_BFR^( float(sxpar(Jupiter_center_header, 'AIRMASS')) / Telluric_AIRMASS_BFR )

      ; Fix any subtle pixel shifts between *known telluric lines* measured in the calibration stars and in Jupiter
        lag                = findgen(41)-20.               ; search range in pixels, at what shift the maximal correlation of target spectrum and telluric spectrum occur? 
        correl_BFR         = C_CORRELATE( Rough_Absorption_BFR[telluric_fitting_ind], jup_center[telluric_fitting_ind], lag)
        correl_A0V         = C_CORRELATE( Rough_Absorption_A0V[telluric_fitting_ind], jup_center[telluric_fitting_ind], lag)
        yfit_BFR           = MPFITPEAK(lag, correl_BFR, A_BFR, nterms = 3, /positive)
        yfit_A0V           = MPFITPEAK(lag, correl_A0V, A_A0V, nterms = 3, /positive)
        aligned_absorp_BFR = interpolate(Rough_Absorption_BFR, findgen(n_elements(Rough_Absorption_BFR)) - A_BFR[1])
        aligned_absorp_A0V = interpolate(Rough_Absorption_A0V, findgen(n_elements(Rough_Absorption_A0V)) - A_A0V[1])

      ; Inspect the alignment of both methods of telluric absorption and Jupiter
        window, 0, title = 'Checking Wavelength Alignment of the Telluric Absorption: Bottom plot looks aligned after the above shift? (Blue = BFR, Orange = A0V) ', ys = 800
        pos = cglayout([1,2])
        cgplot, lag, correl_BFR, pos = pos[*,0], title = 'Correlation coefficient vs shift between target and telluric spectrum', xtitle = 'shift in pixels', ytitle = 'correlation'
        cgplot, lag, yfit_BFR, /overplot, color = 'blue'
        cgplot, lag, correl_A0V, /overplot
        cgplot, lag, yfit_A0V, /overplot, color = 'orange'
  
        cgplot, WL, jup_center, xr = spectral_range, /ynozero, pos = pos[*,1], /noerase
        cgplot, WL, aligned_absorp_BFR, color = 'blue', /overplot
        cgplot, WL, aligned_absorp_A0V, color = 'orange', /overplot

      print, 'USER---> Carefully inspect the Wavelength Alignment of the Telluric Absorption Correction at Jupiter in every Bandpass:', spectral_range

      ; Fit a y = Ax^C + B function to the spectrum, where x is the telluric absorption scaled for 0 to 1 dynamic range, and y is jupiter, scaled for 0 to 1 dynamic range
        p0 = [1.0, 0.0, 1.0] ; guess at initial coefficients
        p = mpfitfun('Match_Telluric_Absorption', aligned_absorp_BFR[telluric_fitting_ind], jup_center[telluric_fitting_ind], 0.1*jup_center_err[telluric_fitting_ind], $
          p0, /NaN, status=status, /quiet)
        telluric_fit_BFR = P[0]*(aligned_absorp_BFR^P[2]) + P[1]
  
        p = mpfitfun('Match_Telluric_Absorption', aligned_absorp_A0V[telluric_fitting_ind], jup_center[telluric_fitting_ind], 0.1*jup_center_err[telluric_fitting_ind], $
          p0, /NaN, status=status, /quiet)
        telluric_fit_A0V = P[0]*(aligned_absorp_A0V^P[2]) + P[1]

      ; Inspect the telluirc fit
        window, 2, xs = 1000, ys = 900, title='Black = Jupiter, Green = Indices where Telluric Absorptions are scaled, Red = Aligned Telluric Absorption, Blue = Scaled absorption to Match Jupiter'
        pos = cgLayout([1,2], OXMargin=[11,3], OYMargin=[9,6], YGap=0)
        sample = where((WL gt spectral_range[0]) and (WL lt spectral_range[1]), /NULL)
  
        yr = minmax(jup_center[sample])
        cgplot, WL, jup_center, xr = spectral_range, yr = [yr[0], 1.1], pos = pos[*,0], xtickformat = '(A1)'
        cgplot, WL[telluric_fitting_ind], jup_center[telluric_fitting_ind], /overplot, color = 'green', psym=14
        cgplot, WL, aligned_absorp_A0V, color = 'red', /overplot
        cgplot, WL, telluric_fit_A0V, color = 'blue', /overplot
        cgtext, spectral_range[0]+1, yr[0]+0.15, 'Jupiter Center FullSpec', charsize = 2
        cgtext, spectral_range[0]+1, yr[0]+0.1,  'Aligned A0V', color = 'red', charsize = 2
        cgtext, spectral_range[0]+1, yr[0]+0.05,  'Fit A0V', color = 'Blue', charsize = 2
        cgtext, spectral_range[1]-10, yr[0]+0.05,  'Fitting indicies', color = 'green', charsize = 2
  
        yr = minmax(jup_center[sample]/telluric_fit_BFR[sample])
        cgplot, WL, jup_center/telluric_fit_BFR, xr = spectral_range, yr = yr, /noerase, pos = pos[*,1], color = 'blue'
        cgplot, WL, jup_center/telluric_fit_A0V, /overplot, color = 'orange'

      ; Keep whichever telluric correction method is closer to unity at every spectral bin
      ; Really you should keep whichever is a better match to Jupiter's reflectance * solar spectrum at every bin, but that's a pain to code up
      ; Interstellar Na absorption is the big problem
        residual_BFR  = abs(1.-jup_center/telluric_fit_BFR)
        residual_A0V  = abs(1.-jup_center/telluric_fit_A0V)
        
          ; For sodium, just use whichever telurc correction gives the lower value.
          ; A correction with Interstellar Na will give high values, so default to lower
          ; Otherwise use whichever telluric spectrum gives closer to unity result
          if ((Io_Airglow[i] eq  Na1) or (Io_Airglow[i] eq  Na2)) then begin 
            junk = min([[jup_center/telluric_fit_BFR],[jup_center/telluric_fit_A0V]], dim = 2, loc)
          endif else junk = min([[residual_BFR],[residual_A0V]], dim = 2, loc) 
          
        combined      = [[jup_center/telluric_fit_BFR], [jup_center/telluric_fit_A0V]]
        tell_combined = [[telluric_fit_BFR], [telluric_fit_A0V]]
        telluric_fit  = tell_combined[loc]
        jup_tell_corr = combined[loc]
        cgplot, WL, jup_tell_corr, /overplot, color = 'red'
        cgtext, spectral_range[0]+1, 0.90, 'Correction using BFR', color = 'blue'
        cgtext, spectral_range[0]+1, 0.85, 'Correction using A0V', color = 'orange'
        cgtext, spectral_range[0]+1, 0.80, 'Correction using a combination', color = 'red'

      ; write in this region of the telluric absorption
        Master_telluric[sample] = telluric_fit[sample]
      
      ;if (Io_Airglow[i] eq  Na1) then stop
    endfor

    ; Apply master telluric absorption correction to the Jupiter center spectrum.
      jup_tell_corr = jup_center / Master_telluric

    ; We've been working with unity-normalized spectra until now. Determine the DN / S from the "Sum" files, which are real DN, not unity-normalized for continuum
      jup_sum        = MRDFITS(reduced_dir+'sum'+Jupiter_Center_File+'.ec.fits', 0, header, /fscale, /silent, /unsigned )
      raw_header     = headfits(dir+''+Jupiter_Center_File+'.fits')
      jup_sum        = jup_sum / sxpar(raw_header, 'EXPTIME')      ; Convert to DN per second
      jup_sum        = jup_sum / Master_telluric                   ; Remove telluric absorption
      jup_Sum_Curve  = smooth(jup_Sum, 2500, /edge_truncate, /nan) ; Smooth the hell out of it

    ; Determine the instrumental sensitivity from the expected versus measured flux at Jupiter Disk Center. This should be a smooth function
      expected_flux = interpol(Rayleighs_per_angstrom, WL_A, WL)   ; move expected flux to the Jovian Doppler-shift, UNITS are R / A
      smoothed_expected_flux = GAUSS_SMOOTH(expected_flux, 1.9)    ; this smoothing looks about right for APO/ARCES, UNITS are R / A
      window, 3
      Sensitivity       = jup_sum / smoothed_expected_flux         ; Sensitivity in (DN / S) / (R / A)
      Sensitivity_Curve = smooth(sensitivity, 5000, /edge_truncate, /nan)
      cgplot, WL, Sensitivity, /xstyle, yr = [0.,0.0025], Ytitle = 'Measured Flux Sensitivity (DN/S) / (R/A)', Xtitle = 'Angstroms'
      fit_Sensitivity   = Sensitivity
      fit_sensitivity[where( (sensitivity gt 0.002) or (sensitivity lt 0.0), /Null)] = !values.F_NaN                                            ; some basic rejection criterion
      Sensitivity_Curve = smooth(fit_sensitivity, 5000, /edge_truncate, /nan)
      fit_sensitivity[where( (Sensitivity_Curve/fit_Sensitivity gt 1.5) or (Sensitivity_Curve/fit_Sensitivity lt 0.5), /Null)] = !values.F_NaN  ; further rejection criterion
      Sensitivity_Curve = smooth(fit_sensitivity, 5000, /edge_truncate, /nan)
      cgplot, WL, Sensitivity_Curve, color = 'red', /overplot

    ; Inspect and write the telluric corrected R/A Calibrated Jupiter disk center spectra to file
      window, 4, Title = 'Inspect Jupiter: Red = Uncorrected, Black = Telluric Corrected'
      Jup_Cal_Tell_Corr       = jup_tell_corr * jup_Sum_Curve / Sensitivity_Curve  ; Telluric Correct and Convert to Rayleighs per Angstrom
      Jup_Cal_No_Tell_Corr    = jup_center * jup_Sum_Curve / Sensitivity_Curve     ; Convert to Rayleighs per Angstrom, no correction
      cgplot, WL, Jup_Cal_Tell_Corr, xr = spectral_range, /ynozero
      cgplot, WL, Jup_Cal_No_Tell_Corr, color = 'red', /overplot
      MWRFITS, Jup_Cal_Tell_Corr, reduced_dir+'R_per_A_Jupiter_Center.ec.fits', header, /CREATE ;/create overwrites
      MWRFITS, SQRT(abs(Jup_Cal_Tell_Corr)), reduced_dir+'sig_R_per_A_Jupiter_Center.ec.fits', sig_header, /CREATE, /silent ; Assume Poisson sig values

  ; **************************** Io Calibration and Absorption Correction **********************************************+
    cspice_UTC2ET, PenUmbra_UTC, PenUmbra_ET
    cspice_UTC2ET, Umbra_UTC, Umbra_ET

  ; Correct tellurics and write the spectra into R/A units...
    for i = 0, n_elements(Eclipse_files)-1 do begin

      Io_ecl            = MRDFITS(reduced_dir+'fullspec'+Eclipse_files[i]+'.ec.fits', 0, Io_header, /fscale, /unsigned, /silent)
      Io_ecl_tell_corr  = Io_ecl / Master_telluric^( float(sxpar(Io_header, 'AIRMASS')) / float(sxpar(Jupiter_center_header, 'AIRMASS')) )
      Io_ecl_sum        = MRDFITS(reduced_dir+'sum'+Eclipse_files[i]+'.ec.fits', 0, Io_header, /fscale, /silent, /unsigned )
      raw_header        = headfits(dir+Eclipse_files[i]+'.fits')
      Io_ecl_sum        = Io_ecl_sum / sxpar(raw_header, 'EXPTIME')      ; Convert to DN per second
      Io_ecl_sum        = Io_ecl_sum / Master_telluric^( float(sxpar(Io_header, 'AIRMASS')) / float(sxpar(Jupiter_center_header, 'AIRMASS')) )
      Io_ecl_Sum_Curve  = smooth(Io_ecl_Sum, 2500, /edge_truncate, /nan)
      Io_ecl_Cal_Tell_Corr       = Io_ecl_tell_corr * Io_ecl_Sum_Curve / Sensitivity_Curve  ; Convert to Rayleighs per Angstrom, Telluric Corrected
      Io_ecl_Cal_No_Tell_Corr    = Io_ecl * Io_ecl_Sum_Curve / Sensitivity_Curve            ; Convert to Rayleighs per Angstrom, no Telluric correction

      ; correct for broad-band extinction between the Jupiter airmass and the Io airmass (this is crude, but better than nothing!)
        K = 0.16 ; APO extinction coefficient in mag / airmass in greenish redish wavelengths (cf. Hogg et al. 2001 Astronomical Journal)
        Scale_for_extinction    = 100. ^ (K * (float(sxpar(Io_header, 'AIRMASS')) - float(sxpar(Jupiter_center_header, 'AIRMASS'))) / 5.)
        Io_ecl_Cal_No_Tell_Corr = Io_ecl_Cal_No_Tell_Corr * Scale_for_extinction
        Io_ecl_Cal_Tell_Corr    = Io_ecl_Cal_Tell_Corr * Scale_for_extinction
        print, 'Io Airmass:', sxpar(Io_header, 'AIRMASS'), ' Jupiter Airmass:', sxpar(Jupiter_center_header, 'AIRMASS'), ' Scaling by:', Scale_for_extinction

      ;      ; Penumbral emissions also need Io reflactance subtracted
      ;        if ((date eq 'UT180320') and (i eq 0)) then begin
      ;          window, 0
      ;          cgplot, WL, Io_ecl_Cal_Tell_Corr, xr = xr, /ynozero
      ;          Io_sunlit = MRDFITS(reduced_dir+'fullspecIo_free_and_clear.0006.ec.fits', 0, Io_header, /fscale, /silent, /unsigned )     ; fullspec, normalized to one
      ;          Io_sunlit = Io_sunlit * ARCES_Correct_Interweave(Io_sunlit, WL)
      ;          Io_sunlit = Io_sunlit / Master_telluric^( float(sxpar(Io_header, 'AIRMASS')) / float(sxpar(Jupiter_center_header, 'AIRMASS')) )
      ;          Io_sunlit = interpolate(Io_sunlit, findgen(n_elements(Io_sunlit)) - .3)
      ;          scale_Io_reflectance = 8.2e4
      ;          ;scale_Io_reflectance = 9.2e4
      ;          cgplot, WL, Io_sunlit*scale_Io_reflectance, /overplot, linestyle = 2
      ;          Io_ecl_Cal_Tell_Corr = Io_ecl_Cal_Tell_Corr - Io_sunlit*scale_Io_reflectance
      ;          cgplot, WL, Io_ecl_Cal_Tell_Corr*2.+8.e4, /overplot, color = 'red'
      ;        endif

      ; find the instantaneous Earth-Io Doppler Shift
        cspice_UTC2ET, sxpar(raw_header, 'DATE-OBS'), ET
        ET_mid_exposure = ET + float(sxpar(raw_header, 'EXPTIME'))/2.
        cspice_spkezr, 'Io', ET_mid_exposure, 'J2000', 'LT+S', 'Earth', Io_Earth_State, ltime        
        ;cspice_spkezr, 'Ganymede', ET_mid_exposure, 'J2000', 'LT+S', 'Earth', Io_Earth_State, ltime

      ; Sodium D1+D2 g-value
        cspice_spkezr, 'Io', ET_mid_exposure - Ltime, 'J2000', 'LT+S', 'Sun', Io_Sun_State, Io_Sun_ltime  
        theta  = cspice_vsep(Io_Sun_state[0:2], Io_Sun_state[3:5])
        Io_heliocentric_vel = cos(theta) * norm(Io_Sun_State[3:5])
        GVALUE, 'Na-D', Io_heliocentric_vel*1.e3, norm(Io_Sun_state[0:2]) / 1.495927e8, WL_Vacuum, flux, g
        SXADDPAR, Io_header, 'Na_G_Val', g, 'D2 + D1 photons/atom/s'

        theta  = cspice_vsep(Io_Earth_state[0:2], Io_Earth_state[3:5])
        Io_wrt_Earth_Dopplershift = cos(theta) * norm(Io_Earth_State[3:5])
        SXADDPAR, Io_header, 'Io_DOPPL', Io_wrt_Earth_Dopplershift, 'Io-Earth V_radial in km/s (mid exposure)'
        SXADDPAR, Io_header, 'T_P_Shad', (ET_mid_exposure-PenUmbra_ET) / 60., 'Minutes since Penumbral ingress'
        SXADDPAR, Io_header, 'T_USHADO', (ET_mid_exposure-Umbra_ET) / 60., 'Minutes since Umbral ingress'

      ; Find the system III longitude and latitude of Io, and write them to the header
        cspice_subpnt, 'Near point: ellipsoid', 'Jupiter', ET_mid_exposure - ltime, 'IAU_Jupiter', 'None', 'Io', Sub_Io, trgepc, srfvec ;get sub solar point in cartesian IAU Jupiter coords
        cspice_bodvrd, 'Jupiter', 'RADII', 3, radii ;get Jupiter shape constants
        re = radii[0]
        rp = radii[2]
        f = (re-rp)/re
        obspos = Sub_Io - srfvec
        cspice_recpgr, 'Jupiter', obspos, re, f, Io_SysIII, Io_SysIII_LATITUDE, opgalt
        ;torus_lat = (2./3.) * 10.31*cos( (196.61 / !radeg) - Io_SysIII )                                  ; Depreciated JRM09 Dipole approximation
        torus_lat = interpol(reverse(torus_deg), torus_lat_out, Io_SysIII*!radeg)
        SXADDPAR, Io_header, 'Sys3_Lat', Io_SysIII, 'Io''s System III Longitude'
        SXADDPAR, Io_header, 'Sys3_Lon', Io_SysIII_LATITUDE, 'Io''s System III Latitude'
        SXADDPAR, Io_header, 'Torus_Lat', torus_lat, 'JRM09 Dipole Approx'

      ; Now do the slit filling factor aperture correction. ***************** ASSUMES EMISSION REGION IS THE SIZE OF IO'S DISK *************
        ARCES_Slit_area = 1.6*3.2                                                                        ; Default ARCES slit size
        Io_Area = !pi * (tan(1821.6 / norm(Io_Earth_State[0:2])) * 206265.)^2                            ; Io's solid angle in square arcseconds
        Io_ecl_Cal_Tell_Corr = Io_ecl_Cal_Tell_Corr * ARCES_Slit_area / Io_Area

      MWRFITS, Io_ecl_Cal_Tell_Corr, reduced_dir+'R_per_A_'+Eclipse_files[i]+'.ec.fits', Io_header, /CREATE ; /create overwrites
      MWRFITS, sqrt(abs(Io_ecl_Cal_Tell_Corr)), reduced_dir+'sig_R_per_A_'+Eclipse_files[i]+'.ec.fits', Io_header, /CREATE ; ASSUME POISSON STATISTICS APPLY, IGNORE DETAILED UNCERTAINTY
      ;MWRFITS, Io_ecl_Cal_No_Tell_Corr, reduced_dir+'R_per_A_No_Tell_Corr'+Eclipse_files[i]+'.ec.fits', Io_header, /CREATE ; /create overwrites
    endfor

    ; **************************** Jovian Scatter (Sky Background) Calibration and Absorption Correction **********************************************+
    for i = 0, n_elements(Jovian_Scatter_files)-1 do begin

      Jovian_Scatter            = MRDFITS(reduced_dir+'fullspec'+Jovian_Scatter_files[i]+'.ec.fits', 0, Jovian_Scatter_header, /fscale, /unsigned, /silent)
      Jovian_Scatter_tell_corr  = Jovian_Scatter / Master_telluric^( float(sxpar(Jovian_Scatter_header, 'AIRMASS')) / float(sxpar(Jupiter_center_header, 'AIRMASS')) )
      Jovian_Scatter_sum        = MRDFITS(reduced_dir+'sum'+Jovian_Scatter_files[i]+'.ec.fits', 0, Jovian_Scatter_header, /fscale, /silent, /unsigned )
      raw_header                = headfits(dir+Eclipse_files[i]+'.fits')
      Jovian_Scatter_sum        = Jovian_Scatter_sum / sxpar(raw_header, 'EXPTIME')                                 ; Convert to DN per second
      Jovian_Scatter_sum        = Jovian_Scatter_sum / Master_telluric^( float(sxpar(Jovian_Scatter_header, 'AIRMASS')) / float(sxpar(Jupiter_center_header, 'AIRMASS')) )
      Jovian_Scatter_Sum_Curve  = smooth(Jovian_Scatter_Sum, 2500, /edge_truncate, /nan)

      Jovian_Scatter_Cal_Tell_Corr       = Jovian_Scatter_tell_corr * Jovian_Scatter_Sum_Curve / Sensitivity_Curve  ; Convert to Rayleighs per Angstrom, Telluric Corrected
      Jovian_Scatter_Cal_No_Tell_Corr    = Jovian_Scatter * Jovian_Scatter_Sum_Curve / Sensitivity_Curve            ; Convert to Rayleighs per Angstrom, no Telluric correction

      ; Now do the same slit filling factor aperture correction for consisteancy with the Io Eclipse frames
        ARCES_Slit_area = 1.6*3.2 ; Default ARCES slit size
        Io_Area = !pi * (tan(1821.6 / norm(Io_Earth_State[0:2])) * 206265.)^2                                       ; Io's solid angle in square arcseconds
        Jovian_Scatter_Cal_Tell_Corr = Jovian_Scatter_Cal_Tell_Corr * ARCES_Slit_area / Io_Area
        Jovian_Scatter_Cal_No_Tell_Corr = Jovian_Scatter_Cal_No_Tell_Corr * ARCES_Slit_area / Io_Area

      ; Jovian Scatter frames are very long exposures with substantial telluric airglow.
      ; Interpolate over the telluric airglow lines, so that Jovian scatter is "just jupiter"
        ARCES_2sigma_bandwidth = 2.*(Telluric_Airglow/31500.) / 2.3548 ; 2 * FWHM in Angstroms / 2sqrt(2ln2) = 99% of flux enclosed within 2 sigma --> Appropriate for bright emissions
        Telluric_AG_Ind = []
        for j = 0, N_elements(Telluric_Airglow)-1 do Telluric_AG_ind = [Telluric_AG_Ind, where(abs(Telluric_Airglow[j] - WL) lt ARCES_2sigma_bandwidth[j], /NULL)] ; indicies within 3 sigma of any Telluric airglow
        Keep_ind   = cgsetdifference( indgen(N_elements(Jovian_Scatter_Cal_Tell_Corr)), Telluric_AG_ind )
        dummy      = Jovian_Scatter_Cal_Tell_Corr
        dummy[Telluric_AG_ind] = !values.F_NaN
        Keep_ind   = where(finite(dummy), /null)
        No_Airglow = INTERPOL( Jovian_Scatter_Cal_Tell_Corr[Keep_ind], WL[Keep_ind], WL)
        cgplot, WL, Jovian_Scatter_Cal_Tell_Corr, /ynozero, xr = [5575., 5579.]
        cgplot, WL, No_airglow, /overplot, color = 'blue'
        Jovian_Scatter_Cal_Tell_Corr = temporary(No_airglow)

      ;cgplot, WL, Jovian_Scatter_Cal_Tell_Corr, /ynozero, xr = [8180, 8197]
      ;cgplot, WL, Jovian_Scatter_Cal_No_Tell_Corr, /overplot, color = 'red'

      MWRFITS, Jovian_Scatter_Cal_Tell_Corr, reduced_dir+'R_per_A_'+Jovian_Scatter_files[i]+'.ec.fits', Jovian_Scatter_header, /CREATE ; /create overwrites
      ;MWRFITS, Jovian_Scatter_Cal_No_Tell_Corr, reduced_dir+'R_per_A_No_Tell_Corr'+Jovian_Scatter_files[i]+'.ec.fits', Jovian_Scatter_header, /CREATE ; /create overwrites
    endfor

    ; write the telluric spectrum since it may be useful for sebsequent analysis
      MWRFITS, Master_telluric, reduced_dir + 'master_telluric_spectrum.fits', /CREATE
  endif ; Part eq 0 ARCES basic redux

  if Part eq 1.00 then begin ; Make waterfall plots

    ;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ;--------------------------------------------------------ARCES Scatter Fit & "Waterfall" Plots----------------------------------------------------------------------------
    ;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ; METHOD:
    ;   Loop though all possible spectra that could be fit to jovian scatter use the minimum residual to determine what works best. Avoid Io and telluric airglow in the fitting.
    ;   Three steps to the fit: 1) Use Mpfit to roughly add & multiply the Jovian scatter spectra to approach the background in the Io exposures
    ;                           2) Use Amoeba to shift and smooth the spectra, maximizing correlation. If issues are encoutered iterate using bigger tolerances.
    ;                              Amoeba is prone to infinite loops if two points in the simplex give the same correlation, and iteration helps get an answer.
    ;                              Since results are somewhat senstive to the guess, run Amoeba a seccond time to assure that the initial simplex covers appropriate territory.
    ;                           3) Use Mpfit again to fine tune the add & multiply of Jovian scatter spectra needed to match the background in the Io exposures.
    ;                              This time, weight the fit by (1/N_pixels)^power, adjusting the power to optimize things visually.
    ;
    ;  Then run the fit again and generate a postscript, this time pegging some options. Use a goodness of fit metric to *MANUALLY* determine which is the best jovian scatter spectrum.
    ;  Use only one *MANUALLY* determined background spectrum throughout. *PEG* the linewidth in the Gaussian fitting.
    ;  Find any systematic shift from expected wavelengths, slightly offset the specctra, and *PEG* the line center in Gaussian fitting.

    ; NOTES:

    ;   It's critical to inspect with plot, WL[correl_indicies], scatter_fit[correl_indicies], uncomment this line to inspect and adjust things
    ;   With this inspection, it is useful to adjust *include_WLs*, and *AG_ind* to cover regions of continuum as close as possible to Io's airglow
    ;   Parameter information is not included in the Gaussian fitting to the line spread function in the residual. I couldn't figgure out why it failed with status = 0

    ; Fit various spectra to the Jovian Scatter, *use* whichever one produces the lowest residual. Note this method will never work if there's much telluric absorption.
      order_Na_D = {WL_range:[5888.,5897.5], waterfall_plot_title:'Io''s Airglow in Eclipse on ' + date, npanels:2, name:'order_Na_D'}
      order_NaIR = {WL_range:[8180.5, 8199.], waterfall_plot_title:'Io''s Airglow in Eclipse on ' + date, npanels:3, name:'order_NaIR'}
      order_7774 = {WL_range:[7770.5, 7777.], waterfall_plot_title:'Io''s Airglow in Eclipse on ' + date, npanels:3, name:'order_7774'}
      order_K_D  = {WL_range:[7660., 7704.], waterfall_plot_title:'Io''s Airglow in Eclipse on ' + date, npanels:3, name:'order_K_D'}
      order_8446 = {WL_range:[8444., 8448.], waterfall_plot_title:'Io''s Airglow in Eclipse on ' + date, npanels:3, name:'order_8446'}
      order_9225 = {WL_range:[9208., 9242.], waterfall_plot_title:'Io''s Airglow in Eclipse on ' + date, npanels:3, name:'order_9225'}
      order_5577 = {WL_range:[5574., 5580.], waterfall_plot_title:'Io''s Airglow in Eclipse on ' + date, npanels:3, name:'order_5577'}
      order_6300 = {WL_range:[6297.6, 6302.6], waterfall_plot_title:'Io''s Airglow in Eclipse on ' + date, npanels:2, name:'order_6300'}
      order_6364 = {WL_range:[6361.53, 6366.03], waterfall_plot_title:'Io''s Airglow in Eclipse on ' + date, npanels:2, name:'order_6364'}
      if keyword_set(ingress) then order_Na_D.waterfall_plot_title = 'Io''s Airglow Response Following ' + date + ' Ingress'
      if keyword_set(ingress) then order_6300.waterfall_plot_title = 'Io''s Airglow Response Following ' + date + ' Ingress'
    
<<<<<<< Updated upstream
    ;orders     = [order_5577, order_6300, order_6364, order_Na_D] ; which lines to extract
    orders     = [order_6300, order_6364, order_Na_D, order_NaIR, order_9225, order_8446, order_7774, order_K_D] ;use this to plot ALL the wavelength regions with potential emission
=======
    orders     = [order_Na_D] ; which lines to extract
    ;orders     = [order_6300, order_6364, order_Na_D, order_NaIR, order_9225, order_8446, order_7774, order_K_D] ;use this to find ALL the regions
>>>>>>> Stashed changes

    MX_plus_B_parinfo          = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0]}, 2)
    MX_plus_B_parinfo[1].fixed = 1        ; peg the additive "B" component of the MX_Plus_B at zero

    ; Define ingrees / egress times
    cspice_UTC2ET, PenUmbra_UTC, PenUmbra_ET
    cspice_UTC2ET, Umbra_UTC, Umbra_ET
    Io_Airglow_params.Umbra_ET = Umbra_ET
    Io_Airglow_params.PenUmbra_ET = PenUmbra_ET

    Case date of
      'UT180320': begin
        Fit_Me_To_The_Scatter    = ['D:\DATA\Apache Point\Echelle\Io Eclipses\UT180320\Reduced\R_per_A_Jovian_Scatter.0002.ec.fits', $
                                    'D:\DATA\Apache Point\Echelle\Io Eclipses\UT180320\Reduced\R_per_A_Jovian_Scatter.0003.ec.fits' ]
      end
      'UT190812': begin
        if sunlit eq 0 then Fit_Me_To_The_Scatter    = ['D:\DATA\Apache Point\Echelle\Io Eclipses\UT190812\Reduced\R_per_A_Jupiter_Center.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200823\Reduced\R_per_A_Jupiter_Scatter.0001.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200823\Reduced\R_per_A_Jupiter_Scatter.0002.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200823\Reduced\R_per_A_Jupiter_Scatter.0003.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200823\Reduced\R_per_A_Jupiter_Scatter.0004.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200908\Reduced\R_per_A_Jupiter_Scatter.0017.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200908\Reduced\R_per_A_Jupiter_Scatter.0030.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200908\Reduced\R_per_A_Jupiter_Scatter.0032.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT201001\Reduced\R_per_A_Jupiter_Scatter.0001.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT201001\Reduced\R_per_A_Jupiter_Scatter.0002.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT201001\Reduced\R_per_A_Jupiter_Scatter.0003.ec.fits' ]
        if sunlit eq 1 then Fit_Me_To_The_Scatter    = ['D:\DATA\Apache Point\Echelle\Io Eclipses\UT190812\Reduced\R_per_A_Ganymede.0001.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200823\Reduced\R_per_A_Ganymede.0001.ec.fits'] ;use this for Sunlight
      end
      'UT200823': begin
         if sunlit eq 0 then Fit_Me_To_The_Scatter    = ['D:\DATA\Apache Point\Echelle\Io Eclipses\UT200823\Reduced\R_per_A_Jupiter_Center.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200823\Reduced\R_per_A_Jupiter_Scatter.0001.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200823\Reduced\R_per_A_Jupiter_Scatter.0002.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200823\Reduced\R_per_A_Jupiter_Scatter.0003.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200823\Reduced\R_per_A_Jupiter_Scatter.0004.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200908\Reduced\R_per_A_Jupiter_Scatter.0017.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200908\Reduced\R_per_A_Jupiter_Scatter.0030.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200908\Reduced\R_per_A_Jupiter_Scatter.0032.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT201001\Reduced\R_per_A_Jupiter_Scatter.0001.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT201001\Reduced\R_per_A_Jupiter_Scatter.0002.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT201001\Reduced\R_per_A_Jupiter_Scatter.0003.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT180320\Reduced\R_per_A_Jovian_Scatter.0001.ec.fits' , $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT201001\Reduced\R_per_A_Jupiter_Scatter.0001.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT201001\Reduced\R_per_A_Jupiter_Scatter.0002.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT201001\Reduced\R_per_A_Jupiter_Scatter.0003.ec.fits' ]
        if sunlit eq 1 then Fit_Me_To_The_Scatter    = ['D:\DATA\Apache Point\Echelle\Io Eclipses\UT200823\Reduced\R_per_A_Ganymede.0001.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT190812\Reduced\R_per_A_Ganymede.0001.ec.fits'] ;use this for Sunlight
      end
      'UT200908': begin
        if sunlit eq 0 then Fit_Me_To_The_Scatter    = ['D:\DATA\Apache Point\Echelle\Io Eclipses\UT200908\Reduced\R_per_A_Jupiter_Center.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200823\Reduced\R_per_A_Jupiter_Scatter.0001.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200823\Reduced\R_per_A_Jupiter_Scatter.0002.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200823\Reduced\R_per_A_Jupiter_Scatter.0003.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200823\Reduced\R_per_A_Jupiter_Scatter.0004.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200908\Reduced\R_per_A_Jupiter_Scatter.0017.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200908\Reduced\R_per_A_Jupiter_Scatter.0030.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200908\Reduced\R_per_A_Jupiter_Scatter.0032.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT201001\Reduced\R_per_A_Jupiter_Scatter.0001.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT201001\Reduced\R_per_A_Jupiter_Scatter.0002.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT201001\Reduced\R_per_A_Jupiter_Scatter.0003.ec.fits' ]
        if sunlit eq 1 then Fit_Me_To_The_Scatter    = ['D:\DATA\Apache Point\Echelle\Io Eclipses\UT200908\Reduced\R_per_A_Ganymede.0035.ec.fits', $
           'D:\DATA\Apache Point\Echelle\Io Eclipses\UT201001\Reduced\R_per_A_Ganymede.0001.ec.fits'] ;use this for Sunlight
      end
      'UT201001': begin
        if sunlit eq 0 then Fit_Me_To_The_Scatter    = ['D:\DATA\Apache Point\Echelle\Io Eclipses\UT201001\Reduced\R_per_A_Jupiter_Center.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200823\Reduced\R_per_A_Jupiter_Scatter.0001.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200823\Reduced\R_per_A_Jupiter_Scatter.0002.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200823\Reduced\R_per_A_Jupiter_Scatter.0003.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200823\Reduced\R_per_A_Jupiter_Scatter.0004.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200908\Reduced\R_per_A_Jupiter_Scatter.0017.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200908\Reduced\R_per_A_Jupiter_Scatter.0030.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200908\Reduced\R_per_A_Jupiter_Scatter.0032.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT201001\Reduced\R_per_A_Jupiter_Scatter.0001.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT201001\Reduced\R_per_A_Jupiter_Scatter.0002.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT201001\Reduced\R_per_A_Jupiter_Scatter.0003.ec.fits' ]
        if sunlit eq 1 then Fit_Me_To_The_Scatter    = ['D:\DATA\Apache Point\Echelle\Io Eclipses\UT201001\Reduced\R_per_A_Ganymede.0001.ec.fits'] ;use this for Sunlight
      end
      'UT201017': begin
        if sunlit eq 0 then Fit_Me_To_The_Scatter    = ['D:\DATA\Apache Point\Echelle\Io Eclipses\UT201017\Reduced\R_per_A_Jupiter_Center.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200823\Reduced\R_per_A_Jupiter_Scatter.0001.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200823\Reduced\R_per_A_Jupiter_Scatter.0002.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200823\Reduced\R_per_A_Jupiter_Scatter.0003.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200823\Reduced\R_per_A_Jupiter_Scatter.0004.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200908\Reduced\R_per_A_Jupiter_Scatter.0017.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200908\Reduced\R_per_A_Jupiter_Scatter.0030.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200908\Reduced\R_per_A_Jupiter_Scatter.0032.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT201001\Reduced\R_per_A_Jupiter_Scatter.0001.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT201001\Reduced\R_per_A_Jupiter_Scatter.0002.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT201001\Reduced\R_per_A_Jupiter_Scatter.0003.ec.fits' ]
        if sunlit eq 1 then Fit_Me_To_The_Scatter    = ['D:\DATA\Apache Point\Echelle\Io Eclipses\UT201017\Reduced\R_per_A_Callisto.0015.ec.fits'] ;use this for Sunlight
      end
      'UT210601': begin
        Fit_Me_To_The_Scatter    = ['D:\DATA\Apache Point\Echelle\Io Eclipses\UT201017\Reduced\R_per_A_Jupiter_Center.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200823\Reduced\R_per_A_Jupiter_Scatter.0001.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200823\Reduced\R_per_A_Jupiter_Scatter.0002.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200823\Reduced\R_per_A_Jupiter_Scatter.0003.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200823\Reduced\R_per_A_Jupiter_Scatter.0004.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200908\Reduced\R_per_A_Jupiter_Scatter.0017.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200908\Reduced\R_per_A_Jupiter_Scatter.0030.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200908\Reduced\R_per_A_Jupiter_Scatter.0032.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT201001\Reduced\R_per_A_Jupiter_Scatter.0001.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT201001\Reduced\R_per_A_Jupiter_Scatter.0002.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT201001\Reduced\R_per_A_Jupiter_Scatter.0003.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT210601\Reduced\R_per_A_Jupiter_Scatter.0019.ec.fits', $         
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT210601\Reduced\R_per_A_Jupiter_Scatter.0020.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT210601\Reduced\R_per_A_Jupiter_Scatter.0025.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT210601\Reduced\R_per_A_Jupiter_Scatter.0026.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT210609\Reduced\R_per_A_Jupiter_Scatter.0023.ec.fits', $
          'D:\DATA\Apache Point\Echelle\Io Eclipses\UT210609\Reduced\R_per_A_Jupiter_Scatter.0024.ec.fits' ]
      end
      'UT210609': begin
          Fit_Me_To_The_Scatter    = ['D:\DATA\Apache Point\Echelle\Io Eclipses\UT201017\Reduced\R_per_A_Jupiter_Center.ec.fits', $
            'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200823\Reduced\R_per_A_Jupiter_Scatter.0001.ec.fits', $
            'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200823\Reduced\R_per_A_Jupiter_Scatter.0002.ec.fits', $
            'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200823\Reduced\R_per_A_Jupiter_Scatter.0003.ec.fits', $
            'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200823\Reduced\R_per_A_Jupiter_Scatter.0004.ec.fits', $
            'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200908\Reduced\R_per_A_Jupiter_Scatter.0017.ec.fits', $
            'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200908\Reduced\R_per_A_Jupiter_Scatter.0030.ec.fits', $
            'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200908\Reduced\R_per_A_Jupiter_Scatter.0032.ec.fits', $
            'D:\DATA\Apache Point\Echelle\Io Eclipses\UT201001\Reduced\R_per_A_Jupiter_Scatter.0001.ec.fits', $
            'D:\DATA\Apache Point\Echelle\Io Eclipses\UT201001\Reduced\R_per_A_Jupiter_Scatter.0002.ec.fits', $
            'D:\DATA\Apache Point\Echelle\Io Eclipses\UT201001\Reduced\R_per_A_Jupiter_Scatter.0003.ec.fits', $
            'D:\DATA\Apache Point\Echelle\Io Eclipses\UT210601\Reduced\R_per_A_Jupiter_Scatter.0019.ec.fits', $
            'D:\DATA\Apache Point\Echelle\Io Eclipses\UT210601\Reduced\R_per_A_Jupiter_Scatter.0020.ec.fits', $
            'D:\DATA\Apache Point\Echelle\Io Eclipses\UT210601\Reduced\R_per_A_Jupiter_Scatter.0025.ec.fits', $
            'D:\DATA\Apache Point\Echelle\Io Eclipses\UT210601\Reduced\R_per_A_Jupiter_Scatter.0026.ec.fits', $
            'D:\DATA\Apache Point\Echelle\Io Eclipses\UT210609\Reduced\R_per_A_Jupiter_Scatter.0023.ec.fits', $
            'D:\DATA\Apache Point\Echelle\Io Eclipses\UT210609\Reduced\R_per_A_Jupiter_Scatter.0024.ec.fits' ]
        end
    endcase

    ; setup the plot axis
      spec      = MRDFITS(reduced_dir+'R_per_A_'+Eclipse_files[1]+'.ec.fits', 0, header, /fscale, /unsigned, /silent)
      WL        = sxpar(header, 'CRVAL1')+findgen(N_elements(spec))*sxpar(header, 'Cdelt1')
      spec      = interpolate(spec, findgen(n_elements(spec)) + 10)
  
      for h = 0, N_elements(orders)-1 do begin
        order           = orders[h]
        xr              = order.WL_range
        Case date of
          'UT180320': begin
            if order.name eq 'order_Na_D' then begin
              yr_residual = [-4, 16]
              YR          = [21,85]
              ;YR          = [30,116]
              adjust_wavelength_solution = []
              offset_waterfall = 20.                    ; for clarity offset the brightness in each spectrum of the waterfall by this amount
            endif
            if order.name eq 'order_6300' then begin
              yr_residual   = [-2, 21]
              adjust_wavelength_solution = 1
            endif
            if order.name eq 'order_6364' then begin
              yr_residual   = [-2, 6]
              adjust_wavelength_solution = 1
            endif
          end
          'UT190812': begin
            if order.name eq 'order_Na_D' then begin
              yr_residual = [-8, 21]
              if sunlit eq 0 then YR          = [40,250]
              if sunlit eq 1 then YR          = [90,500]
              adjust_wavelength_solution = []
            endif
            if order.name eq 'order_6300' then begin
              yr_residual   = [-2, 30]
              adjust_wavelength_solution = 1
            endif
            if order.name eq 'order_6364' then begin
              yr_residual   = [-8, 15]
              adjust_wavelength_solution = 1
            endif
          end
          'UT200823': begin
            if order.name eq 'order_Na_D' then begin
              yr_residual = [-4, 9]
              YR          = [28,170]
              if sunlit eq 1 then yr_residual = [-4, 100]
              if sunlit eq 1 then YR          = [0,2200]
              adjust_wavelength_solution = []
            endif
            if order.name eq 'order_6300' then begin
              yr_residual   = [-2, 15]
              adjust_wavelength_solution = 1
            endif
            if order.name eq 'order_6364' then begin
              yr_residual   = [-3, 6]
              adjust_wavelength_solution = 1
            endif
          end
          'UT200908': begin
            if order.name eq 'order_Na_D' then begin
              yr_residual = [-4, 14]
              YR          = [45,280]
              if sunlit eq 1 then yr_residual = [-50, 350] ;for in sunlight
              if sunlit eq 1 then YR          = [300, 3400] ;for in sunlight
              adjust_wavelength_solution = []
            endif
            if order.name eq 'order_6300' then begin
              yr_residual   = [-2, 25]
              adjust_wavelength_solution = 1
            endif
            if order.name eq 'order_6364' then begin
              yr_residual   = [-2, 8]
              adjust_wavelength_solution = 1
            endif
          end
          'UT201001': begin
            if order.name eq 'order_Na_D' then begin
              yr_residual = [-4, 13]
              YR          = [30, 180]
              if sunlit eq 1 then yr_residual = [-50, 450] ;for in sunlight
              if sunlit eq 1 then YR          = [200, 3500] ;for in sunlight
              adjust_wavelength_solution = 1
            endif
            if order.name eq 'order_6300' then begin
              yr_residual   = [-2, 28]
              adjust_wavelength_solution = 1
            endif
            if order.name eq 'order_6364' then begin
              yr_residual   = [-2, 9]
              adjust_wavelength_solution = 1
            endif
          end
          'UT201017': begin
            if order.name eq 'order_Na_D' then begin
              yr_residual = [-4, 10]
              YR          = [45, 120]
              if sunlit eq 1 then yr_residual = [-50, 200] ;for in sunlight
              if sunlit eq 1 then YR          = [200, 3500] ;for in sunlight
              adjust_wavelength_solution = []
            endif
            if order.name eq 'order_6300' then begin
              yr_residual   = [-2, 35]
              adjust_wavelength_solution = 1
            endif
            if order.name eq 'order_6364' then begin
              yr_residual   = [-3, 12]
              adjust_wavelength_solution = 1
            endif
          end
          'UT210601': begin
            adjust_wavelength_solution = []
            if order.name eq 'order_Na_D' then begin
              yr_residual = [-4, 10]
              YR          = [45, 120]
              adjust_wavelength_solution = []
            endif
            if order.name eq 'order_6300' then begin
              yr_residual   = [-2, 35]
              adjust_wavelength_solution = 1
            endif
            if order.name eq 'order_6364' then begin
              yr_residual   = [-3, 12]
              adjust_wavelength_solution = 1
            endif
          end
          'UT210609': begin
            adjust_wavelength_solution = []
            if order.name eq 'order_Na_D' then begin
              yr_residual = [-4, 10]
              YR          = [45, 120]
              adjust_wavelength_solution = []
            endif
            if order.name eq 'order_6300' then begin
              yr_residual   = [-2, 35]
              adjust_wavelength_solution = 1
            endif
            if order.name eq 'order_6364' then begin
              yr_residual   = [-3, 12]
              adjust_wavelength_solution = 1
            endif
          end
          else: 
        endcase

      ; Define arrays
        include_WLs            = where( (wl gt xr[0]) and (wl lt xr[1]), /NULL )                      ; wavelengths to use in the scattered light fitting
        include_WLs_broad      = where( (wl gt (xr[0]-20.)) and (wl lt (xr[1]+20.)), /NULL )          ; sometime we'll need a broader range of wavelengths to use in the scattered light fitting
        if (ingress and (order.name eq 'order_Na_D')) then include_WLs = where( (abs(wl - Na1+.4) lt 1.0) or (abs(wl - Na2+.4) lt 1.0), /NULL ); wavelength regions over which we'll fit the jovian scatter
        if (ingress and (order.name eq 'order_NaIR')) then include_WLs = where( (abs(wl - Na3+.4) lt .6) or (abs(wl - Na4+.4) lt .6), /NULL ); wavelength regions over which we'll fit the jovian scatter
        if ((not ingress) and (order.name eq 'order_Na_D')) then include_WLs = where( (((abs(wl - Na1-.4) lt 0.51) and abs(wl - Na1-.4) gt 0.2)) or ((abs(wl - Na2-.4) lt 0.51) and abs(wl - Na2-.4) gt 0.2), /NULL ); wavelength regions over which we'll fit the jovian scatter
        if ((not ingress) and (order.name eq 'order_NaIR')) then include_WLs = where( (abs(wl - Na3-.4) lt .6) or (abs(wl - Na4-.4) lt .6), /NULL ); wavelength regions over which we'll fit the jovian scatter

      if order.name eq 'order_Na_D' then MX_plus_B_parinfo   = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0]}, 2)
      plot_WLs               = where( (wl gt xr[0]) and (wl lt xr[1]), /NULL )                      ; wavelengths to actually plot
      residual_array         = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
      plot_spec_array        = fltarr(N_elements(plot_WLs), n_elements(Eclipse_files))
      plot_residual_array    = fltarr(N_elements(plot_WLs), n_elements(Eclipse_files))
      err_plot_residual_array= fltarr(N_elements(plot_WLs), n_elements(Eclipse_files))
      LSF_Fit_array          = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
      WL_array               = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
      T_P_Shadow_array       = fltarr(n_elements(Eclipse_files))
      T_U_Shadow_array       = fltarr(n_elements(Eclipse_files))
      DopplerShift_array     = fltarr(N_elements(Io_airglow), n_elements(Eclipse_files))
      Gof                    = 0.                                                                   ; Goodness of fit (co-addded mean residual)
      GOF_array              = fltarr(N_elements(Fit_Me_To_The_Scatter))
      Use_Scatter            = strarr(N_elements(Eclipse_files))
      Possible_WL_offset     = fltarr(N_elements(Eclipse_files), N_elements(Fit_Me_To_The_Scatter)) ; Save any offsets from Io's expected rest wavelength, then correct any error in the wavelength solution.
      Possible_line_width    = fltarr(N_elements(Eclipse_files), N_elements(Fit_Me_To_The_Scatter)) ; Save the linewidths from all type of scatter subtraction, then average and peg this parameter.
      Line_indices           = where( (Io_airglow gt min(WL[plot_WLs])) and (Io_airglow lt max(WL[plot_WLs])), /Null) ;which lines are in this wavelength range
      Line_index             = Line_indices[0]

      for i = 0, n_elements(Eclipse_files)-1 do begin
        spec        = MRDFITS(reduced_dir+'R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, header, /fscale, /unsigned, /silent)
        spec_err    = MRDFITS(reduced_dir+'sig_R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, sig_header, /fscale, /unsigned, /silent)
        
        WL          = sxpar(header, 'CRVAL1')+findgen(N_elements(spec))*sxpar(header, 'Cdelt1')
        ;window, i
        ;full = mrdfits(reduced_dir+'fullspec' + Eclipse_files[i] + '.ec.fits', 0, /fscale, /silent, /unsigned )
        ;cgplot, wl[include_WLs], full[include_WLs]

        for K = 0, n_elements(Fit_Me_To_The_Scatter)-1 do begin
          ; Which Jovian scatter are we testing out?
            Jup_Cal_Tell_Corr = MRDFITS(Fit_Me_To_The_Scatter[K], 0, Jupiter_Scatter_header, /fscale, /unsigned, /silent)
            Jup_Cal_Tell_Corr      = interpolate(Jup_Cal_Tell_Corr, findgen(n_elements(Jup_Cal_Tell_Corr)) + 15.)
          ; Define fitting indicies we need to avoid because of airglow itself.
            Ios_Airglow            = [Io_airglow + Io_airglow * sxpar(header, 'IO_DOPPL') / cspice_clight()] ; wavelength locations of all airglow lines, Ignore telluric K
            ARCES_1sigma_bandwidth = (mean(xr)/31500.) / 2.3548                                              ; FWHM in Angstroms / 2sqrt(2ln2) = 68% of flux enclosed within 1 sigma --> Appropriate for weak telluric airglow
            Io_AG_Ind = [] & Telluric_AG_Ind = []
            for j = 0, N_elements(Ios_Airglow)-1 do Io_AG_ind = [Io_AG_Ind, where(abs(Ios_Airglow[j] - WL) lt 2.*ARCES_1sigma_bandwidth, /NULL)]                    ; indicies within 2 sigma of an Io airglow line
            for j = 0, N_elements(Telluric_Airglow)-1 do Telluric_AG_ind = [Telluric_AG_Ind, where(abs(Telluric_Airglow[j] - WL) lt ARCES_1sigma_bandwidth, /NULL)] ; indicies within 1 sigma of any Telluric airglow
            AG_ind = [Io_AG_Ind, Telluric_AG_Ind]
            AG_ind = AG_ind[UNIQ(AG_ind, SORT(AG_ind))]

          ; Fit the Jovian scattered light, use Amoeba/Correlation to match Jupiter in terms of smoothing and shifting, then multiply and add until matches the scattered background
            correl_indicies  = cgSetDifference(include_Wls, AG_Ind)
            Mult_and_add     = mpfitfun('MX_Plus_B', Jup_Cal_Tell_Corr[correl_indicies], Spec[correl_indicies], Spec_err[correl_indicies], parinfo = MX_plus_B_parinfo, $
              [mean(Spec[correl_indicies]/Jup_Cal_Tell_Corr[correl_indicies]), 0.], /NaN, status = Scatter_fit_status, PERROR = err_a, /quiet)
            jup              = Mult_and_add[0]*Jup_Cal_Tell_Corr + Mult_and_add[1]

          ; Iterate Amoeba until it gives an answer
            ;correl_indicies  = cgSetDifference(include_Wls, AG_Ind)
            trial_smooth_and_shift = AMOEBAX(1.e-4, 1.e-4, function_name='shift_smooth', SCALE = [1., 1.], P0 = [4., -5.], FUNCTION_VALUE = fval, NMAX = 500, NCalls = NCalls)
            Case 1 of
              ( (fval[1] gt 1.) and (fval[1] lt 1.2) and (N_elements(trial_smooth_and_shift) eq 2) ): smooth_and_shift = trial_smooth_and_shift
              else: begin
                correl_indicies  = cgSetDifference(include_Wls_broad, AG_Ind)                    ; expand the wavelength range that we're fitting smooths & shifts over. 'correl_indicies' is passed in the common block
                trial_smooth_and_shift2 = AMOEBAX(1.e-4, 1.e-4, function_name='shift_smooth', SCALE = [2.5, 10.], P0 = [1.25, 0.], FUNCTION_VALUE = fval, NMAX = 500, NCalls = NCalls)
                if ( (fval[1] gt 1.) and (fval[1] lt 1.2) and (N_elements(trial_smooth_and_shift2) eq 2) ) then smooth_and_shift = trial_smooth_and_shift2 $
                else smooth_and_shift = [0, shift_only]
              end
            endcase
            correl_indicies  = cgSetDifference(include_Wls, AG_Ind)
            ;cgplot, Wl[correl_indicies], Spec[correl_indicies], xr = xr, /ynozero
            ;cgplot, Wl[correl_indicies], Jup[correl_indicies], color = 'red', /overplot
            ;smoothed_shifted = interpolate(gauss_smooth(jup, smooth_and_shift[0]), findgen(n_elements(Jup)) - smooth_and_shift[1])
            ;cgplot, Wl[correl_indicies], smoothed_shifted[correl_indicies], color = 'blue', /overplot
            ;print, smooth_and_shift
            ;stop

          smoothed_shifted = interpolate(gauss_smooth(jup, smooth_and_shift[0]), findgen(n_elements(Jup)) - smooth_and_shift[1])
          Mult_and_add     = mpfitfun('MX_Plus_B', smoothed_shifted[include_WLs_broad ], Spec[include_WLs_broad ], Spec_err[include_WLs_broad ], parinfo = MX_plus_B_parinfo, $
            [mean(Spec[include_WLs_broad ]/smoothed_shifted[include_WLs_broad ]), 0.], /NaN, status = Scatter_fit_status, PERROR = err_a, /quiet)
          scatter_fit      = Mult_and_add[0]*smoothed_shifted + Mult_and_add[1]
          residual         = Spec - scatter_fit
          GOF_array[K]     = GOF + mean(abs(residual[correl_indicies]))

          ; Get the linewidth and the determine any systematic wavelength shift from Io's rest
          ; Ultimately we will peg both of these parameters.
            plot_spec_array[*, i]     = spec[plot_WLs] /1.e3
            LSF_fitting_ind           = cgSetDifference(where( abs(wl - Ios_Airglow[line_Index]) lt 0.3, /NULL), Telluric_AG_Ind)  ; fit region within +/- 0.3A of line center, excluding Telluric airglow
            initial_guess             = [2.5e3, Ios_Airglow[line_Index], 0.077]                                      ; rough values are fine here [height, wavelength, linewidth_sigma]
            fa                        = {x:double(WL[LSF_fitting_ind]), y:double(residual[LSF_fitting_ind]), err:double( sqrt(abs( residual[LSF_fitting_ind] )) )}
            a                         = mpfit('Gaussian_for_MPFIT', initial_guess, PERROR = err_a, funct=fa, maxiter=50, STATUS = Did_it_work, /Quiet, NPEGGED = NPEGGED)
            Possible_WL_offset[i,K]   = a[1] - Ios_Airglow[line_Index]
            Possible_Line_Width[i,K]  = a[2]
        endfor

        ; write which potial background spectrum best matches the science spectrum
          best_GoF       = min(GOF_array, best_fit)
          Use_Scatter[i] = Fit_Me_To_The_Scatter[best_fit]
          print, 'Frame: ', Eclipse_files[i], ' ', order.name,' is best paired with the Jupiter scatter from: ', Use_Scatter[i], ', Goodness of fit =', best_GoF

        ;            Case date of
        ;              'UT180320': begin
        ;                if order.name eq 'order_Na_D' then junk = 0. ;continue ;'R_per_A_Jovian_Scatter.0003.ec.fits' ; Manually input desired scatting spectrum.
        ;                if order.name eq 'order_6300' then junk = 0.; Use_Scatter[i] = 'R_per_A_Jovian_Scatter.0003.ec.fits' ; Manually input desired scatting spectrum.
        ;                if order.name eq 'order_6364' then Use_Scatter[i] = 'R_per_A_Jovian_Scatter.0003.ec.fits' ; Manually input desired scatting spectrum.
        ;                order.waterfall_plot_title = 'Io''s Airglow Response Following ' + date + ' Ingress'
        ;              end
        ;              'UT190812': begin
        ;                best_scatter_spectrum = 'R_per_A_Jupiter_Center.ec.fits' ; Manually input desired scatting spectrum.
        ;              end
        ;              else: junk = 0
        ;            endcase
        print, 'Using: ', Use_Scatter[i], ' with a smooth & shift of' 
      endfor

      ; -------------------------------------------Waterfall Plot Postscript Io and Fit Jovian scatter--------------------------------------------------------------------------------------------------------------
      ; Okay, now that we've established which Jupiter scatter spectrum best optimzes the residual, we can actually use it for scattered light subtraction...
      ; Get color versus ingress time & setup plot positions
        timeColors = BytScl(findgen(N_elements(eclipse_files)), Min=0, Max=N_elements(eclipse_files), Top=N_elements(eclipse_files))
        cgLoadCT, 33, NColors=N_elements(eclipse_files), /reverse
        window, 0
        if order.npanels eq 3 then begin
          pos        = cgLayout([1,3], OXMargin=[10,3], OYMargin=[8,5], YGap=0)
          Pos[1,0]   = Pos[1,0]*.8 & Pos[3,1] = Pos[3,1]*.8
          Pos[1,1]   = Pos[1,1]*.84 & Pos[3,2] = Pos[3,2]*.84
        endif else begin
          pos        = cgLayout([1,2], OXMargin=[10,3], OYMargin=[8,5], YGap=0)
          Pos[1,0]   = Pos[1,0]*.7 & Pos[3,1] = Pos[3,1]*.7
        endelse
        WDELETE, 0

      cgPS_Open, filename = Reduced_Dir + order.name + '_scatterfit.eps', /ENCAPSULATED, xsize = 7.5, ysize = 5
      !P.font = 1
      device, SET_FONT = 'Helvetica Bold', /TT_FONT
      !p.charsize = 1.5

      if not keyword_set(yr) then yr = minmax(plot_spec_array[where(plot_spec_array gt 0., /null)])

      cgplot, WL, spec/1.e3, psym = 0, Ytitle = 'kR / '+cgsymbol('Angstrom'), $
        yr = yr, /xstyle, /ystyle, title = order.waterfall_plot_title, $
        xr = xr, /nodata, pos = pos[*,0], xtickformat = '(A1)', xminor = 10, xticklen = .025

      for i = 0, n_elements(Eclipse_files)-1 do begin
        spec         = MRDFITS(reduced_dir+'R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, header, /fscale, /unsigned, /silent)
        spec_err     = MRDFITS(reduced_dir+'sig_R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, sig_header, /fscale, /unsigned, /silent)
        spec         = spec / 1.e3       ; to KR
        spec_err     = spec_err / 1.e3   ; to KR
        if keyword_set(adjust_wavelength_solution) then $
          WL         = sxpar(header, 'CRVAL1') + findgen(N_elements(spec))*sxpar(header, 'Cdelt1') - median(Possible_WL_offset[1:*,*]) $
        else $
          WL         = sxpar(header, 'CRVAL1') + findgen(N_elements(spec))*sxpar(header, 'Cdelt1')

        ; Which Jovian scatter are we using?
          Jup_Cal_Tell_Corr = MRDFITS(Use_Scatter[i], 0, Jup_Scatter_header, /fscale, /unsigned, /silent)
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;This is hacky. Better solution will be made soon. The wavelength correction depends on the scatter used for the moment guess is required.
          Jup_Cal_Tell_Corr      = interpolate(Jup_Cal_Tell_Corr, findgen(n_elements(Jup_Cal_Tell_Corr)) + 15.)         ;CORRECTION for wavelength shift at 
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

        ; Define fitting indicies we need to avoid because of airglow itself.
          Ios_Airglow    = [Io_airglow + Io_airglow * sxpar(header, 'IO_DOPPL') / cspice_clight()]                                      ; wavelength locations of all airglow lines, Ignore telluric K
          ARCES_1sigma_bandwidth =(mean(xr)/31500.) / 2.3548 ; FWHM in Angstroms / 2sqrt(2ln2) = 68% of flux enclosed within 1 sigma --> Appropriate for weak telluric airglow
          Io_AG_Ind = [] & Telluric_AG_Ind = []
          for j = 0, N_elements(Ios_Airglow)-1 do Io_AG_ind = [Io_AG_Ind, where(abs(Ios_Airglow[j] - WL) lt 2.*ARCES_1sigma_bandwidth, /NULL)] ;  indicies within 2 sigma of an Io airglow line
          for j = 0, N_elements(Telluric_Airglow)-1 do Telluric_AG_ind = [Telluric_AG_Ind, where(abs(Telluric_Airglow[j] - WL) lt ARCES_1sigma_bandwidth, /NULL)] ; indicies within 1 sigma of any Telluric airglow
          AG_ind = [Io_AG_Ind, Telluric_AG_Ind]
          AG_ind = AG_ind[UNIQ(AG_ind, SORT(AG_ind))]

        ; Fit the Jovian scattered light, use Amoeba/Correlation to match Jupiter in terms of smoothing and shifting, then multiply and add until matches the scattered background
          correl_indicies  = cgSetDifference(include_Wls, AG_Ind)
          weights          = 0.
          for ind = 0, N_elements(line_Indices)-1 do weights = weights + 1./abs(Ios_Airglow[line_Indices[ind]] - correl_indicies)^(1.5)
          Mult_and_add     = mpfitfun('MX_Plus_B', Jup_Cal_Tell_Corr[correl_indicies], Spec[correl_indicies], Spec_err[correl_indicies], parinfo = MX_plus_B_parinfo, $
            [mean(Spec[correl_indicies]/Jup_Cal_Tell_Corr[correl_indicies]), 0.], /NaN, status = Scatter_fit_status, PERROR = err_a, /quiet)
          jup              = Mult_and_add[0]*Jup_Cal_Tell_Corr + Mult_and_add[1]

        ; Iterate Amoeba until it gives an answer
          trial_smooth_and_shift = AMOEBAX(1.e-4, 1.e-4, function_name='shift_smooth', SCALE = [1., 1.], P0 = [4., -5.], FUNCTION_VALUE = fval, NMAX = 1000, NCalls = NCalls)
          Case 1 of
            ( (fval[1] gt 1.) and (fval[1] lt 1.2) and (N_elements(trial_smooth_and_shift) eq 2) ): smooth_and_shift = trial_smooth_and_shift
            else: begin
              correl_indicies  = cgSetDifference(include_Wls_broad, AG_Ind)                    ; expand the wavelength range that we're fitting smooths & shifts over. 'correl_indicies' is passed in the common block
              trial_smooth_and_shift2 = AMOEBAX(1.e-4, 1.e-4, function_name='shift_smooth', SCALE = [2.5, 10.], P0 = [1.10, 0.], FUNCTION_VALUE = fval, NMAX = 1000, NCalls = NCalls)
              if ( (fval[1] gt 1.) and (fval[1] lt 1.2) and (N_elements(trial_smooth_and_shift2) eq 2) ) then smooth_and_shift = trial_smooth_and_shift2 $
              else smooth_and_shift = [0, shift_only]
            end
          endcase
          correl_indicies  = cgSetDifference(include_Wls, AG_Ind)
  
          smoothed_shifted = interpolate(gauss_smooth(jup, smooth_and_shift[0]), findgen(n_elements(Jup)) - smooth_and_shift[1])
          Mult_and_add     = mpfitfun('MX_plus_B', smoothed_shifted[correl_indicies], Spec[correl_indicies], Spec_err[correl_indicies], [mean(Spec[correl_indicies]/smoothed_shifted[correl_indicies]), 0.], $
            /NaN, status = Scatter_fit_status, PERROR = err_a, /quiet, weights = weights, parinfo = MX_plus_B_parinfo)
          scatter_fit      = Mult_and_add[0]*smoothed_shifted + Mult_and_add[1]
  
;          if ((order.name eq 'order_Na_D') and (date EQ 'UT200908')) then begin
;            scatter_fit      = (Mult_and_add[0])*(smoothed_shifted*1.01*(-1.0)) + Mult_and_add[1] + 5000
;            residual         = (Spec) - scatter_fit
;            residual_err     = make_array( N_elements(residual), value = stddev(residual[correl_indicies]) )
;          endif 
          residual         = Spec - scatter_fit
          residual_err     = make_array( N_elements(residual), value = stddev(residual[correl_indicies]) )
  
          plot_residual_array[*, i]     = residual[plot_WLs]
          err_plot_residual_array[*, i] = residual_err[plot_WLs]
          DopplerShift_array[*, i]      = Ios_Airglow

        ; plot the spectrum, fit and indices used for the fit (if we're in debug mode)
          if ((order.name eq 'order_6300') and (date EQ 'UT180320')) then begin
            ;cgplot, WL, scatter_fit + ((n_elements(Eclipse_files)-1)*12.) - (i*12.), color = timeColors[i], linestyle = 1, thick = 5, /overplot
            ;cgplot, WL, spec + ((n_elements(Eclipse_files)-1)*12.) - (i*12.), color = timeColors[i], thick = 5, /overplot
            cgplot, WL, scatter_fit, color = timeColors[i], linestyle = 1, thick = 5, /overplot
            cgplot, WL, spec, color = timeColors[i], thick = 5, /overplot
          endif else begin
            ;cgplot, WL[correl_indicies], scatter_fit[correl_indicies], color = timeColors[i], psym=14, /overplot ----> useful to confirm spectral regions being fit for jovian scatter
            cgplot, WL, scatter_fit, color = timeColors[i], linestyle = 1, thick = 5, /overplot
            cgplot, WL, spec, color = timeColors[i], thick = 5, /overplot
          endelse
        T_P_Shadow_array[i]      = sxpar(header, 'T_P_Shad')
        T_U_Shadow_array[i]      = sxpar(header, 'T_USHADO')
      endfor

      ; Annotate w/ legend & text
        IF ((date EQ 'UT180320') and (order.name eq 'order_6300')) THEN BEGIN
          AL_legend, ['Raw Io Spectrum','Jupiter Scattered Light Fit'], Psym = [0,0], linestyle = [0,1], charsize = 1., linsize = 0.5, position = [.67, .9], /normal
          cgtext, 6297.8, yr[1]*0.96, strcompress(string(T_U_Shadow_array[-1], format = '(F10.1)'), /remove_all) +'min Post Umbral Ingress', charsize = 1.4, alignment = 0, color = timeColors[5]
          cgtext, 6297.8, yr[0]*1.07, strcompress(string(T_U_Shadow_array[0], format = '(F10.1)'), /remove_all) +'min Post Umbral Ingress', charsize = 1.4, alignment = 0, color = timeColors[0]
        ENDIF

      ; Residual plot
        lines_to_fit = where(Ios_airglow gt xr[0] and Ios_airglow lt xr[1], /Null, count_lines)
        if (order.name eq 'order_6300') or (order.name eq 'order_6364') or (order.name eq 'order_Na_D') then count_lines = count_lines else count_lines = 0 ;only plot fits for Na and O 6300/6364A
        if count_lines gt 0 then LSF_Fit_array = fltarr(N_elements(Plot_WLs), n_elements(Eclipse_files))
        if not keyword_set(yr_residual) then yr_residual = [-2.*stddev(plot_residual_array), 2.*stddev(plot_residual_array)]

      if order.npanels eq 3 then begin
        cgtext, 0.03, 0.35, 'Residual [kR / '+cgsymbol('Angstrom')+']', orientation = 90, alignment = 0.5, /normal
        cgplot, spec, WL, psym = 0, xtickformat = '(A1)', $
          yr = yr_residual, xr = xr, /nodata, pos = pos[*,1], /xstyle, /ystyle, /noerase, xminor = 10, xticklen = .05
      endif else begin
        cgtext, 0.036, 0.25, 'Residual [kR / '+cgsymbol('Angstrom')+']', orientation = 90, alignment = 0.5, /normal
        cgplot, spec, WL, psym = 0, yr = yr_residual, xr = xr, /nodata, pos = pos[*,1], /xstyle, /ystyle, /noerase, xminor = 10, xticklen = .05, xtitle = 'Wavelength ['+cgsymbol('Angstrom')+']'
      endelse

      for i = 0, n_elements(Eclipse_files)-1 do begin
        for k = 0, N_elements(Io_Airglow)-1 do cgplot, replicate(DopplerShift_array[k,i], 2), [.9*yr_residual[1], yr_residual[1]], COLOR = timeColors[i], /overplot
        cgplot, WL[plot_WLs], smooth(plot_residual_array[*, i], 1), /OVERPLOT, COLOR = timeColors[i], thick = 3, psym=10
        header         = headfits(reduced_dir+'R_per_A_' + Eclipse_files[i] + '.ec.fits')
        fullspec_header= headfits(reduced_dir+'fullspec' + Eclipse_files[i] + '.ec.fits')                     ; not sure why but the EXPTIME header gets corrupted?
        Ios_Airglow    = [Io_airglow + Io_airglow * sxpar(header, 'IO_DOPPL') / cspice_clight()]              ; wavelength locations of all airglow lines

        ; fit any airglow lines and save the brightness
        if count_lines gt 0 then begin
          for j = 0, count_lines-1 do begin

            ; Now do the line fitting to the residual. First define Line Spread Function (LSF) Gaussian fit parameter info
              parinfo               = replicate({value:0., fixed:0, limited:[0,0], limits:[0.,0]}, 3)
              parinfo[1].fixed      = 0
              parinfo[1].value      = Ios_airglow[lines_to_fit[j]]                                            ; Pin the line's wavelength?
              parinfo[2].fixed      = 0
              parinfo[2].value      = median(Possible_line_width[i,*])                                        ; Pin the line's wavelength?
              
            ; Gaussian fit
              LSF_fitting_ind       = where( abs(wl[plot_WLs] - Ios_airglow[lines_to_fit[j]]) lt 0.2, /NULL)
              tsum_integral_ind     = where( (wl[plot_WLs] - Ios_airglow[lines_to_fit[j]] lt 0.4) and (wl[plot_WLs] - Ios_airglow[lines_to_fit[j]] gt -0.3), /NULL) ;HACKED FIX FOR FAULTY GAUSSIAN ONLY FOR DEMONSTRATION OF 201001
              fa                    = { x:wl[plot_WLs[LSF_fitting_ind]], y:plot_residual_array[LSF_fitting_ind, i], err: err_plot_residual_array[LSF_fitting_ind, i] }
              a                     = mpfit('Gaussian_for_MPFIT', [2., parinfo[1].value, parinfo[2].value], funct=fa, maxiter=50, STATUS = Did_it_work, /Quiet, NPEGGED = NPEGGED, parinfo=parinfo)
 
;            ; Check for junk fitting. Force linewidths where necessary.
;              bad_index             = where((a lt 0.) or ~(finite(A)), /null, count_bad_fits)
;              if count_bad_fits gt 0 then begin
;                print, 'Warning: Bad line fit detected:', a[bad_index]
;                if bad_index eq 2 then begin
;                  parinfo[2].fixed = 1
;                  print, 'Forcing line width for', parinfo[1].value 
;                endif
;                if bad_index eq 1 then parinfo[1].fixed = 1
;                a = mpfit('Gaussian_for_MPFIT', [2., parinfo[1].value, parinfo[2].value], funct=fa, maxiter=50, STATUS = Did_it_work, /Quiet, NPEGGED = NPEGGED, parinfo=parinfo)
;              endif  
            ;a = abs(a)

            ; log the results for plotting up later
              LSF_Fit_array[*, i]                                  = LSF_Fit_array[*, i] + gaussian(WL[plot_WLs], a)  ; LSF_Fit
              Io_Airglow_params.brightness[lines_to_fit[j], i]     = A[0]*abs(A[2])*SQRT(2*!DPI)
              ;Io_Airglow_params.brightness[lines_to_fit[j], i]     = tsum(wl[plot_WLs[tsum_integral_ind]], plot_residual_array[tsum_integral_ind, i])
              Io_Airglow_params.err_brightness[lines_to_fit[j], i] = stddev(residual[correl_indicies])*A[2]*SQRT(2*!DPI)
              Io_Airglow_params.linewidth[lines_to_fit[j], i]      = A[2]
              Io_Airglow_params.linecenter[lines_to_fit[j], i]     = A[1]
              Io_Airglow_params.torus_latitude[i]                  = float(sxpar(header, 'Torus_la'))
              Io_Airglow_params.exptime[i]                         = float(sxpar(fullspec_header, 'EXPTIME')) / 60.
              Io_Airglow_params.T_P_Shadow[i]                      = sxpar(header, 'T_P_Shad')
          endfor
          cgplot, WL[plot_WLs], LSF_Fit_array[*, i], /OVERPLOT, COLOR = timeColors[i], thick = 3 ; plot the fit to Io's emission
        endif
      endfor ; loop over eclipse frames
      ;cgplot, [(wl[plot_WLs[tsum_integral_ind]])[0], (wl[plot_WLs[tsum_integral_ind]])[0]], [0, 100000], /overplot, color = 'black'
      ;cgplot, [(wl[plot_WLs[tsum_integral_ind]])[-1], (wl[plot_WLs[tsum_integral_ind]])[-1]], [0, 100000], /overplot, color = 'black'
      IF (date EQ 'UT180320') and (order.name eq 'order_6300') THEN BEGIN
        cgtext, Ios_airglow[line_index] - 0.7, 16, "Io's Doppler Shift", charsize = 1.4, alignment = 0.5
        cgtext, Io_airglow[line_index] + .5, 16, "Telluric [O I]", charsize = 1.4, alignment = 0.5
      endif
      for k = 0, N_elements(Io_Airglow)-1 do cgplot, replicate(Io_Airglow[k], 2), YR_residual, linestyle = 1, /overplot ; plot the telluric line locations

      if order.npanels eq 3 then begin

        ; Co-align to Io's Doppler Shift
          shift_by_array = transpose( DopplerShift_array) - $
            rebin(transpose(mean(DopplerShift_array, dim = 2)), n_elements(eclipse_files), N_elements(Io_Airglow))
          ind = where(io_airglow gt min(WL[plot_WLs]) and Io_airglow lt max(WL[plot_WLs]), /NULL, count)
          if count eq 0 then begin
            shifted_WL = fltarr( n_elements(eclipse_files ) )
            for i = 0, n_elements(eclipse_files )-1 do begin
              header         = headfits(reduced_dir+'R_per_A_' + Eclipse_files[i] + '.ec.fits')
              shifted_WL[i]  = [mean(WL[plot_WLs]) + mean(WL[plot_WLs])* sxpar(header, 'IO_DOPPL') / cspice_clight()]
            endfor
            shift_by = shifted_WL - mean(shifted_WL)
          endif
          if count eq 1 then shift_by = shift_by_array[*,ind]
          if count gt 1 then shift_by = mean(shift_by_array[*,ind], dim = 2)
  
          aligned_residuals = plot_residual_array
          for i = 0, n_elements(eclipse_files )-1 do begin
            aligned_residuals[*,i] = interpol( plot_residual_array[*,i], WL[plot_WLs], WL[plot_WLs] - Shift_By[i])
          endfor
          co_add = median(aligned_residuals, dimension = 2, /even)

        ; plot the aligned & co-added spectra
          yr_co_add = yr_residual ;/ [5., 2.8]
          cgplot, WL[plot_WLs], co_add, pos = pos[*,2], thick=5, /noerase, xr = xr, yr = yr_co_add, /yst, psym=10, $
            xtitle = 'Wavelength ('+cgsymbol('Angstrom')+')', yminor = 2, /nodata, xminor = 10, xticklen = .05
          cgplot, xr, [0.,0.], linestyle = 2, /overplot
  
          airglow_x = congrid( mean(DopplerShift_array, dim = 2), N_elements(Io_Airglow)*2+1)
          airglow_x = airglow_x[0:-2]
          airglow_y = Reform(Rebin([-15, 15], 2, N_elements(Io_Airglow)), 2*N_elements(Io_Airglow))
          for i = 0, N_elements(Io_Airglow)-1 do cgplot, airglow_x[2*i:2*i+1], airglow_y[2*i:2*i+1], /overplot, linestyle = 1, color = 'blue', thick = 1.
          cgplot, WL[plot_WLs], smooth(co_add, 1, /edge_mirror), thick = 5, /overplot
      endif

      ; cleanup temporary definitions written each loop in the order for loop
        junk = temporary(yr) 
        junk = temporary(yr_residual)

      cgps_close
    endfor ; order loop
    save, Io_Airglow_params, filename = Reduced_Dir+'Io_Airglow_params.sav'
    stop
  endif



;  ;part 1 and 12 to do everything APO!
;
;
;  ;==========================================  Combined lightcurve  ======================================================================
;
;  ;print, 'Ajello et al. 2008 ratio [7774:8446:9225]', [3.90e-19, 2.15e-19, 6.63e-19] / 2.15e-19
;  if part eq 2 then begin
;    restore, Reduced_Dir+'Na1Na2_fit_params.sav'
;    restore, Reduced_Dir+'O1_fit_params.sav'
;    restore, Reduced_Dir+'O2_fit_params.sav'
;    restore, Reduced_Dir+'O3_fit_params.sav'
;
;    pos =  [.12,.17,.95,.9]
;    cgPS_Open, filename = Reduced_Dir+'Combined_lightcurve.eps', /ENCAPSULATED, xsize = 6., ysize = 4.
;    !P.font=1
;    device, SET_FONT = 'Helvetica Bold', /TT_FONT
;
;    yr = [0, 5.]
;    cgplot, O2_fit_params.T_P_Shadow, O2_fit_params.Brightness, psym = 15, Ytitle = 'Disk-Averaged Brightness (KR)', xtitle = 'Minutes After Ingress', $
;      title = 'Io''s Response to Ingress', yr = yr, xr = [range1,range2], /nodata, pos = pos,YTICKFORMAT="(A1)",yticks = 1
;
;    cgLoadCT, 33, NColors=8, /reverse
;    if not ingress then x = findgen(Umbra_ET - PenUmbra_ET) / 60.
;    if ingress then x = findgen(Umbra_ET - PenUmbra_ET) / 60.
;    colors = reverse(BytScl(x, MIN=min(x), MAX=2.*max(x)))+128
;    FOR j=0,n_elements(x)-2 DO BEGIN
;      xpoly = [x[j],         x[j], x[j+1],       x[j+1],         x[j]]
;      ypoly = [  0., !Y.CRange[1], !Y.CRange[1], 0., 0.]
;      cgColorFill, xpoly, ypoly, Color=colors[j]
;    ENDFOR
;    if ingress then begin
;      xpoly = [max(x),     max(x), !X.CRange[1],  !X.CRange[1],  max(x)]
;      ypoly = [  0., !Y.CRange[1], !Y.CRange[1], 0., 0.]
;      cgColorFill, xpoly, ypoly, color = 'charcoal', /LINE_FILL, orientation = 45., thick = .1
;      cgColorFill, xpoly, ypoly, color = 'charcoal', /LINE_FILL, orientation = -45., thick = .1
;    endif
;
;    cgaxis, yaxis = 0, yr = yr, /ystyle                        ; repair the axis damage that the Penumbra did
;    cgaxis, xaxis = 0, xr = time_range                         ; repair axis damage
;    cgaxis, xaxis = 1, xr = time_range, xtickformat = '(A1)'   ; repair axis damage
;
;
;    Na1Na2_fit_params[0].Brightness = !Values.F_Nan
;    O1_fit_params[0].Brightness = !Values.F_Nan
;    O2_fit_params[0].Brightness = !Values.F_Nan
;    O3_fit_params[0].Brightness = !Values.F_Nan
;
;    READCOL,'D:\DATA\March20_sun_ec.dat',col1,t_post_eclipse,SO2_val_1,SO2_err_1,SO2_val_2,SO2_err_2,SO2_val_3,SO2_err_3,SO2_val_4,SO2_err_4,SO_val,SO_err, STRINGSKIP = '#', /Silent
;    SO2_val_1in = SO2_val_1[0]
;    SO2_val_2in = SO2_val_2[0]
;    SO2_val_3in = SO2_val_3[0]
;    SO2_val_4in = SO2_val_4[0]
;    SO_valin = SO_val[0]
;
;    SO2_val_1[0] = !Values.F_Nan
;    SO2_val_1[1] = !Values.F_Nan
;    SO2_val_2[0] = !Values.F_Nan
;    SO2_val_2[1] = !Values.F_Nan
;    SO2_val_3[0] = !Values.F_Nan
;    SO2_val_3[1] = !Values.F_Nan
;    SO2_val_4[0] = !Values.F_Nan
;    SO2_val_4[1] = !Values.F_Nan
;    SO_val[0] = !Values.F_Nan
;    SO_val[1] = !Values.F_Nan
;
;    cgplot, O2_fit_params.T_P_Shadow, O2_fit_params.Brightness, psym = 16, /overplot, color = 'red', $
;      ERR_YLOW = O2_fit_params.ERR_Brightness, ERR_YHigh = O2_fit_params.ERR_Brightness, ERR_XLOW = O2_fit_params.exptime/2., ERR_XHigh = O2_fit_params.exptime/2.
;    cgplot, O2_fit_params.T_P_Shadow, O2_fit_params.Brightness, color = 'red', /overplot
;
;    cgplot, O3_fit_params.T_P_Shadow, O3_fit_params.Brightness, psym = 15, /overplot, color = 'red', $
;      ERR_YLOW = O3_fit_params.ERR_Brightness, ERR_YHigh = O3_fit_params.ERR_Brightness, ERR_XLOW = O3_fit_params.exptime/2., ERR_XHigh = O3_fit_params.exptime/2., ERR_CLIP = 1
;    cgplot, O3_fit_params.T_P_Shadow, O3_fit_params.Brightness, color = 'red', /overplot
;
;    ;cgplot, O1_fit_params.T_P_Shadow, O1_fit_params.Brightness, psym = 16, /overplot, color = 'green', /err_clip, $
;    ;ERR_YLOW = O1_fit_params.ERR_Brightness, ERR_YHigh = O1_fit_params.ERR_Brightness, ERR_XLOW = O1_fit_params.exptime/2., ERR_XHigh = O1_fit_params.exptime/2.
;    ;cgplot, O1_fit_params.T_P_Shadow, O1_fit_params.Brightness, color = 'green', /overplot
;
;    cgplot, Na1Na2_fit_params.T_P_Shadow, Na1Na2_fit_params.Brightness, psym = 16, /overplot, color = 'orange', $
;      ERR_YLOW = Na1Na2_fit_params.ERR_Brightness, ERR_YHigh = Na1Na2_fit_params.ERR_Brightness, ERR_XLOW = Na1Na2_fit_params.exptime/2., ERR_XHigh = Na1Na2_fit_params.exptime/2.
;    cgplot, Na1Na2_fit_params.T_P_Shadow, Na1Na2_fit_params.Brightness, color = 'orange', /overplot
;
;    cgAxis, YAxis=1, YRange=[0, 6.5],title= 'Flux Density of Sulfur- Species [Jy]', COLOR = 'blue', /Save
;
;    cgplot, t_post_eclipse, (SO2_val_1 + SO2_val_2)/1000., psym = 17, /overplot, color = 'blue', $
;      ERR_YLOW = (SO2_err_1 + SO2_err_2)/1000., ERR_YHigh = (SO2_err_1 + SO2_err_2)/1000., ERR_XLOW = [0.,0.,1.7,1.84,5.68], ERR_XHigh = [0.,0.,1.3,1.81,4.9]
;    cgplot, t_post_eclipse, (SO2_val_1 + SO2_val_2)/1000., color = 'blue', /overplot
;    cgplot, [0, t_post_eclipse[2]], [(SO2_val_1in + SO2_val_2in)/1000., (SO2_val_1[2] + SO2_val_2[2])/1000.], color = 'blue', LineStyle=2, /overplot
;
;    ;cgplot, t_post_eclipse, SO2_val_2/1000., psym = 18, /overplot, color = 'blue', $
;    ;  ERR_YLOW = SO2_err_1/1000., ERR_YHigh = SO2_err_1/1000., ERR_XLOW = [0.,0.,1.7,1.84,5.68], ERR_XHigh = [0.,0.,1.3,1.81,4.9]
;    ;cgplot, t_post_eclipse, SO2_val_2/1000., color = 'blue', /overplot
;
;    ;cgplot, t_post_eclipse, SO2_val_3/1000., psym = 19, /overplot, color = 'blue', $
;    ;  ERR_YLOW = SO2_err_1/1000., ERR_YHigh = SO2_err_1/1000., ERR_XLOW = [0.,0.,1.7,1.84,5.68], ERR_XHigh = [0.,0.,1.3,1.81,4.9]
;    ;cgplot, t_post_eclipse, SO2_val_3/1000., color = 'blue', /overplot
;
;    ;cgplot, t_post_eclipse, SO2_val_4/1000., psym = 20, /overplot, color = 'blue', $
;    ;  ERR_YLOW = SO2_err_1/1000., ERR_YHigh = SO2_err_1/1000., ERR_XLOW = [0.,0.,1.7,1.84,5.68], ERR_XHigh = [0.,0.,1.3,1.81,4.9]
;    ;cgplot, t_post_eclipse, SO2_val_4/1000., color = 'blue', /overplot
;
;    cgplot, t_post_eclipse, (SO2_val_3 + SO2_val_4)/1000., psym = 19, /overplot, color = 'blue', $
;      ERR_YLOW = (SO2_err_3 + SO2_err_4)/1000., ERR_YHigh = (SO2_err_3 + SO2_err_4)/1000., ERR_XLOW = [0.,0.,1.7,1.84,5.68], ERR_XHigh = [0.,0.,1.3,1.81,4.9]
;    cgplot, t_post_eclipse, (SO2_val_3 + SO2_val_4)/1000., color = 'blue', /overplot
;    cgplot, [0, t_post_eclipse[2]], [(SO2_val_3in + SO2_val_4in)/1000., (SO2_val_3[2] + SO2_val_4[2])/1000.], color = 'blue', LineStyle=2, /overplot
;
;    cgplot, t_post_eclipse, SO_val/1000., psym = 16, /overplot, color = 'blue', $
;      ERR_YLOW = SO2_err_1/1000., ERR_YHigh = SO2_err_1/1000., ERR_XLOW = [0.,0.,1.7,1.84,5.68], ERR_XHigh = [0.,0.,1.3,1.81,4.9]
;    cgplot, t_post_eclipse, SO_val/1000., color = 'blue', /overplot
;
;    cgtext, 2.5, !Y.CRange[1]/4., 'Penumbral Eclipse', orientation = 90., color = 'white'
;
;    Case date of
;      'UT180320': cglegend, title = ['[O I] 6300'+cgsymbol('Angstrom'), '[O I] 6364'+cgsymbol('Angstrom'), 'Na D1 + D2', 'SO2 spw1 + spw2', 'SO2 spw5 + spw6', 'SO spw1 + spw3'], Psym = [16, 15, 16, 17,19, 16], charsize = 0.9, $
;        bg_color = 'white', color = ['red', 'red', 'orange', 'blue', 'blue', 'blue'], Location=[0.67, 0.87];, /Background, /BOX
;      'UT190812': AL_legend, ['K-D in Umbra 190812'], Psym = 4, /right, charsize = 1, /clear
;    endcase
;    cgps_close
;  endif


  ;========================================================================================================================================================
  ;===                                                                                                                                                  ===
  ;===                                                                                                                                                  ===
  ;===                                                        Keck HIRES Analysis                                                                       ===
  ;===                                                                                                                                                  ===
  ;===                                                                                                                                                  ===
  ;========================================================================================================================================================
  ;
  ;-------------------------------------------------------------Basic Reductions---------------------------------------------------------------------------
  if part eq 100.0 then begin
    ;  big_array = fltarr(2139, 4096, 10)
    ;  for i = 0, 9 do begin
    ;    big_array[*,*,i] = [mrdfits(Directory + 'j2970'+strcompress(130+i, /rem)+'.fits', 3, header, /fscale), $
    ;                        mrdfits(Directory + 'j2970'+strcompress(130+i, /rem)+'.fits', 2, header, /fscale), $
    ;                        mrdfits(Directory + 'j2970'+strcompress(130+i, /rem)+'.fits', 1, header, /fscale)]
    ;  endfor
    ;  bias = median(big_array, dim = 3, /even)
    big_array = fltarr(2139, 4096, 5)
    for i = 0, 4 do begin
      big_array[*,*,i] = [mrdfits(Directory + 'j2970'+strcompress(150+i, /rem)+'.fits', 3, header, /fscale), $
        mrdfits(Directory + 'j2970'+strcompress(150+i, /rem)+'.fits', 2, header, /fscale), $
        mrdfits(Directory + 'j2970'+strcompress(150+i, /rem)+'.fits', 1, header, /fscale)]
    endfor
    bias = median(big_array, dim = 3, /even)

    big_array = fltarr(2139, 4096, 5)
    for i = 0, 4 do begin
      big_array[*,*,i] = [mrdfits(Directory + 'j2970'+strcompress(140+i, /rem)+'.fits', 3, header, /fscale), $
        mrdfits(Directory + 'j2970'+strcompress(140+i, /rem)+'.fits', 2, header, /fscale), $
        mrdfits(Directory + 'j2970'+strcompress(140+i, /rem)+'.fits', 1, header, /fscale)]
    endfor
    junk = mrdfits(Directory + 'j2970140.fits', 0, header, /fscale)  ; get the proper header
    flat = median(big_array, dim = 3, /even) - bias                  ; WE WILL NORMALIZE / USE THIS FLAT LATER, ORDER BY ORDER
    MWRFITS, rotate(transpose(FLAT),7), directory+'Reduced\FLAT_BS.fits', header, /create, /silent ;super flat, bias subtracted but ***NOT NORMALIZED***
    flat = flat / mean(flat[where(flat gt 2000.)])                   ; normalize the flat, this seems to work better than order by order flat fielding (2000 is the whole aperture)
    ;     window, 0, xs = 4096, ys = 2200
    ;     tv, bytscl(rotate(transpose(reform(flat)),7), 0, 2)
    ;    stop

    ; ---> LAME!!! This flat mis-aligned. There are no other flats.
    ; There are real pixel level artifacts if we shift it, so there's no way around loosing a few pixels of spatial data on each end of the slit.
    ; Use the 14" slit decker in the future. Because of this misalignement, the 7" slit decker does not leave many pixels fo sky sampling.

    big_array = fltarr(2139, 4096, 5)
    for i = 0, 4 do begin
      big_array[*,*,i] = [mrdfits(Directory + 'j2970'+strcompress(145+i, /rem)+'.fits', 3, header, /fscale), $
        mrdfits(Directory + 'j2970'+strcompress(145+i, /rem)+'.fits', 2, header, /fscale), $
        mrdfits(Directory + 'j2970'+strcompress(145+i, /rem)+'.fits', 1, header, /fscale)]
    endfor
    junk    = mrdfits(Directory + 'j2970145.fits', 0, header, /fscale)
    ThAr    = median(big_array, dim = 3, /even) - bias
    MWRFITS, rotate(transpose(ThAr),7), directory+'Reduced\ThAr.fits', header, /create, /silent

    Jupiter_array  = fltarr(2139, 4096, n_elements(Jupiter_frames))
    for i = 0, n_elements(Jupiter_frames)-1 do begin
      filename = 'j2970'+strcompress(Jupiter_frames[i], /rem)+'.fits'
      Jupiter_array[*,*,i] = [mrdfits(Directory + filename, 3, header, /fscale), $
        mrdfits(Directory + filename, 2, header, /fscale), $
        mrdfits(Directory + filename, 1, header, /fscale)] - bias
      junk                 =  mrdfits(Directory + filename, 0, header, /fscale)
      new_filename         = STRMID(filename, 0, strpos(filename,'.fits'))
      Jupiter_array[*,*,i] = Jupiter_array[*,*,i] / flat
      Jupiter_array[*,*,i] = Jupiter_array[*,*,i] / float(Sxpar(header, 'EXPTIME')) ;normalize to 1 second exposure time
      write_file = rotate(transpose(reform(Jupiter_array[*,*,i])),7)
      MWRFITS, write_file, directory+'Reduced\' + new_filename + '.Cleaned.fits', header, /create, /silent
    endfor

    Ganymede = [mrdfits(Directory + 'j2970302.fits', 3, header, /fscale), $
      mrdfits(Directory + 'j2970302.fits', 2, header, /fscale), $
      mrdfits(Directory + 'j2970302.fits', 1, header, /fscale)]
    junk     =  mrdfits(Directory + 'j2970302.fits', 0, header, /fscale)
    Ganymede = Ganymede - bias ; Don't flat divide, since we'll want to use Ganymede to find the trace in each order, and the misaligned flat would throw off the fit.
    Ganymede = Ganymede / float(Sxpar(header, 'EXPTIME')) ; normalize to 1 second expsosure time
    MWRFITS, rotate(transpose(Ganymede),7), directory+'Reduced\j2970302.Trace.fits', header, /create, /silent

    ;---------------Reduce Io frames --------------------------------------------
    READCOL,'D:\DATA\___Calibration___\Carl_centrifugal_equator.txt', torus_deg, torus_lat_out, skipline = 1, /Silent
    cspice_UTC2ET, penumbra_utc, PenUmbra_ET
    cspice_UTC2ET, umbra_utc, Umbra_ET
    Io_array        = fltarr(2139, 4096, n_elements(Io_frames))
    ET_array        = dblarr(N_elements(Io_frames))
    exptime_array   = fltarr(N_elements(Io_frames))
    torus_lat_array = fltarr(N_elements(Io_frames))
    for i = 0, n_elements(Io_frames)-1 do begin
      filename = 'j2970'+strcompress(Io_frames[i], /rem)+'.fits'
      Io_array[*,*,i] = [mrdfits(Directory + filename, 3, header, /fscale), $
        mrdfits(Directory + filename, 2, header, /fscale), $
        mrdfits(Directory + filename, 1, header, /fscale)] - bias
      junk            =  mrdfits(Directory + filename, 0, header, /fscale)
      new_filename    = STRMID(filename, 0, strpos(filename,'.fits'))

      Io_array[*,*,i] = Io_array[*,*,i] / flat
      Io_array[*,*,i] = Io_array[*,*,i] / float(Sxpar(header, 'EXPTIME')) ; normalize to 1 second expsosure time

      ; find the instantaneous Earth-Io Doppler Shift
        cspice_UTC2ET, sxpar(header, 'DATE_BEG'), ET
        ET_mid_exposure = ET + float(sxpar(header, 'EXPTIME'))/2.
        cspice_spkezr, 'Io', ET_mid_exposure, 'J2000', 'LT+S', 'Earth', Io_Earth_State, ltime
        theta  = cspice_vsep(Io_Earth_state[0:2], Io_Earth_state[3:5])
        Io_wrt_Earth_Dopplershift = cos(theta) * norm(Io_Earth_State[3:5])

      cspice_subpnt, 'Near point: ellipsoid', 'Jupiter', ET_mid_exposure - ltime, 'IAU_Jupiter', 'None', 'Io', Sub_Io, trgepc, srfvec ;get sub solar point in cartesian IAU Jupiter coords
      cspice_bodvrd, 'Jupiter', 'RADII', 3, radii ;get Jupiter shape constants
      re = radii[0]
      rp = radii[2]
      f = (re-rp)/re
      obspos = Sub_Io - srfvec
      cspice_recpgr, 'Jupiter', obspos, re, f, Io_SysIII, Io_SysIII_LATITUDE, opgalt
      torus_lat = interpol(torus_lat_out, reverse(torus_deg), Io_SysIII*!radeg) ; Io's latitude in the torus using Phil Phipp's arrays
      torus_lat_array[i] = torus_lat
      ET_array[i]        = ET_mid_exposure
      exptime_array[i]   = float(sxpar(header, 'EXPTIME'))

      SXADDPAR, header, 'Sys3_Lat', Io_SysIII, 'Io''s System III Latitude'
      SXADDPAR, header, 'Sys3_Lon', Io_SysIII_LATITUDE, 'Io''s System III Longitude'
      SXADDPAR, header, 'Torus_Lat', torus_lat, 'JRM09 Dipole Approx'
      SXADDPAR, header, 'Io_DOPPL', Io_wrt_Earth_Dopplershift, 'Io-Earth V_radial in km/s (mid exposure)'
      SXADDPAR, header, 'T_PSHADO', (ET_mid_exposure-PenUmbra_ET) / 60., 'Minutes since Penumbral ingress'
      SXADDPAR, header, 'T_USHADO', (ET_mid_exposure-Umbra_ET) / 60., 'Minutes since Umbral ingress'

      write_file = rotate(transpose(reform(Io_array[*,*,i])),7)
      MWRFITS, write_file, directory+'Reduced\' + new_filename + '.Cleaned.fits', header, /create, /silent
    endfor
    Io_Airglow_params.penumbra_ET    = PenUmbra_ET
    Io_Airglow_params.T_P_Shadow     = (ET_array - PenUmbra_ET)/60.
    Io_Airglow_params.torus_latitude = torus_lat_array
    Io_Airglow_params.exptime        = exptime_array
    save, Io_Airglow_params, filename = Directory + 'Reduced\Io_Airglow_params.sav'
  endif ; Part 0 (basic bias/flat/exptime reductions)

  order_66 = {guess_coeffs:[5328.86, 0.0260124,-5.03615e-007], low_bound: 64,  hi_bound:130,  WL_range:[5575., 5580], aperture_limit:[12,27], name:'order_66'}  ; best possible
  order_65 = {guess_coeffs:[5410.92, 0.0263225,-4.90198e-007], low_bound:101,  hi_bound:170,  WL_range:[5575., 5580], aperture_limit:[12,27], name:'order_65'}  ; best possible
  order_64 = {guess_coeffs:[5495.41, 0.0267302,-4.83544e-007], low_bound:140,  hi_bound:210,  WL_range:[5575., 5580], aperture_limit:[12,27], name:'order_64'}  ; best possible
  order_63 = {guess_coeffs:[5582.56, 0.0272323,-5.08754e-007], low_bound:179,  hi_bound:247,  WL_range:[5575., 5580], aperture_limit:[12,27], name:'order_63'}  ; best possible
  order_62 = {guess_coeffs:[5672.57, 0.0277664,-5.52980e-007], low_bound:219,  hi_bound:286,  WL_range:[5575., 5580], aperture_limit:[12,27], name:'order_62'}  ; best possible
  order_61 = {guess_coeffs:[5765.50, 0.0281711,-5.24588e-007], low_bound:261,  hi_bound:328,  WL_range:[5575., 5580], aperture_limit:[12,27], name:'order_61'}  ; best possible
  order_60 = {guess_coeffs:[5861.58, 0.0286144,-5.19316e-007], low_bound:308,  hi_bound:377,  WL_range:[5888., 5898], aperture_limit:[12,26], name:'order_60'}  ; best possible, aperture double checked
  order_59 = {guess_coeffs:[5960.78, 0.0292072,-5.44596e-007], low_bound:349,  hi_bound:418,  WL_range:[5888., 5898], aperture_limit:[12,25], name:'order_59'}  ; best possible
  order_58 = {guess_coeffs:[6063.57, 0.0296808,-5.46829e-007], low_bound:400,  hi_bound:471,  WL_range:[5888., 5898], aperture_limit:[12,25], name:'order_58'}  ; best possible
  order_57 = {guess_coeffs:[6169.91, 0.0302184,-5.54848e-007], low_bound:448,  hi_bound:519,  WL_range:[5888., 5898], aperture_limit:[11,25], name:'order_57'}  ; best possible
  order_56 = {guess_coeffs:[6280.03, 0.0307753,-5.60140e-007], low_bound:496,  hi_bound:568,  WL_range:[6300., 6301.6], aperture_limit:[11,25], name:'order_56'}  ; best possible
  order_55 = {guess_coeffs:[6394.13, 0.0314351,-5.94222e-007], low_bound:546,  hi_bound:621,  WL_range:[9508., 9680], aperture_limit:[11,25], name:'order_55'}  ; best possible
  order_54 = {guess_coeffs:[6512.50, 0.0320107,-5.93035e-007], low_bound:600,  hi_bound:675,  WL_range:[9508., 9680], aperture_limit:[11,24], name:'order_54'}  ; best possible
  order_53 = {guess_coeffs:[6635.47, 0.0323705,-5.09817e-007], low_bound:654,  hi_bound:729,  WL_range:[9508., 9680], aperture_limit:[11,24], name:'order_53'}  ; best possible
  order_52 = {guess_coeffs:[6763.03, 0.0332943,-6.24047e-007], low_bound:720,  hi_bound:797,  WL_range:[9508., 9680], aperture_limit:[11,24], name:'order_52'}  ; best possible
  order_51 = {guess_coeffs:[6895.62, 0.0338623,-6.04324e-007], low_bound:786,  hi_bound:864,  WL_range:[9508., 9680], aperture_limit:[11,24], name:'order_51'}  ; best possible
  order_50 = {guess_coeffs:[7033.31, 0.0348218,-6.87466e-007], low_bound:848,  hi_bound:926,  WL_range:[9508., 9680], aperture_limit:[10,24], name:'order_50'}  ; best possible
  order_49 = {guess_coeffs:[7176.80, 0.0355404,-6.96667e-007], low_bound:911,  hi_bound:988,  WL_range:[9508., 9680], aperture_limit:[10,24], name:'order_49'}  ; best possible
  order_48 = {guess_coeffs:[7326.22, 0.0363962,-7.37812e-007], low_bound:977,  hi_bound:1053, WL_range:[9508., 9680], aperture_limit:[10,24], name:'order_48'}  ; best possible
  order_47 = {guess_coeffs:[7482.41, 0.0368767,-6.95045e-007], low_bound:1045, hi_bound:1123, WL_range:[9508., 9680], aperture_limit:[10,24], name:'order_47'} ; Not great
  order_46 = {guess_coeffs:[7645.07, 0.0376087,-6.98756e-007], low_bound:1116, hi_bound:1248, WL_range:[7660., 7705], aperture_limit:[10,24], name:'order_46'} ; best possible, dodgy >7730, aperture double checked
  order_46 = {guess_coeffs:[7645.07, 0.0376087,-6.98756e-007], low_bound:1116, hi_bound:1248, WL_range:[7770., 7778], aperture_limit:[10,24], name:'order_46'} ; best possible, dodgy >7730, aperture double checked
  order_45 = {guess_coeffs:[7814.95, 0.0384298,-6.97203e-007], low_bound:1191, hi_bound:1271, WL_range:[9508., 9680], aperture_limit:[10,24],  name:'order_45'} ; best possible
  order_44 = {guess_coeffs:[7992.28, 0.0396213,-7.90228e-007], low_bound:1267, hi_bound:1350, WL_range:[8070., 8100], aperture_limit:[10,24],  name:'order_44'} ; best possible
  order_43 = {guess_coeffs:[8178.53, 0.0400097,-6.47557e-007], low_bound:1349, hi_bound:1430, WL_range:[8179., 8201], aperture_limit:[10,24],  name:'order_43'} ; best possible,
  order_42 = {guess_coeffs:[8372.70, 0.0413697,-7.96743e-007], low_bound:1449, hi_bound:1532, WL_range:[8440., 8453], aperture_limit:[10,24],  name:'order_42'} ; DODGY WL solution! ok for 8446
  order_41 = {guess_coeffs:[8576.86, 0.0424576,-8.44546e-007], low_bound:1537, hi_bound:1625, WL_range:[9508., 9680], aperture_limit:[9,24],  name:'order_41'} ; best possible
  order_40 = {guess_coeffs:[8791.59, 0.0433164,-8.47361e-007], low_bound:1631, hi_bound:1720, WL_range:[9508., 9680], aperture_limit:[9,24],  name:'order_40'} ; Not great
  order_39 = {guess_coeffs:[9016.92, 0.0444981,-8.89124e-007], low_bound:1731, hi_bound:1821, WL_range:[9508., 9680], aperture_limit:[9,24],  name:'order_39'} ; best possible
  order_38 = {guess_coeffs:[9254.72, 0.0452545,-8.55026e-007], low_bound:1837, hi_bound:1927, WL_range:[9508., 9680], aperture_limit:[9,24],  name:'order_38'} ; best possible
  order_37 = {guess_coeffs:[9504.39, 0.0469215,-9.89839e-007], low_bound:1948, hi_bound:2043, WL_range:[9508., 9680], aperture_limit:[9,24],  name:'order_37'} ; SO Band ~9553? C I 9658.4?

  orders   = [order_37, order_38, order_39, order_40, order_41, order_42, order_43, order_44, order_45, order_46, $
              order_47, order_48, order_49, order_50, order_51, order_52, order_53, order_54, order_55, order_56, $
              order_57, order_58, order_59, order_60, order_61, order_62, order_63, order_64, order_65, order_66]

  if part eq 101. then begin ; flatten, extract and get the wavelength solutions
    READCOL,'C:\IDL\Io\Keck Programs\thar_uves.dat', F='X,F', ThAr_WL, STRINGSKIP = '#', skipline = 1200, numline = 1200

    for h = 0, N_elements(orders) - 1 do begin
      order = orders[h]
      
      ; Use Ganymede to find the trace in within the orders of interest
        G    = mrdfits(Directory + 'Reduced\j2970302.Trace.fits', 0, G_header) ; do not use multiple Ganymede frames here, it's important it's fixed.
        ThAr = mrdfits(Directory + 'Reduced\ThAr.fits', 0, ThAr_header)
        flat = mrdfits(Directory + 'Reduced\FLAT_BS.fits', 0, Flat_header)     ; not yet normalized to unity

      ; Guess and check at a linear wavelength solution
        WL          = poly(findgen(4001), order.guess_coeffs)    ; Only roughly accurate
        WL_0        = WL
        xr          = minmax(WL)
        order_lines = ThAr_WL[where( (ThAr_WL gt xr[0]) and (ThAr_WL lt xr[1]))]
        ID          = make_array(N_elements(order_lines), value = 1.e5)

      cube          = fltarr(4001, 1 + order.aperture_limit[1] - order.aperture_limit[0], N_elements(Io_frames))
      Jupiter_cube  = fltarr(4001, 1 + order.aperture_limit[1] - order.aperture_limit[0], N_elements(Jupiter_frames))

      Frames = 'j2970'+strcompress(Io_frames, /rem)+'.Cleaned.fits'
      for frame = 0, n_elements(Frames)-1 do begin

        Io   = mrdfits(Directory + 'Reduced\' + Frames[frame], 0, header)

        ; -------------------------------------- Wavelength Solution in the Orders of Interest ------------------------------------------------------------------
        ;              Take_a_long_hard_look   = mrdfits(Directory + 'Reduced\' + Frames[-2], 0, header)
        ;              window, 0, xs = 4096, ys = 1000;,ys = 2050
        ;              tv, bytscl(Take_a_long_hard_look, 3, 8 )
        ;              rdpix, Take_a_long_hard_look
        ;              stop

        ; trace each order's footprinnt on the chip using the ecentroid of the sunlit ganymede spectrum
        subframe = G[0:4000, order.low_bound:order.hi_bound]
        s = size(subframe, /dim)
        trace = fltarr(s[0])
        for i = 0, s[0] - 1 do begin
          junk = max(subframe[i,*], loc)
          trace[i] = loc
        endfor

        trace_old = trace
        replace_left = trace[0:2000]
        replace_right= trace[2001:*]
        replace_left[where(replace_left gt 41.)] = !values.F_nan
        replace_right[where( (replace_right lt 30.) or (replace_right gt 100.) )] = !values.F_nan
        trace[0:2000] = replace_left
        trace[2001:*] = replace_right
        keep = where(finite(trace))
        x = findgen(s[0])
        coeffs = poly_fit(x[keep], trace[keep], 3) ;trace ganymede position with a 3rd order polynomial

        window, 1, xs = 1800, ys=800, title = 'WAVELENGTH SOLUTION FOR: ' + ORDER.NAME
        cgplot, findgen(s[0]), trace_old, psym = 4, /ynozero
        cgplot, findgen(s[0]), trace, psym = 4, /overplot, color = 'red'
        cgplot, x, POLY( findgen(s[0]), coeffs), /overplot, color = 'blue'

        ThAr_order_straight = fltarr(s[0], 30)
        Flat_order_straight = fltarr(s[0], 30)
        g_order_straight = fltarr(s[0], 30)
        Io_order_straight = fltarr(s[0], 30)
        for i = 0, s[0] - 1 do begin
          G_order_straight[i,*] = interpolate(G[i,*], order.low_bound + findgen(30) - 15. + POLY( i, coeffs))
          ThAr_order_straight[i,*] = interpolate(ThAr[i,*], order.low_bound + findgen(30) - 15. + POLY( i, coeffs))
          flat_order_straight[i,*] = interpolate(flat[i,*], order.low_bound + findgen(30) - 15. + POLY( i, coeffs))
          Io_order_straight[i,*] = interpolate(Io[i,*], order.low_bound + findgen(30) - 15. + POLY( i, coeffs))
        endfor

        Flat_aperture = Flat_order_straight[*, order.aperture_limit[0]:order.aperture_limit[1]]
        ;nomalize_with = mean(Flat_aperture, dim = 2, /NAN)
        coeffs = poly_fit(findgen(s[0]), mean(Flat_aperture, dim = 2, /NAN), 3)
        nomalize_with = poly(findgen(s[0]), coeffs)
        Flat_aperture = Flat_aperture / rebin(nomalize_with, s[0],order.aperture_limit[1] - order.aperture_limit[0] + 1) ; Normalize the flat field to unity, now its spectral AND spatially normalized
        aperture      = Io_order_straight[*, order.aperture_limit[0]:order.aperture_limit[1]]
        ;aperture      = aperture / Flat_aperture ; don't flatten twice!!!
        ;        if order.name eq 'order_44' then begin
        ;          window, 0, xs  = 1900, ys = 20
        ;          tv, bytscl(aperture[2100:*, *], 15, 25)
        ;          window, 1, xs  = 1900, ys = 20, ypos = 60
        ;          tv, bytscl(Flat_aperture[2100:*, *], 0.9, 1.1)
        ;          window, 4, xs  = 1900, ys = 20, ypos = 120
        ;          test = aperture / Flat_aperture^1.2
        ;          tv, bytscl(test[2100:*, *], 15, 25)
        ;          stop
        ;        endif

        ; Evaluate the guess at the wavelength solution
          ThAr_order_straight[where(ThAr_order_straight lt 0., /null)] = !Values.F_Nan
          ThAr_measured = total(ThAr_order_straight[*, order.aperture_limit[0]:order.aperture_limit[1]], 2)
          P = cglayout([1,2], ygap = 0.)
          cgplot, WL, ThAr_measured, /ylog, yr = [2.e3, 5.e5], /xstyle, pos = p[*,0], xtickformat = '(A1)', ytitle = 'ThAr Lamp Counts'

        ; find the peaks and do a 3rd order wavelength solution in each order
          identwave_1 = [] & identpixl_1 = []
          search    = 50  ; pixel distance from the expected position based of the fits header to search for an arc line
          expected_pixel = round( interpol( findgen(N_elements(WL)), WL, order_lines ) )
          parinfo = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.,0.]}, 4)
          parinfo[2].fixed = 1
          parinfo[2].value = 4.0
          for i = 0, n_elements(order_lines) - 1 do begin
            if ( (expected_pixel[i]-search lt 0) or (expected_pixel[i]+search gt (s[0])-1) ) then continue
            y = ThAr_measured[expected_pixel[i]-search:expected_pixel[i]+search]
            if total(finite(y), /Nan) eq 0. then continue
            result = mpfitpeak(findgen(search*2. + 1), y, a, /POSITIVE, PARINFO = PARINFO, STATUS = STATUS, NFREE = 2, /nan)
            if (a eq !null) then continue
            if (not finite(total(a))) then continue ;else print, a
            if ( (status gt 0) and (a[0] gt 1.e4) and (a[1] gt 0.) and (a[2] gt 2.) and (a[2] lt 6.)) then begin
              identwave_1 = [identwave_1, order_lines[i]]
              identpixl_1 = [identpixl_1, a[1] - search + expected_pixel[i]]
            endif
          endfor
          coeff_1 = ROBUST_POLY_FIT(identpixl_1, identwave_1, 3, yfit_1, SIG)
          WL = poly(findgen(N_elements(WL)), coeff_1)

        ; iterate and try for 5th order
          identwave_2 = [] & identpixl_2 = []
          search    = 10  ; pixel distance from the expected position based of the fits header to search for an arc line
          expected_pixel = round( interpol( findgen(N_elements(WL)), WL, order_lines ) )
          for i = 0, n_elements(order_lines) - 1 do begin
            if ( (expected_pixel[i]-search lt 0) or (expected_pixel[i]+search gt s[0]) ) then continue
            result = mpfitpeak(findgen(search*2. + 1), float(ThAr_measured[expected_pixel[i]-search:expected_pixel[i]+search]), a, /POSITIVE, PARINFO = PARINFO, STATUS = STATUS, NFREE = 2)
            ;print, a
            if not finite(total(a)) then continue
            if ( (status gt 0) and (a[0] gt 1.e3) and (a[1] gt 0.) and (a[2] gt 2.) and (a[2] lt 6.)) then begin
              identwave_2 = [identwave_2, order_lines[i]]
              identpixl_2 = [identpixl_2, a[1] - search + expected_pixel[i]]
            endif
          endfor
          coeff_2 = ROBUST_POLY_FIT(identpixl_2, identwave_2, 5, yfit_2, SIG)

        ; Manual override, when this technique fails
          if order.name eq 'order_42' then coeff_2 = order.GUESS_COEFFS

        cgplot, order_lines, make_array(N_elements(order_lines), value = 1.e5), /overplot, psym=4                           ; plot cataloged lines
        cgplot, interpol(WL_0, findgen(N_elements(WL_0)), identpixl_1), make_array(N_elements(identpixl_1), value = 2.e5), /overplot, psym=5, color = 'red' ; plot which lines were found in iteration 1
        cgplot, interpol(WL_0, findgen(N_elements(WL_0)), identpixl_2), make_array(N_elements(identpixl_2), value = 3.e5), /overplot, psym=5, color = 'blue'; plot which lines were found in iteration 2

        cgplot, identwave_1, identpixl_1, pos = p[*,1], /noerase, XR = XR, /xstyle, psym=5, color = 'red'
        cgplot, identpixl_1, yfit_1, psym = 5, /overplot, color = 'red'
        cgplot, identwave_2, identpixl_2, psym=5, /overplot, color = 'blue'
        cgplot, identpixl_2, yfit_2, psym = 5, /overplot, color = 'blue'
        cgplot, poly(findgen(N_elements(WL)), order.guess_coeffs), findgen(N_elements(WL)), /overplot
        cgplot, poly(findgen(N_elements(WL)), coeff_1), findgen(N_elements(WL)), /overplot, color = 'red'
        cgplot, poly(findgen(N_elements(WL)), coeff_2), findgen(N_elements(WL)), /overplot, color = 'blue'

        Print, 'Fitting Wavelength Solution From'+ strcompress(N_elements(identpixl_2)) + ' ThAr lines.
        Print, 'Guess coefficients for [xstart, dispersion/pixel, 2nd order]', order.guess_coeffs
        Print, 'A better guess would have been:', ROBUST_POLY_FIT(identpixl_2, identwave_2, 2, junk, SIG)
        print, 'Poly fit coefficients:', coeff_2
        WL = poly(findgen(N_elements(WL)), coeff_2)

        ; write the fits files and run the Cosmic Ray Correcctions. Some orders overlap the CCD edge in the extraction so exclude these regions in the CR correction
          SXADDPAR, Header, 'BZERO', 0
          SXADDPAR, Header, 'BSCALE', 0
          MWRFITS, aperture, directory + 'Reduced\Cosmic Rays\CR_' + Frames[frame], header, /create
          Case 1 of
            order.name eq 'order_43': statsec = '[0:2925,*]'
            order.name eq 'order_52': statsec = '[500:4000,*]'
            order.name eq 'order_53': statsec = '[0:2620,*]'
            else: junk = temporary(statsec)
          endcase
          la_cosmic, directory + 'Reduced\Cosmic Rays\CR_' + Frames[frame], outsuff = "CR", sigclip = 4.5, statsec = statsec
          cube[*,*,frame] = mrdfits(Directory + 'Reduced\Cosmic Rays\CR_j2970'+strcompress(Io_frames[frame], /rem)+'.CleanedCR.fits', 0, junk_header)
      endfor ; frame (Io frames)

      ; Now do the same process for Jupiter
      Frames = 'j2970'+strcompress(Jupiter_frames, /rem)+'.Cleaned.fits'
      for frame = 0, n_elements(Frames)-1 do begin

        Jupiter  = mrdfits(Directory + 'Reduced\' + Frames[frame], 0, header)
        subframe = G[0:4000, order.low_bound:order.hi_bound]
        s = size(subframe, /dim)
        trace = fltarr(s[0])
        for i = 0, s[0] - 1 do begin
          junk = max(subframe[i,*], loc)
          trace[i] = loc
        endfor

        trace_old     = trace
        replace_left  = trace[0:2000]
        replace_right = trace[2001:*]
        replace_left[where(replace_left gt 45.)] = !values.F_nan
        replace_right[where( (replace_right lt 30.) or (replace_right gt 100.) )] = !values.F_nan
        trace[0:2000] = replace_left
        trace[2001:*] = replace_right

        keep = where(finite(trace))
        x = findgen(s[0])
        coeffs = poly_fit(x[keep], trace[keep], 2)
        Jupiter_order_straight = fltarr(s[0], 30)
        ;Flat_order_straight = fltarr(s[0], 30)
        for i = 0, s[0] - 1 do begin
          Jupiter_order_straight[i,*] = interpolate(Jupiter[i,*], order.low_bound + findgen(30) - 15. + POLY( i, coeffs))
          ;flat_order_straight[i,*]    = interpolate(flat[i,*], order.low_bound + findgen(30) - 15. + POLY( i, coeffs))
        endfor

        ;        ; prepare the field field
        ;          Flat_aperture = Flat_order_straight[*, order.aperture_limit[0]:order.aperture_limit[1]]
        ;          ;nomalize_with = mean(Flat_aperture, dim = 2, /NAN)
        ;          coeffs = poly_fit(findgen(s[0]), mean(Flat_aperture, dim = 2, /NAN), 3)
        ;          nomalize_with = poly(findgen(s[0]), coeffs)
        ;          Flat_aperture = Flat_aperture / rebin(nomalize_with, s[0],order.aperture_limit[1] - order.aperture_limit[0] + 1) ; Normalize the flat field to unity, now its spectral AND spatially normalized

        Jupiter_cube[*,*,frame] = Jupiter_order_straight[*, order.aperture_limit[0]:order.aperture_limit[1]] ;/ Flat_aperture ; flat field jupiter
      endfor ; frames (Jupiter frame number)
      save, cube, Jupiter_cube, WL, filename = Directory + 'Reduced\' + order.name + '.sav'
    endfor ; h (order number)
  endif

  if Part eq 102. then begin ;

    ;---------------------------------------Find Io's spatial location---------------------------------------
    order           = order_56
    restore, Directory + 'Reduced\' + order.name + '.sav'
    dummy           = Cube
    EXTRACT_ROWS    = intarr(2, N_elements(Io_frames ) )

    loadct, 0

    for f = 0, N_elements(Io_frames )-1 do begin
      Chunk_6300 = Cube[681:691,*,f]
      Chunk_BG = Cube[670:680,*,f]
      Io_profile = total(Chunk_6300 - Chunk_BG, 1)
      fit = mpfitpeak(findgen(N_elements(Io_profile)), Io_profile, a)
      window, 2
      cgplot, io_profile
      cgplot, fit, color = 'red', /overplot

      ind        = sort(fit)
      sky_BG_ind = ind[0:3]   ; largest extraction regions over Io seem to show the best results, but we still need some sky to sample at both ends
      Io_ind     = ind[4:*]
      ;sky_BG_ind = ind[0:6]  ; test for a narrow extraction region ---> noisey
      ;Io_ind     = ind[7:*]

      EXTRACT_ROWS[*,f] = minmax(Io_ind)

      ; Manual override, we want to always have some of both sides of the slitt
        if EXTRACT_ROWS[0,f] eq 0 then EXTRACT_ROWS[0,f] = 1
        ;if EXTRACT_ROWS[1,f] eq max(ind) then EXTRACT_ROWS[1,f] = max(ind) - 1
        if EXTRACT_ROWS[1,f] eq (max(ind)-1) then EXTRACT_ROWS[1,f] = max(ind) - 2 ; test shouldn't need this but does it fix 51-54?
        print, 'Extract Io over rows between:', EXTRACT_ROWS[*,f]

      window, 0, xs = 1901, ys = 300
      cgimage, reform(Cube[0000:1900,*,f]), minvalue = 0.5*mean(reform(Cube[0000:1900,*,f])), 1.5*mean(reform(Cube[0000:1900,*,f]))
      dummy[*,EXTRACT_ROWS[0,f]:EXTRACT_ROWS[1,f],f] = !values.F_Nan

      window, 1, xs = 1900, ys = 300
      cgimage, reform(dummy[0000:1900,*,f]), minvalue = 0.5*mean(reform(Cube[0000:1900,*,f])), 1.5*mean(reform(Cube[0000:1900,*,f]))
      ;      wait, 0.5
    endfor
    save, EXTRACT_ROWS, filename = Directory + 'Reduced\EXTRACT_ROWS.sav'

    ;---------------------------inspect the extraction in 2D----------------------------------------------------------------
    order           = order_43 ; order 40 is broken
    restore, Directory + 'Reduced\' + order.name + '.sav'
    dummy           = Cube
    s               = size(cube, /dim)
    co_add_2d       = fltarr(s[0], s[1])
    for f = 0, N_elements(Io_frames )-1 do begin

      window, 0, xs = 3400, ys = 300
      cgimage, reform(Cube[0000:3400,*,f]), minvalue = 0.5*mean(reform(Cube[0000:3400,*,f])), 1.5*mean(reform(Cube[0000:3400,*,f]))
      dummy[*,EXTRACT_ROWS[0,f]:EXTRACT_ROWS[1,f],f] = !values.F_Nan
      Background = total(dummy[*, *, f], 2, /nan)
      Case 1 of
        order.name eq 'order_43': Background = Background / mean(Background[0:2925])
        order.name eq 'order_52': Background = Background / mean(Background[500:*])
        order.name eq 'order_53': Background = Background / mean(Background[0:2620])
        else: Background = Background / mean(Background)
      endcase
      img = reform(dummy[*,*,f])
      for i = 0, s[0]-1 do begin
        img[i,*] = interpol(img[i,*], findgen(s[1]), findgen(s[1]), /Nan)
      endfor
      dummy[*,*,f] = img
      resid = reform(Cube[*,*,f]) - reform(dummy[*,*,f])
      Case 1 of
        order.name eq 'order_43': statsec = resid[0:2925,*]
        order.name eq 'order_52': statsec = resid[500:*,*]
        order.name eq 'order_53': statsec = resid[0:2620,*]
        else: statsec = resid
      endcase
      slit_profile = mean(statsec, dimension = 1, /NAN)
      fake_profile = rebin(transpose(slit_profile), s[0], s[1]) * rebin(Background, s[0], s[1])

      window, 1, xs = 3400, ys = 300
      cgimage, reform(dummy[0000:3400,*,f]), minvalue = 0.5*mean(reform(Cube[0000:3400,*,f])), 1.5*mean(reform(Cube[0000:3400,*,f]))
      co_add_this_frame_2d = resid[0000:3400,*] - fake_profile[0000:3400,*]
      co_add_2d                       = co_add_2d + co_add_this_frame_2d
      window, 3
      cghistoplot, co_add_this_frame_2d, omin = -10., omax = 10.
      window, 4, xs = 3400, ys = 300, ypos = 350
      cgimage, co_add_2d
      stop
    endfor
  endif

  if Part eq 103. then begin ; Make waterfall plots in each order
    restore, Directory + 'Reduced\Io_Airglow_params.sav'

    ;---------------------------- Determine Sensitivy for flux calibration -------------------------------------

    ; Compare publications of Jupiter's spectral albedo at disk center

    ; absolute brightness: Digitized Plot from Woodman et al. 1979.
      READCOL,'C:\IDL\Io\Woodman_et_al_1979_plot1.txt', F='A,A', WL_1, Albedo_1, STRINGSKIP = '#', /Silent;wavelength in Angstroms I/F unitless
      READCOL,'C:\IDL\Io\Woodman_et_al_1979_plot2_new.txt', F='A,A', WL_2, Albedo_2, STRINGSKIP = '#', /Silent ;wavelength in Angstroms I/F unitless
      Woodman_WL = float([WL_1, WL_2])             ; STITCH THESE TOGETHER
      Woodman_Albedo = Float([albedo_1, albedo_2]) ; STITCH THESE TOGETHER

    ; absolute brightness: from Karkoschka (1998) Icarus on the PDS as ID # ESO-J/S/N/U-SPECTROPHOTOMETER-4-V1.0
      READCOL,'C:\IDL\Io\Karkoschka_1995low.tab', F='X,A,X,A', Karkoschka_WL, Karkoschka_Albedo, STRINGSKIP = '#', /Silent ;wavelength in nm I/F unitless
      READCOL,'C:\IDL\Io\Karkoschka_1995high.tab', F='X,A,X,A', Karkoschka_WL_2, Karkoschka_Albedo_2, STRINGSKIP = '#', /Silent ;wavelength in Angstroms I/F unitless
      Karkoschka_WL = float([Karkoschka_WL, Karkoschka_WL_2])             ; STITCH THESE TOGETHER
      Karkoschka_Albedo = Float([Karkoschka_albedo, Karkoschka_albedo_2]) ; STITCH THESE TOGETHER
      Karkoschka_Albedo = Karkoschka_Albedo[sort(Karkoschka_WL)]
      Karkoschka_WL = Karkoschka_WL[sort(Karkoschka_WL)]

    ; compare the two
      cgplot, Woodman_WL / 10., Woodman_Albedo, color = 'blue', xstyle = 1., psym = 3, Xtitle = 'Wavelength (nm)', ytitle = 'I/F Reflectivity'
      cgplot, Karkoschka_WL, Karkoschka_Albedo*1.35, color = 'red', /overplot
      cgtext, 340, .1, 'EQUATOR AT CENTRAL MERIDIAN (Woodman et al. 1979)', color = 'blue'
      cgtext, 340, .16, 'FULL DISK scaled by 1.35 (Karkoschka 1998)', color = 'red'

    ; make an informed choice via scaling
      Karkoschka_wl = Karkoschka_wl * 10. ;nm to Angstroms
      Karkoschka_Albedo = Karkoschka_Albedo*1.35 ;Scale the "Full disk albedo" to the "Central Meridian Equatorial Absolute Reflectivity"

    ; Get the Solar Spectral Irradiance at 1 AU
      READCOL,'C:\IDL\Io\Kurucz_2005_irradthuwl.dat', F='A,A', WL_nm, flux, STRINGSKIP = '#', /Silent ;flux is in W/m2/nm
      start  = where(WL_nm eq '299.100')
      WL_nm = float(WL_nm[start:*])
      flux = float(flux[start:*])

    ; change flux units from W/m^2/nm to photons/(cm^2 s A)
    ; multiply by ((lambda / hc) / 1 W)  * (1 m^2 / 1e4 cm^2) * (1 nm / 10 A)
      conversion = ((WL_nm*1.e-9)/(6.62606957e-34*299792458.D)) * (1./1.e4) * (1./10.)
      flux = flux * conversion      ; Cross-checked this result against Huebner et al. 1992.
      WL_A = temporary(WL_nm) * 10. ; Wavelength from nm into angstroms
      VACTOAIR, WL_A, WL_A_Air      ; Vacuum to air wavelength conversion
      WL_A = temporary(WL_A_Air)

    Jupiter_center_header = headfits(directory+'Reduced\j2970'+strcompress(Jupiter_frames[0], /rem)+'.Cleaned.fits')
    cspice_UTC2ET, sxpar(Jupiter_center_header, 'DATE-OBS') + ' '+ sxpar(Jupiter_center_header, 'UTC'), ET
    cspice_spkezr, 'Jupiter', ET, 'J2000', 'LT+S', 'Sun', Jupiter_Sun_State, ltime
    cspice_spkezr, 'Jupiter', ET, 'J2000', 'LT+S', 'Earth', Jupiter_Earth_State, ltime
    solar_distance = norm(Jupiter_Sun_State[0:2]) / 149597871.
    flux_at_jupiter = flux / solar_distance^2.

    ; Multiply incident solar irradiance x spectral albedo to determine the theoreitcal brightness of Jupiter at disk center
      Albedo = INTERPOL(Karkoschka_Albedo, Karkoschka_WL, WL_A)
      Rayleighs_per_angstrom = 4.*flux_at_jupiter*albedo / 1.e6
      if keyword_set(debug) then window, Title = 'Instantaneous Rayleighs per Angstrom: Center of Jupiter''s Disk'
      if keyword_set(debug) then plot, WL_A, Rayleighs_per_angstrom, xr = [5885, 5900], charsize = 2 ;compare to 5.5 MR per angstrom by Brown & Schneider, 1981

    ; Adjust the expected absolute flux for Jupiter's Instantaneous Doppler shift
      theta  = cspice_vsep(Jupiter_Earth_State[0:2], Jupiter_Earth_State[3:5])
      Dopplershift = cos(theta) * norm(Jupiter_Earth_State[3:5])   ; scalar projection of the relative velocity along the line of sight
      WL_A = WL_A + WL_A*Dopplershift/cspice_clight()

    restore, Directory + 'Reduced\EXTRACT_ROWS.sav'

    ; get color versus ingress time & setup plot positions
      timeColors = BytScl(findgen(11), Min=0, Max=11, Top=11)
      Color      = timeColors[0]
      cgLoadCT, 33, NColors=12, /reverse

    for h = 0, N_elements(orders)-1 do begin
      order           = orders[h]
      restore, Directory + 'Reduced\' + order.name + '.sav'

      ;if not ((order.name eq 'order_60')) then continue
      case 1 of
        order.name eq 'order_42': BEGIN
          ;xr           = order.WL_range
          xr           = minmax(WL)
        end
        order.name eq 'order_43': BEGIN
          WL           = WL[0:3400]
          cube         = cube[0:3400,*,*]
          Jupiter_cube = Jupiter_cube[0:3400,*,*]
          xr           = order.WL_range
          ;xr           = minmax(WL)
          title        = 'Io''s 2018-08-07 Eclipse: Sodium Emissions'
        end
        order.name eq 'order_46': BEGIN
          xr           = order.WL_range
          title        = 'Io''s 2018-08-07 Eclipse: Potassium D Emissions'
        end
        order.name eq 'order_52': BEGIN
          xr           = minmax(WL)
          WL           = WL[500:*]
          cube         = cube[500:*,*,*]
          Jupiter_cube = Jupiter_cube[500:*,*,*]
        end
        order.name eq 'order_53': BEGIN
          xr           = minmax(WL)
          WL           = WL[0:2620]
          cube         = cube[0:2620,*,*]
          Jupiter_cube = Jupiter_cube[0:2620,*,*]
        end
        order.name eq 'order_56': BEGIN
          ;xr           = order.WL_range
          xr           = minmax(WL)
          yr_residual  = [-1.99,19.99]
        end
        order.name eq 'order_60': BEGIN
          xr           = minmax(WL)
          ;xr           = order.WL_range
          yr_residual  = [-3.99,6.99]
        end
        ;strmatch(order.name, 'order_**b'): xr = order.WL_range
        else: begin
          yr_residual  = !null
          xr           = minmax(WL)
        end
      endcase

    ; ========= Fit & extract the sky ===========================
      dummy  = Cube
      s = size(cube, /dim)
      for f = 0, N_elements(Io_frames )-1 do begin
        dummy[*,EXTRACT_ROWS[0,f]:EXTRACT_ROWS[1,f],f] = !values.F_Nan
        Background = total(dummy[*, *, f], 2, /nan)
        Case 1 of
          order.name eq 'order_43': Background = Background / mean(Background[0:2925])
          order.name eq 'order_52': Background = Background / mean(Background[500:*])
          order.name eq 'order_53': Background = Background / mean(Background[0:2620])
          else: Background = Background / mean(Background)
        endcase
        img = reform(dummy[*,*,f])
        for i = 0, s[0]-1 do begin
          img[i,*] = interpol(img[i,*], findgen(s[1]), findgen(s[1]), /Nan) ; Tests showed this is much better than fitting, or using interpol keywords and takes less time
        endfor
        resid = reform(Cube[*,*,f]) - img
        Case 1 of
          order.name eq 'order_43': statsec = resid[0:2925,*]
          order.name eq 'order_52': statsec = resid[500:*,*]
          order.name eq 'order_53': statsec = resid[0:2620,*]
          else: statsec = resid
        endcase
        slit_profile = mean(statsec, dimension = 1, /NAN)
        fake_profile = rebin(transpose(slit_profile), s[0], s[1]) * rebin(Background, s[0], s[1])
        dummy[*,*,f] = img + fake_profile ; fake profile seems to make suprisingly little difference?

      endfor
      Jupiter_extract = fltarr(N_elements(WL), N_elements(Io_frames))
      Io_extract      = fltarr(N_elements(WL), N_elements(Io_frames))
      Sky_Extract     = fltarr(N_elements(WL), N_elements(Io_frames))
      for f = 0, N_elements(Io_frames)-1 do begin
        Io_extract[*,f]         = total(cube[*,EXTRACT_ROWS[0,f]:EXTRACT_ROWS[1,f],f], 2, /nan)         ; Extract rows over Io to 1D
        Sky_extract[*,f]        = total(dummy[*,EXTRACT_ROWS[0,f]:EXTRACT_ROWS[1,f],f], 2, /nan)        ; Extract the Interpolated background beneath Io over the same regio
        Jupiter_extract[*,f]    = total(Jupiter_cube[*,EXTRACT_ROWS[0,f]:EXTRACT_ROWS[1,f],1], 2, /nan) ; Extract the same region of a spectrum of Jupiter at disk center
        ; The two jupiter frames bracketing the observation do not change ~10% in flux.
      endfor

      ;----------------------Flux Calibrate Using Distance and Absolute Spectral Reflectivity of Jupiter---------------------------

      ; Determine the instrumental sensitivity from the expected versus measured flux at Jupiter Disk Center. This should be a smooth function
        expected_flux          = interpol(Rayleighs_per_angstrom, WL_A, WL)       ; move expected flux to the Jovian Doppler-shift, UNITS are R / A
        smoothed_expected_flux = GAUSS_SMOOTH(expected_flux, 3.6, /edge_truncate) ; this smoothing looks about right for HIRES D3

      ; The aperture of HIRES D3 measures about 20 pixels in every order and it is is 7" long
        Io_Area     = !pi * (tan(1821.6 / norm(Jupiter_Earth_State[0:2])) * 206265.)^2    ; Io's solid angle in square arcseconds
        plate_scale = 7. / 20. ; "/pixel
        slit_width  = 1.722    ; " in D3 mode
        for f = 0, N_elements(Io_frames)-1 do begin
          Sensitivity          = Jupiter_extract[*,f] / smoothed_expected_flux            ; Sensitivity in (DN / S) / (R / A)
          ind                  = reverse(sort(Sensitivity))
          ind                  = ind[0:2000]                                              ; pick only the top 50% most sensitive points
          sens_coeffs          = poly_fit(WL[ind], Sensitivity[ind], 2)
          fit_Sensitivity      = poly(WL,sens_coeffs)
          Slit_area_summed     = slit_width * plate_scale * (EXTRACT_ROWS[1,f] - EXTRACT_ROWS[0,f])     ; square arcsec
          Io_extract[*,f]      = Io_extract[*,f] / (fit_Sensitivity * 1.e3)               ; convert to kR / A
          Sky_extract[*,f]     = Sky_extract[*,f] / (fit_Sensitivity * 1.e3)
          Io_extract[*,f]      = Io_extract[*,f] * Slit_area_summed / Io_Area             ; slit filling factor aperture correction ---> ASSUMES EMISSION REGION IS THE SIZE OF IO'S DISK
          Sky_extract[*,f]     = Sky_extract[*,f] * Slit_area_summed / Io_Area            ; slit filling factor aperture correction ---> ASSUMES EMISSION REGION IS THE SIZE OF IO'S DISK
          ;cgplot, WL, Sensitivity, /xstyle, yr = [0.,0.025], Ytitle = 'Measured Flux Sensitivity (DN/S) / (R/A)', Xtitle = 'Angstroms'
          ;cgplot, WL, fit_Sensitivity, color = 'red', /overplot
        endfor
        
      ; Define arrays & plot scaling
        include_WLs            = where( (wl gt xr[0]) and (wl lt xr[1]), /NULL )          ; wavelength regions over which we'll fit the jovian scatter
        plot_WLs               = where( (wl gt xr[0]) and (wl lt xr[1]), /NULL )          ; wavelengths to actually plot
        yr                     = minmax(Io_extract[plot_WLs,*])
        residual_array         = fltarr(N_elements(include_WLs), n_elements(Io_frames ))
        plot_residual_array    = fltarr(N_elements(plot_WLs), n_elements(Io_frames ))
        err_plot_residual_array= fltarr(N_elements(plot_WLs), n_elements(Io_frames ))
        DopplerShift_array     = fltarr(N_elements(Io_airglow), n_elements(Io_frames ))

      ; setup the plot frames
        window, 0
        pos        = cgLayout([1,3], OXMargin=[10,3], OYMargin=[8,5], YGap=0)
        Pos[1,0]   = Pos[1,0]*.8 & Pos[3,1] = Pos[3,1]*.8
        Pos[1,1]   = Pos[1,1]*.84 & Pos[3,2] = Pos[3,2]*.84

      MX_plus_B_parinfo          = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0]}, 2)
      MX_plus_B_parinfo[1].fixed = 1                                                    ; peg the additive "B" component of the MX_Plus_B at zero

      cgPS_Open, filename = Directory + 'Reduced\'+order.name+'_scatterfit.eps', /ENCAPSULATED, xsize = 7.5, ysize = 5
      !P.font = 1
      device, SET_FONT = 'Helvetica Bold', /TT_FONT
      !p.charsize = 1.5

      cgplot, WL, reform(Io_extract[*,0]), psym = 0, Ytitle = 'kR / '+cgsymbol('Angstrom'), $
        title = title,  yr = yr, /xstyle, /ystyle, $
        xr = xr, /nodata, pos = pos[*,0], xtickformat = '(A1)', xminor = 10, xticklen = .025

      Frames = 'j2970'+strcompress(Io_frames , /rem)+'.Cleaned.fits'
      for i = 0, N_elements(Io_frames )-1 do begin
        header = headfits(Directory + 'Reduced\' + Frames[i])

        ; Define fitting indicies we need to avoid because of airglow itself.
          Ios_Airglow    = [Io_airglow + Io_airglow * sxpar(header, 'IO_DOPPL') / cspice_clight()]              ; wavelength locations of all airglow lines
          ARCES_1sigma_bandwidth = 1. * (mean(xr)/31500.) / 2.3548 ; 1. * FWHM in Angstroms / 2sqrt(2ln2) = 68% of flux enclosed within 1 sigma --> Appropriate for weak telluric airglow
          Io_AG_Ind = [] & Telluric_AG_Ind = []
          for j = 0, N_elements(Ios_Airglow)-1 do Io_AG_ind = [Io_AG_Ind, where(abs(Ios_Airglow[j] - WL) lt ARCES_1sigma_bandwidth, /NULL)] ;  indicies within 2 sigma of an Io airglow line
          for j = 0, N_elements(Telluric_Airglow)-1 do Telluric_AG_ind = [Telluric_AG_Ind, where(abs(Telluric_Airglow[j] - WL) lt ARCES_1sigma_bandwidth, /NULL)] ; indicies within 1 sigma of any Telluric airglow
          AG_ind = [Io_AG_Ind, Telluric_AG_Ind]
          if AG_ind ne !null then AG_ind = AG_ind[UNIQ(AG_ind, SORT(AG_ind))]
          if AG_ind ne !null then correl_indicies  = cgSetDifference(include_Wls, AG_Ind) else correl_indicies = include_Wls

        ; Broadband extinction may occcur if Jupiter calibration and Io have significantly different airmass, scale brightness for that
        ; correct for broad-band extinction between the Jupiter airmass and the Io airmass (this is crude, but better than nothing!)
          K = 0.118 ; Mauna Kea extinction coefficient in mag / airmass in greenish redish wavelengths (cf. Krisciunas et al., PASP 99, 887, 1987)
          Scale_for_extinction    = 100. ^ (K * (float(sxpar(header, 'AIRMASS')) - float(sxpar(Jupiter_center_header, 'AIRMASS'))) / 5.)
          print, 'Io Airmass:', sxpar(header, 'AIRMASS'), ' Jupiter Airmass:', sxpar(Jupiter_center_header, 'AIRMASS'), ' Scaling by:', Scale_for_extinction

        spec                      = reform(Io_extract[*,i]) * Scale_for_extinction
        scatter                   = reform(Sky_extract[*,i]) * Scale_for_extinction
        Mult_and_add              = mpfitfun('MX_Plus_B', scatter[correl_indicies], spec[correl_indicies], replicate(stddev(spec[correl_indicies]), N_elements(correl_indicies)), $
                                             [1., 0], /NaN, status = Scatter_fit_status, parinfo = MX_plus_B_parinfo, /quiet)
        scatter_fit               = scatter * Mult_and_add[0] + Mult_and_add[1]
        residual                  = Spec - scatter_fit
        residual_err              = make_array( N_elements(residual), value = stddev(residual[correl_indicies]) )
        plot_residual_array[*, i] = residual[plot_WLs]
        err_plot_residual_array[*, i] = residual_err[plot_WLs]
        DopplerShift_array[*, i]  = Ios_Airglow
        cgplot, WL, scatter_fit, color = timeColors[i], linestyle = 1, thick = 2, /overplot
        cgplot, WL, spec, color = timeColors[i], thick = 2, /overplot
      endfor

      ; Residual plot
      ; first some basic filtering...
        for i = 0, n_elements(Io_frames)-1 do plot_residual_array[*, i] = med_filter(plot_residual_array[*, i], [7, 2.]) ;[7,2.0] are both carefully tested & conservative parameters

      cgtext, 0.03, 0.35, 'Residual (kR / '+cgsymbol('Angstrom')+')', orientation = 90, alignment = 0.5, /normal
      lines_to_fit = where(Ios_airglow gt xr[0] and Ios_airglow lt xr[1], /Null, count_lines)
      if (order.name eq 'order_60') or (order.name eq 'order_56') then count_lines = count_lines else count_lines = 0 ;only plot fits for Na and O 6300/6364A
      if count_lines gt 0 then LSF_Fit_array = fltarr(N_elements(include_WLs), n_elements(Io_frames))
      if not keyword_set(yr_residual) then yr_residual = [-2.*stddev(plot_residual_array[include_WLs,*]), 2.*stddev(plot_residual_array[include_WLs,*])]

      cgplot, spec, WL, psym = 0, xtickformat = '(A1)', $
        yr = yr_residual, xr = xr, /nodata, pos = pos[*,1], /xstyle, /ystyle, /noerase, xminor = 10, xticklen = .05
      for i = 0, n_elements(Io_frames)-1 do begin
        for k = 0, N_elements(Io_Airglow)-1 do begin
          cgplot, replicate(DopplerShift_array[k,i], 2), [.9*yr_residual[1], yr_residual[1]], COLOR = timeColors[i], /overplot;, linestyle = 2
        endfor

        cgplot, WL[plot_WLs], smooth(plot_residual_array[*, i], 1), /OVERPLOT, COLOR = timeColors[i], thick = 1, psym=10

        header         = headfits(Directory + 'Reduced\' + Frames[i])
        Ios_Airglow    = [Io_airglow + Io_airglow * sxpar(header, 'IO_DOPPL') / cspice_clight()]              ; wavelength locations of all airglow lines



        ; fit any airglow lines and save the brightness
          if count_lines gt 0 then begin
            for j = 0, count_lines-1 do begin
  
              parinfo               = replicate({value:0., fixed:0, limited:[0,0], limits:[0.,0]}, 3)
              parinfo[1].fixed      = 1
              parinfo[1].value      = Ios_airglow[lines_to_fit[j]]                                            ; Pin the line's wavelength
              ;parinfo[2].LIMITED    = [1, 1]
              ;parinfo[2].LIMITS     = [0.05,0.15 ]
              parinfo[2].fixed      = 1
              parinfo[2].value      = 0.1 ;APPLIES TO both NA and 6300A
  
              LSF_fitting_ind       = where( abs(wl[plot_WLs] - Ios_airglow[lines_to_fit[j]]) lt 0.3, /NULL)
              fa                    = { x:wl[plot_WLs[LSF_fitting_ind]], y:plot_residual_array[LSF_fitting_ind, i], err: err_plot_residual_array[LSF_fitting_ind, i] }
              a                     = mpfit('Gaussian_for_MPFIT', [1., parinfo[1].value, 0.1], funct=fa, maxiter=50, STATUS = Did_it_work, /Quiet, NPEGGED = NPEGGED, parinfo=parinfo)
  
              LSF_Fit_array[*, i]     = LSF_Fit_array[*, i] + gaussian(WL[include_WLs], a)  ; LSF_Fit
              ; log the results for plotting up later
              Io_Airglow_params.brightness[lines_to_fit[j], i]     = A[0]*A[2]*SQRT(2*!DPI)
              Io_Airglow_params.err_brightness[lines_to_fit[j], i] = stddev(residual[correl_indicies])*A[2]*SQRT(2*!DPI)
              Io_Airglow_params.linewidth[lines_to_fit[j], i]      = A[2]
              Io_Airglow_params.linecenter[lines_to_fit[j], i]     = A[1]
            endfor
            cgplot, WL[plot_WLs], LSF_Fit_array[*, i], /OVERPLOT, COLOR = timeColors[i] ; plot the fit to Io's emission
          endif
      endfor

      ; Co-align to Io's Doppler Shift
      shift_by_array = transpose( DopplerShift_array) - $
        rebin(transpose(mean(DopplerShift_array, dim = 2)), n_elements(Io_frames ), N_elements(Io_Airglow))
      ind = where(io_airglow gt min(WL) and Io_airglow lt max(WL), /NULL, count)
      if count eq 0 then begin
        shifted_WL = fltarr( n_elements(Io_frames ) )
        for i = 0, n_elements(Io_frames )-1 do begin
          header = headfits(Directory + 'Reduced\' + Frames[i])
          shifted_WL[i]  = [mean(WL) + mean(WL)* sxpar(header, 'IO_DOPPL') / cspice_clight()]
        endfor
        shift_by = shifted_WL - mean(shifted_WL)
      endif
      if count eq 1 then shift_by = shift_by_array[*,ind]
      if count gt 1 then shift_by = mean(shift_by_array[*,ind], dim = 2)

      aligned_residuals = plot_residual_array
      for i = 0, n_elements(Io_frames )-1 do begin
        aligned_residuals[*,i] = interpol( plot_residual_array[*,i], WL[plot_WLs], WL[plot_WLs] - Shift_By[i])
      endfor
      co_add = median(aligned_residuals, dimension = 2, /even)

      ; plot the aligned & co-added spectra
        yr_co_add = yr_residual / [5., 2.8]
        cgplot, WL[plot_WLs], co_add, pos = pos[*,2], thick=2, /noerase, xr = xr, yr = yr_co_add, /yst, psym=10, $
          xtitle = 'Wavelength ('+cgsymbol('Angstrom')+')', yminor = 2, /nodata, xminor = 10, xticklen = .05
        cgplot, xr, [0.,0.], linestyle = 2, /overplot

      airglow_x = congrid( mean(DopplerShift_array, dim = 2), N_elements(Io_Airglow)*2+1)
      airglow_x = airglow_x[0:-2]
      airglow_y = Reform(Rebin([-15, 15], 2, N_elements(Io_Airglow)), 2*N_elements(Io_Airglow))
      for i = 0, N_elements(Io_Airglow)-1 do cgplot, airglow_x[2*i:2*i+1], airglow_y[2*i:2*i+1], /overplot, linestyle = 1, color = 'blue', thick = 1.
      cgplot, WL[plot_WLs], smooth(co_add, 5, /edge_mirror), thick = 1, /overplot
      ;cgplot, congrid(WL[plot_WLs], 500), congrid(co_add, 500), thick = 1, /overplot
      cgps_close
    endfor
    save, Io_Airglow_params, filename = Directory + 'Reduced\Io_Airglow_params.sav'

    ;------------- Searching for the vibronic series in Trafton (2012) ----------------------
    ;WN = 17760. + (findgen(25)-8.)*404.5 ; compare with his Figure 1
    ;cgplot, Findgen(18), WN, /ynozero
    ;print, 1.e8*1./WN
    ;----------------------------------------------------------------------------------------
  endif

  ;========================================================================================================================================================
  ;===                                                                                                                                                  ===
  ;===                                                                                                                                                  ===
  ;===                                                        LBT PEPSI Analysis                                                                        ===
  ;===                                                                                                                                                  ===
  ;===                                                                                                                                                  ===
  ;========================================================================================================================================================

  ; Problems:
  ;   -The S-side cross-disperser does not cover 6300A, so that analysis is one side only.
  ;    Stellar calibrations were only obtained on the D side, so no stellar calibrations redward of 6275A ---> Can't correct telluric absorption for O I 6300.
  ;   -The extraction / blaze normalization of the Io eclipse spectra has distinct issues. Orders/blaze show up in the shape of the spectra.
  ;    Since the methodology relies on using Jupiter for absolute flux calibration, and shapw of the Jupiter 1-D extraction differes strongly from Io in eclipse,
  ;    Absolute calibration is very inacccurate. Moreover, orders/blaze structure differs between the S and D eyes.
  ;

  if part eq 10.0 then begin

    ;the Large Binocular Telescope has two eyes: the S side and the D side, Process the D side first

    READCOL,'D:\DATA\___Calibration___\Carl_centrifugal_equator.txt', torus_lat_out, torus_deg, skipline = 1, /Silent

    ;----------------------; load a unity normalized spectrum of HD 145127, a 6.6mag A0V star
    fits             = MRDFITS(Dir+'pepsib.20190424.057.dxt.ffc.all', 1, header, /fscale, /silent, /unsigned )
    header           = headfits(Dir+'pepsib.20190424.057.dxt.ffc.all')
    Telluric_WL_A0V1 = fits.arg
    A0V_norm1        = fits.fun
    A0V1_airmass     = float(sxpar(header, 'AIRMASS'))

    READCOL,'D:\DATA\___Calibration___\Solar_and_Stellar_Calibration_Spectra\vegallpr25.1000', F='A,A', Model_A0V_WL, Model_A0V_Flux, Skipline = 700000, numline = 500000, /Silent
    Model_A0V_WL    = Model_A0V_WL*10.
    Model_A0V_Flux  = float(Model_A0V_Flux)
    Model_A0V_Flux  = float(shift(Model_A0V_Flux, -170))   ; Align the model stellar lines with the measured stellar lins
    Model_A0V_Flux  = INTERPOL(Model_A0V_Flux, Model_A0V_WL, Telluric_WL_A0V1)
    result          = ROBUST_POLY_FIT(Telluric_WL_A0V1, Model_A0V_flux, 6)
    Model_A0V_norm  = Model_A0V_flux / poly(Telluric_WL_A0V1, result)

    ;window, 0, title = 'Unity Normalized A0V Spectral Model (Blue) vs and Measured A0V (Black: HD 155379)'
    ;cgplot, Telluric_WL_A0V1, Model_A0V_norm, color = 'blue', xr = [7650,7700];, xr = [4820,4900]
    ;cgplot, Telluric_WL_A0V1, A0V_norm1, /overplot
    ;cgplot, Telluric_WL_A0V1, gauss_smooth(Model_A0V_norm, 600), /overplot, color = 'green'

    Model_A0V_norm           = gauss_smooth(Model_A0V_norm, 600)
    Telluric_Absorption_A0V1 = A0V_norm1 / Model_A0V_norm

    ; load a unity normalized spectrum of HD 163336, a 5.9 Vmag A0V star, seemingly with much higher spectral resolution --- WHY??
      fits             = MRDFITS(Dir+'pepsib.20190424.071.dxt.ffc.all', 1, header, /fscale, /silent, /unsigned )
      header           = headfits(Dir+'pepsib.20190424.071.dxt.ffc.all')
      Telluric_WL_A0V2 = fits.arg
      A0V_norm2        = fits.fun
      A0V2_airmass     = float(sxpar(header, 'AIRMASS'))

    READCOL,'D:\DATA\___Calibration___\Solar_and_Stellar_Calibration_Spectra\vegallpr25.50000', F='A,A', Model_A0V_WL, Model_A0V_Flux, Skipline = 700000, numline = 500000, /Silent
    Model_A0V_WL    = Model_A0V_WL*10.
    Model_A0V_Flux  = float(Model_A0V_Flux)
    Model_A0V_Flux  = float(shift(Model_A0V_Flux, -180))   ; Align the model stellar lines with the measured stellar lins
    Model_A0V_Flux  = INTERPOL(Model_A0V_Flux, Model_A0V_WL, Telluric_WL_A0V2)
    result          = ROBUST_POLY_FIT(Telluric_WL_A0V2, Model_A0V_flux, 6)
    Model_A0V_norm  = Model_A0V_flux / poly(Telluric_WL_A0V2, result)

    ;window, 1, title = 'Unity Normalized A0V Spectral Model (Blue) vs and Measured A0V (Black: HD 155379)'
    ;cgplot, Telluric_WL_A0V2, Model_A0V_norm, color = 'blue', xr = [7650,7700]
    ;cgplot, Telluric_WL_A0V2, A0V_norm2, /overplot
    ;cgplot, Telluric_WL_A0V2, gauss_smooth(Model_A0V_norm, 30), /overplot, color = 'green'

    Model_A0V_norm           = gauss_smooth(Model_A0V_norm, 30)
    Telluric_Absorption_A0V2 = A0V_norm2 / Model_A0V_norm

    window, 2, title = 'Comparison of Telluric Absorption Spectra:'
    cgplot, Telluric_WL_A0V1, Telluric_Absorption_A0V1, /xst, color = 'blue', xr = [7650,7700]
    cgplot, Telluric_WL_A0V2, Telluric_Absorption_A0V2, /overplot, color = 'red'

    ; remove interstellar Na absorptions, these will otherwise contaiminate the telluric absorption corrections here
    interstellar_Na_ind1 = where(((Telluric_WL_A0V1 gt 5889.07) and (Telluric_WL_A0V1 lt 5889.32)) or $
      ((Telluric_WL_A0V1 gt 5895.05) and (Telluric_WL_A0V1 lt 5895.28)), /NULL)
    interstellar_Na_ind2 = where(((Telluric_WL_A0V2 gt 5888.92) and (Telluric_WL_A0V2 lt 5889.13)) or $
      ((Telluric_WL_A0V2 gt 5894.88) and (Telluric_WL_A0V2 lt 5895.1)), /NULL)
    dummy1 = Telluric_Absorption_A0V1
    dummy2 = Telluric_Absorption_A0V2
    dummy1[interstellar_Na_ind1] = !values.F_NaN
    dummy2[interstellar_Na_ind2] = !values.F_NaN
    Keep_ind1   = where(finite(dummy1), /null)
    Keep_ind2   = where(finite(dummy2), /null)
    No_interstellar1 = INTERPOL( Telluric_Absorption_A0V1[Keep_ind1], Telluric_WL_A0V1[Keep_ind1], Telluric_WL_A0V1)
    No_interstellar2 = INTERPOL( Telluric_Absorption_A0V2[Keep_ind2], Telluric_WL_A0V1[Keep_ind2], Telluric_WL_A0V2)
    ;    cgplot, Telluric_WL_A0V2, Telluric_Absorption_A0V2, /ynozero, xr = [5888,5897]
    ;    cgplot, Telluric_WL_A0V2, No_interstellar2, /overplot, color = 'blue'
    Telluric_Absorption_A0V1 = temporary(No_interstellar1)
    Telluric_Absorption_A0V2 = temporary(No_interstellar2)

    ;----------------------Flux Calibrate Using Distance and Absolute Spectral Reflectivity of Jupiter---------------------------

    ; Compare publications of Jupiter's spectral albedo at disk center

    ; absolute brightness: Digitized Plot from Woodman et al. 1979.
    READCOL,'C:\IDL\Io\Woodman_et_al_1979_plot1.txt', F='A,A', WL, Albedo, STRINGSKIP = '#', /Silent;wavelength in Angstroms I/F unitless
    READCOL,'C:\IDL\Io\Woodman_et_al_1979_plot2_new.txt', F='A,A', WL_2, Albedo_2, STRINGSKIP = '#', /Silent ;wavelength in Angstroms I/F unitless
    Woodman_WL = float([WL, WL_2])              ; STITCH THESE TOGETHER
    Woodman_Albedo = Float([albedo , albedo_2]) ; STITCH THESE TOGETHER

    ; absolute brightness: from Karkoschka (1998) Icarus on the PDS as ID # ESO-J/S/N/U-SPECTROPHOTOMETER-4-V1.0
    READCOL,'C:\IDL\Io\Karkoschka_1995low.tab', F='X,A,X,A', Karkoschka_WL, Karkoschka_Albedo, STRINGSKIP = '#', /Silent ;wavelength in nm I/F unitless
    READCOL,'C:\IDL\Io\Karkoschka_1995high.tab', F='X,A,X,A', Karkoschka_WL_2, Karkoschka_Albedo_2, STRINGSKIP = '#', /Silent ;wavelength in Angstroms I/F unitless
    Karkoschka_WL = float([Karkoschka_WL, Karkoschka_WL_2])             ; STITCH THESE TOGETHER
    Karkoschka_Albedo = Float([Karkoschka_albedo, Karkoschka_albedo_2]) ; STITCH THESE TOGETHER
    Karkoschka_Albedo = Karkoschka_Albedo[sort(Karkoschka_WL)]
    Karkoschka_WL = Karkoschka_WL[sort(Karkoschka_WL)]

    ; compare the two
    ;cgplot, Woodman_WL / 10., Woodman_Albedo, color = 'blue', xstyle = 1., psym = 3, Xtitle = 'Wavelength (nm)', ytitle = 'I/F Reflectivity'
    ;cgplot, Karkoschka_WL, Karkoschka_Albedo*1.35, color = 'red', /overplot
    ;cgtext, 340, .1, 'EQUATOR AT CENTRAL MERIDIAN (Woodman et al. 1979)', color = 'blue'
    ;cgtext, 340, .16, 'FULL DISK scaled by 1.35 (Karkoschka 1998)', color = 'red'

    ; make an informed choice via scaling
    Karkoschka_wl = Karkoschka_wl * 10. ;nm to Angstroms
    Karkoschka_Albedo = Karkoschka_Albedo*1.35 ;Scale the "Full disk albedo" to the "Central Meridian Equatorial Absolute Reflectivity"

    ; Get the Solar Spectral Irradiance at 1 AU
    READCOL,'C:\IDL\Io\Kurucz_2005_irradthuwl.dat', F='A,A', WL_nm, flux, STRINGSKIP = '#', /Silent ;flux is in W/m2/nm
    start  = where(WL_nm eq '299.100')
    WL_nm = float(WL_nm[start:*])
    flux = float(flux[start:*])

    ; change flux units from W/m^2/nm to photons/(cm^2 s A)
    ; multiply by ((lambda / hc) / 1 W)  * (1 m^2 / 1e4 cm^2) * (1 nm / 10 A)
      conversion = ((WL_nm*1.e-9)/(6.62606957e-34*299792458.D)) * (1./1.e4) * (1./10.)
      flux = flux * conversion      ; Cross-checked this result against Huebner et al. 1992.
      WL_A = temporary(WL_nm) * 10. ; Wavelength from nm into angstroms
      VACTOAIR, WL_A, WL_A_Air      ; Vacuum to air wavelength conversion
      WL_A = temporary(WL_A_Air)

    ; Scale the solar flux to Jupiter's instantaneous distance
      Jupiter_center_header = headfits(Dir + Jupiter_Center_File[0])
      cspice_UTC2ET, '2019 24 April ' + sxpar(Jupiter_center_header, 'UT-OBS'), ET
      cspice_spkezr, 'Jupiter', ET, 'J2000', 'LT+S', 'Sun', Jupiter_Sun_State, ltime
      cspice_spkezr, 'Jupiter', ET, 'J2000', 'LT+S', 'Earth', Jupiter_Earth_State, JE_ltime
      solar_distance = norm(Jupiter_Sun_State[0:2]) / 149597871.
      flux_at_jupiter = flux / solar_distance^2.

    ; Multiply incident solar irradiance x spectral albedo to determine the theoreitcal brightness of Jupiter at disk center
      Albedo = INTERPOL(Karkoschka_Albedo, Karkoschka_WL, WL_A)
      Rayleighs_per_angstrom = 4.*flux_at_jupiter*albedo / 1.e6
      if keyword_set(debug) then window, Title = 'Instantaneous Rayleighs per Angstrom: Center of Jupiter''s Disk'
      if keyword_set(debug) then plot, WL_A, Rayleighs_per_angstrom, xr = [5885, 5900], charsize = 2 ;compare to 5.5 MR per angstrom by Brown & Schneider, 1981

    ; Adjust the expected absolute flux for Jupiter's Instantaneous Doppler shift
      theta  = cspice_vsep(Jupiter_Earth_State[0:2], Jupiter_Earth_State[3:5])
      Dopplershift_E = cos(theta) * norm(Jupiter_Earth_State[3:5])   ; scalar projection of the relative velocity along the line of sight
      cspice_spkezr, 'Sun', ET - JE_ltime, 'J2000', 'LT+S', 'Jupiter', Jupiter_Sun_State, ltime
      theta  = cspice_vsep(Jupiter_Sun_State[0:2], Jupiter_Sun_State[3:5])
      Dopplershift_S = cos(theta) * norm(Jupiter_Sun_State[3:5])   ; scalar projection of the relative velocity along the line of sight

    WL_A = WL_A + WL_A*Dopplershift_E/cspice_clight() + WL_A*Dopplershift_S/cspice_clight()

    ; ***********Correct the Jupiter spectrum for Telluric Absorption**********
    fits_2                  = MRDFITS(Dir + Jupiter_Center_File[0], 1, header, /fscale, /silent, /unsigned )
    fits_3                  = MRDFITS(Dir + Jupiter_Center_File[1], 1, header, /fscale, /silent, /unsigned )
    fits_4                  = MRDFITS(Dir + Jupiter_Center_File[2], 1, header, /fscale, /silent, /unsigned )
    fits_6                  = MRDFITS(Dir + Jupiter_Center_File[3], 1, header, /fscale, /silent, /unsigned )
    Jup_center_header       = headfits(Dir + Jupiter_Center_File[0])
    WL                      = [fits_2.arg, fits_3.arg, fits_4.arg, fits_6.arg]   ; concatenate each cross-disperser spectrum
    Jup_center              = [reverse(fits_2.fun), reverse(fits_3.fun), fits_4.fun, fits_6.fun] / (ten(sxpar(Jup_center_header, 'EXPTIME'))*3600.)
    Jup_center_err          = [reverse(fits_2.var), reverse(fits_3.var), fits_4.var, fits_6.var] / (ten(sxpar(Jup_center_header, 'EXPTIME'))*3600.)
    Jup_center_airmass      = float(sxpar(Jup_center_header, 'AIRMASS'))

    ; Make a master telluric correction that combines the different standard stars
      master_telluric_D       = replicate(1., N_elements(jup_center))
      Master_telluric_airmass = replicate(1., N_elements(jup_center))

    rough_absorption_A0V1 = interpol(Telluric_Absorption_A0V1, Telluric_WL_A0V1, WL)
    rough_absorption_A0V2 = interpol(Telluric_Absorption_A0V2, Telluric_WL_A0V2, WL)

    ; Where telluric absorption corrections are needed, carefully identify wavelength indicies where we can isolate the absorption
    for i = 0, N_elements(Io_Airglow)-1 do begin
      telluric_absorption = 0            ; let's default to no correction for Earth's atmosphere, and handle corrections case by case with wavelength
      Case 1 of
        (Io_Airglow[i] eq  Na1) or (Io_Airglow[i] eq  Na2):                                                                ; no telluric absorption
        (Io_Airglow[i] eq  Na3) or (Io_Airglow[i] eq  Na4): begin
          telluric_absorption = 1
          spectral_range = [8173, 8204]  ; Wavelength region for inspection. Unfortunately the standard star lines here begin to saturate (e.g. 8189A), making W's in the corrected spectra
          telluric_fitting_ind = where( (wl gt 8173.) and (wl lt 8180.), /NULL)                                            ; dominated by telluric lines
          ;telluric_fitting_ind = where( (wl gt 8185.) and (wl lt 8193.5), /NULL)                                          ; dominated by telluric lines
          ; Can't fix this without better standard star data
        end
        (Io_Airglow[i] eq  K1): begin
          telluric_absorption = 1
          spectral_range = [7650, 7680]  ; Wavelength region for inspection
          telluric_fitting_ind = where( ((wl gt 7665.) and (wl lt 7667.)) or ((wl gt 7670.) and (wl lt 7673.)), /NULL)     ; dominated by telluric lines
        end
        (Io_Airglow[i] eq  K2): begin
          telluric_absorption = 1
          spectral_range = [7680, 7720]  ; Wavelength region for inspection
          telluric_fitting_ind = where( (wl gt 7695.5) and (wl lt 7697.5), /NULL)                                          ; dominated by telluric lines
          ;telluric_fitting_ind = where( ((wl gt 7659.) and (wl lt 7661.)) or ((wl gt 7670.) and (wl lt 7673.)), /NULL)    ; dominated by telluric lines
        end
        (Io_Airglow[i] eq O1):                                                                                             ; no telluric absorption
        (Io_Airglow[i] eq O2): begin
          telluric_absorption = 0          ; no telluric data here
          spectral_range = [6292, 6308]  ; Wavelength region for inspection
          telluric_fitting_ind = where( ((wl gt 6295.0) and (wl lt 6296.5)) or ((wl gt 6305.5) and (wl lt 6307.0)), /NULL) ; dominated by telluric lines
        end
        (Io_Airglow[i] eq O3):
        (Io_Airglow[i] eq O4) or (Io_Airglow[i] eq  O5) or (Io_Airglow[i] eq O6):                                          ; no telluric absorption
        (Io_Airglow[i] eq O7) or (Io_Airglow[i] eq  O8) or (Io_Airglow[i] eq O9):                                          ; Certainly some tellurics but just too messy to remove them!
        (Io_Airglow[i] eq S1) or (Io_Airglow[i] eq  S2) or (Io_Airglow[i] eq S3):                                          ; Certainly some tellurics but just too messy to remove them!
        else:
      endcase

      if telluric_absorption eq 0 then continue ; otherwise fix the tellurics

      ; Fix any subtle pixel shifts between *known telluric lines* measured in the calibration stars and in Jupiter
        lag        = findgen(21)-10.
        correl_A0V1 = C_CORRELATE( Rough_Absorption_A0V1[telluric_fitting_ind], jup_center[telluric_fitting_ind], lag)
        correl_A0V2 = C_CORRELATE( Rough_Absorption_A0V2[telluric_fitting_ind], jup_center[telluric_fitting_ind], lag)
        yfit_A0V1   = MPFITPEAK(lag, correl_A0V1, A_A0V1, nterms = 3, /positive)
        yfit_A0V2   = MPFITPEAK(lag, correl_A0V2, A_A0V2, nterms = 3, /positive)
        aligned_absorp_A0V1 = interpolate(Rough_Absorption_A0V1, findgen(n_elements(Rough_Absorption_A0V1)) - A_A0V1[1])
        aligned_absorp_A0V2 = interpolate(Rough_Absorption_A0V2, findgen(n_elements(Rough_Absorption_A0V2)) - A_A0V2[1])

      ; Inspect the alignment of both methods of telluric absorption and Jupiter
        window, 0, title = 'Checking Wavelength Alignment of the Telluric Absorption: Bottom plot looks aligned after the above shift? (Blue = BFR, Oreange = A0V) ', ys = 800
        pos = cglayout([1,2])
        cgplot, lag, correl_A0V1, pos = pos[*,0]
        cgplot, lag, yfit_A0V1, /overplot, color = 'blue'
        cgplot, lag, correl_A0V2, /overplot
        cgplot, lag, yfit_A0V2, /overplot, color = 'orange'
        cgplot, WL, jup_center, xr = spectral_range, /ynozero, pos = pos[*,1], /noerase
        cgplot, WL, jup_center / aligned_absorp_A0V1, color = 'blue', /overplot
        cgplot, WL, jup_center / aligned_absorp_A0V2, color = 'orange', /overplot

      print, 'USER---> Carefully inspect the Wavelength Alignment of the Telluric Absorption Correction at Jupiter in every Bandpass:', spectral_range

      ; Fit the telluric absorption: Multiply P[0], Add P[1] and Exponent P[2] to match the match telluric absorption to the Jupiter disk center spectra
        parinfo            = replicate( {value:0., fixed: 0b, limited: [0b,0b], limits: fltarr(2) }, 3)
        parinfo[1].fixed = 1b                               ; limit additive differences
        parinfo[1].value = 0.                               ; additive adjustments should be near zero
        parinfo[2].limited = 1b                             ; limit the ratio of airmass
        parinfo[2].limits  = [.1, 10.]
        p1 = mpfitfun('Match_Telluric_Absorption', aligned_absorp_A0V1[telluric_fitting_ind], jup_center[telluric_fitting_ind], jup_center_err[telluric_fitting_ind], $
                       [0.1, 0.0, A0V1_airmass/jup_center_Airmass], /NaN, status=status1, /quiet, parinfo = parinfo)
        telluric_fit_A0V1 = P1[0]*(aligned_absorp_A0V1^P1[2]) + P1[1]
        p2 = mpfitfun('Match_Telluric_Absorption', aligned_absorp_A0V2[telluric_fitting_ind], jup_center[telluric_fitting_ind], jup_center_err[telluric_fitting_ind], $
                       [0.1, 0.0, A0V2_airmass/jup_center_Airmass], /NaN, status=status2, /quiet, parinfo = parinfo)
        telluric_fit_A0V2 = P2[0]*(aligned_absorp_A0V2^P2[2]) + P2[1]

      Airmass_solution_1 = P1[2] * jup_center_Airmass ; this is the airmass that provides a best fit for a telluric correction
      Airmass_solution_2 = P2[2] * jup_center_Airmass ;

      window, 2, xs = 1000, ys = 900, title='Black = Jupiter, Green = Indices where Telluric Absorptions are scaled, Red = A0V1, Blue = A0V2'
      pos = cgLayout([1,2], OXMargin=[11,3], OYMargin=[9,6], YGap=0)
      sample = where((WL gt spectral_range[0]) and (WL lt spectral_range[1]), /NULL)

      yr = minmax(jup_center[sample])
      cgplot, WL, jup_center, xr = spectral_range, yr = yr, pos = pos[*,0], xtickformat = '(A1)'
      cgplot, WL[telluric_fitting_ind], jup_center[telluric_fitting_ind], /overplot, color = 'green', psym=14
      cgplot, WL, telluric_fit_A0V1, color = 'red', /overplot
      cgplot, WL, telluric_fit_A0V2, color = 'blue', /overplot

      cgplot, WL, jup_center / aligned_absorp_A0V1^(Airmass_solution_1/jup_center_Airmass), color = 'red', xr = spectral_range, /ynozero, /noerase, pos = pos[*,1]
      cgplot, WL, jup_center / aligned_absorp_A0V2^(Airmass_solution_2/jup_center_Airmass), color = 'blue', /overplot

      ; Keep whichever telluric correction method is closer to unity at every spectral bin
      residual_A0V1  = abs(1.-jup_center/telluric_fit_A0V1)
      residual_A0V2  = abs(1.-jup_center/telluric_fit_A0V2)
      junk           = min([[residual_A0V1],[residual_A0V2]], dim = 2, loc)
      combined       = [[jup_center / aligned_absorp_A0V1^(Airmass_solution_1/jup_center_Airmass)], $
                        [jup_center / aligned_absorp_A0V2^(Airmass_solution_2/jup_center_Airmass)]]
      tell_combined  = [[aligned_absorp_A0V1^(Airmass_solution_1/jup_center_Airmass)], $
                        [aligned_absorp_A0V2^(Airmass_solution_2/jup_center_Airmass)]]
      telluric_fit   = tell_combined[loc]
      cgplot, WL, combined[loc], /overplot, color = 'orange'
      cgtext, spectral_range[0]+1, !y.crange[1]*0.80, 'Correction using A0V1', color = 'red'
      cgtext, spectral_range[0]+1, !y.crange[1]*0.75, 'Correction using A0V2', color = 'blue'
      cgtext, spectral_range[0]+1, !y.crange[1]*0.70, 'Correction using a combination', color = 'orange'

      ; write in this region of the telluric absorption
      master_telluric_D[sample] = telluric_fit[sample] ; master telluric will effectively have absorption depths at Jupiter's airmass, verify this!
    endfor

    ; Apply master telluric absorption correction to the Jupiter center spectrum.
    jup_tell_corr = jup_center / master_telluric_D

    ; Determine the instrumental sensitivity from the expected versus measured flux at Jupiter Disk Center. This should be a smooth function
      expected_flux = interpol(Rayleighs_per_angstrom, WL_A, WL)         ; move expected flux to the Jovian Doppler-shift, UNITS are R / A
      smoothed_expected_flux = shift(GAUSS_SMOOTH(expected_flux, 5), 6)  ; this smoothing looks about right for LBT/PEPSI, UNITS are R / A, the shift is a hack
      window, 3
      Sensitivity       = jup_tell_corr / smoothed_expected_flux         ; Sensitivity in (DN / S) / (R / A)
      Sensitivity_Curve = smooth(sensitivity, 5000, /edge_truncate, /nan)
      cgplot, WL, Sensitivity, /xstyle, Ytitle = 'Measured Flux Sensitivity (DN/S) / (R/A)', Xtitle = 'Angstroms'
      fit_Sensitivity   = Sensitivity
      fit_sensitivity[where( (sensitivity gt 0.002) or (sensitivity lt 0.0), /Null)] = !values.F_NaN                                            ; some basic rejection criterion
      Sensitivity_Curve = smooth(fit_sensitivity, 5000, /edge_truncate, /nan)
      fit_sensitivity[where( (Sensitivity_Curve/fit_Sensitivity gt 1.5) or (Sensitivity_Curve/fit_Sensitivity lt 0.5), /Null)] = !values.F_NaN  ; further rejection criterion
      Sensitivity_Curve_D = smooth(fit_sensitivity, 5000, /edge_truncate, /nan)
      cgplot, WL, Sensitivity_Curve_D, color = 'red', /overplot

    ; Inspect and write the telluric corrected R/A Calibrated Jupiter disk center spectra to file
      window, 4, Title = 'Inspect Jupiter: Red = Uncorrected, Black = Telluric Corrected'
      Jup_Cal_Tell_Corr       = jup_tell_corr / Sensitivity_Curve_D   ; Telluric Correct and Convert to Rayleighs per Angstrom
      Jup_Cal_No_Tell_Corr    = jup_center / Sensitivity_Curve_D      ; Convert to Rayleighs per Angstrom, no correction
      Jup_Cal_Tell_Corr_err   = (Jup_center_err / master_telluric_D) / Sensitivity_Curve_D
      cgplot, WL, Jup_Cal_Tell_Corr, xr = spectral_range, /ynozero
      ;cgplot, WL, Jup_Cal_No_Tell_Corr, color = 'red', /overplot

    Jup_Cen = {arg:Wl, fun:Jup_Cal_Tell_Corr, var:Jup_Cal_Tell_Corr_err, no_tell_corr:Jup_Cal_No_Tell_Corr}
    MWRFITS, [], reduced_dir+'R_per_A_Jupiter_Center_D.fits', Jupiter_center_header, /CREATE ; /create overwrites
    MWRFITS, Jup_cen, reduced_dir+'R_per_A_Jupiter_Center_D.fits'
    WL_D = WL ; wavelength array differs on each side

    ;==================================================== Now repeat this same procedure on the S side spectra ===========================================================================

    ; Stellar spectra only exist for the D side so skip that part, start with the S side jupiter files
    for i = 0, N_elements(Jupiter_Center_File)-1 do Jupiter_Center_File[i] = str_replace(Jupiter_Center_File[i], 'dxt', 'sxt')

    ; ***********Correct the Jupiter spectrum for Telluric Absorption**********
    fits_2                  = MRDFITS(Dir + Jupiter_Center_File[0], 1, header, /fscale, /silent, /unsigned )
    fits_3                  = MRDFITS(Dir + Jupiter_Center_File[1], 1, header, /fscale, /silent, /unsigned )
    fits_4                  = MRDFITS(Dir + Jupiter_Center_File[2], 1, header, /fscale, /silent, /unsigned )
    fits_6                  = MRDFITS(Dir + Jupiter_Center_File[3], 1, header, /fscale, /silent, /unsigned )
    Jup_center_header       = headfits(Dir + Jupiter_Center_File[0])
    WL                      = [fits_2.arg, fits_3.arg, fits_4.arg, fits_6.arg]   ; concatenate each cross-disperser spectrum
    Jup_center              = [reverse(fits_2.fun), reverse(fits_3.fun), fits_4.fun, fits_6.fun] / (ten(sxpar(Jup_center_header, 'EXPTIME'))*3600.)
    Jup_center_err          = [reverse(fits_2.var), reverse(fits_3.var), fits_4.var, fits_6.var] / (ten(sxpar(Jup_center_header, 'EXPTIME'))*3600.)
    Jup_center_airmass      = float(sxpar(Jup_center_header, 'AIRMASS'))

    ; Make a master telluric correction that combines the different standard stars
    master_telluric_S       = replicate(1., N_elements(jup_center))
    Master_telluric_airmass = replicate(1., N_elements(jup_center))

    rough_absorption_A0V1 = interpol(Telluric_Absorption_A0V1, Telluric_WL_A0V1, WL)
    rough_absorption_A0V2 = interpol(Telluric_Absorption_A0V2, Telluric_WL_A0V2, WL)

    ; Where telluric absorption corrections are needed, carefully identify wavelength indicies where we can isolate the absorption
    for i = 0, N_elements(Io_Airglow)-1 do begin
      telluric_absorption = 0            ; let's default to no correction for Earth's atmosphere, and handle corrections case by case with wavelength
      Case 1 of
        (Io_Airglow[i] eq  Na1) or (Io_Airglow[i] eq  Na2): begin                                                          ; mild telluric (H2O) absorption
          telluric_absorption = 1
          spectral_range = [5880, 5910]  ; Wavelength region for inspection. Unfortunately the standard star lines here begin to saturate (e.g. 8189A), making W's in the corrected spectra
          telluric_fitting_ind = where( (wl gt 5885.) and (wl lt 5888.), /NULL)                                            ; regions dominated by telluric lines
        end
        (Io_Airglow[i] eq  Na3) or (Io_Airglow[i] eq  Na4): begin
          telluric_absorption = 1
          spectral_range = [8173, 8204]  ; Wavelength region for inspection. Unfortunately the standard star lines here begin to saturate (e.g. 8189A), making W's in the corrected spectra
          telluric_fitting_ind = where( (wl gt 8173.) and (wl lt 8180.), /NULL)                                            ; dominated by telluric lines
          ;telluric_fitting_ind = where( (wl gt 8185.) and (wl lt 8193.5), /NULL)                                          ; dominated by telluric lines
          ; Can't fix this without better standard star data
        end
        (Io_Airglow[i] eq  K1): begin
          telluric_absorption = 1
          spectral_range = [7650, 7680]  ; Wavelength region for inspection
          telluric_fitting_ind = where( ((wl gt 7665.) and (wl lt 7667.)) or ((wl gt 7670.) and (wl lt 7673.)), /NULL)     ; dominated by telluric lines
        end
        (Io_Airglow[i] eq  K2): begin
          telluric_absorption = 1
          spectral_range = [7680, 7720]  ; Wavelength region for inspection
          telluric_fitting_ind = where( (wl gt 7695.5) and (wl lt 7697.5), /NULL)                                          ; dominated by telluric lines
          ;telluric_fitting_ind = where( ((wl gt 7659.) and (wl lt 7661.)) or ((wl gt 7670.) and (wl lt 7673.)), /NULL)    ; dominated by telluric lines
        end
        (Io_Airglow[i] eq O1):                                                                                             ; no telluric absorption
        (Io_Airglow[i] eq O2): begin
          telluric_absorption = 0          ; no telluric data here
          spectral_range = [6292, 6308]  ; Wavelength region for inspection
          telluric_fitting_ind = where( ((wl gt 6295.0) and (wl lt 6296.5)) or ((wl gt 6305.5) and (wl lt 6307.0)), /NULL) ; dominated by telluric lines
        end
        (Io_Airglow[i] eq O3):
        (Io_Airglow[i] eq O4) or (Io_Airglow[i] eq  O5) or (Io_Airglow[i] eq O6):                                          ; no telluric absorption
        (Io_Airglow[i] eq O7) or (Io_Airglow[i] eq  O8) or (Io_Airglow[i] eq O9):                                          ; Certainly some tellurics but just too messy to remove them!
        (Io_Airglow[i] eq S1) or (Io_Airglow[i] eq  S2) or (Io_Airglow[i] eq S3):                                          ; Certainly some tellurics but just too messy to remove them!
        else:
      endcase

      if telluric_absorption eq 0 then continue ; otherwise fix the tellurics

      ; Fix any subtle pixel shifts between *known telluric lines* measured in the calibration stars and in Jupiter
        lag        = findgen(21)-10.
        correl_A0V1 = C_CORRELATE( Rough_Absorption_A0V1[telluric_fitting_ind], jup_center[telluric_fitting_ind], lag)
        correl_A0V2 = C_CORRELATE( Rough_Absorption_A0V2[telluric_fitting_ind], jup_center[telluric_fitting_ind], lag)
        yfit_A0V1   = MPFITPEAK(lag, correl_A0V1, A_A0V1, nterms = 3, /positive)
        yfit_A0V2   = MPFITPEAK(lag, correl_A0V2, A_A0V2, nterms = 3, /positive)
        aligned_absorp_A0V1 = interpolate(Rough_Absorption_A0V1, findgen(n_elements(Rough_Absorption_A0V1)) - A_A0V1[1])
        aligned_absorp_A0V2 = interpolate(Rough_Absorption_A0V2, findgen(n_elements(Rough_Absorption_A0V2)) - A_A0V2[1])

      ; Inspect the alignment of both methods of telluric absorption and Jupiter
        window, 0, title = 'Checking Wavelength Alignment of the Telluric Absorption: Bottom plot looks aligned after the above shift? (Blue = BFR, Oreange = A0V) ', ys = 800
        pos = cglayout([1,2])
        cgplot, lag, correl_A0V1, pos = pos[*,0]
        cgplot, lag, yfit_A0V1, /overplot, color = 'blue'
        cgplot, lag, correl_A0V2, /overplot
        cgplot, lag, yfit_A0V2, /overplot, color = 'orange'
  
        cgplot, WL, jup_center, xr = spectral_range, /ynozero, pos = pos[*,1], /noerase
        cgplot, WL, jup_center / aligned_absorp_A0V1, color = 'blue', /overplot
        cgplot, WL, jup_center / aligned_absorp_A0V2, color = 'orange', /overplot
  
        print, 'USER---> Carefully inspect the Wavelength Alignment of the Telluric Absorption Correction at Jupiter in every Bandpass:', spectral_range
        parinfo            = replicate( {value:0., fixed: 0b, limited: [0b,0b], limits: fltarr(2) }, 3)
        parinfo[1].fixed = 1b                               ; limit additive differences
        parinfo[1].value = 0.                               ; additive adjustments should be near zero
        parinfo[2].limited = 1b                             ; limit the ratio of airmass
        parinfo[2].limits  = [.1, 10.]

      ; Fit the telluric absorption: Multiply P[0], Add P[1] and Exponent P[2] to match the match telluric absorption to the Jupiter disk center spectra
      p1 = mpfitfun('Match_Telluric_Absorption', aligned_absorp_A0V1[telluric_fitting_ind], jup_center[telluric_fitting_ind], jup_center_err[telluric_fitting_ind], $
        [0.1, 0.0, A0V1_airmass/jup_center_Airmass], /NaN, status=status1, /quiet, parinfo = parinfo)
      telluric_fit_A0V1 = P1[0]*(aligned_absorp_A0V1^P1[2]) + P1[1]
      p2 = mpfitfun('Match_Telluric_Absorption', aligned_absorp_A0V2[telluric_fitting_ind], jup_center[telluric_fitting_ind], jup_center_err[telluric_fitting_ind], $
        [0.1, 0.0, A0V2_airmass/jup_center_Airmass], /NaN, status=status2, /quiet, parinfo = parinfo)
      telluric_fit_A0V2 = P2[0]*(aligned_absorp_A0V2^P2[2]) + P2[1]

      Airmass_solution_1 = P1[2] * jup_center_Airmass ; this is the airmass that provides a best fit for a telluric correction
      Airmass_solution_2 = P2[2] * jup_center_Airmass ;

      window, 2, xs = 1000, ys = 900, title='Black = Jupiter, Green = Indices where Telluric Absorptions are scaled, Red = A0V1, Blue = A0V2'
      pos = cgLayout([1,2], OXMargin=[11,3], OYMargin=[9,6], YGap=0)
      sample = where((WL gt spectral_range[0]) and (WL lt spectral_range[1]), /NULL)

      yr = minmax(jup_center[sample])
      cgplot, WL, jup_center, xr = spectral_range, yr = yr, pos = pos[*,0], xtickformat = '(A1)'
      cgplot, WL[telluric_fitting_ind], jup_center[telluric_fitting_ind], /overplot, color = 'green', psym=14
      cgplot, WL, telluric_fit_A0V1, color = 'red', /overplot
      cgplot, WL, telluric_fit_A0V2, color = 'blue', /overplot

      cgplot, WL, jup_center / aligned_absorp_A0V1^(Airmass_solution_1/jup_center_Airmass), color = 'red', xr = spectral_range, /ynozero, /noerase, pos = pos[*,1]
      cgplot, WL, jup_center / aligned_absorp_A0V2^(Airmass_solution_2/jup_center_Airmass), color = 'blue', /overplot

      ; Keep whichever telluric correction method is closer to unity at every spectral bin
        residual_A0V1  = abs(1.-jup_center/telluric_fit_A0V1)
        residual_A0V2  = abs(1.-jup_center/telluric_fit_A0V2)
        junk           = min([[residual_A0V1],[residual_A0V2]], dim = 2, loc)
        combined       = [[jup_center / aligned_absorp_A0V1^(Airmass_solution_1/jup_center_Airmass)], $
                          [jup_center / aligned_absorp_A0V2^(Airmass_solution_2/jup_center_Airmass)]]
        tell_combined  = [[aligned_absorp_A0V1^(Airmass_solution_1/jup_center_Airmass)], $
                          [aligned_absorp_A0V2^(Airmass_solution_2/jup_center_Airmass)]]
        telluric_fit   = tell_combined[loc]

      If ((Io_Airglow[i] eq  Na1) or (Io_Airglow[i] eq  Na2)) then telluric_fit = aligned_absorp_A0V1^(Airmass_solution_1/jup_center_Airmass) ; A0V 2 has major issues near Na-D

      cgplot, WL, combined[loc], /overplot, color = 'orange'
      cgtext, spectral_range[0]+1, !y.crange[1]*0.80, 'Correction using A0V1', color = 'red'
      cgtext, spectral_range[0]+1, !y.crange[1]*0.75, 'Correction using A0V2', color = 'blue'
      cgtext, spectral_range[0]+1, !y.crange[1]*0.70, 'Correction using a combination', color = 'orange'

      ; This is a good stopping place for inspection
      ; STOP STOP  STOP  STOP  STOP  STOP  STOP  STOP  STOP  STOP  STOP  STOP  STOP  STOP

      ; write in this region of the telluric absorption
      master_telluric_S[sample] = telluric_fit[sample] ; master telluric will effectively have absorption depths at Jupiter's airmass, verify this!
    endfor

    ; Apply master telluric absorption correction to the Jupiter center spectrum.
      jup_tell_corr = jup_center / master_telluric_S

    ; Determine the instrumental sensitivity from the expected versus measured flux at Jupiter Disk Center. This should be a smooth function
      expected_flux = interpol(Rayleighs_per_angstrom, WL_A, WL)         ; move expected flux to the Jovian Doppler-shift, UNITS are R / A
      smoothed_expected_flux = shift(GAUSS_SMOOTH(expected_flux, 5), 6)  ; this smoothing looks about right for LBT/PEPSI, UNITS are R / A, the shift is a hack
      window, 3
      Sensitivity       = jup_tell_corr / smoothed_expected_flux         ; Sensitivity in (DN / S) / (R / A)
      Sensitivity_Curve = smooth(sensitivity, 5000, /edge_truncate, /nan)
      cgplot, WL, Sensitivity, /xstyle, Ytitle = 'Measured Flux Sensitivity (DN/S) / (R/A)', Xtitle = 'Angstroms'
      fit_Sensitivity   = Sensitivity
      fit_sensitivity[where( (sensitivity gt 0.002) or (sensitivity lt 0.0), /Null)] = !values.F_NaN                                            ; some basic rejection criterion
      Sensitivity_Curve = smooth(fit_sensitivity, 5000, /edge_truncate, /nan)
      fit_sensitivity[where( (Sensitivity_Curve/fit_Sensitivity gt 1.5) or (Sensitivity_Curve/fit_Sensitivity lt 0.5), /Null)] = !values.F_NaN  ; further rejection criterion
      Sensitivity_Curve_S = smooth(fit_sensitivity, 5000, /edge_truncate, /nan)
      cgplot, WL, Sensitivity_Curve_S, color = 'red', /overplot

    ; Inspect and write the telluric corrected R/A Calibrated Jupiter disk center spectra to file
      window, 4, Title = 'Inspect Jupiter: Red = Uncorrected, Black = Telluric Corrected'
      Jup_Cal_Tell_Corr       = jup_tell_corr / Sensitivity_Curve_S   ; Telluric Correct and Convert to Rayleighs per Angstrom
      Jup_Cal_No_Tell_Corr    = jup_center / Sensitivity_Curve_S      ; Convert to Rayleighs per Angstrom, no correction
      Jup_Cal_Tell_Corr_err   = (Jup_center_err / master_telluric_S) / Sensitivity_Curve_S
      cgplot, WL, Jup_Cal_Tell_Corr, xr = spectral_range, /ynozero
  
      Jup_Cen = {arg:Wl, fun:Jup_Cal_Tell_Corr, var:Jup_Cal_Tell_Corr_err, no_tell_corr:Jup_Cal_No_Tell_Corr}
      MWRFITS, [], reduced_dir+'R_per_A_Jupiter_Center_S.fits', Jupiter_center_header, /CREATE ; /create overwrites
      MWRFITS, Jup_cen, reduced_dir+'R_per_A_Jupiter_Center_S.fits'
      WL_S = WL ; wavelength array differs on each side

    ; **************************************************************** Io Calibration and Absorption Correction ************************************************************
    cspice_UTC2ET, PenUmbra_UTC, PenUmbra_ET
    cspice_UTC2ET, Umbra_UTC, Umbra_ET

    ; the binocular telescope has two "eyes" the s side and d side, process the Io frames from both"
      sunlit_Files_s = sunlit_Files_d
      for i = 0, N_elements(sunlit_Files_d)-1 do sunlit_Files_s[i] = str_replace(sunlit_Files_d[i], 'dxt', 'sxt')

    ; We'll want to process spectra from the sky fibers the same way...
      sky_Files_d = sunlit_Files_d
      sky_Files_s = sunlit_Files_s
      for i = 0, N_elements(sunlit_Files_d)-1 do sky_Files_d[i] = str_replace(sunlit_Files_d[i], 'dxt', 'dxs')
      for i = 0, N_elements(sunlit_Files_s)-1 do sky_Files_s[i] = str_replace(sunlit_Files_s[i], 'sxt', 'sxs')

    ; the binocular telescope has two "eyes" the s side and d side, process the Io frames from both"
      Eclipse_Files_s = Eclipse_Files_d
      for i = 0, N_elements(Eclipse_Files_d)-1 do Eclipse_Files_s[i] = str_replace(Eclipse_Files_d[i], 'dxt', 'sxt')

    ; We'll want to process spectra from the sky fibers the same way...
      sky_Files_d = Eclipse_Files_d
      sky_Files_s = Eclipse_Files_s
      for i = 0, N_elements(Eclipse_Files_d)-1 do sky_Files_d[i] = str_replace(Eclipse_Files_d[i], 'dxt', 'dxs')
      for i = 0, N_elements(Eclipse_Files_s)-1 do sky_Files_s[i] = str_replace(Eclipse_Files_s[i], 'sxt', 'sxs')

    ; ============================================== Sunlit Spectra: Correct tellurics and write the spectra into R/A units ================================================

    for i = 0, n_elements(sunlit_Files_d[0,*])-1 do begin

      fits_b            = MRDFITS(dir+sunlit_Files_d[0,i], 1, header, /fscale, /silent, /unsigned )
      fits_r            = MRDFITS(dir+sunlit_Files_d[1,i], 1, header, /fscale, /silent, /unsigned )
      fits_b_Sky        = MRDFITS(dir+sky_Files_d[0,i], 1, header, /fscale, /silent, /unsigned )
      fits_r_Sky        = MRDFITS(dir+sky_Files_d[1,i], 1, header, /fscale, /silent, /unsigned )
      Io_header         = headfits(dir+sunlit_Files_d[0,i])
      Io_WL             = [fits_b.arg, fits_r.arg]   ; concatenate each cross-disperser spectrum
      Io_sun            = [reverse(fits_b.fun), fits_r.fun] / (ten(sxpar(Io_header, 'EXPTIME'))*3600.)
      Io_sun_err        = [reverse(fits_b.var), fits_r.var] / (ten(sxpar(Io_header, 'EXPTIME'))*3600.)
      Io_Sky_WL         = [fits_b_Sky.arg, fits_r_sky.arg]   ; concatenate each cross-disperser spectrum
      Io_sky            = [reverse(fits_b_Sky.fun), fits_r_Sky.fun] / (ten(sxpar(Io_header, 'EXPTIME'))*3600.)
      Io_sky_err        = [reverse(fits_b_Sky.var), fits_r_Sky.var] / (ten(sxpar(Io_header, 'EXPTIME'))*3600.)

      ; Do the telluric absorption corrections
        telluric_abs_sun         = interpol(master_telluric_D, WL_D, Io_Wl)
        Io_sun_tell_corr         = Io_sun / telluric_abs_sun^( float(sxpar(Io_header, 'airmass')) / float(sxpar(Jup_center_header, 'airmass')) )
        telluric_abs_sky         = interpol(master_telluric_D, WL_D, Io_sky_WL)
        Io_sky_tell_corr         = Io_sky / telluric_abs_sky^( float(sxpar(Io_header, 'airmass')) / float(sxpar(Jup_center_header, 'airmass')) )
        Io_sun_Cal_tell_corr_err = Io_sun_err / telluric_abs_sun^( float(sxpar(Io_header, 'airmass')) / float(sxpar(Jup_center_header, 'airmass')) )
        Io_sky_Cal_tell_corr_err = Io_sky_err / telluric_abs_sky^( float(sxpar(Io_header, 'airmass')) / float(sxpar(Jup_center_header, 'airmass')) )

        K = 0.15 ; Estimate of the broadband extinction coefficient in mag / airmass 
        Scale_for_extinction    = 100. ^ (K * (float(sxpar(Io_header, 'AIRMASS')) - float(sxpar(Jup_center_header, 'AIRMASS'))) / 5.)
        
      ; Convert everything into R / A units      
        Io_sun_Cal_Tell_Corr     = Io_sun_tell_corr * Scale_for_extinction / Sensitivity_Curve_D
        Io_sky_Cal_Tell_Corr     = Io_sky_tell_corr * Scale_for_extinction / Sensitivity_Curve_D
        Io_sun_Cal_tell_corr_err = Io_sun_Cal_tell_corr_err * Scale_for_extinction / Sensitivity_Curve_D
        Io_sky_Cal_tell_corr_err = Io_sky_Cal_tell_corr_err * Scale_for_extinction / Sensitivity_Curve_D
        Io_sun_Cal_No_Tell_Corr  = Io_sun * Scale_for_extinction / Sensitivity_Curve_D
        Io_sky_Cal_No_Tell_Corr  = Io_sky * Scale_for_extinction / Sensitivity_Curve_D

      ; Find the instantaneous Earth-Io Doppler Shift
        cspice_UTC2ET, '2019 24 April ' + sxpar(Io_header, 'UT-OBS'), ET
        ET_mid_exposure = ET + float(sxpar(Io_header, 'EXPTIME'))/2.
        cspice_spkezr, 'Io', ET_mid_exposure, 'J2000', 'LT+S', 'Earth', Io_Earth_State, ltime
        theta  = cspice_vsep(Io_Earth_state[0:2], Io_Earth_state[3:5])
        Io_wrt_Earth_Dopplershift = cos(theta) * norm(Io_Earth_State[3:5])
        SXADDPAR, Io_header, 'Io_DOPPL', Io_wrt_Earth_Dopplershift, 'Io-Earth V_radial in km/s (mid exposure)'
        SXADDPAR, Io_header, 'T_P_Shad', (ET_mid_exposure-PenUmbra_ET) / 60., 'Minutes since Penumbral ingress'
        SXADDPAR, Io_header, 'T_USHADO', (ET_mid_exposure-Umbra_ET) / 60., 'Minutes since Uumbral ingress'

      ; Find the system III longitude and latitude of Io, and write them to the header
        cspice_subpnt, 'Near point: ellipsoid', 'Jupiter', ET_mid_exposure - ltime, 'IAU_Jupiter', 'None', 'Io', Sub_Io, trgepc, srfvec ;get sub solar point in cartesian IAU Jupiter coords
        cspice_bodvrd, 'Jupiter', 'RADII', 3, radii ;get Jupiter shape constants
        re = radii[0]
        rp = radii[2]
        f = (re-rp)/re
        obspos = Sub_Io - srfvec
        cspice_recpgr, 'Jupiter', obspos, re, f, Io_SysIII, Io_SysIII_LATITUDE, opgalt
        ;torus_lat = (2./3.) * 10.31*cos( (196.61 / !radeg) - Io_SysIII )    ; JRM09 Dipole approximation
        torus_lat = interpol(reverse(torus_deg), torus_lat_out, Io_SysIII*!radeg)
        SXADDPAR, Io_header, 'Sys3_Lat', Io_SysIII, 'Io''s System III Latitude'
        SXADDPAR, Io_header, 'Sys3_Lon', Io_SysIII_LATITUDE, 'Io''s System III Longitude'
        SXADDPAR, Io_header, 'Torus_Lat', torus_lat, 'JRM09 Dipole Approx'

      ; Now do the slit filling factor aperture correction. ***************** ASSUMES EMISSION REGION IS THE SIZE OF IO'S DISK *************
        Slit_area = !pi*(2.3/2.)^2                                                                 ; Default ARCES slit size in square arcsec
        Io_Area   = !pi * (tan(1821.6 / norm(Io_Earth_State[0:2])) * 206265.)^2                    ; Io's solid angle in square arcseconds
        Io_sun_Cal_Tell_Corr     = Io_sun_Cal_Tell_Corr * Slit_area / Io_Area
        Io_sky_Cal_Tell_Corr     = Io_sky_Cal_Tell_Corr * Slit_area / Io_Area
        Io_sun_Cal_Tell_Corr_err = Io_sun_Cal_Tell_Corr_err * Slit_area / Io_Area
        Io_sky_Cal_Tell_Corr_err = Io_sky_Cal_Tell_Corr_err * Slit_area / Io_Area
        Io_sun_Cal_No_Tell_Corr  = Io_sun_Cal_No_Tell_Corr * Slit_area / Io_Area
        Io_sky_Cal_No_Tell_Corr  = Io_sky_Cal_No_Tell_Corr * Slit_area / Io_Area

      ; write an output structure formatted after the 1D extraction '.rec' files
        Io_sun = {arg:Io_Wl, fun:Io_sun_Cal_Tell_Corr, var:Io_sun_Cal_Tell_Corr_err, no_tell_corr:Io_sun_Cal_No_Tell_Corr, $
          arg_sky:Io_Sky_WL, fun_sky:Io_sky_Cal_Tell_Corr, var_sky:Io_sky_Cal_Tell_Corr_err, no_tell_corr_sky:Io_sky_Cal_No_Tell_Corr }
        MWRFITS, [], reduced_dir+'R_per_A_Io_Sunlit_'+strcompress(string(i),/remove_all)+'_D.fits', Io_header, /CREATE ; /create overwrites
        MWRFITS, Io_sun, reduced_dir+'R_per_A_Io_Sunlit_'+strcompress(string(i),/remove_all)+'_D.fits'
    endfor
    for i = 0, n_elements(sunlit_Files_s[0,*])-1 do begin ; now the S-side

      fits_b            = MRDFITS(dir+sunlit_Files_s[0,i], 1, header, /fscale, /silent, /unsigned )
      fits_r            = MRDFITS(dir+sunlit_Files_s[1,i], 1, header, /fscale, /silent, /unsigned )
      fits_b_Sky        = MRDFITS(dir+sky_Files_s[0,i], 1, header, /fscale, /silent, /unsigned )
      fits_r_Sky        = MRDFITS(dir+sky_Files_s[1,i], 1, header, /fscale, /silent, /unsigned )
      Io_header         = headfits(dir+sunlit_Files_s[0,i])
      Io_WL             = [fits_b.arg, fits_r.arg]   ; concatenate each cross-disperser spectrum
      Io_sun            = [reverse(fits_b.fun), fits_r.fun] / (ten(sxpar(Io_header, 'EXPTIME'))*3600.)
      Io_sun_err        = [reverse(fits_b.var), fits_r.var] / (ten(sxpar(Io_header, 'EXPTIME'))*3600.)
      Io_Sky_WL         = [fits_b_Sky.arg, fits_r_sky.arg]   ; concatenate each cross-disperser spectrum
      Io_sky            = [reverse(fits_b_Sky.fun), fits_r_Sky.fun] / (ten(sxpar(Io_header, 'EXPTIME'))*3600.)
      Io_sky_err        = [reverse(fits_b_Sky.var), fits_r_Sky.var] / (ten(sxpar(Io_header, 'EXPTIME'))*3600.)

      ; Estimate the
        Io_sun_airmass    = float(sxpar(Io_header, 'AIRMASS'))
        telluric_abs_sun  = interpol(master_telluric_s, WL_s, Io_Wl)
        Io_sun_tell_corr  = Io_sun / telluric_abs_sun^( Io_sun_airmass / Jup_center_airmass )
        telluric_abs_sky  = interpol(master_telluric_s, WL_s, Io_sky_WL)
        Io_sky_tell_corr  = Io_sky / telluric_abs_sky^( Io_sun_airmass / Jup_center_airmass )

        K = 0.15 ; Estimate of the broadband extinction coefficient in mag / airmass
        Scale_for_extinction    = 100. ^ (K * (float(sxpar(Io_header, 'AIRMASS')) - float(sxpar(Jup_center_header, 'AIRMASS'))) / 5.)

      ; do the telluric correction and place things into R/A units
        Io_sun_Cal_Tell_Corr     = Io_sun_tell_corr * Scale_for_extinction / Sensitivity_Curve_s   ; Convert to Rayleighs per Angstrom, Telluric Corrected
        Io_sky_Cal_Tell_Corr     = Io_sky_tell_corr * Scale_for_extinction / Sensitivity_Curve_s   ; Convert to Rayleighs per Angstrom, Telluric Corrected
        Io_sun_Cal_tell_corr_err = (Io_sun_err / telluric_abs_sun^( Io_sun_airmass / Jup_center_airmass ) ) * Scale_for_extinction / Sensitivity_Curve_s
        Io_sky_Cal_tell_corr_err = (Io_sky_err / telluric_abs_sky^( Io_sun_airmass / Jup_center_airmass ) ) * Scale_for_extinction / Sensitivity_Curve_s
        Io_sun_Cal_No_Tell_Corr  = Io_sun * Scale_for_extinction / Sensitivity_Curve_s             ; Convert to Rayleighs per Angstrom, no Telluric correction
        Io_sky_Cal_No_Tell_Corr  = Io_sky * Scale_for_extinction / Sensitivity_Curve_s             ; Convert to Rayleighs per Angstrom, no Telluric correction

      ; find the instantaneous Earth-Io Doppler Shift
        cspice_UTC2ET, '2019 24 April ' + sxpar(Io_header, 'UT-OBS'), ET
        ET_mid_exposure = ET + float(sxpar(Io_header, 'EXPTIME'))/2.
        cspice_spkezr, 'Io', ET_mid_exposure, 'J2000', 'LT+S', 'Earth', Io_Earth_State, ltime
        theta  = cspice_vsep(Io_Earth_state[0:2], Io_Earth_state[3:5])
        Io_wrt_Earth_Dopplershift = cos(theta) * norm(Io_Earth_State[3:5])
        SXADDPAR, Io_header, 'Io_DOPPL', Io_wrt_Earth_Dopplershift, 'Io-Earth V_radial in km/s (mid exposure)'
        SXADDPAR, Io_header, 'T_P_Shad', (ET_mid_exposure-PenUmbra_ET) / 60., 'Minutes since Penumbral ingress'
        SXADDPAR, Io_header, 'T_USHADO', (ET_mid_exposure-Umbra_ET) / 60., 'Minutes since Uumbral ingress'

      ; Find the system III longitude and latitude of Io, and write them to the header
        cspice_subpnt, 'Near point: ellipsoid', 'Jupiter', ET_mid_exposure - ltime, 'IAU_Jupiter', 'None', 'Io', Sub_Io, trgepc, srfvec ;get sub solar point in cartesian IAU Jupiter coords
        cspice_bodvrd, 'Jupiter', 'RADII', 3, radii ;get Jupiter shape constants
        re = radii[0]
        rp = radii[2]
        f = (re-rp)/re
        obspos = Sub_Io - srfvec
        cspice_recpgr, 'Jupiter', obspos, re, f, Io_SysIII, Io_SysIII_LATITUDE, opgalt
        ;torus_lat = (2./3.) * 10.31*cos( (196.61 / !radeg) - Io_SysIII )    ; JRM09 Dipole approximation
        torus_lat = interpol(reverse(torus_deg), torus_lat_out, Io_SysIII*!radeg)
        SXADDPAR, Io_header, 'Sys3_Lat', Io_SysIII, 'Io''s System III Latitude'
        SXADDPAR, Io_header, 'Sys3_Lon', Io_SysIII_LATITUDE, 'Io''s System III Longitude'
        SXADDPAR, Io_header, 'Torus_Lat', torus_lat, 'JRM09 Dipole Approx'

      ; Now do the slit filling factor aperture correction. ***************** ASSUMES EMISSION REGION IS THE SIZE OF IO'S DISK *************
        Slit_area = !pi*(2.3/2.)^2                                                                 ; Default ARCES slit size in square arcsec
        Io_Area   = !pi * (tan(1821.6 / norm(Io_Earth_State[0:2])) * 206265.)^2                    ; Io's solid angle in square arcseconds
        Io_sun_Cal_Tell_Corr     = Io_sun_Cal_Tell_Corr * Slit_area / Io_Area
        Io_sky_Cal_Tell_Corr     = Io_sky_Cal_Tell_Corr * Slit_area / Io_Area
        Io_sun_Cal_Tell_Corr_err = Io_sun_Cal_Tell_Corr_err * Slit_area / Io_Area
        Io_sky_Cal_Tell_Corr_err = Io_sky_Cal_Tell_Corr_err * Slit_area / Io_Area
        Io_sun_Cal_No_Tell_Corr  = Io_sun_Cal_No_Tell_Corr * Slit_area / Io_Area
        Io_sky_Cal_No_Tell_Corr  = Io_sky_Cal_No_Tell_Corr * Slit_area / Io_Area

      ; write an output structure formatted after the 1D extraction '.rec' files
        Io_sun = {arg:Io_Wl, fun:Io_sun_Cal_Tell_Corr, var:Io_sun_Cal_Tell_Corr_err, no_tell_corr:Io_sun_Cal_No_Tell_Corr, $
          arg_sky:Io_Sky_WL, fun_sky:Io_sky_Cal_Tell_Corr, var_sky:Io_sky_Cal_Tell_Corr_err, no_tell_corr_sky:Io_sky_Cal_No_Tell_Corr }
        MWRFITS, [], reduced_dir+'R_per_A_Io_Sunlit_'+strcompress(string(i),/remove_all)+'_S.fits', Io_header, /CREATE ; /create overwrites
        MWRFITS, Io_sun, reduced_dir+'R_per_A_Io_Sunlit_'+strcompress(string(i),/remove_all)+'_S.fits'
    endfor

    ; Carrying around both Eyes all these different wavelength arrays is tedious. Combine eyes and interpolate everything to a "master" wavelength grid
      spec_D     = MRDFITS(reduced_dir + 'R_per_A_Io_sunlit_0_D.fits', 1, junk_header, /Fscale, /silent )
      MASTER_WL  = spec_D.arg
      SD         = Spec_D

    ; Combine the S and D eyes
    frames              = strcompress([0,2], /remove_all)
    for i = 0, n_elements(frames)-1 do begin
      D                   = MRDFITS(reduced_dir + 'R_per_A_Io_sunlit_' + frames[i] + '_D.fits', 1, junk_header, /Fscale, /silent )
      S                   = MRDFITS(reduced_dir + 'R_per_A_Io_sunlit_' + frames[i] + '_S.fits', 1, junk_header, /Fscale, /silent )
      header              = headfits(reduced_dir + 'R_per_A_Io_sunlit_' + frames[i] + '_D.fits')
      SD.fun              = mean([transpose(interpol(S.fun, S.arg, MASTER_WL)), transpose(interpol(D.fun, D.arg, MASTER_WL))], dim=1)
      SD.NO_TELL_CORR     = mean([transpose(interpol(S.NO_TELL_CORR, S.arg, MASTER_WL)), transpose(interpol(D.NO_TELL_CORR, D.arg, MASTER_WL))], dim=1)
      SD.var              = sqrt( interpol(S.var, S.arg, MASTER_WL)^2 + interpol(D.var, D.arg, MASTER_WL)^2 )
      SD.fun_sky          = mean([transpose(interpol(S.fun_sky, S.arg_sky, MASTER_WL)), transpose(interpol(D.fun_sky, D.arg_sky, MASTER_WL))], dim=1)
      SD.NO_TELL_CORR_sky = mean([transpose(interpol(S.NO_TELL_CORR_SKY, S.arg_sky, MASTER_WL)), transpose(interpol(D.NO_TELL_CORR_SKY, D.arg_sky, MASTER_WL))], dim=1)
      SD.var_sky          = sqrt( interpol(S.var_sky, S.arg_sky, MASTER_WL)^2 + interpol(D.var_sky, D.arg_sky, MASTER_WL)^2 )
      SD.arg              = MASTER_WL
      SD.arg_sky          = MASTER_WL
      MWRFITS, [], reduced_dir+'R_per_A_Io_sunlit_'+strcompress(string(frames[i]),/remove_all)+'_SD.fits', Header, /CREATE ; /create overwrites
      MWRFITS, SD, reduced_dir+'R_per_A_Io_sunlit_'+strcompress(string(frames[i]),/remove_all)+'_SD.fits'
    endfor

    ; ============================================== Eclipse Spectra: Correct tellurics and write D-Side the spectra into R/A units ================================================

    for i = 0, n_elements(Eclipse_Files_d[0,*])-1 do begin

      fits_b            = MRDFITS(dir+Eclipse_Files_d[0,i], 1, header, /fscale, /silent, /unsigned )
      fits_r            = MRDFITS(dir+Eclipse_Files_d[1,i], 1, header, /fscale, /silent, /unsigned )
      fits_b_Sky        = MRDFITS(dir+sky_Files_d[0,i], 1, header, /fscale, /silent, /unsigned )
      fits_r_Sky        = MRDFITS(dir+sky_Files_d[1,i], 1, header, /fscale, /silent, /unsigned )
      Io_header         = headfits(dir+Eclipse_Files_d[0,i])
      Io_WL             = [fits_b.arg, fits_r.arg]   ; concatenate each cross-disperser spectrum
      Io_ecl            = [reverse(fits_b.fun), fits_r.fun] / (ten(sxpar(Io_header, 'EXPTIME'))*3600.)
      Io_ecl_err        = [reverse(fits_b.var), fits_r.var] / (ten(sxpar(Io_header, 'EXPTIME'))*3600.)
      Io_Sky_WL         = [fits_b_Sky.arg, fits_r_sky.arg]   ; concatenate each cross-disperser spectrum
      Io_sky            = [reverse(fits_b_Sky.fun), fits_r_Sky.fun] / (ten(sxpar(Io_header, 'EXPTIME'))*3600.)
      Io_sky_err        = [reverse(fits_b_Sky.var), fits_r_Sky.var] / (ten(sxpar(Io_header, 'EXPTIME'))*3600.)

      ; Do the telluric absorption corrections
        Io_ecl_airmass           = float(sxpar(Io_header, 'AIRMASS'))
        telluric_abs_ecl         = interpol(master_telluric_S, WL_S, Io_Wl)
        Io_ecl_tell_corr         = Io_ecl / telluric_abs_ecl^( Io_ecl_airmass / Jup_center_airmass )
        telluric_abs_sky         = interpol(master_telluric_S, WL_S, Io_sky_WL)
        Io_sky_tell_corr         = Io_sky / telluric_abs_sky^( Io_ecl_airmass / Jup_center_airmass )
        Io_ecl_Cal_tell_corr_err = Io_ecl_err / telluric_abs_ecl^( Io_ecl_airmass / Jup_center_airmass )
        Io_sky_Cal_tell_corr_err = Io_sky_err / telluric_abs_sky^( Io_ecl_airmass / Jup_center_airmass )
        
        K = 0.15 ; Estimate of the broadband extinction coefficient in mag / airmass
        Scale_for_extinction    = 100. ^ (K * (float(sxpar(Io_header, 'AIRMASS')) - float(sxpar(Jup_center_header, 'AIRMASS'))) / 5.)
        
      ; Convert everything into R / A units
        Io_ecl_Cal_Tell_Corr     = Io_ecl_tell_corr * Scale_for_extinction / Sensitivity_Curve_S
        Io_sky_Cal_Tell_Corr     = Io_sky_tell_corr * Scale_for_extinction / Sensitivity_Curve_S
        Io_ecl_Cal_tell_corr_err = Io_ecl_Cal_tell_corr_err * Scale_for_extinction / Sensitivity_Curve_S
        Io_sky_Cal_tell_corr_err = Io_sky_Cal_tell_corr_err * Scale_for_extinction / Sensitivity_Curve_S
        Io_ecl_Cal_No_Tell_Corr  = Io_ecl * Scale_for_extinction / Sensitivity_Curve_S
        Io_sky_Cal_No_Tell_Corr  = Io_sky * Scale_for_extinction / Sensitivity_Curve_S

      ; find the instantaneous Earth-Io Doppler Shift
        cspice_UTC2ET, '2019 24 April ' + sxpar(Io_header, 'UT-OBS'), ET
        ET_mid_exposure = ET + float(sxpar(Io_header, 'EXPTIME'))/2.
        cspice_spkezr, 'Io', ET_mid_exposure, 'J2000', 'LT+S', 'Earth', Io_Earth_State, ltime
        theta  = cspice_vsep(Io_Earth_state[0:2], Io_Earth_state[3:5])
        Io_wrt_Earth_Dopplershift = cos(theta) * norm(Io_Earth_State[3:5])
        SXADDPAR, Io_header, 'Io_DOPPL', Io_wrt_Earth_Dopplershift, 'Io-Earth V_radial in km/s (mid exposure)'
        SXADDPAR, Io_header, 'T_P_Shad', (ET_mid_exposure-PenUmbra_ET) / 60., 'Minutes since Penumbral ingress'
        SXADDPAR, Io_header, 'T_USHADO', (ET_mid_exposure-Umbra_ET) / 60., 'Minutes since Uumbral ingress'

      ; Find the system III longitude and latitude of Io, and write them to the header
        cspice_subpnt, 'Near point: ellipsoid', 'Jupiter', ET_mid_exposure - ltime, 'IAU_Jupiter', 'None', 'Io', Sub_Io, trgepc, srfvec ;get sub solar point in cartesian IAU Jupiter coords
        cspice_bodvrd, 'Jupiter', 'RADII', 3, radii ;get Jupiter shape constants
        re = radii[0]
        rp = radii[2]
        f = (re-rp)/re
        obspos = Sub_Io - srfvec
        cspice_recpgr, 'Jupiter', obspos, re, f, Io_SysIII, Io_SysIII_LATITUDE, opgalt
        ;torus_lat = (2./3.) * 10.31*cos( (196.61 / !radeg) - Io_SysIII )    ; JRM09 Dipole approximation
        torus_lat = interpol(reverse(torus_deg), torus_lat_out, Io_SysIII*!radeg)
        SXADDPAR, Io_header, 'Sys3_Lat', Io_SysIII, 'Io''s System III Latitude'
        SXADDPAR, Io_header, 'Sys3_Lon', Io_SysIII_LATITUDE, 'Io''s System III Longitude'
        SXADDPAR, Io_header, 'Torus_Lat', torus_lat, 'JRM09 Dipole Approx'

      ; Now do the slit filling factor aperture correction. ***************** ASSUMES EMISSION REGION IS THE SIZE OF IO'S DISK *************
        Slit_area = !pi*(2.3/2.)^2                                                                 ; Default ARCES slit size in square arcsec
        Io_Area   = !pi * (tan(1821.6 / norm(Io_Earth_State[0:2])) * 206265.)^2                    ; Io's solid angle in square arcseconds
        Io_ecl_Cal_Tell_Corr     = Io_ecl_Cal_Tell_Corr * Slit_area / Io_Area
        Io_sky_Cal_Tell_Corr     = Io_sky_Cal_Tell_Corr * Slit_area / Io_Area
        Io_ecl_Cal_Tell_Corr_err = Io_ecl_Cal_Tell_Corr_err * Slit_area / Io_Area
        Io_sky_Cal_Tell_Corr_err = Io_sky_Cal_Tell_Corr_err * Slit_area / Io_Area
        Io_ecl_Cal_No_Tell_Corr  = Io_ecl_Cal_No_Tell_Corr * Slit_area / Io_Area
        Io_sky_Cal_No_Tell_Corr  = Io_sky_Cal_No_Tell_Corr * Slit_area / Io_Area

      ; write an output structure formatted after the 1D extraction '.rec' files
      Io_ecl = {arg:Io_Wl, fun:Io_ecl_Cal_Tell_Corr, var:Io_ecl_Cal_Tell_Corr_err, no_tell_corr:Io_ecl_Cal_No_Tell_Corr, $
        arg_sky:Io_Sky_WL, fun_sky:Io_sky_Cal_Tell_Corr, var_sky:Io_sky_Cal_Tell_Corr_err, no_tell_corr_sky:Io_sky_Cal_No_Tell_Corr }
      MWRFITS, [], reduced_dir+'R_per_A_Io_Eclipsed_'+strcompress(string(i),/remove_all)+'_D.fits', Io_header, /CREATE ; /create overwrites
      MWRFITS, Io_ecl, reduced_dir+'R_per_A_Io_Eclipsed_'+strcompress(string(i),/remove_all)+'_D.fits'
    endfor

    ; ============================================== Correct tellurics and write the S side Io spectra into R/A units ================================================

    for i = 0, n_elements(Eclipse_Files_s[0,*])-1 do begin

      fits_b            = MRDFITS(dir+Eclipse_Files_s[0,i], 1, header, /fscale, /silent, /unsigned )
      fits_r            = MRDFITS(dir+Eclipse_Files_s[1,i], 1, header, /fscale, /silent, /unsigned )
      fits_b_Sky        = MRDFITS(dir+sky_Files_s[0,i], 1, header, /fscale, /silent, /unsigned )
      fits_r_Sky        = MRDFITS(dir+sky_Files_s[1,i], 1, header, /fscale, /silent, /unsigned )
      Io_header         = headfits(dir+Eclipse_Files_s[0,i])
      Io_WL             = [fits_b.arg, fits_r.arg]   ; concatenate each cross-disperser spectrum
      Io_ecl            = [reverse(fits_b.fun), fits_r.fun] / (ten(sxpar(Io_header, 'EXPTIME'))*3600.)
      Io_ecl_err        = [reverse(fits_b.var), fits_r.var] / (ten(sxpar(Io_header, 'EXPTIME'))*3600.)
      Io_Sky_WL         = [fits_b_Sky.arg, fits_r_sky.arg]   ; concatenate each cross-disperser spectrum
      Io_sky            = [reverse(fits_b_Sky.fun), fits_r_Sky.fun] / (ten(sxpar(Io_header, 'EXPTIME'))*3600.)
      Io_sky_err        = [reverse(fits_b_Sky.var), fits_r_Sky.var] / (ten(sxpar(Io_header, 'EXPTIME'))*3600.)

      Io_ecl_airmass    = float(sxpar(Io_header, 'AIRMASS'))
      telluric_abs_ecl  = interpol(master_telluric_S, WL_S, Io_Wl)
      Io_ecl_tell_corr  = Io_ecl / telluric_abs_ecl^( Io_ecl_airmass / Jup_center_airmass )
      telluric_abs_sky  = interpol(master_telluric_S, WL_S, Io_sky_WL)
      Io_sky_tell_corr  = Io_sky / telluric_abs_sky^( Io_ecl_airmass / Jup_center_airmass )

      K = 0.15 ; Estimate of the broadband extinction coefficient in mag / airmass
      Scale_for_extinction    = 100. ^ (K * (float(sxpar(Io_header, 'AIRMASS')) - float(sxpar(Jup_center_header, 'AIRMASS'))) / 5.)

      Io_ecl_Cal_Tell_Corr     = Io_ecl_tell_corr * Scale_for_extinction / Sensitivity_Curve_s   ; Convert to Rayleighs per Angstrom, Telluric Corrected
      Io_sky_Cal_Tell_Corr     = Io_sky_tell_corr * Scale_for_extinction / Sensitivity_Curve_s   ; Convert to Rayleighs per Angstrom, Telluric Corrected
      Io_ecl_Cal_tell_corr_err = (Io_ecl_err / telluric_abs_ecl^( Io_ecl_airmass / Jup_center_airmass ) ) * Scale_for_extinction / Sensitivity_Curve_s
      Io_sky_Cal_tell_corr_err = (Io_sky_err / telluric_abs_sky^( Io_ecl_airmass / Jup_center_airmass ) ) * Scale_for_extinction / Sensitivity_Curve_s
      Io_ecl_Cal_No_Tell_Corr  = Io_ecl * Scale_for_extinction / Sensitivity_Curve_s             ; Convert to Rayleighs per Angstrom, no Telluric correction
      Io_sky_Cal_No_Tell_Corr  = Io_sky * Scale_for_extinction / Sensitivity_Curve_s             ; Convert to Rayleighs per Angstrom, no Telluric correction

      ; find the instantaneous Earth-Io Doppler Shift
        cspice_UTC2ET, '2019 24 April ' + sxpar(Io_header, 'UT-OBS'), ET
        ET_mid_exposure = ET + float(sxpar(Io_header, 'EXPTIME'))/2.
        cspice_spkezr, 'Io', ET_mid_exposure, 'J2000', 'LT+S', 'Earth', Io_Earth_State, ltime
        theta  = cspice_vsep(Io_Earth_state[0:2], Io_Earth_state[3:5])
        Io_wrt_Earth_Dopplershift = cos(theta) * norm(Io_Earth_State[3:5])
        SXADDPAR, Io_header, 'Io_DOPPL', Io_wrt_Earth_Dopplershift, 'Io-Earth V_radial in km/s (mid exposure)'
        SXADDPAR, Io_header, 'T_P_Shad', (ET_mid_exposure-PenUmbra_ET) / 60., 'Minutes since Penumbral ingress'
        SXADDPAR, Io_header, 'T_USHADO', (ET_mid_exposure-Umbra_ET) / 60., 'Minutes since Uumbral ingress'

      ; Find the system III longitude and latitude of Io, and write them to the header
        cspice_subpnt, 'Near point: ellipsoid', 'Jupiter', ET_mid_exposure - ltime, 'IAU_Jupiter', 'None', 'Io', Sub_Io, trgepc, srfvec ;get sub solar point in cartesian IAU Jupiter coords
        cspice_bodvrd, 'Jupiter', 'RADII', 3, radii ;get Jupiter shape constants
        re = radii[0]
        rp = radii[2]
        f = (re-rp)/re
        obspos = Sub_Io - srfvec
        cspice_recpgr, 'Jupiter', obspos, re, f, Io_SysIII, Io_SysIII_LATITUDE, opgalt
        ;torus_lat = (2./3.) * 10.31*cos( (196.61 / !radeg) - Io_SysIII )    ; JRM09 Dipole approximation
        torus_lat = interpol(reverse(torus_deg), torus_lat_out, Io_SysIII*!radeg)
        SXADDPAR, Io_header, 'Sys3_Lat', Io_SysIII, 'Io''s System III Latitude'
        SXADDPAR, Io_header, 'Sys3_Lon', Io_SysIII_LATITUDE, 'Io''s System III Longitude'
        SXADDPAR, Io_header, 'Torus_Lat', torus_lat, 'JRM09 Dipole Approx'

      ; Now do the slit filling factor aperture correction. ***************** ASSUMES EMISSION REGION IS THE SIZE OF IO'S DISK *************
        Slit_area = !pi*(2.3/2.)^2                                                                 ; Default ARCES slit size in square arcsec
        Io_Area   = !pi * (tan(1821.6 / norm(Io_Earth_State[0:2])) * 206265.)^2                    ; Io's solid angle in square arcseconds
        Io_ecl_Cal_Tell_Corr     = Io_ecl_Cal_Tell_Corr * Slit_area / Io_Area
        Io_sky_Cal_Tell_Corr     = Io_sky_Cal_Tell_Corr * Slit_area / Io_Area
        Io_ecl_Cal_Tell_Corr_err = Io_ecl_Cal_Tell_Corr_err * Slit_area / Io_Area
        Io_sky_Cal_Tell_Corr_err = Io_sky_Cal_Tell_Corr_err * Slit_area / Io_Area
        Io_ecl_Cal_No_Tell_Corr  = Io_ecl_Cal_No_Tell_Corr * Slit_area / Io_Area
        Io_sky_Cal_No_Tell_Corr  = Io_sky_Cal_No_Tell_Corr * Slit_area / Io_Area

      ; write an output structure formatted after the 1D extraction '.rec' files
        Io_ecl = {arg:Io_Wl, fun:Io_ecl_Cal_Tell_Corr, var:Io_ecl_Cal_Tell_Corr_err, no_tell_corr:Io_ecl_Cal_No_Tell_Corr, $
          arg_sky:Io_Sky_WL, fun_sky:Io_sky_Cal_Tell_Corr, var_sky:Io_sky_Cal_Tell_Corr_err, no_tell_corr_sky:Io_sky_Cal_No_Tell_Corr }
        MWRFITS, [], reduced_dir+'R_per_A_Io_Eclipsed_'+strcompress(string(i),/remove_all)+'_S.fits', Io_header, /CREATE ; /create overwrites
        MWRFITS, Io_ecl, reduced_dir+'R_per_A_Io_Eclipsed_'+strcompress(string(i),/remove_all)+'_S.fits'
    endfor

    ; Carrying around both Eyes all these different wavelength arrays is tedious. Combine eyes and interpolate everything to a "master" wavelength grid
      spec_D     = MRDFITS(reduced_dir + 'R_per_A_Io_Eclipsed_0_D.fits', 1, junk_header, /Fscale, /silent )
      MASTER_WL  = spec_D.arg
      SD         = Spec_D

    ; Combine the S and D eyes
      frames              = strcompress([0,2,4,6], /remove_all)
      filename_matrix     = [['R_per_A_Io_Eclipsed_'+ frames +'_S.fits'], ['R_per_A_Io_Eclipsed_'+ frames +'_D.fits']]
      dims                = size(filename_matrix, /dim)

    for i = 0, n_elements(frames)-1 do begin
      D                   = MRDFITS(reduced_dir + 'R_per_A_Io_Eclipsed_' + frames[i] + '_D.fits', 1, junk_header, /Fscale, /silent )
      S                   = MRDFITS(reduced_dir + 'R_per_A_Io_Eclipsed_' + frames[i] + '_S.fits', 1, junk_header, /Fscale, /silent )
      header              = headfits(reduced_dir + 'R_per_A_Io_Eclipsed_' + frames[i] + '_D.fits')
      SD.fun              = mean([transpose(interpol(S.fun, S.arg, MASTER_WL)), transpose(interpol(D.fun, D.arg, MASTER_WL))], dim=1)
      SD.NO_TELL_CORR     = mean([transpose(interpol(S.NO_TELL_CORR, S.arg, MASTER_WL)), transpose(interpol(D.NO_TELL_CORR, D.arg, MASTER_WL))], dim=1)
      SD.var              = sqrt( interpol(S.var, S.arg, MASTER_WL)^2 + interpol(D.var, D.arg, MASTER_WL)^2 )
      SD.fun_sky          = mean([transpose(interpol(S.fun_sky, S.arg_sky, MASTER_WL)), transpose(interpol(D.fun_sky, D.arg_sky, MASTER_WL))], dim=1)

      ;          window, 0
      ;          cgplot, S.arg_sky, S.fun_sky, xr = [5800, 6000], /ynozero
      ;          cgplot, D.arg_sky, D.fun_sky, /overplot, color = 'red'
      ;
      ;          window, 2
      ;          cgplot, S.arg, S.fun, xr = [5800, 6000], /ynozero
      ;          cgplot, D.arg, D.fun, /overplot, color = 'red'
      ;          wait, 5

      SD.NO_TELL_CORR_sky = mean([transpose(interpol(S.NO_TELL_CORR_SKY, S.arg_sky, MASTER_WL)), transpose(interpol(D.NO_TELL_CORR_SKY, D.arg_sky, MASTER_WL))], dim=1)
      SD.var_sky          = sqrt( interpol(S.var_sky, S.arg_sky, MASTER_WL)^2 + interpol(D.var_sky, D.arg_sky, MASTER_WL)^2 )
      SD.arg              = MASTER_WL
      SD.arg_sky          = MASTER_WL
      MWRFITS, [], reduced_dir+'R_per_A_Io_Eclipsed_'+strcompress(string(frames[i]),/remove_all)+'_SD.fits', Header, /CREATE ; /create overwrites
      MWRFITS, SD, reduced_dir+'R_per_A_Io_Eclipsed_'+strcompress(string(frames[i]),/remove_all)+'_SD.fits'
    endfor

    ; And to the same for the odd frames with the other cross-disperser wavelengths
      spec_D     = MRDFITS(reduced_dir + 'R_per_A_Io_Eclipsed_1_D.fits', 1, junk_header, /Fscale, /silent )
      MASTER_WL  = spec_D.arg
      SD         = Spec_D

    ; Combine the S and D eyes
      frames              = strcompress([1,3,5,7], /remove_all)
      filename_matrix     = [['R_per_A_Io_Eclipsed_'+ frames +'_S.fits'], ['R_per_A_Io_Eclipsed_'+ frames +'_D.fits']]
      dims                = size(filename_matrix, /dim)

    for i = 0, n_elements(frames)-1 do begin
      D                   = MRDFITS(reduced_dir + 'R_per_A_Io_Eclipsed_' + frames[i] + '_D.fits', 1, junk_header, /Fscale, /silent )
      S                   = MRDFITS(reduced_dir + 'R_per_A_Io_Eclipsed_' + frames[i] + '_S.fits', 1, junk_header, /Fscale, /silent )
      header              = headfits(reduced_dir + 'R_per_A_Io_Eclipsed_' + frames[i] + '_D.fits')
      SD.fun              = mean([transpose(interpol(S.fun, S.arg, MASTER_WL)), transpose(interpol(D.fun, D.arg, MASTER_WL))], dim=1)
      SD.NO_TELL_CORR     = mean([transpose(interpol(S.NO_TELL_CORR, S.arg, MASTER_WL)), transpose(interpol(D.NO_TELL_CORR, D.arg, MASTER_WL))], dim=1)
      SD.var              = sqrt( interpol(S.var, S.arg, MASTER_WL)^2 + interpol(D.var, D.arg, MASTER_WL)^2 )
      SD.fun_sky          = mean([transpose(interpol(S.fun_sky, S.arg_sky, MASTER_WL)), transpose(interpol(D.fun_sky, D.arg_sky, MASTER_WL))], dim=1)
      SD.NO_TELL_CORR_sky = mean([transpose(interpol(S.NO_TELL_CORR_SKY, S.arg_sky, MASTER_WL)), transpose(interpol(D.NO_TELL_CORR_SKY, D.arg_sky, MASTER_WL))], dim=1)
      SD.var_sky          = sqrt( interpol(S.var_sky, S.arg_sky, MASTER_WL)^2 + interpol(D.var_sky, D.arg_sky, MASTER_WL)^2 )
      SD.arg              = MASTER_WL
      SD.arg_sky          = MASTER_WL
      MWRFITS, [], reduced_dir+'R_per_A_Io_Eclipsed_'+strcompress(string(frames[i]),/remove_all)+'_SD.fits', Header, /CREATE ; /create overwrites
      MWRFITS, SD, reduced_dir+'R_per_A_Io_Eclipsed_'+strcompress(string(frames[i]),/remove_all)+'_SD.fits'
    endfor
  endif

  ;================================================================= LBT 6300A Waterfall Plot =================================================================================
  if part eq 11.0 then begin

    ; setup the plot axis
      spec   = MRDFITS(reduced_dir + 'R_per_A_Io_Eclipsed_0_D.fits', 1, junk_header, /Fscale, /silent )
      header = headfits(reduced_dir+'R_per_A_Io_Eclipsed_0_D.fits')
      WL     = spec.arg
      bandwidth = 2.5                                                                               ; half with of the plot's x axis in Angstroms
      YR       = [20, 150]
      YR_resid = [-10, 35]
      XR       = [O2 + O2 * sxpar(header, 'IO_DOPPL') / cspice_clight() - Bandwidth, O2 + O2 * sxpar(header, 'IO_DOPPL') / cspice_clight() + Bandwidth]

    ; Which files to analyze
      frames              = strcompress([0,2,4,6], /remove_all)
      filename_matrix     = ['R_per_A_Io_Eclipsed_'+ frames +'_D.fits']
      correl_coeff_matrix = fltarr(size(filename_matrix, /dim))

    ; Define arrays
    include_WLs            = where( (wl gt xr[0]) and (wl lt xr[1]), /NULL )
    plot_WLs               = where( (wl gt xr[0]) and (wl lt xr[1]), /NULL )
    residual_array         = fltarr(N_elements(include_WLs), n_elements(frames))
    plot_residual_array    = fltarr(N_elements(plot_WLs), n_elements(frames))
    LSF_Fit_array          = fltarr(N_elements(include_WLs), n_elements(frames))
    T_P_Shadow_array       = fltarr(n_elements(frames))
    T_U_Shadow_array       = fltarr(n_elements(frames))
    EXPTime_array          = fltarr(n_elements(frames))
    DopplerShift_array1    = fltarr(n_elements(frames))
    Torus_lat              = fltarr(n_elements(frames))
    O2_Fit_params          = REPLICATE({params, brightness:0.0, err_brightness:0.0, mult:0.0, smooth:0.0, shift:0.0, background:'', T_P_shadow:0.0, exptime:0.0}, n_elements(frames))
    Torus_params           = REPLICATE({params, brightness:0.0, err_brightness:0.0, mult:0.0, smooth:0.0, shift:0.0, background:'', T_P_shadow:0.0, exptime:0.0}, n_elements(frames))

    ; get color versus ingress time & setup plot positions
    timeColors = BytScl(findgen(11), Min=0, Max=11, Top=11)
    Color      = timeColors[0]
    cgLoadCT, 33, NColors=5, /reverse
    pos        = cgLayout([1,2], OXMargin=[10,3], OYMargin=[8,5], YGap=0)
    Pos[1,0]   = Pos[1,0]*.7 & Pos[3,1] = Pos[3,1]*.7

    cgPS_Open, filename = Reduced_Dir+'LBT_O6300_scatterfit.eps', /ENCAPSULATED, xsize = 7.5, ysize = 5.
    !P.font = 1
    device, SET_FONT = 'Helvetica Bold', /TT_FONT
    !p.charsize = 1.5

    cgplot, spec.arg, spec.fun/1.e3, psym = 0, Ytitle = 'KR / '+cgsymbol('Angstrom'), $             ; Rayleighs to KiloRayleighs
      title = 'Io''s [O I] 6300'+cgsymbol('Angstrom')+' Response to Eclipse', $
      xr = xr, /nodata, pos = pos[*,0], xtickformat = '(A1)', yr = yr

    for i = 0, N_elements(frames) -1 do begin
      spec     = MRDFITS(reduced_dir + 'R_per_A_Io_Eclipsed_' + frames[i] + '_D.fits', 1, junk_header, /Fscale, /silent )
      WL       = spec.arg
      header   = headfits(reduced_dir+'R_per_A_Io_Eclipsed_' + frames[i] + '_D.fits')
      O2_Io_frame = O2 + O2 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      Torus_lat[i] = sxpar(header, 'Torus_la')

      ; Define fitting indicies we need to avoid because of airglow itself.
      Ios_Airglow    = [Io_airglow + Io_airglow * sxpar(header, 'IO_DOPPL') / cspice_clight()]                                      ; wavelength locations of all airglow lines, Ignore telluric K
      PEPSI_1sigma_bandwidth = 1. * (mean(xr)/43000.) / 2.3548 ; 1. * FWHM in Angstroms / 2sqrt(2ln2) = 68% of flux enclosed within 1 sigma --> Appropriate for weak telluric airglow
      Io_AG_Ind = [] & Telluric_AG_Ind = []
      for j = 0, N_elements(Ios_Airglow)-1 do Io_AG_ind = [Io_AG_Ind, where(abs(Ios_Airglow[j] - WL) lt 2.5*PEPSI_1sigma_bandwidth, /NULL)] ;  indicies within 2 sigma of an Io airglow line
      for j = 0, N_elements(Telluric_Airglow)-1 do Telluric_AG_ind = [Telluric_AG_Ind, where(abs(Telluric_Airglow[j] - WL) lt 1.5*PEPSI_1sigma_bandwidth, /NULL)] ; indicies within 1.5 sigma of any Telluric airglow
      AG_ind = [Io_AG_Ind, Telluric_AG_Ind]
      AG_ind = AG_ind[UNIQ(AG_ind, SORT(AG_ind))]

      ; Find which sky spectra is best correlated with the data
      correl_indicies  = cgSetDifference(include_Wls, AG_Ind)
      for n = 0,3 do begin
        trial_sky              = MRDFITS(reduced_dir + 'R_per_A_Io_Eclipsed_' + frames[n] +'_D.fits', 1, junk_header, /Fscale, /silent )
        aligned_trial_sky      = interpol(trial_sky.fun_sky, trial_sky.arg_sky, WL)
        correl_coeff_matrix[n] = CORRELATE(spec.fun[correl_indicies], aligned_trial_sky[correl_indicies])
      endfor
      junk             = max(correl_coeff_matrix, loc)
      sky_spec         = MRDFITS(reduced_dir + filename_matrix[loc], 1, Junk_header, /Fscale, /silent )
      aligned_sky_spec = interpol(sky_spec.fun_sky, sky_spec.arg_sky, WL)
      scale_sky        = median(spec.fun[correl_indicies] / aligned_sky_spec[correl_indicies])
      print, 'Best sky is: ', filename_matrix[loc], junk
      cgplot, WL, spec.fun/1.e3, Color=timeColors[i], thick = 5, /overplot
      cgplot, WL, scale_sky*aligned_sky_spec/1.e3, COLOR = timeColors[i], linestyle = 1, thick = 5, /overplot

      residual                = (spec.fun - scale_sky*aligned_sky_spec) / 1.e3

      ; Now do the line fitting to the residual. First define Line Spread Function (LSF) Gaussian fit parameter info
      parinfo               = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0]}, 3)
      parinfo[1].fixed      = 1
      parinfo[1].value      = double(O2_Io_frame)                                                           ; Pin the line's wavelength
      LSF_fitting_ind1      = cgSetDifference(where( abs(wl - O2_Io_frame) lt 0.3, /NULL), Telluric_AG_Ind) ; fit region within +/- 0.3A of line center, excluding Telluric airglow

      ; Gaussian fit
      initial_guess = [5.D, parinfo[1].value, 0.1]
      fa            = {x:double(WL[LSF_fitting_ind1]), y:double(residual[LSF_fitting_ind1]), err:double( sqrt(abs( spec.fun[LSF_fitting_ind1] / 1.e3 ) ) )}
      a             = abs(mpfit('Gaussian_for_MPFIT', initial_guess, PERROR = err_a, funct=fa, parinfo = parinfo, STATUS = Did_it_work, /Quiet, NPEGGED = NPEGGED)) ;ABS needed? Looks good, but getting negative linewidths?

      print, 'O 6300 FWHM linewidths =' + string(a[2]*2.3548) + 'A or' + string(a[2]/O2*cspice_clight()) + ' km/s'

      ; write the results
      plot_residual_array[*, i] = residual[plot_WLs]
      LSF_Fit_array[*, i]       = gaussian(WL[include_WLs], a)                         ; LSF_Fit
      Torus_params[i]          = {params, A[0]*A[2]*SQRT(2*!DPI), stddev(residual[correl_indicies])*A[2]*SQRT(2*!DPI), $
        0, 0, float(sxpar(header, 'Torus_la')), filename_matrix[loc], float(sxpar(header, 'T_P_Shad')), ten(sxpar(header, 'EXPTIME'))*60.}
      O2_fit_params[i]          = {params, A[0]*A[2]*SQRT(2*!DPI), stddev(residual[correl_indicies])*A[2]*SQRT(2*!DPI), $
        0, 0, 0, filename_matrix[loc], float(sxpar(header, 'T_P_Shad')), ten(sxpar(header, 'EXPTIME'))*60.}

      DopplerShift_array1[i]    = O2_Io_frame
    endfor

    ; Plot the residual and Gaussian fits
    cgplot, WL[plot_WLs], plot_residual_array[*, 0], psym = 0, Ytitle = 'Residual (kR / '+cgsymbol('Angstrom')+')', xtitle = 'Wavelength ('+cgsymbol('Angstrom')+')', $
      yr = YR_resid, xr = XR, /nodata, pos = pos[*,1], /noerase
    cgtext, O2_Io_frame - 0.7, 17, "Io's Doppler Shift", charsize = 1.4, alignment = 0.5
    cgtext, O2 + .5, 17, "Telluric [O I]", charsize = 1.4, alignment = 0.5
    for i = 0, n_elements(frames)-1 do begin
      spec     = MRDFITS(reduced_dir + 'R_per_A_Io_Eclipsed_' + frames[i] +'_D.fits', 1, junk_header, /Fscale, /silent )
      WL       = spec.arg
      cgplot, WL[plot_WLs], plot_residual_array[*, i], /OVERPLOT, COLOR = timeColors[i], psym=10
      cgplot, WL[plot_WLs], LSF_Fit_array[*, i], /OVERPLOT, COLOR = timeColors[i]
      cgplot, [DopplerShift_array1[i], DopplerShift_array1[i]], [YR_resid[1], 31], COLOR = timeColors[i], /overplot
      cgplot, [O2, O2], [YR_resid[1], 3.], linestyle = 1, /overplot
    endfor
    cgps_close
    save, O2_fit_params, filename = Reduced_Dir+'O2_fit_params.sav'
    save, Torus_params, filename = Reduced_Dir+'O2_Torus_params.sav'
  endif

  ;================================================================= LBT 5577A No Detection =================================================================================
  if part eq 11.02 then begin
    ; Which side of the binocular are we plotting up a waterfall for?
    Eye = 'SD'

    ; setup the plot axis
    spec   = MRDFITS(reduced_dir + 'R_per_A_Io_Eclipsed_0_' + Eye + '.fits', 1, junk_header, /Fscale, /silent )
    header = headfits(reduced_dir+'R_per_A_Io_Eclipsed_0_' + Eye + '.fits')
    WL     = spec.arg
    Bandwidth= 2.5
    XR       = [O1 + O1 * sxpar(header, 'IO_DOPPL') / cspice_clight() - Bandwidth, O1 + O1 * sxpar(header, 'IO_DOPPL') / cspice_clight() + Bandwidth]
    YR       = [30, 200]
    YR_resid = [-5, 5]

    ; Which files to analyze
    frames              = strcompress([0,2,4,6], /remove_all)
    filename_matrix     = ['R_per_A_Io_Eclipsed_'+ frames +'_SD.fits']
    dims                = size(filename_matrix, /dim)
    correl_coeff_matrix = fltarr(dims)

    ; Define arrays
    include_WLs            = where( (wl gt xr[0]) and (wl lt xr[1]), /NULL )
    plot_WLs               = where( (wl gt xr[0]) and (wl lt xr[1]), /NULL )
    residual_array         = fltarr(N_elements(include_WLs), n_elements(frames))
    plot_residual_array    = fltarr(N_elements(plot_WLs), n_elements(frames))
    LSF_Fit_array          = fltarr(N_elements(plot_WLs), n_elements(frames))
    T_P_Shadow_array       = fltarr(n_elements(frames))
    T_U_Shadow_array       = fltarr(n_elements(frames))
    EXPTime_array          = fltarr(n_elements(frames))
    DopplerShift_array1    = fltarr(n_elements(frames))
    DopplerShift_array2    = fltarr(n_elements(frames))
    O1_Fit_params      = REPLICATE({params, brightness:0.0, err_brightness:0.0, mult:0.0, smooth:0.0, shift:0.0, background:'', T_P_shadow:0.0, exptime:0.0}, n_elements(frames))

    ; get color versus ingress time & setup plot positions
    timeColors = BytScl(findgen(11), Min=0, Max=11, Top=11)
    Color      = timeColors[0]
    cgLoadCT, 33, NColors=5, /reverse
    pos        = cgLayout([1,2], OXMargin=[10,3], OYMargin=[8,5], YGap=0)
    Pos[1,0]   = Pos[1,0]*.7 & Pos[3,1] = Pos[3,1]*.7

    cgPS_Open, filename = Reduced_Dir+'LBT_O5577_scatterfit_' + Eye + '.eps', /ENCAPSULATED, xsize = 7.5, ysize = 5.
    !P.font = 1
    device, SET_FONT = 'Helvetica Bold', /TT_FONT
    !p.charsize = 1.5

    cgplot, spec.arg, spec.fun/1.e3, psym = 0, Ytitle = 'kR / '+cgsymbol('Angstrom'), $             ; Rayleighs to KiloRayleighs
      title = 'Io''s [O I] 5577'+cgsymbol('Angstrom')+' Response to Eclipse: ' + Date, $
      xr = xr, /nodata, pos = pos[*,0], xtickformat = '(A1)', yr = yr

    for i = 0, N_elements(frames) -1 do begin
      spec     = MRDFITS(reduced_dir + 'R_per_A_Io_Eclipsed_' + frames[i] + '_' + Eye + '.fits', 1, junk_header, /Fscale, /silent )
      WL       = spec.arg
      header   = headfits(reduced_dir+'R_per_A_Io_Eclipsed_' + frames[i] + '_' + Eye + '.fits')
      O1_Io_frame = O1 + O1 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      junk         = min( abs(WL - O1_Io_frame), O1_Io_Ind)

      ; Define fitting indicies we need to avoid because of airglow itself.
      Ios_Airglow    = [Io_airglow + Io_airglow * sxpar(header, 'IO_DOPPL') / cspice_clight()]                                      ; wavelength locations of all airglow lines, Ignore telluric K
      PEPSI_1sigma_bandwidth = 1. * (mean(xr)/43000.) / 2.3548 ; 1. * FWHM in Angstroms / 2sqrt(2ln2) = 68% of flux enclosed within 1 sigma --> Appropriate for weak telluric airglow
      Io_AG_Ind = [] & Telluric_AG_Ind = []
      for j = 0, N_elements(Ios_Airglow)-1 do Io_AG_ind = [Io_AG_Ind, where(abs(Ios_Airglow[j] - WL) lt 2.5*PEPSI_1sigma_bandwidth, /NULL)] ;  indicies within 2 sigma of an Io airglow line
      for j = 0, N_elements(Telluric_Airglow)-1 do Telluric_AG_ind = [Telluric_AG_Ind, where(abs(Telluric_Airglow[j] - WL) lt 1.5*PEPSI_1sigma_bandwidth, /NULL)] ; indicies within 1.5 sigma of any Telluric airglow
      AG_ind = [Io_AG_Ind, Telluric_AG_Ind]
      AG_ind = AG_ind[UNIQ(AG_ind, SORT(AG_ind))]

      ; Find which sky spectra is best correlated with the data
      correl_indicies  = cgSetDifference(include_Wls, AG_Ind)
      for n = 0, dims[0]-1 do begin
        trial_sky              = MRDFITS(reduced_dir + filename_matrix[n], 1, junk_header, /Fscale, /silent )
        correl_coeff_matrix[n] = CORRELATE(spec.fun[correl_indicies], trial_sky.fun_sky[correl_indicies])
      endfor
      junk      = max(correl_coeff_matrix, loc)
      sky_spec  = MRDFITS(reduced_dir + filename_matrix[loc], 1, Junk_header, /Fscale, /silent )
      scale_sky = median(spec.fun[correl_indicies] / sky_spec.fun_sky[correl_indicies])

      Spec = spec.fun / 1.e3
      Jup  = scale_sky*sky_spec.fun_sky/1.e3

      Ftol = 1.e-4 &  Xtol = 1.e-4 & count = 0                                          ; Define initial tolerances for Amoeba
      smooth_and_shift = [4., 0.] & scale  = [20., 5.]                                  ; Define intial guess and simplex for Amoeba
      trial_smooth_and_shift = AMOEBAX(Ftol, xtol, function_name='shift_smooth', SCALE = Scale, P0 = smooth_and_shift, FUNCTION_VALUE = fval)
      smooth_and_shift = trial_smooth_and_shift
      smoothed_shifted = interpolate(gauss_smooth(jup, smooth_and_shift[0]), findgen(n_elements(Jup)) - smooth_and_shift[1])
      print, 'Best sky is: ', filename_matrix[loc], ', Scaling:', scale_sky, ', Correlation:', junk

      cgplot, WL, Spec, Color=timeColors[i], thick = 5, /overplot
      cgplot, WL, smoothed_shifted, COLOR = timeColors[i], linestyle = 1, thick = 5, /overplot
      ;cgplot, WL[correl_indicies], spec[correl_indicies], psym = 4, /overplot
      residual                 = Spec - smoothed_shifted

      ; Now do the line fitting to the residual. First define Line Spread Function (LSF) Gaussian fit parameter info
      parinfo1               = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0]}, 3)
      parinfo1[1].fixed      = 1
      parinfo1[1].value      = double(O1_Io_frame)                                                        ; Pin the line's wavelength
      LSF_fitting_ind1       = cgSetDifference(where( abs(wl - O1_Io_frame) lt 0.3, /NULL), Telluric_AG_Ind) ; fit region within +/- 0.3A of line center, excluding Telluric airglow

      ; Gaussian fit
      initial_guess = [5.D, parinfo1[1].value, .05]
      fa            = {x:double(WL[LSF_fitting_ind1]), y:double(residual[LSF_fitting_ind1]), err:double( sqrt(abs( spec[LSF_fitting_ind1] )) )}
      a             = mpfit('Gaussian_for_MPFIT', initial_guess, PERROR = err_a, funct=fa, parinfo = parinfo1, STATUS = Did_it_work1, /Quiet, NPEGGED = NPEGGED)

      ; Log the results
      print, 'O1 FWHM linewidths =' + string(a[2]*2.3548) + 'A or' + string(a[2]/Na1*cspice_clight()) + ' km/s'
      plot_residual_array[*, i] = residual[plot_WLs]
      LSF_Fit_array[*, i]       = gaussian(WL[plot_WLs], a)                 ; LSF_Fit
      O1_fit_params[i]      = {params, A[0]*A[2]*SQRT(2*!DPI), stddev(residual[correl_indicies])*A[2]*SQRT(2*!DPI), $
        0, 0, 0, filename_matrix[loc], float(sxpar(header, 'T_P_Shad')), ten(sxpar(header, 'EXPTIME'))*60.}
      plot_residual_array[*, i] = residual[plot_WLs]
      DopplerShift_array1[i]    = O1_Io_frame
    endfor

    ; Plot the residual
    shift_By = mean( [transpose( DopplerShift_array1 - Mean(DopplerShift_array1) )], dimension = 1 )
    aligned_residuals = plot_residual_array

    cgplot, WL[plot_WLs], plot_residual_array[*, 0], psym = 0, Ytitle = 'Residual (kR / '+cgsymbol('Angstrom')+')', xtitle = 'Wavelength ('+cgsymbol('Angstrom')+')', $
      yr = YR_resid, xr = XR, /nodata, pos = pos[*,1], /noerase
    for i = 0, n_elements(frames)-1 do begin
      spec     = MRDFITS(reduced_dir + 'R_per_A_Io_Eclipsed_' + frames[i] +'_' + Eye + '.fits', 1, junk_header, /Fscale, /silent )
      WL       = spec.arg
      cgplot, [DopplerShift_array1[i], DopplerShift_array1[i]], [-5, YR_resid[1]], COLOR = timeColors[i], /overplot
      cgplot, WL[plot_WLs], LSF_Fit_array[*, i], /OVERPLOT, COLOR = timeColors[i]
      cgplot, [O1, O1], YR_resid, linestyle = 1, /overplot
      cgplot, WL[plot_WLs], plot_residual_array[*, i], /OVERPLOT, COLOR = timeColors[i];, psym=10
      aligned_residuals[*,i] = interpol( plot_residual_array[*,i], WL[plot_WLs], WL[plot_WLs] - Shift_By[i])
    endfor
    cgplot, WL[plot_WLs], total(aligned_residuals[*,0:2], 2) / n_elements(frames), /overplot
    cgps_close
    save, O1_fit_params, filename = Reduced_Dir+'O1_fit_params.sav'
  endif

  ;================================================================= LBT Na D =================================================================================
  if part eq 11.2 then begin

    ; Which side of the binocular are we plotting up a waterfall for?
    Eye = 'SD'

    ; setup the plot axis
      spec     = MRDFITS(reduced_dir + 'R_per_A_Io_Eclipsed_0_' + Eye + '.fits', 1, junk_header, /Fscale, /silent )
      header   = headfits(reduced_dir+'R_per_A_Io_Eclipsed_0_' + Eye + '.fits')
      WL       = spec.arg
      XR       = [5888.,5897.5]
      YR       = [15, 125]
      YR_resid = [-4.9, 35]
      Fix_linewidths = 0
      linewidths = fltarr(2,4)

    ; Which files to analyze
      frames              = strcompress([0,2,4,6], /remove_all)
      filename_matrix     = ['R_per_A_Io_Eclipsed_'+ frames +'_SD.fits']
      dims                = size(filename_matrix, /dim)
      correl_coeff_matrix = fltarr(dims)

    ; Define arrays
      include_WLs            = where( (abs(wl - Na1+.4) lt 1.1) or (abs(wl - Na2+.4) lt 1.1), /NULL ); wavelength regions over which we'll fit the jovian scatter
      plot_WLs               = where( (wl gt xr[0]) and (wl lt xr[1]), /NULL )
      residual_array         = fltarr(N_elements(include_WLs), n_elements(frames))
      plot_residual_array    = fltarr(N_elements(plot_WLs), n_elements(frames))
      LSF_Fit_array          = fltarr(N_elements(plot_WLs), n_elements(frames))
      T_P_Shadow_array       = fltarr(n_elements(frames))
      T_U_Shadow_array       = fltarr(n_elements(frames))
      EXPTime_array          = fltarr(n_elements(frames))
      DopplerShift_array1    = fltarr(n_elements(frames))
      DopplerShift_array2    = fltarr(n_elements(frames))
      Na1Na2_Fit_params      = REPLICATE({params, brightness:0.0, err_brightness:0.0, mult:0.0, smooth:0.0, shift:0.0, background:'', T_P_shadow:0.0, exptime:0.0}, n_elements(frames))

    ; get color versus ingress time & setup plot positions
      timeColors = BytScl(findgen(11), Min=0, Max=11, Top=11)
      Color      = timeColors[0]
      cgLoadCT, 33, NColors=5, /reverse
      pos        = cgLayout([1,2], OXMargin=[10,3], OYMargin=[8,5], YGap=0)
      Pos[1,0]   = Pos[1,0]*.7 & Pos[3,1] = Pos[3,1]*.7

    cgPS_Open, filename = Reduced_Dir+'LBT_Na-D_scatterfit_' + Eye + '.eps', /ENCAPSULATED, xsize = 7.5, ysize = 5.
    !P.font = 1
    device, SET_FONT = 'Helvetica Bold', /TT_FONT
    !p.charsize = 1.5

    cgplot, spec.arg, spec.fun/1.e3, psym = 0, Ytitle = 'kR / '+cgsymbol('Angstrom'), $             ; Rayleighs to KiloRayleighs
      title = 'Io''s Na D Response to Eclipse: ' + Date, $
      xr = xr, /nodata, pos = pos[*,0], xtickformat = '(A1)', yr = yr

    for i = 0, N_elements(frames) -1 do begin
      spec     = MRDFITS(reduced_dir + 'R_per_A_Io_Eclipsed_' + frames[i] + '_' + Eye + '.fits', 1, junk_header, /Fscale, /silent )
      WL       = spec.arg
      header   = headfits(reduced_dir+'R_per_A_Io_Eclipsed_' + frames[i] + '_' + Eye + '.fits')
      Na1_Io_frame = Na1 + Na1 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      Na2_Io_frame = Na2 + Na2 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      junk         = min( abs(WL - Na1_Io_frame), Na1_Io_Ind)
      junk         = min( abs(WL - Na2_Io_frame), Na2_Io_Ind)

      ; Define fitting indicies we need to avoid because of airglow itself.
      Ios_Airglow    = [Io_airglow + Io_airglow * sxpar(header, 'IO_DOPPL') / cspice_clight()]                                      ; wavelength locations of all airglow lines, Ignore telluric K
      PEPSI_1sigma_bandwidth = 1. * (mean(xr)/43000.) / 2.3548 ; 1. * FWHM in Angstroms / 2sqrt(2ln2) = 68% of flux enclosed within 1 sigma --> Appropriate for weak telluric airglow
      Io_AG_Ind = [] & Telluric_AG_Ind = []
      for j = 0, N_elements(Ios_Airglow)-1 do Io_AG_ind = [Io_AG_Ind, where(abs(Ios_Airglow[j] - WL) lt 2.5*PEPSI_1sigma_bandwidth, /NULL)] ;  indicies within 2 sigma of an Io airglow line
      for j = 0, N_elements(Telluric_Airglow)-1 do Telluric_AG_ind = [Telluric_AG_Ind, where(abs(Telluric_Airglow[j] - WL) lt 1.5*PEPSI_1sigma_bandwidth, /NULL)] ; indicies within 1.5 sigma of any Telluric airglow
      AG_ind = [Io_AG_Ind, Telluric_AG_Ind]
      AG_ind = AG_ind[UNIQ(AG_ind, SORT(AG_ind))]

      ; Find which sky spectra is best correlated with the data
      correl_indicies  = cgSetDifference(include_Wls, AG_Ind)
      for n = 0, dims[0]-1 do begin
        trial_sky              = MRDFITS(reduced_dir + filename_matrix[n], 1, junk_header, /Fscale, /silent )
        correl_coeff_matrix[n] = CORRELATE(spec.fun[correl_indicies], trial_sky.fun_sky[correl_indicies])
      endfor
      junk      = max(correl_coeff_matrix, loc)
      sky_spec  = MRDFITS(reduced_dir + filename_matrix[loc], 1, Junk_header, /Fscale, /silent )
      scale_sky = median(spec.fun[correl_indicies] / sky_spec.fun_sky[correl_indicies])

      Spec = spec.fun / 1.e3
      Jup  = scale_sky*sky_spec.fun_sky/1.e3

      Ftol = 1.e-4 &  Xtol = 1.e-4 & count = 0                                          ; Define initial tolerances for Amoeba
      smooth_and_shift = [4., 0.] & scale  = [20., 5.]                                  ; Define intial guess and simplex for Amoeba
      trial_smooth_and_shift = AMOEBAX(Ftol, xtol, function_name='shift_smooth', SCALE = Scale, P0 = smooth_and_shift, FUNCTION_VALUE = fval)
      smooth_and_shift = trial_smooth_and_shift
      smoothed_shifted = interpolate(gauss_smooth(jup, smooth_and_shift[0]), findgen(n_elements(Jup)) - smooth_and_shift[1])
      print, 'Best sky is: ', filename_matrix[loc], ', Scaling:', scale_sky, ', Correlation:', junk

      cgplot, WL, Spec, Color=timeColors[i], thick = 5, /overplot
      cgplot, WL, smoothed_shifted, COLOR = timeColors[i], linestyle = 1, thick = 5, /overplot
      ;cgplot, WL[correl_indicies], spec[correl_indicies], psym = 4, /overplot
      residual                 = Spec - smoothed_shifted

      ; Fit the Io emissions
      LSF_fitting_ind1       = where( abs(wl - Na1_Io_frame) lt 0.3, /NULL)
      LSF_fitting_ind2       = where( abs(wl - Na2_Io_frame) lt 0.3, /NULL)
      LSF_Fit1               = mpfitpeak(WL[LSF_fitting_ind1], residual[LSF_fitting_ind1], a, STATUS = Line_fit_Status1, /POSITIVE, nterms = 3)
      LSF_Fit2               = mpfitpeak(WL[LSF_fitting_ind2], residual[LSF_fitting_ind2], b, STATUS = Line_fit_Status2, /POSITIVE, nterms = 3)

      ; Now do the line fitting to the residual. First define Line Spread Function (LSF) Gaussian fit parameter info
      parinfo1               = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0]}, 3)
      parinfo1[1].fixed      = 1
      parinfo1[1].value      = double(Na1_Io_frame)

      parinfo2               = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0]}, 3)
      parinfo2[1].fixed      = 1
      parinfo2[1].value      = double(Na2_Io_frame)

      if fix_linewidths then begin         ; some of the fits get hairy near Jupiter so this is best for estimating brightness
        parinfo1[2].fixed    = 1
        parinfo2[2].fixed    = 1
        parinfo2[2].value    = 0.0521963   ; the median of the fit linewidths,
        parinfo2[2].value    = 0.0521963
      endif

      LSF_fitting_ind1       = cgSetDifference(where( abs(wl - Na1_Io_frame) lt 0.3, /NULL), Telluric_AG_Ind) ; fit region within +/- 0.3A of line center, excluding Telluric airglow
      LSF_fitting_ind2       = cgSetDifference(where( abs(wl - Na2_Io_frame) lt 0.3, /NULL), Telluric_AG_Ind) ; fit region within +/- 0.3A of line center, excluding Telluric airglow

      ; Gaussian fit
      initial_guess = [5.D, parinfo1[1].value, 0.0521963]
      fa            = {x:double(WL[LSF_fitting_ind1]), y:double(residual[LSF_fitting_ind1]), err:double( sqrt(abs( spec[LSF_fitting_ind1] )) )}
      a             = mpfit('Gaussian_for_MPFIT', initial_guess, PERROR = err_a, funct=fa, parinfo = parinfo1, STATUS = Did_it_work1, /Quiet, NPEGGED = NPEGGED)

      initial_guess = [5.D, parinfo2[1].value, 0.0521963]
      fb            = {x:double(WL[LSF_fitting_ind2]), y:double(residual[LSF_fitting_ind2]), err:double( sqrt(abs( spec[LSF_fitting_ind2] )) )}
      b             = mpfit('Gaussian_for_MPFIT', initial_guess, PERROR = err_b, funct=fb, parinfo = parinfo2, STATUS = Did_it_work2, /Quiet, NPEGGED = NPEGGED)

      ; Log the results
      if not fix_linewidths then linewidths[*,i] = [a[2], b[2]]
      D2_over_D1     = (A[0]*A[2]*SQRT(2*!DPI)) / (B[0]*B[2]*SQRT(2*!DPI))
      D2_over_D1_err = D2_over_D1 * sqrt( ((stddev(residual[correl_indicies])*A[2]*SQRT(2*!DPI)) / (A[0]*A[2]*SQRT(2*!DPI)))^2 + $
        ((stddev(residual[correl_indicies])*B[2]*SQRT(2*!DPI)) / (B[0]*B[2]*SQRT(2*!DPI)))^2 )
      print, 'Na D2 FWHM linewidths =' + string(a[2]*2.3548) + 'A or' + string(a[2]/Na1*cspice_clight()) + ' km/s'
      print, 'Na D1 FWHM linewidths =' + string(b[2]*2.3548) + 'A or' + string(b[2]/Na2*cspice_clight()) + ' km/s'
      print, 'D2 ='+string(A[0]*A[2]*SQRT(2*!DPI))+', D1 ='+string(B[0]*B[2]*SQRT(2*!DPI))+', D2/D1 line ratio ='+string(D2_over_D1)+'/-'+string(D2_over_D1_err)

      plot_residual_array[*, i] = residual[plot_WLs]
      LSF_Fit_array[*, i]       = gaussian(WL[plot_WLs], a) + gaussian(WL[plot_WLs], b)                ; LSF_Fit
      Na1Na2_fit_params[i]      = {params, A[0]*A[2]*SQRT(2*!DPI) + B[0]*B[2]*SQRT(2*!DPI), stddev(residual[correl_indicies])*A[2]*SQRT(2*!DPI) + stddev(residual[correl_indicies])*B[2]*SQRT(2*!DPI), $
        0, 0, 0, filename_matrix[loc], float(sxpar(header, 'T_P_Shad')), ten(sxpar(header, 'EXPTIME'))*60.}
      plot_residual_array[*, i] = residual[plot_WLs]
      DopplerShift_array1[i]    = Na1_Io_frame
      DopplerShift_array2[i]    = Na2_Io_frame
    endfor

    ; Plot the residual
    shift_By = mean( [transpose( DopplerShift_array1 - Mean(DopplerShift_array1[1:*]) )], dimension = 1 )
    aligned_residuals = plot_residual_array

    cgplot, WL[plot_WLs], plot_residual_array[*, 0], psym = 0, Ytitle = 'Residual (kR / '+cgsymbol('Angstrom')+')', xtitle = 'Wavelength ('+cgsymbol('Angstrom')+')', $
      yr = YR_resid, xr = XR, /nodata, pos = pos[*,1], /noerase
    for i = 0, n_elements(frames)-1 do begin
      spec     = MRDFITS(reduced_dir + 'R_per_A_Io_Eclipsed_' + frames[i] +'_' + Eye + '.fits', 1, junk_header, /Fscale, /silent )
      WL       = spec.arg
      cgplot, [DopplerShift_array1[i], DopplerShift_array1[i]], [33, YR_resid[1]], COLOR = timeColors[i], /overplot
      cgplot, [DopplerShift_array2[i], DopplerShift_array2[i]], [20, YR_resid[1]], COLOR = timeColors[i], /overplot
      cgplot, WL[plot_WLs], LSF_Fit_array[*, i], /OVERPLOT, COLOR = timeColors[i]
      cgplot, [Na1, Na1], YR_resid, linestyle = 1, /overplot
      cgplot, [Na2, Na2], YR_resid, linestyle = 1, /overplot
      cgplot, WL[plot_WLs], plot_residual_array[*, i], /OVERPLOT, COLOR = timeColors[i], psym=10
      aligned_residuals[*,i] = interpol( plot_residual_array[*,i], WL[plot_WLs], WL[plot_WLs] - Shift_By[i])
    endfor
    cgps_close
    save, Na1Na2_fit_params, filename = Reduced_Dir+'Na1Na2_fit_params.sav'
;    
;    
;    
;    
;    
;    
;    
;    ; fit any airglow lines and save the brightness
;    if count_lines gt 0 then begin
;      for j = 0, count_lines-1 do begin
;
;        ; Now do the line fitting to the residual. First define Line Spread Function (LSF) Gaussian fit parameter info
;        parinfo               = replicate({value:0., fixed:0, limited:[0,0], limits:[0.,0]}, 3)
;        parinfo[1].fixed      = 0
;        parinfo[1].value      = Ios_airglow[lines_to_fit[j]]                                            ; Pin the line's wavelength
;        parinfo[2].fixed      = 0
;        parinfo[2].value      = median(Possible_line_width[i,*])
;
;        ; Gaussian fit
;        LSF_fitting_ind       = where( abs(wl[plot_WLs] - Ios_airglow[lines_to_fit[j]]) lt 0.25, /NULL)
;        tsum_integral_ind     = where( (wl[plot_WLs] - Ios_airglow[lines_to_fit[j]] lt 0.4) and (wl[plot_WLs] - Ios_airglow[lines_to_fit[j]] gt -0.3), /NULL) ;HACKED FIX FOR FAULTY GAUSSIAN ONLY FOR DEMONSTRATION OF 201001
;        fa                    = { x:wl[plot_WLs[LSF_fitting_ind]], y:plot_residual_array[LSF_fitting_ind, i], err: err_plot_residual_array[LSF_fitting_ind, i] }
;        a                     = mpfit('Gaussian_for_MPFIT', [2., parinfo[1].value, parinfo[2].value], funct=fa, maxiter=50, STATUS = Did_it_work, /Quiet, NPEGGED = NPEGGED, parinfo=parinfo)
;        LSF_Fit_array[*, i]   = LSF_Fit_array[*, i] + gaussian(WL[plot_WLs], a)  ; LSF_Fit
;        print, a
;
;        ; log the results for plotting up later
;        Io_Airglow_params.brightness[lines_to_fit[j], i]     = A[0]*abs(A[2])*SQRT(2*!DPI)
;        ;Io_Airglow_params.brightness[lines_to_fit[j], i]     = tsum(wl[plot_WLs[tsum_integral_ind]], plot_residual_array[tsum_integral_ind, i])
;        Io_Airglow_params.err_brightness[lines_to_fit[j], i] = stddev(residual[correl_indicies])*A[2]*SQRT(2*!DPI)
;        Io_Airglow_params.linewidth[lines_to_fit[j], i]      = A[2]
;        Io_Airglow_params.linecenter[lines_to_fit[j], i]     = A[1]
;        Io_Airglow_params.torus_latitude[i]                  = float(sxpar(header, 'Torus_la'))
;        Io_Airglow_params.exptime[i]                         = float(sxpar(fullspec_header, 'EXPTIME')) / 60.
;        Io_Airglow_params.T_P_Shadow[i]                      = sxpar(header, 'T_P_Shad')
;
;      endfor

;    
  endif

  ;================================================================= LBT K-D Waterfall Plot =================================================================================
  if part eq 11.3 then begin

    ; setup the plot axis
    spec   = MRDFITS(reduced_dir + 'R_per_A_Io_Eclipsed_1_S.fits', 1, junk_header, /Fscale, /silent )
    header = headfits(reduced_dir+'R_per_A_Io_Eclipsed_1_S.fits')
    WL     = spec.arg
    ;stop
    bandwidth = 12.5                                                                               ; half with of the plot's x axis in Angstroms
    YR       = [20, 160]
    YR_resid = [-10, 35]
    ;XR       = [7991.80, 8138.06]
    XR       = [8440,9300]


    ; Which files to analyze
    ;frames              = strcompress([1,3,5,7], /remove_all)
    frames              = strcompress([1,3,5], /remove_all)
    Eclipse_files = frames
    filename_matrix     = ['R_per_A_Io_Eclipsed_'+ frames +'_S.fits']
    correl_coeff_matrix = fltarr(size(filename_matrix, /dim))

    ; Define arrays
    include_WLs            = where( (wl gt xr[0]) and (wl lt xr[1]), /NULL )
    plot_WLs               = where( (wl gt xr[0]) and (wl lt xr[1]), /NULL )
    residual_array         = fltarr(N_elements(include_WLs), n_elements(frames))
    plot_residual_array    = fltarr(N_elements(plot_WLs), n_elements(frames))
    LSF_Fit_array          = fltarr(N_elements(include_WLs), n_elements(frames))
    T_P_Shadow_array       = fltarr(n_elements(frames))
    T_U_Shadow_array       = fltarr(n_elements(frames))
    EXPTime_array          = fltarr(n_elements(frames))
    DopplerShift_array     = fltarr(N_elements(Io_airglow), n_elements(frames))
    Torus_lat              = fltarr(n_elements(frames))

    K1_Fit_params          = REPLICATE({params, brightness:0.0, err_brightness:0.0, mult:0.0, smooth:0.0, shift:0.0, background:'', T_P_shadow:0.0, exptime:0.0}, n_elements(frames))
    Torus_params          = REPLICATE({params, brightness:0.0, err_brightness:0.0, mult:0.0, smooth:0.0, shift:0.0, background:'', T_P_shadow:0.0, exptime:0.0}, n_elements(frames))

    ; get color versus ingress time & setup plot positions
    timeColors = BytScl(findgen(11), Min=0, Max=11, Top=11)
    Color      = timeColors[0]
    cgLoadCT, 33, NColors=5, /reverse
    pos        = cgLayout([1,3], OXMargin=[10,3], OYMargin=[8,5], YGap=0)
    Pos[1,0]   = Pos[1,0]*.8 & Pos[3,1] = Pos[3,1]*.8
    Pos[1,1]   = Pos[1,1]*.84 & Pos[3,2] = Pos[3,2]*.84

    cgPS_Open, filename = Reduced_Dir+'LBT_K_scatterfit.eps', /ENCAPSULATED, xsize = 7.5, ysize = 5.
    !P.font = 1
    device, SET_FONT = 'Helvetica Bold', /TT_FONT
    !p.charsize = 1.5

    cgplot, spec.arg, spec.fun/1.e3, psym = 0, Ytitle = 'kR / '+cgsymbol('Angstrom'), $             ; Rayleighs to KiloRayleighs
      title = 'Io''s Airglow in Eclipse', $
      xr = xr, /nodata, pos = pos[*,0], xtickformat = '(A1)', yr = yr

    for i = 0, N_elements(frames) -1 do begin
      spec     = MRDFITS(reduced_dir + 'R_per_A_Io_Eclipsed_' + frames[i] + '_D.fits', 1, junk_header, /Fscale, /silent )
      WL       = spec.arg
      header   = headfits(reduced_dir+'R_per_A_Io_Eclipsed_' + frames[i] + '_D.fits')
      K1_Io_frame = K1 + K1 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      Torus_lat[i] = sxpar(header, 'Torus_la')

      ; Define fitting indicies we need to avoid because of airglow itself.
      Ios_Airglow    = [Io_airglow + Io_airglow * sxpar(header, 'IO_DOPPL') / cspice_clight()]                                      ; wavelength locations of all airglow lines, Ignore telluric K
      PEPSI_1sigma_bandwidth = 1. * (mean(xr)/43000.) / 2.3548 ; 1. * FWHM in Angstroms / 2sqrt(2ln2) = 68% of flux enclosed within 1 sigma --> Appropriate for weak telluric airglow
      Io_AG_Ind = [] & Telluric_AG_Ind = []
      for j = 0, N_elements(Ios_Airglow)-1 do Io_AG_ind = [Io_AG_Ind, where(abs(Ios_Airglow[j] - WL) lt 2.5*PEPSI_1sigma_bandwidth, /NULL)] ;  indicies within 2 sigma of an Io airglow line
      for j = 0, N_elements(Telluric_Airglow)-1 do Telluric_AG_ind = [Telluric_AG_Ind, where(abs(Telluric_Airglow[j] - WL) lt 1.5*PEPSI_1sigma_bandwidth, /NULL)] ; indicies within 1.5 sigma of any Telluric airglow
      AG_ind = [Io_AG_Ind, Telluric_AG_Ind]
      AG_ind = AG_ind[UNIQ(AG_ind, SORT(AG_ind))]

      ; Find which sky spectra is best correlated with the data
      correl_indicies  = cgSetDifference(include_Wls, AG_Ind)
      for n = 0, N_elements(frames)-1 do begin
        trial_sky              = MRDFITS(reduced_dir + 'R_per_A_Io_Eclipsed_' + frames[n] +'_D.fits', 1, junk_header, /Fscale, /silent )
        aligned_trial_sky      = interpol(trial_sky.fun_sky, trial_sky.arg_sky, WL)
        correl_coeff_matrix[n] = CORRELATE(spec.fun[correl_indicies], aligned_trial_sky[correl_indicies])
      endfor
      junk             = max(correl_coeff_matrix, loc)
      sky_spec         = MRDFITS(reduced_dir + filename_matrix[loc], 1, Junk_header, /Fscale, /silent )
      aligned_sky_spec = interpol(sky_spec.fun_sky, sky_spec.arg_sky, WL)
      scale_sky        = median(spec.fun[correl_indicies] / aligned_sky_spec[correl_indicies])
      print, 'Best sky is: ', filename_matrix[loc], junk
      cgplot, WL, spec.fun/1.e3, Color=timeColors[i], thick = 2, /overplot
      cgplot, WL, scale_sky*aligned_sky_spec/1.e3, COLOR = timeColors[i], linestyle = 1, thick = 2, /overplot
      residual                = (spec.fun - scale_sky*aligned_sky_spec) / 1.e3

      ; Now do the line fitting to the residual. First define Line Spread Function (LSF) Gaussian fit parameter info
      parinfo               = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0]}, 3)
      parinfo[1].fixed      = 1
      parinfo[1].value      = double(K1_Io_frame)                                                           ; Pin the line's wavelength
      LSF_fitting_ind1      = cgSetDifference(where( abs(wl - K1_Io_frame) lt 0.3, /NULL), Telluric_AG_Ind) ; fit region within +/- 0.3A of line center, excluding Telluric airglow

      ; Gaussian fit
      initial_guess = [5.D, parinfo[1].value, 0.1]
      fa            = {x:double(WL[LSF_fitting_ind1]), y:double(residual[LSF_fitting_ind1]), err:double( sqrt(abs( spec.fun[LSF_fitting_ind1] / 1.e3 ) ) )}
      a             = abs(mpfit('Gaussian_for_MPFIT', initial_guess, PERROR = err_a, funct=fa, parinfo = parinfo, STATUS = Did_it_work, /Quiet, NPEGGED = NPEGGED)) ;ABS needed? Looks good, but getting negative linewidths?

      print, 'K D1 FWHM linewidths =' + string(a[2]*2.3548) + 'A or' + string(a[2]/K1*cspice_clight()) + ' km/s'

      ; write the results
      plot_residual_array[*, i] = residual[plot_WLs]
      LSF_Fit_array[*, i]       = gaussian(WL[include_WLs], a)                         ; LSF_Fit
      Torus_params[i]          = {params, A[0]*A[2]*SQRT(2*!DPI), stddev(residual[correl_indicies])*A[2]*SQRT(2*!DPI), $
        0, 0, float(sxpar(header, 'Torus_la')), filename_matrix[loc], float(sxpar(header, 'T_P_Shad')), ten(sxpar(header, 'EXPTIME'))*60.}
      K1_fit_params[i]          = {params, A[0]*A[2]*SQRT(2*!DPI), stddev(residual[correl_indicies])*A[2]*SQRT(2*!DPI), $
        0, 0, 0, filename_matrix[loc], float(sxpar(header, 'T_P_Shad')), ten(sxpar(header, 'EXPTIME'))*60.}
      DopplerShift_array[*, i]  = Ios_Airglow
    endfor

    ; Plot the residual and Gaussian fits
    cgplot, WL[plot_WLs], plot_residual_array[*, 0], psym = 0, Ytitle = 'Residual (kR / '+cgsymbol('Angstrom')+')', xtickformat = '(A1)', $
      yr = YR_resid, xr = XR, /nodata, pos = pos[*,1], /noerase
    cgtext, K1_Io_frame - 0.7, 17, "Io's Doppler Shift", charsize = 1.4, alignment = 0.5
    cgtext, K1 + .5, 17, "Telluric [O I]", charsize = 1.4, alignment = 0.5
    for i = 0, n_elements(frames)-1 do begin
      spec     = MRDFITS(reduced_dir + 'R_per_A_Io_Eclipsed_' + frames[i] +'_D.fits', 1, junk_header, /Fscale, /silent )
      WL       = spec.arg
      cgplot, WL[plot_WLs], plot_residual_array[*, i], /OVERPLOT, COLOR = timeColors[i], psym=10
      cgplot, WL[plot_WLs], LSF_Fit_array[*, i], /OVERPLOT, COLOR = timeColors[i]
    endfor

    ; Co-align to Io's Doppler Shift
    shift_by_array = transpose( DopplerShift_array) - $
      rebin(transpose(mean(DopplerShift_array, dim = 2)), n_elements(Eclipse_files), N_elements(Io_Airglow))
    ind = where(io_airglow gt min(WL) and Io_airglow lt max(WL), /NULL, count)
    if count eq 0 then begin
      shifted_WL = fltarr( n_elements(Eclipse_files) )
      for i = 0, n_elements(Eclipse_files)-1 do begin
        header = headfits(Directory + 'Processed\' + Frames[i])
        shifted_WL[i]  = [mean(WL) + mean(WL)* sxpar(header, 'IO_DOPPL') / cspice_clight()]
      endfor
      shift_by = shifted_WL - mean(shifted_WL)
    endif
    if count eq 1 then shift_by = shift_by_array[*,ind]
    if count gt 1 then shift_by = mean(shift_by_array[*,ind], dim = 2)

    aligned_residuals = plot_residual_array
    for i = 0, n_elements(Eclipse_files)-1 do begin
      aligned_residuals[*,i] = interpol( plot_residual_array[*,i], WL[plot_WLs], WL[plot_WLs] - Shift_By[i])
    endfor
    S = size(aligned_residuals, /dim)
    co_add = median(aligned_residuals, dimension = 2, /even)

    ; plot the aligned & co-added spectra
    cgplot, WL[plot_WLs], co_add, pos = pos[*,2], thick=2, /noerase, xr = xr, yr = [-10.99,10.99], /yst, psym=10, $
      xtitle = 'Wavelength ('+cgsymbol('Angstrom')+')', yminor = 2, /nodata, xminor = 10, xticklen = .1
    cgplot, xr, [0.,0.], linestyle = 2, /overplot

    airglow_x = congrid( mean(DopplerShift_array, dim = 2), N_elements(Io_Airglow)*2+1)
    airglow_x = airglow_x[0:-1]
    airglow_y = Reform(Rebin([-15, 15], 2, N_elements(Io_Airglow)), 2*N_elements(Io_Airglow))
    for i =0, N_elements(Io_Airglow)-1 do cgplot, airglow_x[2*i:2*i+1], airglow_y[2*i:2*i+1], /overplot, linestyle = 2, color = 'blue'
    cgplot, WL[plot_WLs], smooth(co_add, 5), thick = 1, /overplot

    cgps_close
  endif

  ; ================================================ Get the Sunlit Na-D and [O I] 6300A Emissions ==================================================================
  if part eq 11.9 then begin

    ; Get the Solar Spectral Irradiance at 1 AU
    READCOL,'C:\IDL\Io\Kurucz_2005_irradthuwl.dat', F='A,A', WL_nm, flux, STRINGSKIP = '#', /Silent ;flux is in W/m2/nm
    start  = where(WL_nm eq '299.100')
    WL_nm = float(WL_nm[start:*])
    flux = float(flux[start:*])

    ; change flux units from W/m^2/nm to photons/(cm^2 s A)
    ; multiply by ((lambda / hc) / 1 W)  * (1 m^2 / 1e4 cm^2) * (1 nm / 10 A)
    conversion = ((WL_nm*1.e-9)/(6.62606957e-34*299792458.D)) * (1./1.e4) * (1./10.)
    flux = flux * conversion      ; Cross-checked this result against Huebner et al. 1992.
    WL_A = temporary(WL_nm) * 10. ; Wavelength from nm into angstroms
    VACTOAIR, WL_A, WL_A_Air      ; Vacuum to air wavelength conversion
    WL_A = temporary(WL_A_Air)

    Frames      = strcompress(string([0,2]), /remove_all)
    Eye         = 'SD'
    spec        = MRDFITS(reduced_dir + 'R_per_A_Io_Sunlit_' + frames[0] + '_' + Eye + '.fits', 1, junk_header, /Fscale, /silent )
    WL          = spec.arg
    XR          = [5888.,5897.5]
    include_WLs = where( (wl gt xr[0]) and (wl lt xr[1]), /NULL )
    plot_WLs               = where( (wl gt xr[0]) and (wl lt xr[1]), /NULL )
    plot_residual_array    = fltarr(N_elements(plot_WLs), n_elements(frames))
    LSF_Fit_array          = fltarr(N_elements(plot_WLs), n_elements(frames))
    Na1Na2_Fit_params      = REPLICATE({params, brightness:0.0, err_brightness:0.0, mult:0.0, smooth:0.0, shift:0.0, background:'', T_P_shadow:0.0, exptime:0.0}, n_elements(frames))

    for i = 0, N_elements(frames) -1 do begin
      spec     = MRDFITS(reduced_dir + 'R_per_A_Io_Sunlit_' + frames[i] + '_' + Eye + '.fits', 1, junk_header, /Fscale, /silent )
      WL       = spec.arg
      header   = headfits(reduced_dir+'R_per_A_Io_Sunlit_' + frames[i] + '_' + Eye + '.fits')
      Na1_Io_frame = Na1 + Na1 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      Na2_Io_frame = Na2 + Na2 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      junk         = min( abs(WL - Na1_Io_frame), Na1_Io_Ind)
      junk         = min( abs(WL - Na2_Io_frame), Na2_Io_Ind)

      ; Define fitting indicies we need to avoid because of airglow itself.
      Ios_Airglow    = [Io_airglow + Io_airglow * sxpar(header, 'IO_DOPPL') / cspice_clight()]                                      ; wavelength locations of all airglow lines, Ignore telluric K
      PEPSI_1sigma_bandwidth = 1. * (mean(xr)/43000.) / 2.3548 ; 1. * FWHM in Angstroms / 2sqrt(2ln2) = 68% of flux enclosed within 1 sigma --> Appropriate for weak telluric airglow
      Io_AG_Ind = [] & Telluric_AG_Ind = []
      for j = 0, N_elements(Ios_Airglow)-1 do Io_AG_ind = [Io_AG_Ind, where(abs(Ios_Airglow[j] - WL) lt 2.5*PEPSI_1sigma_bandwidth, /NULL)] ;  indicies within 2 sigma of an Io airglow line
      for j = 0, N_elements(Telluric_Airglow)-1 do Telluric_AG_ind = [Telluric_AG_Ind, where(abs(Telluric_Airglow[j] - WL) lt 1.5*PEPSI_1sigma_bandwidth, /NULL)] ; indicies within 1.5 sigma of any Telluric airglow
      AG_ind = [Io_AG_Ind, Telluric_AG_Ind]
      AG_ind = AG_ind[UNIQ(AG_ind, SORT(AG_ind))]

      ; Find which sky spectra is best correlated with the data
      correl_indicies  = cgSetDifference(include_Wls, AG_Ind)                       ; Define wavelength indices to use for Amoeba correlation

      ; shift and smooth a solar spectrum until it matches the sunlit Io data
      spec = spec.fun                                                               ; Define spectrum to correlate
      Jup  = interpol(Flux, WL_A, WL)                                               ; Define spectrum to smooth and shift

      Jup = Jup * median(spec[correl_indicies] / Jup[correl_indicies])
      Ftol = 1.e-4 &  Xtol = 1.e-4 & count = 0                                      ; Define initial tolerances for Amoeba
      smooth_and_shift = [4., -20.] & scale  = [2., 5.]                             ; Define initial guess and simplex for Amoeba
      trial_smooth_and_shift = AMOEBAX(Ftol, xtol, function_name='shift_smooth', SCALE = Scale, P0 = smooth_and_shift, FUNCTION_VALUE = fval)
      smooth_and_shift = trial_smooth_and_shift
      smoothed_shifted = interpolate(gauss_smooth(jup, smooth_and_shift[0]), findgen(n_elements(Jup)) - smooth_and_shift[1])

      p1 = mpfitfun('MX_plus_B', smoothed_shifted[correl_indicies], spec[correl_indicies], 0.1*Spec[correl_indicies], $
        [median(spec[correl_indicies] / smoothed_shifted[correl_indicies]), 0.], /NaN, status=status1, /quiet)

      smoothed_shifted = smoothed_shifted*p1[0] + p1[1]
      residual         = (spec - smoothed_shifted) / 1.e3

      window, 0
      cgplot, WL, spec, xr = [5888, 5897], /ynozero, ytitle = 'R / A'
      cgplot, WL, smoothed_shifted, /overplot, color = 'red'
      cgplot, WL[correl_indicies], spec[correl_indicies], /overplot, psym=4

      ; Gaussian fits
      LSF_fitting_ind1       = cgSetDifference(where( abs(wl - Na1_Io_frame) lt 0.3, /NULL), Telluric_AG_Ind) ; fit region within +/- 0.3A of line center, excluding Telluric airglow
      LSF_fitting_ind2       = cgSetDifference(where( abs(wl - Na2_Io_frame) lt 0.3, /NULL), Telluric_AG_Ind) ; fit region within +/- 0.3A of line center, excluding Telluric airglow
      fa            = {x:double(WL[LSF_fitting_ind1]), y:double(residual[LSF_fitting_ind1]), err:double( sqrt(abs( spec[LSF_fitting_ind1] )) )}
      a             = mpfit('Gaussian_for_MPFIT', [50.D, Na1_Io_frame, 0.0521963], PERROR = err_a, funct=fa, STATUS = Did_it_work1, /Quiet, NPEGGED = NPEGGED)
      fb            = {x:double(WL[LSF_fitting_ind2]), y:double(residual[LSF_fitting_ind2]), err:double( sqrt(abs( spec[LSF_fitting_ind2] )) )}
      b             = mpfit('Gaussian_for_MPFIT', [50.D, Na2_Io_frame, 0.0521963], PERROR = err_b, funct=fb, STATUS = Did_it_work2, /Quiet, NPEGGED = NPEGGED)
      plot_residual_array[*, i] = residual[plot_WLs]

      LSF_Fit_array[*, i]       = gaussian(WL[plot_WLs], a) + gaussian(WL[plot_WLs], b)                ; LSF_Fit
      Na1Na2_fit_params[i]      = {params, A[0]*A[2]*SQRT(2*!DPI) + B[0]*B[2]*SQRT(2*!DPI), stddev(residual[correl_indicies])*A[2]*SQRT(2*!DPI) + stddev(residual[correl_indicies])*B[2]*SQRT(2*!DPI), $
                                   0, 0, 0, '', float(sxpar(header, 'T_P_Shad')), ten(sxpar(header, 'EXPTIME'))*60.}

      D2_over_D1     = (A[0]*A[2]*SQRT(2*!DPI)) / (B[0]*B[2]*SQRT(2*!DPI))
      D2_over_D1_err = D2_over_D1 * sqrt( ((stddev(residual[correl_indicies])*A[2]*SQRT(2*!DPI)) / (A[0]*A[2]*SQRT(2*!DPI)))^2 + $
        ((stddev(residual[correl_indicies])*B[2]*SQRT(2*!DPI)) / (B[0]*B[2]*SQRT(2*!DPI)))^2 )
      print, 'Na D2 FWHM linewidths =' + string(a[2]*2.3548) + 'A or' + string(a[2]/Na1*cspice_clight()) + ' km/s'
      print, 'Na D1 FWHM linewidths =' + string(b[2]*2.3548) + 'A or' + string(b[2]/Na2*cspice_clight()) + ' km/s'
      print, 'D2 ='+string(A[0]*A[2]*SQRT(2*!DPI))+', D1 ='+string(B[0]*B[2]*SQRT(2*!DPI))+', D2/D1 line ratio ='+string(D2_over_D1)+' +/-'+string(D2_over_D1_err)

      window, 2
      cgplot, WL, residual, xr = [5888, 5897], /ynozero
      cgplot, WL[plot_WLs], LSF_Fit_array[*, i], /overplot, color = 'red'
    endfor
  endif

  ;==========================================================================================================================================================================
  ;===                                                                                                                                                                  =====
  ;===                                                PLOT TIME SERIES AND BRIGHTNESS WITH RESPECT TO THE TORUS GEOMETRY                                                =====
  ;===                                                                                                                                                                  =====
  ;==========================================================================================================================================================================
  plot_combined:
  if not Keyword_set(date) then begin

    ; Plot temporal lightcurves of the UT 180320 and LBT ingress O and Na D results
    Reduced_Dirs = ['D:\DATA\Apache Point\Echelle\Io Eclipses\UT180320\Reduced\', 'D:\DATA\LBT\Reduced\']
    ;Dates        = ['20 March 2018', '19 April 24']
    Dates        = ['UT180320', 'UT190424']

    for dir_index = 0, N_elements(Reduced_Dirs)-1 do begin
      Reduced_Dir = Reduced_Dirs[dir_index]
      Date        = Dates[dir_index]
      cgPS_Open, filename = Reduced_Dir+'Combined_lightcurve.eps', /ENCAPSULATED, xsize = 5., ysize = 4.
      !P.font=1
      device, SET_FONT = 'Helvetica Bold', /TT_FONT
      pos =  [.12,.19,.9,.9]

      if Reduced_Dir eq 'D:\DATA\LBT\Reduced\' then begin
        yr = [0, 7.]
        restore, Reduced_Dir+'Na1Na2_fit_params.sav'
        restore, Reduced_Dir+'O2_fit_params.sav'
        restore, Reduced_Dir+'O2_Torus_params.sav'
        Penumbra_UTC          = '2019-Apr-24 10:15:12'
        Umbra_UTC             = '2019-Apr-24 10:18:52'
        cspice_UTC2ET, PenUmbra_UTC, PenUmbra_ET
        cspice_UTC2ET, Umbra_UTC, Umbra_ET

        cgplot, O2_fit_params.T_P_Shadow, O2_fit_params.Brightness, psym = 15, Ytitle = 'Disk-Averaged Brightness (kR)', xtitle = 'Minutes After Ingress', $
          title = 'Io''s Response to Ingress: '+Date, yr = yr, xr = [0.,51.], /nodata, pos = pos

        cgLoadCT, 33, NColors=8, /reverse
        ;if ingress then x = findgen(Umbra_ET - PenUmbra_ET) / 60.
        x = findgen(Umbra_ET - PenUmbra_ET) / 60.
        colors = reverse(BytScl(x, MIN=min(x), MAX=2.*max(x)))+128
        FOR j=0,n_elements(x)-2 DO BEGIN
          xpoly = [x[j],         x[j], x[j+1],       x[j+1],         x[j]]
          ypoly = [  0., !Y.CRange[1], !Y.CRange[1], 0., 0.]
          cgColorFill, xpoly, ypoly, Color=colors[j]
        ENDFOR
        ;if ingress then begin
        xpoly = [max(x),     max(x), !X.CRange[1],  !X.CRange[1],  max(x)]
        ypoly = [  0., !Y.CRange[1], !Y.CRange[1], 0., 0.]
        cgColorFill, xpoly, ypoly, color = 'charcoal', /LINE_FILL, orientation = 45., thick = .1
        cgColorFill, xpoly, ypoly, color = 'charcoal', /LINE_FILL, orientation = -45., thick = .1
        ; endif
        cgtext, 2.5, !Y.CRange[1]/4., 'Penumbral Eclipse', orientation = 90., color = 'white'
        cgaxis, yaxis = 0, yr = yr, /ystyle                                      ; repair the axis damage that the Penumbra did
        cgaxis, xaxis = 0, xr = time_range, /xstyle                              ; repair axis damage
        cgaxis, xaxis = 1, xr = time_range, xtickformat = '(A1)', /xstyle        ; repair axis damage

        cgplot, O2_fit_params.T_P_Shadow, O2_fit_params.Brightness, psym = 16, /overplot, color = 'red', /err_clip, $
          ERR_YLOW = O2_fit_params.ERR_Brightness, ERR_YHigh = O2_fit_params.ERR_Brightness, ERR_XLOW = O2_fit_params.exptime/2., ERR_XHigh = O2_fit_params.exptime/2.
        cgplot, O2_fit_params.T_P_Shadow, O2_fit_params.Brightness, color = 'red', /overplot
        cgplot, Na1Na2_fit_params.T_P_Shadow, Na1Na2_fit_params.Brightness, psym = 16, /overplot, color = 'orange', $
          ERR_YLOW = Na1Na2_fit_params.ERR_Brightness, ERR_YHigh = Na1Na2_fit_params.ERR_Brightness, ERR_XLOW = Na1Na2_fit_params.exptime/2., ERR_XHigh = Na1Na2_fit_params.exptime/2.
        cgplot, Na1Na2_fit_params.T_P_Shadow, Na1Na2_fit_params.Brightness, color = 'orange', /overplot


      endif else begin
        yr = [0, 7.1]
        restore, Reduced_Dir+'Io_Airglow_params.sav'
        Ind_6300 = where(Io_Airglow_params.line eq 6300.3042)
        Ind_6364 = where(Io_Airglow_params.line eq 6363.7759)
        Ind_NaD1 = where(Io_Airglow_params.line eq 5895.9243)
        Ind_NaD2 = where(Io_Airglow_params.line eq 5889.9512)

        cgplot, Io_Airglow_params.T_P_Shadow, Io_Airglow_params.Brightness[Ind_6300,*], psym = 15, Ytitle = 'Disk-Averaged Brightness (kR)', xtitle = 'Minutes After Ingress', $
          title = 'Io''s Airglow Response to Ingress', yr = yr, xr = [0.,50.], /nodata, pos = pos, ystyle = 9

        cgLoadCT, 33, NColors=8, /reverse
        x = findgen(Io_Airglow_params.Umbra_ET - Io_Airglow_params.PenUmbra_ET) / 60. ;else x = findgen(PenUmbra_ET - Umbra_ET) / 60.

        colors = reverse(BytScl(x, MIN=min(x), MAX=2.*max(x)))+128
        FOR j=0,n_elements(x)-2 DO BEGIN
          xpoly = [x[j],         x[j], x[j+1],       x[j+1],         x[j]]
          ypoly = [  0., !Y.CRange[1], !Y.CRange[1], 0., 0.]
          cgColorFill, xpoly, ypoly, Color=colors[j]
        ENDFOR
        ;if ingress then begin
        xpoly = [max(x),     max(x), !X.CRange[1],  !X.CRange[1],  max(x)]
        ypoly = [  0., !Y.CRange[1], !Y.CRange[1], 0., 0.]
        cgColorFill, xpoly, ypoly, color = 'charcoal', /LINE_FILL, orientation = 45., thick = .1
        cgColorFill, xpoly, ypoly, color = 'charcoal', /LINE_FILL, orientation = -45., thick = .1
        ;endif

        cgplot, Io_Airglow_params.T_P_Shadow, Io_Airglow_params.Brightness[Ind_6300,*], psym = 16, symsize = 1.1, /overplot, color = 'red', /err_clip, $
          ERR_YLOW  = sqrt(Io_Airglow_params.ERR_Brightness[Ind_6300,*]^2), $
          ERR_YHigh = sqrt(Io_Airglow_params.ERR_Brightness[Ind_6300,*]^2), $
          ERR_XLOW = Io_Airglow_params.exptime/2., ERR_XHigh = Io_Airglow_params.exptime/2., ERR_WIDTH = .005
        cgplot, Io_Airglow_params.T_P_Shadow, Io_Airglow_params.Brightness[Ind_6300,*] , color = 'red', /overplot

        ;cgplot, Io_Airglow_params.T_P_Shadow, Io_Airglow_params.Brightness[Ind_6364,*], psym = 16, /overplot, color = 'red', /err_clip, $
        ;  ERR_YLOW = Io_Airglow_params.ERR_Brightness[Ind_6364,*], ERR_YHigh = Io_Airglow_params.ERR_Brightness[Ind_6364,*], ERR_XLOW = Io_Airglow_params.exptime/2., ERR_XHigh = Io_Airglow_params.exptime/2.
        ;cgplot, Io_Airglow_params.T_P_Shadow, Io_Airglow_params.Brightness[Ind_6364,*], color = 'red', /overplot

        cgplot, Io_Airglow_params.T_P_Shadow, Io_Airglow_params.Brightness[Ind_NaD1,*] + Io_Airglow_params.Brightness[Ind_NaD2,*], psym = 16, symsize = 1.1, /overplot, color = 'orange', /err_clip, $
          ERR_YLOW  = sqrt(Io_Airglow_params.ERR_Brightness[Ind_NaD1,*]^2 + Io_Airglow_params.ERR_Brightness[Ind_NaD2,*]^2), $
          ERR_YHigh = sqrt(Io_Airglow_params.ERR_Brightness[Ind_NaD1,*]^2 + Io_Airglow_params.ERR_Brightness[Ind_NaD2,*]^2), $
          ERR_XLOW = Io_Airglow_params.exptime/2., ERR_XHigh = Io_Airglow_params.exptime/2., ERR_WIDTH = .005
        cgplot, Io_Airglow_params.T_P_Shadow, Io_Airglow_params.Brightness[Ind_NaD1,*] + Io_Airglow_params.Brightness[Ind_NaD2,*], color = 'orange', /overplot

        restore, 'D:\DATA\LBT\Reduced\Na1Na2_fit_params.sav'
        restore, 'D:\DATA\LBT\Reduced\O2_fit_params.sav'
        restore, 'D:\DATA\LBT\Reduced\O2_Torus_params.sav'

        cgplot, O2_fit_params.T_P_Shadow, O2_fit_params.Brightness, psym = 14, symsize = 1.2, /overplot, color = 'red', /err_clip, $
          ERR_YLOW = O2_fit_params.ERR_Brightness, ERR_YHigh = O2_fit_params.ERR_Brightness, ERR_XLOW = O2_fit_params.exptime/2., ERR_XHigh = O2_fit_params.exptime/2.
        cgplot, O2_fit_params.T_P_Shadow, O2_fit_params.Brightness, color = 'red', /overplot
        cgplot, Na1Na2_fit_params.T_P_Shadow, Na1Na2_fit_params.Brightness, psym = 14, symsize = 1.2, /overplot, color = 'orange', /err_clip, $
          ERR_YLOW = Na1Na2_fit_params.ERR_Brightness, ERR_YHigh = Na1Na2_fit_params.ERR_Brightness, ERR_XLOW = Na1Na2_fit_params.exptime/2., ERR_XHigh = Na1Na2_fit_params.exptime/2.
        cgplot, Na1Na2_fit_params.T_P_Shadow, Na1Na2_fit_params.Brightness, color = 'orange', /overplot

        restore, Reduced_Dir+'Io_Airglow_params.sav'

        ;------------------ Mikhail's code ---------------------------------------
        READCOL, Reduced_Dir + 'March20_sun_ec.dat',col1,t_post_eclipse,SO2_val_1,SO2_err_1,SO2_val_2,SO2_err_2,SO2_val_3,SO2_err_3,SO2_val_4,SO2_err_4,SO_val,SO_err, STRINGSKIP = '#', /Silent
        SO2_val_1in = SO2_val_1[0]
        SO2_val_2in = SO2_val_2[0]
        SO2_val_3in = SO2_val_3[0]
        SO2_val_4in = SO2_val_4[0]
        SO_valin = SO_val[0]

        SO2_val_1[0] = !Values.F_Nan
        SO2_val_1[1] = !Values.F_Nan
        SO2_val_2[0] = !Values.F_Nan
        SO2_val_2[1] = !Values.F_Nan
        SO2_val_3[0] = !Values.F_Nan
        SO2_val_3[1] = !Values.F_Nan
        SO2_val_4[0] = !Values.F_Nan
        SO2_val_4[1] = !Values.F_Nan
        SO_val[0] = !Values.F_Nan
        SO_val[1] = !Values.F_Nan

        cgAxis, YAxis=1, YRange=[0, 8.7], title= 'Flux Density [Jy]', COLOR = 'blue', /Save, /ystyle

        cgplot, t_post_eclipse, (SO2_val_1 + SO2_val_2+ SO2_val_3 + SO2_val_4)/1000., psym = 17, symsize = 1.2, /overplot, color = 'blue', $
          ERR_YLOW = (SO2_err_1 + SO2_err_2 + SO2_err_3 + SO2_err_4)/1000., ERR_YHigh = (SO2_err_1 + SO2_err_2 + SO2_err_3 + SO2_err_4)/1000., ERR_XLOW = [0.,0.,1.5,1.825,5.29], ERR_XHigh = [0.,0.,1.5,1.825,5.29], /ystyle
        cgplot, t_post_eclipse, (SO2_val_1 + SO2_val_2 + SO2_val_3 + SO2_val_4)/1000., color = 'blue', /overplot, /ystyle
        cgplot, [0, t_post_eclipse[2]], [(SO2_val_1in + SO2_val_2in + SO2_val_3in + SO2_val_4in)/1000., (SO2_val_1[2] + SO2_val_2[2] + SO2_val_3[2] + SO2_val_4[2])/1000.], color = 'blue', LineStyle=2, /overplot, /ystyle

        cgplot, t_post_eclipse, SO_val/1000., psym = 18, symsize = 1.2, /overplot, color = 'blue', $
          ERR_YLOW = SO2_err_1/1000., ERR_YHigh = SO2_err_1/1000., ERR_XLOW = [0.,0.,1.5,1.825,5.29], ERR_XHigh = [0.,0.,1.5,1.825,5.29], /ystyle
        cgplot, t_post_eclipse, SO_val/1000., color = 'blue', /overplot, /ystyle
        cgplot, [0, t_post_eclipse[2]], [SO_valin/1000., SO_val[2]/1000.], color = 'blue', LineStyle=2, /overplot, /ystyle
        ;------------------ end Mikhail's code ---------------------------------------


        cgtext, 2.5, !Y.CRange[1]/4., 'Penumbral Eclipse', orientation = 90., color = 'white'
        cgaxis, yaxis = 0, yr = yr, /ystyle                                      ; repair the axis damage that the Penumbra did
        cgaxis, xaxis = 0, xr = time_range, /xstyle                              ; repair axis damage
        cgaxis, xaxis = 1, xr = time_range, xtickformat = '(A1)', /xstyle        ; repair axis damage
      endelse

      cglegend, title = ['[O I] 6300'+cgsymbol('Angstrom'), 'Na D!D1!N + D!D2!N', 'SO!D2!N 0.865 + 0.900mm ',  'SO 0.865 + 0.871mm', 'UT190424'], Psym = [16, 16, 17, 18, 14], charsize = 1.2, $
        bg_color = 'white', color = ['red', 'orange', 'blue', 'blue', 'black'], alignment = 1, Location=[1.4, 0.88], /Background, /BOX, vspace = 1.2, length = 0
      ;'UT190424': cglegend, title = ['[O I] 6300'+cgsymbol('Angstrom'), 'Na D1 + D2'], Psym = [16, 16], charsize = 1.2, $
      ; bg_color = 'white', color = ['red', 'orange'], Location=[0.6, 0.88], /Background, /BOX
      ;'UT190812': AL_legend, ['K-D in Umbra 190812'], Psym = 4, /right, charsize = 1.2, /clear

      cgps_close
    endfor ;dir_index

  ; Plot temporal lightcurves of the APO UT180320, LBT ingress O and Na D results
    Reduced_Dirs = ['D:\DATA\Apache Point\Echelle\Io Eclipses\UT180320\Reduced\', 'D:\DATA\LBT\Reduced\']
    Dates        = ['UT180320', 'UT190424']
    for dir_index = 0, N_elements(Reduced_Dirs)-1 do begin
      Reduced_Dir = Reduced_Dirs[dir_index]
      Date        = Dates[dir_index]
      cgPS_Open, filename = Reduced_Dir+'Combined_lightcurve_egress.eps', /ENCAPSULATED, xsize = 6., ysize = 4.
      !P.font=1
      device, SET_FONT = 'Helvetica Bold', /TT_FONT
      pos =  [.12,.19,.9,.9]

      if Reduced_Dir eq 'D:\DATA\LBT\Reduced\' then begin
        yr = [0, 9.]
        restore, Reduced_Dir+'Na1Na2_fit_params.sav'
        restore, Reduced_Dir+'O2_fit_params.sav'
        restore, Reduced_Dir+'O2_Torus_params.sav'
        Penumbra_UTC          = '2019-Apr-24 10:15:12'
        Umbra_UTC             = '2019-Apr-24 10:18:52'
        cspice_UTC2ET, PenUmbra_UTC, PenUmbra_ET
        cspice_UTC2ET, Umbra_UTC, Umbra_ET

        cgplot, O2_fit_params.T_P_Shadow, O2_fit_params.Brightness, psym = 15, Ytitle = 'Disk-Averaged Brightness (kR)', xtitle = 'Minutes After Ingress', $
          title = 'Io''s Airglow Response to Ingress and Egress", yr = yr, xr = [0.,51.], /nodata, pos = pos

        cgLoadCT, 33, NColors=8, /reverse
        ;if ingress then x = findgen(Umbra_ET - PenUmbra_ET) / 60.
        x = findgen(Umbra_ET - PenUmbra_ET) / 60.
        colors = reverse(BytScl(x, MIN=min(x), MAX=2.*max(x)))+128
        FOR j=0,n_elements(x)-2 DO BEGIN
          xpoly = [x[j],         x[j], x[j+1],       x[j+1],         x[j]]
          ypoly = [  0., !Y.CRange[1], !Y.CRange[1], 0., 0.]
          cgColorFill, xpoly, ypoly, Color=colors[j]
        ENDFOR
        ;if ingress then begin
        xpoly = [max(x),     max(x), !X.CRange[1],  !X.CRange[1],  max(x)]
        ypoly = [  0., !Y.CRange[1], !Y.CRange[1], 0., 0.]
        cgColorFill, xpoly, ypoly, color = 'charcoal', /LINE_FILL, orientation = 45., thick = .1
        cgColorFill, xpoly, ypoly, color = 'charcoal', /LINE_FILL, orientation = -45., thick = .1
        ; endif
        cgtext, 2.5, !Y.CRange[1]/4., 'Penumbral Eclipse', orientation = 90., color = 'white'
        cgaxis, yaxis = 0, yr = yr, /ystyle                                      ; repair the axis damage that the Penumbra did
        cgaxis, xaxis = 0, xr = time_range, /xstyle                              ; repair axis damage
        cgaxis, xaxis = 1, xr = time_range, xtickformat = '(A1)', /xstyle        ; repair axis damage

        cgplot, O2_fit_params.T_P_Shadow, O2_fit_params.Brightness, psym = 16, /overplot, color = 'red', /err_clip, $
          ERR_YLOW = O2_fit_params.ERR_Brightness, ERR_YHigh = O2_fit_params.ERR_Brightness, ERR_XLOW = O2_fit_params.exptime/2., ERR_XHigh = O2_fit_params.exptime/2.
        cgplot, O2_fit_params.T_P_Shadow, O2_fit_params.Brightness, color = 'red', /overplot
        cgplot, Na1Na2_fit_params.T_P_Shadow, Na1Na2_fit_params.Brightness, psym = 16, /overplot, color = 'orange', $
          ERR_YLOW = Na1Na2_fit_params.ERR_Brightness, ERR_YHigh = Na1Na2_fit_params.ERR_Brightness, ERR_XLOW = Na1Na2_fit_params.exptime/2., ERR_XHigh = Na1Na2_fit_params.exptime/2.
        cgplot, Na1Na2_fit_params.T_P_Shadow, Na1Na2_fit_params.Brightness, color = 'orange', /overplot
      endif else begin
        yr = [0, 9.]
        restore, Reduced_Dir+'Io_Airglow_params.sav'
        Ind_6300 = where(Io_Airglow_params.line eq 6300.3042)
        Ind_6364 = where(Io_Airglow_params.line eq 6363.7759)
        Ind_NaD1 = where(Io_Airglow_params.line eq 5895.9243)
        Ind_NaD2 = where(Io_Airglow_params.line eq 5889.9512)

        cgplot, Io_Airglow_params.T_P_Shadow, Io_Airglow_params.Brightness[Ind_6300,*], psym = 15, Ytitle = 'Disk-Averaged Brightness (kR)', xtitle = 'Minutes After Ingress', $
          title = 'Io''s Aiglow Response to Ingress and Egress', yr = yr, xr = [0.,140.], /nodata, pos = pos, ystyle = 9

        cgLoadCT, 33, NColors=8, /reverse
        x = findgen(Io_Airglow_params.Umbra_ET - Io_Airglow_params.PenUmbra_ET) / 60. ;else x = findgen(PenUmbra_ET - Umbra_ET) / 60.

        colors = reverse(BytScl(x, MIN=min(x), MAX=2.*max(x)))+128
        FOR j=0,n_elements(x)-2 DO BEGIN
          xpoly = [x[j],         x[j], x[j+1],       x[j+1],         x[j]]
          ypoly = [  0., !Y.CRange[1], !Y.CRange[1], 0., 0.]
          cgColorFill, xpoly, ypoly, Color=colors[j]
        ENDFOR
        ;if ingress then begin
        xpoly = [max(x),     max(x), !X.CRange[1],  !X.CRange[1],  max(x)]
        ypoly = [  0., !Y.CRange[1], !Y.CRange[1], 0., 0.]
        cgColorFill, xpoly, ypoly, color = 'charcoal', /LINE_FILL, orientation = 45., thick = .1
        cgColorFill, xpoly, ypoly, color = 'charcoal', /LINE_FILL, orientation = -45., thick = .1
        ;endif

        cgplot, Io_Airglow_params.T_P_Shadow, Io_Airglow_params.Brightness[Ind_6300,*], psym = 16, /overplot, color = 'red', /err_clip, $
          ERR_YLOW  = sqrt(Io_Airglow_params.ERR_Brightness[Ind_6300,*]^2), $
          ERR_YHigh = sqrt(Io_Airglow_params.ERR_Brightness[Ind_6300,*]^2), $
          ERR_XLOW = Io_Airglow_params.exptime/2., ERR_XHigh = Io_Airglow_params.exptime/2., ERR_WIDTH = .005
        cgplot, Io_Airglow_params.T_P_Shadow, Io_Airglow_params.Brightness[Ind_6300,*] , color = 'red', /overplot

        ;cgplot, Io_Airglow_params.T_P_Shadow, Io_Airglow_params.Brightness[Ind_6364,*], psym = 16, /overplot, color = 'red', /err_clip, $
        ;  ERR_YLOW = Io_Airglow_params.ERR_Brightness[Ind_6364,*], ERR_YHigh = Io_Airglow_params.ERR_Brightness[Ind_6364,*], ERR_XLOW = Io_Airglow_params.exptime/2., ERR_XHigh = Io_Airglow_params.exptime/2.
        ;cgplot, Io_Airglow_params.T_P_Shadow, Io_Airglow_params.Brightness[Ind_6364,*], color = 'red', /overplot

        cgplot, Io_Airglow_params.T_P_Shadow, Io_Airglow_params.Brightness[Ind_NaD1,*] + Io_Airglow_params.Brightness[Ind_NaD2,*], psym = 16, /overplot, color = 'orange', $
          ERR_YLOW  = sqrt(Io_Airglow_params.ERR_Brightness[Ind_NaD1,*]^2 + Io_Airglow_params.ERR_Brightness[Ind_NaD2,*]^2), $
          ERR_YHigh = sqrt(Io_Airglow_params.ERR_Brightness[Ind_NaD1,*]^2 + Io_Airglow_params.ERR_Brightness[Ind_NaD2,*]^2), $
          ERR_XLOW = Io_Airglow_params.exptime/2., ERR_XHigh = Io_Airglow_params.exptime/2., ERR_WIDTH = .005
        cgplot, Io_Airglow_params.T_P_Shadow, Io_Airglow_params.Brightness[Ind_NaD1,*] + Io_Airglow_params.Brightness[Ind_NaD2,*], color = 'orange', /overplot

        cgplot, [130.2, 130.2], [0, 10], color = 'black', LineStyle=2, /overplot

;        ;------------------ Mikhail's code ---------------------------------------
;        restore, 'D:\DATA\Apache Point\Echelle\Io Eclipses\UT190812\Reduced\Na1Na2_fit_params.sav'
;        restore, 'D:\DATA\Apache Point\Echelle\Io Eclipses\UT190812\Reduced\O2_fit_params.sav'
;        restore, 'D:\DATA\Apache Point\Echelle\Io Eclipses\UT190812\Reduced\O2_Torus_params.sav'
;
;        cgplot, O2_fit_params[where(O2_fit_params.T_P_Shadow GT 0)].T_P_Shadow, O2_fit_params[where(O2_fit_params.T_P_Shadow GT 0)].Brightness, psym = 15, /overplot, color = 'red', /err_clip, $
;          ERR_YLOW = O2_fit_params[where(O2_fit_params.T_P_Shadow GT 0)].ERR_Brightness, ERR_YHigh = O2_fit_params[where(O2_fit_params.T_P_Shadow GT 0)].ERR_Brightness, ERR_XLOW = O2_fit_params[where(O2_fit_params.T_P_Shadow GT 0)].exptime/2., ERR_XHigh = O2_fit_params[where(O2_fit_params.T_P_Shadow GT 0)].exptime/2.
;        cgplot, O2_fit_params[where(O2_fit_params.T_P_Shadow GT 0)].T_P_Shadow, O2_fit_params[where(O2_fit_params.T_P_Shadow GT 0)].Brightness, color = 'red', /overplot

        restore, 'D:\DATA\LBT\Reduced\Na1Na2_fit_params.sav'
        restore, 'D:\DATA\LBT\Reduced\O2_fit_params.sav'
        restore, 'D:\DATA\LBT\Reduced\O2_Torus_params.sav'

        cgplot, O2_fit_params.T_P_Shadow, O2_fit_params.Brightness, psym = 14, /overplot, color = 'red', /err_clip, $
          ERR_YLOW = O2_fit_params.ERR_Brightness, ERR_YHigh = O2_fit_params.ERR_Brightness, ERR_XLOW = O2_fit_params.exptime/2., ERR_XHigh = O2_fit_params.exptime/2.
        cgplot, O2_fit_params.T_P_Shadow, O2_fit_params.Brightness, color = 'red', /overplot
        cgplot, Na1Na2_fit_params.T_P_Shadow, Na1Na2_fit_params.Brightness, psym = 14, /overplot, color = 'orange', $
          ERR_YLOW = Na1Na2_fit_params.ERR_Brightness, ERR_YHigh = Na1Na2_fit_params.ERR_Brightness, ERR_XLOW = Na1Na2_fit_params.exptime/2., ERR_XHigh = Na1Na2_fit_params.exptime/2.
        cgplot, Na1Na2_fit_params.T_P_Shadow, Na1Na2_fit_params.Brightness, color = 'orange', /overplot

;        restore, 'D:\DATA\Keck\Io_HIRES_Orders\Katherine\Processed\Io_Airglow_params.sav'
;        Ind_6300 = where(Io_Airglow_params.line eq 6300.3042)
;        Ind_6364 = where(Io_Airglow_params.line eq 6363.7759)
;        Ind_NaD1 = where(Io_Airglow_params.line eq 5895.9243)
;        Ind_NaD2 = where(Io_Airglow_params.line eq 5889.9512)
;
;        cgplot, Io_Airglow_params.T_P_Shadow, Io_Airglow_params.Brightness[Ind_6300,*], psym = 4, /overplot, color = 'red', /err_clip, $
;          ERR_YLOW  = sqrt(Io_Airglow_params.ERR_Brightness[Ind_6300,*]^2), $
;          ERR_YHigh = sqrt(Io_Airglow_params.ERR_Brightness[Ind_6300,*]^2), $
;          ERR_XLOW = Io_Airglow_params.exptime/120., ERR_XHigh = Io_Airglow_params.exptime/120., ERR_WIDTH = .005
;        cgplot, Io_Airglow_params.T_P_Shadow, Io_Airglow_params.Brightness[Ind_6300,*], color = 'red', /overplot
;
;        cgplot, (Io_Airglow_params.T_P_Shadow)[1:*], (Io_Airglow_params.Brightness[Ind_NaD1,*] + Io_Airglow_params.Brightness[Ind_NaD2,*])[1:*], psym = 4, /overplot, color = 'orange', $
;          ERR_YLOW  = (sqrt(Io_Airglow_params.ERR_Brightness[Ind_NaD1,*]^2 + Io_Airglow_params.ERR_Brightness[Ind_NaD2,*]^2))[1:*], $
;          ERR_YHigh = (sqrt(Io_Airglow_params.ERR_Brightness[Ind_NaD1,*]^2 + Io_Airglow_params.ERR_Brightness[Ind_NaD2,*]^2))[1:*], $
;          ERR_XLOW = Io_Airglow_params.exptime/120., ERR_XHigh = Io_Airglow_params.exptime/120., ERR_WIDTH = .005
;        cgplot, (Io_Airglow_params.T_P_Shadow)[1:*], (Io_Airglow_params.Brightness[Ind_NaD1,*] + Io_Airglow_params.Brightness[Ind_NaD2,*])[1:*], color = 'orange', /overplot

        restore, 'D:\DATA\Apache Point\Echelle\Io Eclipses\UT201001\Reduced\Io_Airglow_params.sav'
        Ind_6300 = where(Io_Airglow_params.line eq 6300.3042)
        Ind_6364 = where(Io_Airglow_params.line eq 6363.7759)
        Ind_NaD1 = where(Io_Airglow_params.line eq 5895.9243)
        Ind_NaD2 = where(Io_Airglow_params.line eq 5889.9512)

        cgplot, Io_Airglow_params.T_P_Shadow, (Io_Airglow_params.Brightness[Ind_6300,*]) , psym = 34, /overplot, color = 'red', /err_clip, $
          ERR_YLOW  = (sqrt(Io_Airglow_params.ERR_Brightness[Ind_6300,*]^2)), $
          ERR_YHigh = (sqrt(Io_Airglow_params.ERR_Brightness[Ind_6300,*]^2)) , $
          ERR_XLOW = Io_Airglow_params.exptime /2., ERR_XHigh = Io_Airglow_params.exptime /2., ERR_WIDTH = .005
        cgplot, Io_Airglow_params.T_P_Shadow, (Io_Airglow_params.Brightness[Ind_6300,*]) , color = 'red', /overplot

        cgplot, Io_Airglow_params.T_P_Shadow, Io_Airglow_params.Brightness[Ind_NaD1,*] + Io_Airglow_params.Brightness[Ind_NaD2,*], psym = 34, /overplot, color = 'orange', $
          ERR_YLOW  = sqrt(Io_Airglow_params.ERR_Brightness[Ind_NaD1,*]^2 + Io_Airglow_params.ERR_Brightness[Ind_NaD2,*]^2), $
          ERR_YHigh = sqrt(Io_Airglow_params.ERR_Brightness[Ind_NaD1,*]^2 + Io_Airglow_params.ERR_Brightness[Ind_NaD2,*]^2), $
          ERR_XLOW = Io_Airglow_params.exptime/2., ERR_XHigh = Io_Airglow_params.exptime/2., ERR_WIDTH = .005
        cgplot, Io_Airglow_params.T_P_Shadow, Io_Airglow_params.Brightness[Ind_NaD1,*] + Io_Airglow_params.Brightness[Ind_NaD2,*], color = 'orange', /overplot


        ;restore, 'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200823\Reduced\Io_Airglow_params.sav'
        ;Ind_6300 = where(Io_Airglow_params.line eq 6300.3042)
        ;Ind_6364 = where(Io_Airglow_params.line eq 6363.7759)
        ;Ind_NaD1 = where(Io_Airglow_params.line eq 5895.9243)
        ;Ind_NaD2 = where(Io_Airglow_params.line eq 5889.9512)

        ;cgplot, Io_Airglow_params.T_P_Shadow, Io_Airglow_params.Brightness[Ind_6300,*], psym = 35, /overplot, color = 'red', /err_clip, $
        ; ERR_YLOW  = sqrt(Io_Airglow_params.ERR_Brightness[Ind_6300,*]^2), $
        ; ERR_YHigh = sqrt(Io_Airglow_params.ERR_Brightness[Ind_6300,*]^2), $
        ;ERR_XLOW = Io_Airglow_params.exptime/2., ERR_XHigh = Io_Airglow_params.exptime/2., ERR_WIDTH = .005
        ;cgplot, Io_Airglow_params.T_P_Shadow, Io_Airglow_params.Brightness[Ind_6300,*], color = 'red', /overplot

        ;cgplot, Io_Airglow_params.T_P_Shadow[0:2:2], (Io_Airglow_params.Brightness[Ind_NaD1,*])[0:2:2] + (Io_Airglow_params.Brightness[Ind_NaD2,*])[0:2:2], psym = 35, /overplot, color = 'orange', $
        ; ERR_YLOW  = (sqrt(Io_Airglow_params.ERR_Brightness[Ind_NaD1,*]^2 + Io_Airglow_params.ERR_Brightness[Ind_NaD2,*]^2))[0:2:2], $
        ;ERR_YHigh = (sqrt(Io_Airglow_params.ERR_Brightness[Ind_NaD1,*]^2 + Io_Airglow_params.ERR_Brightness[Ind_NaD2,*]^2))[0:2:2], $
        ; ERR_XLOW = Io_Airglow_params.exptime/2., ERR_XHigh = Io_Airglow_params.exptime/2., ERR_WIDTH = .005
        ;cgplot, Io_Airglow_params.T_P_Shadow[0:2:2], (Io_Airglow_params.Brightness[Ind_NaD1,*])[0:2:2] + (Io_Airglow_params.Brightness[Ind_NaD2,*])[0:2:2], color = 'orange', /overplot

        restore, 'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200908\Reduced\Io_Airglow_params.sav'
        Ind_6300 = where(Io_Airglow_params.line eq 6300.3042)
        Ind_6364 = where(Io_Airglow_params.line eq 6363.7759)
        Ind_NaD1 = where(Io_Airglow_params.line eq 5895.9243)
        Ind_NaD2 = where(Io_Airglow_params.line eq 5889.9512)

        cgplot, (Io_Airglow_params.T_P_Shadow)[0:2], (Io_Airglow_params.Brightness[Ind_6300,*])[0:2], psym = 36, /overplot, color = 'red', /err_clip, $
          ERR_YLOW  = sqrt(Io_Airglow_params.ERR_Brightness[Ind_6300,*]^2), $
          ERR_YHigh = sqrt(Io_Airglow_params.ERR_Brightness[Ind_6300,*]^2), $
          ERR_XLOW = Io_Airglow_params.exptime/2., ERR_XHigh = Io_Airglow_params.exptime/2., ERR_WIDTH = .005
        cgplot, (Io_Airglow_params.T_P_Shadow[0:2]), (Io_Airglow_params.Brightness[Ind_6300,*])[0:2], color = 'red', /overplot

        cgplot, Io_Airglow_params.T_P_Shadow[0:2], (Io_Airglow_params.Brightness[Ind_NaD1,*])[0:2] + (Io_Airglow_params.Brightness[Ind_NaD2,*])[0:2], psym = 36, /overplot, color = 'orange', $
          ERR_YLOW  = (sqrt(Io_Airglow_params.ERR_Brightness[Ind_NaD1,*]^2 + Io_Airglow_params.ERR_Brightness[Ind_NaD2,*]^2))[0:2], $
          ERR_YHigh = (sqrt(Io_Airglow_params.ERR_Brightness[Ind_NaD1,*]^2 + Io_Airglow_params.ERR_Brightness[Ind_NaD2,*]^2))[0:2], $
          ERR_XLOW = Io_Airglow_params.exptime/2., ERR_XHigh = Io_Airglow_params.exptime/2., ERR_WIDTH = .005, err_clip = 1
        cgplot, Io_Airglow_params.T_P_Shadow[0:2], (Io_Airglow_params.Brightness[Ind_NaD1,*])[0:2] + (Io_Airglow_params.Brightness[Ind_NaD2,*])[0:2], color = 'orange', /overplot

        restore, 'D:\DATA\Apache Point\Echelle\Io Eclipses\UT201017\Reduced\Io_Airglow_params.sav'
        Ind_6300 = where(Io_Airglow_params.line eq 6300.3042)
        Ind_6364 = where(Io_Airglow_params.line eq 6363.7759)
        Ind_NaD1 = where(Io_Airglow_params.line eq 5895.9243)
        Ind_NaD2 = where(Io_Airglow_params.line eq 5889.9512)

        cgplot, Io_Airglow_params.T_P_Shadow, (Io_Airglow_params.Brightness[Ind_6300,*]), psym = 23, /overplot, color = 'red', /err_clip, $
          ERR_YLOW  = sqrt(Io_Airglow_params.ERR_Brightness[Ind_6300,*]^2), $
          ERR_YHigh = sqrt(Io_Airglow_params.ERR_Brightness[Ind_6300,*]^2), $
          ERR_XLOW = Io_Airglow_params.exptime/2., ERR_XHigh = Io_Airglow_params.exptime/2., ERR_WIDTH = .005
        cgplot, Io_Airglow_params.T_P_Shadow, (Io_Airglow_params.Brightness[Ind_6300,*]), color = 'red', /overplot

        cgplot, Io_Airglow_params.T_P_Shadow, (Io_Airglow_params.Brightness[Ind_NaD1,*] + Io_Airglow_params.Brightness[Ind_NaD2,*]), psym = 23, /overplot, color = 'orange', $
          ERR_YLOW  = (sqrt(Io_Airglow_params.ERR_Brightness[Ind_NaD1,*]^2 + Io_Airglow_params.ERR_Brightness[Ind_NaD2,*]^2)), $
          ERR_YHigh = (sqrt(Io_Airglow_params.ERR_Brightness[Ind_NaD1,*]^2 + Io_Airglow_params.ERR_Brightness[Ind_NaD2,*]^2)), $
          ERR_XLOW = Io_Airglow_params.exptime/2., ERR_XHigh = Io_Airglow_params.exptime/2., ERR_WIDTH = .005
        cgplot, Io_Airglow_params.T_P_Shadow, (Io_Airglow_params.Brightness[Ind_NaD1,*] + Io_Airglow_params.Brightness[Ind_NaD2,*]), color = 'orange', /overplot

        READCOL, Reduced_Dir + 'March20_sun_ec.dat',col1,t_post_eclipse,SO2_val_1,SO2_err_1,SO2_val_2,SO2_err_2,SO2_val_3,SO2_err_3,SO2_val_4,SO2_err_4,SO_val,SO_err, STRINGSKIP = '#', /Silent
        SO2_val_1in = SO2_val_1[0]
        SO2_val_2in = SO2_val_2[0]
        SO2_val_3in = SO2_val_3[0]
        SO2_val_4in = SO2_val_4[0]
        SO_valin = SO_val[0]

        SO2_val_1[0] = !Values.F_Nan
        SO2_val_1[1] = !Values.F_Nan
        SO2_val_2[0] = !Values.F_Nan
        SO2_val_2[1] = !Values.F_Nan
        SO2_val_3[0] = !Values.F_Nan
        SO2_val_3[1] = !Values.F_Nan
        SO2_val_4[0] = !Values.F_Nan
        SO2_val_4[1] = !Values.F_Nan
        SO_val[0] = !Values.F_Nan
        SO_val[1] = !Values.F_Nan

        cgAxis, YAxis=1, YRange=[0, 8.1], title= 'Flux Density [Jy]', COLOR = 'blue', /Save

        cgplot, t_post_eclipse, (SO2_val_1 + SO2_val_2+ SO2_val_3 + SO2_val_4)/1000., psym = 17, symsize = 0.9, /overplot, color = 'blue', $
          ERR_YLOW = (SO2_err_1 + SO2_err_2 + SO2_err_3 + SO2_err_4)/1000., ERR_YHigh = (SO2_err_1 + SO2_err_2 + SO2_err_3 + SO2_err_4)/1000., ERR_XLOW = [0.,0.,1.5,1.825,5.29], ERR_XHigh = [0.,0.,1.5,1.825,5.29], /ystyle
        cgplot, t_post_eclipse, (SO2_val_1 + SO2_val_2 + SO2_val_3 + SO2_val_4)/1000., color = 'blue', /overplot, /ystyle
        cgplot, [0, t_post_eclipse[2]], [(SO2_val_1in + SO2_val_2in + SO2_val_3in + SO2_val_4in)/1000., (SO2_val_1[2] + SO2_val_2[2] + SO2_val_3[2] + SO2_val_4[2])/1000.], color = 'blue', LineStyle=2, /overplot, /ystyle

        cgplot, t_post_eclipse, SO_val/1000., psym = 16, symsize = 0.9, /overplot, color = 'blue', $
          ERR_YLOW = SO2_err_1/1000., ERR_YHigh = SO2_err_1/1000., ERR_XLOW = [0.,0.,1.5,1.825,5.29], ERR_XHigh = [0.,0.,1.5,1.825,5.29], /ystyle
        cgplot, t_post_eclipse, SO_val/1000., color = 'blue', /overplot, /ystyle
        cgplot, [0, t_post_eclipse[2]], [SO_valin/1000., SO_val[2]/1000.], color = 'blue', LineStyle=2, /overplot, /ystyle
        ;------------------ end Mikhail's code ---------------------------------------

        cgaxis, yaxis = 0, yr = yr, /ystyle                                      ; repair the axis damage that the Penumbra did
        cgaxis, xaxis = 0, xr = time_range, /xstyle                              ; repair axis damage
        cgaxis, xaxis = 1, xr = time_range, xtickformat = '(A1)', /xstyle        ; repair axis damage
      endelse

      ;UT180320': cglegend, title = ['[O I] 6300'+cgsymbol('Angstrom') + ' + 6364'+cgsymbol('Angstrom'), 'Na D!D1!N + D!D2!N', 'SO!D2!N 346.52+346.65Ghz', 'SO!D2!N 332.09+332.50Ghz', 'SO 346.52+344.31Ghz'], Psym = [16, 16, 17, 19, 16], charsize = 1.5, $
      cglegend, title = ['[O I] 6300'+cgsymbol('Angstrom'), 'Na D!D1!N + D!D2!N', 'SO!D2!N spw1+2+5+6',  'SO spw1+3', 'UT190812','UT190424','UT201001','UT200823','UT200908','UT201017', 'UT180807'], Psym = [16, 16, 19, 16,15,14,34,35,36,23, 4], charsize = 1.2, $
        bg_color = 'white', color = ['red', 'orange', 'blue', 'blue', 'red','red','red','red','red','red', 'red'], alignment = 1, Location=[1.3, 0.88], /Background, /BOX, vspace = 1.2, length = 0
      ;'UT190424': cglegend, title = ['[O I] 6300'+cgsymbol('Angstrom'), 'Na D1 + D2'], Psym = [16, 16], charsize = 1.2, $
      ; bg_color = 'white', color = ['red', 'orange'], Location=[0.6, 0.88], /Background, /BOX
      ;'UT190812': AL_legend, ['K-D in Umbra 190812'], Psym = 4, /right, charsize = 1.2, /clear
      cgps_close
    endfor ;dir_index

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;                     plotting [O I] versus torus latitude                       ;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;    cgPS_Open, filename = 'D:\DATA\Apache Point\Echelle\Io Eclipses\UT180320\Reduced\6300_vs_Torus_latitude.eps', /ENCAPSULATED, xsize = 6., ysize = 4.
;      !P.font=1
;      device, SET_FONT = 'Helvetica Bold', /TT_FONT
;      
;      param_files = [ 'D:\DATA\Apache Point\Echelle\Io Eclipses\UT180320\Reduced\Io_Airglow_params.sav', $
;                      'D:\DATA\Keck\Io Eclipse HIRES\Katherine\Reduced\Io_Airglow_params.sav', $
;                      'D:\DATA\LBT\Reduced\O2_fit_params.sav', $
;                      'D:\DATA\Apache Point\Echelle\Io Eclipses\UT190812\Reduced\Io_Airglow_params.sav', $
;                      'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200823\Reduced\Io_Airglow_params.sav', $
;                      'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200908\Reduced\Io_Airglow_params.sav', $
;                      'D:\DATA\Apache Point\Echelle\Io Eclipses\UT201001\Reduced\Io_Airglow_params.sav', $
;                      'D:\DATA\Apache Point\Echelle\Io Eclipses\UT201017\Reduced\Io_Airglow_params.sav' ]
;      
;      Legend_dates = ['180320', '180807', '190424', '190812', '200823', '200908', '201001','201017']
;      Plot_symbols = [16, 21, 14, 6, 4, 36, 23, 35]
;      
;      cgplot, [0,1], [0,1], /nodata, title = '[O I] Red Lines vs Io''s Latitude in the Torus', $
;              xtitle = 'Io''s Latitude from the Torus Centrifugal Equator [deg]', xr = [0.,7.], $
;              ytitle = 'Disk-Averaged 6300 + 6364'+cgsymbol('Angstrom')+' [kR]', yr = [0., 12.], pos = [0.12, 0.16, 0.98, 0.9]
;                                       
;      for i = 0, n_elements(Legend_dates)-1 do begin
;        ; LBT's results are formatted a bit differently than APO and Keck HIRES. the 6364A line is cutoff, so assume a 1 to 3 line ratio
;          if i eq 2 then begin
;            cgplot, abs(Torus_params.shift), O2_fit_params.Brightness * 1.333, psym = plot_symbols[i], /overplot, color = 'red', /err_clip, $
;              ERR_YLOW = O2_fit_params.ERR_Brightness, ERR_YHigh = O2_fit_params.ERR_Brightness, symsize = 1.5
;            continue
;          endif  
;        
;        ; APO and Keck   
;          restore, param_files[i]
;          Ind_6300 = where(Io_Airglow_params.line eq 6300.3042)
;          Ind_6364 = where(Io_Airglow_params.line eq 6363.7759)
;          Red_Sum = Io_Airglow_params.brightness[Ind_6300,*] + Io_Airglow_params.brightness[Ind_6364,*]
;          err_Red_Sum = sqrt(Io_Airglow_params.ERR_Brightness[Ind_6300,*]^2 + Io_Airglow_params.ERR_Brightness[Ind_6364,*]^2)
;          cgplot, abs(Io_Airglow_params.torus_latitude), Red_Sum, psym = plot_symbols[i], /overplot, color = 'red', /err_clip, $
;            ERR_YLOW = err_Red_Sum, ERR_YHIGH = err_Red_Sum , symsize = 1.5
;
;          ;cgplot, abs(Io_Airglow_params.torus_latitude), Io_Airglow_params.brightness[ind_6300,*], psym = plot_symbols[i], /overplot, color = 'red', /err_clip, $
;          ;  ERR_YLOW = Io_Airglow_params.ERR_Brightness[ind_6300,*], ERR_YHigh = Io_Airglow_params.ERR_Brightness[ind_6300,*], symsize = 1.5
;      endfor
;      
;      cglegend, title = Legend_dates, Psym = plot_symbols, charsize = 1.2, bg_color = 'white', color = ['red', 'red', 'red', 'red', 'red', 'red', 'red', 'red'], $
;        alignment = 1, Location=[0.31, 0.6], /Background, /BOX, vspace = 1.2, length = 0, symsize = 1.5
;    cgps_close

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;                     plotting  Na D versus torus latitude                       ;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
;    cgPS_Open, filename = 'D:\DATA\Apache Point\Echelle\Io Eclipses\UT180320\Reduced\Na_D_vs_Torus_latitude.eps', /ENCAPSULATED, xsize = 6., ysize = 4.
;      !P.font=1
;      device, SET_FONT = 'Helvetica Bold', /TT_FONT
;  
;      param_files = [ 'D:\DATA\Apache Point\Echelle\Io Eclipses\UT180320\Reduced\Io_Airglow_params.sav', $
;                      'D:\DATA\Keck\Io Eclipse HIRES\Katherine\Reduced\Io_Airglow_params.sav', $
;                      'D:\DATA\LBT\Reduced\Na1Na2_fit_params.sav', $
;                      'D:\DATA\Apache Point\Echelle\Io Eclipses\UT190812\Reduced\Io_Airglow_params.sav', $
;                      'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200823\Reduced\Io_Airglow_params.sav', $
;                      'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200908\Reduced\Io_Airglow_params.sav', $
;                      'D:\DATA\Apache Point\Echelle\Io Eclipses\UT201001\Reduced\Io_Airglow_params.sav', $
;                      'D:\DATA\Apache Point\Echelle\Io Eclipses\UT201017\Reduced\Io_Airglow_params.sav' ]
;  
;      Legend_dates = ['180320', '180807', '190424', '190812', '200823', '200908', '201001','201017']
;      Plot_symbols = [16, 21, 14, 6, 4, 36, 23, 35]
;  
;      cgplot, [0,1], [0,1], /nodata, title = 'Sodium Aurora vs Io''s Latitude in the Torus', $
;        xtitle = 'Io''s Latitude from the Torus Centrifugal Equator [deg]', xr = [0.,7.], $
;        ytitle = 'Disk-Averaged Na D!D1!N + D!D2!N [kR]', yr = [0., 10.], pos = [0.12, 0.16, 0.98, 0.9]
;  
;      for i = 0, n_elements(Legend_dates)-1 do begin
;        ; LBT's results are formatted a bit differently than APO and Keck HIRES
;          if i eq 2 then begin
;            cgplot, abs(Torus_params.shift), Na1Na2_fit_params.Brightness, psym = plot_symbols[i], /overplot, color = 'orange', /err_clip, $
;              ERR_YLOW = Na1Na2_fit_params.ERR_Brightness, ERR_YHigh = Na1Na2_fit_params.ERR_Brightness, symsize = 1.5
;            continue
;          endif
;  
;        ; APO and Keck
;          restore, param_files[i]
;          Ind_NaD1 = where(Io_Airglow_params.line eq 5895.9243)
;          Ind_NaD2 = where(Io_Airglow_params.line eq 5889.9512)
;          D_sum = Io_Airglow_params.brightness[Ind_NaD2,*] + Io_Airglow_params.brightness[Ind_NaD1,*]
;          err_D_sum = sqrt(Io_Airglow_params.ERR_Brightness[Ind_NaD2,*]^2 + Io_Airglow_params.ERR_Brightness[Ind_NaD1,*]^2)
;          cgplot, abs(Io_Airglow_params.torus_latitude), D_sum, psym = plot_symbols[i], /overplot, color = 'orange', /err_clip, $
;                ERR_YLOW = err_D_sum, ERR_YHIGH = err_D_sum, symsize = 1.5
;      endfor
;  
;      cglegend, title = Legend_dates, Psym = plot_symbols, charsize = 1.2, bg_color = 'white', color = ['orange','orange','orange','orange','orange','orange','orange','orange'], $
;        alignment = 1, Location=[0.31, 0.6], /Background, /BOX, vspace = 1.2, length = 0, symsize = 1.5
;    cgps_close

stop

;    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;    ;plotting Na versus Torus lattitude
;    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;    cgPS_Open, filename = 'D:\DATA\Apache Point\Echelle\Io Eclipses\UT180320\Reduced\NaD_vs_Torus_latitude_Ingress.eps', /ENCAPSULATED, xsize = 6., ysize = 4.
;    !P.font=1
;    device, SET_FONT = 'Helvetica Bold', /TT_FONT
;    yr = [0, 7]
;    pos = [0.120000,0.170000,0.900000,0.900000]
;
;    restore, 'D:\DATA\LBT\Reduced\O2_fit_params.sav'
;    restore, 'D:\DATA\LBT\Reduced\O2_Torus_params.sav'
;    restore, 'D:\DATA\LBT\Reduced\Na1Na2_fit_params.sav'
;    cgplot, abs(Torus_params.shift), O2_fit_params.Brightness, psym = 14, xtitle = 'Latitude from Centrifugal Equator [deg]', $
;      title = 'Na D!D1!N + D!D2!N'  +' vs Io''s Torus Latitude', yr = yr, xr = [1.5,6.7],YTICKFORMAT="(A1)",yticks = 1, /nodata, pos = pos
;    cgAxis, YAxis=0, YRange=yr, title= 'Disk-Averaged Brightness [kR]', COLOR = 'black', /Save, /ystyle
;    cgplot, abs(Torus_params.shift), Na1Na2_fit_params.Brightness, psym = 16, /overplot, color = 'orange', $
;      ERR_YLOW = Na1Na2_fit_params.ERR_Brightness, ERR_YHigh = Na1Na2_fit_params.ERR_Brightness
;    cgplot,abs(Torus_params.shift), Na1Na2_fit_params.Brightness, color = 'orange', /overplot
;
;    restore, 'D:\DATA\Apache Point\Echelle\Io Eclipses\UT180320\Reduced\Io_Airglow_params.sav'
;    Ind_6300 = where(Io_Airglow_params.line eq 6300.3042)
;    Ind_6364 = where(Io_Airglow_params.line eq 6363.7759)
;    Ind_NaD1 = where(Io_Airglow_params.line eq 5895.9243)
;    Ind_NaD2 = where(Io_Airglow_params.line eq 5889.9512)
;    Io_Airglow_params.brightness[ind_6300,0] = !Values.F_Nan ;skip penumbral frame
;    cgplot, abs(Io_Airglow_params.torus_latitude), (Io_Airglow_params.Brightness[Ind_NaD1,*] + Io_Airglow_params.Brightness[Ind_NaD2,*]), psym = 16, /overplot, color = 'orange', /err_clip, $
;      ERR_YLOW = Io_Airglow_params.ERR_Brightness[ind_6300,*], ERR_YHigh = Io_Airglow_params.ERR_Brightness[ind_6300,*]
;    cgplot, abs(Io_Airglow_params.torus_latitude), (Io_Airglow_params.Brightness[Ind_NaD1,*] + Io_Airglow_params.Brightness[Ind_NaD2,*]), color = 'orange', /overplot
;
;    restore, 'D:\DATA\Apache Point\Echelle\Io Eclipses\UT201001\Reduced\Io_Airglow_params.sav'
;    Ind_6300 = where(Io_Airglow_params.line eq 6300.3042)
;    Ind_6364 = where(Io_Airglow_params.line eq 6363.7759)
;    Ind_NaD1 = where(Io_Airglow_params.line eq 5895.9243)
;    Ind_NaD2 = where(Io_Airglow_params.line eq 5889.9512)
;    cgplot, (abs(Io_Airglow_params.torus_latitude)), (Io_Airglow_params.Brightness[Ind_NaD1,*] + Io_Airglow_params.Brightness[Ind_NaD2,*]), psym = 21, /overplot, color = 'orange', /err_clip, $
;      ERR_YLOW = Io_Airglow_params.ERR_Brightness[ind_6300,*], ERR_YHigh = Io_Airglow_params.ERR_Brightness[ind_6300,*]
;    cgplot, (abs(Io_Airglow_params.torus_latitude)), (Io_Airglow_params.Brightness[Ind_NaD1,*] + Io_Airglow_params.Brightness[Ind_NaD2,*]), color = 'orange', /overplot
;
;    restore, 'D:\DATA\Keck\Io Eclipse HIRES\Katherine\Reduced\Io_Airglow_params.sav'
;    Ind_6300 = where(Io_Airglow_params.line eq 6300.3042)
;    Ind_6364 = where(Io_Airglow_params.line eq 6363.7759)
;    Ind_NaD1 = where(Io_Airglow_params.line eq 5895.9243)
;    Ind_NaD2 = where(Io_Airglow_params.line eq 5889.9512)
;    Io_Airglow_params.brightness[ind_6300,0] = !Values.F_Nan ;skip penumbral frame
;    cgplot, (abs(Io_Airglow_params.torus_latitude))[1:*], (Io_Airglow_params.Brightness[Ind_NaD1,*] + Io_Airglow_params.Brightness[Ind_NaD2,*])[1:*], psym = 4, /overplot, color = 'orange', /err_clip, $
;      ERR_YLOW = Io_Airglow_params.ERR_Brightness[ind_6300,*], ERR_YHigh = (Io_Airglow_params.ERR_Brightness[ind_6300,*])[1:*]
;    cgplot, (abs(Io_Airglow_params.torus_latitude))[1:*], (Io_Airglow_params.Brightness[Ind_NaD1,*] + Io_Airglow_params.Brightness[Ind_NaD2,*])[1:*], color = 'orange', /overplot
;
;    restore, 'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200908\Reduced\Io_Airglow_params.sav'
;    Ind_6300 = where(Io_Airglow_params.line eq 6300.3042)
;    Ind_6364 = where(Io_Airglow_params.line eq 6363.7759)
;    Ind_NaD1 = where(Io_Airglow_params.line eq 5895.9243)
;    Ind_NaD2 = where(Io_Airglow_params.line eq 5889.9512)
;
;    cgplot, (abs(Io_Airglow_params.torus_latitude))[0:2], (Io_Airglow_params.Brightness[Ind_NaD1,*])[0:2] + (Io_Airglow_params.Brightness[Ind_NaD2,*])[0:2], psym = 36, /overplot, color = 'orange', /err_clip, $
;      ERR_YLOW  = sqrt(Io_Airglow_params.ERR_Brightness[Ind_6300,*]^2), ERR_YHigh = sqrt(Io_Airglow_params.ERR_Brightness[Ind_6300,*]^2)
;    cgplot, (abs(Io_Airglow_params.torus_latitude))[0:2], (Io_Airglow_params.Brightness[Ind_NaD1,*])[0:2] + (Io_Airglow_params.Brightness[Ind_NaD2,*])[0:2], color = 'orange', /overplot
;
;    restore, 'D:\DATA\Apache Point\Echelle\Io Eclipses\UT201017\Reduced\Io_Airglow_params.sav'
;    Ind_6300 = where(Io_Airglow_params.line eq 6300.3042)
;    Ind_6364 = where(Io_Airglow_params.line eq 6363.7759)
;    Ind_NaD1 = where(Io_Airglow_params.line eq 5895.9243)
;    Ind_NaD2 = where(Io_Airglow_params.line eq 5889.9512)
;
;
;    cgplot, abs(Io_Airglow_params.torus_latitude), (Io_Airglow_params.Brightness[Ind_NaD1,*] + Io_Airglow_params.Brightness[Ind_NaD2,*]), psym = 23, /overplot, color = 'orange', $
;      ERR_YLOW  = (sqrt(Io_Airglow_params.ERR_Brightness[Ind_NaD1,*]^2 + Io_Airglow_params.ERR_Brightness[Ind_NaD2,*]^2)), ERR_YHigh = (sqrt(Io_Airglow_params.ERR_Brightness[Ind_NaD1,*]^2 + Io_Airglow_params.ERR_Brightness[Ind_NaD2,*]^2))
;    cgplot, abs(Io_Airglow_params.torus_latitude), (Io_Airglow_params.Brightness[Ind_NaD1,*] + Io_Airglow_params.Brightness[Ind_NaD2,*]), color = 'orange', /overplot
;
;
;
;    cglegend, title = ['190424', '180320', '201001', '180807', '200908', '201017'], Psym = [14, 16, 21, 4, 36, 23], charsize = 1.2, $
;      bg_color = 'white', color = ['orange', 'orange', 'orange', 'orange', 'orange', 'orange'], alignment = 1, Location=[0.30, 0.88], /Background, /BOX, vspace = 1.2, length = 0
;    cgps_close

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;Plotting Na D in sunlight: currently only plotting 201001
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    cgPS_Open, filename = 'D:\DATA\Apache Point\Echelle\Io Eclipses\UT201001\Reduced\NaD_Sunlit.eps', /ENCAPSULATED, xsize = 6., ysize = 4.
    !P.font=1
    device, SET_FONT = 'Helvetica Bold', /TT_FONT
    yr = [0, 0.16]
    pos = [0.120000,0.170000,0.900000,0.900000]

    restore, 'D:\DATA\Apache Point\Echelle\Io Eclipses\UT201001\Reduced\Io_Airglow_params.sav'
    Ind_NaD1 = where(Io_Airglow_params.line eq 5895.9243)
    Ind_NaD2 = where(Io_Airglow_params.line eq 5889.9512)

    cgplot, Io_Airglow_params.T_P_Shadow, Io_Airglow_params.Brightness[Ind_NaD1,*] + Io_Airglow_params.Brightness[Ind_NaD2,*], psym = 14, xtitle = 'Time after Egress', $
      title = 'Na D!D1!N + D!D2!N'  +' Emerging into Sunlight', yr = yr, xr = [0, 120.],YTICKFORMAT="(A1)",yticks = 1, /nodata, pos = pos
    cgAxis, YAxis=0, YRange=yr, title= 'Disk-Averaged Brightness [MR]', COLOR = 'black', /Save, /ystyle

    cgplot, Io_Airglow_params.T_P_Shadow[3:*] - 133.32, (Io_Airglow_params.Brightness[Ind_NaD1,*] + Io_Airglow_params.Brightness[Ind_NaD2,*])[3:*]/1000., psym = 34, /overplot, color = 'orange', $
      ERR_YLOW  = (sqrt(Io_Airglow_params.ERR_Brightness[Ind_NaD1,*]^2 + Io_Airglow_params.ERR_Brightness[Ind_NaD2,*]^2))[3:*]/1000., $
      ERR_YHigh = (sqrt(Io_Airglow_params.ERR_Brightness[Ind_NaD1,*]^2 + Io_Airglow_params.ERR_Brightness[Ind_NaD2,*]^2))[3:*]/1000., $
      ERR_XLOW = (Io_Airglow_params.exptime/2.)[3:*], ERR_XHigh = (Io_Airglow_params.exptime/2.)[2:*], ERR_WIDTH = .005
    cgplot, Io_Airglow_params.T_P_Shadow[3:*] - 133.32, (Io_Airglow_params.Brightness[Ind_NaD1,*] + Io_Airglow_params.Brightness[Ind_NaD2,*])[3:*]/1000., color = 'orange', /overplot

    cgplot, Io_Airglow_params.T_P_Shadow[0:3] - 133.32, [0,0,0, (Io_Airglow_params.Brightness[Ind_NaD1,*] + Io_Airglow_params.Brightness[Ind_NaD2,*])[3]]/1000., psym = 34, /overplot, color = 'orange', $
      ERR_YLOW  = [(sqrt(Io_Airglow_params.ERR_Brightness[Ind_NaD1,*]^2 + Io_Airglow_params.ERR_Brightness[Ind_NaD2,*]^2))]/1000., $
      ERR_YHigh = [(sqrt(Io_Airglow_params.ERR_Brightness[Ind_NaD1,*]^2 + Io_Airglow_params.ERR_Brightness[Ind_NaD2,*]^2))]/1000., $
      ERR_XLOW = (Io_Airglow_params.exptime/2.)[0:3], ERR_XHigh = (Io_Airglow_params.exptime/2.)[0:3], ERR_WIDTH = .005, /err_clip
    cgplot, Io_Airglow_params.T_P_Shadow[0:3] - 133.32, [0,0,0, (Io_Airglow_params.Brightness[Ind_NaD1,*] + Io_Airglow_params.Brightness[Ind_NaD2,*])[3]]/1000., color = 'orange', /overplot
    
    
    restore, 'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200908\Reduced\Io_Airglow_params.sav'
    Ind_NaD1 = where(Io_Airglow_params.line eq 5895.9243)
    Ind_NaD2 = where(Io_Airglow_params.line eq 5889.9512)
    
    cgplot, Io_Airglow_params.T_P_Shadow - 133, (Io_Airglow_params.Brightness[Ind_NaD1,*] + Io_Airglow_params.Brightness[Ind_NaD2,*])/1000., psym = 34, /overplot, color = 'purple', $
      ERR_YLOW  = (sqrt(Io_Airglow_params.ERR_Brightness[Ind_NaD1,*]^2 + Io_Airglow_params.ERR_Brightness[Ind_NaD2,*]^2))/1000., $
      ERR_YHigh = (sqrt(Io_Airglow_params.ERR_Brightness[Ind_NaD1,*]^2 + Io_Airglow_params.ERR_Brightness[Ind_NaD2,*]^2))/1000., $
      ERR_XLOW = (Io_Airglow_params.exptime/2.), ERR_XHigh = (Io_Airglow_params.exptime/2.), ERR_WIDTH = .005
    cgplot, Io_Airglow_params.T_P_Shadow - 133, (Io_Airglow_params.Brightness[Ind_NaD1,*] + Io_Airglow_params.Brightness[Ind_NaD2,*])/1000., color = 'purple', /overplot
    
    restore, 'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200823\Reduced\Io_Airglow_params.sav'
    Ind_NaD1 = where(Io_Airglow_params.line eq 5895.9243)
    Ind_NaD2 = where(Io_Airglow_params.line eq 5889.9512)

    cgplot, Io_Airglow_params.T_P_Shadow - 132.75, (Io_Airglow_params.Brightness[Ind_NaD1,*] + Io_Airglow_params.Brightness[Ind_NaD2,*])/1000., psym = 34, /overplot, color = 'green', $
      ERR_YLOW  = (sqrt(Io_Airglow_params.ERR_Brightness[Ind_NaD1,*]^2 + Io_Airglow_params.ERR_Brightness[Ind_NaD2,*]^2))/1000., $
      ERR_YHigh = (sqrt(Io_Airglow_params.ERR_Brightness[Ind_NaD1,*]^2 + Io_Airglow_params.ERR_Brightness[Ind_NaD2,*]^2))/1000., $
      ERR_XLOW = (Io_Airglow_params.exptime/2.), ERR_XHigh = (Io_Airglow_params.exptime/2.), ERR_WIDTH = .005
    cgplot, Io_Airglow_params.T_P_Shadow - 132.75, (Io_Airglow_params.Brightness[Ind_NaD1,*] + Io_Airglow_params.Brightness[Ind_NaD2,*])/1000., color = 'green', /overplot
    
    restore, 'D:\DATA\Apache Point\Echelle\Io Eclipses\UT201017\Reduced\Io_Airglow_params.sav'
    Ind_NaD1 = where(Io_Airglow_params.line eq 5895.9243)
    Ind_NaD2 = where(Io_Airglow_params.line eq 5889.9512)

    cgplot, Io_Airglow_params.T_P_Shadow - 133.33, (Io_Airglow_params.Brightness[Ind_NaD1,*] + Io_Airglow_params.Brightness[Ind_NaD2,*])/1000., psym = 34, /overplot, color = 'green', $
      ERR_YLOW  = (sqrt(Io_Airglow_params.ERR_Brightness[Ind_NaD1,*]^2 + Io_Airglow_params.ERR_Brightness[Ind_NaD2,*]^2))/1000., $
      ERR_YHigh = (sqrt(Io_Airglow_params.ERR_Brightness[Ind_NaD1,*]^2 + Io_Airglow_params.ERR_Brightness[Ind_NaD2,*]^2))/1000., $
      ERR_XLOW = (Io_Airglow_params.exptime/2.), ERR_XHigh = (Io_Airglow_params.exptime/2.), ERR_WIDTH = .005
    cgplot, Io_Airglow_params.T_P_Shadow - 132.33, (Io_Airglow_params.Brightness[Ind_NaD1,*] + Io_Airglow_params.Brightness[Ind_NaD2,*])/1000., color = 'green', /overplot


    cgps_close
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;Plotting Na D ratios in sunlight: currently only plotting 201001
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    cgPS_Open, filename = 'D:\DATA\Apache Point\Echelle\Io Eclipses\UT201001\Reduced\NaD_ratios.eps', /ENCAPSULATED, xsize = 6., ysize = 4.
    !P.font=1
    device, SET_FONT = 'Helvetica Bold', /TT_FONT
    yr = [0, 0.16]
    pos = [0.120000,0.170000,0.900000,0.900000]

    restore, 'D:\DATA\Apache Point\Echelle\Io Eclipses\UT201001\Reduced\Io_Airglow_params.sav'
    Ind_6300 = where(Io_Airglow_params.line eq 6300.3042)
    Ind_6364 = where(Io_Airglow_params.line eq 6363.7759)
    Ind_NaD1 = where(Io_Airglow_params.line eq 5895.9243)
    Ind_NaD2 = where(Io_Airglow_params.line eq 5889.9512)

    cgplot, Io_Airglow_params.T_P_Shadow, Io_Airglow_params.Brightness[Ind_NaD1,*] + Io_Airglow_params.Brightness[Ind_NaD2,*], psym = 14, xtitle = 'Ratio of Lines D2/D1', $
      title = 'Na D!D1!N + D!D2!N'  +' Emerging into Sunlight', yr = yr, xr = [1, 2],YTICKFORMAT="(A1)",yticks = 1, /nodata, pos = pos
    cgAxis, YAxis=0, YRange=yr, title= 'Disk-Averaged Brightness [MR]', COLOR = 'black', /Save, /ystyle

    cgplot, (Io_Airglow_params.Brightness[Ind_NaD2,*] / Io_Airglow_params.Brightness[Ind_NaD1,*])[3:*], (Io_Airglow_params.Brightness[Ind_NaD1,*] + Io_Airglow_params.Brightness[Ind_NaD2,*])[3:*]/1000., psym = 34, /overplot, color = 'orange', $
      ERR_YLOW  = (sqrt(Io_Airglow_params.ERR_Brightness[Ind_NaD1,*]^2 + Io_Airglow_params.ERR_Brightness[Ind_NaD2,*]^2))[3:*]/1000., $
      ERR_YHigh = (sqrt(Io_Airglow_params.ERR_Brightness[Ind_NaD1,*]^2 + Io_Airglow_params.ERR_Brightness[Ind_NaD2,*]^2))[3:*]/1000., $
      ERR_WIDTH = .005
    cgplot, (Io_Airglow_params.Brightness[Ind_NaD2,*] / Io_Airglow_params.Brightness[Ind_NaD1,*])[3:*], (Io_Airglow_params.Brightness[Ind_NaD1,*] + Io_Airglow_params.Brightness[Ind_NaD2,*])[3:*]/1000., color = 'orange', /overplot

    restore, 'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200908\Reduced\Io_Airglow_params.sav'
    Ind_NaD1 = where(Io_Airglow_params.line eq 5895.9243)
    Ind_NaD2 = where(Io_Airglow_params.line eq 5889.9512)

    cgplot, (Io_Airglow_params.Brightness[Ind_NaD2,*] / Io_Airglow_params.Brightness[Ind_NaD1,*]), (Io_Airglow_params.Brightness[Ind_NaD1,*] + Io_Airglow_params.Brightness[Ind_NaD2,*])/1000., psym = 34, /overplot, color = 'orange', $
      ERR_YLOW  = (sqrt(Io_Airglow_params.ERR_Brightness[Ind_NaD1,*]^2 + Io_Airglow_params.ERR_Brightness[Ind_NaD2,*]^2))/1000., $
      ERR_YHigh = (sqrt(Io_Airglow_params.ERR_Brightness[Ind_NaD1,*]^2 + Io_Airglow_params.ERR_Brightness[Ind_NaD2,*]^2))/1000., $
      ERR_WIDTH = .005
    cgplot, (Io_Airglow_params.Brightness[Ind_NaD2,*] / Io_Airglow_params.Brightness[Ind_NaD1,*]), (Io_Airglow_params.Brightness[Ind_NaD1,*] + Io_Airglow_params.Brightness[Ind_NaD2,*])/1000., color = 'orange', /overplot

    restore, 'D:\DATA\Apache Point\Echelle\Io Eclipses\UT200823\Reduced\Io_Airglow_params.sav'
    Ind_NaD1 = where(Io_Airglow_params.line eq 5895.9243)
    Ind_NaD2 = where(Io_Airglow_params.line eq 5889.9512)

    cgplot, (Io_Airglow_params.Brightness[Ind_NaD2,*] / Io_Airglow_params.Brightness[Ind_NaD1,*]), (Io_Airglow_params.Brightness[Ind_NaD1,*] + Io_Airglow_params.Brightness[Ind_NaD2,*])/1000., psym = 34, /overplot, color = 'orange', $
      ERR_YLOW  = (sqrt(Io_Airglow_params.ERR_Brightness[Ind_NaD1,*]^2 + Io_Airglow_params.ERR_Brightness[Ind_NaD2,*]^2))/1000., $
      ERR_YHigh = (sqrt(Io_Airglow_params.ERR_Brightness[Ind_NaD1,*]^2 + Io_Airglow_params.ERR_Brightness[Ind_NaD2,*]^2))/1000., $
      ERR_WIDTH = .005
    cgplot, (Io_Airglow_params.Brightness[Ind_NaD2,*] / Io_Airglow_params.Brightness[Ind_NaD1,*]), (Io_Airglow_params.Brightness[Ind_NaD1,*] + Io_Airglow_params.Brightness[Ind_NaD2,*])/1000., color = 'orange', /overplot

   cgps_close

  endif
  stop
end
