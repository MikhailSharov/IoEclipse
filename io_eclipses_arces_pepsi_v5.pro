; Analysis of Io's optical spectra in eclipse from datasets with Apache Point and Large Binocular telescopes
; Written by C. Schmidt, M. Sharov. BU Center for Space Physics, 2019
Function ARCES_Correct_Interweave, Input_spectrum, Input_wavelength, show_spectral_range = spectral_range

  ; Perform correcction this ONLY ON UNITY NORMALIZED SPECTRA
  ; This is the few % correction factor to force the continuum to exactly 1. It should be applied to *all* science data

  Airglow_threshold  = 2.e2 ; we don't care much about faint telluric airglow
  step               = 25   ; how many indices to go through before

  ; First, when fitting continuum, one wants to avoid telluric emissions which could throw this procedure off
  ; Load a telluric emission line list

  Files = file_search('D:\DATA\Solar and Telluric Spectra\Telluric_airglow\UVES_Sky\*.tfits', count = n_files)
  WL = [] & flux = [] & FWHM = []
  for i = 0, N_files-1 do begin
    F    = MRDFITS(files[i], 1, header, /USE_COLNUM )
    WL   = [WL,f.c1]
    flux = [flux,f.c2]
    FWHM = [FWHM,f.c3]
  endfor

  keep   = where(flux gt Airglow_threshold, /Null)
  WL     = WL[Keep]
  Flux   = Flux[Keep]
  AG_ind = []
  for i = 0, N_elements(WL)-1 do AG_ind = [AG_Ind, where(abs(wl[i] - Input_wavelength) lt 0.1, /NULL)] ; grab indicies that are within 0.1 A of any telluric airglow line
  AG_ind = AG_ind[UNIQ(AG_ind, SORT(AG_ind))]

  ; Correct the interweave errors to get a better unity normalization of the continuum
  continuum_ind = []
  for i = 0, N_elements(Input_spectrum) - (N_elements(Input_spectrum) mod step) -1, step do begin
    ind = reverse(sort(Input_spectrum[i:i+step]))
    continuum_ind = [continuum_ind, i+ind[1:5]]
  endfor
  continuum_ind    = continuum_ind[UNIQ(continuum_ind, SORT(continuum_ind))]
  continuum_ind    = cgSetDifference(continuum_ind, AG_ind) ; don't include any indicies
  Correct_IW       = interpol(1./Input_spectrum[continuum_ind], Input_wavelength[continuum_ind], Input_wavelength)

  if keyword_set(spectral_range) then begin
    yr = minmax(Input_spectrum[where((Input_wavelength gt spectral_range[0]) and (Input_wavelength lt spectral_range[1]), /NULL)])
    Window, Title = 'Black = original ''FullSpec'', Grey = Continuum Indices, Blue = Echelle Order Interweave Normalized, Red = Omitted Region of Telluric airglow'
    cgplot, Input_wavelength, Input_spectrum, xr = spectral_range, thick = 2, yr = [yr[0],1.05]
    cgplot, Input_wavelength[continuum_ind], Input_spectrum[continuum_ind], color = 'grey', Psym = 16, /overplot
    cgplot, Input_wavelength, Input_spectrum*correct_IW , color = 'blue', thick = 2, /overplot
    for i = 0, n_elements(WL)-1 do cgplot, [WL[i],WL[i]], [0.,1.1], color = 'Green', linestyle = 2, /overplot
    cgplot, Input_wavelength[AG_ind], Input_spectrum[AG_ind], Color = 'red', psym=16, /overplot
  endif
  return, Correct_IW ; The array for correction. Multiplying Correct_IW should fix interweave...
end

FUNCTION Match_Reference_Spectrum, X, P
  ; Multiply P[0], Add P[1] and Smooth P[2] a reference spectrum (X) until it best matches Y
  return, P[0]*gauss_smooth(X, P[2], /EDGE_TRUNCATE) + P[1]
  ;return, P[0]*X + P[1]
end

FUNCTION Match_Telluric_Absorption, X, P
  ; Multiply P[0], Add P[1] and exponent P[2] a reference spectrum (X) until it best matches Y
  return, P[0]*(X^P[2]) + P[1]
end

FUNCTION Jup_Shift_Smooth, X, S
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
  ;RETURN, P[0] + GAUSS1(X, P[1:3])
end

FUNCTION MX_plus_B, X, P
  ; multipliy P[0] and offset P[1]
  return, P[0]*X+ P[1]
end

function shift_smooth, shift_and_smooth_by
  Common shift_smooth_common, spec, jup, correl_indicies, debug

  ; shift_and_smooth_by[0] = # of pixels to smooth
  ; shift_and_smooth_by[1] = # of pixels to shift

    A                = spec[correl_indicies]
    match_to_scatter = interpolate( GAUSS_SMOOTH(jup, shift_and_smooth_by[0]), findgen(n_elements(jup)) - shift_and_smooth_by[1], Missing = !values.F_NaN)
    B                = match_to_scatter[correl_indicies]
    correl           = correlate( A, B, /double )
    ;print, 'Smooth:', shift_and_smooth_by[0], ', Shift:', shift_and_smooth_by[1], ', Correlation:', correl
    ;if keyword_set(debug) then print, 'Smooth:', shift_and_smooth_by[0], ', Shift:', shift_and_smooth_by[1], ', Correlation:', correl
  return, 1.D / correl
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
                 
Pro Io_Eclipses_ARCES_PEPSI_V5, Part=Part, Date=Date
  Common shift_smooth_common, spec, jup, correl_indicies, debug
  
  ;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ;:Part 0 = loading all spectra + functions [this is used to change from tellurically corrected to non corrected]                                                            ::
  ;:Part 1 = O6300 Waterfall and Time Series                                                                                                                                  ::
  ;:Part 2 = Na D Lines Waterfall and Time Series                                                                                                                             ::
  ;:Part 3 = Combined Light Curve [still fairly manual will work on it]                                                                                                       ::
  ;:Part 4 = OI 10.7 - 9.1                                                                                                                                                    ::
  ;:Part 5 = OI 11.0 - 9.5                                                                                                                                                    ::
  ;:Part 6 = SI 7.9 - 6.5                                                                                                                                                     ::
  ;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  Case date of
    'UT190812': begin
      Penumbra_UTC          = '2019-Aug-12 03:44:40'
      Umbra_UTC             = '2019-Aug-12 03:48:19'
      dir                   = 'D:\DATA\Apache Point\Echelle\Io Eclipses\'+date+'\'
      reduced_dir           = 'D:\DATA\Apache Point\Echelle\Io Eclipses\'+date+'\Reduced\'
      calibration_dir       = 'D:\DATA\Apache Point\Echelle\Io Eclipses\'+date+'\Reduced\'  ; Directory with Stellar Spectra for Telluric Corrections
      Eclipse_files         = ['Io_eclipsed.000'+strcompress(indgen(7)+3, /remove_all), 'Io_eclipsed.00'+strcompress(indgen(2)+10, /remove_all) ]
      Jupiter_Center_File   = 'Jupiter_Disk_Center.0002'
    end
    'UT180320': begin
      ingress               = 1
      Penumbra_UTC          = '2018-Mar-20 10:46:46'
      Umbra_UTC             = '2018-Mar-20 10:50:22'
      dir                   = 'D:\DATA\Apache Point\Echelle\Io Eclipses\'+date+'\'
      reduced_dir           = 'D:\DATA\Apache Point\Echelle\Io Eclipses\'+date+'\Reduced\'
      calibration_dir       = 'D:\DATA\Apache Point\Echelle\Io Eclipses\UT190812\Reduced\'  ; Directory with Stellar Spectra for Telluric Corrections
      Jovian_Scatter_files  = [ 'Jovian_Scatter.0001', 'Jovian_Scatter.0002', 'Jovian_Scatter.0003' ]
      Eclipse_files         = [ 'Io_penumbra.0001', 'Io_eclipsed.000'+strcompress(indgen(6)+1, /remove_all) ]
      Jupiter_Center_File   = 'Jupiter_Center_Spectrum.0001'
    end
    'UT190424': begin
      dir                   = 'D:\DATA\LBT\'
      reduced_dir           = 'D:\DATA\LBT\Reduced\'
      calibration_dir       = 'D:\DATA\LBT\Reduced\'
      Penumbra_UTC          = '2019-Apr-24 10:15:12'
      Umbra_UTC             = '2019-Apr-24 10:18:52'
      Jupiter_Center_File   = ['pepsib.20190424.069.dxt.ffc.rec', 'pepsib.20190424.070.dxt.ffc.rec',  $   ; taken over 4 cross-disperser settings (2,3,4,6) in two exposures 3 minutes apart
                               'pepsir.20190424.049.dxt.ffc.rec', 'pepsir.20190424.050.dxt.ffc.rec']      ; [4265-4800A, 4800-5441A, 5441-6278A, 7419-9067A ]  
                               
      Sunlit_files_D        = [['pepsib.20190424.055.dxt.ffc.rec', 'pepsir.20190424.035.dxt.ffc.rec'], $   ; [2,4]  09:41:44.7 Sunlit
                               ['pepsib.20190424.056.dxt.ffc.rec', 'pepsir.20190424.036.dxt.ffc.rec'], $   ; [3,6]  09:46:57.1 Sunlit
                               ['pepsib.20190424.059.dxt.ffc.rec', 'pepsir.20190424.039.dxt.ffc.rec'], $   ; [2,4]  10:10:09.0 Sunlit
                               ['pepsib.20190424.060.dxt.ffc.rec', 'pepsir.20190424.040.dxt.ffc.rec']]     ; [3,6]  10:15:24.5 Penumbra
                               
      Eclipse_files_D      = [['pepsib.20190424.061.dxt.ffc.rec', 'pepsir.20190424.041.dxt.ffc.rec'], $   ; [2,4]  10:20:40.1
                              ['pepsib.20190424.062.dxt.ffc.rec', 'pepsir.20190424.042.dxt.ffc.rec'], $   ; [3,6]  10:25:55.6
                              ['pepsib.20190424.063.dxt.ffc.rec', 'pepsir.20190424.043.dxt.ffc.rec'], $   ; [2,4]  10:32:53.1
                              ['pepsib.20190424.064.dxt.ffc.rec', 'pepsir.20190424.044.dxt.ffc.rec'], $   ; [3,6]  10:40:08.3
                              ['pepsib.20190424.065.dxt.ffc.rec', 'pepsir.20190424.045.dxt.ffc.rec'], $   ; [2,4]  10:47:36.4
                              ['pepsib.20190424.066.dxt.ffc.rec', 'pepsir.20190424.046.dxt.ffc.rec'], $   ; [3,6]  10:54:52.7
                              ['pepsib.20190424.067.dxt.ffc.rec', 'pepsir.20190424.047.dxt.ffc.rec'], $   ; [2,4]  11:02:49.6
                              ['pepsib.20190424.068.dxt.ffc.rec', 'pepsir.20190424.048.dxt.ffc.rec']]     ; [3,6]  11:10:04.7               
    end
  endcase

  ; Desired time range to plot (minutes since Penumbral ingress)
    CASE DATE OF
      'UT180320': BEGIN
        ingress = 1
        range1 = 0
        range2 = 40
      end
      'UT190424': BEGIN
        ingress = 1
        range1 = 0
        range2 = 51
      end
      'UT190812': begin
        ingress = 0
        range1 = 70
        range2 = 135
      end
    ENDCASE

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

  ; Define ingrees / egress times
    cspice_UTC2ET, PenUmbra_UTC, PenUmbra_ET
    cspice_UTC2ET, Umbra_UTC, Umbra_ET
    
  ; Get line lists for individual emission components
    Io_Airglow         = [Na1, Na2, Na3, Na4, K1, K2, O1, O2, O3, O4, O5, O6, O7, O8, O9, S1, S2, S3]            ; Load an Io emission line list

    Airglow_threshold  = 2.e2 ; we don't care much about faint telluric airglow, so threshold it here 
    Files = file_search('D:\DATA\Solar and Telluric Spectra\Telluric_airglow\UVES_Sky\*.tfits', count = n_files) ; Load a telluric emission line list
    AG_WL = [] & AG_flux = [] & AG_FWHM = []
    for i = 0, N_files-1 do begin
      UVES    = MRDFITS(files[i], 1, header, /USE_COLNUM )
      AG_WL   = [AG_WL,UVES.c1]
      AG_flux = [AG_flux,UVES.c2]
      AG_FWHM = [AG_FWHM,UVES.c3]
    endfor
    keep   = where(AG_flux gt Airglow_threshold, /Null)
    Telluric_Airglow = AG_WL[Keep] & Telluric_Airglow_Flux = AG_Flux[Keep]

    ;read in torus latitude data
    READCOL,'D:\DATA\___Calibration___\Carl_centrifugal_equator.txt', torus_lat_out, torus_deg, skipline = 1, /Silent
    

  ;================================Part 0 Generate Telluric Absorption Spectrum==================================================================================
  if part eq 0 then begin
    
;;; test some PyRaf fixes....
;    window, xsize =1000, ys = 800
;    t = MRDFITS(Calibration_Dir+'\fullspecHD_155379_A0V.0001.ec.fits', 0, header, /fscale, /silent, /unsigned )
;    W = sxpar(header, 'CRVAL1')+findgen(N_elements(t))*sxpar(header, 'Cdelt1')
;    cgplot, w, t, /ynozero, xr = [4500, 5700]
;    
;    tf = MRDFITS(Calibration_Dir+'\sumHD_155379_A0V.0001.ec.fits', 0, header, /fscale, /silent, /unsigned )
;    Wf = sxpar(header, 'CRVAL1')+findgen(N_elements(t))*sxpar(header, 'Cdelt1')
;    ;cgplot, wf, tf/t, /ynozero, xr = [4500, 5700]

    ;----------------------Compare 2 methods to get a decent telluric absorption spectrum---------------------------

    ; Method 1: Assume absorption *is* a unity-normalized blue fast rotator stellar spectrum
    ;           Does pretty well for Oxygen 6300
    ;           Beware that interstellar sodium absorption will screw this up
    Telluric_Absorption_BFR = MRDFITS(Calibration_Dir+'\fullspecHD_159975_BFR.0001.ec.fits', 0, header, /fscale, /silent, /unsigned )
    Telluric_WL_BFR         = sxpar(header, 'CRVAL1')+findgen(N_elements(Telluric_Absorption_BFR))*sxpar(header, 'Cdelt1')
    Telluric_AIRMASS_BFR    = sxpar(header, 'AIRMASS')

    ; Method 2: Divide an A0V stellar type model emission spectrum with a measured A0V
    Telluric_Absorption_A0V = MRDFITS(Calibration_Dir+'\sumfitcontHD_155379_A0V.0001.ec.fits', 0, header, /fscale, /silent, /unsigned )
    Telluric_WL_A0V         = sxpar(header, 'CRVAL1')+findgen(N_elements(Telluric_Absorption_A0V))*sxpar(header, 'Cdelt1')
    Telluric_AIRMASS_A0V    = sxpar(header, 'AIRMASS')

    window, 0, xs = 1900, ys = 1000
    test = Telluric_Absorption_BFR*ARCES_Correct_Interweave(Telluric_Absorption_BFR, Telluric_WL_BFR)      ; Test the correction for interweave
    cgplot, Telluric_WL_BFR, Telluric_Absorption_BFR,  xr = [5400, 6900], yr = [0.8, 1.1], color = 'grey'
    cgplot, Telluric_WL_BFR, Test, /overplot
    cgplot, Telluric_WL_BFR, ARCES_Correct_Interweave(Telluric_Absorption_BFR, Telluric_WL_BFR), /overplot, color = 'blue'

    READCOL,'D:\DATA\___Calibration___\Solar_and_Stellar_Calibration_Spectra\vegallpr25.30000', F='A,A', Model_A0V_WL, Model_A0V_Flux, Skipline = 700000, numline = 500000, /Silent
    Model_A0V_WL = Model_A0V_WL*10.
    Model_A0V_Flux = float(shift(Model_A0V_Flux, -100))
    Model_A0V_Flux = INTERPOL(Model_A0V_Flux, Model_A0V_WL, Telluric_WL_A0V)

    result = ROBUST_POLY_FIT(Telluric_WL_A0V , Telluric_Absorption_A0V, 9)
    Echelle_Order_Interweave = Telluric_Absorption_A0V / poly(Telluric_WL_A0V, result)
    result = ROBUST_POLY_FIT(Telluric_WL_A0V, Model_A0V_flux, 6)
    Model_A0V_norm = Model_A0V_flux / poly(Telluric_WL_A0V, result)

    A0V_norm = MRDFITS(Calibration_Dir+'\FullspecHD_155379_A0V.0001.ec.fits', 0, header, /fscale, /silent, /unsigned )
    Telluric_Absorption_A0V = A0V_norm / Model_A0V_norm

    window, 0, title = 'Unity Normalized A0V Spectral Model (Blue) vs and Measured A0V (Black: HD 155379)'
    cgplot, Telluric_WL_A0V, Model_A0V_norm, color = 'blue', xr = [5570,9300]
    cgplot, Telluric_WL_A0V, A0V_norm, color = 'green', /overplot

    window, 1, title = 'Comparison of Telluric Absorption Spectra: Blue fast rotator (black) vs. A0V Stellar Model (Red)'
    cgplot, Telluric_WL_BFR, Telluric_Absorption_BFR, xr = [5000, 6270]
    cgplot, Telluric_WL_A0V, Telluric_Absorption_A0V, /overplot, color = 'red'
    cgplot, Telluric_WL_A0V, Echelle_Order_Interweave - 0.5, /overplot, color = 'green'

    ;----------------------Flux Calibrate Using Distance and Absolute Spectral Reflectivity of Jupiter---------------------------

    ; Comapare publications of Jupiter's spectral albedo at disk center

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
              telluric_fitting_ind = where( (wl gt 7695.5) and (wl lt 7697.5), /NULL)                                            ; dominated by telluric lines
              ;telluric_fitting_ind = where( ((wl gt 7659.) and (wl lt 7661.)) or ((wl gt 7670.) and (wl lt 7673.)), /NULL)     ; dominated by telluric lines
            end
          (Io_Airglow[i] eq O1):                                                                                             ; no telluric absorption 
          (Io_Airglow[i] eq O2): begin
            telluric_absorption = 1
            spectral_range = [6292, 6308]  ; Wavelength region for inspection
            telluric_fitting_ind = where( ((wl gt 6295.0) and (wl lt 6296.5)) or ((wl gt 6305.5) and (wl lt 6307.0)), /NULL) ; dominated by telluric lines
            end
          (Io_Airglow[i] eq O3):  
          (Io_Airglow[i] eq O4) or (Io_Airglow[i] eq  O5) or (Io_Airglow[i] eq O6):                                          ; no telluric absorption  
          (Io_Airglow[i] eq O7) or (Io_Airglow[i] eq  O8) or (Io_Airglow[i] eq O9):                                          ; Certainly some tellurics but just too messy to remove them!    
          (Io_Airglow[i] eq S1) or (Io_Airglow[i] eq  S2) or (Io_Airglow[i] eq S3):                                          ; Certainly some tellurics but just too messy to remove them!
        endcase
      
      if telluric_absorption eq 0 then continue ; otherwise fix the tellurics
      
      ; Fix any subtle pixel shifts between *known telluric lines* measured in the calibration stars and in Jupiter
        lag        = findgen(21)-10.
        correl_BFR = C_CORRELATE( Rough_Absorption_BFR[telluric_fitting_ind], jup_center[telluric_fitting_ind], lag)
        correl_A0V = C_CORRELATE( Rough_Absorption_A0V[telluric_fitting_ind], jup_center[telluric_fitting_ind], lag)
        yfit_BFR   = MPFITPEAK(lag, correl_BFR, A_BFR, nterms = 3, /positive)
        yfit_A0V   = MPFITPEAK(lag, correl_A0V, A_A0V, nterms = 3, /positive)
        aligned_absorp_BFR = interpolate(Rough_Absorption_BFR, findgen(n_elements(Rough_Absorption_BFR)) - A_BFR[1])
        aligned_absorp_A0V = interpolate(Rough_Absorption_A0V, findgen(n_elements(Rough_Absorption_A0V)) - A_A0V[1])
  
      ; Inspect the alignment of both methods of telluric absorption and Jupiter
        window, 0, title = 'Checking Wavelength Alignment of the Telluric Absorption: Bottom plot looks aligned after the above shift? (Blue = BFR, Oreange = A0V) ', ys = 800
        pos = cglayout([1,2])
        cgplot, lag, correl_BFR, pos = pos[*,0]
        cgplot, lag, yfit_BFR, /overplot, color = 'blue'
        cgplot, lag, correl_A0V, /overplot
        cgplot, lag, yfit_A0V, /overplot, color = 'orange'
        
        cgplot, WL, jup_center, xr = spectral_range, /ynozero, pos = pos[*,1], /noerase
        cgplot, WL, aligned_absorp_BFR, color = 'blue', /overplot
        cgplot, WL, aligned_absorp_A0V, color = 'orange', /overplot

        print, 'USER---> Carefully inspect the Wavelength Alignment of the Telluric Absorption Correction at Jupiter in every Bandpass:', spectral_range

    ; Correct the interweave errors in Jupiter and get a better unity normalization of the continuum
      jup_center_N     = jup_center*ARCES_Correct_Interweave(jup_center, WL)      ; Correct for interweave
      jup_center_err_N = jup_center_err*ARCES_Correct_Interweave(jup_center, WL)  ; Correct for interweave

    ; And do the same echelle order interweave normalization for the telluric (BFR and A0V) spectrum
      aligned_absorp_BFR_n = aligned_absorp_BFR*ARCES_Correct_Interweave(aligned_absorp_BFR, WL)
      aligned_absorp_A0V_n = aligned_absorp_A0V*ARCES_Correct_Interweave(aligned_absorp_A0V, WL)

;    ; Set the region of the interstellar absorption equal to zero
;      na_absorb_ind = where( ((wl gt 5889.75) and (wl lt 5890.7)) or ((wl gt 5895.5) and (wl lt 5896.65)), /NULL)

      parinfo            = replicate( {value:0., fixed: 0b, limited: [0b,0b], limits: fltarr(2) }, 3)
      parinfo[0].limited = 1b                               ; limit differences in multiplier amplitude
      parinfo[0].limits  = [0.2, 5.]                        ; multiplicative range limits
      parinfo[1].limited = 1b                               ; limit additive differences
      parinfo[1].limits  = [-0.9, 0.9]                      ; additive adjustments should be small
      parinfo[2].Fixed   = 1b                               ; Fix spectral smoothing
      parinfo[2].value   = 0.                               ; No spectral smoothing
      p0 = [0.8, 0.2, 0.] ; guess at initial Multiply P[0], Add P[1] and Smooth P[2] to

    ; Fit a y = mx + b function to the spectrum, where x is the telluric absorption scaled for 0 to 1 dynamic range, and y is jupiter, scaled for 0 to 1 dynamic range
      p = mpfitfun('Match_Reference_Spectrum', aligned_absorp_BFR[telluric_fitting_ind], jup_center[telluric_fitting_ind], 0.1*jup_center_err[telluric_fitting_ind], $
        p0, /NaN, status=status, parinfo = parinfo, /quiet)
      telluric_fit_BFR = P[0]*aligned_absorp_BFR + P[1]

      p = mpfitfun('Match_Reference_Spectrum', aligned_absorp_A0V[telluric_fitting_ind], jup_center[telluric_fitting_ind], 0.1*jup_center_err[telluric_fitting_ind], $
        p0, /NaN, status=status, parinfo = parinfo, /quiet) ; Fit a y = mx + b function to the spectrum
      telluric_fit_A0V = P[0]*aligned_absorp_A0V + P[1]
  
      window, 2, xs = 1000, ys = 900, title='Black = Jupiter, Green = Indices where Telluric Absorptions are scaled, Red = Aligned Telluric Absorption, Blue = Scaled absorption to Match Jupiter'
      pos = cgLayout([1,2], OXMargin=[11,3], OYMargin=[9,6], YGap=0)
      sample = where((WL gt spectral_range[0]) and (WL lt spectral_range[1]), /NULL)
      
      yr = minmax(jup_center[sample])  
      cgplot, WL, jup_center_n, xr = spectral_range, yr = [yr[0], 1.1], pos = pos[*,0], xtickformat = '(A1)'
      cgplot, WL, jup_center, /overplot, color = 'grey'
      cgplot, WL[telluric_fitting_ind], jup_center[telluric_fitting_ind], /overplot, color = 'green', psym=14
      cgplot, WL, aligned_absorp_A0V, color = 'red', /overplot
      cgplot, WL, telluric_fit_A0V, color = 'blue', /overplot
      cgtext, spectral_range[0]+1, yr[0]+0.05, 'FullSpec', color = 'grey', charsize = 2
      cgtext, spectral_range[0]+1, yr[0]+0.1,  'Blaze Corrected', color = 'black', charsize = 2
      
      yr = minmax(jup_center[sample]/telluric_fit_BFR[sample])  
      cgplot, WL, jup_center/telluric_fit_BFR, xr = spectral_range, yr = yr, /noerase, pos = pos[*,1], color = 'blue'
      cgplot, WL, jup_center/telluric_fit_A0V, /overplot, color = 'orange'

    ; Keep whichever telluric correction method is closer to unity at every spectral bin
      residual_BFR  = abs(1.-jup_center/telluric_fit_BFR)
      residual_A0V  = abs(1.-jup_center/telluric_fit_A0V)
      junk          = min([[residual_BFR],[residual_A0V]], dim = 2, loc)
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
      ;if (Io_Airglow[i] eq  K2) then stop
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
      Jup_Cal_No_Tell_Corr    = jup_center_N * jup_Sum_Curve / Sensitivity_Curve   ; Convert to Rayleighs per Angstrom, no correction
      cgplot, WL, Jup_Cal_Tell_Corr, xr = spectral_range, /ynozero
      cgplot, WL, Jup_Cal_No_Tell_Corr, color = 'red', /overplot
      MWRFITS, Jup_Cal_Tell_Corr, reduced_dir+'R_per_A_Jupiter_Center.ec.fits', header, /CREATE ;/create overwrites
      ;MWRFITS, Jup_Cal_Tell_Corr_err, reduced_dir+'sig_R_per_A_Jupiter_Center.ec.fits', sig_header, /CREATE, /silent ;/create overwrites for sig values (not fancy mate :( will be fixed later)

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

      ; Penumbral emissions also need Io reflactance subtracted
        if date eq 'UT180320' and i eq 0 then begin
          window, 0
          cgplot, WL, Io_ecl_Cal_Tell_Corr, xr = xr, /ynozero
          Io_sunlit = MRDFITS(reduced_dir+'fullspecIo_free_and_clear.0006.ec.fits', 0, Io_header, /fscale, /silent, /unsigned )     ; fullspec, normalized to one
          Io_sunlit = Io_sunlit * ARCES_Correct_Interweave(Io_sunlit, WL)   
          Io_sunlit = Io_sunlit / Master_telluric^( float(sxpar(Io_header, 'AIRMASS')) / float(sxpar(Jupiter_center_header, 'AIRMASS')) ) 
          Io_sunlit = interpolate(Io_sunlit, findgen(n_elements(Io_sunlit)) - .3)
          scale_Io_reflectance = 8.2e4
          ;scale_Io_reflectance = 9.2e4
          cgplot, WL, Io_sunlit*scale_Io_reflectance, /overplot, linestyle = 2
          Io_ecl_Cal_Tell_Corr = Io_ecl_Cal_Tell_Corr - Io_sunlit*scale_Io_reflectance
          cgplot, WL, Io_ecl_Cal_Tell_Corr*2.+8.e4, /overplot, color = 'red'
        endif

      ; find the instantaneous Earth-Io Doppler Shift
        cspice_UTC2ET, sxpar(raw_header, 'DATE-OBS'), ET
        ET_mid_exposure = ET + float(sxpar(raw_header, 'EXPTIME'))/2.
        cspice_spkezr, 'Io', ET_mid_exposure, 'J2000', 'LT+S', 'Earth', Io_Earth_State, ltime
        theta  = cspice_vsep(Io_Earth_state[0:2], Io_Earth_state[3:5])
        Io_wrt_Earth_Dopplershift = cos(theta) * norm(Io_Earth_State[3:5])
        SXADDPAR, Io_header, 'Io_DOPPL', Io_wrt_Earth_Dopplershift, 'Io-Earth V_radial in km/s (mid exposure)'
        SXADDPAR, Io_header, 'T_PSHADO', (ET_mid_exposure-PenUmbra_ET) / 60., 'Minutes since Penumbral ingress'
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
        SXADDPAR, Io_header, 'Sys3_Lat', Io_SysIII, 'Io''s System III Longitude'
        SXADDPAR, Io_header, 'Sys3_Lon', Io_SysIII_LATITUDE, 'Io''s System III Latitude'
        SXADDPAR, Io_header, 'Torus_Lat', torus_lat, 'JRM09 Dipole Approx'
  
      ; Now do the slit filling factor aperture correction. ***************** ASSUMES EMISSION REGION IS THE SIZE OF IO'S DISK *************
        ARCES_Slit_area = 1.6*3.2                                                                        ; Default ARCES slit size
        Io_Area = !pi * (tan(1821.6 / norm(Io_Earth_State[0:2])) * 206265.)^2                            ; Io's solid angle in square arcseconds
        Io_ecl_Cal_Tell_Corr = Io_ecl_Cal_Tell_Corr * ARCES_Slit_area / Io_Area
        
      MWRFITS, Io_ecl_Cal_Tell_Corr, reduced_dir+'R_per_A_'+Eclipse_files[i]+'.ec.fits', Io_header, /CREATE ; /create overwrites 
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

endif ; Part eq 0

MX_plus_B_parinfo          = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0]}, 2)
MX_plus_B_parinfo[1].fixed = 1                                                                        ; peg the additive "B" component of the MX_Plus_B at zero

if Part eq 1.00 then begin ; Make waterfall plots, O 6300 requires telluric corrections

    ;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ;--------------------------------------------------------Get the O6300 Time Series and Scatter Fit-------------------------------------------------------------------------------
    ;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ; METHOD:
    ;   Loop though all possible spectra that could be fit to jovian scatter and *MANUALLY* choose what works best. Avoid Io and telluric airglow in the fitting.
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
    Case date of
      'UT180320': begin
        Fit_Me_To_The_Scatter    = ['R_per_A_Jupiter_Center.ec.fits', $
                                  'R_per_A_Jovian_Scatter.0001.ec.fits', $
                                  'R_per_A_Jovian_Scatter.0002.ec.fits', $
                                  'R_per_A_Jovian_Scatter.0003.ec.fits' ]
                  end
      'UT190812': begin
        Fit_Me_To_The_Scatter    = ['R_per_A_Jupiter_Center.ec.fits']
                  end
    endcase

    ; setup the plot axis
      spec      = MRDFITS(reduced_dir+'R_per_A_'+Eclipse_files[1]+'.ec.fits', 0, header, /fscale, /unsigned, /silent)
      WL        = sxpar(header, 'CRVAL1')+findgen(N_elements(spec))*sxpar(header, 'Cdelt1')      
      bandwidth = 2.5                                                                               ; half with of the plot's x axis in Angstroms

      Case date of
        'UT180320': begin
                      YR       = [38, 86]
                      YR_resid = [-2, 20]
                      XR       = [O2 + O2 * sxpar(header, 'IO_DOPPL') / cspice_clight() - Bandwidth, O2 + O2 * sxpar(header, 'IO_DOPPL') / cspice_clight() + Bandwidth]
                    end
        'UT190812': begin
                      yr   = [5.5e1, 2e2]
                      YR_resid = [-2, 20]
                      XR       = [O2 + O2 * sxpar(header, 'IO_DOPPL') / cspice_clight() - Bandwidth, O2 + O2 * sxpar(header, 'IO_DOPPL') / cspice_clight() + Bandwidth]
                    end
      endcase

    ; Define arrays
      include_WLs            = where( (wl gt xr[0]) and (wl lt xr[1]), /NULL )                      ; wavelengths to use in the scattered light fitting
      plot_WLs               = where( (wl gt xr[0]) and (wl lt xr[1]), /NULL )                      ; wavelengths to actually plot
      residual_array         = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
      plot_residual_array    = fltarr(N_elements(plot_WLs), n_elements(Eclipse_files))
      LSF_Fit_array          = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
      WL_array               = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
      Brightness_array       = fltarr(n_elements(Eclipse_files))
      err_Brightness_array   = fltarr(n_elements(Eclipse_files))
      T_P_Shadow_array       = fltarr(n_elements(Eclipse_files))
      T_U_Shadow_array       = fltarr(n_elements(Eclipse_files))
      EXPTime_array          = fltarr(n_elements(Eclipse_files))
      DopplerShift_array1    = fltarr(n_elements(Eclipse_files))                                    ; For the D2 line
      Torus_lat              = fltarr(n_elements(Eclipse_files))
      Gof                    = 0.                                                                   ; Goodness of fit (co-addded mean residual)
      GOF_array              = fltarr(N_elements(Fit_Me_To_The_Scatter))
      Use_Scatter            = strarr(N_elements(Eclipse_files))
      Possible_WL_offset     = fltarr(N_elements(Eclipse_files), N_elements(Fit_Me_To_The_Scatter)) ; Save any offsets from Io's expected rest wavelength, then correct any error in the wavelength solution.
      Possible_line_width    = fltarr(N_elements(Eclipse_files), N_elements(Fit_Me_To_The_Scatter)) ; Save the linewidths from all type of scatter subtraction, then average and peg this parameter.
      O2_Fit_params          = REPLICATE({params, brightness:0.0, err_brightness:0.0, mult:0.0, smooth:0.0, shift:0.0, background:'', T_P_shadow:0.0, exptime:0.0}, n_elements(Eclipse_files))
      Torus_params          = REPLICATE({params, brightness:0.0, err_brightness:0.0, mult:0.0, smooth:0.0, shift:0.0, background:'', T_P_shadow:0.0, exptime:0.0}, n_elements(Eclipse_files))
    
    for i = 0, n_elements(Eclipse_files)-1 do begin
      if i eq 0 then continue ; skip penumbra
      spec        = MRDFITS(reduced_dir+'R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, header, /fscale, /unsigned, /silent)
      spec_err    = MRDFITS(reduced_dir+'sig_R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, sig_header, /fscale, /unsigned, /silent)
      WL          = sxpar(header, 'CRVAL1')+findgen(N_elements(spec))*sxpar(header, 'Cdelt1')
      

      for h = 0, n_elements(Fit_Me_To_The_Scatter)-1 do begin
        ; Which Jovian scatter are we using?
          Jup_Cal_Tell_Corr = MRDFITS(reduced_dir+Fit_Me_To_The_Scatter[h], 0, Jupiter_Scatter_header, /fscale, /unsigned, /silent)

        ; Define fitting indicies we need to avoid because of airglow itself.
          Ios_Airglow    = [Io_airglow + Io_airglow * sxpar(header, 'IO_DOPPL') / cspice_clight()]                                      ; wavelength locations of all airglow lines, Ignore telluric K
          ARCES_2sigma_bandwidth = 2. * (mean(xr)/31500.) / 2.3548 ; 2. * FWHM in Angstroms / 2sqrt(2ln2) = 95% of flux enclosed within 2 sigma --> Appropriate for Bright Io emission
          ARCES_1sigma_bandwidth = 1. * (mean(xr)/31500.) / 2.3548 ; 1. * FWHM in Angstroms / 2sqrt(2ln2) = 68% of flux enclosed within 1 sigma --> Appropriate for weak telluric airglow
          Io_AG_Ind = [] & Telluric_AG_Ind = []
          for j = 0, N_elements(Ios_Airglow)-1 do Io_AG_ind = [Io_AG_Ind, where(abs(Ios_Airglow[j] - WL) lt ARCES_2sigma_bandwidth, /NULL)] ;  indicies within 2 sigma of an Io airglow line
          for j = 0, N_elements(Telluric_Airglow)-1 do Telluric_AG_ind = [Telluric_AG_Ind, where(abs(Telluric_Airglow[j] - WL) lt ARCES_1sigma_bandwidth, /NULL)] ; indicies within 1 sigma of any Telluric airglow
          AG_ind = [Io_AG_Ind, Telluric_AG_Ind]
          AG_ind = AG_ind[UNIQ(AG_ind, SORT(AG_ind))]

        ; Fit the Jovian scattered light, use Amoeba/Correlation to match Jupiter in terms of smoothing and shifting, then multiply and add until matches the scattered background
          correl_indicies  = cgSetDifference(include_Wls, AG_Ind)
          Mult_and_add     = mpfitfun('MX_Plus_B', Jup_Cal_Tell_Corr[correl_indicies], Spec[correl_indicies], Spec_err[correl_indicies], parinfo = MX_plus_B_parinfo, $
                                      [mean(Spec[correl_indicies]/Jup_Cal_Tell_Corr[correl_indicies]), 0.], /NaN, status = Scatter_fit_status, PERROR = err_a, /quiet)
          jup              = Mult_and_add[0]*Jup_Cal_Tell_Corr + Mult_and_add[1]

          ; Iterate Amoeba until it gives an answer
            Ftol = 1.e-4 &  Xtol = 1.e-4 & count = 0                                          ; Define initial tolerances for Amoeba
            smooth_and_shift = [2., 0.] & scale  = [2., 2.]                                   ; Define intial guess and simplex for Amoeba
            REPEAT BEGIN
              trial_smooth_and_shift = AMOEBAX(Ftol, xtol, function_name='shift_smooth', SCALE = Scale, P0 = smooth_and_shift, FUNCTION_VALUE = fval)
              if N_elements (N_elements(trial_smooth_and_shift) lt 2) then begin              ; If Amoeba crashes
                xtol = xtol * 1.5 & ftol = ftol * 1.5                                         ; Grow the tolerances
                smooth_and_shift = smooth_and_shift + [RANDOMN(seed), RANDOMN(seed)]          ; Change the guess
                SCALE            = SCALE + [RANDOMN(seed), RANDOMN(seed)]                     ; Change the simplex
                count = count+1
              endif
            ENDREP UNTIL ( (N_elements(trial_smooth_and_shift) eq 2) or (count gt 500) )
            smooth_and_shift = AMOEBAX(Ftol, xtol, function_name='shift_smooth', SCALE = [2., 2.], P0 = smooth_and_shift, FUNCTION_VALUE = fval)
            if N_elements(smooth_and_shift) eq 1 then stop ; check for bugs

          smoothed_shifted = interpolate(gauss_smooth(jup, smooth_and_shift[0]), findgen(n_elements(Jup)) - smooth_and_shift[1])
          Mult_and_add     = mpfitfun('MX_Plus_B', smoothed_shifted[correl_indicies], Spec[correl_indicies], Spec_err[correl_indicies], parinfo = MX_plus_B_parinfo, $
                            [mean(Spec[correl_indicies]/smoothed_shifted[correl_indicies]), 0.], /NaN, status = Scatter_fit_status, PERROR = err_a, /quiet)
          scatter_fit      = Mult_and_add[0]*smoothed_shifted + Mult_and_add[1]
          residual         = Spec - scatter_fit
          GOF_array[h]     = GOF + median(abs(residual[correl_indicies]))
          
        ; Get the linewidth and the determine any systematic wavelength shift from Io's rest
        ; Ultimately we will peg both of these parameters.
          O2_Io_frame               = O2 + O2 * sxpar(header, 'IO_DOPPL') / cspice_clight()
          LSF_fitting_ind1          = cgSetDifference(where( abs(wl - O2_Io_frame) lt 0.3, /NULL), Telluric_AG_Ind)  ; fit region within +/- 0.3A of line center, excluding Telluric airglow    
          initial_guess             = [2.e3, O2_Io_frame, 0.077]  ; rough values are fine here [height, wavelength, linewidth_sigma]
          fa                        = {x:double(WL[LSF_fitting_ind1]), y:double(residual[LSF_fitting_ind1]), err:double( sqrt(abs( residual[LSF_fitting_ind1] )) )}
          a                         = mpfit('Gaussian_for_MPFIT', initial_guess, PERROR = err_a, funct=fa, maxiter=50, STATUS = Did_it_work, /Quiet, NPEGGED = NPEGGED) 
          Possible_WL_offset[i,h]   = a[1] - O2_Io_frame
          Possible_Line_Width[i,h]  = a[2]        
      endfor
      best_GoF = min(GOF_array, best_fit)
      Use_Scatter[i] = Fit_Me_To_The_Scatter[best_fit]
      print, 'Frame: ', Eclipse_files[i], ' is best paired with the Jupiter scatter from: ', Use_Scatter[i], ', Goodness of fit = ', best_GoF
    endfor
    Case date of
      'UT180320': begin
          best_scatter_spectrum = 'R_per_A_Jovian_Scatter.0003.ec.fits' ; Manually input desired scatting spectrum.
      end
      'UT190812': begin
          best_scatter_spectrum = 'R_per_A_Jupiter_Center.ec.fits' ; Manually input desired scatting spectrum.
      end
    endcase
    print, 'Upon inspection, best choice seems to be: ', best_scatter_spectrum
    
    ; -------------------------------------------Waterfall Plot Postscript Io and Fit Jovian scatter--------------------------------------------------------------------------------------------------------------

    ; Okay, now that we've established which Jupiter scatter spectrum best optimzes the residual, we can actually use it for scattered light subtraction

    ; get color versus ingress time & setup plot positions
      timeColors = BytScl(findgen(11), Min=0, Max=11, Top=11)
      Color      = timeColors[0]
      cgLoadCT, 33, NColors=8, /reverse
      pos        = cgLayout([1,2], OXMargin=[10,3], OYMargin=[8,5], YGap=0)
      Pos[1,0]   = Pos[1,0]*.7 & Pos[3,1] = Pos[3,1]*.7

    cgPS_Open, filename = Reduced_Dir+'O6300_scatterfit.eps', /ENCAPSULATED, xsize = 7.5, ysize = 5.
      !P.font = 1
      device, SET_FONT = 'Helvetica Bold', /TT_FONT
      !p.charsize = 1.5
    
      cgplot, WL, spec/1.e3, psym = 0, Ytitle = 'KR / '+cgsymbol('Angstrom'), $             ; Rayleighs to KiloRayleighs
        title = 'Io''s [O I] 6300'+cgsymbol('Angstrom')+' Response to Eclipse',  $
        xr = xr, /nodata, pos = pos[*,0], xtickformat = '(A1)', yr = yr
  
      for i = 0, n_elements(Eclipse_files)-1 do begin
        if i eq 0 then continue ; skip penumbra
        spec         = MRDFITS(reduced_dir+'R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, header, /fscale, /unsigned, /silent)
        spec_err     = MRDFITS(reduced_dir+'sig_R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, sig_header, /fscale, /unsigned, /silent)
        spec         = spec / 1.e3       ; to KR
        spec_err     = spec_err / 1.e3   ; to KR
        WL           = sxpar(header, 'CRVAL1') + findgen(N_elements(spec))*sxpar(header, 'Cdelt1') - median(Possible_WL_offset[1:*,*])
        O2_Io_frame  = O2 + O2 * sxpar(header, 'IO_DOPPL') / cspice_clight()
        junk         = min( abs(WL - O2_Io_frame), O2_Io_Ind)
  
        ; Which Jovian scatter are we using?
          Jup_Cal_Tell_Corr = MRDFITS(reduced_dir+best_scatter_spectrum, 0, Jup_Scatter_header, /fscale, /unsigned, /silent)
  
        ; Define fitting indicies we need to avoid because of airglow itself.
          Ios_Airglow    = [Io_airglow + Io_airglow * sxpar(header, 'IO_DOPPL') / cspice_clight()]                                      ; wavelength locations of all airglow lines, Ignore telluric K
          ARCES_2sigma_bandwidth = 2. * (mean(xr)/31500.) / 2.3548 ; 2. * FWHM in Angstroms / 2sqrt(2ln2) = 95% of flux enclosed within 2 sigma --> Appropriate for Bright Io emission
          ARCES_1sigma_bandwidth = 1. * (mean(xr)/31500.) / 2.3548 ; 1. * FWHM in Angstroms / 2sqrt(2ln2) = 68% of flux enclosed within 1 sigma --> Appropriate for weak telluric airglow
          Io_AG_Ind = [] & Telluric_AG_Ind = []
          for j = 0, N_elements(Ios_Airglow)-1 do Io_AG_ind = [Io_AG_Ind, where(abs(Ios_Airglow[j] - WL) lt ARCES_2sigma_bandwidth, /NULL)] ;  indicies within 2 sigma of an Io airglow line
          for j = 0, N_elements(Telluric_Airglow)-1 do Telluric_AG_ind = [Telluric_AG_Ind, where(abs(Telluric_Airglow[j] - WL) lt ARCES_1sigma_bandwidth, /NULL)] ; indicies within 1 sigma of any Telluric airglow
          AG_ind = [Io_AG_Ind, Telluric_AG_Ind]
          AG_ind = AG_ind[UNIQ(AG_ind, SORT(AG_ind))]
    
        ; Fit the Jovian scattered light, use Amoeba/Correlation to match Jupiter in terms of smoothing and shifting, then multiply and add until matches the scattered background
          correl_indicies  = cgSetDifference(include_Wls, AG_Ind)
          Mult_and_add     = mpfitfun('MX_Plus_B', Jup_Cal_Tell_Corr[correl_indicies], Spec[correl_indicies], Spec_err[correl_indicies], parinfo = MX_plus_B_parinfo, $
                                      [mean(Spec[correl_indicies]/Jup_Cal_Tell_Corr[correl_indicies]), 0.], /NaN, status = Scatter_fit_status, PERROR = err_a, /quiet)
          jup              = Mult_and_add[0]*Jup_Cal_Tell_Corr + Mult_and_add[1]
          
          ; Iterate Amoeba until it gives an answer
            Ftol = 1.e-4 &  Xtol = 1.e-4 & count = 0                                          ; Define initial tolerances for Amoeba
            smooth_and_shift = [2., 0.] & scale  = [2., 2.]                                   ; Define intial guess and simplex for Amoeba
            REPEAT BEGIN
              trial_smooth_and_shift = AMOEBAX(Ftol, xtol, function_name='shift_smooth', SCALE = Scale, P0 = smooth_and_shift, FUNCTION_VALUE = fval)
              if N_elements (N_elements(trial_smooth_and_shift) lt 2) then begin              ; If Amoeba crashes
                xtol = xtol * 1.5 & ftol = ftol * 1.5                                         ; Grow the tolerances
                smooth_and_shift = smooth_and_shift + [RANDOMN(seed), RANDOMN(seed)]          ; Change the guess
                SCALE            = SCALE + [RANDOMN(seed), RANDOMN(seed)]                     ; Change the simplex
                count = count+1
              endif
            ENDREP UNTIL ( (N_elements(trial_smooth_and_shift) eq 2) or (count gt 500) )
            smooth_and_shift = AMOEBAX(Ftol, xtol, function_name='shift_smooth', SCALE = [2., 2.], P0 = smooth_and_shift, FUNCTION_VALUE = fval)
            if N_elements(smooth_and_shift) eq 1 then stop ; check for bugs
          
          smoothed_shifted = interpolate(gauss_smooth(jup, smooth_and_shift[0]), findgen(n_elements(Jup)) - smooth_and_shift[1])
          Mult_and_add     = mpfitfun('MX_Plus_B', smoothed_shifted[correl_indicies], Spec[correl_indicies], Spec_err[correl_indicies], [mean(Spec[correl_indicies]/smoothed_shifted[correl_indicies]), 0.], $
                                      /NaN, status = Scatter_fit_status, PERROR = err_a, /quiet, weights = 1./abs(O2_Io_Ind - correl_indicies)^(1.5), parinfo = MX_plus_B_parinfo)                                  
          scatter_fit      = Mult_and_add[0]*smoothed_shifted + Mult_and_add[1]
          residual         = Spec - scatter_fit
          residual_err     = make_array( N_elements(residual), value = stddev(residual[correl_indicies]) )
  
        ; plot the spectrum, fit and indices used for the fit (if we're in debug mode)
          ;cgplot, WL[correl_indicies], scatter_fit[correl_indicies], psym=14, /overplot
          cgplot, WL, scatter_fit, color = timeColors[i], linestyle = 1, thick = 5, /overplot
          cgplot, WL, spec, color = timeColors[i], thick = 5, /overplot
  
        ; Now do the line fitting to the residual. First define Line Spread Function (LSF) Gaussian fit parameter info
          parinfo               = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0]}, 3)
          parinfo[1].fixed      = 1
          parinfo[1].value      = double(O2_Io_frame)                                                           ; Pin the line's wavelength
          parinfo[2].fixed      = 1
          parinfo[2].value      = median(Possible_line_width[1:*,*])
          LSF_fitting_ind1      = cgSetDifference(where( abs(wl - O2_Io_frame) lt 0.3, /NULL), Telluric_AG_Ind) ; fit region within +/- 0.3A of line center, excluding Telluric airglow  
  
        ; Gaussian fit
          initial_guess = [5.D, parinfo[1].value, parinfo[2].value]
          fa            = {x:double(WL[LSF_fitting_ind1]), y:double(residual[LSF_fitting_ind1]), err:double( sqrt(abs( spec[LSF_fitting_ind1] )) )}
          a             = mpfit('Gaussian_for_MPFIT', initial_guess, PERROR = err_a, funct=fa, parinfo = parinfo, STATUS = Did_it_work, /Quiet, NPEGGED = NPEGGED) 

        ; write the results
          WL_array[*, i]           = WL[include_WLs]
          plot_residual_array[*, i]= residual[plot_WLs]
          LSF_Fit_array[*, i]      = gaussian(WL[include_WLs], a)                         ; LSF_Fit
          O2_fit_params[i]         = {params, A[0]*A[2]*SQRT(2*!DPI), stddev(residual[correl_indicies])*A[2]*SQRT(2*!DPI), $
                                              Mult_and_add[0], smooth_and_shift[0], smooth_and_shift[1], Use_Scatter[i], float(sxpar(header, 'T_PSHADO')), float(sxpar(sig_header, 'EXPTIME')) / 60.}
          Torus_params[i]          = {params, A[0]*A[2]*SQRT(2*!DPI), stddev(residual[correl_indicies])*A[2]*SQRT(2*!DPI), $
                                              0, 0, float(sxpar(header, 'Torus_la')), 0, float(sxpar(header, 'T_PSHADO')), ten(sxpar(header, 'EXPTIME'))*60.}          
          Brightness_array[i]      = A[0]*A[2]*SQRT(2*!DPI)                               ; Area under the Gaussian
          err_Brightness_array[i]  = stddev(residual[correl_indicies])*A[2]*SQRT(2*!DPI)  ; use the std dev of the residual and the width of the line
          T_P_Shadow_array[i]      = sxpar(header, 'T_PSHADO')
          T_U_Shadow_array[i]      = sxpar(header, 'T_USHADO')
          EXPTime_array[i]         = sxpar(sig_header, 'EXPTIME') / 60.                   ; for some reason, this keyword is bad in the regular headers, so use sig_header as a workaround
          DopplerShift_array1[i]   = O2_Io_frame
      endfor
      
      ; Annotate w/ legend & text
        AL_legend, ['Raw Io Spectrum','Jupiter Scattered Light Fit'], Psym = [0,0], linestyle = [0,1], charsize = 1., linsize = 0.5, position = [6300.5, yr[1]*.99]
        IF (date EQ 'UT180320') THEN BEGIN
          cgtext, 6297.8, 81, strcompress(string(T_U_Shadow_array[-1], format = '(F10.1)'), /remove_all) +'min Post Umbral Ingress', charsize = 1.4, alignment = 0, color = timeColors[6]
          cgtext, 6297.8, 41, strcompress(string(T_U_Shadow_array[1], format = '(F10.1)'), /remove_all) +'min Post Umbral Ingress', charsize = 1.4, alignment = 0, color = timeColors[1]
        ENDIF

      ; Plot the residual and Gaussian fits
        cgplot, spec, WL, psym = 0, Ytitle = 'Residual (KR / '+cgsymbol('Angstrom')+')', xtitle = 'Wavelength ('+cgsymbol('Angstrom')+')', $
          yr = YR_resid, xr = XR, /nodata, pos = pos[*,1], /noerase
        cgtext, O2_Io_frame - 0.7, 17, "Io's Doppler Shift", charsize = 1.4, alignment = 0.5
        cgtext, O2 + .5, 17, "Telluric [O I]", charsize = 1.4, alignment = 0.5
        for i = 1, n_elements(Eclipse_files)-1 do begin
          cgplot, WL[plot_WLs], plot_residual_array[*, i], /OVERPLOT, COLOR = timeColors[i], psym=10
          cgplot, WL_array[*, i], LSF_Fit_array[*, i], /OVERPLOT, COLOR = timeColors[i]
          cgplot, [DopplerShift_array1[i], DopplerShift_array1[i]], [YR_resid[1], 18], COLOR = timeColors[i], /overplot
          cgplot, [O2, O2], [YR_resid[1], 3.], linestyle = 1, /overplot
        endfor
    cgps_close
    save, O2_fit_params, filename = Reduced_Dir+'O2_fit_params.sav'
    save, Torus_params, filename = Reduced_Dir+'O2_Torus_params.sav'
  endif
  if Part eq 1.01 then begin ; Make waterfall plots, O 6364 requires telluric corrections

    ;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ;--------------------------------------------------------Get the O6364 Time Series and Scatter Fit-------------------------------------------------------------------------------
    ;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ; METHOD:
    ;   Loop though all possible spectra that could be fit to jovian scatter and *MANUALLY* choose what works best. Avoid Io and telluric airglow in the fitting.
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
    Case date of
      'UT180320': begin
        Fit_Me_To_The_Scatter    = ['R_per_A_Jupiter_Center.ec.fits', $
                                  'R_per_A_Jovian_Scatter.0001.ec.fits', $
                                  'R_per_A_Jovian_Scatter.0002.ec.fits', $
                                  'R_per_A_Jovian_Scatter.0003.ec.fits' ]
                  end
      'UT190812': begin
        Fit_Me_To_The_Scatter    = ['R_per_A_Jupiter_Center.ec.fits']
                  end
    endcase

    ; setup the plot axis
      spec      = MRDFITS(reduced_dir+'R_per_A_'+Eclipse_files[1]+'.ec.fits', 0, header, /fscale, /unsigned, /silent)
      WL        = sxpar(header, 'CRVAL1')+findgen(N_elements(spec))*sxpar(header, 'Cdelt1')
      bandwidth = 2.5                                                                               ; half with of the plot's x axis in Angstroms

    Case date of
      'UT180320': begin
        YR       = [41, 76]
        YR_resid = [-2, 6]
        XR       = [O3 + O3 * sxpar(header, 'IO_DOPPL') / cspice_clight() - Bandwidth, O3 + O3 * sxpar(header, 'IO_DOPPL') / cspice_clight() + Bandwidth]
      end
      'UT190812': begin
        yr   = [5.5e1, 2e2]
        YR_resid = [-2, 6]
        XR       = [O3 + O3 * sxpar(header, 'IO_DOPPL') / cspice_clight() - Bandwidth, O3 + O3 * sxpar(header, 'IO_DOPPL') / cspice_clight() + Bandwidth]
      end
    endcase

    ; Define arrays
    include_WLs            = where( (wl gt xr[0]) and (wl lt xr[1]), /NULL )                      ; wavelengths to use in the scattered light fitting
    plot_WLs               = where( (wl gt xr[0]) and (wl lt xr[1]), /NULL )                      ; wavelengths to actually plot
    residual_array         = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
    plot_residual_array    = fltarr(N_elements(plot_WLs), n_elements(Eclipse_files))
    LSF_Fit_array          = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
    WL_array               = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
    Brightness_array       = fltarr(n_elements(Eclipse_files))
    err_Brightness_array   = fltarr(n_elements(Eclipse_files))
    T_P_Shadow_array       = fltarr(n_elements(Eclipse_files))
    T_U_Shadow_array       = fltarr(n_elements(Eclipse_files))
    EXPTime_array          = fltarr(n_elements(Eclipse_files))
    DopplerShift_array1    = fltarr(n_elements(Eclipse_files))                                    ; For the D2 line
    Gof                    = 0.                                                                   ; Goodness of fit (co-addded mean residual)
    GOF_array              = fltarr(N_elements(Fit_Me_To_The_Scatter))
    Use_Scatter            = strarr(N_elements(Eclipse_files))
    Possible_WL_offset     = fltarr(N_elements(Eclipse_files), N_elements(Fit_Me_To_The_Scatter)) ; Save any offsets from Io's expected rest wavelength, then correct any error in the wavelength solution.
    Possible_line_width    = fltarr(N_elements(Eclipse_files), N_elements(Fit_Me_To_The_Scatter)) ; Save the linewidths from all type of scatter subtraction, then average and peg this parameter.
    O3_Fit_params          = REPLICATE({params, brightness:0.0, err_brightness:0.0, mult:0.0, smooth:0.0, shift:0.0, background:'', T_P_shadow:0.0, exptime:0.0}, n_elements(Eclipse_files))

    for i = 0, n_elements(Eclipse_files)-1 do begin
      if i eq 0 then continue ; skip penumbra
      spec        = MRDFITS(reduced_dir+'R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, header, /fscale, /unsigned, /silent)
      spec_err    = MRDFITS(reduced_dir+'sig_R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, sig_header, /fscale, /unsigned, /silent)
      WL          = sxpar(header, 'CRVAL1')+findgen(N_elements(spec))*sxpar(header, 'Cdelt1')

      for h = 0, n_elements(Fit_Me_To_The_Scatter)-1 do begin
        ; Which Jovian scatter are we using?
        Jup_Cal_Tell_Corr = MRDFITS(reduced_dir+Fit_Me_To_The_Scatter[h], 0, Jupiter_Scatter_header, /fscale, /unsigned, /silent)

        ; Define fitting indicies we need to avoid because of airglow itself.
        Ios_Airglow    = [Io_airglow + Io_airglow * sxpar(header, 'IO_DOPPL') / cspice_clight()]                                      ; wavelength locations of all airglow lines, Ignore telluric K
        ARCES_2sigma_bandwidth = 2. * (mean(xr)/31500.) / 2.3548 ; 2. * FWHM in Angstroms / 2sqrt(2ln2) = 95% of flux enclosed within 2 sigma --> Appropriate for Bright Io emission
        ARCES_1sigma_bandwidth = 1. * (mean(xr)/31500.) / 2.3548 ; 1. * FWHM in Angstroms / 2sqrt(2ln2) = 68% of flux enclosed within 1 sigma --> Appropriate for weak telluric airglow
        Io_AG_Ind = [] & Telluric_AG_Ind = []
        for j = 0, N_elements(Ios_Airglow)-1 do Io_AG_ind = [Io_AG_Ind, where(abs(Ios_Airglow[j] - WL) lt ARCES_2sigma_bandwidth, /NULL)] ;  indicies within 2 sigma of an Io airglow line
        for j = 0, N_elements(Telluric_Airglow)-1 do Telluric_AG_ind = [Telluric_AG_Ind, where(abs(Telluric_Airglow[j] - WL) lt ARCES_1sigma_bandwidth, /NULL)] ; indicies within 1 sigma of any Telluric airglow
        AG_ind = [Io_AG_Ind, Telluric_AG_Ind]
        AG_ind = AG_ind[UNIQ(AG_ind, SORT(AG_ind))]

        ; Fit the Jovian scattered light, use Amoeba/Correlation to match Jupiter in terms of smoothing and shifting, then multiply and add until matches the scattered background
        correl_indicies  = cgSetDifference(include_Wls, AG_Ind)
        Mult_and_add     = mpfitfun('MX_Plus_B', Jup_Cal_Tell_Corr[correl_indicies], Spec[correl_indicies], Spec_err[correl_indicies], parinfo = MX_plus_B_parinfo, $
          [mean(Spec[correl_indicies]/Jup_Cal_Tell_Corr[correl_indicies]), 0.], /NaN, status = Scatter_fit_status, PERROR = err_a, /quiet)
        jup              = Mult_and_add[0]*Jup_Cal_Tell_Corr + Mult_and_add[1]

        ; Iterate Amoeba until it gives an answer
        Ftol = 1.e-4 &  Xtol = 1.e-4 & count = 0                                          ; Define initial tolerances for Amoeba
        smooth_and_shift = [2., 0.] & scale  = [2., 2.]                                   ; Define intial guess and simplex for Amoeba
        REPEAT BEGIN
          trial_smooth_and_shift = AMOEBAX(Ftol, xtol, function_name='shift_smooth', SCALE = Scale, P0 = smooth_and_shift, FUNCTION_VALUE = fval)
          if N_elements (N_elements(trial_smooth_and_shift) lt 2) then begin              ; If Amoeba crashes
            xtol = xtol * 1.5 & ftol = ftol * 1.5                                         ; Grow the tolerances
            smooth_and_shift = smooth_and_shift + [RANDOMN(seed), RANDOMN(seed)]          ; Change the guess
            SCALE            = SCALE + [RANDOMN(seed), RANDOMN(seed)]                     ; Change the simplex
            count = count+1
          endif
        ENDREP UNTIL ( (N_elements(trial_smooth_and_shift) eq 2) or (count gt 500) )
        smooth_and_shift = AMOEBAX(Ftol, xtol, function_name='shift_smooth', SCALE = [2., 2.], P0 = smooth_and_shift, FUNCTION_VALUE = fval)
        if N_elements(smooth_and_shift) eq 1 then stop ; check for bugs

        smoothed_shifted = interpolate(gauss_smooth(jup, smooth_and_shift[0]), findgen(n_elements(Jup)) - smooth_and_shift[1])
        Mult_and_add     = mpfitfun('MX_Plus_B', smoothed_shifted[correl_indicies], Spec[correl_indicies], Spec_err[correl_indicies], parinfo = MX_plus_B_parinfo, $
          [mean(Spec[correl_indicies]/smoothed_shifted[correl_indicies]), 0.], /NaN, status = Scatter_fit_status, PERROR = err_a, /quiet)
        scatter_fit      = Mult_and_add[0]*smoothed_shifted + Mult_and_add[1]
        residual         = Spec - scatter_fit
        GOF_array[h]     = GOF + median(abs(residual[correl_indicies]))

        ; Get the linewidth and the determine any systematic wavelength shift from Io's rest
        ; Ultimately we will peg both of these parameters.
        O3_Io_frame               = O3 + O3 * sxpar(header, 'IO_DOPPL') / cspice_clight()
        LSF_fitting_ind1          = cgSetDifference(where( abs(wl - O3_Io_frame) lt 0.3, /NULL), Telluric_AG_Ind)  ; fit region within +/- 0.3A of line center, excluding Telluric airglow
        initial_guess             = [2.e3, O3_Io_frame, 0.077]  ; rough values are fine here [height, wavelength, linewidth_sigma]
        fa                        = {x:double(WL[LSF_fitting_ind1]), y:double(residual[LSF_fitting_ind1]), err:double( sqrt(abs( residual[LSF_fitting_ind1] )) )}
        a                         = mpfit('Gaussian_for_MPFIT', initial_guess, PERROR = err_a, funct=fa, maxiter=50, STATUS = Did_it_work, /Quiet, NPEGGED = NPEGGED)
        Possible_WL_offset[i,h]   = a[1] - O3_Io_frame
        Possible_Line_Width[i,h]  = a[2]
      endfor
      best_GoF = min(GOF_array, best_fit)
      Use_Scatter[i] = Fit_Me_To_The_Scatter[best_fit]
      print, 'Frame: ', Eclipse_files[i], ' is best paired with the Jupiter scatter from: ', Use_Scatter[i], ', Goodness of fit = ', best_GoF
    endfor
    Case date of
      'UT180320': begin
          best_scatter_spectrum = 'R_per_A_Jovian_Scatter.0003.ec.fits' ; Manually input desired scatting spectrum.
      end
      'UT190812': begin
          best_scatter_spectrum = 'R_per_A_Jupiter_Center.ec.fits' ; Manually input desired scatting spectrum.
      end
    endcase
    print, 'Upon inspection, best choice seems to be: ', best_scatter_spectrum

    ; -------------------------------------------Waterfall Plot Postscript Io and Fit Jovian scatter--------------------------------------------------------------------------------------------------------------

    ; Okay, now that we've established which Jupiter scatter spectrum best optimzes the residual, we can actually use it for scattered light subtraction

    ; get color versus ingress time & setup plot positions
    timeColors = BytScl(findgen(11), Min=0, Max=11, Top=11)
    Color      = timeColors[0]
    cgLoadCT, 33, NColors=8, /reverse
    pos        = cgLayout([1,2], OXMargin=[10,3], OYMargin=[8,5], YGap=0)
    Pos[1,0]   = Pos[1,0]*.7 & Pos[3,1] = Pos[3,1]*.7

    cgPS_Open, filename = Reduced_Dir+'O6364_scatterfit.eps', /ENCAPSULATED, xsize = 7.5, ysize = 5.
    !P.font = 1
    device, SET_FONT = 'Helvetica Bold', /TT_FONT
    !p.charsize = 1.5

    cgplot, WL, spec/1.e3, psym = 0, Ytitle = 'KR / '+cgsymbol('Angstrom'), $             ; Rayleighs to KiloRayleighs
      title = 'Io''s [O I] 6364'+cgsymbol('Angstrom')+' Response to Eclipse',  $
      xr = xr, /nodata, pos = pos[*,0], xtickformat = '(A1)', yr = yr

    for i = 0, n_elements(Eclipse_files)-1 do begin
      if i eq 0 then continue ; skip penumbra
      spec         = MRDFITS(reduced_dir+'R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, header, /fscale, /unsigned, /silent)
      spec_err     = MRDFITS(reduced_dir+'sig_R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, sig_header, /fscale, /unsigned, /silent)
      spec         = spec / 1.e3       ; to KR
      spec_err     = spec_err / 1.e3   ; to KR
      WL           = sxpar(header, 'CRVAL1') + findgen(N_elements(spec))*sxpar(header, 'Cdelt1') - median(Possible_WL_offset[1:*,*])
      O3_Io_frame  = O3 + O3 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      junk         = min( abs(WL - O3_Io_frame), O3_Io_Ind)

      ; Which Jovian scatter are we using?
      Jup_Cal_Tell_Corr = MRDFITS(reduced_dir+best_scatter_spectrum, 0, Jup_Scatter_header, /fscale, /unsigned, /silent)

      ; Define fitting indicies we need to avoid because of airglow itself.
      Ios_Airglow    = [Io_airglow + Io_airglow * sxpar(header, 'IO_DOPPL') / cspice_clight()]                                      ; wavelength locations of all airglow lines, Ignore telluric K
      ARCES_2sigma_bandwidth = 2. * (mean(xr)/31500.) / 2.3548 ; 2. * FWHM in Angstroms / 2sqrt(2ln2) = 95% of flux enclosed within 2 sigma --> Appropriate for Bright Io emission
      ARCES_1sigma_bandwidth = 1. * (mean(xr)/31500.) / 2.3548 ; 1. * FWHM in Angstroms / 2sqrt(2ln2) = 68% of flux enclosed within 1 sigma --> Appropriate for weak telluric airglow
      Io_AG_Ind = [] & Telluric_AG_Ind = []
      for j = 0, N_elements(Ios_Airglow)-1 do Io_AG_ind = [Io_AG_Ind, where(abs(Ios_Airglow[j] - WL) lt ARCES_2sigma_bandwidth, /NULL)] ;  indicies within 2 sigma of an Io airglow line
      for j = 0, N_elements(Telluric_Airglow)-1 do Telluric_AG_ind = [Telluric_AG_Ind, where(abs(Telluric_Airglow[j] - WL) lt ARCES_1sigma_bandwidth, /NULL)] ; indicies within 1 sigma of any Telluric airglow
      AG_ind = [Io_AG_Ind, Telluric_AG_Ind]
      AG_ind = AG_ind[UNIQ(AG_ind, SORT(AG_ind))]

      ; Fit the Jovian scattered light, use Amoeba/Correlation to match Jupiter in terms of smoothing and shifting, then multiply and add until matches the scattered background
      correl_indicies  = cgSetDifference(include_Wls, AG_Ind)
      Mult_and_add     = mpfitfun('MX_Plus_B', Jup_Cal_Tell_Corr[correl_indicies], Spec[correl_indicies], Spec_err[correl_indicies], parinfo = MX_plus_B_parinfo, $
        [mean(Spec[correl_indicies]/Jup_Cal_Tell_Corr[correl_indicies]), 0.], /NaN, status = Scatter_fit_status, PERROR = err_a, /quiet)
      jup              = Mult_and_add[0]*Jup_Cal_Tell_Corr + Mult_and_add[1]

      ; Iterate Amoeba until it gives an answer
      Ftol = 1.e-4 &  Xtol = 1.e-4 & count = 0                                          ; Define initial tolerances for Amoeba
      smooth_and_shift = [2., 0.] & scale  = [2., 2.]                                   ; Define intial guess and simplex for Amoeba
      REPEAT BEGIN
        trial_smooth_and_shift = AMOEBAX(Ftol, xtol, function_name='shift_smooth', SCALE = Scale, P0 = smooth_and_shift, FUNCTION_VALUE = fval)
        if N_elements (N_elements(trial_smooth_and_shift) lt 2) then begin              ; If Amoeba crashes
          xtol = xtol * 1.5 & ftol = ftol * 1.5                                         ; Grow the tolerances
          smooth_and_shift = smooth_and_shift + [RANDOMN(seed), RANDOMN(seed)]          ; Change the guess
          SCALE            = SCALE + [RANDOMN(seed), RANDOMN(seed)]                     ; Change the simplex
          count = count+1
        endif
      ENDREP UNTIL ( (N_elements(trial_smooth_and_shift) eq 2) or (count gt 500) )
      smooth_and_shift = AMOEBAX(Ftol, xtol, function_name='shift_smooth', SCALE = [2., 2.], P0 = smooth_and_shift, FUNCTION_VALUE = fval)
      if N_elements(smooth_and_shift) eq 1 then stop ; check for bugs

      smoothed_shifted = interpolate(gauss_smooth(jup, smooth_and_shift[0]), findgen(n_elements(Jup)) - smooth_and_shift[1])
      Mult_and_add     = mpfitfun('MX_Plus_B', smoothed_shifted[correl_indicies], Spec[correl_indicies], Spec_err[correl_indicies], [mean(Spec[correl_indicies]/smoothed_shifted[correl_indicies]), 0.], $
        /NaN, status = Scatter_fit_status, PERROR = err_a, /quiet, weights = 1./abs(O3_Io_Ind - correl_indicies)^(1.5), parinfo = MX_plus_B_parinfo)
      scatter_fit      = Mult_and_add[0]*smoothed_shifted + Mult_and_add[1]
      residual         = Spec - scatter_fit
      residual_err     = make_array( N_elements(residual), value = stddev(residual[correl_indicies]) )

      ; plot the spectrum, fit and indices used for the fit (if we're in debug mode)
      ;cgplot, WL[correl_indicies], scatter_fit[correl_indicies], psym=14, /overplot
      cgplot, WL, scatter_fit, color = timeColors[i], linestyle = 1, thick = 7, /overplot
      cgplot, WL, spec, color = timeColors[i], thick = 5, /overplot

      ; Now do the line fitting to the residual. First define Line Spread Function (LSF) Gaussian fit parameter info
        parinfo               = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0]}, 3)
        parinfo[1].fixed      = 1
        parinfo[1].value      = double(O3_Io_frame)                                                           ; Pin the line's wavelength
        parinfo[2].fixed      = 1
        parinfo[2].value      = median(Possible_line_width[1:*,*])
        LSF_fitting_ind1      = cgSetDifference(where( abs(wl - O3_Io_frame) lt 0.3, /NULL), Telluric_AG_Ind) ; fit region within +/- 0.3A of line center, excluding Telluric airglow

      ; Gaussian fit
        initial_guess = [5.D, parinfo[1].value, parinfo[2].value]
        fa            = {x:double(WL[LSF_fitting_ind1]), y:double(residual[LSF_fitting_ind1]), err:double( sqrt(abs( spec[LSF_fitting_ind1] )) )}
        a             = mpfit('Gaussian_for_MPFIT', initial_guess, PERROR = err_a, funct=fa, parinfo = parinfo, STATUS = Did_it_work, /Quiet, NPEGGED = NPEGGED)

      ; write the results
        WL_array[*, i]           = WL[include_WLs]
        plot_residual_array[*, i]= residual[plot_WLs]
        LSF_Fit_array[*, i]      = gaussian(WL[include_WLs], a)                         ; LSF_Fit
        O3_fit_params[i]         = {params, A[0]*A[2]*SQRT(2*!DPI), stddev(residual[correl_indicies])*A[2]*SQRT(2*!DPI), $
                                            Mult_and_add[0], smooth_and_shift[0], smooth_and_shift[1], Use_Scatter[i], float(sxpar(header, 'T_PSHADO')), float(sxpar(sig_header, 'EXPTIME')) / 60.}
        Brightness_array[i]      = A[0]*A[2]*SQRT(2*!DPI)                               ; Area under the Gaussian
        err_Brightness_array[i]  = stddev(residual[correl_indicies])*A[2]*SQRT(2*!DPI)  ; use the std dev of the residual and the width of the line
        T_P_Shadow_array[i]      = sxpar(header, 'T_PSHADO')
        T_U_Shadow_array[i]      = sxpar(header, 'T_USHADO')
        EXPTime_array[i]         = sxpar(sig_header, 'EXPTIME') / 60.                   ; for some reason, this keyword is bad in the regular headers, so use sig_header as a workaround
        DopplerShift_array1[i]   = O3_Io_frame
    endfor

      ; Plot the residual and Gaussian fits
        cgplot, spec, WL, psym = 0, Ytitle = 'Residual (KR / '+cgsymbol('Angstrom')+')', xtitle = 'Wavelength ('+cgsymbol('Angstrom')+')', $
          yr = YR_resid, xr = XR, /nodata, pos = pos[*,1], /noerase
        cgtext, O3_Io_frame - 0.7, 5, "Io's Doppler Shift", charsize = 1.4, alignment = 0.5
        cgtext, O3 + .5, 5, "Telluric [O I]", charsize = 1.4, alignment = 0.5
        for i = 1, n_elements(Eclipse_files)-1 do begin
          cgplot, WL[plot_WLs], plot_residual_array[*, i], /OVERPLOT, COLOR = timeColors[i], psym=10
          cgplot, WL_array[*, i], LSF_Fit_array[*, i], /OVERPLOT, COLOR = timeColors[i]
          cgplot, [DopplerShift_array1[i], DopplerShift_array1[i]], [YR_resid[1], 5], COLOR = timeColors[i], /overplot
          cgplot, [O3, O3], [YR_resid[1], 3.], linestyle = 1, /overplot
        endfor
    cgps_close

    ; *****************************************************************O I 6364A LIGHT CURVE***********************************************************************************

    pos =  [.12,.17,.95,.9]
    cgPS_Open, filename = Reduced_Dir+'O6364_lightcurve.eps', /ENCAPSULATED, xsize = 6., ysize = 4.
    !P.font=1
    device, SET_FONT = 'Helvetica Bold', /TT_FONT

    loadct, 0
    yr = [0., 5.]
    cgplot, T_P_Shadow_array, Brightness_array, psym = 15, Ytitle = 'Disk-Averaged Brightness (KR)', xtitle = 'Minutes After Ingress', $
      title = 'Io''s Oxygen 6364A Line Response to Eclipse', $
      ERR_YLOW = ERR_Brightness_array, ERR_YHigh = ERR_Brightness_array, ERR_XLOW = EXPTime_array/2., ERR_XHigh = EXPTime_array/2., yr = yr, xr = time_range, /nodata, pos = pos

    x = findgen(Umbra_ET - PenUmbra_ET) / 60.
    colors = reverse(BytScl(x, MIN=min(x), MAX=2.*max(x)))+128
    FOR j=0,n_elements(x)-2 DO BEGIN
      xpoly = [x[j],         x[j], x[j+1],       x[j+1],         x[j]]
      ypoly = [  0., !Y.CRange[1], !Y.CRange[1], 0., 0.]
      cgColorFill, xpoly, ypoly, Color=colors[j]
    ENDFOR

    xpoly = [max(x), max(x), !X.CRange[1],  !X.CRange[1],  max(x)]
    ypoly = [  0., !Y.CRange[1], !Y.CRange[1], 0., 0.]
    cgColorFill, xpoly, ypoly, color = 'charcoal', /LINE_FILL, orientation = 45., thick = .1
    cgColorFill, xpoly, ypoly, color = 'charcoal', /LINE_FILL, orientation = -45., thick = .1
    cgplot, T_P_Shadow_array, Brightness_array, psym = 4, /noerase, symsize = 2, yr = yr, xr = time_range, pos = pos, $
      ERR_YLOW = ERR_Brightness_array, ERR_YHigh = ERR_Brightness_array, ERR_XLOW = EXPTime_array/2., ERR_XHigh = EXPTime_array/2.
    cgLoadCT, 33, NColors=8, /reverse

    if ingress then x = findgen(Umbra_ET - PenUmbra_ET) / 60.
    FOR j=0,n_elements(x)-2 DO BEGIN
      xpoly = [x[j],         x[j], x[j+1],       x[j+1],         x[j]]
      ypoly = [  0., !Y.CRange[1], !Y.CRange[1], 0., 0.]
      cgColorFill, xpoly, ypoly, Color=colors[j]
    ENDFOR
    if ingress then begin
      xpoly = [max(x),     max(x), !X.CRange[1],  !X.CRange[1],  max(x)]
      ypoly = [  0., !Y.CRange[1], !Y.CRange[1], 0., 0.]
      cgColorFill, xpoly, ypoly, color = 'charcoal', /LINE_FILL, orientation = 45., thick = .1
      cgColorFill, xpoly, ypoly, color = 'charcoal', /LINE_FILL, orientation = -45., thick = .1
    endif
    cgtext, 2.5, !Y.CRange[1]/4., 'Penumbral Eclipse', orientation = 90., color = 'white'
    cgaxis, yaxis = 0, yr = yr, /ystyle                        ; repair the axis damage that the Penumbra did
    cgaxis, xaxis = 0, xr = time_range                         ; repair axis damage
    cgaxis, xaxis = 1, xr = time_range, xtickformat = '(A1)'   ; repair axis damage

    Case date of
      'UT180320': AL_legend, ['O 6364 in Ingress 180320'], Psym = 4, /right, charsize = 1.5, /clear
      'UT190812': AL_legend, ['O 6364 in Umbra 190812'], Psym = 4, /right, charsize = 1.5, /clear
    endcase
    cgps_close
    save, O3_fit_params, filename = Reduced_Dir+'O3_fit_params.sav'
endif
if Part eq 1.02 then begin ; Make waterfall plots, O 5577 requires telluric corrections

    ;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ;--------------------------------------------------------Get the O5577 Time Series and Scatter Fit-------------------------------------------------------------------------------
    ;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ; METHOD:
    ;   Loop though all possible spectra that could be fit to jovian scatter and *MANUALLY* choose what works best. Avoid Io and telluric airglow in the fitting.
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
    ;  
    ; NOTES:
    ;   It's critical to inspect with plot, WL[correl_indicies], scatter_fit[correl_indicies], uncomment this line to inspect and adjust things
    ;   With this inspection, it is useful to adjust *include_WLs*, and *AG_ind* to cover regions of continuum as close as possible to Io's airglow
    ;   Parameter information is not included in the Gaussian fitting to the line spread function in the residual. I couldn't figgure out why it failed with status = 0

    ;!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ;!
    Unexpected_Doppler_Shift = -0.06   ; The seperation between Telluric and Io lines does not seem to be what it should! 
    ;                                    This term accounts for this in the indicies that are fit.
    ;                                    Allow the Gaussian centroid to be a free parameter in the line fitting
    ;!
    ;!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ; Fit various spectra to the Jovian Scatter, *use* whichever one produces the lowest residual. Note this method will never work if there's much telluric absorption.
    Case date of
      'UT180320': begin
        Fit_Me_To_The_Scatter    = ['R_per_A_Jupiter_Center.ec.fits', $
                                  'R_per_A_Jovian_Scatter.0001.ec.fits', $
                                  'R_per_A_Jovian_Scatter.0002.ec.fits', $
                                  'R_per_A_Jovian_Scatter.0003.ec.fits' ]
                  end
      'UT190812': begin
        Fit_Me_To_The_Scatter    = ['R_per_A_Jupiter_Center.ec.fits']
                  end
    endcase

    ; setup the plot axis
      spec      = MRDFITS(reduced_dir+'R_per_A_'+Eclipse_files[1]+'.ec.fits', 0, header, /fscale, /unsigned, /silent)
      WL        = sxpar(header, 'CRVAL1')+findgen(N_elements(spec))*sxpar(header, 'Cdelt1')
      bandwidth = 2.                                                                              ; half with of the plot's x axis in Angstroms

    Case date of
      'UT180320': begin
        YR   = [47, 102]
        XR   = [O1 + O1 * sxpar(header, 'IO_DOPPL') / cspice_clight() - Bandwidth, O1 + O1 * sxpar(header, 'IO_DOPPL') / cspice_clight() + Bandwidth]
      end
      'UT190812': begin
        yr   = [5.5e1, 2e2]
        YR_resid = [-2, 20]
        XR       = [O1 + O1 * sxpar(header, 'IO_DOPPL') / cspice_clight() - Bandwidth, O1 + O1 * sxpar(header, 'IO_DOPPL') / cspice_clight() + Bandwidth]
      end
    endcase

    ; Define arrays
      include_WLs            = where( (wl gt xr[0]) and (wl lt xr[1]), /NULL )                      ; wavelengths to use in the scattered light fitting
      plot_WLs               = where( (wl gt xr[0]) and (wl lt xr[1]), /NULL )                      ; wavelengths to actually plot
      residual_array         = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
      plot_residual_array    = fltarr(N_elements(plot_WLs), n_elements(Eclipse_files))
      LSF_Fit_array          = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
      WL_array               = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
      Brightness_array       = fltarr(n_elements(Eclipse_files))
      err_Brightness_array   = fltarr(n_elements(Eclipse_files))
      T_P_Shadow_array       = fltarr(n_elements(Eclipse_files))
      T_U_Shadow_array       = fltarr(n_elements(Eclipse_files))
      EXPTime_array          = fltarr(n_elements(Eclipse_files))
      DopplerShift_array1    = fltarr(n_elements(Eclipse_files))                                    ; For the O1 line
      Gof                    = 0.                                                                   ; Goodness of fit (co-addded mean residual)
      GOF_array              = fltarr(N_elements(Fit_Me_To_The_Scatter))
      Use_Scatter            = strarr(N_elements(Eclipse_files))
      Possible_WL_offset     = fltarr(N_elements(Eclipse_files), N_elements(Fit_Me_To_The_Scatter)) ; Save any offsets from Io's expected rest wavelength, then correct any error in the wavelength solution.
      Possible_line_width    = fltarr(N_elements(Eclipse_files), N_elements(Fit_Me_To_The_Scatter)) ; Save the linewidths from all type of scatter subtraction, then average and peg this parameter.
      O1_Fit_params          = REPLICATE({params, brightness:0.0, err_brightness:0.0, mult:0.0, smooth:0.0, shift:0.0, background:'', T_P_shadow:0.0, exptime:0.0}, n_elements(Eclipse_files))

    for i = 0, n_elements(Eclipse_files)-1 do begin
      if i eq 0 then continue ; skip penumbra
      spec        = MRDFITS(reduced_dir+'R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, header, /fscale, /unsigned, /silent)
      spec_err    = MRDFITS(reduced_dir+'sig_R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, sig_header, /fscale, /unsigned, /silent)
      WL          = sxpar(header, 'CRVAL1')+findgen(N_elements(spec))*sxpar(header, 'Cdelt1')

      for h = 0, n_elements(Fit_Me_To_The_Scatter)-1 do begin
        ; Which Jovian scatter are we using?
          Jup_Cal_Tell_Corr = MRDFITS(reduced_dir+Fit_Me_To_The_Scatter[h], 0, Jupiter_Scatter_header, /fscale, /unsigned, /silent)

        ; Define fitting indicies we need to avoid because of airglow itself.
          Ios_Airglow    = [Io_airglow + Io_airglow * sxpar(header, 'IO_DOPPL') / cspice_clight()]                                      ; wavelength locations of all airglow lines, Ignore telluric K
          ARCES_2sigma_bandwidth = 2. * (mean(xr)/31500.) / 2.3548 ; 2. * FWHM in Angstroms / 2sqrt(2ln2) = 95% of flux enclosed within 2 sigma --> Appropriate for Bright Io emission
          ARCES_1sigma_bandwidth = 1. * (mean(xr)/31500.) / 2.3548 ; 1. * FWHM in Angstroms / 2sqrt(2ln2) = 68% of flux enclosed within 1 sigma --> Appropriate for weak telluric airglow
          Io_AG_Ind = [] & Telluric_AG_Ind = []
          for j = 0, N_elements(Ios_Airglow)-1 do Io_AG_ind = [Io_AG_Ind, where(abs(Ios_Airglow[j] - WL) lt ARCES_2sigma_bandwidth, /NULL)] ;  indicies within 2 sigma of an Io airglow line
          for j = 0, N_elements(Telluric_Airglow)-1 do Telluric_AG_ind = [Telluric_AG_Ind, where(abs(Telluric_Airglow[j] - WL) lt ARCES_1sigma_bandwidth, /NULL)] ; indicies within 1 sigma of any Telluric airglow
          AG_ind = [Io_AG_Ind, Telluric_AG_Ind]
          AG_ind = AG_ind[UNIQ(AG_ind, SORT(AG_ind))]

        ; Fit the Jovian scattered light, use Amoeba/Correlation to match Jupiter in terms of smoothing and shifting, then multiply and add until matches the scattered background
          correl_indicies  = cgSetDifference(include_Wls, AG_Ind)
          Mult_and_add     = mpfitfun('MX_Plus_B', Jup_Cal_Tell_Corr[correl_indicies], Spec[correl_indicies], Spec_err[correl_indicies], $
            [mean(Spec[correl_indicies]/Jup_Cal_Tell_Corr[correl_indicies]), 0.], /NaN, status = Scatter_fit_status, PERROR = err_a, /quiet)
          jup              = Mult_and_add[0]*Jup_Cal_Tell_Corr + Mult_and_add[1]
  
          ; Iterate Amoeba until it gives an answer
          Ftol = 1.e-4 &  Xtol = 1.e-4 & count = 0                                          ; Define initial tolerances for Amoeba
          smooth_and_shift = [2., 0.] & scale  = [2., 2.]                                   ; Define intial guess and simplex for Amoeba
          REPEAT BEGIN
            trial_smooth_and_shift = AMOEBAX(Ftol, xtol, function_name='shift_smooth', SCALE = Scale, P0 = smooth_and_shift, FUNCTION_VALUE = fval)
            if N_elements (N_elements(trial_smooth_and_shift) lt 2) then begin              ; If Amoeba crashes
              xtol = xtol * 1.5 & ftol = ftol * 1.5                                         ; Grow the tolerances
              smooth_and_shift = smooth_and_shift + [RANDOMN(seed), RANDOMN(seed)]          ; Change the guess
              SCALE            = SCALE + [RANDOMN(seed), RANDOMN(seed)]                     ; Change the simplex
              count = count+1
            endif
          ENDREP UNTIL ( (N_elements(trial_smooth_and_shift) eq 2) or (count gt 500) )
          smooth_and_shift = AMOEBAX(Ftol, xtol, function_name='shift_smooth', SCALE = [2., 2.], P0 = smooth_and_shift, FUNCTION_VALUE = fval)
          if N_elements(smooth_and_shift) eq 1 then stop ; check for bugs
  
          smoothed_shifted = interpolate(gauss_smooth(jup, smooth_and_shift[0]), findgen(n_elements(Jup)) - smooth_and_shift[1])
          Mult_and_add     = mpfitfun('MX_Plus_B', smoothed_shifted[correl_indicies], Spec[correl_indicies], Spec_err[correl_indicies], $
            [mean(Spec[correl_indicies]/smoothed_shifted[correl_indicies]), 0.], /NaN, status = Scatter_fit_status, PERROR = err_a, /quiet)
          scatter_fit      = Mult_and_add[0]*smoothed_shifted + Mult_and_add[1]
          residual         = Spec - scatter_fit
          GOF_array[h]     = GOF + median(abs(residual[correl_indicies]))

        ; Get the linewidth and the determine any systematic wavelength shift from Io's rest
        ; Ultimately we will peg both of these parameters.
          O1_Io_frame               = O1 + O1 * sxpar(header, 'IO_DOPPL') / cspice_clight() + Unexpected_Doppler_Shift
          LSF_fitting_ind1          = cgSetDifference(where( abs(wl - O1_Io_frame) lt 0.3, /NULL), Telluric_AG_Ind)  ; fit region within +/- 0.3A of line center, excluding Telluric airglow
          initial_guess             = [2.e3, O1_Io_frame, 0.077]  ; rough values are fine here [height, wavelength, linewidth_sigma]
          fa                        = {x:double(WL[LSF_fitting_ind1]), y:double(residual[LSF_fitting_ind1]), err:double( sqrt(abs( residual[LSF_fitting_ind1] )) )}
          a                         = mpfit('Gaussian_for_MPFIT', initial_guess, PERROR = err_a, funct=fa, maxiter=50, STATUS = Did_it_work, /Quiet, NPEGGED = NPEGGED)
          Possible_WL_offset[i,h]   = a[1] - O1_Io_frame
          Possible_Line_Width[i,h]  = a[2]
      endfor
      best_GoF = min(GOF_array, best_fit)
      Use_Scatter[i] = Fit_Me_To_The_Scatter[best_fit]
      print, 'Frame: ', Eclipse_files[i], ' is best paired with the Jupiter scatter from: ', Use_Scatter[i], ', Goodness of fit = ', best_GoF
    endfor
        Case date of
      'UT180320': begin
          best_scatter_spectrum = 'R_per_A_Jovian_Scatter.0003.ec.fits' ; Manually input desired scatting spectrum.
      end
      'UT190812': begin
          best_scatter_spectrum = 'R_per_A_Jupiter_Center.ec.fits' ; Manually input desired scatting spectrum.
      end
    endcase
    print, 'Upon inspection, best choice seems to be: ', best_scatter_spectrum

    ; -------------------------------------------Waterfall Plot Postscript Io and Fit Jovian scatter--------------------------------------------------------------------------------------------------------------

    ; Okay, now that we've established which Jupiter scatter spectrum best optimzes the residual, we can actually use it for scattered light subtraction

    ; get color versus ingress time & setup plot positions
    timeColors = BytScl(findgen(11), Min=0, Max=11, Top=11)
    Color      = timeColors[0]
    cgLoadCT, 33, NColors=8, /reverse
    pos        = cgLayout([1,2], OXMargin=[10,3], OYMargin=[8,5], YGap=0)
    Pos[1,0]   = Pos[1,0]*.7 & Pos[3,1] = Pos[3,1]*.7

    cgPS_Open, filename = Reduced_Dir+'O5577_scatterfit.eps', /ENCAPSULATED, xsize = 7.5, ysize = 5.
    !P.font = 1
    device, SET_FONT = 'Helvetica Bold', /TT_FONT
    !p.charsize = 1.5

    cgplot, WL, spec/1.e3, psym = 0, Ytitle = 'KR / '+cgsymbol('Angstrom'), $             ; Rayleighs to KiloRayleighs
      title = 'Io''s [O I] 5577'+cgsymbol('Angstrom')+' Response to Eclipse',  $
      xr = xr, /nodata, pos = pos[*,0], xtickformat = '(A1)', yr = yr

    for i = 0, n_elements(Eclipse_files)-1 do begin
      if i eq 0 then continue ; skip penumbra
      spec         = MRDFITS(reduced_dir+'R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, header, /fscale, /unsigned, /silent)
      spec_err     = MRDFITS(reduced_dir+'sig_R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, sig_header, /fscale, /unsigned, /silent)
      spec         = spec / 1.e3       ; to KR
      spec_err     = spec_err / 1.e3   ; to KR
      WL           = sxpar(header, 'CRVAL1') + findgen(N_elements(spec))*sxpar(header, 'Cdelt1')  
      O1_Io_frame  = O1 + O1 * sxpar(header, 'IO_DOPPL') / cspice_clight()               
      junk         = min( abs(WL - O1_Io_frame), O1_Io_Ind)

      ; Which Jovian scatter are we using?
      Jup_Cal_Tell_Corr = MRDFITS(reduced_dir+best_scatter_spectrum, 0, Jup_Scatter_header, /fscale, /unsigned, /silent)

      ; Define fitting indicies we need to avoid because of airglow itself.
      Ios_Airglow    = [Io_airglow + Io_airglow * sxpar(header, 'IO_DOPPL') / cspice_clight()] + Unexpected_Doppler_Shift     ; wavelength locations of Io's 5577 lines seems blue shifted
      ARCES_2sigma_bandwidth = 2. * (mean(xr)/31500.) / 2.3548 ; 2. * FWHM in Angstroms / 2sqrt(2ln2) = 95% of flux enclosed within 2 sigma --> Appropriate for faint Io emission
      ARCES_1sigma_bandwidth = 1. * (mean(xr)/31500.) / 2.3548 ; 1. * FWHM in Angstroms / 2sqrt(2ln2) = 68% of flux enclosed within 1 sigma --> Appropriate for Bright 5577 telluric airglow
      Io_AG_Ind = [] & Telluric_AG_Ind = []
      for j = 0, N_elements(Ios_Airglow)-1 do Io_AG_ind = [Io_AG_Ind, where(abs(Ios_Airglow[j] - WL) lt ARCES_1sigma_bandwidth, /NULL)] ;  indicies within 1.5 sigma of an Io airglow line
      for j = 0, N_elements(Telluric_Airglow)-1 do Telluric_AG_ind = [Telluric_AG_Ind, where(abs(Telluric_Airglow[j] - WL) lt ARCES_2sigma_bandwidth, /NULL)] ; indicies within 2 sigma of any Telluric airglow
      AG_ind = [Io_AG_Ind, Telluric_AG_Ind]
      AG_ind = AG_ind[UNIQ(AG_ind, SORT(AG_ind))]

      ; Fit the Jovian scattered light, use Amoeba/Correlation to match Jupiter in terms of smoothing and shifting, then multiply and add until matches the scattered background
        correl_indicies  = cgSetDifference(include_Wls, AG_Ind)
        Mult_and_add     = mpfitfun('MX_Plus_B', Jup_Cal_Tell_Corr[correl_indicies], Spec[correl_indicies], Spec_err[correl_indicies], $
                                    [mean(Spec[correl_indicies]/Jup_Cal_Tell_Corr[correl_indicies]), 0.], /NaN, status = Scatter_fit_status, PERROR = err_a, /quiet)
        jup              = Mult_and_add[0]*Jup_Cal_Tell_Corr + Mult_and_add[1]
  
        ; Iterate Amoeba until it gives an answer
        Ftol = 1.e-4 &  Xtol = 1.e-4 & count = 0                                          ; Define initial tolerances for Amoeba
        smooth_and_shift = [2., 0.] & scale  = [2., 2.]                                   ; Define intial guess and simplex for Amoeba
        REPEAT BEGIN
          trial_smooth_and_shift = AMOEBAX(Ftol, xtol, function_name='shift_smooth', SCALE = Scale, P0 = smooth_and_shift, FUNCTION_VALUE = fval)
          if N_elements (N_elements(trial_smooth_and_shift) lt 2) then begin              ; If Amoeba crashes
            xtol = xtol * 1.5 & ftol = ftol * 1.5                                         ; Grow the tolerances
            smooth_and_shift = smooth_and_shift + [RANDOMN(seed), RANDOMN(seed)]          ; Change the guess
            SCALE            = SCALE + [RANDOMN(seed), RANDOMN(seed)]                     ; Change the simplex
            count = count+1
          endif
        ENDREP UNTIL ( (N_elements(trial_smooth_and_shift) eq 2) or (count gt 500) )
        smooth_and_shift = AMOEBAX(Ftol, xtol, function_name='shift_smooth', SCALE = [2., 2.], P0 = smooth_and_shift, FUNCTION_VALUE = fval)
        if N_elements(smooth_and_shift) eq 1 then stop ; check for bugs
  
        smoothed_shifted = interpolate(gauss_smooth(jup, smooth_and_shift[0]), findgen(n_elements(Jup)) - smooth_and_shift[1])
        Mult_and_add     = mpfitfun('MX_Plus_B', smoothed_shifted[correl_indicies], Spec[correl_indicies], Spec_err[correl_indicies], [mean(Spec[correl_indicies]/smoothed_shifted[correl_indicies]), 0.], $
                                    /NaN, status = Scatter_fit_status, PERROR = err_a, /quiet, weights = 1./abs(O1_Io_Ind - correl_indicies)^(1.5))
        scatter_fit      = Mult_and_add[0]*smoothed_shifted + Mult_and_add[1]
        residual         = Spec - scatter_fit
        residual_err     = make_array( N_elements(residual), value = stddev(residual[correl_indicies]) )

      ; plot the spectrum, fit and indices used for the fit (if we're in debug mode)
        ;cgplot, WL[correl_indicies], scatter_fit[correl_indicies], psym=14, /overplot
        cgplot, WL, scatter_fit, color = timeColors[i], linestyle = 1, thick = 7, /overplot
        cgplot, WL, spec, color = timeColors[i], thick = 5, /overplot

      ; Now do the line fitting to the residual. First define Line Spread Function (LSF) Gaussian fit parameter info
        parinfo               = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0]}, 3)
        parinfo[1].limited    = [1, 1]
        parinfo[1].limits     = double(O1_Io_frame) + Unexpected_Doppler_shift + [-.07, 0.07]                 ; Give some room for the lines cetner to drift 
        parinfo[2].limited    = [1, 1]
        parinfo[2].limits     = median(Possible_line_width[1:*,*]) + [-0.005, 0.005]                          ; Give some room for line widths to drift as well
        LSF_fitting_ind1      = cgSetDifference(where( abs(wl - O1_Io_frame) lt 0.3, /NULL), Telluric_AG_Ind) ; fit region within +/- 0.3A of line center, excluding Telluric airglow

      ; Gaussian fit
        initial_guess = [2.D, mean(parinfo[1].limits), mean(parinfo[2].limits)]
        fa            = {x:double(WL[LSF_fitting_ind1]), y:double(residual[LSF_fitting_ind1]), err:double( sqrt(abs( spec[LSF_fitting_ind1] )) )}
        a             = mpfit('Gaussian_for_MPFIT', initial_guess, PERROR = err_a, funct=fa, parinfo = parinfo, STATUS = Did_it_work, /Quiet, NPEGGED = NPEGGED)
  
        ; write the results
        WL_array[*, i]           = WL[include_WLs]
        plot_residual_array[*, i]= residual[plot_WLs]
        LSF_Fit_array[*, i]      = gaussian(WL[include_WLs], a)                         ; LSF_Fit
        O1_fit_params[i]         = {params, A[0]*A[2]*SQRT(2*!DPI), stddev(residual[correl_indicies])*A[2]*SQRT(2*!DPI), $
                                            Mult_and_add[0], smooth_and_shift[0], smooth_and_shift[1], Use_Scatter[i], float(sxpar(header, 'T_PSHADO')), float(sxpar(sig_header, 'EXPTIME')) / 60.}
        Brightness_array[i]      = A[0]*A[2]*SQRT(2*!DPI)                               ; Area under the Gaussian
        err_Brightness_array[i]  = stddev(residual[correl_indicies])*A[2]*SQRT(2*!DPI)  ; use the std dev of the residual and the width of the line
        T_P_Shadow_array[i]      = sxpar(header, 'T_PSHADO')
        T_U_Shadow_array[i]      = sxpar(header, 'T_USHADO')
        EXPTime_array[i]         = sxpar(sig_header, 'EXPTIME') / 60.                   ; for some reason, this keyword is bad in the regular headers, so use sig_header as a workaround
        DopplerShift_array1[i]   = O1_Io_frame
    endfor

      ; Plot the residual and Gaussian fits
        cgplot, spec, WL, psym = 0, Ytitle = 'Residual (KR / '+cgsymbol('Angstrom')+')', xtitle = 'Wavelength ('+cgsymbol('Angstrom')+')', $
          yr = [-2., 14], xr = XR, /nodata, pos = pos[*,1], /noerase
        cgtext, O1_Io_frame - 0.7, 5, "Io's Doppler Shift", charsize = 1.4, alignment = 0.5
        cgtext, O1 + .5, 5, "Telluric [O I]", charsize = 1.4, alignment = 0.5
        for i = 1, n_elements(Eclipse_files)-1 do begin
          cgplot, WL[plot_WLs], plot_residual_array[*, i], /OVERPLOT, COLOR = timeColors[i];, psym=10
          cgplot, WL_array[*, i], LSF_Fit_array[*, i], /OVERPLOT, COLOR = timeColors[i]
          cgplot, [DopplerShift_array1[i], DopplerShift_array1[i]], [14, 3], COLOR = timeColors[i], /overplot
          cgplot, [O1, O1], [-2, 14], linestyle = 1, /overplot
        endfor
    cgps_close

    ; *****************************************************************O I 5577A LIGHT CURVE***********************************************************************************

    pos =  [.12,.17,.95,.9]
    cgPS_Open, filename = Reduced_Dir+'O5577_lightcurve.eps', /ENCAPSULATED, xsize = 6., ysize = 4.
    !P.font=1
    device, SET_FONT = 'Helvetica Bold', /TT_FONT

    loadct, 0
    yr = [0., 1.4]
    cgplot, T_P_Shadow_array, Brightness_array, psym = 15, Ytitle = 'Disk-Averaged Brightness (KR)', xtitle = 'Minutes After Ingress', $
      title = 'Io''s Oxygen 5577A Line Response to Eclipse', $
      ERR_YLOW = ERR_Brightness_array, ERR_YHigh = ERR_Brightness_array, ERR_XLOW = EXPTime_array/2., ERR_XHigh = EXPTime_array/2., yr = yr, xr = time_range, /nodata, pos = pos

    x = findgen(Umbra_ET - PenUmbra_ET) / 60.
    colors = reverse(BytScl(x, MIN=min(x), MAX=2.*max(x)))+128
    FOR j=0,n_elements(x)-2 DO BEGIN
      xpoly = [x[j],         x[j], x[j+1],       x[j+1],         x[j]]
      ypoly = [  0., !Y.CRange[1], !Y.CRange[1], 0., 0.]
      cgColorFill, xpoly, ypoly, Color=colors[j]
    ENDFOR

    xpoly = [max(x), max(x), !X.CRange[1],  !X.CRange[1],  max(x)]
    ypoly = [  0., !Y.CRange[1], !Y.CRange[1], 0., 0.]
    cgColorFill, xpoly, ypoly, color = 'charcoal', /LINE_FILL, orientation = 45., thick = .1
    cgColorFill, xpoly, ypoly, color = 'charcoal', /LINE_FILL, orientation = -45., thick = .1

    cgplot, T_P_Shadow_array, Brightness_array, psym = 4, /noerase, symsize = 2, yr = yr, xr = time_range, pos = pos, $
      ERR_YLOW = ERR_Brightness_array, ERR_YHigh = ERR_Brightness_array, ERR_XLOW = EXPTime_array/2., ERR_XHigh = EXPTime_array/2.
    cgLoadCT, 33, NColors=8, /reverse

    if ingress then x = findgen(Umbra_ET - PenUmbra_ET) / 60.
    FOR j=0,n_elements(x)-2 DO BEGIN
      xpoly = [x[j],         x[j], x[j+1],       x[j+1],         x[j]]
      ypoly = [  0., !Y.CRange[1], !Y.CRange[1], 0., 0.]
      cgColorFill, xpoly, ypoly, Color=colors[j]
    ENDFOR
    if ingress then begin
      xpoly = [max(x),     max(x), !X.CRange[1],  !X.CRange[1],  max(x)]
      ypoly = [  0., !Y.CRange[1], !Y.CRange[1], 0., 0.]
      cgColorFill, xpoly, ypoly, color = 'charcoal', /LINE_FILL, orientation = 45., thick = .1
      cgColorFill, xpoly, ypoly, color = 'charcoal', /LINE_FILL, orientation = -45., thick = .1
    endif
    cgtext, 2.5, !Y.CRange[1]/4., 'Penumbral Eclipse', orientation = 90., color = 'white'
    cgaxis, yaxis = 0, yr = yr, /ystyle                        ; repair the axis damage that the Penumbra did
    cgaxis, xaxis = 0, xr = time_range                         ; repair axis damage
    cgaxis, xaxis = 1, xr = time_range, xtickformat = '(A1)'   ; repair axis damage

    Case date of
      'UT180320': AL_legend, ['O 5577 in Ingress 180320'], Psym = 4, /right, charsize = 1.5, /clear
      'UT190812': AL_legend, ['O 5577 in Umbra 190812'], Psym = 4, /right, charsize = 1.5, /clear
    endcase
    cgps_close
    save, O1_fit_params, filename = Reduced_Dir+'O1_fit_params.sav'
    print, 'Can it be true? O 5577A seems well blue-shifted from Io''s disk? 
  endif
  if Part eq 1.03 then begin ; Make waterfall plots, O 7774 requires telluric corrections

    ;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ;--------------------------------------------------------Get the O7774 Time Series and Scatter Fit-------------------------------------------------------------------------------
    ;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ; METHOD:
    ;   Loop though all possible spectra that could be fit to jovian scatter and *MANUALLY* choose what works best. Avoid Io and telluric airglow in the fitting.
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
    Case date of
      'UT180320': begin
        Fit_Me_To_The_Scatter    = ['R_per_A_Jupiter_Center.ec.fits', $
                                  'R_per_A_Jovian_Scatter.0001.ec.fits', $
                                  'R_per_A_Jovian_Scatter.0002.ec.fits', $
                                  'R_per_A_Jovian_Scatter.0003.ec.fits' ]
                  end
      'UT190812': begin
        Fit_Me_To_The_Scatter    = ['R_per_A_Jupiter_Center.ec.fits']
                  end
    endcase

    ; setup the plot axis
    spec      = MRDFITS(reduced_dir+'R_per_A_'+Eclipse_files[1]+'.ec.fits', 0, header, /fscale, /unsigned, /silent)
    WL        = sxpar(header, 'CRVAL1')+findgen(N_elements(spec))*sxpar(header, 'Cdelt1')
    bandwidth = 6.                                                                               ; half with of the plot's x axis in Angstroms

    Case date of
      'UT180320': begin
        YR   = [21, 41.6]
        XR   = [O5 + O5 * sxpar(header, 'IO_DOPPL') / cspice_clight() - Bandwidth, O5 + O5 * sxpar(header, 'IO_DOPPL') / cspice_clight() + Bandwidth] - 0.5
      end
      'UT190812': begin
        yr   = [5.5e1, 2e2]
        YR_resid = [-2, 20]
        XR   = [O5 + O5 * sxpar(header, 'IO_DOPPL') / cspice_clight() - Bandwidth, O5 + O5 * sxpar(header, 'IO_DOPPL') / cspice_clight() + Bandwidth] - 0.5
      end
    endcase

    ; Define arrays
      include_WLs            = where( (wl gt xr[0]) and (wl lt xr[1]), /NULL )                      ; wavelengths to use in the scattered light fitting
      plot_WLs               = where( (wl gt xr[0]) and (wl lt xr[1]), /NULL )                      ; wavelengths to actually plot
      residual_array         = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
      plot_residual_array    = fltarr(N_elements(plot_WLs), n_elements(Eclipse_files))
      LSF_Fit_array          = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
      WL_array               = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
      Brightness_array       = fltarr(n_elements(Eclipse_files))
      err_Brightness_array   = fltarr(n_elements(Eclipse_files))
      T_P_Shadow_array       = fltarr(n_elements(Eclipse_files))
      T_U_Shadow_array       = fltarr(n_elements(Eclipse_files))
      EXPTime_array          = fltarr(n_elements(Eclipse_files))
      DopplerShift_array1    = fltarr(n_elements(Eclipse_files))                                    ; For the first line in the triplet
      DopplerShift_array2    = fltarr(n_elements(Eclipse_files))
      DopplerShift_array3    = fltarr(n_elements(Eclipse_files))                                      
      Gof                    = 0.                                                                   ; Goodness of fit (co-addded mean residual)
      GOF_array              = fltarr(N_elements(Fit_Me_To_The_Scatter))
      Use_Scatter            = strarr(N_elements(Eclipse_files))
      Possible_WL_offset     = fltarr(N_elements(Eclipse_files), N_elements(Fit_Me_To_The_Scatter)) ; Save any offsets from Io's expected rest wavelength, then correct any error in the wavelength solution.
      Possible_line_width    = fltarr(N_elements(Eclipse_files), N_elements(Fit_Me_To_The_Scatter)) ; Save the linewidths from all type of scatter subtraction, then average and peg this parameter.

    for i = 0, n_elements(Eclipse_files)-1 do begin
      if i eq 0 then continue ; skip penumbra
      spec        = MRDFITS(reduced_dir+'R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, header, /fscale, /unsigned, /silent)
      spec_err    = MRDFITS(reduced_dir+'sig_R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, sig_header, /fscale, /unsigned, /silent)
      WL          = sxpar(header, 'CRVAL1')+findgen(N_elements(spec))*sxpar(header, 'Cdelt1')

      for h = 0, n_elements(Fit_Me_To_The_Scatter)-1 do begin
        ; Which Jovian scatter are we using?
          Jup_Cal_Tell_Corr = MRDFITS(reduced_dir+Fit_Me_To_The_Scatter[h], 0, Jupiter_Scatter_header, /fscale, /unsigned, /silent)

        ; Define fitting indicies we need to avoid because of airglow itself.
          Ios_Airglow    = [Io_airglow + Io_airglow * sxpar(header, 'IO_DOPPL') / cspice_clight()]                                      ; wavelength locations of all airglow lines, Ignore telluric K
          ARCES_1sigma_bandwidth = (mean(xr)/31500.) / 2.3548 ; 1. * FWHM in Angstroms / 2sqrt(2ln2) = 68% of flux enclosed within 1 sigma --> Appropriate for weak telluric airglow
          Io_AG_Ind = [] & Telluric_AG_Ind = []
          for j = 0, N_elements(Ios_Airglow)-1 do Io_AG_ind = [Io_AG_Ind, where(abs(Ios_Airglow[j] - WL) lt 2.5*ARCES_1sigma_bandwidth, /NULL)] ;  indicies within 2.5 sigma of an Io airglow line
          for j = 0, N_elements(Telluric_Airglow)-1 do Telluric_AG_ind = [Telluric_AG_Ind, where(abs(Telluric_Airglow[j] - WL) lt ARCES_1sigma_bandwidth, /NULL)] ; indicies within 1 sigma of any Telluric airglow
          AG_ind = [Io_AG_Ind, Telluric_AG_Ind]
          AG_ind = AG_ind[UNIQ(AG_ind, SORT(AG_ind))]

        ; Fit the Jovian scattered light, use Amoeba/Correlation to match Jupiter in terms of smoothing and shifting, then multiply and add until matches the scattered background
          correl_indicies  = cgSetDifference(include_Wls, AG_Ind)
          Mult_and_add     = mpfitfun('MX_Plus_B', Jup_Cal_Tell_Corr[correl_indicies], Spec[correl_indicies], Spec_err[correl_indicies], $
                                      [mean(Spec[correl_indicies]/Jup_Cal_Tell_Corr[correl_indicies]), 0.], /NaN, status = Scatter_fit_status, PERROR = err_a, /quiet)
          jup              = Mult_and_add[0]*Jup_Cal_Tell_Corr + Mult_and_add[1]

          ; Iterate Amoeba until it gives an answer
          Ftol = 1.e-4 &  Xtol = 1.e-4 & count = 0                                          ; Define initial tolerances for Amoeba
          smooth_and_shift = [2., 0.] & scale  = [2., 2.]                                   ; Define intial guess and simplex for Amoeba
          REPEAT BEGIN
            trial_smooth_and_shift = AMOEBAX(Ftol, xtol, function_name='shift_smooth', SCALE = Scale, P0 = smooth_and_shift, FUNCTION_VALUE = fval)
            if N_elements (N_elements(trial_smooth_and_shift) lt 2) then begin              ; If Amoeba crashes
              xtol = xtol * 1.5 & ftol = ftol * 1.5                                         ; Grow the tolerances
              smooth_and_shift = smooth_and_shift + [RANDOMN(seed), RANDOMN(seed)]          ; Change the guess
              SCALE            = SCALE + [RANDOMN(seed), RANDOMN(seed)]                     ; Change the simplex
              count = count+1
            endif
          ENDREP UNTIL ( (N_elements(trial_smooth_and_shift) eq 2) or (count gt 500) )
          smooth_and_shift = AMOEBAX(Ftol, xtol, function_name='shift_smooth', SCALE = [2., 2.], P0 = smooth_and_shift, FUNCTION_VALUE = fval)
          if N_elements(smooth_and_shift) eq 1 then stop ; check for bugs

        smoothed_shifted = interpolate(gauss_smooth(jup, smooth_and_shift[0]), findgen(n_elements(Jup)) - smooth_and_shift[1])
        Mult_and_add     = mpfitfun('MX_Plus_B', smoothed_shifted[correl_indicies], Spec[correl_indicies], Spec_err[correl_indicies], $
          [mean(Spec[correl_indicies]/smoothed_shifted[correl_indicies]), 0.], /NaN, status = Scatter_fit_status, PERROR = err_a, /quiet)
        scatter_fit      = Mult_and_add[0]*smoothed_shifted + Mult_and_add[1]
        residual         = Spec - scatter_fit
        GOF_array[h]     = GOF + median(abs(residual[correl_indicies]))

        ; Get the linewidth and the determine any systematic wavelength shift from Io's rest
        ; Ultimately we will peg both of these parameters.
        O6_Io_frame               = O6 + O6 * sxpar(header, 'IO_DOPPL') / cspice_clight()
        LSF_fitting_ind1          = cgSetDifference(where( abs(wl - O6_Io_frame) lt 0.3, /NULL), Telluric_AG_Ind)  ; fit region within +/- 0.3A of line center, excluding Telluric airglow
        initial_guess             = [2.e3, O6_Io_frame, 0.077]  ; rough values are fine here [height, wavelength, linewidth_sigma]
        fa                        = {x:double(WL[LSF_fitting_ind1]), y:double(residual[LSF_fitting_ind1]), err:double( sqrt(abs( residual[LSF_fitting_ind1] )) )}
        a                         = mpfit('Gaussian_for_MPFIT', initial_guess, PERROR = err_a, funct=fa, maxiter=50, STATUS = Did_it_work, /Quiet, NPEGGED = NPEGGED)
        Possible_WL_offset[i,h]   = a[1] - O6_Io_frame
        Possible_Line_Width[i,h]  = a[2]
      endfor ; h
      best_GoF = min(GOF_array, best_fit)
      Use_Scatter[i] = Fit_Me_To_The_Scatter[best_fit]
      print, 'Frame: ', Eclipse_files[i], ' is best paired with the Jupiter scatter from: ', Use_Scatter[i], ', Goodness of fit = ', best_GoF
    endfor
  
    Case date of
      'UT180320': begin
          best_scatter_spectrum = 'R_per_A_Jupiter_Center.ec.fits' ; Manually input desired scatting spectrum.
      end
      'UT190812': begin
          best_scatter_spectrum = 'R_per_A_Jupiter_Center.ec.fits' ; Manually input desired scatting spectrum.
      end
    endcase
    print, 'UPON INSPECTION, best choice seems to be: ', best_scatter_spectrum

    ; -------------------------------------------Waterfall Plot Postscript Io and Fit Jovian scatter--------------------------------------------------------------------------------------------------------------

    ; Okay, now that we've established which Jupiter scatter spectrum best optimzes the residual, we can actually use it for scattered light subtraction

    ; get color versus ingress time & setup plot positions
      timeColors = BytScl(findgen(11), Min=0, Max=11, Top=11)
      Color      = timeColors[0]
      cgLoadCT, 33, NColors=8, /reverse
      pos        = cgLayout([1,3], OXMargin=[10,3], OYMargin=[8,5], YGap=0)
      Pos[1,0]   = Pos[1,0]*.8 & Pos[3,1] = Pos[3,1]*.8 
      Pos[1,1]   = Pos[1,1]*.84 & Pos[3,2] = Pos[3,2]*.84 

    cgPS_Open, filename = Reduced_Dir+'O7774_scatterfit.eps', /ENCAPSULATED, xsize = 7.5, ysize = 5.
    !P.font = 1
    device, SET_FONT = 'Helvetica Bold', /TT_FONT
    !p.charsize = 1.5

    cgplot, WL, spec/1.e3, psym = 0, Ytitle = 'kR / '+cgsymbol('Angstrom'), $             ; Rayleighs to KiloRayleighs
      title = 'Io''s O I 7774'+cgsymbol('Angstrom')+' Airglow in Eclipse',  $
      xr = xr, /nodata, pos = pos[*,0], xtickformat = '(A1)', yr = yr

    for i = 0, n_elements(Eclipse_files)-1 do begin
      if i eq 0 then continue ; skip penumbra
      spec         = MRDFITS(reduced_dir+'R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, header, /fscale, /unsigned, /silent)
      spec_err     = MRDFITS(reduced_dir+'sig_R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, sig_header, /fscale, /unsigned, /silent)
      spec         = spec / 1.e3       ; to KR
      spec_err     = spec_err / 1.e3   ; to KR
      WL           = sxpar(header, 'CRVAL1') + findgen(N_elements(spec))*sxpar(header, 'Cdelt1')
      O4_Io_frame  = O4 + O4 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      O5_Io_frame  = O5 + O5 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      O6_Io_frame  = O6 + O6 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      junk         = min( abs(WL - O4_Io_frame), O4_Io_Ind)
      junk         = min( abs(WL - O5_Io_frame), O5_Io_Ind)
      junk         = min( abs(WL - O6_Io_frame), O6_Io_Ind)

      ; Which Jovian scatter are we using?
        ;Jup_Cal_Tell_Corr = MRDFITS(reduced_dir+best_scatter_spectrum, 0, Jup_Scatter_header, /fscale, /unsigned, /silent)
        Jup_Cal_Tell_Corr = MRDFITS(reduced_dir+Use_Scatter[i], 0, Jup_Scatter_header, /fscale, /unsigned, /silent)

      ; Define fitting indicies we need to avoid because of airglow itself.
        Ios_Airglow    = [Io_airglow + Io_airglow * sxpar(header, 'IO_DOPPL') / cspice_clight()]                                      ; wavelength locations of all airglow lines, Ignore telluric K
        ARCES_1sigma_bandwidth = (mean(xr)/31500.) / 2.3548 ; FWHM in Angstroms / 2sqrt(2ln2) = 68% of flux enclosed within 1 sigma --> Appropriate for weak telluric airglow
        Io_AG_Ind = [] & Telluric_AG_Ind = []
        for j = 0, N_elements(Ios_Airglow)-1 do Io_AG_ind = [Io_AG_Ind, where(abs(Ios_Airglow[j] - WL) lt 2.5*ARCES_1sigma_bandwidth, /NULL)] ;  indicies within 2.5 sigma of an Io airglow line
        for j = 0, N_elements(Telluric_Airglow)-1 do Telluric_AG_ind = [Telluric_AG_Ind, where(abs(Telluric_Airglow[j] - WL) lt ARCES_1sigma_bandwidth, /NULL)] ; indicies within 1 sigma of any Telluric airglow
        AG_ind = [Io_AG_Ind, Telluric_AG_Ind]
        AG_ind = AG_ind[UNIQ(AG_ind, SORT(AG_ind))]

      ; Fit the Jovian scattered light, use Amoeba/Correlation to match Jupiter in terms of smoothing and shifting, then multiply and add until matches the scattered background
        correl_indicies  = cgSetDifference(include_Wls, AG_Ind)
        weights_7774     = 1./abs(O4_Io_Ind - correl_indicies)^(1.5) + 1./abs(O5_Io_Ind - correl_indicies)^(1.5) + 1./abs(O6_Io_Ind - correl_indicies)^(1.5)
        Mult_and_add     = mpfitfun('MX_Plus_B', Jup_Cal_Tell_Corr[correl_indicies], Spec[correl_indicies], Spec_err[correl_indicies], parinfo = MX_plus_B_parinfo, $
                                    [mean(Spec[correl_indicies]/Jup_Cal_Tell_Corr[correl_indicies]), 0.], /NaN, status = Scatter_fit_status, PERROR = err_a, /quiet)
        jup              = Mult_and_add[0]*Jup_Cal_Tell_Corr + Mult_and_add[1]

        ; Iterate Amoeba until it gives an answer
        Ftol = 1.e-4 &  Xtol = 1.e-4 & count = 0                                          ; Define initial tolerances for Amoeba
        smooth_and_shift = [2., 0.] & scale  = [2., 2.]                                   ; Define intial guess and simplex for Amoeba
        REPEAT BEGIN
          trial_smooth_and_shift = AMOEBAX(Ftol, xtol, function_name='shift_smooth', SCALE = Scale, P0 = smooth_and_shift, FUNCTION_VALUE = fval)
          if N_elements (N_elements(trial_smooth_and_shift) lt 2) then begin              ; If Amoeba crashes
            xtol = xtol * 1.5 & ftol = ftol * 1.5                                         ; Grow the tolerances
            smooth_and_shift = smooth_and_shift + [RANDOMN(seed), RANDOMN(seed)]          ; Change the guess
            SCALE            = SCALE + [RANDOMN(seed), RANDOMN(seed)]                     ; Change the simplex
            count = count+1
          endif
        ENDREP UNTIL ( (N_elements(trial_smooth_and_shift) eq 2) or (count gt 500) )
        smooth_and_shift = AMOEBAX(Ftol, xtol, function_name='shift_smooth', SCALE = [2., 2.], P0 = smooth_and_shift, FUNCTION_VALUE = fval)
        if N_elements(smooth_and_shift) eq 1 then stop ; check for bugs

      smoothed_shifted = interpolate(gauss_smooth(jup, smooth_and_shift[0]), findgen(n_elements(Jup)) - smooth_and_shift[1])
      
      Mult_and_add     = mpfitfun('MX_Plus_B', smoothed_shifted[correl_indicies], Spec[correl_indicies], Spec_err[correl_indicies], [mean(Spec[correl_indicies]/smoothed_shifted[correl_indicies]), 0.], $
                                  /NaN, status = Scatter_fit_status, PERROR = err_a, /quiet, weights = weights_7774, parinfo = MX_plus_B_parinfo) ; weights override error, local indices weighted 
      
      scatter_fit      = Mult_and_add[0]*smoothed_shifted + Mult_and_add[1]
      residual         = Spec - scatter_fit
      residual_err     = make_array( N_elements(residual), value = stddev(residual[correl_indicies]) )

      ; plot the spectrum, fit and indices used for the fit (if we're in debug mode)
        ;cgplot, WL[correl_indicies], scatter_fit[correl_indicies], psym=14, /overplot
        cgplot, WL, scatter_fit, color = timeColors[i], linestyle = 1, thick = 5, /overplot
        cgplot, WL, spec, color = timeColors[i], thick = 5, /overplot

      ; Now do the line fitting to the residual. First define Line Spread Function (LSF) Gaussian fit parameter info
        parinfo               = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0]}, 3)
        parinfo[1].fixed      = 1
        parinfo[1].value      = double(O6_Io_frame)                                                           ; Pin the line's wavelength
        parinfo[2].fixed      = 1
        parinfo[2].value      = median(Possible_line_width[1:*,*])
        LSF_fitting_ind1      = cgSetDifference(where( abs(wl - O6_Io_frame) lt 0.3, /NULL), Telluric_AG_Ind) ; fit region within +/- 0.3A of line center, excluding Telluric airglow

      ; Gaussian fit
        initial_guess = [5.D, parinfo[1].value, parinfo[2].value]
        fa            = {x:double(WL[LSF_fitting_ind1]), y:double(residual[LSF_fitting_ind1]), err:double( sqrt(abs( spec[LSF_fitting_ind1] )) )}
        a             = mpfit('Gaussian_for_MPFIT', initial_guess, PERROR = err_a, funct=fa, parinfo = parinfo, STATUS = Did_it_work, /Quiet, NPEGGED = NPEGGED)

      ; write the results
        WL_array[*, i]           = WL[include_WLs]
        plot_residual_array[*, i]= residual[plot_WLs]
        LSF_Fit_array[*, i]      = gaussian(WL[include_WLs], a)                         ; LSF_Fit
        Brightness_array[i]      = A[0]*A[2]*SQRT(2*!DPI)                               ; Area under the Gaussian
        err_Brightness_array[i]  = stddev(residual[correl_indicies])*A[2]*SQRT(2*!DPI)  ; use the std dev of the residual and the width of the line
        T_P_Shadow_array[i]      = sxpar(header, 'T_PSHADO')
        T_U_Shadow_array[i]      = sxpar(header, 'T_USHADO')
        EXPTime_array[i]         = sxpar(sig_header, 'EXPTIME') / 60.                   ; for some reason, this keyword is bad in the regular headers, so use sig_header as a workaround
        DopplerShift_array1[i]   = O4_Io_frame
        DopplerShift_array2[i]   = O5_Io_frame
        DopplerShift_array3[i]   = O6_Io_frame
    endfor

    ; Plot the residual and Gaussian fits
      N_frames            = total(Brightness_array gt 0.)
      cgplot, spec, WL, psym = 0, yminor = 2, $
        yr = [-1.99, 4.9], xr = XR, /nodata, pos = pos[*,1], xtickformat= '(A1)', /noerase
      cgtext, xr[0] - 0.8, -2., 'Residual (kR / '+cgsymbol('Angstrom')+')', orientation = 90, alignment = 0.5
      ; Co-align to Io's Doppler Shift ???
        shift_By = mean( [transpose( DopplerShift_array1 - Mean(DopplerShift_array1[1:*]) ), $
                          transpose( DopplerShift_array2 - Mean(DopplerShift_array2[1:*]) ), $
                          transpose( DopplerShift_array3 - Mean(DopplerShift_array3[1:*]) )], dimension = 1 )
        aligned_residuals = plot_residual_array
  
      for i = 1, n_elements(Eclipse_files)-1 do begin
        cgplot, WL[plot_WLs], plot_residual_array[*, i], /OVERPLOT, COLOR = timeColors[i], psym=10
        cgplot, [DopplerShift_array1[i], DopplerShift_array1[i]], [4.9, 2.5], COLOR = timeColors[i], /overplot
        cgplot, [DopplerShift_array2[i], DopplerShift_array2[i]], [4.9, 4.5], COLOR = timeColors[i], /overplot
        cgplot, [DopplerShift_array3[i], DopplerShift_array3[i]], [4.9, 2.5], COLOR = timeColors[i], /overplot
        aligned_residuals[*,i] = interpol( plot_residual_array[*,i], WL[plot_WLs], WL[plot_WLs] - Shift_By[i]) 
      endfor
      co_add = total(aligned_residuals[*,1:*], 2) / N_frames
 
      ; plot the aligned & co-added spectra      
        cgplot, WL[plot_WLs], co_add, pos = pos[*,2], thick=5, /noerase, xr = xr, yr = [-0.6, 1.6], /yst,  $
          xtitle = 'Wavelength ('+cgsymbol('Angstrom')+')'                
        cgplot, xr, [0.,0.], linestyle = 2, /overplot
        cgplot, [Mean(DopplerShift_array1[1:*]), Mean(DopplerShift_array1[1:*])], [-0.6, 1.6], linestyle = 1, /overplot
        cgplot, [Mean(DopplerShift_array2[1:*]), Mean(DopplerShift_array2[1:*])], [-0.6, 1.6], linestyle = 1, /overplot
        cgplot, [Mean(DopplerShift_array3[1:*]), Mean(DopplerShift_array3[1:*])], [-0.6, 1.6], linestyle = 1, /overplot
    cgps_close
    
    ; Estimate the total time-averaged brightness of the emission feature or features
      integrate_multiplet = cgSetIntersection(Io_AG_Ind, plot_WLs) - min(plot_WLs) ;indices of the multiplet over
      print, 'Disk- and Time-Averaged 7774A Brightness =', total(co_add[integrate_multiplet]) * N_elements(integrate_multiplet)*sxpar(header, 'Cdelt1') / N_frames, ' KiloRayleighs'      
      window, 0
      x = WL[plot_WLs]
      cgplot, x[integrate_multiplet], co_add[integrate_multiplet] / N_frames, /xst, psym=16, ytitle = 'Disk- and Time-Averaged 7774A Brightness kR/A)'
      cgplot, x[integrate_multiplet], co_add[integrate_multiplet] / N_frames, /overplot, psym=10
  endif  

  if Part eq 1.04 then begin ; Make waterfall plots, O 8446 requires telluric corrections

    ;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ;--------------------------------------------------------Get the O8446 Time Series and Scatter Fit-------------------------------------------------------------------------------
    ;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ; METHOD:
    ;   Loop though all possible spectra that could be fit to jovian scatter and *MANUALLY* choose what works best. Avoid Io and telluric airglow in the fitting.
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

    ; The is a local error in the Wavelength solution near the 8446A lines, likely because ThAr calibration lines saturate in this region during the 200s APO recommended exopsure time
    ; Estimate the error in the wavelength solution so that it can be corrected. 
      Io_Sunlit        = MRDFITS(reduced_dir+'fullspecIo_Sunlit.0001.ec.fits', 0, header, /fscale, /unsigned )
      WL               = sxpar(header, 'CRVAL1')+findgen(N_elements(spec))*sxpar(header, 'Cdelt1')
      line             = 8446.382   ; Bass2000 wavelength of the Solar O I line core
      WL_Solution_Line = 8446.310   ; Location best matching the line core by visual inspection
  
      ; find and plot the intantaneous Earth-Io Doppler Shift and light time
        cspice_UTC2ET, sxpar(header, 'DATE-OBS'), ET
        cspice_spkezr, 'Io', ET + float(sxpar(header, 'EXPTIME'))/2., 'J2000', 'LT', 'Earth', Io_Earth_State, Earth2Io_ltime
        theta  = cspice_vsep(Io_Earth_state[0:2], Io_Earth_state[3:5])
        Io_wrt_Earth_Dopplershift = line * cos(theta) * norm(Io_Earth_State[3:5]) / cspice_clight()
  
      ; find Sun-Io Doppler Shift and light time
        cspice_spkezr, 'Io', ET + float(sxpar(header, 'EXPTIME'))/2. - Earth2Io_ltime, 'J2000', 'LT', 'Sun', Io_Sun_State, Io2Sun_ltime
        theta  = cspice_vsep(Io_Sun_state[0:2], Io_Sun_state[3:5])
        Io_wrt_Sun_Dopplershift = line * cos(theta) * norm(Io_Sun_State[3:5]) / cspice_clight()
  
      ; Apply the light time corrected 2-way Doppler Shifts
        reflected_line             = line + Io_wrt_Earth_Dopplershift + Io_wrt_Sun_Dopplershift
        reflected_WL_Solution_Line = WL_Solution_Line + Io_wrt_Earth_Dopplershift + Io_wrt_Sun_Dopplershift
  
      ; Visual Inspection: Red line is in the core?
        window, 0
        cgplot, wl, Io_Sunlit, xr = [reflected_line-1., reflected_line+1.], xstyle = 1, /ynozero
        cgplot, [reflected_line, reflected_line], [0., 1.], /overplot
        cgplot, [reflected_WL_Solution_Line, reflected_WL_Solution_Line], [0., 1.], /overplot, color = 'red'
        WS_Correction = reflected_line - reflected_WL_Solution_Line
        print, 'Wavelengths here are off by: ', WS_Correction*1000., +'mA. ---> shifting all ARCES spectra to the red near 8446A.

    ;-------------------------------------------------------------------------------------------------------------------------------------------

    ; Fit various spectra to the Jovian Scatter, *use* whichever one produces the lowest residual. Note this method will never work if there's much telluric absorption.
    Case date of
      'UT180320': begin
        Fit_Me_To_The_Scatter    = ['R_per_A_Jupiter_Center.ec.fits', $
                                  'R_per_A_Jovian_Scatter.0001.ec.fits', $
                                  'R_per_A_Jovian_Scatter.0002.ec.fits', $
                                  'R_per_A_Jovian_Scatter.0003.ec.fits' ]
                  end
      'UT190812': begin
        Fit_Me_To_The_Scatter    = ['R_per_A_Jupiter_Center.ec.fits']
                  end
    endcase

    ; setup the plot axis
      spec      = MRDFITS(reduced_dir+'R_per_A_'+Eclipse_files[1]+'.ec.fits', 0, header, /fscale, /unsigned, /silent)
      WL        = sxpar(header, 'CRVAL1')+findgen(N_elements(spec))*sxpar(header, 'Cdelt1')
      bandwidth = 6.                                                                               ; half with of the plot's x axis in Angstroms
  
      Case date of
        'UT180320': begin
          YR   = [10.3, 30.3]
          XR   = [O8 + O8 * sxpar(header, 'IO_DOPPL') / cspice_clight() - Bandwidth, O8 + O8 * sxpar(header, 'IO_DOPPL') / cspice_clight() + Bandwidth] 
        end
        'UT190812': begin
          yr   = [5.5e1, 2e2]
          YR_resid = [-2, 20]
          XR   = [O8 + O8 * sxpar(header, 'IO_DOPPL') / cspice_clight() - Bandwidth, O8 + O8 * sxpar(header, 'IO_DOPPL') / cspice_clight() + Bandwidth]
          end
     endcase

    ; Define arrays
      include_WLs            = where( (wl gt xr[0]) and (wl lt xr[1]), /NULL )                      ; wavelengths to use in the scattered light fitting
      plot_WLs               = where( (wl gt xr[0]) and (wl lt xr[1]), /NULL )                      ; wavelengths to actually plot
      residual_array         = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
      plot_residual_array    = fltarr(N_elements(plot_WLs), n_elements(Eclipse_files))
      LSF_Fit_array          = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
      WL_array               = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
      Brightness_array       = fltarr(n_elements(Eclipse_files))
      err_Brightness_array   = fltarr(n_elements(Eclipse_files))
      T_P_Shadow_array       = fltarr(n_elements(Eclipse_files))
      T_U_Shadow_array       = fltarr(n_elements(Eclipse_files))
      EXPTime_array          = fltarr(n_elements(Eclipse_files))
      DopplerShift_array1    = fltarr(n_elements(Eclipse_files))                                    ; For the first line in the triplet
      DopplerShift_array2    = fltarr(n_elements(Eclipse_files))
      DopplerShift_array3    = fltarr(n_elements(Eclipse_files))
      Gof                    = 0.                                                                   ; Goodness of fit (co-addded mean residual)
      GOF_array              = fltarr(N_elements(Fit_Me_To_The_Scatter))
      Use_Scatter            = strarr(N_elements(Eclipse_files))
      Possible_WL_offset     = fltarr(N_elements(Eclipse_files), N_elements(Fit_Me_To_The_Scatter)) ; Save any offsets from Io's expected rest wavelength, then correct any error in the wavelength solution.
      Possible_line_width    = fltarr(N_elements(Eclipse_files), N_elements(Fit_Me_To_The_Scatter)) ; Save the linewidths from all type of scatter subtraction, then average and peg this parameter.

    for i = 0, n_elements(Eclipse_files)-1 do begin
      if i eq 0 then continue ; skip penumbra
      spec        = MRDFITS(reduced_dir+'R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, header, /fscale, /unsigned, /silent)
      spec_err    = MRDFITS(reduced_dir+'sig_R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, sig_header, /fscale, /unsigned, /silent)
      WL          = sxpar(header, 'CRVAL1')+findgen(N_elements(spec))*sxpar(header, 'Cdelt1') + WS_Correction

      for h = 0, n_elements(Fit_Me_To_The_Scatter)-1 do begin
        ; Which Jovian scatter are we using?
        Jup_Cal_Tell_Corr = MRDFITS(reduced_dir+Fit_Me_To_The_Scatter[h], 0, Jupiter_Scatter_header, /fscale, /unsigned, /silent)

        ; Define fitting indicies we need to avoid because of airglow itself.
        Ios_Airglow    = [Io_airglow + Io_airglow * sxpar(header, 'IO_DOPPL') / cspice_clight()]                                      ; wavelength locations of all airglow lines, Ignore telluric K
        ARCES_1sigma_bandwidth = (mean(xr)/31500.) / 2.3548 ; 1. * FWHM in Angstroms / 2sqrt(2ln2) = 68% of flux enclosed within 1 sigma --> Appropriate for weak telluric airglow
        Io_AG_Ind = [] & Telluric_AG_Ind = []
        for j = 0, N_elements(Ios_Airglow)-1 do Io_AG_ind = [Io_AG_Ind, where(abs(Ios_Airglow[j] - WL) lt ARCES_1sigma_bandwidth, /NULL)] ;  indicies within 2 sigma of an Io airglow line
        for j = 0, N_elements(Telluric_Airglow)-1 do Telluric_AG_ind = [Telluric_AG_Ind, where(abs(Telluric_Airglow[j] - WL) lt ARCES_1sigma_bandwidth, /NULL)] ; indicies within 1 sigma of any Telluric airglow
        AG_ind = [Io_AG_Ind, Telluric_AG_Ind]
        AG_ind = AG_ind[UNIQ(AG_ind, SORT(AG_ind))]

        ; Fit the Jovian scattered light, use Amoeba/Correlation to match Jupiter in terms of smoothing and shifting, then multiply and add until matches the scattered background
        correl_indicies  = cgSetDifference(include_Wls, AG_Ind)
        Mult_and_add     = mpfitfun('MX_Plus_B', Jup_Cal_Tell_Corr[correl_indicies], Spec[correl_indicies], Spec_err[correl_indicies], parinfo = MX_plus_B_parinfo,$
                                     [mean(Spec[correl_indicies]/Jup_Cal_Tell_Corr[correl_indicies]), 0.], /NaN, status = Scatter_fit_status, PERROR = err_a, /quiet)
        jup              = Mult_and_add[0]*Jup_Cal_Tell_Corr + Mult_and_add[1]

          ; Iterate Amoeba until it gives an answer
          Ftol = 1.e-4 &  Xtol = 1.e-4 & count = 0                                          ; Define initial tolerances for Amoeba
          smooth_and_shift = [2., 0.] & scale  = [2., 2.]                                   ; Define intial guess and simplex for Amoeba
          REPEAT BEGIN
            trial_smooth_and_shift = AMOEBAX(Ftol, xtol, function_name='shift_smooth', SCALE = Scale, P0 = smooth_and_shift, FUNCTION_VALUE = fval)
            if N_elements (N_elements(trial_smooth_and_shift) lt 2) then begin              ; If Amoeba crashes
              xtol = xtol * 1.5 & ftol = ftol * 1.5                                         ; Grow the tolerances
              smooth_and_shift = smooth_and_shift + [RANDOMN(seed), RANDOMN(seed)]          ; Change the guess
              SCALE            = SCALE + [RANDOMN(seed), RANDOMN(seed)]                     ; Change the simplex
              count = count+1
            endif
          ENDREP UNTIL ( (N_elements(trial_smooth_and_shift) eq 2) or (count gt 500) )
          smooth_and_shift = AMOEBAX(Ftol, xtol, function_name='shift_smooth', SCALE = [2., 2.], P0 = smooth_and_shift, FUNCTION_VALUE = fval)
          if N_elements(smooth_and_shift) eq 1 then stop ; check for bugs

        smoothed_shifted = interpolate(gauss_smooth(jup, smooth_and_shift[0]), findgen(n_elements(Jup)) - smooth_and_shift[1])
        Mult_and_add     = mpfitfun('MX_Plus_B', smoothed_shifted[correl_indicies], Spec[correl_indicies], Spec_err[correl_indicies], parinfo = MX_plus_B_parinfo, $
                                    [mean(Spec[correl_indicies]/smoothed_shifted[correl_indicies]), 0.], /NaN, status = Scatter_fit_status, PERROR = err_a, /quiet)
        scatter_fit      = Mult_and_add[0]*smoothed_shifted + Mult_and_add[1]
        residual         = Spec - scatter_fit
        GOF_array[h]     = GOF + median(abs(residual[correl_indicies]))

        ; Get the linewidth and the determine any systematic wavelength shift from Io's rest
        ; Ultimately we will peg both of these parameters.
        O9_Io_frame               = O9 + O9 * sxpar(header, 'IO_DOPPL') / cspice_clight()
        LSF_fitting_ind1          = cgSetDifference(where( abs(wl - O9_Io_frame) lt 0.3, /NULL), Telluric_AG_Ind)  ; fit region within +/- 0.3A of line center, excluding Telluric airglow
        initial_guess             = [2.e3, O9_Io_frame, 0.077]  ; rough values are fine here [height, wavelength, linewidth_sigma]
        fa                        = {x:double(WL[LSF_fitting_ind1]), y:double(residual[LSF_fitting_ind1]), err:double( sqrt(abs( residual[LSF_fitting_ind1] )) )}
        a                         = mpfit('Gaussian_for_MPFIT', initial_guess, PERROR = err_a, funct=fa, maxiter=50, STATUS = Did_it_work, /Quiet, NPEGGED = NPEGGED)
        Possible_WL_offset[i,h]   = a[1] - O9_Io_frame
        Possible_Line_Width[i,h]  = a[2]
      endfor ; h
      best_GoF = min(GOF_array, best_fit)
      Use_Scatter[i] = Fit_Me_To_The_Scatter[best_fit]
      print, 'Frame: ', Eclipse_files[i], ' is best paired with the Jupiter scatter from: ', Use_Scatter[i], ', Goodness of fit = ', best_GoF
    endfor
    Case date of
      'UT180320': begin
          best_scatter_spectrum = 'R_per_A_Jupiter_Center.ec.fits' ; Manually input desired scatting spectrum.
      end
      'UT190812': begin
          best_scatter_spectrum = 'R_per_A_Jupiter_Center.ec.fits' ; Manually input desired scatting spectrum.
      end
    endcase
    print, 'UPON INSPECTION, best choice seems to be: ', best_scatter_spectrum

    ; -------------------------------------------Waterfall Plot Postscript Io and Fit Jovian scatter--------------------------------------------------------------------------------------------------------------

    ; Okay, now that we've established which Jupiter scatter spectrum best optimzes the residual, we can actually use it for scattered light subtraction

    ; get color versus ingress time & setup plot positions
      timeColors = BytScl(findgen(11), Min=0, Max=11, Top=11)
      Color      = timeColors[0]
      cgLoadCT, 33, NColors=8, /reverse
      pos        = cgLayout([1,3], OXMargin=[10,3], OYMargin=[8,5], YGap=0)
      Pos[1,0]   = Pos[1,0]*.8 & Pos[3,1] = Pos[3,1]*.8 
      Pos[1,1]   = Pos[1,1]*.84 & Pos[3,2] = Pos[3,2]*.84 

    cgPS_Open, filename = Reduced_Dir+'O8446_scatterfit.eps', /ENCAPSULATED, xsize = 7.5, ysize = 5.
    !P.font = 1
    device, SET_FONT = 'Helvetica Bold', /TT_FONT
    !p.charsize = 1.5

    cgplot, WL, spec/1.e3, psym = 0, Ytitle = 'kR / '+cgsymbol('Angstrom'), $             ; Rayleighs to KiloRayleighs
      title = 'Io''s O I 8446'+cgsymbol('Angstrom')+' Airglow in Eclipse',  $
      xr = xr, /nodata, pos = pos[*,0], xtickformat = '(A1)', yr = yr

    for i = 0, n_elements(Eclipse_files)-1 do begin
      if i eq 0 then continue ; skip penumbra
      spec         = MRDFITS(reduced_dir+'R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, header, /fscale, /unsigned, /silent)
      spec_err     = MRDFITS(reduced_dir+'sig_R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, sig_header, /fscale, /unsigned, /silent)
      spec         = spec / 1.e3       ; to KR
      spec_err     = spec_err / 1.e3   ; to KR
      WL           = sxpar(header, 'CRVAL1') + findgen(N_elements(spec))*sxpar(header, 'Cdelt1') + WS_Correction
      O7_Io_frame  = O7 + O7 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      O8_Io_frame  = O8 + O8 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      O9_Io_frame  = O9 + O9 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      junk         = min( abs(WL - O7_Io_frame), O7_Io_Ind)
      junk         = min( abs(WL - O8_Io_frame), O8_Io_Ind)
      junk         = min( abs(WL - O9_Io_frame), O9_Io_Ind)

      ; Which Jovian scatter are we using?
        Jup_Cal_Tell_Corr = MRDFITS(reduced_dir+best_scatter_spectrum, 0, Jup_Scatter_header, /fscale, /unsigned, /silent)

      ; Define fitting indicies we need to avoid because of airglow itself.
      Ios_Airglow    = [Io_airglow + Io_airglow * sxpar(header, 'IO_DOPPL') / cspice_clight()]                                      ; wavelength locations of all airglow lines, Ignore telluric K
      ARCES_1sigma_bandwidth = (mean(xr)/31500.) / 2.3548 ; FWHM in Angstroms / 2sqrt(2ln2) = 68% of flux enclosed within 1 sigma --> Appropriate for weak telluric airglow
      Io_AG_Ind = [] & Telluric_AG_Ind = []
      for j = 0, N_elements(Ios_Airglow)-1 do Io_AG_ind = [Io_AG_Ind, where(abs(Ios_Airglow[j] - WL) lt 2.*ARCES_1sigma_bandwidth, /NULL)] ;  indicies within 2 sigma of an Io airglow line
      for j = 0, N_elements(Telluric_Airglow)-1 do Telluric_AG_ind = [Telluric_AG_Ind, where(abs(Telluric_Airglow[j] - WL) lt ARCES_1sigma_bandwidth, /NULL)] ; indicies within 1 sigma of any Telluric airglow
      AG_ind = [Io_AG_Ind, Telluric_AG_Ind]
      AG_ind = AG_ind[UNIQ(AG_ind, SORT(AG_ind))]

      ; Fit the Jovian scattered light, use Amoeba/Correlation to match Jupiter in terms of smoothing and shifting, then multiply and add until matches the scattered background
      correl_indicies  = cgSetDifference(include_Wls, AG_Ind)
      weights_8446     = 1./abs(O7_Io_Ind - correl_indicies)^(1.5) + 1./abs(O8_Io_Ind - correl_indicies)^(1.5) + 1./abs(O9_Io_Ind - correl_indicies)^(1.5)
      Mult_and_add     = mpfitfun('MX_Plus_B', Jup_Cal_Tell_Corr[correl_indicies], Spec[correl_indicies], Spec_err[correl_indicies], parinfo = MX_plus_B_parinfo, $
                                  [mean(Spec[correl_indicies]/Jup_Cal_Tell_Corr[correl_indicies]), 0.], /NaN, status = Scatter_fit_status, PERROR = err_a, /quiet)
      jup              = Mult_and_add[0]*Jup_Cal_Tell_Corr + Mult_and_add[1]

          ; Iterate Amoeba until it gives an answer
          Ftol = 1.e-4 &  Xtol = 1.e-4 & count = 0                                          ; Define initial tolerances for Amoeba
          smooth_and_shift = [2., 0.] & scale  = [2., 2.]                                   ; Define intial guess and simplex for Amoeba
          REPEAT BEGIN
            trial_smooth_and_shift = AMOEBAX(Ftol, xtol, function_name='shift_smooth', SCALE = Scale, P0 = smooth_and_shift, FUNCTION_VALUE = fval)
            if N_elements (N_elements(trial_smooth_and_shift) lt 2) then begin              ; If Amoeba crashes
              xtol = xtol * 1.5 & ftol = ftol * 1.5                                         ; Grow the tolerances
              smooth_and_shift = smooth_and_shift + [RANDOMN(seed), RANDOMN(seed)]          ; Change the guess
              SCALE            = SCALE + [RANDOMN(seed), RANDOMN(seed)]                     ; Change the simplex
              count = count+1
            endif
          ENDREP UNTIL ( (N_elements(trial_smooth_and_shift) eq 2) or (count gt 500) )
          smooth_and_shift = AMOEBAX(Ftol, xtol, function_name='shift_smooth', SCALE = [2., 2.], P0 = smooth_and_shift, FUNCTION_VALUE = fval)
          if N_elements(smooth_and_shift) eq 1 then stop ; check for bugs

      smoothed_shifted = interpolate(gauss_smooth(jup, smooth_and_shift[0]), findgen(n_elements(Jup)) - smooth_and_shift[1])

      Mult_and_add     = mpfitfun('MX_Plus_B', smoothed_shifted[correl_indicies], Spec[correl_indicies], Spec_err[correl_indicies], [mean(Spec[correl_indicies]/smoothed_shifted[correl_indicies]), 0.], $
        /NaN, status = Scatter_fit_status, PERROR = err_a, /quiet, weights = weights_8446, parinfo = MX_plus_B_parinfo) ; weights override error, local indices weighted

      scatter_fit      = Mult_and_add[0]*smoothed_shifted + Mult_and_add[1]
      residual         = Spec - scatter_fit
      residual_err     = make_array( N_elements(residual), value = stddev(residual[correl_indicies]) )

      ; plot the spectrum, fit and indices used for the fit (if we're in debug mode)
        ;cgplot, WL[correl_indicies], scatter_fit[correl_indicies], psym=14, /overplot
        cgplot, WL, scatter_fit, color = timeColors[i], linestyle = 1, thick = 7, /overplot
        cgplot, WL, spec, color = timeColors[i], thick = 5, /overplot

      ; Now do the line fitting to the residual. First define Line Spread Function (LSF) Gaussian fit parameter info
        parinfo               = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0]}, 3)
        parinfo[1].fixed      = 1
        parinfo[1].value      = double(O9_Io_frame)                                                           ; Pin the line's wavelength
        parinfo[2].fixed      = 1
        parinfo[2].value      = median(Possible_line_width[1:*,*])
        LSF_fitting_ind1      = cgSetDifference(where( abs(wl - O9_Io_frame) lt 0.3, /NULL), Telluric_AG_Ind) ; fit region within +/- 0.3A of line center, excluding Telluric airglow

      ; Gaussian fit
        initial_guess = [5.D, parinfo[1].value, parinfo[2].value]
        fa            = {x:double(WL[LSF_fitting_ind1]), y:double(residual[LSF_fitting_ind1]), err:double( sqrt(abs( spec[LSF_fitting_ind1] )) )}
        a             = mpfit('Gaussian_for_MPFIT', initial_guess, PERROR = err_a, funct=fa, parinfo = parinfo, STATUS = Did_it_work, /Quiet, NPEGGED = NPEGGED)

      ; write the results
        WL_array[*, i]           = WL[include_WLs]
        plot_residual_array[*, i]= residual[plot_WLs]
        LSF_Fit_array[*, i]      = gaussian(WL[include_WLs], a)                         ; LSF_Fit
        Brightness_array[i]      = A[0]*A[2]*SQRT(2*!DPI)                               ; Area under the Gaussian
        err_Brightness_array[i]  = stddev(residual[correl_indicies])*A[2]*SQRT(2*!DPI)  ; use the std dev of the residual and the width of the line
        T_P_Shadow_array[i]      = sxpar(header, 'T_PSHADO')
        T_U_Shadow_array[i]      = sxpar(header, 'T_USHADO')
        EXPTime_array[i]         = sxpar(sig_header, 'EXPTIME') / 60.                   ; for some reason, this keyword is bad in the regular headers, so use sig_header as a workaround
        DopplerShift_array1[i]   = O7_Io_frame
        DopplerShift_array2[i]   = O8_Io_frame
        DopplerShift_array3[i]   = O9_Io_frame
    endfor

    ; Plot the residual and Gaussian fits
      N_frames            = total(Brightness_array gt 0.)
      cgplot, spec, WL, psym = 0, xtickformat = '(A1)', yminor = 2, $
        yr = [-1.99,3.5 ], xr = XR, /nodata, pos = pos[*,1], /noerase      
      cgtext, xr[0] - 0.8, -2., 'Residual (kR / '+cgsymbol('Angstrom')+')', orientation = 90, alignment = 0.5
      
      ; Co-align to Io's Doppler Shift 
        shift_By = mean( [transpose( DopplerShift_array1 - Mean(DopplerShift_array1[1:*]) ), $
          transpose( DopplerShift_array2 - Mean(DopplerShift_array2[1:*]) ), $
          transpose( DopplerShift_array3 - Mean(DopplerShift_array3[1:*]) )], dimension = 1 )
        aligned_residuals = plot_residual_array
  
      for i = 1, n_elements(Eclipse_files)-1 do begin
        cgplot, WL[plot_WLs], plot_residual_array[*, i], /OVERPLOT, COLOR = timeColors[i], psym=10
        cgplot, [DopplerShift_array1[i], DopplerShift_array1[i]], [3.5, 2.5], COLOR = timeColors[i], /overplot
        cgplot, [DopplerShift_array2[i], DopplerShift_array2[i]], [3.5, 3.0], COLOR = timeColors[i], /overplot
        cgplot, [DopplerShift_array3[i], DopplerShift_array3[i]], [3.5, 2.0], COLOR = timeColors[i], /overplot
        aligned_residuals[*,i] = interpol( plot_residual_array[*,i], WL[plot_WLs], WL[plot_WLs] - Shift_By[i])
      endfor
      co_add = total(aligned_residuals[*,1:*], 2) / n_frames         
      
      ; plot the aligned & co-added spectra
        cgplot, WL[plot_WLs], co_add, pos = pos[*,2], thick=5, /noerase, xr = xr, yr = [-0.6, 1.6], /yst,  $
          xtitle = 'Wavelength ('+cgsymbol('Angstrom')+')'                ; plot the aligned & co-added spectra
        cgplot, xr, [0.,0.], linestyle = 2, /overplot
        cgplot, [Mean(DopplerShift_array1[1:*]), Mean(DopplerShift_array1[1:*])], [-0.6, 1.6], linestyle = 1, /overplot
        cgplot, [Mean(DopplerShift_array2[1:*]), Mean(DopplerShift_array2[1:*])], [-0.6, 1.6], linestyle = 1, /overplot
        cgplot, [Mean(DopplerShift_array3[1:*]), Mean(DopplerShift_array3[1:*])], [-0.6, 1.6], linestyle = 1, /overplot
    cgps_close
    
    ; Estimate the total time-averaged brightness of the emission feature or features
      integrate_multiplet = cgSetIntersection(Io_AG_Ind, plot_WLs) - min(plot_WLs) ;indices of the multiplet over
      print, 'Disk- and Time-Averaged 8446A Brightness =', total(co_add[integrate_multiplet]) * N_elements(integrate_multiplet)*sxpar(header, 'Cdelt1') / N_frames, ' KiloRayleighs'  
      window, 0
      x = WL[plot_WLs]
      cgplot, x[integrate_multiplet], co_add[integrate_multiplet] / N_frames, /xst, psym=16, ytitle = 'Disk- and Time-Averaged 8446A Brightness kR/A)''
      cgplot, x[integrate_multiplet], co_add[integrate_multiplet] / N_frames, /overplot, psym=10
  endif

  if Part eq 1.1 then begin ; Make waterfall plots, S 9225 (requires telluric corrections?)

    ;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ;--------------------------------------------------------Get the O9225 Time Series and Scatter Fit-------------------------------------------------------------------------------
    ;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ; METHOD:
    ;   Loop though all possible spectra that could be fit to jovian scatter and *MANUALLY* choose what works best. Avoid Io and telluric airglow in the fitting.
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
    Case date of
      'UT180320': begin
        Fit_Me_To_The_Scatter    = ['R_per_A_Jupiter_Center.ec.fits', $
                                  'R_per_A_Jovian_Scatter.0001.ec.fits', $
                                  'R_per_A_Jovian_Scatter.0002.ec.fits', $
                                  'R_per_A_Jovian_Scatter.0003.ec.fits' ]
                  end
      'UT190812': begin
        Fit_Me_To_The_Scatter    = ['R_per_A_Jupiter_Center.ec.fits']
                  end
    endcase

    ; setup the plot axis
      spec      = MRDFITS(reduced_dir+'R_per_A_'+Eclipse_files[1]+'.ec.fits', 0, header, /fscale, /unsigned, /silent)
      WL        = sxpar(header, 'CRVAL1')+findgen(N_elements(spec))*sxpar(header, 'Cdelt1')
      bandwidth = 16.                                                                               ; half with of the plot's x axis in Angstroms

    Case date of
      'UT180320': begin
        YR   = [9, 26]
        XR   = [S2 + S2 * sxpar(header, 'IO_DOPPL') / cspice_clight() - Bandwidth, S2 + S2 * sxpar(header, 'IO_DOPPL') / cspice_clight() + Bandwidth] - 3.
      end
      'UT190812': begin
        yr   = [5.5e1, 2e2]
        YR_resid = [-2, 20]
        XR   = [S2 + S2 * sxpar(header, 'IO_DOPPL') / cspice_clight() - Bandwidth, S2 + S2 * sxpar(header, 'IO_DOPPL') / cspice_clight() + Bandwidth] - 3.
      end
    endcase

    ; Define arrays
      include_WLs            = where( (wl gt xr[0]) and (wl lt xr[1]), /NULL )                      ; wavelengths to use in the scattered light fitting
      plot_WLs               = where( (wl gt xr[0]) and (wl lt xr[1]), /NULL )                      ; wavelengths to actually plot
      residual_array         = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
      plot_residual_array    = fltarr(N_elements(plot_WLs), n_elements(Eclipse_files))
      LSF_Fit_array          = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
      WL_array               = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
      Brightness_array       = fltarr(n_elements(Eclipse_files))
      err_Brightness_array   = fltarr(n_elements(Eclipse_files))
      T_P_Shadow_array       = fltarr(n_elements(Eclipse_files))
      T_U_Shadow_array       = fltarr(n_elements(Eclipse_files))
      EXPTime_array          = fltarr(n_elements(Eclipse_files))
      DopplerShift_array1    = fltarr(n_elements(Eclipse_files))                                    ; For the first line in the triplet
      DopplerShift_array2    = fltarr(n_elements(Eclipse_files))
      DopplerShift_array3    = fltarr(n_elements(Eclipse_files))
      Gof                    = 0.                                                                   ; Goodness of fit (co-addded mean residual)
      GOF_array              = fltarr(N_elements(Fit_Me_To_The_Scatter))
      Use_Scatter            = strarr(N_elements(Eclipse_files))
      Possible_WL_offset     = fltarr(N_elements(Eclipse_files), N_elements(Fit_Me_To_The_Scatter)) ; Save any offsets from Io's expected rest wavelength, then correct any error in the wavelength solution.
      Possible_line_width    = fltarr(N_elements(Eclipse_files), N_elements(Fit_Me_To_The_Scatter)) ; Save the linewidths from all type of scatter subtraction, then average and peg this parameter.

    for i = 0, n_elements(Eclipse_files)-1 do begin
      if i eq 0 then continue ; skip penumbra
      spec        = MRDFITS(reduced_dir+'R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, header, /fscale, /unsigned, /silent)
      spec_err    = MRDFITS(reduced_dir+'sig_R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, sig_header, /fscale, /unsigned, /silent)
      WL          = sxpar(header, 'CRVAL1')+findgen(N_elements(spec))*sxpar(header, 'Cdelt1')

      for h = 0, n_elements(Fit_Me_To_The_Scatter)-1 do begin
        ; Which Jovian scatter are we using?
        Jup_Cal_Tell_Corr = MRDFITS(reduced_dir+Fit_Me_To_The_Scatter[h], 0, Jupiter_Scatter_header, /fscale, /unsigned, /silent)

        ; Define fitting indicies we need to avoid because of airglow itself.
        Ios_Airglow    = [Io_airglow + Io_airglow * sxpar(header, 'IO_DOPPL') / cspice_clight()]                                      ; wavelength locations of all airglow lines, Ignore telluric K
        ARCES_1sigma_bandwidth = (mean(xr)/31500.) / 2.3548 ; 1. * FWHM in Angstroms / 2sqrt(2ln2) = 68% of flux enclosed within 1 sigma --> Appropriate for weak telluric airglow
        Io_AG_Ind = [] & Telluric_AG_Ind = []
        for j = 0, N_elements(Ios_Airglow)-1 do Io_AG_ind = [Io_AG_Ind, where(abs(Ios_Airglow[j] - WL) lt ARCES_1sigma_bandwidth, /NULL)] ;  indicies within 2 sigma of an Io airglow line
        for j = 0, N_elements(Telluric_Airglow)-1 do Telluric_AG_ind = [Telluric_AG_Ind, where(abs(Telluric_Airglow[j] - WL) lt ARCES_1sigma_bandwidth, /NULL)] ; indicies within 1 sigma of any Telluric airglow
        AG_ind = [Io_AG_Ind, Telluric_AG_Ind]
        AG_ind = AG_ind[UNIQ(AG_ind, SORT(AG_ind))]

        ; Fit the Jovian scattered light, use Amoeba/Correlation to match Jupiter in terms of smoothing and shifting, then multiply and add until matches the scattered background
        correl_indicies  = cgSetDifference(include_Wls, AG_Ind)
        Mult_and_add     = mpfitfun('MX_Plus_B', Jup_Cal_Tell_Corr[correl_indicies], Spec[correl_indicies], Spec_err[correl_indicies], parinfo = MX_plus_B_parinfo,$
                                    [mean(Spec[correl_indicies]/Jup_Cal_Tell_Corr[correl_indicies]), 0.], /NaN, status = Scatter_fit_status, PERROR = err_a, /quiet)
        jup              = Mult_and_add[0]*Jup_Cal_Tell_Corr + Mult_and_add[1]

          ; Iterate Amoeba until it gives an answer
          Ftol = 1.e-4 &  Xtol = 1.e-4 & count = 0                                          ; Define initial tolerances for Amoeba
          smooth_and_shift = [2., 0.] & scale  = [2., 2.]                                   ; Define intial guess and simplex for Amoeba
          REPEAT BEGIN
            trial_smooth_and_shift = AMOEBAX(Ftol, xtol, function_name='shift_smooth', SCALE = Scale, P0 = smooth_and_shift, FUNCTION_VALUE = fval)
            if N_elements (N_elements(trial_smooth_and_shift) lt 2) then begin              ; If Amoeba crashes
              xtol = xtol * 1.5 & ftol = ftol * 1.5                                         ; Grow the tolerances
              smooth_and_shift = smooth_and_shift + [RANDOMN(seed), RANDOMN(seed)]          ; Change the guess
              SCALE            = SCALE + [RANDOMN(seed), RANDOMN(seed)]                     ; Change the simplex
              count = count+1
            endif
          ENDREP UNTIL ( (N_elements(trial_smooth_and_shift) eq 2) or (count gt 500) )
          smooth_and_shift = AMOEBAX(Ftol, xtol, function_name='shift_smooth', SCALE = [2., 2.], P0 = smooth_and_shift, FUNCTION_VALUE = fval)
          if N_elements(smooth_and_shift) eq 1 then stop ; check for bugs

        smoothed_shifted = interpolate(gauss_smooth(jup, smooth_and_shift[0]), findgen(n_elements(Jup)) - smooth_and_shift[1])
        Mult_and_add     = mpfitfun('MX_Plus_B', smoothed_shifted[correl_indicies], Spec[correl_indicies], Spec_err[correl_indicies], parinfo = MX_plus_B_parinfo, $
                                    [mean(Spec[correl_indicies]/smoothed_shifted[correl_indicies]), 0.], /NaN, status = Scatter_fit_status, PERROR = err_a, /quiet)
        scatter_fit      = Mult_and_add[0]*smoothed_shifted + Mult_and_add[1]
        residual         = Spec - scatter_fit
        GOF_array[h]     = GOF + median(abs(residual[correl_indicies]))

        ; Get the linewidth and the determine any systematic wavelength shift from Io's rest
        ; Ultimately we will peg both of these parameters.
        S3_Io_frame               = S3 + S3 * sxpar(header, 'IO_DOPPL') / cspice_clight()
        LSF_fitting_ind1          = cgSetDifference(where( abs(wl - S3_Io_frame) lt 0.3, /NULL), Telluric_AG_Ind)  ; fit region within +/- 0.3A of line center, excluding Telluric airglow
        initial_guess             = [2.e3, S3_Io_frame, 0.077]  ; rough values are fine here [height, wavelength, linewidth_sigma]
        fa                        = {x:double(WL[LSF_fitting_ind1]), y:double(residual[LSF_fitting_ind1]), err:double( sqrt(abs( residual[LSF_fitting_ind1] )) )}
        a                         = mpfit('Gaussian_for_MPFIT', initial_guess, PERROR = err_a, funct=fa, maxiter=50, STATUS = Did_it_work, /Quiet, NPEGGED = NPEGGED)
        Possible_WL_offset[i,h]   = a[1] - S3_Io_frame
        Possible_Line_Width[i,h]  = a[2]
      endfor ; h
      best_GoF = min(GOF_array, best_fit)
      Use_Scatter[i] = Fit_Me_To_The_Scatter[best_fit]
      print, 'Frame: ', Eclipse_files[i], ' is best paired with the Jupiter scatter from: ', Use_Scatter[i], ', Goodness of fit = ', best_GoF
    endfor
    Case date of
      'UT180320': begin
          best_scatter_spectrum = 'R_per_A_Jupiter_Center.ec.fits' ; Manually input desired scatting spectrum.
      end
      'UT190812': begin
          best_scatter_spectrum = 'R_per_A_Jupiter_Center.ec.fits' ; Manually input desired scatting spectrum.
      end
    endcase
    print, 'UPON INSPECTION, best choice seems to be: ', best_scatter_spectrum, 'Unnerving how Jovian scatter cannot verify a smoothed jovian center spectrum'

    ; -------------------------------------------Waterfall Plot Postscript Io and Fit Jovian scatter--------------------------------------------------------------------------------------------------------------

    ; Okay, now that we've established which Jupiter scatter spectrum best optimzes the residual, we can actually use it for scattered light subtraction

    ; get color versus ingress time & setup plot positions
      timeColors = BytScl(findgen(11), Min=0, Max=11, Top=11)
      Color      = timeColors[0]
      cgLoadCT, 33, NColors=8, /reverse
      pos        = cgLayout([1,3], OXMargin=[10,3], OYMargin=[8,5], YGap=0)
      Pos[1,0]   = Pos[1,0]*.8 & Pos[3,1] = Pos[3,1]*.8 
      Pos[1,1]   = Pos[1,1]*.84 & Pos[3,2] = Pos[3,2]*.84 

    cgPS_Open, filename = Reduced_Dir+'S9225_scatterfit.eps', /ENCAPSULATED, xsize = 7.5, ysize = 5.
    !P.font = 1
    device, SET_FONT = 'Helvetica Bold', /TT_FONT
    !p.charsize = 1.5

    cgplot, WL, spec/1.e3, psym = 0, Ytitle = 'kR / '+cgsymbol('Angstrom'), $             ; Rayleighs to KiloRayleighs
      title = 'Io''s Sulfur I 9225'+cgsymbol('Angstrom')+' Airglow in Eclipse',  $
      xr = xr, /nodata, pos = pos[*,0], xtickformat = '(A1)', yr = yr

    for i = 0, n_elements(Eclipse_files)-1 do begin
      if i eq 0 then continue ; skip penumbra
      spec         = MRDFITS(reduced_dir+'R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, header, /fscale, /unsigned, /silent)
      spec_err     = MRDFITS(reduced_dir+'sig_R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, sig_header, /fscale, /unsigned, /silent)
      spec         = spec / 1.e3       ; to KR
      spec_err     = spec_err / 1.e3   ; to KR
      WL           = sxpar(header, 'CRVAL1') + findgen(N_elements(spec))*sxpar(header, 'Cdelt1')
      S1_Io_frame  = S1 + S1 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      S2_Io_frame  = S2 + S2 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      S3_Io_frame  = S3 + S3 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      junk         = min( abs(WL - S1_Io_frame), S1_Io_Ind)
      junk         = min( abs(WL - S2_Io_frame), S2_Io_Ind)
      junk         = min( abs(WL - S3_Io_frame), S3_Io_Ind)

      ; Which Jovian scatter are we using?
        Jup_Cal_Tell_Corr = MRDFITS(reduced_dir+best_scatter_spectrum, 0, Jup_Scatter_header, /fscale, /unsigned, /silent)
        ;Jup_Cal_Tell_Corr = MRDFITS(reduced_dir+Use_Scatter[i], 0, Jup_Scatter_header, /fscale, /unsigned, /silent)

      ; Define fitting indicies we need to avoid because of airglow itself.
        Ios_Airglow    = [Io_airglow + Io_airglow * sxpar(header, 'IO_DOPPL') / cspice_clight()]                                      ; wavelength locations of all airglow lines, Ignore telluric K
        ARCES_1sigma_bandwidth = (mean(xr)/31500.) / 2.3548 ; FWHM in Angstroms / 2sqrt(2ln2) = 68% of flux enclosed within 1 sigma --> Appropriate for weak telluric airglow
        Io_AG_Ind = [] & Telluric_AG_Ind = []
        for j = 0, N_elements(Ios_Airglow)-1 do Io_AG_ind = [Io_AG_Ind, where(abs(Ios_Airglow[j] - WL) lt 2.*ARCES_1sigma_bandwidth, /NULL)] ;  indicies within 2 sigma of an Io airglow line
        for j = 0, N_elements(Telluric_Airglow)-1 do Telluric_AG_ind = [Telluric_AG_Ind, where(abs(Telluric_Airglow[j] - WL) lt ARCES_1sigma_bandwidth, /NULL)] ; indicies within 1 sigma of any Telluric airglow
        AG_ind = [Io_AG_Ind, Telluric_AG_Ind]
        AG_ind = AG_ind[UNIQ(AG_ind, SORT(AG_ind))]

      ; Fit the Jovian scattered light, use Amoeba/Correlation to match Jupiter in terms of smoothing and shifting, then multiply and add until matches the scattered background                                                                    ; peg the additive "B" component of the MX_Plus_B at zero
        correl_indicies  = cgSetDifference(include_Wls, AG_Ind)
        weights_9225     = 1./abs(S1_Io_Ind - correl_indicies)^(1.5) + 1./abs(S2_Io_Ind - correl_indicies)^(1.5) + 1./abs(S3_Io_Ind - correl_indicies)^(1.5)
        Mult_and_add     = mpfitfun('MX_Plus_B', Jup_Cal_Tell_Corr[correl_indicies], Spec[correl_indicies], Spec_err[correl_indicies], parinfo = MX_plus_B_parinfo, $
                                    [mean(Spec[correl_indicies]/Jup_Cal_Tell_Corr[correl_indicies]), 0.], /NaN, status = Scatter_fit_status, PERROR = err_a, /quiet)
        jup              = Mult_and_add[0]*Jup_Cal_Tell_Corr + Mult_and_add[1]

          ; Iterate Amoeba until it gives an answer
          Ftol = 1.e-4 &  Xtol = 1.e-4 & count = 0                                          ; Define initial tolerances for Amoeba
          smooth_and_shift = [2., 0.] & scale  = [2., 2.]                                   ; Define intial guess and simplex for Amoeba
          REPEAT BEGIN
            trial_smooth_and_shift = AMOEBAX(Ftol, xtol, function_name='shift_smooth', SCALE = Scale, P0 = smooth_and_shift, FUNCTION_VALUE = fval)
            if N_elements (N_elements(trial_smooth_and_shift) lt 2) then begin              ; If Amoeba crashes
              xtol = xtol * 1.5 & ftol = ftol * 1.5                                         ; Grow the tolerances
              smooth_and_shift = smooth_and_shift + [RANDOMN(seed), RANDOMN(seed)]          ; Change the guess
              SCALE            = SCALE + [RANDOMN(seed), RANDOMN(seed)]                     ; Change the simplex
              count = count+1
            endif
          ENDREP UNTIL ( (N_elements(trial_smooth_and_shift) eq 2) or (count gt 500) )
          smooth_and_shift = AMOEBAX(Ftol, xtol, function_name='shift_smooth', SCALE = [2., 2.], P0 = smooth_and_shift, FUNCTION_VALUE = fval)
          if N_elements(smooth_and_shift) eq 1 then stop ; check for bugs
  
        smoothed_shifted = interpolate(gauss_smooth(jup, smooth_and_shift[0]), findgen(n_elements(Jup)) - smooth_and_shift[1])
  
        Mult_and_add     = mpfitfun('MX_Plus_B', smoothed_shifted[correl_indicies], Spec[correl_indicies], Spec_err[correl_indicies], [mean(Spec[correl_indicies]/smoothed_shifted[correl_indicies]), 0.], $
          /NaN, status = Scatter_fit_status, PERROR = err_a, /quiet, weights = weights_9225, parinfo = MX_plus_B_parinfo) ; weights override error, local indices weighted
  
        scatter_fit      = Mult_and_add[0]*smoothed_shifted + Mult_and_add[1]
        residual         = Spec - scatter_fit
        residual_err     = make_array( N_elements(residual), value = stddev(residual[correl_indicies]) )

      ; plot the spectrum, fit and indices used for the fit (if we're in debug mode)
      ;cgplot, WL[correl_indicies], scatter_fit[correl_indicies], psym=14, /overplot
      cgplot, WL, scatter_fit, color = timeColors[i], linestyle = 1, thick = 7, /overplot
      cgplot, WL, spec, color = timeColors[i], thick = 5, /overplot

      ; Now do the line fitting to the residual. First define Line Spread Function (LSF) Gaussian fit parameter info
      parinfo               = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0]}, 3)
      parinfo[1].fixed      = 1
      parinfo[1].value      = double(S3_Io_frame)                                                           ; Pin the line's wavelength
      parinfo[2].fixed      = 1
      parinfo[2].value      = median(Possible_line_width[1:*,*])
      LSF_fitting_ind1      = cgSetDifference(where( abs(wl - S3_Io_frame) lt 0.3, /NULL), Telluric_AG_Ind) ; fit region within +/- 0.3A of line center, excluding Telluric airglow

      ; Gaussian fit
      initial_guess = [5.D, parinfo[1].value, parinfo[2].value]
      fa            = {x:double(WL[LSF_fitting_ind1]), y:double(residual[LSF_fitting_ind1]), err:double( sqrt(abs( spec[LSF_fitting_ind1] )) )}
      a             = mpfit('Gaussian_for_MPFIT', initial_guess, PERROR = err_a, funct=fa, parinfo = parinfo, STATUS = Did_it_work, /Quiet, NPEGGED = NPEGGED)

      ; write the results
      WL_array[*, i]           = WL[include_WLs]
      plot_residual_array[*, i]= residual[plot_WLs]
      LSF_Fit_array[*, i]      = gaussian(WL[include_WLs], a)                         ; LSF_Fit
      Brightness_array[i]      = A[0]*A[2]*SQRT(2*!DPI)                               ; Area under the Gaussian
      err_Brightness_array[i]  = stddev(residual[correl_indicies])*A[2]*SQRT(2*!DPI)  ; use the std dev of the residual and the width of the line
      T_P_Shadow_array[i]      = sxpar(header, 'T_PSHADO')
      T_U_Shadow_array[i]      = sxpar(header, 'T_USHADO')
      EXPTime_array[i]         = sxpar(sig_header, 'EXPTIME') / 60.                   ; for some reason, this keyword is bad in the regular headers, so use sig_header as a workaround
      DopplerShift_array1[i]   = S1_Io_frame
      DopplerShift_array2[i]   = S2_Io_frame
      DopplerShift_array3[i]   = S3_Io_frame
    endfor

    ; Plot the residual and Gaussian fits
      N_frames            = total(Brightness_array gt 0.)
      cgplot, spec, WL, psym = 0, xtickformat = '(A1)', yminor = 2, $
        yr = [-2.99, 5.5 ], xr = XR, /nodata, pos = pos[*,1], /noerase      
      cgtext, xr[0] - 2., -2., 'Residual (kR / '+cgsymbol('Angstrom')+')', orientation = 90, alignment = 0.5
      
    ; Co-align to Io's Doppler Shift ???
      shift_By = mean( [transpose( DopplerShift_array1 - Mean(DopplerShift_array1[1:*]) ), $
        transpose( DopplerShift_array2 - Mean(DopplerShift_array2[1:*]) ), $
        transpose( DopplerShift_array3 - Mean(DopplerShift_array3[1:*]) )], dimension = 1 )
      aligned_residuals = plot_residual_array

    for i = 1, n_elements(Eclipse_files)-1 do begin
      cgplot, WL[plot_WLs], plot_residual_array[*, i], /OVERPLOT, COLOR = timeColors[i], psym=10
      cgplot, [DopplerShift_array1[i], DopplerShift_array1[i]], [5.5, 3.9], COLOR = timeColors[i], /overplot
      cgplot, [DopplerShift_array2[i], DopplerShift_array2[i]], [5.5, 4.7], COLOR = timeColors[i], /overplot
      cgplot, [DopplerShift_array3[i], DopplerShift_array3[i]], [5.5, 4.5], COLOR = timeColors[i], /overplot
      aligned_residuals[*,i] = interpol( plot_residual_array[*,i], WL[plot_WLs], WL[plot_WLs] - Shift_By[i])
    endfor
    co_add = total(aligned_residuals[*,1:*], 2) / n_frames         
      
      ; plot the aligned & co-added spectra
        cgplot, WL[plot_WLs], co_add, pos = pos[*,2], thick=5, /noerase, xr = xr, yr = [-1.9, 3.2], /yst,  $
          xtitle = 'Wavelength ('+cgsymbol('Angstrom')+')', yminor = 2                
        cgplot, xr, [0.,0.], linestyle = 2, /overplot
        cgplot, [Mean(DopplerShift_array1[1:*]), Mean(DopplerShift_array1[1:*])], [-1.9, 3.2], linestyle = 1, /overplot
        cgplot, [Mean(DopplerShift_array2[1:*]), Mean(DopplerShift_array2[1:*])], [-1.9, 3.2], linestyle = 1, /overplot
        cgplot, [Mean(DopplerShift_array3[1:*]), Mean(DopplerShift_array3[1:*])], [-1.9, 3.2], linestyle = 1, /overplot
    cgps_close
    
  ; Estimate the total time-averaged brightness of the emission feature or features
    N_frames            = total(Brightness_array gt 0.)
    integrate_multiplet = cgSetIntersection(Io_AG_Ind, plot_WLs) - min(plot_WLs) ;indices of the multiplet over
    print, 'Disk- and Time-Averaged 9225A Brightness =', total(co_add[integrate_multiplet]) * N_elements(integrate_multiplet)*sxpar(header, 'Cdelt1') / N_frames, ' KiloRayleighs'  
    window, 0
    x = WL[plot_WLs]
    cgplot, x[integrate_multiplet], co_add[integrate_multiplet] / N_frames, /xst, psym=16, ytitle = 'Disk- and Time-Averaged 9225A Brightness (kR/A)'
    cgplot, x[integrate_multiplet], co_add[integrate_multiplet] / N_frames, /overplot, psym=10
  endif

  if Part eq 1.2 then begin ; Make waterfall plots, Na-D does not need telluric corrected spectra.

    ;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ;--------------------------------------------------------Get the Na Time Series and Scatter Fit-------------------------------------------------------------------------------
    ;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ; METHOD:
    ;   Loop though all possible spectra that could be fit to jovian scatter and choose what works best. Avoid Io and telluric airglow in the fitting.
    ;   Three steps to the fit: 1) Use Mpfit to roughly add & multiply the Jovian scatter spectra to approach the background in the Io exposures
    ;                           2) Use Amoeba to shift and smooth the spectra, maximizing correlation
    ;                           3) Use Mpfit again to fine tune the add & multiply of Jovian scatter spectra needed to match the background in the Io exposures
    ; NOTES:
    ;   This section can be a little fiddly. Amoeba is prone to infinite loops if two points in the simplex give the same correlation.
    ;   Therefore, include_WLs needs to be careful tuned to keep Amoeba stable.
    ;   It's critical to inspect with plot, WL[correl_indicies], scatter_fit[correl_indicies], uncomment this line to inspect and adjust things
    ;   With this inspection, it is useful to adjust *include_WLs*, and *AG_ind* to cover regions of continuum as close as possible to Io's airglow
    ;   Parameter information is not included in the Gaussian fitting to the line spread function in the residual. I couldn't figgure out why it failed with status = 0


    ; Multiplicative and additive DC offset allowed
      MX_plus_B_parinfo          = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0]}, 2)

    ; Fit various spectra to the Jovian Scatter, *use* whichever one produces the lowest residual. Note this method will never work if there's much telluric absorption.
    Case date of
      'UT180320': begin
        Fit_Me_To_The_Scatter    = ['R_per_A_Jupiter_Center.ec.fits', $
                                  'R_per_A_Jovian_Scatter.0001.ec.fits', $
                                  'R_per_A_Jovian_Scatter.0002.ec.fits', $
                                  'R_per_A_Jovian_Scatter.0003.ec.fits' ]
                  end
      'UT190812': begin
        Fit_Me_To_The_Scatter    = ['R_per_A_Jupiter_Center.ec.fits']
                  end
    endcase

    ; setup the plot axis
      spec = MRDFITS(reduced_dir+'R_per_A_'+Eclipse_files[1]+'.ec.fits', 0, header, /fscale, /unsigned, /silent)
      WL   = sxpar(header, 'CRVAL1')+findgen(N_elements(spec))*sxpar(header, 'Cdelt1')

      Case date of
        'UT180320': begin
          YR   = [18,104]
          XR   = [5888.,5897.5]
        end
        'UT190812': begin
          yr   = [5.5e1, 2e2]
          YR_resid = [-2, 20]
          XR   = [5888.,5897.5]
        end
      endcase

    ; Define arrays
      include_WLs            = where( (abs(wl - Na1+.4) lt 1.1) or (abs(wl - Na2+.4) lt 1.1), /NULL ); wavelength regions over which we'll fit the jovian scatter
      plot_WLs               = where( (wl gt xr[0]) and (wl lt xr[1]), /NULL )                       ; wavelengths to actually plot
      residual_array         = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
      plot_residual_array    = fltarr(N_elements(plot_WLs), n_elements(Eclipse_files))
      LSF_Fit_array          = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
      WL_array               = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
      Brightness_array       = fltarr(n_elements(Eclipse_files))
      err_Brightness_array   = fltarr(n_elements(Eclipse_files))
      T_P_Shadow_array       = fltarr(n_elements(Eclipse_files))
      T_U_Shadow_array       = fltarr(n_elements(Eclipse_files))
      EXPTime_array          = fltarr(n_elements(Eclipse_files))
      DopplerShift_array1    = fltarr(n_elements(Eclipse_files))                          ; For the D2 line
      DopplerShift_array2    = fltarr(n_elements(Eclipse_files))                          ; For the D1 line
      Gof                    = 0.                                                         ; Goodness of fit (co-addded mean residual)
      GOF_array              = fltarr(N_elements(Fit_Me_To_The_Scatter))
      Use_Scatter            = strarr(N_elements(Eclipse_files))
      Na1Na2_Fit_params      = REPLICATE({params, brightness:0.0, err_brightness:0.0, mult:0.0, smooth:0.0, shift:0.0, background:'', T_P_shadow:0.0, exptime:0.0}, n_elements(Eclipse_files))

    for i = 0, n_elements(Eclipse_files)-1 do begin
      if i eq 0 then continue ; skip penumbra
      spec        = MRDFITS(reduced_dir+'R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, header, /fscale, /unsigned, /silent)
      spec_err    = MRDFITS(reduced_dir+'sig_R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, sig_header, /fscale, /unsigned, /silent)
      WL          = sxpar(header, 'CRVAL1')+findgen(N_elements(spec))*sxpar(header, 'Cdelt1')

      Na1_Io_frame = Na1 + Na1 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      Na2_Io_frame = Na2 + Na2 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      junk         = min( abs(WL - Na1_Io_frame), Na1_Io_Ind)
      junk         = min( abs(WL - Na2_Io_frame), Na2_Io_Ind)

      for h = 0, n_elements(Fit_Me_To_The_Scatter)-1 do begin
        ; Which Jovian scatter are we using?
          Jup_Cal_Tell_Corr = MRDFITS(reduced_dir+Fit_Me_To_The_Scatter[h], 0, header, /fscale, /unsigned, /silent)

        ; Define fitting indicies we need to avoid because of airglow itself.
          Ios_Airglow    = [Io_airglow + Io_airglow * sxpar(header, 'IO_DOPPL') / cspice_clight()]                                      ; wavelength locations of all airglow lines, Ignore telluric Na
          ARCES_1sigma_bandwidth = 1. * (mean(xr)/31500.) / 2.3548 ; 1. * FWHM in Angstroms / 2sqrt(2ln2) = 68% of flux enclosed within 1 sigma --> Appropriate for weak telluric airglow
          Io_AG_Ind = [] & Telluric_AG_Ind = []
          for j = 0, N_elements(Ios_Airglow)-1 do Io_AG_ind = [Io_AG_Ind, where(abs(Ios_Airglow[j] - WL) lt 2.*ARCES_1sigma_bandwidth, /NULL)] ;  indicies within 2 sigma of an Io airglow line
          for j = 0, N_elements(Telluric_Airglow)-1 do Telluric_AG_ind = [Telluric_AG_Ind, where(abs(Telluric_Airglow[j] - WL) lt ARCES_1sigma_bandwidth, /NULL)] ; indicies within 1 sigma of any Telluric airglow
          AG_ind = [Io_AG_Ind, Telluric_AG_Ind]
          AG_ind = AG_ind[UNIQ(AG_ind, SORT(AG_ind))]

        ; Fit the Jovian scattered light, use Amoeba/Correlation to match Jupiter in terms of smoothing and shifting, then multiply and add until matches the scattered background
          correl_indicies  = cgSetDifference(include_Wls, AG_Ind)
          Mult_and_add     = mpfitfun('MX_Plus_B', Jup_Cal_Tell_Corr[correl_indicies], Spec[correl_indicies], Spec_err[correl_indicies], parinfo = MX_plus_B_parinfo, $
                                      [mean(Spec[correl_indicies]/Jup_Cal_Tell_Corr[correl_indicies]), 0.], /NaN, status = Scatter_fit_status, PERROR = err_a, /quiet)
          jup              = Mult_and_add[0]*Jup_Cal_Tell_Corr + Mult_and_add[1]
  
          ; Iterate Amoeba until it gives an answer
            Ftol = 1.e-4 &  Xtol = 1.e-4 & count = 0        ; Define initial tolerances for Amoeba 
            smooth_and_shift = [2., 0.] & scale  = [2., 2.] ; Define nitial guess and simplex for Amoeba 
            REPEAT BEGIN
              trial_smooth_and_shift = AMOEBAX(Ftol, xtol, function_name='shift_smooth', SCALE = Scale, P0 = smooth_and_shift, FUNCTION_VALUE = fval)
              if N_elements (N_elements(trial_smooth_and_shift) lt 2) then begin              ; If Amoeba crashes
                xtol = xtol * 1.5 & ftol = ftol * 1.5                                         ; Grow the tolerances
                smooth_and_shift = smooth_and_shift + [RANDOMN(seed), RANDOMN(seed)]          ; Change the guess
                SCALE            = SCALE + [RANDOMN(seed), RANDOMN(seed)]                     ; Change the simplex
                count = count+1
              endif  
            ENDREP UNTIL ( (N_elements(trial_smooth_and_shift) eq 2) or (count gt 500) )
            smooth_and_shift = AMOEBAX(Ftol, xtol, function_name='shift_smooth', SCALE = [2., 2.], P0 = trial_smooth_and_shift, FUNCTION_VALUE = fval)
            if N_elements(smooth_and_shift) eq 1 then stop ; check for bugs

          smoothed_shifted = interpolate(gauss_smooth(jup, smooth_and_shift[0]), findgen(n_elements(Jup)) - smooth_and_shift[1])
          Mult_and_add     = mpfitfun('MX_Plus_B', smoothed_shifted[correl_indicies], Spec[correl_indicies], Spec_err[correl_indicies], parinfo = MX_plus_B_parinfo, $
                                      [mean(Spec[correl_indicies]/smoothed_shifted[correl_indicies]), 0.], /NaN, status = Scatter_fit_status, PERROR = err_a, /quiet)
                                      
          scatter_fit      = Mult_and_add[0]*smoothed_shifted + Mult_and_add[1]
          residual         = Spec - scatter_fit
          GOF_array[h]     = GOF + mean(abs(residual[correl_indicies]))
        ;print, 'Na D Jupiter Scatter Goodness of Fit = ', GOF_array[h]
      endfor
      best_GoF = min(GOF_array, best_fit)
      Use_Scatter[i] = Fit_Me_To_The_Scatter[best_fit]
      print, 'Frame: ', Eclipse_files[i], ' should take the Jupiter scatter from: ', Use_Scatter[i], ', Goodness of fit = ', best_GoF
    endfor
    ; -------------------------------------------Waterfall Plot Postscript Io and Fit Jovian scatter--------------------------------------------------------------------------------------------------------------

    ; Okay, now that we've established which Jupiter scatter spectrum best optimzes the residual, we can actually user it and do this scattered light subtraction

    ; get color versus ingress time
      timeColors = BytScl(findgen(11), Min=0, Max=11, Top=11)
      Color=timeColors[0]
      cgLoadCT, 33, NColors=8, /reverse

      pos = cgLayout([1,2], OXMargin=[10,3], OYMargin=[8,5], YGap=0)
      Pos[1,0] = Pos[1,0]*.7 & Pos[3,1] = Pos[3,1]*.7
      cgPS_Open, filename = Reduced_Dir+'Na-D_scatterfit.eps', /ENCAPSULATED, xsize = 7.5, ysize = 5.
      !P.font = 1
      device, SET_FONT = 'Helvetica Bold', /TT_FONT
      !p.charsize = 1.5

    cgplot, WL, spec/1.e3, psym = 0, Ytitle = 'kR / '+cgsymbol('Angstrom'), $             ; Rayleighs to KiloRayleighs
      title = 'Io''s Sodium Response to Eclipse',  $
      xr = xr, /nodata, pos = pos[*,0], xtickformat = '(A1)', yr = yr

    for i = 0, n_elements(Eclipse_files)-1 do begin
      if i eq 0 then continue ; skip penumbra
      spec         = MRDFITS(reduced_dir+'R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, header, /fscale, /unsigned, /silent)
      spec_err     = MRDFITS(reduced_dir+'sig_R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, sig_header, /fscale, /unsigned, /silent)
      spec         = spec / 1.e3       ; to KR
      spec_err     = spec_err / 1.e3   ; to KR
      WL           = sxpar(header, 'CRVAL1')+findgen(N_elements(spec))*sxpar(header, 'Cdelt1')
      Na1_Io_frame = Na1 + Na1 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      Na2_Io_frame = Na2 + Na2 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      junk         = min( abs(WL - Na1_Io_frame), Na1_Io_Ind)
      junk         = min( abs(WL - Na2_Io_frame), Na2_Io_Ind)
      
      ; Which Jovian scatter are we using?
        Jup_Cal_Tell_Corr = MRDFITS(reduced_dir+Use_Scatter[i], 0, Jup_Scatter_header, /fscale, /unsigned, /silent)

      ; Define fitting indicies we need to avoid because of airglow itself.
        Ios_Airglow    = [Io_airglow + Io_airglow * sxpar(header, 'IO_DOPPL') / cspice_clight()]                                      ; wavelength locations of all airglow lines, Ignore telluric Na
        ARCES_1sigma_bandwidth = 1. * (mean(xr)/31500.) / 2.3548 ; 1. * FWHM in Angstroms / 2sqrt(2ln2) = 68% of flux enclosed within 1 sigma --> Appropriate for weak telluric airglow
        Io_AG_Ind = [] & Telluric_AG_Ind = []
        for j = 0, N_elements(Ios_Airglow)-1 do Io_AG_ind = [Io_AG_Ind, where(abs(Ios_Airglow[j] - WL) lt 2.*ARCES_1sigma_bandwidth, /NULL)] ;  indicies within 2 sigma of an Io airglow line
        for j = 0, N_elements(Telluric_Airglow)-1 do Telluric_AG_ind = [Telluric_AG_Ind, where(abs(Telluric_Airglow[j] - WL) lt ARCES_1sigma_bandwidth, /NULL)] ; indicies within 1 sigma of any Telluric airglow
        AG_ind = [Io_AG_Ind, Telluric_AG_Ind]
        AG_ind = AG_ind[UNIQ(AG_ind, SORT(AG_ind))]

      ; Fit the Jovian scattered light, use Amoeba/Correlation to match Jupiter in terms of smoothing and shifting, then multiply and add until matches the scattered background
        correl_indicies  = cgSetDifference(include_Wls, AG_Ind)
        weights          = 1./abs(Na1_Io_Ind - correl_indicies)^(1.5) + 1./abs(Na2_Io_Ind - correl_indicies)^(1.5)
        Mult_and_add     = mpfitfun('MX_Plus_B', Jup_Cal_Tell_Corr[correl_indicies], Spec[correl_indicies], Spec_err[correl_indicies], parinfo = MX_plus_B_parinfo, $
                            [mean(Spec[correl_indicies]/Jup_Cal_Tell_Corr[correl_indicies]), 0.], /NaN, status = Scatter_fit_status, PERROR = err_a, /quiet)
        jup              = Mult_and_add[0]*Jup_Cal_Tell_Corr + Mult_and_add[1]

        ; Iterate Amoeba until it gives an answer
          Ftol = 1.e-4 &  Xtol = 1.e-4 & count = 0        ; Define initial tolerances for Amoeba
          smooth_and_shift = [3., 0.] & scale  = [3., 2.] ; Define nitial guess and simplex for Amoeba
          REPEAT BEGIN
            trial_smooth_and_shift = AMOEBAX(Ftol, xtol, function_name='shift_smooth', SCALE = Scale, P0 = smooth_and_shift, FUNCTION_VALUE = fval)
            if N_elements (N_elements(trial_smooth_and_shift) lt 2) then begin              ; If Amoeba crashes
              xtol = xtol * 1.5 & ftol = ftol * 1.5                                         ; Grow the tolerances
              smooth_and_shift = smooth_and_shift + [RANDOMN(seed), RANDOMN(seed)]          ; Change the guess
              SCALE            = SCALE + [RANDOMN(seed), RANDOMN(seed)]                     ; Change the simplex
              count = count+1
            endif
          ENDREP UNTIL ( (N_elements(trial_smooth_and_shift) eq 2) or (count gt 500) )
          smooth_and_shift = AMOEBAX(Ftol, xtol, function_name='shift_smooth', SCALE = [2., 2.], P0 = trial_smooth_and_shift, FUNCTION_VALUE = fval)
          if N_elements(smooth_and_shift) eq 1 then stop ; check for bugs
        
        smoothed_shifted = interpolate(gauss_smooth(jup, smooth_and_shift[0]), findgen(n_elements(Jup)) - smooth_and_shift[1])
        Mult_and_add     = mpfitfun('MX_Plus_B', smoothed_shifted[correl_indicies], Spec[correl_indicies], Spec_err[correl_indicies], [mean(Spec[correl_indicies]/smoothed_shifted[correl_indicies]), 0.], $
                                    /NaN, status = Scatter_fit_status, PERROR = err_a, /quiet, weights = weights, parinfo = MX_plus_B_parinfo)
        scatter_fit      = Mult_and_add[0]*smoothed_shifted + Mult_and_add[1]
        residual         = Spec - scatter_fit
        residual_err     = make_array( N_elements(residual), value = stddev(residual[correl_indicies]) )

      ; plot the spectrum, fit and indices used for the fit (if we're in debug mode)
        ;cgplot, WL[correl_indicies], scatter_fit[correl_indicies] + 60.-(i*12.), psym=14, /overplot
        cgplot, WL, scatter_fit + 60.-(i*12.), color = timeColors[i], linestyle = 1, thick = 5, /overplot
        cgplot, WL, spec + 60. -(i*12.), color = timeColors[i], thick = 5, /overplot

      LSF_fitting_ind1         = where( abs(wl - Na1_Io_frame) lt 0.3, /NULL) 
      LSF_fitting_ind2         = where( abs(wl - Na2_Io_frame) lt 0.3, /NULL) 
      LSF_Fit1                 = mpfitpeak(WL[LSF_fitting_ind1], residual[LSF_fitting_ind1], a, STATUS = Line_fit_Status1, /POSITIVE, nterms = 3)
      LSF_Fit2                 = mpfitpeak(WL[LSF_fitting_ind2], residual[LSF_fitting_ind2], b, STATUS = Line_fit_Status2, /POSITIVE, nterms = 3)
      WL_array[*, i]           = WL[include_WLs]
      plot_residual_array[*, i]= residual[plot_WLs]
      LSF_Fit_array[*, i]      = gaussian(WL[include_WLs], a) + gaussian(WL[include_WLs], b)                ; LSF_Fit
      Na1Na2_fit_params[i]     = {params, A[0]*A[2]*SQRT(2*!DPI) + B[0]*B[2]*SQRT(2*!DPI), stddev(residual[correl_indicies])*A[2]*SQRT(2*!DPI) + stddev(residual[correl_indicies])*B[2]*SQRT(2*!DPI), $
                                          Mult_and_add[0], smooth_and_shift[0], smooth_and_shift[1], Use_Scatter[i], float(sxpar(header, 'T_PSHADO')), float(sxpar(sig_header, 'EXPTIME')) / 60.}      
      Brightness_array[i]      = A[0]*A[2]*SQRT(2*!DPI) + B[0]*B[2]*SQRT(2*!DPI)                            ; Area under the Gaussian
      err_Brightness_array[i]  = stddev(residual[correl_indicies])*A[2]*SQRT(2*!DPI) + stddev(residual[correl_indicies])*B[2]*SQRT(2*!DPI)  ; use the std dev of the residual and the width of the line
      T_P_Shadow_array[i]      = sxpar(header, 'T_PSHADO')
      T_U_Shadow_array[i]      = sxpar(header, 'T_USHADO')
      EXPTime_array[i]         = sxpar(sig_header, 'EXPTIME') / 60.                                         ; for some reason, this keyword is bad in the regular headers, so use sig_header as a workaround
      DopplerShift_array1[i]   = Na1_Io_frame
      DopplerShift_array2[i]   = Na2_Io_frame
    
      print, 'Na D1 FWHM linewidths =' + string(a[2]*2.3548) + 'A or' + string(a[2]/Na1*cspice_clight()) + ' km/s'
      print, 'Na D2 FWHM linewidths =' + string(b[2]*2.3548) + 'A or' + string(b[2]/Na2*cspice_clight()) + ' km/s'
      print, 'D2 ='+string(A[0]*A[2]*SQRT(2*!DPI))+', D1 ='+string(B[0]*B[2]*SQRT(2*!DPI))+', D2/D1 line ratio =', (A[0]*A[2]*SQRT(2*!DPI)) / (B[0]*B[2]*SQRT(2*!DPI))
    endfor

    ; Residual plot
      cgplot, spec, WL, psym = 0, Ytitle = 'Residual (kR / '+cgsymbol('Angstrom')+')', xtitle = 'Wavelength ('+cgsymbol('Angstrom')+')', $
        yr = [-4., 16], xr = xr, /nodata, pos = pos[*,1], /noerase
      for i = 1, n_elements(Eclipse_files)-1 do begin
        cgplot, WL[plot_WLs], plot_residual_array[*, i], /OVERPLOT, COLOR = timeColors[i], psym=10
        cgplot, WL_array[*, i], LSF_Fit_array[*, i], /OVERPLOT, COLOR = timeColors[i]
        cgplot, [Na1, Na1], [16, 4.], linestyle = 1, /overplot
        cgplot, [Na2, Na2], [16, 3.], linestyle = 1, /overplot
        cgplot, [DopplerShift_array1[i], DopplerShift_array1[i]], [16, 14.5], COLOR = timeColors[i], /overplot
        cgplot, [DopplerShift_array2[i], DopplerShift_array2[i]], [16, 10.5], COLOR = timeColors[i], /overplot
      endfor
    cgps_close
    save, Na1Na2_fit_params, filename = Reduced_Dir+'Na1Na2_fit_params.sav'
endif
if Part eq 1.21 then begin ; Make waterfall plots, Na in the near-IR is crippled by telluric absorption---beware!

    ;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ;--------------------------------------------------------Get the Near-IR Na Time Series and Scatter Fit-------------------------------------------------------------------------------
    ;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ; METHOD:
    ;   Loop though all possible spectra that could be fit to jovian scatter and choose what works best. Avoid Io and telluric airglow in the fitting.
    ;   Three steps to the fit: 1) Use Mpfit to roughly add & multiply the Jovian scatter spectra to approach the background in the Io exposures
    ;                           2) Use Amoeba to shift and smooth the spectra, maximizing correlation
    ;                           3) Use Mpfit again to fine tune the add & multiply of Jovian scatter spectra needed to match the background in the Io exposures
    ; NOTES:
    ;   This section can be a little fiddly. Amoeba is prone to infinite loops if two points in the simplex give the same correlation.
    ;   Therefore, include_WLs needs to be careful tuned to keep Amoeba stable.
    ;   It's critical to inspect with plot, WL[correl_indicies], scatter_fit[correl_indicies], uncomment this line to inspect and adjust things
    ;   With this inspection, it is useful to adjust *include_WLs*, and *AG_ind* to cover regions of continuum as close as possible to Io's airglow
    ;   Parameter information is not included in the Gaussian fitting to the line spread function in the residual. I couldn't figgure out why it failed with status = 0

      restore, Reduced_Dir+'Na1Na2_fit_params.sav'      
      ;restore, Reduced_Dir+'O2_fit_params.sav'   
      ;Na1Na2_fit_params = O2_fit_params

    ; load telluric absorption
      telluric_absorp       = MRDFITS(reduced_dir + 'master_telluric_spectrum.fits', 0, Junkheader, /fscale, /silent)

    ; Define local weighting
      local  = 1.2 ; 0.0 to 2.0 bigger = weight the fits more locally to the Io emission's wavelength

    ; Fit various spectra to the Jovian Scatter, *use* whichever one produces the lowest residual. Note this method will never work if there's much telluric absorption.
      Fit_Me_To_The_Scatter = [ 'R_per_A_Jupiter_Center.ec.fits', $
                                'R_per_A_Jovian_Scatter.0001.ec.fits', $
                                'R_per_A_Jovian_Scatter.0002.ec.fits', $
                                'R_per_A_Jovian_Scatter.0003.ec.fits' ]

    ; setup the plot axis
      spec = MRDFITS(reduced_dir+'R_per_A_'+Eclipse_files[1]+'.ec.fits', 0, header, /fscale, /unsigned, /silent)
      WL   = sxpar(header, 'CRVAL1')+findgen(N_elements(spec))*sxpar(header, 'Cdelt1')
      Case date of
        'UT180320': begin
          YR   = [15,41]
          XR   = [8181.,8196.]
        end
        'UT190812': yr   = [1.1e1, 4e1]
      endcase

    ; Define arrays
      include_WLs            = where( ((wl gt 8193.75) and (wl lt 8195.25)) or $                     ; best yet for no telluric correction 
                                      ((wl gt 8182.25) and (wl lt 8183.75)) , /NULL )                ; best yet for no telluric correction 
      plot_WLs               = where( (wl gt xr[0]) and (wl lt xr[1]), /NULL )                       ; wavelengths to actually plot
      residual_array         = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
      plot_residual_array    = fltarr(N_elements(plot_WLs), n_elements(Eclipse_files))
      LSF_Fit_array          = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
      WL_array               = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
      Brightness_array       = fltarr(n_elements(Eclipse_files))
      err_Brightness_array   = fltarr(n_elements(Eclipse_files))
      T_P_Shadow_array       = fltarr(n_elements(Eclipse_files))
      T_U_Shadow_array       = fltarr(n_elements(Eclipse_files))
      EXPTime_array          = fltarr(n_elements(Eclipse_files))
      DopplerShift_array1    = fltarr(n_elements(Eclipse_files))                          ; For the D2 line
      DopplerShift_array2    = fltarr(n_elements(Eclipse_files))                          ; For the D1 line
      Gof                    = 0.                                                         ; Goodness of fit (co-addded mean residual)
      GOF_array              = fltarr(N_elements(Fit_Me_To_The_Scatter))
      Use_Scatter            = strarr(N_elements(Eclipse_files))

    ; determine indices where telluric absorption is too severe to recover
      telluric_absorp_ind = where(telluric_absorp lt 0.94)                                ; don't use regions of 10% absorption or more when fitting jovian scatter

    for i = 0, n_elements(Eclipse_files)-1 do begin
      if i eq 0 then continue ; skip penumbra
      spec        = MRDFITS(reduced_dir+'R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, header, /fscale, /unsigned, /silent)
      spec_err    = MRDFITS(reduced_dir+'sig_R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, sig_header, /fscale, /unsigned, /silent)
      WL          = sxpar(header, 'CRVAL1')+findgen(N_elements(spec))*sxpar(header, 'Cdelt1')
      Na3_Io_frame = Na3 + Na3 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      Na4_Io_frame = Na4 + Na4 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      junk         = min( abs(WL - Na3_Io_frame), Na3_Io_Ind)
      junk         = min( abs(WL - Na4_Io_frame), Na4_Io_Ind)

      for h = 0, n_elements(Fit_Me_To_The_Scatter)-1 do begin
        ; Which Jovian scatter are we using?
          Jup_Cal_Tell_Corr = MRDFITS(reduced_dir+Fit_Me_To_The_Scatter[h], 0, header, /fscale, /unsigned, /silent)

        ; Define fitting indicies we need to avoid because of airglow itself.
          Ios_Airglow    = [Io_airglow + Io_airglow * sxpar(header, 'IO_DOPPL') / cspice_clight()]                                      ; wavelength locations of all airglow lines, Ignore telluric Na
          ARCES_1sigma_bandwidth = 1. * (mean(xr)/31500.) / 2.3548 ; 1. * FWHM in Angstroms / 2sqrt(2ln2) = 68% of flux enclosed within 1 sigma --> Appropriate for weak telluric airglow
          Io_AG_Ind = [] & Telluric_AG_Ind = []
          for j = 0, N_elements(Ios_Airglow)-1 do Io_AG_ind = [Io_AG_Ind, where(abs(Ios_Airglow[j] - WL) lt 1.5*ARCES_1sigma_bandwidth, /NULL)] ;  indicies within 2 sigma of an Io airglow line
          for j = 0, N_elements(Telluric_Airglow)-1 do Telluric_AG_ind = [Telluric_AG_Ind, where(abs(Telluric_Airglow[j] - WL) lt ARCES_1sigma_bandwidth, /NULL)] ; indicies within 1 sigma of any Telluric airglow
          AG_ind = [Io_AG_Ind, Telluric_AG_Ind]
          AG_ind = AG_ind[UNIQ(AG_ind, SORT(AG_ind))]

        ; Fit the Jovian scattered light, use Amoeba/Correlation to match Jupiter in terms of smoothing and shifting, then multiply and add until matches the scattered background
          correl_indicies  = cgSetDifference(include_Wls, AG_Ind)
          correl_indicies  = cgSetDifference(correl_indicies, telluric_absorp_ind)
          parinfo          = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0]}, 2)
          parinfo[1].fixed = 1
          weights          = 1./abs(Na3_Io_Ind - correl_indicies)^(local) + 1./abs(Na4_Io_Ind - correl_indicies)^(local)                                                                        ; peg the additive "B" component of the MX_Plus_B at zero
          Mult_and_add     = mpfitfun('MX_Plus_B', Jup_Cal_Tell_Corr[correl_indicies], Spec[correl_indicies], Spec_err[correl_indicies], parinfo = parinfo, $
                                      [mean(Spec[correl_indicies]/Jup_Cal_Tell_Corr[correl_indicies]), 0.], /NaN, status = Scatter_fit_status, weights = weights, PERROR = err_a, /quiet)
          jup              = Mult_and_add[0]*Jup_Cal_Tell_Corr + Mult_and_add[1]
  
          ; Iterate Amoeba until it gives an answer
            Ftol = 1.e-4 &  Xtol = 1.e-4 & count = 0                                          ; Define initial tolerances for Amoeba 
            smooth_and_shift = [2., 0.] & scale  = [2., 2.]                                   ; Define intial guess and simplex for Amoeba 
            REPEAT BEGIN
              trial_smooth_and_shift = AMOEBAX(Ftol, xtol, function_name='shift_smooth', SCALE = Scale, P0 = smooth_and_shift, FUNCTION_VALUE = fval)
              if N_elements (N_elements(trial_smooth_and_shift) lt 2) then begin              ; If Amoeba crashes
                xtol = xtol * 1.5 & ftol = ftol * 1.5                                         ; Grow the tolerances
                smooth_and_shift = smooth_and_shift + [RANDOMN(seed), RANDOMN(seed)]          ; Change the guess
                SCALE            = SCALE + [RANDOMN(seed), RANDOMN(seed)]                     ; Change the simplex
                count = count+1
              endif  
            ENDREP UNTIL ( (N_elements(trial_smooth_and_shift) eq 2) or (count gt 500) )
            smooth_and_shift = AMOEBAX(Ftol, xtol, function_name='shift_smooth', SCALE = [2., 2.], P0 = smooth_and_shift, FUNCTION_VALUE = fval)
            if N_elements(smooth_and_shift) eq 1 then stop ; check for bugs
         
          smoothed_shifted = interpolate(gauss_smooth(jup, smooth_and_shift[0]), findgen(n_elements(Jup)) - smooth_and_shift[1])
          Mult_and_add     = mpfitfun('MX_Plus_B', smoothed_shifted[correl_indicies], Spec[correl_indicies], Spec_err[correl_indicies], parinfo = parinfo, $
                                      [mean(Spec[correl_indicies]/smoothed_shifted[correl_indicies]), 0.], /NaN, weights = weights, status = Scatter_fit_status, PERROR = err_a, /quiet)
        scatter_fit      = Mult_and_add[0]*smoothed_shifted + Mult_and_add[1]
        residual         = Spec - scatter_fit
        GOF_array[h]     = GOF + mean(abs(residual[correl_indicies]))

      endfor
      best_GoF = min(GOF_array, best_fit)
      Use_Scatter[i] = Fit_Me_To_The_Scatter[best_fit]
      print, 'Frame: ', Eclipse_files[i], ' should take the Jupiter scatter from: ', Use_Scatter[i], ', Goodness of fit = ', best_GoF
    endfor
    ;best_scatter_spectrum = 'R_per_A_Jupiter_Center.ec.fits' ; Manually input desired scattering spectrum.
    best_scatter_spectrum = 'R_per_A_Jovian_Scatter.0001.ec.fits' ; abit unnerving how much this removes the emissions
    print, 'UPON INSPECTION, best choice seems to be: ', best_scatter_spectrum
    
    ; -------------------------------------------Waterfall Plot Postscript Io and Fit Jovian scatter--------------------------------------------------------------------------------------------------------------

    ; Okay, now that we've established which Jupiter scatter spectrum best optimzes the residual, we can actually user it and do this scattered light subtraction

    ; get color versus ingress time
      timeColors = BytScl(findgen(11), Min=0, Max=11, Top=11)
      Color=timeColors[0]
      cgLoadCT, 33, NColors=8, /reverse
      pos        = cgLayout([1,3], OXMargin=[10,3], OYMargin=[8,5], YGap=0)
      Pos[1,0]   = Pos[1,0]*.8 & Pos[3,1] = Pos[3,1]*.8
      Pos[1,1]   = Pos[1,1]*.84 & Pos[3,2] = Pos[3,2]*.84
      
    cgPS_Open, filename = Reduced_Dir+'Na near-IR_scatterfit.eps', /ENCAPSULATED, xsize = 7.5, ysize = 5.
    !P.font = 1
    device, SET_FONT = 'Helvetica Bold', /TT_FONT
    !p.charsize = 1.5

    cgplot, WL, spec/1.e3, psym = 0, Ytitle = 'kR / '+cgsymbol('Angstrom'), $             ; Rayleighs to KiloRayleighs
      title = 'Io''s Sodium Response to Eclipse',  $
      xr = xr, /nodata, pos = pos[*,0], xtickformat = '(A1)', yr = yr

    for i = 0, n_elements(Eclipse_files)-1 do begin
      if i eq 0 then continue ; skip penumbra
      spec         = MRDFITS(reduced_dir+'R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, header, /fscale, /unsigned, /silent)
      spec_err     = MRDFITS(reduced_dir+'sig_R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, sig_header, /fscale, /unsigned, /silent)
      spec         = spec / 1.e3       ; to KR
      spec_err     = spec_err / 1.e3   ; to KR
      WL           = sxpar(header, 'CRVAL1')+findgen(N_elements(spec))*sxpar(header, 'Cdelt1')
      Na3_Io_frame = Na3 + Na3 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      Na4_Io_frame = Na4 + Na4 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      junk         = min( abs(WL - Na3_Io_frame), Na3_Io_Ind)
      junk         = min( abs(WL - Na4_Io_frame), Na4_Io_Ind)

      ; Which Jovian scatter are we using?
        ;Jup_Cal_Tell_Corr = MRDFITS(reduced_dir+best_scatter_spectrum, 0, Jup_Scatter_header, /fscale, /unsigned, /silent)
        Jup_Cal_Tell_Corr = MRDFITS(reduced_dir + Na1Na2_fit_params[i].background, 0, Jup_Scatter_header, /fscale, /unsigned, /silent)

      ; Define fitting indicies we need to avoid because of airglow itself.
        Ios_Airglow    = [Io_airglow + Io_airglow * sxpar(header, 'IO_DOPPL') / cspice_clight()]                                      ; wavelength locations of all airglow lines, Ignore telluric Na
        ARCES_1sigma_bandwidth = 1. * (mean(xr)/31500.) / 2.3548 ; 1. * FWHM in Angstroms / 2sqrt(2ln2) = 68% of flux enclosed within 1 sigma --> Appropriate for weak telluric airglow
        Io_AG_Ind = [] & Telluric_AG_Ind = []
        for j = 0, N_elements(Ios_Airglow)-1 do Io_AG_ind = [Io_AG_Ind, where(abs(Ios_Airglow[j] - WL) lt 1.2*ARCES_1sigma_bandwidth, /NULL)] ;  indicies within 1.5 sigma of an Io airglow line
        for j = 0, N_elements(Telluric_Airglow)-1 do Telluric_AG_ind = [Telluric_AG_Ind, where(abs(Telluric_Airglow[j] - WL) lt 1.2*ARCES_1sigma_bandwidth, /NULL)] ; indicies within 1 sigma of any Telluric airglow
        AG_ind = [Io_AG_Ind, Telluric_AG_Ind]
        AG_ind = AG_ind[UNIQ(AG_ind, SORT(AG_ind))]

      ; Fit the Jovian scattered light, use Amoeba/Correlation to match Jupiter in terms of smoothing and shifting, then multiply and add until matches the scattered background
        parinfo          = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0]}, 2)
        parinfo[1].fixed = 1                                                                        ; peg the additive "B" component of the MX_Plus_B at zero
        correl_indicies  = cgSetDifference(include_Wls, AG_Ind)
        correl_indicies  = cgSetDifference(correl_indicies, telluric_absorp_ind)
        weights          = 1./abs(Na3_Io_Ind - correl_indicies)^(local) + 1./abs(Na4_Io_Ind - correl_indicies)^(local)
        Mult_and_add     = mpfitfun('MX_Plus_B', Jup_Cal_Tell_Corr[correl_indicies], Spec[correl_indicies], Spec_err[correl_indicies], parinfo = parinfo,$
                                    [mean(Spec[correl_indicies]/Jup_Cal_Tell_Corr[correl_indicies]), 0.], /NaN, weights = weights, status = Scatter_fit_status, PERROR = err_a, /quiet)
        jup              = Mult_and_add[0]*Jup_Cal_Tell_Corr + Mult_and_add[1]
  
          ; Iterate Amoeba until it gives an answer
            Ftol = 1.e-4 &  Xtol = 1.e-4 & count = 0        ; Define initial tolerances for Amoeba 
            smooth_and_shift = [2., 0.] & scale  = [2., 2.] ; Define intial guess and simplex for Amoeba 
            REPEAT BEGIN
              trial_smooth_and_shift = AMOEBAX(Ftol, xtol, function_name='shift_smooth', SCALE = Scale, P0 = smooth_and_shift, FUNCTION_VALUE = fval)
              if N_elements (N_elements(trial_smooth_and_shift) lt 2) then begin              ; If Amoeba crashes
                xtol = xtol * 1.5 & ftol = ftol * 1.5                                         ; Grow the tolerances
                smooth_and_shift = smooth_and_shift + [RANDOMN(seed), RANDOMN(seed)]          ; Change the guess
                SCALE            = SCALE + [RANDOMN(seed), RANDOMN(seed)]                     ; Change the simplex
                count = count+1
              endif  
            ENDREP UNTIL ( (N_elements(trial_smooth_and_shift) eq 2) or (count gt 500) )
            smooth_and_shift = AMOEBAX(Ftol, xtol, function_name='shift_smooth', SCALE = [2., 2.], P0 = smooth_and_shift, FUNCTION_VALUE = fval)
            if N_elements(smooth_and_shift) eq 1 then stop ; check for bugs
          
        smooth_and_shift = [Na1Na2_fit_params[i].smooth, Na1Na2_fit_params[i].shift]

        smoothed_shifted = interpolate(gauss_smooth(jup, smooth_and_shift[0]), findgen(n_elements(Jup)) - smooth_and_shift[1])
        Mult_and_add     = mpfitfun('MX_Plus_B', smoothed_shifted[correl_indicies], Spec[correl_indicies], Spec_err[correl_indicies], [mean(Spec[correl_indicies]/smoothed_shifted[correl_indicies]), 0.], $
                                    /NaN, status = Scatter_fit_status, PERROR = err_a, /quiet, weights = weights, parinfo = parinfo)
        scatter_fit      = Mult_and_add[0]*smoothed_shifted + Mult_and_add[1]
        residual         = Spec - scatter_fit
        residual_err     = make_array( N_elements(residual), value = stddev(residual[correl_indicies]) )

      ; plot the spectrum, fit and indices used for the fit (if we're in debug mode)
        ;cgplot, WL[correl_indicies], scatter_fit[correl_indicies], psym=14, /overplot
        cgplot, WL, scatter_fit, color = timeColors[i], linestyle = 1, thick = 5, /overplot
        cgplot, WL, spec, color = timeColors[i], thick = 5, /overplot

      LSF_fitting_ind1         = where( abs(wl - Na3_Io_frame) lt 0.3, /NULL)
      LSF_fitting_ind2         = where( abs(wl - Na4_Io_frame) lt 0.3, /NULL)
      LSF_Fit1                 = mpfitpeak(WL[LSF_fitting_ind1], residual[LSF_fitting_ind1], a, STATUS = Line_fit_Status1, /POSITIVE, nterms = 3)
      LSF_Fit2                 = mpfitpeak(WL[LSF_fitting_ind2], residual[LSF_fitting_ind2], b, STATUS = Line_fit_Status2, /POSITIVE, nterms = 3)
      WL_array[*, i]           = WL[include_WLs]
      plot_residual_array[*, i]= residual[plot_WLs]
      LSF_Fit_array[*, i]      = gaussian(WL[include_WLs], a) + gaussian(WL[include_WLs], b)                ; LSF_Fit
      Brightness_array[i]      = A[0]*A[2]*SQRT(2*!DPI) + B[0]*B[2]*SQRT(2*!DPI)                            ; Area under the Gaussian
      err_Brightness_array[i]  = stddev(residual[correl_indicies])*A[2]*SQRT(2*!DPI) + stddev(residual[correl_indicies])*B[2]*SQRT(2*!DPI)  ; use the std dev of the residual and the width of the line
      T_P_Shadow_array[i]      = sxpar(header, 'T_PSHADO')
      T_U_Shadow_array[i]      = sxpar(header, 'T_USHADO')
      EXPTime_array[i]         = sxpar(sig_header, 'EXPTIME') / 60.                                         ; for some reason, this keyword is bad in the regular headers, so use sig_header as a workaround
      DopplerShift_array1[i]   = Na3_Io_frame
      DopplerShift_array2[i]   = Na4_Io_frame
    endfor

      ; Residual plot
        N_frames            = total(Brightness_array gt 0.)
        cgplot, spec, WL, psym = 0, xtickformat = '(A1)', yminor = 2, $
          yr = [-2.99, 5.5 ], xr = XR, /nodata, pos = pos[*,1], /noerase
        cgtext, xr[0] - 2., -2., 'Residual (kR / '+cgsymbol('Angstrom')+')', orientation = 90, alignment = 0.5

        ; Co-align to Io's Doppler Shift ???
          shift_By = mean( [transpose( DopplerShift_array1 - Mean(DopplerShift_array1[1:*]) ), $
                            transpose( DopplerShift_array2 - Mean(DopplerShift_array2[1:*]) )], dimension = 1 )
          aligned_residuals = plot_residual_array

        for i = 1, n_elements(Eclipse_files)-1 do begin
          cgplot, WL[plot_WLs], plot_residual_array[*, i], /OVERPLOT, COLOR = timeColors[i], psym=10
          cgplot, [DopplerShift_array1[i], DopplerShift_array1[i]], [5.5, 3.9], COLOR = timeColors[i], /overplot
          cgplot, [DopplerShift_array2[i], DopplerShift_array2[i]], [5.5, 4.7], COLOR = timeColors[i], /overplot
          aligned_residuals[*,i] = interpol( plot_residual_array[*,i], WL[plot_WLs], WL[plot_WLs] - Shift_By[i])
        endfor
        ;cgplot, WL[correl_indicies], 5.5*weights/max(weights), psym = 1, /overplot
        co_add = total(aligned_residuals[*,1:*], 2) / n_frames

        ; plot the aligned & co-added spectra         
          cgplot, WL[plot_WLs], co_add, pos = pos[*,2], thick=5, /noerase, xr = xr, yr = [-1.9, 3.2], /yst,  $
            xtitle = 'Wavelength ('+cgsymbol('Angstrom')+')', yminor = 2
          cgplot, xr, [0.,0.], linestyle = 2, /overplot
          cgplot, [Mean(DopplerShift_array1[1:*]), Mean(DopplerShift_array1[1:*])], [-1.9, 3.2], linestyle = 1, /overplot
          cgplot, [Mean(DopplerShift_array2[1:*]), Mean(DopplerShift_array2[1:*])], [-1.9, 3.2], linestyle = 1, /overplot      
  cgps_close
endif
if Part eq 1.3 then begin ; Make waterfall plots, K-D does not need telluric corrected spectra.

    ;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ;--------------------------------------------------------Get the K Time Series and Scatter Fit-------------------------------------------------------------------------------
    ;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ; METHOD:
    ;   Loop though all possible spectra that could be fit to jovian scatter and choose what works best. Avoid Io and telluric airglow in the fitting.
    ;   Three steps to the fit: 1) Use Mpfit to roughly add & multiply the Jovian scatter spectra to approach the background in the Io exposures
    ;                           2) Use Amoeba to shift and smooth the spectra, maximizing correlation
    ;                           3) Use Mpfit again to fine tune the add & multiply of Jovian scatter spectra needed to match the background in the Io exposures
    ; NOTES:
    ;   This section can be a little fiddly. Amoeba is prone to infinite loops if two points in the simplex give the same correlation.
    ;   Therefore, include_WLs needs to be careful tuned to keep Amoeba stable.
    ;   It's critical to inspect with plot, WL[correl_indicies], scatter_fit[correl_indicies], uncomment this line to inspect and adjust things
    ;   With this inspection, it is useful to adjust *include_WLs*, and *AG_ind* to cover regions of continuum as close as possible to Io's airglow
    ;   Parameter information is not included in the Gaussian fitting to the line spread function in the residual. I couldn't figgure out why it failed with status = 0

    ; Fit various spectra to the Jovian Scatter, *use* whichever one produces the lowest residual. Note this method will never work if there's much telluric absorption.
      Fit_Me_To_The_Scatter = [ 'R_per_A_Jupiter_Center.ec.fits', $
                                'R_per_A_Jovian_Scatter.0001.ec.fits', $
                                'R_per_A_Jovian_Scatter.0002.ec.fits', $
                                'R_per_A_Jovian_Scatter.0003.ec.fits' ]

    ; setup the plot axis
      spec = MRDFITS(reduced_dir+'R_per_A_'+Eclipse_files[1]+'.ec.fits', 0, header, /fscale, /unsigned, /silent)
      WL   = sxpar(header, 'CRVAL1')+findgen(N_elements(spec))*sxpar(header, 'Cdelt1')

      Case date of
        'UT180320': begin
          YR       = [15,78]
          YR_resid = [-4, 4]  
          XR       = [7661,7703]
        end
        'UT190812': yr   = [1.1e1, 4e1]
      endcase

    ; Define arrays
      include_WLs            = where( (abs(wl - K1+.4) lt 1.1) or (abs(wl - K2+.4) lt 1.1), /NULL ); wavelength regions over which we'll fit the jovian scatter
      plot_WLs               = where( (wl gt xr[0]) and (wl lt xr[1]), /NULL )                       ; wavelengths to actually plot
      residual_array         = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
      plot_residual_array    = fltarr(N_elements(plot_WLs), n_elements(Eclipse_files))
      LSF_Fit_array          = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
      WL_array               = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
      Brightness_array       = fltarr(n_elements(Eclipse_files))
      err_Brightness_array   = fltarr(n_elements(Eclipse_files))
      T_P_Shadow_array       = fltarr(n_elements(Eclipse_files))
      T_U_Shadow_array       = fltarr(n_elements(Eclipse_files))
      EXPTime_array          = fltarr(n_elements(Eclipse_files))
      DopplerShift_array1    = fltarr(n_elements(Eclipse_files))                          ; For the D2 line
      DopplerShift_array2    = fltarr(n_elements(Eclipse_files))                          ; For the D1 line
      Gof                    = 0.                                                         ; Goodness of fit (co-addded mean residual)
      GOF_array              = fltarr(N_elements(Fit_Me_To_The_Scatter))
      Use_Scatter            = strarr(N_elements(Eclipse_files))
      K1K2_Fit_params        = REPLICATE({params, brightness:0.0, err_brightness:0.0, mult:0.0, smooth:0.0, shift:0.0, background:'', T_P_shadow:0.0, exptime:0.0}, n_elements(Eclipse_files))

    for i = 0, n_elements(Eclipse_files)-1 do begin
      if i eq 0 then continue ; skip penumbra
      spec        = MRDFITS(reduced_dir+'R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, header, /fscale, /unsigned, /silent)
      spec_err    = MRDFITS(reduced_dir+'sig_R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, sig_header, /fscale, /unsigned, /silent)
      WL          = sxpar(header, 'CRVAL1')+findgen(N_elements(spec))*sxpar(header, 'Cdelt1')
      K1_Io_frame = K1 + K1 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      K2_Io_frame = K2 + K2 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      junk         = min( abs(WL - K1_Io_frame), K1_Io_Ind)
      junk         = min( abs(WL - K2_Io_frame), K2_Io_Ind)

      for h = 0, n_elements(Fit_Me_To_The_Scatter)-1 do begin
        ; Which Jovian scatter are we using?
          Jup_Cal_Tell_Corr = MRDFITS(reduced_dir+Fit_Me_To_The_Scatter[h], 0, header, /fscale, /unsigned, /silent)

        ; Define fitting indicies we need to avoid because of airglow itself.
          Ios_Airglow    = [Io_airglow + Io_airglow * sxpar(header, 'IO_DOPPL') / cspice_clight()]                                      ; wavelength locations of all airglow lines, Ignore telluric K
          ARCES_1sigma_bandwidth = 1. * (mean(xr)/31500.) / 2.3548 ; 1. * FWHM in Angstroms / 2sqrt(2ln2) = 68% of flux enclosed within 1 sigma --> Appropriate for weak telluric airglow
          Io_AG_Ind = [] & Telluric_AG_Ind = []
          for j = 0, N_elements(Ios_Airglow)-1 do Io_AG_ind = [Io_AG_Ind, where(abs(Ios_Airglow[j] - WL) lt 2.*ARCES_1sigma_bandwidth, /NULL)] ;  indicies within 2 sigma of an Io airglow line
          for j = 0, N_elements(Telluric_Airglow)-1 do Telluric_AG_ind = [Telluric_AG_Ind, where(abs(Telluric_Airglow[j] - WL) lt ARCES_1sigma_bandwidth, /NULL)] ; indicies within 1 sigma of any Telluric airglow
          AG_ind = [Io_AG_Ind, Telluric_AG_Ind]
          AG_ind = AG_ind[UNIQ(AG_ind, SORT(AG_ind))]

        ; Fit the Jovian scattered light, use Amoeba/Correlation to match Jupiter in terms of smoothing and shifting, then multiply and add until matches the scattered background
          correl_indicies  = cgSetDifference(include_Wls, AG_Ind)
          Mult_and_add     = mpfitfun('MX_Plus_B', Jup_Cal_Tell_Corr[correl_indicies], Spec[correl_indicies], Spec_err[correl_indicies], parinfo = MX_plus_B_parinfo, $
                                      [mean(Spec[correl_indicies]/Jup_Cal_Tell_Corr[correl_indicies]), 0.], /KN, status = Scatter_fit_status, PERROR = err_a, /quiet)
          jup              = Mult_and_add[0]*Jup_Cal_Tell_Corr + Mult_and_add[1]
  
          ; Iterate Amoeba until it gives an answer
            Ftol = 1.e-4 &  Xtol = 1.e-4 & count = 0        ; Define initial tolerances for Amoeba 
            smooth_and_shift = [2., 0.] & scale  = [2., 2.] ; Define nitial guess and simplex for Amoeba 
            REPEAT BEGIN
              trial_smooth_and_shift = AMOEBAX(Ftol, xtol, Function_Name='shift_smooth', SCALE = Scale, P0 = smooth_and_shift, FUNCTION_VALUE = fval)
              if N_elements (N_elements(trial_smooth_and_shift) lt 2) then begin              ; If Amoeba crashes
                xtol = xtol * 1.5 & ftol = ftol * 1.5                                         ; Grow the tolerances
                smooth_and_shift = smooth_and_shift + [RANDOMN(seed), RANDOMN(seed)]          ; Change the guess
                SCALE            = SCALE + [RANDOMN(seed), RANDOMN(seed)]                     ; Change the simplex
                count = count+1
              endif  
            ENDREP UNTIL ( (N_elements(trial_smooth_and_shift) eq 2) or (count gt 500) )
            smooth_and_shift = AMOEBAX(Ftol, xtol, Function_Name='shift_smooth', SCALE = [2., 2.], P0 = smooth_and_shift, FUNCTION_VALUE = fval)
            if N_elements(smooth_and_shift) eq 1 then stop ; check for bugs

          smoothed_shifted = interpolate(gauss_smooth(jup, smooth_and_shift[0]), findgen(n_elements(Jup)) - smooth_and_shift[1])
          Mult_and_add     = mpfitfun('MX_Plus_B', smoothed_shifted[correl_indicies], Spec[correl_indicies], Spec_err[correl_indicies], parinfo = MX_plus_B_parinfo, $
                                      [mean(Spec[correl_indicies]/smoothed_shifted[correl_indicies]), 0.], /KN, status = Scatter_fit_status, PERROR = err_a, /quiet)
                                      
          scatter_fit      = Mult_and_add[0]*smoothed_shifted + Mult_and_add[1]
          residual         = Spec - scatter_fit
          GOF_array[h]     = GOF + mean(abs(residual[correl_indicies]))
        ;print, 'K D Jupiter Scatter Goodness of Fit = ', GOF_array[h]
      endfor
      best_GoF = min(GOF_array, best_fit)
      Use_Scatter[i] = Fit_Me_To_The_Scatter[best_fit]
      print, 'Frame: ', Eclipse_files[i], ' should take the Jupiter scatter from: ', Use_Scatter[i], ', Goodness of fit = ', best_GoF
    endfor
    ; -------------------------------------------Waterfall Plot Postscript Io and Fit Jovian scatter--------------------------------------------------------------------------------------------------------------

    ; Okay, now that we've established which Jupiter scatter spectrum best optimzes the residual, we can actually user it and do this scattered light subtraction

    ; get color versus ingress time
      timeColors = BytScl(findgen(11), Min=0, Max=11, Top=11)
      Color=timeColors[0]
      cgLoadCT, 33, NColors=8, /reverse
      pos        = cgLayout([1,3], OXMargin=[10,3], OYMargin=[8,5], YGap=0)
      Pos[1,0]   = Pos[1,0]*.8 & Pos[3,1] = Pos[3,1]*.8
      Pos[1,1]   = Pos[1,1]*.84 & Pos[3,2] = Pos[3,2]*.84
      
    cgPS_Open, Filename = Reduced_Dir+'K-D_scatterfit.eps', /ENCAPSULATED, xsize = 7.5, ysize = 5.
      !P.font = 1
      device, SET_FONT = 'Helvetica Bold', /TT_FONT
      !p.charsize = 1.5

    cgplot, WL, spec/1.e3, psym = 0, Ytitle = 'kR / '+cgsymbol('Angstrom'), $             ; Rayleighs to KiloRayleighs
      title = 'Io''s Potassium Response to Eclipse',  $
      xr = xr, /nodata, pos = pos[*,0], xtickformat = '(A1)', yr = yr

    for i = 0, n_elements(Eclipse_files)-1 do begin
      if i eq 0 then continue ; skip penumbra
      spec         = MRDFITS(reduced_dir+'R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, header, /fscale, /unsigned, /silent)
      spec_err     = MRDFITS(reduced_dir+'sig_R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, sig_header, /fscale, /unsigned, /silent)
      spec         = spec / 1.e3       ; to KR
      spec_err     = spec_err / 1.e3   ; to KR
      WL           = sxpar(header, 'CRVAL1')+findgen(N_elements(spec))*sxpar(header, 'Cdelt1')
      K1_Io_frame = K1 + K1 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      K2_Io_frame = K2 + K2 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      junk         = min( abs(WL - K1_Io_frame), K1_Io_Ind)
      junk         = min( abs(WL - K2_Io_frame), K2_Io_Ind)
      
      ; Which Jovian scatter are we using?
        Jup_Cal_Tell_Corr = MRDFITS(reduced_dir+Use_Scatter[i], 0, Jup_Scatter_header, /fscale, /unsigned, /silent)

      ; Define fitting indicies we need to avoid because of airglow itself.
        Ios_Airglow    = [Io_airglow + Io_airglow * sxpar(header, 'IO_DOPPL') / cspice_clight()]                                      ; wavelength locations of all airglow lines, Ignore telluric K
        ARCES_1sigma_bandwidth = 1. * (mean(xr)/31500.) / 2.3548 ; 1. * FWHM in Angstroms / 2sqrt(2ln2) = 68% of flux enclosed within 1 sigma --> Appropriate for weak telluric airglow
        Io_AG_Ind = [] & Telluric_AG_Ind = []
        for j = 0, N_elements(Ios_Airglow)-1 do Io_AG_ind = [Io_AG_Ind, where(abs(Ios_Airglow[j] - WL) lt 2.*ARCES_1sigma_bandwidth, /NULL)] ;  indicies within 2 sigma of an Io airglow line
        for j = 0, N_elements(Telluric_Airglow)-1 do Telluric_AG_ind = [Telluric_AG_Ind, where(abs(Telluric_Airglow[j] - WL) lt ARCES_1sigma_bandwidth, /NULL)] ; indicies within 1 sigma of any Telluric airglow
        AG_ind = [Io_AG_Ind, Telluric_AG_Ind]
        AG_ind = AG_ind[UNIQ(AG_ind, SORT(AG_ind))]

      ; Fit the Jovian scattered light, use Amoeba/Correlation to match Jupiter in terms of smoothing and shifting, then multiply and add until matches the scattered background
        correl_indicies  = cgSetDifference(include_Wls, AG_Ind)
        weights          = 1./abs(K1_Io_Ind - correl_indicies)^(1.5) + 1./abs(K2_Io_Ind - correl_indicies)^(1.5)
        Mult_and_add     = mpfitfun('MX_Plus_B', Jup_Cal_Tell_Corr[correl_indicies], Spec[correl_indicies], Spec_err[correl_indicies], parinfo = MX_plus_B_parinfo, $
                            [mean(Spec[correl_indicies]/Jup_Cal_Tell_Corr[correl_indicies]), 0.], /KN, status = Scatter_fit_status, PERROR = err_a, /quiet)
        jup              = Mult_and_add[0]*Jup_Cal_Tell_Corr + Mult_and_add[1]

        ; Iterate Amoeba until it gives an answer
          Ftol = 1.e-4 &  Xtol = 1.e-4 & count = 0        ; Define initial tolerances for Amoeba
          smooth_and_shift = [2., 0.] & scale  = [2., 2.] ; Define nitial guess and simplex for Amoeba
          REPEAT BEGIN
            trial_smooth_and_shift = AMOEBAX(Ftol, xtol, Function_Name='shift_smooth', SCALE = Scale, P0 = smooth_and_shift, FUNCTION_VALUE = fval)
            if N_elements (N_elements(trial_smooth_and_shift) lt 2) then begin              ; If Amoeba crashes
              xtol = xtol * 1.5 & ftol = ftol * 1.5                                         ; Grow the tolerances
              smooth_and_shift = smooth_and_shift + [RANDOMN(seed), RANDOMN(seed)]          ; Change the guess
              SCALE            = SCALE + [RANDOMN(seed), RANDOMN(seed)]                     ; Change the simplex
              count = count+1
            endif
          ENDREP UNTIL ( (N_elements(trial_smooth_and_shift) eq 2) or (count gt 500) )
          smooth_and_shift = AMOEBAX(Ftol, xtol, Function_Name='shift_smooth', SCALE = [2., 2.], P0 = trial_smooth_and_shift, FUNCTION_VALUE = fval)
          if N_elements(smooth_and_shift) eq 1 then stop ; check for bugs
        
        smoothed_shifted = interpolate(gauss_smooth(jup, smooth_and_shift[0]), findgen(n_elements(Jup)) - smooth_and_shift[1])
        Mult_and_add     = mpfitfun('MX_Plus_B', smoothed_shifted[correl_indicies], Spec[correl_indicies], Spec_err[correl_indicies], [mean(Spec[correl_indicies]/smoothed_shifted[correl_indicies]), 0.], $
                                    /KN, status = Scatter_fit_status, PERROR = err_a, /quiet, weights = weights, parinfo = MX_plus_B_parinfo)
        scatter_fit      = Mult_and_add[0]*smoothed_shifted + Mult_and_add[1]
        residual         = Spec - scatter_fit
        residual_err     = make_array( N_elements(residual), value = stddev(residual[correl_indicies]) )

      ; plot the spectrum, fit and indices used for the fit (if we're in debug mode)
        ;cgplot, WL[correl_indicies], scatter_fit[correl_indicies] + 60.-(i*12.), psym=14, /overplot
        cgplot, WL, scatter_fit + 60.-(i*12.), color = timeColors[i], linestyle = 1, thick = 5, /overplot
        cgplot, WL, spec + 60. -(i*12.), color = timeColors[i], thick = 5, /overplot

      LSF_fitting_ind1         = where( abs(wl - K1_Io_frame) lt 0.3, /NULL) 
      LSF_fitting_ind2         = where( abs(wl - K2_Io_frame) lt 0.3, /NULL) 
      LSF_Fit1                 = mpfitpeak(WL[LSF_fitting_ind1], residual[LSF_fitting_ind1], a, STATUS = Line_fit_Status1, /POSITIVE, nterms = 3)
      LSF_Fit2                 = mpfitpeak(WL[LSF_fitting_ind2], residual[LSF_fitting_ind2], b, STATUS = Line_fit_Status2, /POSITIVE, nterms = 3)
      WL_array[*, i]           = WL[include_WLs]
      plot_residual_array[*, i]= residual[plot_WLs]
      LSF_Fit_array[*, i]      = gaussian(WL[include_WLs], a) + gaussian(WL[include_WLs], b)                ; LSF_Fit
      K1K2_fit_params[i]       = {params, A[0]*A[2]*SQRT(2*!DPI) + B[0]*B[2]*SQRT(2*!DPI), stddev(residual[correl_indicies])*A[2]*SQRT(2*!DPI) + stddev(residual[correl_indicies])*B[2]*SQRT(2*!DPI), $
                                          Mult_and_add[0], smooth_and_shift[0], smooth_and_shift[1], Use_Scatter[i], float(sxpar(header, 'T_PSHADO')), float(sxpar(sig_header, 'EXPTIME')) / 60.}      
      Brightness_array[i]      = A[0]*A[2]*SQRT(2*!DPI) + B[0]*B[2]*SQRT(2*!DPI)                            ; Area under the Gaussian
      err_Brightness_array[i]  = stddev(residual[correl_indicies])*A[2]*SQRT(2*!DPI) + stddev(residual[correl_indicies])*B[2]*SQRT(2*!DPI)  ; use the std dev of the residual and the width of the line
      T_P_Shadow_array[i]      = sxpar(header, 'T_PSHADO')
      T_U_Shadow_array[i]      = sxpar(header, 'T_USHADO')
      EXPTime_array[i]         = sxpar(sig_header, 'EXPTIME') / 60.                                         ; for some reason, this keyword is bad in the regular headers, so use sig_header as a workaround
      DopplerShift_array1[i]   = K1_Io_frame
      DopplerShift_array2[i]   = K2_Io_frame
    endfor

    ; Residual plot
      N_frames            = total(Brightness_array gt 0.)
      cgplot, spec, WL, psym = 0, xtickformat = '(A1)', yminor = 2, yr = YR_resid, xr = XR, /nodata, pos = pos[*,1], /noerase
      cgtext, xr[0] - 2., -2., 'Residual (kR / '+cgsymbol('Angstrom')+')', orientation = 90, alignment = 0.5

      ; Co-align to Io's Doppler Shift ???
      shift_By = mean( [transpose( DopplerShift_array1 - Mean(DopplerShift_array1[1:*]) ), $
        transpose( DopplerShift_array2 - Mean(DopplerShift_array2[1:*]) )], dimension = 1 )
      aligned_residuals = plot_residual_array

      for i = 1, n_elements(Eclipse_files)-1 do begin
        cgplot, WL[plot_WLs], plot_residual_array[*, i], /OVERPLOT, COLOR = timeColors[i], psym=10
        cgplot, [DopplerShift_array1[i], DopplerShift_array1[i]], [YR_resid[1], 4], COLOR = timeColors[i], /overplot
        cgplot, [DopplerShift_array2[i], DopplerShift_array2[i]], [YR_resid[1], 4], COLOR = timeColors[i], /overplot
        aligned_residuals[*,i] = interpol( plot_residual_array[*,i], WL[plot_WLs], WL[plot_WLs] - Shift_By[i])
      endfor
      co_add = total(aligned_residuals[*,1:*], 2) / n_frames

      ; plot the aligned & co-added spectra
      cgplot, WL[plot_WLs], co_add, pos = pos[*,2], thick=5, /noerase, xr = xr, yr = YR_resid/2., /yst,  $
              xtitle = 'Wavelength ('+cgsymbol('Angstrom')+')', yminor = 2
      cgplot, xr, [0.,0.], linestyle = 2, /overplot
      cgplot, [Mean(DopplerShift_array1[1:*]), Mean(DopplerShift_array1[1:*])], YR_resid, linestyle = 1, /overplot
      cgplot, [Mean(DopplerShift_array2[1:*]), Mean(DopplerShift_array2[1:*])], YR_resid, linestyle = 1, /overplot  
    cgps_close
endif

;==========================================  Combined lightcurve  ======================================================================

;print, 'Ajello et al. 2008 ratio [7774:8446:9225]', [3.90e-19, 2.15e-19, 6.63e-19] / 2.15e-19  
if part eq 2 then begin
  restore, Reduced_Dir+'Na1Na2_fit_params.sav' 
  restore, Reduced_Dir+'O1_fit_params.sav' 
  restore, Reduced_Dir+'O2_fit_params.sav' 
  restore, Reduced_Dir+'O3_fit_params.sav'

  pos =  [.12,.17,.95,.9]
  cgPS_Open, filename = Reduced_Dir+'Combined_lightcurve.eps', /ENCAPSULATED, xsize = 6., ysize = 4.
    !P.font=1
    device, SET_FONT = 'Helvetica Bold', /TT_FONT

    yr = [0, 5.]
    cgplot, O2_fit_params.T_P_Shadow, O2_fit_params.Brightness, psym = 15, Ytitle = 'Disk-Averaged Brightness (KR)', xtitle = 'Minutes After Ingress', $
            title = 'Io''s Response to Ingress', yr = yr, xr = [range1,range2], /nodata, pos = pos,YTICKFORMAT="(A1)",yticks = 1
            
    cgLoadCT, 33, NColors=8, /reverse
    if not ingress then x = findgen(Umbra_ET - PenUmbra_ET) / 60.
    if ingress then x = findgen(Umbra_ET - PenUmbra_ET) / 60.
    colors = reverse(BytScl(x, MIN=min(x), MAX=2.*max(x)))+128
    FOR j=0,n_elements(x)-2 DO BEGIN
      xpoly = [x[j],         x[j], x[j+1],       x[j+1],         x[j]]
      ypoly = [  0., !Y.CRange[1], !Y.CRange[1], 0., 0.]
      cgColorFill, xpoly, ypoly, Color=colors[j]
    ENDFOR
    if ingress then begin
      xpoly = [max(x),     max(x), !X.CRange[1],  !X.CRange[1],  max(x)]
      ypoly = [  0., !Y.CRange[1], !Y.CRange[1], 0., 0.]
      cgColorFill, xpoly, ypoly, color = 'charcoal', /LINE_FILL, orientation = 45., thick = .1
      cgColorFill, xpoly, ypoly, color = 'charcoal', /LINE_FILL, orientation = -45., thick = .1
    endif

    cgaxis, yaxis = 0, yr = yr, /ystyle                        ; repair the axis damage that the Penumbra did
    cgaxis, xaxis = 0, xr = time_range                         ; repair axis damage
    cgaxis, xaxis = 1, xr = time_range, xtickformat = '(A1)'   ; repair axis damage
    
    
    Na1Na2_fit_params[0].Brightness = !Values.F_Nan 
    O1_fit_params[0].Brightness = !Values.F_Nan 
    O2_fit_params[0].Brightness = !Values.F_Nan 
    O3_fit_params[0].Brightness = !Values.F_Nan 
    
    READCOL,'D:\DATA\March20_sun_ec.dat',col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12, STRINGSKIP = '#', /Silent
    col3in = col3[0]
    col5in = col5[0]
    col7in = col7[0]
    col9in = col9[0]
    col11in = col11[0]
    
    col3[0] = !Values.F_Nan
    col3[1] = !Values.F_Nan
    col5[0] = !Values.F_Nan
    col5[1] = !Values.F_Nan
    col7[0] = !Values.F_Nan
    col7[1] = !Values.F_Nan
    col9[0] = !Values.F_Nan
    col9[1] = !Values.F_Nan
    col11[0] = !Values.F_Nan
    col11[1] = !Values.F_Nan

    cgplot, O2_fit_params.T_P_Shadow, O2_fit_params.Brightness, psym = 16, /overplot, color = 'red', $
      ERR_YLOW = O2_fit_params.ERR_Brightness, ERR_YHigh = O2_fit_params.ERR_Brightness, ERR_XLOW = O2_fit_params.exptime/2., ERR_XHigh = O2_fit_params.exptime/2.
    cgplot, O2_fit_params.T_P_Shadow, O2_fit_params.Brightness, color = 'red', /overplot
    
    cgplot, O3_fit_params.T_P_Shadow, O3_fit_params.Brightness, psym = 15, /overplot, color = 'red', $
      ERR_YLOW = O3_fit_params.ERR_Brightness, ERR_YHigh = O3_fit_params.ERR_Brightness, ERR_XLOW = O3_fit_params.exptime/2., ERR_XHigh = O3_fit_params.exptime/2., ERR_CLIP = 1
    cgplot, O3_fit_params.T_P_Shadow, O3_fit_params.Brightness, color = 'red', /overplot
    
    ;cgplot, O1_fit_params.T_P_Shadow, O1_fit_params.Brightness, psym = 16, /overplot, color = 'green', /err_clip, $
      ;ERR_YLOW = O1_fit_params.ERR_Brightness, ERR_YHigh = O1_fit_params.ERR_Brightness, ERR_XLOW = O1_fit_params.exptime/2., ERR_XHigh = O1_fit_params.exptime/2.
    ;cgplot, O1_fit_params.T_P_Shadow, O1_fit_params.Brightness, color = 'green', /overplot
    
    cgplot, Na1Na2_fit_params.T_P_Shadow, Na1Na2_fit_params.Brightness, psym = 16, /overplot, color = 'orange', $
      ERR_YLOW = Na1Na2_fit_params.ERR_Brightness, ERR_YHigh = Na1Na2_fit_params.ERR_Brightness, ERR_XLOW = Na1Na2_fit_params.exptime/2., ERR_XHigh = Na1Na2_fit_params.exptime/2.
    cgplot, Na1Na2_fit_params.T_P_Shadow, Na1Na2_fit_params.Brightness, color = 'orange', /overplot
    
    cgAxis, YAxis=1, YRange=[0, 6.5],title= 'Flux Density of Sulfur- Species [Jy]', COLOR = 'blue', /Save
    
    cgplot, col2, (col3 + col5)/1000., psym = 17, /overplot, color = 'blue', $
      ERR_YLOW = (col4 + col6)/1000., ERR_YHigh = (col4 + col6)/1000., ERR_XLOW = [0.,0.,1.7,1.84,5.68], ERR_XHigh = [0.,0.,1.3,1.81,4.9]
    cgplot, col2, (col3 + col5)/1000., color = 'blue', /overplot
    cgplot, [0, col2[2]], [(col3in + col5in)/1000., (col3[2] + col5[2])/1000.], color = 'blue', LineStyle=2, /overplot
    
    ;cgplot, col2, col5/1000., psym = 18, /overplot, color = 'blue', $
    ;  ERR_YLOW = col4/1000., ERR_YHigh = col4/1000., ERR_XLOW = [0.,0.,1.7,1.84,5.68], ERR_XHigh = [0.,0.,1.3,1.81,4.9]
    ;cgplot, col2, col5/1000., color = 'blue', /overplot
    
    ;cgplot, col2, col7/1000., psym = 19, /overplot, color = 'blue', $
    ;  ERR_YLOW = col4/1000., ERR_YHigh = col4/1000., ERR_XLOW = [0.,0.,1.7,1.84,5.68], ERR_XHigh = [0.,0.,1.3,1.81,4.9]
    ;cgplot, col2, col7/1000., color = 'blue', /overplot
    
    ;cgplot, col2, col9/1000., psym = 20, /overplot, color = 'blue', $
    ;  ERR_YLOW = col4/1000., ERR_YHigh = col4/1000., ERR_XLOW = [0.,0.,1.7,1.84,5.68], ERR_XHigh = [0.,0.,1.3,1.81,4.9]
    ;cgplot, col2, col9/1000., color = 'blue', /overplot
    
    cgplot, col2, (col7 + col9)/1000., psym = 19, /overplot, color = 'blue', $
      ERR_YLOW = (col8 + col10)/1000., ERR_YHigh = (col8 + col10)/1000., ERR_XLOW = [0.,0.,1.7,1.84,5.68], ERR_XHigh = [0.,0.,1.3,1.81,4.9]
    cgplot, col2, (col7 + col9)/1000., color = 'blue', /overplot
    cgplot, [0, col2[2]], [(col7in + col9in)/1000., (col7[2] + col9[2])/1000.], color = 'blue', LineStyle=2, /overplot
    
    cgplot, col2, col11/1000., psym = 16, /overplot, color = 'blue', $
      ERR_YLOW = col4/1000., ERR_YHigh = col4/1000., ERR_XLOW = [0.,0.,1.7,1.84,5.68], ERR_XHigh = [0.,0.,1.3,1.81,4.9]
    cgplot, col2, col11/1000., color = 'blue', /overplot
    
    cgtext, 2.5, !Y.CRange[1]/4., 'Penumbral Eclipse', orientation = 90., color = 'white'

    Case date of
      'UT180320': cglegend, title = ['[O I] 6300'+cgsymbol('Angstrom'), '[O I] 6364'+cgsymbol('Angstrom'), 'Na D1 + D2', 'SO2 spw1 + spw2', 'SO2 spw5 + spw6', 'SO spw1 + spw3'], Psym = [16, 15, 16, 17,19, 16], charsize = 0.9, $
                            bg_color = 'white', color = ['red', 'red', 'orange', 'blue', 'blue', 'blue'], Location=[0.67, 0.87];, /Background, /BOX
      'UT190812': AL_legend, ['K-D in Umbra 190812'], Psym = 4, /right, charsize = 1, /clear
    endcase
    cgps_close
endif
stop

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
        A0V1_airmass     = 1./cos( (90. - ten(sxpar(header, 'TELAL'))) / !radeg )
      
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
          A0V2_airmass     = 1./cos( (90. - ten(sxpar(header, 'TELAL'))) / !radeg )
      
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
      
        ; Comapare publications of Jupiter's spectral albedo at disk center
      
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
            Jup_center_airmass     = 1./cos( (90. - ten(sxpar(Jup_center_header, 'TELAL'))) / !radeg )
            
          ; Make a master telluric correction that combines the different standard stars 
            master_telluric_D         = replicate(1., N_elements(jup_center))
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
            ;cgplot, WL, jup_tell_corr, xr = [5888, 5898]
            ;cgplot, WL, shift(GAUSS_SMOOTH(expected_flux, 4.9)*3.75e-9, 6), /overplot, color = 'blue'
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
      Jup_center_airmass      = 1./cos( (90. - ten(sxpar(Jup_center_header, 'TELAL'))) / !radeg )

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
      ;cgplot, WL, jup_tell_corr, xr = [5888, 5898]
      ;cgplot, WL, shift(GAUSS_SMOOTH(expected_flux, 4.9)*3.75e-9, 6), /overplot, color = 'blue'
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
      ;cgplot, WL, Jup_Cal_No_Tell_Corr, color = 'red', /overplot
  
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
          Io_sun_airmass           = 1./cos( (90. - ten(sxpar(Io_header, 'TELAL'))) / !radeg )
          telluric_abs_sun         = interpol(master_telluric_D, WL_D, Io_Wl)
          Io_sun_tell_corr         = Io_sun / telluric_abs_sun^( Io_sun_airmass / Jup_center_airmass )
          telluric_abs_sky         = interpol(master_telluric_D, WL_D, Io_sky_WL)
          Io_sky_tell_corr         = Io_sky / telluric_abs_sky^( Io_sun_airmass / Jup_center_airmass )
          Io_sun_Cal_tell_corr_err = Io_sun_err / telluric_abs_sun^( Io_sun_airmass / Jup_center_airmass ) 
          Io_sky_Cal_tell_corr_err = Io_sky_err / telluric_abs_sky^( Io_sun_airmass / Jup_center_airmass )

        ; Convert everything into R / A units
          Io_sun_Cal_Tell_Corr     = Io_sun_tell_corr / Sensitivity_Curve_D  
          Io_sky_Cal_Tell_Corr     = Io_sky_tell_corr / Sensitivity_Curve_D  
          Io_sun_Cal_tell_corr_err = Io_sun_Cal_tell_corr_err / Sensitivity_Curve_D 
          Io_sky_Cal_tell_corr_err = Io_sky_Cal_tell_corr_err / Sensitivity_Curve_D 
          Io_sun_Cal_No_Tell_Corr  = Io_sun / Sensitivity_Curve_D             
          Io_sky_Cal_No_Tell_Corr  = Io_sky / Sensitivity_Curve_D             

        ; Find the instantaneous Earth-Io Doppler Shift
          cspice_UTC2ET, '2019 24 April ' + sxpar(Io_header, 'UT-OBS'), ET
          ET_mid_exposure = ET + float(sxpar(Io_header, 'EXPTIME'))/2.
          cspice_spkezr, 'Io', ET_mid_exposure, 'J2000', 'LT+S', 'Earth', Io_Earth_State, ltime
          theta  = cspice_vsep(Io_Earth_state[0:2], Io_Earth_state[3:5])
          Io_wrt_Earth_Dopplershift = cos(theta) * norm(Io_Earth_State[3:5])
          SXADDPAR, Io_header, 'Io_DOPPL', Io_wrt_Earth_Dopplershift, 'Io-Earth V_radial in km/s (mid exposure)'
          SXADDPAR, Io_header, 'T_PSHADO', (ET_mid_exposure-PenUmbra_ET) / 60., 'Minutes since Penumbral ingress'
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
          Io_sun_airmass    = 1./cos( (90. - ten(sxpar(Io_header, 'TELAL'))) / !radeg )
          telluric_abs_sun  = interpol(master_telluric_s, WL_s, Io_Wl)
          Io_sun_tell_corr  = Io_sun / telluric_abs_sun^( Io_sun_airmass / Jup_center_airmass )
          telluric_abs_sky  = interpol(master_telluric_s, WL_s, Io_sky_WL)
          Io_sky_tell_corr  = Io_sky / telluric_abs_sky^( Io_sun_airmass / Jup_center_airmass )

        ; do the telluric correction and place things into R/A units  
          Io_sun_Cal_Tell_Corr     = Io_sun_tell_corr / Sensitivity_Curve_s   ; Convert to Rayleighs per Angstrom, Telluric Corrected
          Io_sky_Cal_Tell_Corr     = Io_sky_tell_corr / Sensitivity_Curve_s   ; Convert to Rayleighs per Angstrom, Telluric Corrected
          Io_sun_Cal_tell_corr_err = (Io_sun_err / telluric_abs_sun^( Io_sun_airmass / Jup_center_airmass ) ) / Sensitivity_Curve_s
          Io_sky_Cal_tell_corr_err = (Io_sky_err / telluric_abs_sky^( Io_sun_airmass / Jup_center_airmass ) ) / Sensitivity_Curve_s
          Io_sun_Cal_No_Tell_Corr  = Io_sun / Sensitivity_Curve_s             ; Convert to Rayleighs per Angstrom, no Telluric correction
          Io_sky_Cal_No_Tell_Corr  = Io_sky / Sensitivity_Curve_s             ; Convert to Rayleighs per Angstrom, no Telluric correction

        ; find the instantaneous Earth-Io Doppler Shift
          cspice_UTC2ET, '2019 24 April ' + sxpar(Io_header, 'UT-OBS'), ET
          ET_mid_exposure = ET + float(sxpar(Io_header, 'EXPTIME'))/2.
          cspice_spkezr, 'Io', ET_mid_exposure, 'J2000', 'LT+S', 'Earth', Io_Earth_State, ltime
          theta  = cspice_vsep(Io_Earth_state[0:2], Io_Earth_state[3:5])
          Io_wrt_Earth_Dopplershift = cos(theta) * norm(Io_Earth_State[3:5])
          SXADDPAR, Io_header, 'Io_DOPPL', Io_wrt_Earth_Dopplershift, 'Io-Earth V_radial in km/s (mid exposure)'
          SXADDPAR, Io_header, 'T_PSHADO', (ET_mid_exposure-PenUmbra_ET) / 60., 'Minutes since Penumbral ingress'
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
          Io_ecl_airmass           = 1./cos( (90. - ten(sxpar(Io_header, 'TELAL'))) / !radeg )
          telluric_abs_ecl         = interpol(master_telluric_S, WL_S, Io_Wl)
          Io_ecl_tell_corr         = Io_ecl / telluric_abs_ecl^( Io_ecl_airmass / Jup_center_airmass )
          telluric_abs_sky         = interpol(master_telluric_S, WL_S, Io_sky_WL)
          Io_sky_tell_corr         = Io_sky / telluric_abs_sky^( Io_ecl_airmass / Jup_center_airmass )
          Io_ecl_Cal_tell_corr_err = Io_ecl_err / telluric_abs_ecl^( Io_ecl_airmass / Jup_center_airmass ) 
          Io_sky_Cal_tell_corr_err = Io_sky_err / telluric_abs_sky^( Io_ecl_airmass / Jup_center_airmass )

        ; Convert everything into R / A units
          Io_ecl_Cal_Tell_Corr     = Io_ecl_tell_corr / Sensitivity_Curve_S  
          Io_sky_Cal_Tell_Corr     = Io_sky_tell_corr / Sensitivity_Curve_S  
          Io_ecl_Cal_tell_corr_err = Io_ecl_Cal_tell_corr_err / Sensitivity_Curve_S 
          Io_sky_Cal_tell_corr_err = Io_sky_Cal_tell_corr_err / Sensitivity_Curve_S 
          Io_ecl_Cal_No_Tell_Corr  = Io_ecl / Sensitivity_Curve_S             
          Io_sky_Cal_No_Tell_Corr  = Io_sky / Sensitivity_Curve_S        
          
        ; find the instantaneous Earth-Io Doppler Shift
          cspice_UTC2ET, '2019 24 April ' + sxpar(Io_header, 'UT-OBS'), ET
          ET_mid_exposure = ET + float(sxpar(Io_header, 'EXPTIME'))/2.
          cspice_spkezr, 'Io', ET_mid_exposure, 'J2000', 'LT+S', 'Earth', Io_Earth_State, ltime
          theta  = cspice_vsep(Io_Earth_state[0:2], Io_Earth_state[3:5])
          Io_wrt_Earth_Dopplershift = cos(theta) * norm(Io_Earth_State[3:5])
          SXADDPAR, Io_header, 'Io_DOPPL', Io_wrt_Earth_Dopplershift, 'Io-Earth V_radial in km/s (mid exposure)'
          SXADDPAR, Io_header, 'T_PSHADO', (ET_mid_exposure-PenUmbra_ET) / 60., 'Minutes since Penumbral ingress'
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

        Io_ecl_airmass    = 1./cos( (90. - ten(sxpar(Io_header, 'TELAL'))) / !radeg )
        telluric_abs_ecl  = interpol(master_telluric_S, WL_S, Io_Wl)
        Io_ecl_tell_corr  = Io_ecl / telluric_abs_ecl^( Io_ecl_airmass / Jup_center_airmass )
        telluric_abs_sky  = interpol(master_telluric_S, WL_S, Io_sky_WL)
        Io_sky_tell_corr  = Io_sky / telluric_abs_sky^( Io_ecl_airmass / Jup_center_airmass )

        Io_ecl_Cal_Tell_Corr     = Io_ecl_tell_corr / Sensitivity_Curve_s   ; Convert to Rayleighs per Angstrom, Telluric Corrected
        Io_sky_Cal_Tell_Corr     = Io_sky_tell_corr / Sensitivity_Curve_s   ; Convert to Rayleighs per Angstrom, Telluric Corrected
        Io_ecl_Cal_tell_corr_err = (Io_ecl_err / telluric_abs_ecl^( Io_ecl_airmass / Jup_center_airmass ) ) / Sensitivity_Curve_s
        Io_sky_Cal_tell_corr_err = (Io_sky_err / telluric_abs_sky^( Io_ecl_airmass / Jup_center_airmass ) ) / Sensitivity_Curve_s
        Io_ecl_Cal_No_Tell_Corr  = Io_ecl / Sensitivity_Curve_s             ; Convert to Rayleighs per Angstrom, no Telluric correction
        Io_sky_Cal_No_Tell_Corr  = Io_sky / Sensitivity_Curve_s             ; Convert to Rayleighs per Angstrom, no Telluric correction

        ; find the instantaneous Earth-Io Doppler Shift
          cspice_UTC2ET, '2019 24 April ' + sxpar(Io_header, 'UT-OBS'), ET
          ET_mid_exposure = ET + float(sxpar(Io_header, 'EXPTIME'))/2.
          cspice_spkezr, 'Io', ET_mid_exposure, 'J2000', 'LT+S', 'Earth', Io_Earth_State, ltime
          theta  = cspice_vsep(Io_Earth_state[0:2], Io_Earth_state[3:5])
          Io_wrt_Earth_Dopplershift = cos(theta) * norm(Io_Earth_State[3:5])
          SXADDPAR, Io_header, 'Io_DOPPL', Io_wrt_Earth_Dopplershift, 'Io-Earth V_radial in km/s (mid exposure)'
          SXADDPAR, Io_header, 'T_PSHADO', (ET_mid_exposure-PenUmbra_ET) / 60., 'Minutes since Penumbral ingress'
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
    Torus_params          = REPLICATE({params, brightness:0.0, err_brightness:0.0, mult:0.0, smooth:0.0, shift:0.0, background:'', T_P_shadow:0.0, exptime:0.0}, n_elements(frames))
    
    
    
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
            0, 0, float(sxpar(header, 'Torus_la')), filename_matrix[loc], float(sxpar(header, 'T_PSHADO')), ten(sxpar(header, 'EXPTIME'))*60.}
          O2_fit_params[i]          = {params, A[0]*A[2]*SQRT(2*!DPI), stddev(residual[correl_indicies])*A[2]*SQRT(2*!DPI), $
                                              0, 0, 0, filename_matrix[loc], float(sxpar(header, 'T_PSHADO')), ten(sxpar(header, 'EXPTIME'))*60.}

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
                                     0, 0, 0, filename_matrix[loc], float(sxpar(header, 'T_PSHADO')), ten(sxpar(header, 'EXPTIME'))*60.}
        plot_residual_array[*, i] = residual[plot_WLs]
        DopplerShift_array1[i]    = O1_Io_frame
    endfor

    ; Plot the residual
      shift_By = mean( [transpose( DopplerShift_array1 - Mean(DopplerShift_array1) )], dimension = 1 )
      aligned_residuals = plot_residual_array

      cgplot, WL[plot_WLs], plot_residual_array[*, 0], psym = 0, Ytitle = 'Residual (kR / '+cgsymbol('Angstrom')+')', xtitle = 'Wavelength ('+cgsymbol('Angstrom')+')', $
        yr = YR_resid, xr = XR, /nodata, pos = pos[*,1], /noerase
      ;cgtext, O1_Io_frame - 0.7, YR_resid[1]*0.75, "Io's Doppler Shift", charsize = 1.4, alignment = 0.5
      ;cgtext, O1 + .5, YR_resid[1]*0.75, "Telluric Na", charsize = 1.4, alignment = 0.5
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
    ;YR       = [11, 100]
    YR       = [15, 125]
    ;YR_resid = [-4.9, 29]
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
                                     0, 0, 0, filename_matrix[loc], float(sxpar(header, 'T_PSHADO')), ten(sxpar(header, 'EXPTIME'))*60.}
        plot_residual_array[*, i] = residual[plot_WLs]
        DopplerShift_array1[i]    = Na1_Io_frame
        DopplerShift_array2[i]    = Na2_Io_frame
    endfor

    ; Plot the residual
      shift_By = mean( [transpose( DopplerShift_array1 - Mean(DopplerShift_array1[1:*]) )], dimension = 1 )
      aligned_residuals = plot_residual_array

      cgplot, WL[plot_WLs], plot_residual_array[*, 0], psym = 0, Ytitle = 'Residual (kR / '+cgsymbol('Angstrom')+')', xtitle = 'Wavelength ('+cgsymbol('Angstrom')+')', $
        yr = YR_resid, xr = XR, /nodata, pos = pos[*,1], /noerase
      ;cgtext, Na1_Io_frame - 0.7, YR_resid[1]*0.75, "Io's Doppler Shift", charsize = 1.4, alignment = 0.5
      ;cgtext, Na1 + .5, YR_resid[1]*0.75, "Telluric Na", charsize = 1.4, alignment = 0.5
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
                                     0, 0, 0, '', float(sxpar(header, 'T_PSHADO')), ten(sxpar(header, 'EXPTIME'))*60.}

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

if part eq 12 then begin
  restore, Reduced_Dir+'Na1Na2_fit_params.sav'
  restore, Reduced_Dir+'O2_fit_params.sav'
  restore, Reduced_Dir+'O2_Torus_params.sav'
  

  pos =  [.12,.17,.95,.9]
  cgPS_Open, filename = Reduced_Dir+'Combined_lightcurve.eps', /ENCAPSULATED, xsize = 6., ysize = 4.
  !P.font=1
  device, SET_FONT = 'Helvetica Bold', /TT_FONT

  yr = [0, 7.]
  cgplot, O2_fit_params.T_P_Shadow, O2_fit_params.Brightness, psym = 15, Ytitle = 'Disk-Averaged Brightness (kR)', xtitle = 'Minutes After Ingress', $
    title = 'Io''s Response to Ingress: '+Date, yr = yr, xr = [range1,range2], /nodata, pos = pos

  cgLoadCT, 33, NColors=8, /reverse
  if ingress then x = findgen(Umbra_ET - PenUmbra_ET) / 60.
  colors = reverse(BytScl(x, MIN=min(x), MAX=2.*max(x)))+128
  FOR j=0,n_elements(x)-2 DO BEGIN
    xpoly = [x[j],         x[j], x[j+1],       x[j+1],         x[j]]
    ypoly = [  0., !Y.CRange[1], !Y.CRange[1], 0., 0.]
    cgColorFill, xpoly, ypoly, Color=colors[j]
  ENDFOR
  if ingress then begin
    xpoly = [max(x),     max(x), !X.CRange[1],  !X.CRange[1],  max(x)]
    ypoly = [  0., !Y.CRange[1], !Y.CRange[1], 0., 0.]
    cgColorFill, xpoly, ypoly, color = 'charcoal', /LINE_FILL, orientation = 45., thick = .1
    cgColorFill, xpoly, ypoly, color = 'charcoal', /LINE_FILL, orientation = -45., thick = .1
  endif
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

  Case date of
    'UT180320': cglegend, title = ['[O I] 6300'+cgsymbol('Angstrom'), '[O I] 6364'+cgsymbol('Angstrom'), '[O I] 5577'+cgsymbol('Angstrom'), 'Na D1 + D2'], Psym = [16, 15, 16, 16], charsize = 1.5, $
      bg_color = 'white', color = ['red', 'red', 'green', 'orange'], Location=[0.62, 0.88], /Background, /BOX
    'UT190424': cglegend, title = ['[O I] 6300'+cgsymbol('Angstrom'), 'Na D1 + D2'], Psym = [16, 16], charsize = 1.5, $
      bg_color = 'white', color = ['red', 'orange'], Location=[0.6, 0.88], /Background, /BOX  
    'UT190812': AL_legend, ['K-D in Umbra 190812'], Psym = 4, /right, charsize = 1.5, /clear
  endcase
  cgps_close
  
  ;begin plotting procedure of torus curves in comparison to oxygen 6300 value
  cgPS_Open, filename = Reduced_Dir+'torus_curve.eps', /ENCAPSULATED, xsize = 6., ysize = 4.
  !P.font=1
  device, SET_FONT = 'Helvetica Bold', /TT_FONT
  yr = [0, 7]
  
  
  cgplot, O2_fit_params.T_P_Shadow, O2_fit_params.Brightness, psym = 15, xtitle = 'Minutes After Ingress', $
    title = '[O I] 6300'+cgsymbol('Angstrom') +' Response vs. Plasma Torus Latitude', yr = yr, xr = [range1,55],YTICKFORMAT="(A1)",yticks = 1, /nodata, pos = pos
  cgAxis, YAxis=0, YRange=yr, title= 'Disk-Averaged Brightness [kR]', COLOR = 'red', /Save
  cgplot, O2_fit_params.T_P_Shadow, O2_fit_params.Brightness, psym = 16, /overplot, color = 'red', /err_clip, $
    ERR_YLOW = O2_fit_params.ERR_Brightness, ERR_YHigh = O2_fit_params.ERR_Brightness, ERR_XLOW = O2_fit_params.exptime/2., ERR_XHigh = O2_fit_params.exptime/2.
  cgplot, O2_fit_params.T_P_Shadow, O2_fit_params.Brightness, color = 'red', /overplot
  
  
  restore, 'D:\DATA\Apache Point\Echelle\Io Eclipses\UT180320\Reduced\O2_fit_params.sav'
  restore, 'D:\DATA\Apache Point\Echelle\Io Eclipses\UT180320\Reduced\O2_Torus_params.sav'
  O2_fit_params[0].Brightness = !Values.F_Nan
  Torus_params[0].shift = !Values.F_Nan
  cgplot, O2_fit_params.T_P_Shadow, O2_fit_params.Brightness, psym = 15, /overplot, color = 'red', /err_clip, $
    ERR_YLOW = O2_fit_params.ERR_Brightness, ERR_YHigh = O2_fit_params.ERR_Brightness, ERR_XLOW = O2_fit_params.exptime/2., ERR_XHigh = O2_fit_params.exptime/2.
  cgplot, O2_fit_params.T_P_Shadow, O2_fit_params.Brightness, color = 'red', /overplot
  
  restore, 'D:\DATA\Apache Point\Echelle\Io Eclipses\UT190812\Reduced\O2_fit_params.sav'
  restore, 'D:\DATA\Apache Point\Echelle\Io Eclipses\UT190812\Reduced\O2_Torus_params.sav'
  O2_fit_params[0].Brightness = !Values.F_Nan
  Torus_params[0].shift = !Values.F_Nan
  cgplot, O2_fit_params.T_P_Shadow, O2_fit_params.Brightness, psym = 17, /overplot, color = 'red', /err_clip, $
    ERR_YLOW = O2_fit_params.ERR_Brightness, ERR_YHigh = O2_fit_params.ERR_Brightness, ERR_XLOW = O2_fit_params.exptime/2., ERR_XHigh = O2_fit_params.exptime/2.
  cgplot, O2_fit_params.T_P_Shadow, O2_fit_params.Brightness, color = 'red', /overplot
  
  
  
  restore, 'D:\DATA\Apache Point\Echelle\Io Eclipses\UT180320\Reduced\O2_fit_params.sav'
  restore, 'D:\DATA\Apache Point\Echelle\Io Eclipses\UT180320\Reduced\O2_Torus_params.sav'
  O2_fit_params[0].Brightness = !Values.F_Nan
  Torus_params[0].shift = !Values.F_Nan
  cgAxis, YAxis=1, YRange=[4, 7],title= 'Latitude from Centrifugal Equator [deg]', COLOR = 'purple', /Save
  cgplot, O2_fit_params.T_P_Shadow, Torus_params.shift, psym = 15, /overplot, color = 'purple'
  cgplot, O2_fit_params.T_P_Shadow, Torus_params.shift, /overplot, color = 'purple' 
  
  restore, 'D:\DATA\Apache Point\Echelle\Io Eclipses\UT190812\Reduced\O2_fit_params.sav'
  restore, 'D:\DATA\Apache Point\Echelle\Io Eclipses\UT190812\Reduced\O2_Torus_params.sav'
  O2_fit_params[0].Brightness = !Values.F_Nan
  Torus_params[0].shift = !Values.F_Nan
  cgplot, O2_fit_params.T_P_Shadow, Torus_params.shift, psym = 17, /overplot, color = 'purple'
  cgplot, O2_fit_params.T_P_Shadow, Torus_params.shift, /overplot, color = 'purple'
  
  restore, Reduced_Dir+'O2_fit_params.sav'
  restore, Reduced_Dir+'O2_Torus_params.sav'
  cgplot, O2_fit_params.T_P_Shadow, Torus_params.shift, psym = 16, /overplot, color = 'purple'
  cgplot, O2_fit_params.T_P_Shadow, Torus_params.shift, /overplot, color = 'purple'


  cgps_close
  
  
  
  
  
endif
stop
end