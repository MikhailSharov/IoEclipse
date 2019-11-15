FUNCTION Match_Reference_Spectrum, X, P
    ; Multiply P[0], Add P[1] and Smooth P[2] a reference spectrum (X) until it best matches Y
    return, P[0]*gauss_smooth(X, P[2], /EDGE_TRUNCATE) + P[1]
end

FUNCTION Jup_Shift_Smooth, X, S
    ; Smooth S[0], then Shift S[1] a jupiter scatter spectrum (X) until it best matches Y
    return, interpolate(gauss_smooth(X, S[0], /EDGE_TRUNCATE), findgen(n_elements(X)) - S[1])
end

FUNCTION Match_Jup_Combo, X, P
    ; Smooth P[2], the Shift P[3] some transformed (X) with multiplicative P[0] and offset P[1]
    ; This is a combination of the previous two functions such that they occur at the same time
    return, interpolate(gauss_smooth(P[0]*X + P[1], P[2], /EDGE_TRUNCATE), findgen(n_elements(P[0]*X + P[1])) - P[3])
end

Pro ARCES_redux_v4_190812set

;Written by C. Schmidt, BU Center for Space Physics, 2018
 
;Dir = 'D:\DATA\Apache Point Data\UT180320\'
dir             = 'D:\DATA\Apache Point\Echelle\Io Eclipses\UT190812\
reduced_dir     = 'D:\DATA\Apache Point\Echelle\Io Eclipses\UT190812\Reduced\'
calibration_dir = 'D:\DATA\Apache Point\Echelle\Io Eclipses\UT190812\Reduced\'
  
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

  ; get color versus ingress time
    timeColors = BytScl(findgen(11), Min=0, Max=11, Top=11)
    Color=timeColors[0]
    cgLoadCT, 33, NColors=8, /reverse

  ; Define Rest wavelengths
    O1 = 5577.330
    O2 = 6300.304
    O3 = 6363.776 
    
    O4 = 7771.944  ; Directly produced by e- smashing of SO2, See Ajello et al. (2008)
    O5 = 7774.166
    O6 = 7775.388
    
    O7 = 8446.25   ; Directly produced by e- smashing of SO2, See Ajello et al. (2008)
    O8 = 8446.36 
    O9 = 8446.76
    
    S1 = 9212.865  ; Directly produced by e- smashing of SO2, See Ajello et al. (2008)
    S2 = 9228.092
    S3 = 9237.538

line = 6562.84    ; Fraunhofer H alpha
line = 5270.39    ; Fraunhofer Fe I
line = 5892.7     ; Na


;Eclipse_files = [ 'Io_penumbra.0001' , 'Io_eclipsed.000'+strcompress(indgen(8)+2, /remove_all), 'Io_eclipsed.00'+strcompress(indgen(2)+10, /remove_all) ]
Eclipse_files = ['Io_eclipsed.000'+strcompress(indgen(7)+3, /remove_all), 'Io_eclipsed.00'+strcompress(indgen(2)+10, /remove_all) ]


;----------------------Compare 2 methods to get a decent telluric absorption spectrum---------------------------
  
  ; Method 1, Blue fast rotator
    Telluric_Absorption1 = MRDFITS(Calibration_Dir+'\fullspecHD_159975_BFR.0001.ec.fits', 0, header, /fscale, /silent, /unsigned )
    Telluric_WL          = sxpar(header, 'CRVAL1')+findgen(N_elements(Telluric_Absorption1))*sxpar(header, 'Cdelt1')
    Telluric_AIRMASS     = sxpar(header, 'AIRMASS')

  ; Method 2, A0V star with model spectral reference
    A0V = MRDFITS(Calibration_Dir+'\sumfitcontHD_155379_A0V.0001.ec.fits', 0, header, /fscale, /silent, /unsigned )
    WL  = sxpar(header, 'CRVAL1')+findgen(N_elements(A0V))*sxpar(header, 'Cdelt1')
    READCOL,'D:\DATA\___Calibration___\Solar_and_Stellar_Calibration_Spectra\vegallpr25.30000', F='A,A', Model_A0V_WL, Model_A0V_Flux, Skipline = 700000, numline = 500000, /Silent
    Model_A0V_WL = Model_A0V_WL*10.
    Model_A0V_Flux = float(shift(Model_A0V_Flux, -100))
    Model_A0V_Flux = INTERPOL(Model_A0V_Flux, Model_A0V_WL, WL)
    
    result = ROBUST_POLY_FIT(WL, A0V, 9)
    IRAF_wiggles = A0V / poly(WL, result)
    result = ROBUST_POLY_FIT(WL, Model_A0V_flux, 6)
    Model_A0V_norm = Model_A0V_flux / poly(WL, result)    
    A0V_norm = MRDFITS(Calibration_Dir+'\FullspecHD_155379_A0V.0001.ec.fits', 0, header, /fscale, /silent, /unsigned )
    Telluric_Absorption2 = A0V_norm / (Model_A0V_norm)
    
    window, 0, title = 'Unity Normalized A0V Spectral Model and Measured A0V (HD 155379)'
    cgplot, WL, Model_A0V_norm, color = 'blue', xr = [5570,9300]
    cgplot, WL, A0V_norm, color = 'green', /overplot

    window, 1, title = '2 Telluric absorption spectra methods: Blue fast rotator (black) versus A0V'
    cgplot, Telluric_WL, Telluric_Absorption1, xr = [5000, 6270]
    cgplot, WL, Telluric_Absorption2, /overplot, color = 'red'
    cgplot, WL, IRAF_wiggles, /overplot, color = 'green'
stop
;----------------------get the absolute spectral reflectivity of Jupiter---------------------------

  ; absolute brightness: Data Thief'ed from Woodman et al. 1979.
    READCOL,'C:\IDL\Io\Woodman_et_al_1979_plot1.txt', F='A,A', WL, Albedo, STRINGSKIP = '#', /Silent;wavelength in Angstroms I/F unitless
    READCOL,'C:\IDL\Io\Woodman_et_al_1979_plot2_new.txt', F='A,A', WL_2, Albedo_2, STRINGSKIP = '#', /Silent ;wavelength in Angstroms I/F unitless
    Woodman_WL = float([WL, WL_2]) ;STITCH THESE TOGETHER
    Woodman_Albedo = Float([albedo , albedo_2]);STITCH THESE TOGETHER
  
  ; absolute brightness: from Karkoschka (1998) Icarus on the PDS as ID # ESO-J/S/N/U-SPECTROPHOTOMETER-4-V1.0
    READCOL,'C:\IDL\Io\Karkoschka_1995low.tab', F='X,A,X,A', Karkoschka_WL, Karkoschka_Albedo, STRINGSKIP = '#', /Silent ;wavelength in nm I/F unitless
    READCOL,'C:\IDL\Io\Karkoschka_1995high.tab', F='X,A,X,A', Karkoschka_WL_2, Karkoschka_Albedo_2, STRINGSKIP = '#', /Silent ;wavelength in Angstroms I/F unitless
    Karkoschka_WL = float([Karkoschka_WL, Karkoschka_WL_2]) ;STITCH THESE TOGETHER
    Karkoschka_Albedo = Float([Karkoschka_albedo, Karkoschka_albedo_2]) ;STITCH THESE TOGETHER
    Karkoschka_Albedo = Karkoschka_Albedo[sort(Karkoschka_WL)]
    Karkoschka_WL = Karkoschka_WL[sort(Karkoschka_WL)]
  
  ; compare the two
    cgplot, Woodman_WL / 10., Woodman_Albedo, color = 'blue', xstyle = 1., psym = 3, Xtitle = 'Wavelength (nm)', $
      ytitle = 'I/F Reflectivity'
    cgplot, Karkoschka_WL, Karkoschka_Albedo*1.35, color = 'red', /overplot
    cgtext, 340, .1, 'EQUATOR AT CENTRAL MERIDIAN (Woodman et al. 1979)', color = 'blue'
    cgtext, 340, .16, 'FULL DISK scaled by 1.35 (Karkoschka 1998)', color = 'red'

  ; make an informed choice via scaling
    Karkoschka_wl = Karkoschka_wl * 10. ;nm to Angstroms
    Karkoschka_Albedo = Karkoschka_Albedo*1.35 ;Scale the "Full disk albedo" to the "Central Meridian Equatorial Absolute Reflectivity"
  
  ; Woodman only gives I/F to get the output light from Jupiter, multiply this by a solar spectrum
    READCOL,'C:\IDL\Io\Kurucz_2005_irradthuwl.dat', F='A,A', WL_nm, flux, STRINGSKIP = '#', /Silent ;flux is in W/m2/nm
    start  = where(WL_nm eq '299.100')
    WL_nm = float(WL_nm[start:*])
    flux = float(flux[start:*])
    ; change flux units from W/m^2/nm to photons/(cm^2 s A)
    ; multiply by ((lambda / hc) / 1 W)  * (1 m^2 / 1e4 cm^2) * (1 nm / 10 A)
    conversion = ((WL_nm*1.e-9)/(6.62606957e-34*299792458.D)) * (1./1.e4) * (1./10.)
    flux = flux * conversion
    WL_A = temporary(WL_nm) * 10. ;wavelength from nm into angstroms
    VACTOAIR, WL_A, WL_A_Air ;Vacuum to air wavelength conversion
    WL_A = temporary(WL_A_Air)
    ;cross-checked this result against Huebner et al. 1992.
    
    Jupiter_center_header = headfits(reduced_dir+'sumfitcontJupiter_Disk_Center.0001.ec.fits')
    cspice_UTC2ET, sxpar(Jupiter_center_header, 'DATE-OBS'), ET

  ; scale solar flux to Jupiter's instantaneous distance
    cspice_spkezr, 'Jupiter', ET, 'J2000', 'LT+S', 'Sun', Jupiter_Sun_State, ltime
    cspice_spkezr, 'Jupiter', ET, 'J2000', 'LT+S', 'Earth', Jupiter_Earth_State, ltime
    solar_distance = norm(Jupiter_Sun_State[0:2]) / 149597871.
    flux_at_jupiter = flux / solar_distance^2.
    Albedo = INTERPOL(Karkoschka_Albedo, Karkoschka_WL, WL_A)
    Rayleighs_per_angstrom = 4.*flux_at_jupiter*albedo / 1.e6
    window, 0, Title = 'Instantaneous Rayleighs per Angstrom: Center of Jupiter''s Disk'
    plot, WL_A, Rayleighs_per_angstrom, xr = [5885, 5900], charsize = 2 ;compare to 5.5 MR per angstrom by Brown & Schneider, 1981
  
  ; Adjust the expected absolute flux for jupiter's instantaneous Doppler shift
    theta  = cspice_vsep(Jupiter_Earth_State[0:2], Jupiter_Earth_State[3:5])
    Dopplershift = cos(theta) * norm(Jupiter_Earth_State[3:5])  ; scalar projection of the relative velocity along the line of sight
    WL_A = WL_A + WL_A*Dopplershift/cspice_clight()

  ; get the counts per pixel on Jupiter's disk
    cont_fit    = MRDFITS(reduced_dir+'sumfitcontJupiter_Disk_Center.0001.ec.fits', 0, header, /fscale, /silent, /unsigned )
    jup_center_forfit  = MRDFITS(reduced_dir+'sumJupiter_Disk_Center.0001.ec.fits', 0, header, /fscale, /silent, /unsigned ) ;sum of jupitercenter, used for poly fit
    jup_center = MRDFITS(reduced_dir+'fullspecJupiter_Disk_Center.0001.ec.fits', 0, header, /fscale, /silent, /unsigned )    ;fullspec, normalized to one
    WL          = sxpar(header, 'CRVAL1')+findgen(N_elements(jup_center))*sxpar(header, 'Cdelt1')
    jup_center_err =  MRDFITS(reduced_dir+'sig_fullspecJupiter_Disk_Center.0001.ec.fits', 0, header, /fscale, /silent, /unsigned ) ;errors, normalized to one
    jup_center= jup_center* poly(WL, ROBUST_POLY_FIT(WL, jup_center_forfit, 11))                                                       ;fullspec unnormalized using poly fit
    jup_center_err= jup_center_err* poly(WL, ROBUST_POLY_FIT(WL, jup_center_forfit, 11))                                               ;fullspec errors unnormalized using poly fit
    raw_header  = headfits(dir+'Jupiter_Disk_Center.0001.fits')
    cont_fit    = cont_fit / sxpar(raw_header, 'EXPTIME')
    jup_center  = jup_center / sxpar(raw_header, 'EXPTIME')
    jup_center_err = jup_center_err / sxpar(raw_header, 'EXPTIME')
    
    expected_flux = interpol(Rayleighs_per_angstrom, WL_A, WL)   ; move expected flux to the jovian doppler-shift
    smoothed_expected_flux = GAUSS_SMOOTH(expected_flux, 1.9)    ; this smoothing looks about right for APO/ARCES
    
    Sensitivity = jup_center / smoothed_expected_flux            ; Sensetivity in (DN / S) / (R / A)
    Cont_Sensitivity = cont_fit / smoothed_expected_flux         ; Sensetivity in (DN / S) / (R / A)
 
    aligned_telluric_absorption = interpol(Telluric_Absorption1, Telluric_WL, WL) 
    telluric_fitting_ind        = where( (aligned_telluric_absorption lt 0.97) and (((wl gt 6302.4) and (wl lt 6308.)) or ((wl gt 6287.) and (wl lt 6290.))), /NULL)

        parinfo            = replicate( {value:0., fixed: 0b, limited: [0b,0b], limits: fltarr(2) }, 4)
        parinfo[0].limited = 1b                                 ; limit differences in multiplier amplitude
        parinfo[0].limits  = [-1.e3, 1.e7]
        parinfo[1].limited = 1b                                 ; limit additive differences
        parinfo[1].limits  = [-1.e7, 1.e7]
        parinfo[2].limited   = 1b                                 ; fix differences in smooth width
        parinfo[2].limits  = [-1.e7, 1.e7]
        parinfo[3].limited   = 1b                                 ; fix differences in smooth width
        parinfo[3].limits  = [-1.e7, 1.e7]
    
        X = aligned_telluric_absorption
        Y = jup_center / Sensitivity
        Err_Y = jup_center_err / Sensitivity ;propogated from previous import

    p0 = [1.e6, 0., 0., 0.] ; guess at initial Multiply P[0], Add P[1] and Smooth P[2] to
    p = mpfitfun('Match_Reference_Spectrum', x[telluric_fitting_ind], y[telluric_fitting_ind], err_y[telluric_fitting_ind], p0, /NaN, status=status, parinfo = parinfo) ; Fit a y = mx + b function to the spectrum at every spatial bin... SLOW

    telluric_fit = P[0]*X + P[1]
    Jup_Cal_Tell_Corr = Y / (telluric_fit / median(telluric_fit)) ; Jupiter Spectrum Calibrated into R / A and corrected for tellutic absorption
    Jup_Cal_Tell_Corr_err = Y / (telluric_fit / median(telluric_fit))
    window, 1, xs = 1800, ys = 800
    cgplot, WL, Y, xr = [6275,6325], color = 'blue', ytitle = 'Rayleighs per Angstrom'
    cgplot, WL[telluric_fitting_ind], Y[telluric_fitting_ind], /overplot, color = 'green', psym=14
    cgplot, WL, telluric_fit, /overplot, color = 'red'
    cgplot, WL, Jup_Cal_Tell_Corr, /overplot, color = 'black'

    Penumbra_UTC = '2019-Aug-12 03:44:40'
    Umbra_UTC    = '2019-Aug-12 03:48:19'
    cspice_UTC2ET, PenUmbra_UTC, PenUmbra_ET
    cspice_UTC2ET, Umbra_UTC, Umbra_ET
    
    window, 3
  ; Correct tellurics and write the spectra into R/A units...   
    for i = 0, n_elements(Eclipse_files)-1 do begin 
      spec_forfit = MRDFITS(reduced_dir+'sum'+Eclipse_files[i]+'.ec.fits', 0, header, /fscale, /unsigned )
   
      spec_fullspec = MRDFITS(reduced_dir+'fullspec'+Eclipse_files[i]+'.ec.fits', 0, full_header, /fscale, /unsigned)
      spec = spec_fullspec*poly(WL, ROBUST_POLY_FIT(WL, spec_forfit, 11))
      
      spec_err = MRDFITS(reduced_dir+'sig_fullspec'+Eclipse_files[i]+'.ec.fits', 0, sig_header, /fscale, /unsigned )
      spec_err = spec_err*poly(WL, ROBUST_POLY_FIT(WL, spec_forfit, 11))
      
      raw_header  = headfits(dir+Eclipse_files[i]+'.fits')
      spec = spec / sxpar(raw_header, 'EXPTIME')
      spec_err = spec_err / sxpar(raw_header, 'EXPTIME')
      WL   = sxpar(header, 'CRVAL1')+findgen(N_elements(spec))*sxpar(header, 'Cdelt1')
 
      Y = spec / Sensitivity
      Err_Y = spec_err / Sensitivity  ;propogated from sig_ files
      p = mpfitfun('Match_Reference_Spectrum', x[telluric_fitting_ind], y[telluric_fitting_ind], err_y[telluric_fitting_ind], p0, /NaN, status=status, parinfo = parinfo, /quiet) ; Fit a y = mx + b function to the spectrum at every spatial bin... SLOW
      telluric_fit = P[0]*X + P[1]
      
      Spec_Cal_Tell_Corr = Y / (telluric_fit / median(telluric_fit)) ; Jupiter Spectrum Calibrated into R / A and corrected for tellutic absorption
      Spec_Cal_Tell_Corr_err = Err_Y / (telluric_fit / median(telluric_fit)) ; Propogation of errors into R/A and corrected for tellucric absorption
      ; Inspect
        cgplot, WL, Y, xr = [6298,6303], color = 'blue', ytitle = 'Rayleighs per Angstrom'
        cgplot, WL[telluric_fitting_ind], Y[telluric_fitting_ind], /overplot, color = 'green', psym=14
        cgplot, WL, telluric_fit, /overplot, color = 'red'
        cgplot, WL, Spec_Cal_Tell_Corr, /overplot, color = 'black'  
      ; find the instantaneous Earth-Io Doppler Shift
        cspice_UTC2ET, sxpar(raw_header, 'DATE-OBS'), ET
        ET_mid_exposure = ET + float(sxpar(raw_header, 'EXPTIME'))/2.
        cspice_spkezr, 'Io', ET_mid_exposure, 'J2000', 'LT+S', 'Earth', Io_Earth_State, ltime
        theta  = cspice_vsep(Io_Earth_state[0:2], Io_Earth_state[3:5])
        Io_wrt_Earth_Dopplershift = cos(theta) * norm(Io_Earth_State[3:5]) 
        SXADDPAR, Header, 'Io_DOPPL', Io_wrt_Earth_Dopplershift, 'Io-Earth V_radial in km/s (mid exposure)'
        SXADDPAR, Header, 'T_SHADOW', (ET_mid_exposure-PenUmbra_ET) / 60., 'Minutes since Penumbral ingress'

      print, sxpar(header, 'NAXIS2')
      MWRFITS, Spec_Cal_Tell_Corr, reduced_dir+'R_per_A_'+Eclipse_files[i]+'.ec.fits', header, /CREATE ;/create overwrites
      MWRFITS, Spec_Cal_Tell_Corr_err, reduced_dir+'sig_R_per_A_'+Eclipse_files[i]+'.ec.fits', sig_header, /CREATE, /silent ;/create overwrites for sig values (not fancy mate :( will be fixed later)
    endfor
  ;--------------------------------------------------------Get the O 6300 Time Series---------------------------------------------------------------------------------------------

  ; Correct Jovian Scattered Light
    shift_array   = [2.9, 0, 0, 0, 0, 0, 0]  
    Smooth_array  = [0, 0, 0, 0, 0, 0, 0]  
    
    parinfo[0].limits      = [0., 1.e7]
    Gaussian_parinfo       = replicate({value:0., fixed:0, limited:[0,0], limits:[0.,0.]}, 3)    
    
    include_WLs            = where( ((wl gt 6293.) and (wl lt 6308.)), /NULL)
    exclude_WLs            = where( ((wl gt 6300.3) and (wl lt 6300.9)), /NULL)
    scatter_fitting_ind    = cgSetDifference(include_WLs, exclude_WLs)
    LSF_fitting_ind        =  where( ((wl gt 6300.4) and (wl lt 6300.89)), /NULL)
    
    
    residual_array         = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
    LSF_Fit_array          = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
    WL_array               = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
    Brightness_array       = fltarr(n_elements(Eclipse_files))
    err_Brightness_array   = fltarr(n_elements(Eclipse_files))
    T_Shadow_array         = fltarr(n_elements(Eclipse_files))
    EXPTime_array          = fltarr(n_elements(Eclipse_files))
    DopplerShift_array     = fltarr(n_elements(Eclipse_files))
    window, 0, xs = 800, ys = 1000
    cgplot, WL, spec, xr = [6296,6302], yr = [2.e4, 5.e4], ytitle = 'R/A', /nodata
    ;cgplot, WL, spec, xr = [5888,5898], yr = [1.e3, 1.e5], ytitle = 'R/A', /nodata
    for i = 0, n_elements(Eclipse_files)-1 do begin
      spec = MRDFITS(reduced_dir+'R_per_A_'+Eclipse_files[i]+'.ec.fits', 0, header, /fscale, /unsigned, /silent)
      spec_err = MRDFITS(reduced_dir+'sig_R_per_A_' +Eclipse_files[i]+'.ec.fits', 0, sig_header, /fscale, /unsigned, /silent)
      WL   = sxpar(header, 'CRVAL1')+findgen(N_elements(spec))*sxpar(header, 'Cdelt1')
      raw_header  = headfits(dir+Eclipse_files[i]+'.fits')
      cgplot, WL, spec+i*4.e3, Color=timeColors[i], ytitle = 'R/A', /overplot, thick =3
      
      X = Jup_Cal_Tell_Corr
      X_err = Jup_Cal_Tell_Corr_err
      Y = spec
      Err_Y = spec_err
      
      s0 = [10., 10., 1.1, 1.2]
      
      p = mpfitfun('Match_Jup_Combo', X[scatter_fitting_ind], Y[scatter_fitting_ind], spec_err[scatter_fitting_ind], s0, /NaN, status=status_smooth, parinfo = parinfo, /verbose)
      scatter_fit = interpolate(gauss_smooth(P[0]*X + P[1], P[2], /EDGE_TRUNCATE), findgen(n_elements(P[0]*X + P[1])) - P[3])
  
      ;p = mpfitfun('Match_Reference_Spectrum', x[scatter_fitting_ind], y[scatter_fitting_ind], err_y[scatter_fitting_ind], p0, /NaN, status=status, parinfo = parinfo, /quiet) 
        ;scatter_fit = P[0]*X + P[1]
      cgplot, WL, scatter_fit + i*4.e3, /OVERPLOT, COLOR = timeColors[i], linestyle =3 ; offset each plot

      residual = y - scatter_fit
      residual_err = err_y ;real replacement for hacked_errors, not sure if anything from scatter_fit should be subtracted.

      O2_Io_frame                       = O2 + O2 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      ;Gaussian_parinfo[1].fixed        = 1
      ;Gaussian_parinfo[1].value        = O2_Io_frame
      Gaussian_parinfo[1].limited       = 1b
      Gaussian_parinfo[1].limits        = [6300.4, 6300.89]
      ;Gaussian_parinfo[2].limited       = 1b
      ;Gaussian_parinfo[2].limits        = [0.07, 0.15]
      ;Gaussian_parinfo[2].fixed        = 1b
      ;Gaussian_parinfo[2].value        = 0.09

      ;if i eq 0 then hacked_errors = replicate(1000., n_elements(LSF_fitting_ind)) else $
                     ;hacked_errors = replicate(400., n_elements(LSF_fitting_ind))  ;hack
      LSF_Fit = mpfitpeak(WL[LSF_fitting_ind], residual[LSF_fitting_ind], a, PERROR = err_a, /POSITIVE, PARINFO = Gaussian_parinfo, $
                          STATUS = STATUS, ERRORS = abs(residual_err[LSF_fitting_ind]), nterms = 3, /quiet)
        WL_array[*, i]              = WL[include_WLs]
        residual_array[*, i]        = residual[include_WLs]
        LSF_Fit_array[*, i]         = gaussian(WL[include_WLs], a) ;LSF_Fit
        Brightness_array[i]         = $   ;  total(residual[LSF_fitting_ind])
                                      A[0]*0.09*SQRT(2*!DPI)                                                      ; HACK HAck HAck fix the linewidth used in the fitting
                                          ;  A[0]*A[2]*SQRT(2*!DPI)
        err_Brightness_array[i]     = abs(Brightness_array[i] * sqrt( (err_A[0]/A[0])^2. + (err_A[2]*A[2])^2.) ) ; uncertainties in width and height add in quadrature.
        ;err_Brightness_array = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        T_Shadow_array[i]           = sxpar(header, 'T_SHADOW')
        EXPTime_array[i]            = sxpar(raw_header, 'EXPTIME') / 60.
        DopplerShift_array[i]       = O2_Io_frame
      print, 'linewidth = ', a[2]
    endfor

stop
    window, 2, xs = 800, ys = 1200
    cgplot, WL, spec, xr = [6299,6301], yr = [-1.e3, 7.e3], ytitle = 'Residual Jovian scattered light subtraction R/A', /nodata
    for i = 0, n_elements(Eclipse_files)-1 do begin
      cgplot, WL_array[*, i], residual_array[*, i], /OVERPLOT, COLOR = timeColors[i], psym=10
      cgplot, WL_array[*, i], LSF_Fit_array[*, i], /OVERPLOT, COLOR = timeColors[i]
      cgplot, [6300.4, 6300.4], [-10000., 10000], linestyle = 1, /overplot
      cgplot, [6300.89, 6300.89], [-10000., 10000], linestyle = 1, /overplot
      cgplot, [DopplerShift_array[i],DopplerShift_array[i]], [-10000., 10000], COLOR = timeColors[i], /overplot
    endfor
    
    pos =  [.17,.2,.98,.9]
    cgPS_Open, filename = Reduced_Dir+'O6300_lightcurve_V3.eps', /ENCAPSULATED, xsize = 10., ysize = 6.5
    
      !P.font=1
      device, SET_FONT = 'Helvetica Bold', /TT_FONT
      !p.charsize = 3.
      loadct, 0
      cgplot, T_Shadow_array, Brightness_array, psym = 15, Ytitle = 'Rayleighs (5.1 sq arcsec aperture)', xtitle = 'Minutes After Ingress', title = 'Io''s Oxygen 6300'+cgsymbol('Angstrom') + ' Response to Eclipse',  $
              ERR_YLOW = ERR_Brightness_array, ERR_YHigh = ERR_Brightness_array, ERR_XLOW = EXPTime_array/2., ERR_XHigh = EXPTime_array/2., yr = [0, 1700], xr = [70,140], /nodata, pos = pos
  
      x = findgen(Umbra_ET - PenUmbra_ET) / 60.
      colors = reverse(BytScl(x, MIN=min(x), MAX=2.*max(x)))+128 
               ;reverse(BytScl(x, MIN=min(x), MAX=max(x))) 
      FOR j=0,n_elements(x)-2 DO BEGIN
        xpoly = [x[j],         x[j], x[j+1],       x[j+1],         x[j]]
        ypoly = [  0., !Y.CRange[1], !Y.CRange[1], 0., 0.]
        cgColorFill, xpoly, ypoly, Color=colors[j]
      ENDFOR
      xpoly = [70,     70, !X.CRange[1],  !X.CRange[1],  70]
      ypoly = [  0., !Y.CRange[1], !Y.CRange[1], 0., 0.]
      cgColorFill, xpoly, ypoly, color = 'charcoal', /LINE_FILL, orientation = 45., thick = .1
      cgColorFill, xpoly, ypoly, color = 'charcoal', /LINE_FILL, orientation = -45., thick = .1              
      cgtext, 2.5, !Y.CRange[1]/6., 'Penumbral Eclipse', orientation = 90., color = 'white' 
      cgplot, T_Shadow_array, Brightness_array, psym = 4, /noerase, symsize = 2, yr = [0, 1700], xr = [70,140], pos = pos, $
        ERR_YLOW = ERR_Brightness_array, ERR_YHigh = ERR_Brightness_array, ERR_XLOW = EXPTime_array/2., ERR_XHigh = EXPTime_array/2. 
    cgps_close
    cgLoadCT, 33, NColors=8, /reverse
    
stop
    ;--------------------------------------------------------Get the Na Time Series---------------------------------------------------------------------------------------------



    spec = MRDFITS(reduced_dir+'fullspecIo_eclipsed.0002.ec.fits', 0, header, /fscale, /silent, /unsigned )
    WL   = sxpar(header, 'CRVAL1')+findgen(170105)*sxpar(header, 'Cdelt1')
    
    include_WLs            = where( ((wl gt 5888.) and (wl lt 5898.)), /NULL)
    exclude_WLs            = where( ((wl gt 5889.25) and (wl lt 5889.8)), /NULL)
    exclude_WLs2            = where( ((wl gt 5895.25) and (wl lt 5895.8)), /NULL)
    scatter_fitting_ind    = cgSetDifference(include_WLs, exclude_WLs)
    scatter_fitting_ind    = cgSetDifference(scatter_fitting_ind, exclude_WLs2)
    LSF_fitting_ind        =  where( ((wl gt 5889) and (wl lt 5890)), /NULL)
    
    parinfo[0].limits  = [0., 1.e7]
    
    wl_range      = 9.5 ;angstroms
    shift_array   = [0, 1.5, 1.75, 2., 2.25, 2.5, 2.75,1,1,1,1,1]
    Smooth_array  = [0,  15,  14,  14,  14,  14,  14,1,1,1,1,1]  ; Good for 5893
    Smooth_array  = [0,  15,  14,  14,  14,  12,  12,1,1,1,1,1]  ; Good for 5893
    

    ;pos =  [.17,.2,.98,.9]
  cgPS_Open, filename = Reduced_Dir+'Na_Reduction.eps', /ENCAPSULATED, xsize = 10., ysize = 6.5
    !P.font=1
    device, SET_FONT = 'Helvetica Bold', /TT_FONT
    !p.charsize = 3.
      
    pos = cgLayout([1,2], YGap=0, OXMargin=[5,5], OyMargin=[10,10])
    jup = MRDFITS(reduced_dir+'fullspecJupiter_Disk_Center.0001.ec.fits', 0, header, /silent, /fscale, /unsigned )
    junk = [min(abs(wl - (line-wl_range/2)), low_ind), min(abs(wl - (line+wl_range/2)), hi_ind)]
    cgplot, wl, spec, xr = [Line-wl_range/2, line+wl_range/2], /nodata, pos = pos[*,0], yr = [-.4, 1.], xtickformat = '(A1)', title = 'Removing Jupiter Contaimination - Sodium', ytickformat = '(A1)' 
    for i = 1, 6 do begin
      spec = MRDFITS(reduced_dir+'fullspecIo_eclipsed.000'+strcompress(i, /remove_all)+'.ec.fits', 0, header, /silent, /fscale, /unsigned )
      cgplot, wl, spec-.15*(i-1), /overplot, Color=timeColors[i], thick = 2, pos = pos[*,0]
      
      s0 = [5., 2.5, 0., 0.]
      s = mpfitfun('Jup_Shift_Smooth', jup[scatter_fitting_ind], spec[scatter_fitting_ind], err_y[scatter_fitting_ind], s0, /NaN, status=status, parinfo = parinfo, /verbose)
      shifted = interpolate(gauss_smooth(jup, S[0], /EDGE_TRUNACATE), findgen(n_elements(Jup)) - S[1]) - .15*(i-1)
      
      ;shifted = interpolate(smooth(jup, smooth_array[i]), findgen(n_elements(Jup)) - shift_array[i]) - .15*(i-1)
      cgplot, wl, shifted, /overplot, Color=timeColors[i], thick = 2, pos = pos[*,0]  ; sub-pixel shifting
      ;cgplot, wl, smart_shift(smooth(jup, smooth_array[i]), shift_array[i], /Interp)-.15*i, /overplot, Color=timeColors[i], thick = 2, pos = pos[*,0]
    endfor
    
    spec_forfit = MRDFITS(reduced_dir+'sum'+Eclipse_files[0]+'.ec.fits', 0, header, /fscale, /unsigned )
    spec_fullspec = MRDFITS(reduced_dir+'fullspec'+Eclipse_files[0]+'.ec.fits', 0, header, /fscale, /unsigned )
    spec = spec_fullspec*poly(WL, ROBUST_POLY_FIT(WL, spec_forfit, 10))
    
    
;    raw_header  = headfits(dir+Eclipse_files[i]+'.fits')
;    spec = spec / sxpar(raw_header, 'EXPTIME')
;    WL   = sxpar(header, 'CRVAL1')+findgen(N_elements(spec))*sxpar(header, 'Cdelt1')
;    Y = spec / Cont_Sensitivity
    
    cgplot, wl, spec, xr = [Line-wl_range/2, line+wl_range/2], /nodata, yr = [-.1, .3], pos = pos[*,1], ytickformat = '(A1)', /noerase, ytitle = cgsymbol('Angstrom') ;, yr = [-.6, 3.e6]
    cgplot, [5889.95095, 5889.95095], [-2., 2], linestyle = 1, /overplot
    cgplot, [5895.92424, 5895.92424], [-2., 2], linestyle = 1, /overplot
    for i = 2, n_elements(Eclipse_files)-1 do begin
;      spec = MRDFITS(reduced_dir+'sum'+Eclipse_files[i]+'.ec.fits', 0, header, /fscale, /unsigned )
;      cont = MRDFITS(reduced_dir+'sumfitcont'+Eclipse_files[i]+'.ec.fits', 0, header, /fscale, /unsigned )
;      ;spec = spec/cont
;      spec = spec/Cont_Sensitivity
;      WL   = sxpar(header, 'CRVAL1')+findgen(N_elements(spec))*sxpar(header, 'Cdelt1')
;      cgplot, wl, Spec-0.*i, /overplot, Color=timeColors[i], thick = 2, pos = pos[*,1]     
      
      
      ; find and plot the intantaneous Earth-Io Doppler Shift
      cspice_UTC2ET, sxpar(header, 'DATE-OBS'), ET
      cspice_spkezr, 'Io', ET + float(sxpar(header, 'EXPTIME'))/2. , 'J2000', 'LT+S', 'Earth', Io_Earth_State, ltime
      theta  = cspice_vsep(Io_Earth_state[0:2], Io_Earth_state[3:5])
      Io_wrt_Earth_Dopplershift = line * cos(theta) * norm(Io_Earth_State[3:5]) / cspice_clight()
      Io_wrt_Earth_Dopplershift = line * cos(theta) * sqrt(Io_Earth_State[3]^2.+Io_Earth_State[4]^2.+Io_Earth_State[5]^2.) / cspice_clight()
      cgplot, [5889.95095, 5889.95095]+Io_wrt_Earth_Dopplershift, [-2., 2], Color=timeColors[i], linestyle = 2, /overplot
      cgplot, [5895.92424, 5895.92424]+Io_wrt_Earth_Dopplershift, [-2., 2], Color=timeColors[i], linestyle = 2, /overplot
      
      spec     = MRDFITS(reduced_dir+'fullspecIo_eclipsed.000'+strcompress(i, /remove_all)+'.ec.fits', 0, header, /silent, /fscale, /unsigned )
      residual = spec - interpolate(smooth(jup, smooth_array[i]), findgen(n_elements(Jup)) - shift_array[i])
      cgplot, wl, residual, /overplot, Color=timeColors[i], thick = 2, pos = pos[*,1]
      
      cgtext, 5895.1, .1,  'Io',charsize = 2.
      cgtext, 5896.1, .2, 'Telluric',charsize = 2.
   endfor 
 cgps_close
   stop

    ; Correct Jovian Scattered Light
    ;shift_array   = [1.5, 1, 1.5, 1.5, 1.5, 1.5, 1.6]
    ;Smooth_array  = [10,  14,  14,  14,  14,  14,  14]

    parinfo[0].limits      = [0., 1.e6]
    Gaussian_parinfo       = replicate({value:0., fixed:0, limited:[0,0], limits:[0.,0.]}, 3)

    include_WLs            = where( ((wl gt 5888.) and (wl lt 5898.)), /NULL)
    exclude_WLs            = where( ((wl gt 5890.05) and (wl lt 5890.5)), /NULL)
    exclude_WLs2            = where( ((wl gt 5896.) and (wl lt 5896.52)), /NULL)
    scatter_fitting_ind    = cgSetDifference(include_WLs, exclude_WLs)
    scatter_fitting_ind    = cgSetDifference(scatter_fitting_ind, exclude_WLs2)
    LSF_fitting_ind        =  where( ((wl gt 5890.05) and (wl lt 5890.5)), /NULL)

    residual_array         = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
    LSF_Fit_array          = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
    WL_array               = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
    Brightness_array       = fltarr(n_elements(Eclipse_files))
    err_Brightness_array   = fltarr(n_elements(Eclipse_files))
    T_Shadow_array         = fltarr(n_elements(Eclipse_files))
    EXPTime_array          = fltarr(n_elements(Eclipse_files))
    DopplerShift_array     = fltarr(n_elements(Eclipse_files))
    window, 0, xs = 800, ys = 1000
    cgplot, WL, spec, xr = [5882.,5902.], yr = [5.e3, 4.e4], ytitle = 'R/A', /nodata
    ;cgplot, WL, spec, xr = [6556.,6568.], yr = [5.e3, 4.e4], ytitle = 'R/A', /nodata
    for i = 0, n_elements(Eclipse_files)-1 do begin
      spec = MRDFITS(reduced_dir+'R_per_A_'+Eclipse_files[i]+'.ec.fits', 0, header, /fscale, /unsigned, /silent)
      spec_err = MRDFITS(reduced_dir+'sig_R_per_A_'+Eclipse_files[i]+'.ec.fits', 0, sig_header, /fscale, /unsigned, /silent)
      WL   = sxpar(header, 'CRVAL1')+findgen(N_elements(spec))*sxpar(header, 'Cdelt1')
      raw_header  = headfits(dir+Eclipse_files[i]+'.fits')
      cgplot, WL, spec+i*5.e3, Color=timeColors[i], ytitle = 'R/A', /overplot, thick =3
      
      X = Jup_Cal_Tell_Corr
      Y = spec
      
      ;p = mpfitfun('Match_Reference_Spectrum', x[scatter_fitting_ind], y[scatter_fitting_ind], err_y[scatter_fitting_ind], p0, /NaN, status=status_specfit, parinfo = parinfo, /quiet)
      ;scatter_fit = P[0]*X + P[1]
      
      s0 = [0.003, -300., 2.0, 2.1]
      ;s = mpfitfun('Jup_Shift_Smooth', scatter_fit[scatter_fitting_ind], spec[scatter_fitting_ind], spec_err[scatter_fitting_ind], s0, /NaN, status=status_smooth, parinfo = parinfo, /verbose)
      ;scatter_fit = interpolate(gauss_smooth(scatter_fit, S[0]), findgen(n_elements(scatter_fit)) - S[1])
      ;X = interpolate(smooth(Jup_Cal_Tell_Corr, smooth_array[i]), findgen(n_elements(Jup_Cal_Tell_Corr)) - shift_array[i])   using guessed scatter
      Err_Y = spec_err
      p = mpfitfun('Match_Jup_Combo', X[scatter_fitting_ind], Y[scatter_fitting_ind], spec_err[scatter_fitting_ind], s0, /NaN, status=status_smooth, parinfo = parinfo, /verbose)
      scatter_fit = interpolate(gauss_smooth(P[0]*X + P[1], P[2], /EDGE_TRUNCATE), findgen(n_elements(P[0]*X + P[1])) - P[3])
      ;p = mpfitfun('Match_Jup_Combo', X[include_WLs], Y[include_WLs], spec_err[include_WLs], s0, /NaN, status=status_smooth, parinfo = parinfo, /verbose)
      ;scatter_fit_2 = interpolate(gauss_smooth(P[0]*X + P[1], P[2], /EDGE_TRUNCATE), findgen(n_elements(P[0]*X + P[1])) - P[3])
      
      ;scatter_fit = (scatter_fit + i/2.5*scatter_fit_2)/(1.+i/2.5)
      
      ;s = mpfitfun('Jup_Shift_Smooth', scatter_fit[scatter_fitting_ind], Y[scatter_fitting_ind], spec_err[scatter_fitting_ind], s0, /NaN, status=status_smooth, parinfo = parinfo, /verbose)
      ;scatter_fit = interpolate(gauss_smooth(scatter_fit, S[0], /EDGE_TRUNCATE), findgen(n_elements(scatter_fit)) - S[1])
      
      cgplot, WL, scatter_fit + i*5.e3, /OVERPLOT, COLOR = timeColors[i], linestyle =3 ; offset each plot
      cgplot, [5890.05, 5890.05] , [0, 1e6], /overplot ; plotting the edge cases of the scatter fitting indicies
      cgplot, [5890.5, 5890.5] , [0, 1e6], /overplot
      cgplot, [5896., 5896.] , [0, 1e6], /overplot
      cgplot, [5896.52, 5896.52] , [0, 1e6], /overplot
      residual = y - scatter_fit
      residual_err = err_y

      O2_Io_frame                       = O2 + O2 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      ;Gaussian_parinfo[1].fixed        = 1
      ;Gaussian_parinfo[1].value        = O2_Io_frame
      Gaussian_parinfo[1].limited       = 1b
      Gaussian_parinfo[1].limits        = [O2_Io_frame-0.05, O2_Io_frame+0.05]
      ;Gaussian_parinfo[2].limited       = 1b
      ;Gaussian_parinfo[2].limits        = [0.07, 0.15]
      ;Gaussian_parinfo[2].fixed        = 1b
      ;Gaussian_parinfo[2].value        = 0.09

      ;if i eq 0 then hacked_errors = replicate(1000., n_elements(LSF_fitting_ind)) else $
      ;  hacked_errors = replicate(400., n_elements(LSF_fitting_ind))  ;hack
      LSF_Fit = mpfitpeak(WL[LSF_fitting_ind], residual[LSF_fitting_ind], a, PERROR = err_a, /POSITIVE, PARINFO = Gaussian_parinfo, $
                          STATUS = STATUS, MEASURE_ERRORS = residual_err[LSF_fitting_ind], nterms = 3, /quiet)

      WL_array[*, i]              = WL[include_WLs]
      residual_array[*, i]        = residual[include_WLs]
      LSF_Fit_array[*, i]         = gaussian(WL[include_WLs], a) ;LSF_Fit
      Brightness_array[i]         = $   ;  total(residual[LSF_fitting_ind])
        A[0]*0.09*SQRT(2*!DPI)                                                      ; fix the linewidth used in the fitting
      ;  A[0]*A[2]*SQRT(2*!DPI)
      err_Brightness_array[i]     = abs(Brightness_array[i] * sqrt( (err_A[0]/A[0])^2. + (err_A[2]*A[2])^2.) ) ; uncertainties in width and height add in quadrature.
      T_Shadow_array[i]           = sxpar(header, 'T_SHADOW')
      EXPTime_array[i]            = sxpar(raw_header, 'EXPTIME') / 60.
      DopplerShift_array[i]       = O2_Io_frame
      print, 'linewidth = ', a[2]
      stop
    endfor

    window, 2, xs = 800, ys = 1000
    cgplot, WL, spec, xr = [5888,5891], yr = [-1.e3, 7.e3], ytitle = 'Residual Jovian scattered light subtraction R/A', /nodata
    for i = 0, n_elements(Eclipse_files)-1 do begin
      cgplot, WL_array[*, i], residual_array[*, i], /OVERPLOT, COLOR = timeColors[i], psym=10
      cgplot, WL_array[*, i], LSF_Fit_array[*, i], /OVERPLOT, COLOR = timeColors[i]
      cgplot, [5890.05, 5890.05], [-10000., 10000], linestyle = 1, /overplot
      cgplot, [5890.5, 5890.5], [-10000., 10000], linestyle = 1, /overplot
      cgplot, [DopplerShift_array[i],DopplerShift_array[i]], [-10000., 10000], COLOR = timeColors[i], /overplot
    endfor
stop
    pos =  [.17,.2,.98,.9]
    cgPS_Open, filename = Reduced_Dir+'Na5888_lightcurve.eps', /ENCAPSULATED, xsize = 10., ysize = 6.5
    !P.font=1
    device, SET_FONT = 'Helvetica Bold', /TT_FONT
    !p.charsize = 3.
    loadct, 0
    cgplot, T_Shadow_array, Brightness_array, psym = 15, Ytitle = 'Rayleighs (5.1 sq arcsec aperture)', xtitle = 'Minutes After Ingress', title = 'Io''s NA 5888'+cgsymbol('Angstrom') + ' Response to Eclipse',  $
      ERR_YLOW = ERR_Brightness_array, ERR_YHigh = ERR_Brightness_array, ERR_XLOW = EXPTime_array/2., ERR_XHigh = EXPTime_array/2., yr = [0, 3500], xr = [70,150], /nodata, pos = pos

    x = findgen(Umbra_ET - PenUmbra_ET) / 60.
    colors = reverse(BytScl(x, MIN=min(x), MAX=2.*max(x)))+128
    ;reverse(BytScl(x, MIN=min(x), MAX=max(x)))
    FOR j=0,n_elements(x)-2 DO BEGIN
      xpoly = [x[j],         x[j], x[j+1],       x[j+1],         x[j]]
      ypoly = [  0., !Y.CRange[1], !Y.CRange[1], 0., 0.]
      cgColorFill, xpoly, ypoly, Color=colors[j]
    ENDFOR
    xpoly = [70,     70, !X.CRange[1],  !X.CRange[1],  70]
    ypoly = [  0., !Y.CRange[1], !Y.CRange[1], 0., 0.]
    cgColorFill, xpoly, ypoly, color = 'charcoal', /LINE_FILL, orientation = 45., thick = .1
    cgColorFill, xpoly, ypoly, color = 'charcoal', /LINE_FILL, orientation = -45., thick = .1
    cgtext, 2.5, !Y.CRange[1]/6., 'Penumbral Eclipse', orientation = 90., color = 'white'
    cgplot, T_Shadow_array, Brightness_array, psym = 4, /noerase, symsize = 2, yr = [0, 3500], xr = [70,150], pos = pos, $
      ERR_YLOW = ERR_Brightness_array, ERR_YHigh = ERR_Brightness_array, ERR_XLOW = EXPTime_array/2., ERR_XHigh = EXPTime_array/2.
    cgps_close











;    window, 1
;    for i = 1, 6 do begin
;      spec = MRDFITS( reduced_dir+'sumIo_free_and_clear.000'+strcompress(i, /remove_all)+'.ec.fits', 0, header, /fscale, /unsigned )
;      raw_header  = headfits(reduced_dir.Substring(0,-9)+'Io_free_and_clear.000'+strcompress(i, /remove_all)+'.fits')
;      spec = spec / sxpar(raw_header, 'EXPTIME')
;      MWRFITS, spec / Cont_Sensitivity, reduced_dir+'R_per_A_Io_free_and_clear.000'+strcompress(i, /remove_all)+'.ec.fits', header, /CREATE, /silent ;/create overwrites
;    endfor
;    window, 0
;    cgplot, WL, spec, xr = [5888,5898], yr = [0, 1.e6], ytitle = 'R/A', /nodata
;    for i = 1, 6 do begin
;      spec = MRDFITS( reduced_dir+'R_per_A_Io_free_and_clear.000'+strcompress(i, /remove_all)+'.ec.fits', 0, header, /fscale, /unsigned )
;      WL   = sxpar(header, 'CRVAL1')+findgen(N_elements(jup_center))*sxpar(header, 'Cdelt1')
;      cgplot, WL, spec-8.e2*i, xr = [5888,5898], Color=timeColors[i], yr = [0, 1.5e4], ytitle = 'R/A', /overplot
;    endfor
;    
   
    
    stop
    
;;    
;;    cgplot, WL_A, Rayleighs_per_angstrom-3.e5, color = 'blue', /overplot
;;    LSF = LSF_ROTATE(1., 20)
;;    cgplot, WL_A, convol(Rayleighs_per_angstrom, GAUSSIAN_FUNCTION(5), /normalize)-3.e5, color = 'orange', /overplot
;;    
;    window, 1
;    Telluric_Absorption_Spectrum = test_cal / smoothed_expected_flux
;    
;    cgplot, wl, test_cal, xr = [6290, 6310], psym =10
;    cgplot, wl, expected_flux, /overplot, color = 'red'
;stop
;
;
;    parinfo            = replicate( {value:0., fixed: 0b, limited: [0b,0b], limits: fltarr(2) }, 3)
;    parinfo[0].limited = 1b                                 ; limit differences in multiplier amplitude
;    parinfo[0].limits  = [0.5, 2.]
;    parinfo[1].limited = 1b                                 ; limit additive differences                 
;    parinfo[1].limits  = [-2.e6, 2.e6]     
;    parinfo[2].limited = 1b                                 ; limit differences in smooth width 
;    parinfo[2].limits  = [0.01, 5.]                              
;
;    fit_ind = where((wl gt 5880.) and (wl lt 5910.), /NULL)
;    
; 
;    X = expected_flux
;    Y = test_cal
;    Err_Y = sqrt(test) / Cont_Sensitivity
;    
;    dr_s         = minmax(Y[fit_ind], /NaN)
;    dr_r         = minmax(x[fit_ind], /NaN)
;    multiplier   = (dr_s[1]-dr_s[0]) / (dr_r[1]-dr_r[0])
;    X = x * multiplier
;    x = x + min(Y[fit_ind]) - min(x[fit_ind])
;    cgplot, wl, x, /overplot, color = 'blue'
;      
;    p0 = [1.11, -1.e6, 2.5] ; guess at initial Multiply P[0], Add P[1] and Smooth P[2] to 
;    p = mpfitfun('Match_Reference_Spectrum', x[fit_ind], y[fit_ind], err_y[fit_ind], p0, /NaN, status=status, parinfo = parinfo) ; Fit a y = mx + b function to the spectrum at every spatial bin... SLOW
;    
;    ; Multiply P[0], Add P[1] and Smooth P[2] a reference spectrum (X) until it best matches Y
;    hacked_expected_flux = P[0]*gauss_smooth(X, P[2], /EDGE_TRUNCATE) + P[1]
;    cgplot, wl, hacked_expected_flux, /overplot, color = 'green'
;    cgplot, wl, .95*gauss_smooth(X, 2.6, /EDGE_TRUNCATE) -10000., /overplot, color = 'orange'
;    airmass = sxpar(raw_header, 'Airmass')
stop

; ------------------------UT180320 Eclipse ---------------------------------
  dir  = 'D:\DATA\Apache Point\Echelle\Io Eclipses\UT180320\Reduced\'



spec = MRDFITS(dir+'fullspecIo_eclipsed.0001.ec.fits', 0, header, /fscale, /silent, /unsigned )
WL   = sxpar(header, 'CRVAL1')+findgen(170105)*sxpar(header, 'Cdelt1')

wl_range      = 8. ;angstroms
shift_array   = [0, 1.5, 1.75, 2., 2.25, 2.5, 2.75]
;Smooth_array  = [0,  18,  18,  18,  18,  18,  18] ; Good for H alpha 6563A
;Smooth_array  = [0,  14,  14,  14,  14,  14,  14] ; Good for Fe I    5270A
Smooth_array  = [0,  14,  14,  14,  14,  14,  14]  ; Good for 5893

window, 0, xs = 1000, ys = 800
pos = cgLayout([1,2], YGap=0)
jup = MRDFITS(dir+'fullspecJupiter_Disk_Center.0001.ec.fits', 0, header, /silent, /fscale, /unsigned )
junk = [min(abs(wl - (line-wl_range/2)), low_ind), min(abs(wl - (line+wl_range/2)), hi_ind)]
cgplot, wl, spec, xr = [Line-wl_range/2, line+wl_range/2], /nodata, pos = pos[*,0], yr = [-.6, 1.], xtickformat = '(A1)'
for i = 1, 6 do begin
  spec = MRDFITS(dir+'fullspecIo_eclipsed.000'+strcompress(i, /remove_all)+'.ec.fits', 0, header, /silent, /fscale, /unsigned )
  cgplot, wl, spec-.15*i, /overplot, Color=timeColors[i], thick = 2, pos = pos[*,0] 
  
  shifted = interpolate(smooth(jup, smooth_array[i]), findgen(n_elements(Jup)) - shift_array[i]) - .15*i
  cgplot, wl, shifted, /overplot, Color=timeColors[i], thick = 2, pos = pos[*,0]  ; sub-pixel shifting
  ;cgplot, wl, smart_shift(smooth(jup, smooth_array[i]), shift_array[i], /Interp)-.15*i, /overplot, Color=timeColors[i], thick = 2, pos = pos[*,0]  
endfor
cgplot, wl, spec, xr = [Line-wl_range/2, line+wl_range/2], /nodata, yr = [-.3, .5], pos = pos[*,1], /noerase 
cgplot, [5889.95095, 5889.95095], [-2., 2], linestyle = 1, /overplot
cgplot, [5895.92424, 5895.92424], [-2., 2], linestyle = 1, /overplot
cgplot, [O1, O1], [-2., 2], linestyle = 1, /overplot
cgplot, [O2, O2], [-2., 2], linestyle = 1, /overplot
cgplot, [O3, O3], [-2., 2], linestyle = 1, /overplot
cgplot, [O4, O4], [-2., 2], linestyle = 1, /overplot
cgplot, [O5, O5], [-2., 2], linestyle = 1, /overplot
cgplot, [O6, O6], [-2., 2], linestyle = 1, /overplot
cgplot, [O7, O7], [-2., 2], linestyle = 1, /overplot
cgplot, [O8, O8], [-2., 2], linestyle = 1, /overplot
cgplot, [O9, O9], [-2., 2], linestyle = 1, /overplot
cgplot, [S1, S1], [-2., 2], linestyle = 1, /overplot
cgplot, [S2, S2], [-2., 2], linestyle = 1, /overplot
cgplot, [S3, S3], [-2., 2], linestyle = 1, /overplot

tot_stat = 0
co_add = fltarr(N_elements(WL))
for i = 1, 6 do begin
  spec     = MRDFITS(dir+'fullspecIo_eclipsed.000'+strcompress(i, /remove_all)+'.ec.fits', 0, header, /silent, /fscale, /unsigned )
  ;residual = spec - smart_shift(smooth(jup, smooth_array[i]), shift_array[i], /Interp)
  
  residual = spec - interpolate(smooth(jup, smooth_array[i]), findgen(n_elements(Jup)) - shift_array[i])
  
  co_add = co_add+residual
  cgplot, wl, residual, /overplot, Color=timeColors[i], thick = 2, pos = pos[*,1]
  
  ; find and plot the intantaneous Earth-Io Doppler Shift
    cspice_UTC2ET, sxpar(header, 'DATE-OBS'), ET
    cspice_spkezr, 'Io', ET + float(sxpar(header, 'EXPTIME'))/2. , 'J2000', 'LT+S', 'Earth', Io_Earth_State, ltime 
    theta  = cspice_vsep(Io_Earth_state[0:2], Io_Earth_state[3:5])
    Io_wrt_Earth_Dopplershift = line * cos(theta) * norm(Io_Earth_State[3:5]) / cspice_clight()
    Io_wrt_Earth_Dopplershift = line * cos(theta) * sqrt(Io_Earth_State[3]^2.+Io_Earth_State[4]^2.+Io_Earth_State[5]^2.) / cspice_clight()
    cgplot, [5889.95095, 5889.95095]+Io_wrt_Earth_Dopplershift, [-2., 2], Color=timeColors[i], linestyle = 2, /overplot
    cgplot, [5895.92424, 5895.92424]+Io_wrt_Earth_Dopplershift, [-2., 2], Color=timeColors[i], linestyle = 2, /overplot
    
    cgplot, [O1, O1]+Io_wrt_Earth_Dopplershift, [-2., 2], Color=timeColors[i], linestyle = 2, /overplot
    cgplot, [O2, O2]+Io_wrt_Earth_Dopplershift, [-2., 2], Color=timeColors[i], linestyle = 2, /overplot
    cgplot, [O3, O3]+Io_wrt_Earth_Dopplershift, [-2., 2], Color=timeColors[i], linestyle = 2, /overplot
    cgplot, [O4, O4]+Io_wrt_Earth_Dopplershift, [-2., 2], Color=timeColors[i], linestyle = 2, /overplot
    cgplot, [O5, O5]+Io_wrt_Earth_Dopplershift, [-2., 2], Color=timeColors[i], linestyle = 2, /overplot
    cgplot, [O6, O6]+Io_wrt_Earth_Dopplershift, [-2., 2], Color=timeColors[i], linestyle = 2, /overplot
    cgplot, [O7, O7]+Io_wrt_Earth_Dopplershift, [-2., 2], Color=timeColors[i], linestyle = 2, /overplot
    cgplot, [O8, O8]+Io_wrt_Earth_Dopplershift, [-2., 2], Color=timeColors[i], linestyle = 2, /overplot
    cgplot, [O9, O9]+Io_wrt_Earth_Dopplershift, [-2., 2], Color=timeColors[i], linestyle = 2, /overplot
    cgplot, [S1, S1]+Io_wrt_Earth_Dopplershift, [-2., 2], Color=timeColors[i], linestyle = 2, /overplot
    cgplot, [S2, S2]+Io_wrt_Earth_Dopplershift, [-2., 2], Color=timeColors[i], linestyle = 2, /overplot
    cgplot, [S3, S3]+Io_wrt_Earth_Dopplershift, [-2., 2], Color=timeColors[i], linestyle = 2, /overplot
    
  stat = total(abs(residual[low_ind:hi_ind]))
  tot_stat = tot_stat + stat
  print, i, shift_array[i], smooth_array[i], stat
endfor
stop




cgplot, wl, co_add, thick = 2, /overplot

print, tot_stat



window, 1
Line = 5577.330
pos = cgLayout([1,2], YGap=0)
cgplot, wl, spec, xr = [line-1.5, line+1.5], /nodata, yr = [.9, 1.25], pos = pos[*,0]
for i = 1, 6 do begin
  spec = MRDFITS(dir+'fullspecIo_eclipsed.000'+strcompress(i, /remove_all)+'.ec.fits', 0, header, /fscale, /unsigned )
  cgplot, wl, spec, /overplot, Color=timeColors[i], thick = 2
endfor
cgplot, wl, spec, xr = [line-1.5, line+1.5], /nodata, yr = [.9, 1.25], pos = pos[*,1], /noerase
for i = 1, 3 do begin
  spec = MRDFITS(dir+'fullspecJovian_Scatter.000'+strcompress(i, /remove_all)+'.ec.fits', 0, header, /fscale, /unsigned )
  cgplot, wl, spec, /overplot, thick = i
endfor

window, 2
Line = 6300.304
cgplot, wl, spec, xr = [line-1.5, line+1.5], /nodata, yr = [.75, 1.4], pos = pos[*,0]
for i = 1, 6 do begin
  spec = MRDFITS(dir+'fullspecIo_eclipsed.000'+strcompress(i, /remove_all)+'.ec.fits', 0, header, /fscale, /unsigned )
  cgplot, wl, spec, /overplot, Color=timeColors[i], thick = 2
endfor
cgplot, wl, spec, xr = [line-1.5, line+1.5], /nodata, yr = [.8, 1.3], pos = pos[*,1], /noerase
for i = 1, 3 do begin
  spec = MRDFITS(dir+'fullspecJovian_Scatter.000'+strcompress(i, /remove_all)+'.ec.fits', 0, header, /fscale, /unsigned )
  cgplot, wl, spec, /overplot, thick = i
endfor

window, 3
Line = 6363.776
cgplot, wl, spec, xr = [line-1.5, line+1.5], /nodata, yr = [.9, 1.1], pos = pos[*,0]
for i = 1, 6 do begin
  spec = MRDFITS(dir+'fullspecIo_eclipsed.000'+strcompress(i, /remove_all)+'.ec.fits', 0, header, /fscale, /unsigned )
  cgplot, wl, spec, /overplot, Color=timeColors[i], thick = 2
endfor
cgplot, wl, spec, xr = [line-1.5, line+1.5], /nodata, yr = [.9, 1.1], pos = pos[*,1], /noerase
for i = 1, 3 do begin
  spec = MRDFITS(dir+'fullspecJovian_Scatter.000'+strcompress(i, /remove_all)+'.ec.fits', 0, header, /fscale, /unsigned )
  cgplot, wl, spec, /overplot, thick = i, pos = pos[*,1]
endfor


window, 3
Line = 6363.776
cgplot, wl, spec, xr = [line-1.5, line+1.5], /nodata, yr = [.9, 1.1], pos = pos[*,0]
for i = 1, 6 do begin
  spec = MRDFITS(dir+'fullspecIo_eclipsed.000'+strcompress(i, /remove_all)+'.ec.fits', 0, header, /fscale, /unsigned )
  cgplot, wl, spec, /overplot, Color=timeColors[i], thick = 2
endfor
cgplot, wl, spec, xr = [line-1.5, line+1.5], /nodata, yr = [.9, 1.1], pos = pos[*,1], /noerase
for i = 1, 3 do begin
  spec = MRDFITS(dir+'fullspecJovian_Scatter.000'+strcompress(i, /remove_all)+'.ec.fits', 0, header, /fscale, /unsigned )
  cgplot, wl, spec, /overplot, thick = i, pos = pos[*,1]
endfor

window, 4
Line = 7664.8969 ; K D2
;Line = 7698.9600  ; K D1
;Line = 7682.9600 ; K both
;Line = 8375.94   ; Cl I   
;Line = 4227.     ; Ca I
;Line = 8446.38   ; O I
;Line = 6731.      ; S II
;Line = 6716.      ; S II
jup1 = MRDFITS(dir+'fullspecJovian_Scatter.0001.ec.fits', 0, header, /fscale, /unsigned )
jup2 = MRDFITS(dir+'fullspecJovian_Scatter.0002.ec.fits', 0, header, /fscale, /unsigned )
jup3 = MRDFITS(dir+'fullspecJovian_Scatter.0003.ec.fits', 0, header, /fscale, /unsigned )
jup =mean([[jup1],[jup2],[jup3]], dimension=2)

cgplot, wl, spec, xr = [line-1.5, line+1.5], /nodata, yr = [-.2, .2], pos = pos[*,0]
for i = 1, 6 do begin
  spec = MRDFITS(dir+'fullspecIo_eclipsed.000'+strcompress(i, /remove_all)+'.ec.fits', 0, header, /fscale, /unsigned )
  ;if i le 3 then jup = MRDFITS(dir+'fullspecJovian_Scatter.000'+strcompress(i, /remove_all)+'.ec.fits', 0, header, /fscale, /unsigned )
  cgplot, wl, spec-jup, /overplot, Color=timeColors[i], thick = 2
endfor
cgplot, wl, spec, xr = [line-1.5, line+1.5], /nodata, yr = [.7, 1.1], pos = pos[*,1], /noerase
for i = 1, 3 do begin
  spec = MRDFITS(dir+'fullspecJovian_Scatter.000'+strcompress(i, /remove_all)+'.ec.fits', 0, header, /fscale, /unsigned )
  cgplot, wl, spec, /overplot, thick = i, pos = pos[*,1]
endfor

stop





end