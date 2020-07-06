; Analysis of Io's optical spectra in eclipse from datasets with Apache Point and Large Binocular telescopes
; Written by C. Schmidt, M. Sharov. BU Center for Space Physics, 2019

FUNCTION Match_Reference_Spectrum, X, P
  ; Multiply P[0], Add P[1] and Smooth P[2] a reference spectrum (X) until it best matches Y
  return, P[0]*gauss_smooth(X, P[2], /EDGE_TRUNCATE) + P[1]
  ;return, P[0]*X + P[1]
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

Pro Io_Eclipses_ARCES_PEPSI_V3, Part=Part, Date=Date


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
      Penumbra_UTC    = '2019-Aug-12 03:44:40'
      Umbra_UTC       = '2019-Aug-12 03:48:19'
      dir             = 'D:\DATA\Apache Point\Echelle\Io Eclipses\'+date+'\'
      reduced_dir     = 'D:\DATA\Apache Point\Echelle\Io Eclipses\'+date+'\Reduced\'
      calibration_dir = 'D:\DATA\Apache Point\Echelle\Io Eclipses\'+date+'\Reduced\'  ; Directory with Stellar Spectra for Telluric Corrections
      Eclipse_files   = ['Io_eclipsed.000'+strcompress(indgen(7)+3, /remove_all), 'Io_eclipsed.00'+strcompress(indgen(2)+10, /remove_all) ]
      Jupiter_Center_File = 'Jupiter_Disk_Center.0002'
      ;Jupiter_Center_File = 'Jupiter_Disk_Center.0003'
    end
    'UT180320': begin
      Penumbra_UTC    = '2018-Mar-20 10:46:46'
      Umbra_UTC       = '2018-Mar-20 10:50:22'
      dir             = 'D:\DATA\Apache Point\Echelle\Io Eclipses\'+date+'\'
      reduced_dir     = 'D:\DATA\Apache Point\Echelle\Io Eclipses\'+date+'\Reduced\'
      calibration_dir = 'D:\DATA\Apache Point\Echelle\Io Eclipses\UT190812\Reduced\'  ; Directory with Stellar Spectra for Telluric Corrections
      ;Eclipse_files   = [ 'Io_penumbra.0001' , 'Io_eclipsed.000'+strcompress(indgen(8)+2, /remove_all), 'Io_eclipsed.00'+strcompress(indgen(2)+10, /remove_all) ]
      Eclipse_files = [ 'Io_penumbra.0001' , 'Io_eclipsed.000'+strcompress(indgen(6)+1, /remove_all) ]
      Jupiter_Center_File = 'Jupiter_Center_Spectrum.0001'
    end
    'UT190424': begin
      dir             = 'D:\DATA\LBT\'
      reduced_dir     = 'D:\DATA\LBT\Reduced\'
      calibration_dir = 'D:\DATA\LBT\Reduced\'
    end
  endcase

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

  ; Define ingrees / egress times
  cspice_UTC2ET, PenUmbra_UTC, PenUmbra_ET
  cspice_UTC2ET, Umbra_UTC, Umbra_ET

  ;================================Part 0 Generate of Telluric Absorption Spectrum==================================================================================
  if part eq 0 then begin

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

    READCOL,'D:\DATA\___Calibration___\Solar_and_Stellar_Calibration_Spectra\vegallpr25.30000', F='A,A', Model_A0V_WL, Model_A0V_Flux, Skipline = 700000, numline = 500000, /Silent
    Model_A0V_WL = Model_A0V_WL*10.
    Model_A0V_Flux = float(shift(Model_A0V_Flux, -100))
    Model_A0V_Flux = INTERPOL(Model_A0V_Flux, Model_A0V_WL, Telluric_WL_A0V)

    result = ROBUST_POLY_FIT(Telluric_WL_A0V , Telluric_Absorption_A0V, 9)
    Echelle_Order_Interweave = Telluric_Absorption_A0V / poly(Telluric_WL_A0V, result)
    result = ROBUST_POLY_FIT(Telluric_WL_A0V, Model_A0V_flux, 6)
    Model_A0V_norm = Model_A0V_flux / poly(Telluric_WL_A0V, result)

    A0V_norm = MRDFITS(Calibration_Dir+'\FullspecHD_155379_A0V.0001.ec.fits', 0, header, /fscale, /silent, /unsigned )
    Telluric_Absorption_A0V = A0V_norm / (Model_A0V_norm)

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
    ;Jupiter_center_header = headfits(reduced_dir+'sumfitcontJupiter_Disk_Center.0001.ec.fits')
    Jupiter_center_header = headfits(reduced_dir+'sumfitcont'+Jupiter_Center_File+'.ec.fits')
    cspice_UTC2ET, sxpar(Jupiter_center_header, 'DATE-OBS'), ET
    cspice_spkezr, 'Jupiter', ET, 'J2000', 'LT+S', 'Sun', Jupiter_Sun_State, ltime
    cspice_spkezr, 'Jupiter', ET, 'J2000', 'LT+S', 'Earth', Jupiter_Earth_State, ltime
    solar_distance = norm(Jupiter_Sun_State[0:2]) / 149597871.
    flux_at_jupiter = flux / solar_distance^2.

    ; Multiply incident solar irradiance x spectral albedo to determine the theoreitcal brightness of Jupiter at disk center
    Albedo = INTERPOL(Karkoschka_Albedo, Karkoschka_WL, WL_A)
    Rayleighs_per_angstrom = 4.*flux_at_jupiter*albedo / 1.e6
    window, 0, Title = 'Instantaneous Rayleighs per Angstrom: Center of Jupiter''s Disk'
    plot, WL_A, Rayleighs_per_angstrom, xr = [5885, 5900], charsize = 2 ;compare to 5.5 MR per angstrom by Brown & Schneider, 1981

    ; Adjust the expected absolute flux for Jupiter's Instantaneous Doppler shift
    theta  = cspice_vsep(Jupiter_Earth_State[0:2], Jupiter_Earth_State[3:5])
    Dopplershift = cos(theta) * norm(Jupiter_Earth_State[3:5])   ; scalar projection of the relative velocity along the line of sight
    WL_A = WL_A + WL_A*Dopplershift/cspice_clight()

    ; ***********Correct the Unity-Normalized Jupiter spectrum for Telluric Absorption**********
    ; Try both the BFR and the A0V star methods.
    jup_center           = MRDFITS(reduced_dir+'fullspec'+Jupiter_Center_File+'.ec.fits', 0, header, /fscale, /silent, /unsigned )     ; fullspec, normalized to one
    jup_center_err       = MRDFITS(reduced_dir+'sig_fullspec'+Jupiter_Center_File+'.ec.fits', 0, header, /fscale, /silent, /unsigned ) ; 1 sig uncertainty in fullspec
    WL                   = sxpar(header, 'CRVAL1')+findgen(N_elements(jup_center))*sxpar(header, 'Cdelt1')
    rough_absorption_A0V = interpol(Telluric_Absorption_A0V, Telluric_WL_A0V, WL) ;using A0V star correction instead of BFR
    rough_absorption_BFR = interpol(Telluric_Absorption_BFR, Telluric_WL_BFR, WL) ;using BFR star correction instead of AOV

    ; Fix any subtle pixel shifts between *known telluric lines* measured in the calibration stars and in Jupiter
    xr = [6292,6308]  ; Wavelength region for inspection
    telluric_fitting_ind = where( ((wl gt 6295.0) and (wl lt 6296.5)) or ((wl gt 6305.5) and (wl lt 6307.0)), /NULL)                 ; *known easily isolated telluric lines
    lag        = findgen(21)-10.
    correl_BFR = C_CORRELATE( Rough_Absorption_BFR[telluric_fitting_ind], jup_center[telluric_fitting_ind], lag)
    correl_A0V = C_CORRELATE( Rough_Absorption_A0V[telluric_fitting_ind], jup_center[telluric_fitting_ind], lag)
    yfit_BFR   = MPFITPEAK(lag, correl_BFR, A_BFR, nterms = 3, /positive)
    yfit_A0V   = MPFITPEAK(lag, correl_A0V, A_A0V, nterms = 3, /positive)
    aligned_absorp_BFR = interpolate(Rough_Absorption_BFR, findgen(n_elements(Rough_Absorption_BFR)) - A_BFR[1])
    aligned_absorp_A0V = interpolate(Rough_Absorption_A0V, findgen(n_elements(Rough_Absorption_A0V)) - A_A0V[1])

    ; Inspect the alignment of both methods of telluric absorption and Jupiter
    window, 0, title = 'Checking Wavelength Alignment of the Telluric Absorption: Bottom plot looks aligned after the above shift? ', ys = 800
    pos = cglayout([1,2])
    cgplot, lag, correl_BFR, pos = pos[*,0]
    cgplot, lag, yfit_BFR, /overplot, color = 'blue'
    cgplot, lag, correl_A0V, /overplot
    cgplot, lag, yfit_A0V, /overplot, color = 'orange'
    cgplot, WL, jup_center, xr = xr, /ynozero, pos = pos[*,1], /noerase
    cgplot, WL, aligned_absorp_BFR, color = 'blue', /overplot
    cgplot, WL, aligned_absorp_A0V, color = 'orange', /overplot

    print, 'USER---> Carefully inspect the Telluric Absorption Correction at Jupiter in every Bandpass:', XR

    ; Correct the interweave errors in Jupiter and get a better unity normalization of the continuum
    continuum_ind = []
    step = 25
    for i = 0, N_elements(jup_center) - (N_elements(jup_center) mod step) -1, step do begin
      ind = reverse(sort(jup_center[i:i+step]))
      continuum_ind = [continuum_ind, i+ind[1:5]]
    endfor
    continuum_ind    = continuum_ind[UNIQ(continuum_ind, SORT(continuum_ind))]
    Correct_IW       = interpol(1./jup_center[continuum_ind], WL[continuum_ind], WL) ; This is the few % correction factor to force the continuum to exactly 1. It is applied to *all* science data
    jup_center_N     = jup_center*correct_IW                                         ; Correct for interweave
    jup_center_err_N = jup_center_err*correct_IW                                     ; Correct for interweave

    Window, 0, Title = 'Jupiter Center: Black = FullSpec, Grey = Continuum Indices, Blue = Echelle Order Interweave Normalized', xpos = -1000, ypos=100
    cgplot, WL, jup_center, xr = xr, /ynozero
    cgplot, WL[continuum_ind], jup_center[continuum_ind], color = 'grey', /overplot
    cgplot, WL, jup_center_N, color = 'blue', /overplot

    ; And now do the same echelle order interweave normalization for the telluric (BFR and A0V) spectrum
    continuum_ind = []
    step = 25
    for i = 0, N_elements(aligned_absorp_BFR) - (N_elements(aligned_absorp_BFR) mod step) -1, step do begin
      ind = reverse(sort(aligned_absorp_BFR[i:i+step]))
      continuum_ind = [continuum_ind, i+ind[1:5]]
    endfor
    continuum_ind    = continuum_ind[UNIQ(continuum_ind, SORT(continuum_ind))]
    correct_IW_star  = interpol(1./aligned_absorp_BFR[continuum_ind], WL[continuum_ind], WL)
    aligned_absorp_BFR_n = aligned_absorp_BFR * correct_IW_star

    continuum_ind = []
    step = 25
    for i = 0, N_elements(aligned_absorp_A0V) - (N_elements(aligned_absorp_A0V) mod step) -1, step do begin
      ind = reverse(sort(aligned_absorp_A0V[i:i+step]))
      continuum_ind = [continuum_ind, i+ind[1:5]]
    endfor
    continuum_ind    = continuum_ind[UNIQ(continuum_ind, SORT(continuum_ind))]
    correct_IW_star  = interpol(1./aligned_absorp_A0V[continuum_ind], WL[continuum_ind], WL)
    aligned_absorp_A0V_n = aligned_absorp_A0V * correct_IW_star

    ; Set the region of the interstellar absorption equal to zero
    na_absorb_ind = where( ((wl gt 5889.75) and (wl lt 5890.7)) or ((wl gt 5895.5) and (wl lt 5896.65)), /NULL)
    ;aligned_absorp_A0V_n[na_absorb_ind] = 1.0
    ;aligned_absorp_BFR_n[na_absorb_ind] = 1.0

    Window, 1, Title = 'Blue Fast Rotator: Black = FullSpec, Grey = Continuum Indices, Blue = Echelle Order Interweave Normalized', xpos = -1000, ypos=700
    cgplot, WL, aligned_absorp_BFR, xr = xr, /ynozero
    cgplot, WL[continuum_ind], aligned_absorp_BFR[continuum_ind], color = 'grey', /overplot
    cgplot, WL, aligned_absorp_BFR_n, color = 'blue', /overplot

    ; Match the depth of the telluric absorption to the Jupiter Disk-Center measurement.
    ;        telluric_fitting_ind = where( ((wl gt 6276.4) and (wl lt 6280.)) or ((wl gt 6281.3) and (wl lt 6282.7)) or $
    ;                                      ((wl gt 6287.0) and (wl lt 6291.)) or ((wl gt 6292.0) and (wl lt 6293.0)), /NULL)

    parinfo            = replicate( {value:0., fixed: 0b, limited: [0b,0b], limits: fltarr(2) }, 3)
    parinfo[0].limited = 1b                               ; limit differences in multiplier amplitude
    parinfo[0].limits  = [0.2, 5.]                        ; multiplicative range limits
    parinfo[1].limited = 1b                               ; limit additive differences
    parinfo[1].limits  = [-0.9, 0.9]                      ; additive adjustments should be small
    parinfo[2].Fixed   = 1b                               ; Fix spectral smoothing
    parinfo[2].value   = 0.                               ; No spectral smoothing
    p0 = [0.8, 0.2, 0.] ; guess at initial Multiply P[0], Add P[1] and Smooth P[2] to

    ; Fit a y = mx + b function to the spectrum, where x is the tellutic absorbtion scaled for 0 to 1 dynamic range, and y is jupiter, scaled for 0 to 1 dyanmic range
    p = mpfitfun('Match_Reference_Spectrum', aligned_absorp_BFR_n[telluric_fitting_ind], jup_center_n[telluric_fitting_ind], 0.1*jup_center_err_n[telluric_fitting_ind], $
      p0, /NaN, status=status, parinfo = parinfo, /quiet)
    telluric_fit_BFR = P[0]*aligned_absorp_BFR_n + P[1]

    p = mpfitfun('Match_Reference_Spectrum', aligned_absorp_A0V_n[telluric_fitting_ind], jup_center_n[telluric_fitting_ind], 0.1*jup_center_err_n[telluric_fitting_ind], $
      p0, /NaN, status=status, parinfo = parinfo, /quiet) ; Fit a y = mx + b function to the spectrum
    telluric_fit_A0V = P[0]*aligned_absorp_A0V_n + P[1]

    window, 2, xs = 1000, ys = 900, title='Black = Jupiter, Green = Indices where Telluric Absorptions are scaled, Red = Aligned Telluric Absorption, Blue = Scaled absorption to Match Jupiter'
    pos = cgLayout([1,2], OXMargin=[11,3], OYMargin=[9,6], YGap=0)
    cgplot, WL, jup_center_n, xr = xr, /ynozero, pos = pos[*,0], xtickformat = '(A1)'
    cgplot, WL[telluric_fitting_ind], jup_center_n[telluric_fitting_ind], /overplot, color = 'green', psym=14
    cgplot, WL, aligned_absorp_A0V_n, color = 'red', /overplot
    cgplot, WL, telluric_fit_A0V, color = 'blue', /overplot
    cgplot, WL, jup_center_n/telluric_fit_BFR, xr = xr, /ynozero, /noerase, pos = pos[*,1], color = 'blue'
    cgplot, WL, jup_center_n/telluric_fit_A0V, /overplot, color = 'orange'
    cgtext, xr[0]+1, 0.90, 'Correction using BFR', color = 'blue'
    cgtext, xr[0]+1, 0.85, 'Correction using A0V', color = 'orange'

    ; Keep whichever telluric correction method is closer to unity at every spectral bin
    residual_BFR  = abs(1.-jup_center_n/telluric_fit_BFR)
    residual_A0V  = abs(1.-jup_center_n/telluric_fit_A0V)
    junk          = min([[residual_BFR],[residual_A0V]], dim = 2, loc)
    combined      = [[jup_center_n/telluric_fit_BFR], [jup_center_n/telluric_fit_A0V]]
    tell_combined = [[telluric_fit_BFR], [telluric_fit_A0V]]
    telluric_fit  = tell_combined[loc]
    jup_tell_corr = combined[loc]
    cgplot, WL, jup_tell_corr, /overplot, color = 'red'
    cgtext, xr[0]+1, 0.80, 'Correction using a combination', color = 'red'

    ; We've been working with unity-normalized spectra until now. Determine the DN / S from the "Sum" files, which are real DN, not unity-normalized for continuum
    jup_sum        = MRDFITS(reduced_dir+'sum'+Jupiter_Center_File+'.ec.fits', 0, header, /fscale, /silent, /unsigned )
    raw_header     = headfits(dir+''+Jupiter_Center_File+'.fits')
    jup_sum        = jup_sum / sxpar(raw_header, 'EXPTIME')    ; Convert to DN per second
    jup_sum        = jup_sum / telluric_fit                    ; Remove telluric absorption
    jup_Sum_Curve  = smooth(jup_Sum, 2500, /edge_truncate, /nan)

    ; Determine the instrumental sensitivity from the expected versus measured flux at Jupiter Disk Center. This should be a smooth function
    expected_flux = interpol(Rayleighs_per_angstrom, WL_A, WL)   ; move expected flux to the Jovian Doppler-shift, UNITS are R / A
    smoothed_expected_flux = GAUSS_SMOOTH(expected_flux, 1.9)    ; this smoothing looks about right for APO/ARCES, UNITS are R / A
    window, 3
    Sensitivity       = jup_sum / smoothed_expected_flux         ; Sensitivity in (DN / S) / (R / A)
    Sensitivity_Curve = smooth(sensitivity, 5000, /edge_truncate, /nan)
    cgplot, WL, Sensitivity, xr = [3700, 10000], yr = [0.,0.0025], Ytitle = 'Measured Flux Sensitivity (DN/S) / (R/A)', Xtitle = 'Angstroms'
    fit_Sensitivity   = Sensitivity
    fit_sensitivity[where( (sensitivity gt 0.002) or (sensitivity lt 0.0), /Null)] = !values.F_NaN                                            ; some basic rejection criterion
    Sensitivity_Curve = smooth(fit_sensitivity, 5000, /edge_truncate, /nan)
    fit_sensitivity[where( (Sensitivity_Curve/fit_Sensitivity gt 1.5) or (Sensitivity_Curve/fit_Sensitivity lt 0.5), /Null)] = !values.F_NaN  ; further rejection criterion
    Sensitivity_Curve = smooth(fit_sensitivity, 5000, /edge_truncate, /nan)
    cgplot, WL, Sensitivity_Curve, color = 'red', /overplot

    ; Inspect and write the telluric corrected R/A Calibrated Jupiter disk center spectra to file
    window, 4, Title = 'Inspect Jupiter: Red = Uncorrected, Black = Telluric Corrected'
    Jup_Cal_Tell_Corr       = jup_tell_corr * jup_Sum_Curve / (Sensitivity_Curve)  ; Telluric Correct and Convert to Rayleighs per Angstrom
    Jup_Cal_No_Tell_Corr    = jup_center_N * jup_Sum_Curve / (Sensitivity_Curve)   ; Convert to Rayleighs per Angstrom, no correction
    cgplot, WL, Jup_Cal_Tell_Corr, xr = xr, /ynozero
    cgplot, WL, Jup_Cal_No_Tell_Corr, color = 'red', /overplot
    MWRFITS, Jup_Cal_Tell_Corr, reduced_dir+'R_per_A_Jupiter_Center.ec.fits', header, /CREATE ;/create overwrites
    ;MWRFITS, Jup_Cal_Tell_Corr_err, reduced_dir+'sig_R_per_A_Jupiter_Center.ec.fits', sig_header, /CREATE, /silent ;/create overwrites for sig values (not fancy mate :( will be fixed later)

    ; **************************** Finally... onto Io Calibration and Absorption Correction **********************************************
    window, 5, title = 'Io Eclipse Spectrum: Black = Interweave normalized, Blue = Telluric Abs, Red = Telluric corrected', xs = 1900
    cspice_UTC2ET, PenUmbra_UTC, PenUmbra_ET
    cspice_UTC2ET, Umbra_UTC, Umbra_ET

    ; Correct tellurics and write the spectra into R/A units...
    for i = 0, n_elements(Eclipse_files)-1 do begin

      Io_ecl     = MRDFITS(reduced_dir+'fullspec'+Eclipse_files[i]+'.ec.fits', 0, full_header, /fscale, /unsigned, /silent)
      Io_ecl_err = MRDFITS(reduced_dir+'sig_fullspec'+Eclipse_files[i]+'.ec.fits', 0, sig_header, /fscale, /unsigned, /silent )
      Io_WL      = sxpar(header, 'CRVAL1')+findgen(N_elements(Io_ecl))*sxpar(header, 'Cdelt1')

      ; Correct the echell order interweave errors in Io. Assume that Io's interweave correcction is the same as measured at Jupiter
      Io_ecl_N         = Io_ecl*correct_IW                                         ; Correct for interweave
      Io_ecl_err_N     = Io_ecl_err*correct_IW                                     ; Correct for interweave

      correl     = C_CORRELATE( telluric_fit[telluric_fitting_ind], Io_ecl_N[telluric_fitting_ind], lag )
      yfit_BFR   = MPFITPEAK(lag, correl, A, nterms = 3, /positive)
      aligned_telluric_fit = interpolate(telluric_fit, findgen(n_elements(telluric_fit)) - A[1])

      p = mpfitfun('Match_Reference_Spectrum', aligned_telluric_fit[telluric_fitting_ind], Io_Ecl[telluric_fitting_ind], Io_ecl_err[telluric_fitting_ind], $
        p0, /NaN, status=status, parinfo = parinfo, /quiet) ; Fit a y = mx + b function to the spectrum
      telluric_fit_for_Io  = P[0]*aligned_telluric_fit + P[1]

      ; Inspect
      window, 2, xs = 1000, ys = 900, title='Black = Io, Green = Indices where Telluric Absorptions are scaled, Red = Aligned Telluric Absorption, Blue = Scaled absorption to Match Io'
      pos = cgLayout([1,2], OXMargin=[11,3], OYMargin=[9,6], YGap=0)
      cgplot, WL, Io_ecl_n, xr = xr, /ynozero, pos = pos[*,0], xtickformat = '(A1)'
      cgplot, WL[telluric_fitting_ind], Io_ecl_n[telluric_fitting_ind], /overplot, color = 'green', psym=14
      cgplot, WL, aligned_absorp_A0V_n, color = 'red', /overplot
      cgplot, WL, telluric_fit_A0V, color = 'blue', /overplot
      cgplot, WL, Io_ecl_n/telluric_fit_BFR, xr = xr, /ynozero, /noerase, pos = pos[*,1], color = 'blue'
      cgplot, WL, Io_ecl_n/telluric_fit_A0V, /overplot, color = 'orange'
      cgplot, WL, Io_ecl_n/telluric_fit_for_Io, /overplot, color = 'red'
      cgtext, xr[0]+1, 0.90, 'Correction using BFR', color = 'blue'
      cgtext, xr[0]+1, 0.85, 'Correction using A0V', color = 'orange'

      Io_ecl_tell_corr = Io_ecl_n / telluric_fit_for_Io

      Io_ecl_sum        = MRDFITS(reduced_dir+'sum'+Eclipse_files[i]+'.ec.fits', 0, header, /fscale, /silent, /unsigned )
      raw_header        = headfits(dir+Eclipse_files[i]+'.fits')
      Io_ecl_sum        = Io_ecl_sum / sxpar(raw_header, 'EXPTIME')      ; Convert to DN per second
      Io_ecl_sum        = Io_ecl_sum / telluric_fit_for_Io
      Io_ecl_Sum_Curve  = smooth(Io_ecl_Sum, 2500, /edge_truncate, /nan)

      Io_ecl_Cal_Tell_Corr       = Io_ecl_tell_corr * Io_ecl_Sum_Curve / (Sensitivity_Curve)  ; Convert to Rayleighs per Angstrom, Telluric Corrected
      Io_ecl_Cal_No_Tell_Corr    = Io_ecl_N * Io_ecl_Sum_Curve / (Sensitivity_Curve)          ; Convert to Rayleighs per Angstrom, no correction
      ;window, 0
      ;cgplot, WL, Io_ecl_Cal_Tell_Corr, xr = xr, /ynozero
      ;cgplot, WL, Io_ecl_Cal_No_Tell_Corr, color = 'red', /overplot

      ;Penumbral emissions also need Io reflactance subtracted
      if date eq 'UT180320' and i eq 0 then begin
        window, 0
        cgplot, WL, Io_ecl_Cal_Tell_Corr, xr = xr, /ynozero
        Io_sunlit = MRDFITS(reduced_dir+'fullspecIo_free_and_clear.0006.ec.fits', 0, header, /fscale, /silent, /unsigned )     ; fullspec, normalized to one
        Io_sunlit = Io_sunlit * correct_IW
        Io_sunlit = Io_sunlit / telluric_fit_for_Io
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
      SXADDPAR, Header, 'Io_DOPPL', Io_wrt_Earth_Dopplershift, 'Io-Earth V_radial in km/s (mid exposure)'
      SXADDPAR, Header, 'T_SHADOW', (ET_mid_exposure-PenUmbra_ET) / 60., 'Minutes since Penumbral ingress'
      
      MWRFITS, Io_ecl_Cal_Tell_Corr, reduced_dir+'R_per_A_'+Eclipse_files[i]+'.ec.fits', header, /CREATE ;/create overwrites
      
    endfor
  endif ; Part eq 0

  if Part eq 1 then begin ; Make waterfall plots

    ; get color versus ingress time
    timeColors = BytScl(findgen(11), Min=0, Max=11, Top=11)
    Color=timeColors[0]
    cgLoadCT, 33, NColors=8, /reverse
    
    ;O6300 demands telluric corrected spectra:
    
    ;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ;--------------------------------------------------------Get the O 6300 Time Series and Scatter Fitting-------------------------------------------------------------------------
    ;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    IF (date EQ 'UT180320') THEN BEGIN
      range1 = 0
      range2 = 40
    ENDIF ELSE BEGIN
      range1 = 70
      range2 = 135
    ENDELSE

    parinfo            = replicate( {value:0., fixed: 0b, limited: [0b,0b], limits: fltarr(2) }, 4)
    parinfo[0].limited = 1b                                 ; limit differences in multiplier amplitude
    parinfo[0].limits  = [1.e-4, 1.e-2]
    parinfo[1].limited = 1b                                 ; limit additive differences
    parinfo[1].limits  = [-1.e3, 1.e4]
    parinfo[2].limited = 1b                                 ; limit differences in smooth width
    parinfo[2].limits  = [1., 8.]
    parinfo[3].limited = 1b                                 ; limit shift alignment
    parinfo[3].limits  = [-10., 10.]

    ;window, 0, xs = 800, ys = 1000
    ;cgplot, WL, WL, xr = [6296,6302], yr = [0.9e4, 8.e4], ytitle = 'R/A', /nodata
    ;cgplot, WL, spec, xr = [5888,5898], yr = [1.e3, 1.e5], ytitle = 'R/A', /nodata

    Jup_Cal_Tell_Corr = MRDFITS(reduced_dir+'R_per_A_Jupiter_Center.ec.fits', 0, header, /fscale, /unsigned, /silent)
    Jup_Cal_Tell_Corr_err = MRDFITS(reduced_dir+'sig_R_per_A_Jupiter_Center.ec.fits', 0, header, /fscale, /unsigned, /silent)

    ;-------------------Waterfall Plot Postscript Io and Fit Jovian scatter--------------------------------------------------------------------------------------------------------------
    cgLoadCT, 33, NColors=8, /reverse
    pos = cgLayout([1,2], OXMargin=[44,11], OYMargin=[10,8], YGap=0)
    Pos[1,0] = Pos[1,0]*.7 & Pos[3,1] = Pos[3,1]*.7
    cgPS_Open, filename = Reduced_Dir+'O6300_scatterfit_V1.eps', /ENCAPSULATED, xsize = 10, ysize = 10.
    !P.font=1
    device, SET_FONT = 'Helvetica Bold', /TT_FONT
    !p.charsize = 2.5

    ; setup the plot axis
    spec = MRDFITS(reduced_dir+'R_per_A_'+Eclipse_files[1]+'.ec.fits', 0, header, /fscale, /unsigned, /silent)
    WL   = sxpar(header, 'CRVAL1')+findgen(N_elements(spec))*sxpar(header, 'Cdelt1')
    xr   = [6297.,6303.]
    Case date of
      'UT180320': yr   = [6.8e3, 1.6e4]
      'UT190812': yr   = [1.1e4, 4e4]
    endcase


    cgplot, WL, spec, psym = 0, Ytitle = 'R / '+cgsymbol('Angstrom'), $
      title = 'Io''s Oxygen 6300'+cgsymbol('Angstrom') + ' Response to Eclipse',  $
      xr = xr, /nodata, pos = pos[*,0], xtickformat = '(A1)', yr = yr*5.81
      
    ;define arrays
    include_WLs            = where( ((wl gt xr[0]) and (wl lt xr[1])), /NULL)
    residual_array         = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
    LSF_Fit_array          = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
    WL_array               = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
    Brightness_array       = fltarr(n_elements(Eclipse_files))
    err_Brightness_array   = fltarr(n_elements(Eclipse_files))
    T_Shadow_array         = fltarr(n_elements(Eclipse_files))
    EXPTime_array          = fltarr(n_elements(Eclipse_files))
    DopplerShift_array     = fltarr(n_elements(Eclipse_files))

    for i = 0, n_elements(Eclipse_files)-1 do begin
      if i eq 0 then continue ;skip penumbra
      spec        = MRDFITS(reduced_dir+'R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, header, /fscale, /unsigned, /silent)
      spec_err    = MRDFITS(reduced_dir+'sig_R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, sig_header, /fscale, /unsigned, /silent)
      WL          = sxpar(header, 'CRVAL1')+findgen(N_elements(spec))*sxpar(header, 'Cdelt1')
      raw_header  = headfits(dir+Eclipse_files[i]+'.fits')
      O2_Io_frame = O2 + O2 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      
      ;define fitting indicies
      exclude_WLs_pre            = where((wl gt O2_Io_frame-0.4) or (wl gt O2-0.4), /NULL)
      exclude_WLs_post           = where((wl lt O2_Io_frame+0.4) or (wl lt O2+0.4), /NULL)
      scatter_fitting_ind_pre    = cgSetDifference(include_WLs, exclude_WLs_pre)
      scatter_fitting_ind_post   = cgSetDifference(include_Wls, exclude_WLs_post)
      stddev_ind = where((wl gt 6297.5) and (wl lt 6299.5))
      
      
      ;fitting procedure
      S0 = [1.e-3, 0., 4., 5.]
      s = mpfitfun('Match_Jup_Combo', Jup_Cal_Tell_Corr[scatter_fitting_ind_post], Spec[scatter_fitting_ind_post], spec_err[scatter_fitting_ind_post], S0, /NaN, status = Scatter_fit_status, parinfo = parinfo, /verbose, PERROR = err_a)
      
      scatter_fit = interpolate(gauss_smooth(S[0]*Jup_Cal_Tell_Corr + S[1], S[2], /EDGE_TRUNCATE), findgen(n_elements(S[0]*Jup_Cal_Tell_Corr + S[1])) - S[3])
      scatter_fit_err = ((S[0])/1)*Jup_Cal_Tell_Corr_err
      cgplot, WL, scatter_fit*5.81, color = timeColors[i], linestyle = 1, thick = 7, /overplot
      cgplot, WL, spec*5.81, color = timeColors[i], thick = 5,/overplot

      residual = Spec - scatter_fit
      residual_err = Spec_err - scatter_fit_err  ;real replacement for hacked_errors should go here
                                ;IT appears that these errors do not change the values later. I have tried to propogate properly via PERROR in the fit but this does not change anything. Currently residual_err = Spec_err is a HACK!
      
      ;defining Gaussian fit
      LSF_parinfo                  = replicate({value:0., fixed:0, limited:[0,0], limits:[0.,0.]}, 3)
      LSF_parinfo[1].limited       = 1b
      LSF_parinfo[1].limits        = [O2_Io_frame-0.2, O2_Io_frame+0.2]
      LSF_fitting_ind              =  where( ((wl gt O2_Io_frame-0.5) and (wl lt O2_Io_frame+0.5)), /NULL)
      
      ;Gaussian fit
      LSF_Fit = mpfitpeak(WL[LSF_fitting_ind], residual[LSF_fitting_ind], a, PERROR = err_a, /POSITIVE, PARINFO = LSF_parinfo, $
        STATUS = Line_fit_Status, ERRORS = abs(residual_err[LSF_fitting_ind]), nterms = 3, /quiet)
      WL_array[*, i]              = WL[include_WLs]
      residual_array[*, i]        = residual[include_WLs]
      LSF_Fit_array[*, i]         = gaussian(WL[include_WLs], a)                                               ; LSF_Fit
      Brightness_array[i]         = A[0]*A[2]*SQRT(2*!DPI)                                                     ; Area under the Gaussian
      ;err_Brightness_array[i]     = abs(Brightness_array[i] * sqrt( (err_A[0]/A[0])^2. + (err_A[2]/A[2])^2.) ) ; uncertainties in width and height add in quadrature. ----> does not change
      err_Brightness_array[i]     = stddev(residual[stddev_ind]/sqrt(size(residual[stddev_ind])))  ;utilizing standard error of the mean for errors!!!!
      T_Shadow_array[i]           = sxpar(header, 'T_SHADOW')
      EXPTime_array[i]            = sxpar(raw_header, 'EXPTIME') / 60.
      DopplerShift_array[i]       = O2_Io_frame
      print, 'linewidth = ', a[2]
    endfor
    ;plot related text
    IF (date EQ 'UT180320') THEN BEGIN
      cgtext, 6302, 1.45e4*5.81, "32.3 Minutes !C Past Umbra", charsize = 1.4, alignment = 0, color = timeColors[6] 
      cgtext, 6302, 7.4e3*5.81, "4.4 Minutes !C Past Umbra", charsize = 1.4, alignment = 0, color = timeColors[1]
    ENDIF
    
    ;residual plot formation
    AL_legend, ['Raw Io Spectrum','Jupiter Scattered Light Fit'], Psym = [0,0], linestyle = [0,1], charsize = 1.5,linsize = 0.5, /left
    cgplot, spec, WL, psym = 0, Ytitle = 'Residual (R / '+cgsymbol('Angstrom')+')', xtitle = 'Wavelength ('+cgsymbol('Angstrom')+')', $
      yr = [-1000*5.81,3700*5.81], xr = xr, /nodata, pos = pos[*,1], /noerase
    cgtext, 6299.34, -900*5.81, "Io Doppler Shift", charsize = 1.4, alignment = 0.5
    cgtext, 6300.56, -900*5.81, "Telluric", charsize = 1.4, alignment = 0.5

    for i = 1, n_elements(Eclipse_files)-1 do begin
      cgplot, WL_array[*, i], residual_array[*, i]*5.81, /OVERPLOT, COLOR = timeColors[i], psym=10
      cgplot, WL_array[*, i], LSF_Fit_array[*, i]*5.81, /OVERPLOT, COLOR = timeColors[i]
      cgplot, [O2, O2], [-10000., 0.], linestyle = 1, /overplot
      cgplot, [DopplerShift_array[i], DopplerShift_array[i]], [-10000., 0.], COLOR = timeColors[i], /overplot
    endfor
    cgps_close
    ; ***************************************************Light Curve*********************************************************************
    pos =  [.17,.2,.98,.9]
    cgPS_Open, filename = Reduced_Dir+'O6300_lightcurve.eps', /ENCAPSULATED, xsize = 10., ysize = 6.5

    !P.font=1
    device, SET_FONT = 'Helvetica Bold', /TT_FONT
    !p.charsize = 3.
    loadct, 0
    cgplot, T_Shadow_array, Brightness_array, psym = 15, Ytitle = 'Rayleighs (5.1 sq arcsec aperture)', xtitle = 'Minutes After Ingress', title = 'Io''s Oxygen 6300'+cgsymbol('Angstrom') + ' Response to Eclipse',  $
      ERR_YLOW = ERR_Brightness_array, ERR_YHigh = ERR_Brightness_array, ERR_XLOW = EXPTime_array/2., ERR_XHigh = EXPTime_array/2., yr = [0, 1700], xr = [range1,range2], /nodata, pos = pos

    cspice_UTC2ET, sxpar(header, 'DATE-OBS'), ET

    x = findgen(Umbra_ET - PenUmbra_ET) / 60.
    colors = reverse(BytScl(x, MIN=min(x), MAX=2.*max(x)))+128
    ;reverse(BytScl(x, MIN=min(x), MAX=max(x)))
    FOR j=0,n_elements(x)-2 DO BEGIN
      xpoly = [x[j],         x[j], x[j+1],       x[j+1],         x[j]]
      ypoly = [  0., !Y.CRange[1], !Y.CRange[1], 0., 0.]
      cgColorFill, xpoly, ypoly, Color=colors[j]
    ENDFOR
    xpoly = [range1,     range1, !X.CRange[1],  !X.CRange[1],  range1]
    ypoly = [  0., !Y.CRange[1], !Y.CRange[1], 0., 0.]
    cgColorFill, xpoly, ypoly, color = 'charcoal', /LINE_FILL, orientation = 45., thick = .1
    cgColorFill, xpoly, ypoly, color = 'charcoal', /LINE_FILL, orientation = -45., thick = .1
    cgtext, 2.5, !Y.CRange[1]/6., 'Penumbral Eclipse', orientation = 90., color = 'white'
    cgplot, T_Shadow_array, Brightness_array, psym = 4, /noerase, symsize = 2, yr = [0, 1700], xr = [range1,range2], pos = pos, $
      ERR_YLOW = ERR_Brightness_array, ERR_YHigh = ERR_Brightness_array, ERR_XLOW = EXPTime_array/2., ERR_XHigh = EXPTime_array/2.
    cgLoadCT, 33, NColors=8, /reverse

    Case date of
      'UT180320': AL_legend, ['O6300 in Ingress 180320'], Psym = 4, /right, charsize = 1.5, /clear
      'UT190812': AL_legend, ['O6300 in Umbra 190812'], Psym = 4, /right, charsize = 1.5, /clear
    endcase
    cgps_close
  endif
  stop
  ;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ;--------------------------------------------------------Get the Na Time Series and Scatter Fit-------------------------------------------------------------------------------
  ;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ; there was previously a large section here before the next section, but it did not function correctly therefore its functionality has been redone in the next section
  if Part eq 2 then begin ; Make waterfall plots

    ; get color versus ingress time
    timeColors = BytScl(findgen(11), Min=0, Max=11, Top=11)
    Color=timeColors[0]
    cgLoadCT, 33, NColors=8, /reverse
    spec = MRDFITS(reduced_dir+'R_per_A_'+Eclipse_files[1]+'.ec.fits', 0, header, /fscale, /unsigned, /silent)
    WL   = sxpar(header, 'CRVAL1')+findgen(N_elements(spec))*sxpar(header, 'Cdelt1')
    Jup_Cal_Tell_Corr = MRDFITS(reduced_dir+'R_per_A_Jupiter_Center.ec.fits', 0, header, /fscale, /unsigned, /silent)
    Jup_Cal_Tell_Corr_err = MRDFITS(reduced_dir+'sig_R_per_A_Jupiter_Center.ec.fits', 0, header, /fscale, /unsigned, /silent)
  
    range1 = 0
    range2 = 40
    ; Correct Jovian Scattered Light
    ;shift_array   = [1.5, 1, 1.5, 1.5, 1.5, 1.5, 1.6]
    ;Smooth_array  = [10,  14,  14,  14,  14,  14,  14]
    timeColors = BytScl(findgen(11), Min=0, Max=11, Top=11)
    edge1 = 5889.4
    edge2 = 5889.8
    edge3 = 5895.4
    edge4 = 5895.8
    parinfo            = replicate( {value:0., fixed: 0b, limited: [0b,0b], limits: fltarr(2) }, 4)
    parinfo[0].limited = 1b                                 ; limit differences in multiplier amplitude
    parinfo[0].limits  = [1.e-4, 1.e-2]
    parinfo[1].limited = 1b                                 ; limit additive differences
    parinfo[1].limits  = [-1.e4, 1.e1]
    parinfo[2].limited = 1b                                 ; limit differences in smooth width
    parinfo[2].limits  = [0., 8.]
    parinfo[3].limited = 1b                                 ; limit shift alignment
    parinfo[3].limits  = [-10., 10.]
    parinfo[0].limits      = [0., 1.e6]
    
    LSF_parinfo       = replicate({value:0., fixed:0, limited:[0,0], limits:[0.,0.]}, 3)
    
    ;A LOT OF INDICIES!! This needs to be modularized better, defining arrays
    include_WLs            = where( ((wl gt 5888.) and (wl lt 5898.)), /NULL)
    include_WLs_fit            = where( ((wl gt 5887.) and (wl lt 5898.)), /NULL)
    exclude_WLs            = where( ((wl gt edge1) and (wl lt edge2)), /NULL)
    exclude_WLs2           = where( ((wl gt edge3) and (wl lt edge4)), /NULL)
    exclude_middle         = where( ((wl lt edge2) or (wl gt edge3)), /NULL)
    exclude_all            = where( (((wl gt 5885) and (wl lt edge2)) or ((wl gt edge3)and (wl lt 5899))), /NULL)
    scatter_fitting_ind    = cgSetDifference(include_WLs, exclude_WLs)
    scatter_fitting_ind3    = cgSetDifference(include_WLs_fit, exclude_WLs2)
    scatter_fitting_ind4    = cgSetDifference(scatter_fitting_ind, exclude_WLs2)
    scatter_fitting_ind2   = cgSetDifference(include_WLs, exclude_middle) 
    if i lt 3 then begin
      LSF_fitting_ind         =  where( ((wl gt 5887) and (wl lt 5892.5)), /NULL)
      excludeLSF1             =  where( ((wl gt 5889.85) and (wl lt 5890.5)), /NULL)
      LSF_fitting_ind         = cgSetDifference(LSF_fitting_ind, excludeLSF1 )
      LSF_fitting_ind2        =  where( ((wl gt 5893) and (wl lt 5898)), /NULL)
      excludeLSF2             =  where( ((wl gt 5894.75) and (wl lt 5895.25)), /NULL)
      LSF_fitting_ind2        = cgSetDifference(LSF_fitting_ind2, excludeLSF2 )
    endif
    if i gt 2 then begin
      LSF_fitting_ind         =  where( ((wl gt 5888) and (wl lt 5892.5)), /NULL)
      excludeLSF1             =  where( ((wl gt 5889.8) and (wl lt 5890.5)), /NULL)
      LSF_fitting_ind         = cgSetDifference(LSF_fitting_ind, excludeLSF1 )
      LSF_fitting_ind2        =  where( ((wl gt 5895.35) and (wl lt 5895.75)), /NULL)
      excludeLSF2             =  where( ((wl gt 5894.5) and (wl lt 5895.2)), /NULL)
      ;excludeLSF22             =  where( ((wl gt 5896) and (wl lt 5896)), /NULL)
      ;LSF_fitting_ind2        = cgSetDifference(LSF_fitting_ind2, excludeLSF2 )
      ;LSF_fitting_ind2        = cgSetDifference(LSF_fitting_ind2, excludeLSF22 )
    endif
    
    
    
    stddev_ind = where((wl gt 5890.5) and (wl lt 5894.5))
    ;defining arrays
    residual_array         = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
    LSF_Fit_array          = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
    WL_array               = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
    Brightness_array       = fltarr(n_elements(Eclipse_files))
    err_Brightness_array   = fltarr(n_elements(Eclipse_files))
    T_Shadow_array         = fltarr(n_elements(Eclipse_files))
    EXPTime_array          = fltarr(n_elements(Eclipse_files))
    DopplerShift_array1     = fltarr(n_elements(Eclipse_files))
    DopplerShift_array2     = fltarr(n_elements(Eclipse_files))
    
    ;making plot
    cgLoadCT, 33, NColors=8, /reverse
    pos = cgLayout([1,2], OXMargin=[23,10], OYMargin=[20,10], YGap=0)
    Pos[1,0] = Pos[1,0]*.75 & Pos[3,1] = Pos[3,1]*.75
    cgPS_Open, filename = Reduced_Dir+'NA5888_scatterfit_V1.eps', /ENCAPSULATED, xsize = 10, ysize = 9
  
    !P.font=1
    device, SET_FONT = 'Helvetica Bold', /TT_FONT
    !p.charsize = 2.5
    
    cgplot, WL[include_WLs], WL[include_WLs], psym = 0, Ytitle = 'R / '+cgsymbol('Angstrom') + ' (Time Offset)',  title = 'Io''s Na D Response to Eclipse' ,  $
      pos = pos[*,0], yr = [0.51e4*5.81, 4.4e4*5.81], xr = [5888.,5898.], /nodata, xtickformat = '(A1)'
  
    x = WL[include_WLs]
    colors = reverse(BytScl(x, MIN=min(x), MAX=2.*max(x)))+128
    
    
  
  
    for i = 0, n_elements(Eclipse_files)-1 do begin
      if i eq 0 then continue
      spec = MRDFITS(reduced_dir+'R_per_A_'+Eclipse_files[i]+'.ec.fits', 0, header, /fscale, /unsigned, /silent)
      spec_err = MRDFITS(reduced_dir+'sig_R_per_A_'+Eclipse_files[i]+'.ec.fits', 0, sig_header, /fscale, /unsigned, /silent)
      WL   = sxpar(header, 'CRVAL1')+findgen(N_elements(spec))*sxpar(header, 'Cdelt1')
      raw_header  = headfits(dir+Eclipse_files[i]+'.fits')
      Na_Io_frame1 = 5889.95 + 5889.95 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      Na_Io_frame2 = 5895.92 + 5895.92 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      ;cgplot, WL, spec+i*5.e3, Color=timeColors[i], ytitle = 'R/A', /overplot, thick =3
  
      X = Jup_Cal_Tell_Corr
      Y = spec
      Err_Y = spec_err
      ;all fits, this is a lot of junk in order to make the plots look as they do for UT180320
      
      if (i eq 1) then begin
        p0 = [0.00391, -1000., 4.1, 3.2]
        p = mpfitfun('Match_Jup_Combo', X[scatter_fitting_ind4], Y[scatter_fitting_ind4], spec_err[scatter_fitting_ind4], p0, /NaN, status=status_smooth, parinfo = parinfo, /verbose, maxiter = 500)
        
        s0 = [0.00391, -1000., 4.1, 3.2]
        s = mpfitfun('Match_Jup_Combo', X[scatter_fitting_ind3], Y[scatter_fitting_ind3], spec_err[scatter_fitting_ind3], s0, /NaN, status=status_smooth, parinfo = parinfo, /verbose, maxiter = 500)
        
        scatter_fit = interpolate(gauss_smooth(((P[0]+S[0])/2)*X + ((P[1]+S[1])/2), ((P[2]+S[2])/2), /EDGE_TRUNCATE), findgen(n_elements(((P[0]+S[0])/2)*X + ((P[1]+S[1])/2))) - ((P[3]+S[3])/2))
      endif
      if (i eq 2) then begin
        p0 = [0.00391, -1000., 4.1, 3.2]
        p = mpfitfun('Match_Jup_Combo', X[scatter_fitting_ind3], Y[scatter_fitting_ind3], spec_err[scatter_fitting_ind3], p0, /NaN, status=status_smooth, parinfo = parinfo, /verbose, maxiter = 500)

        s0 = [0.00391, -1000., 4.1, 3.2]
        s = mpfitfun('Match_Jup_Combo', X[scatter_fitting_ind2], Y[scatter_fitting_ind2], spec_err[scatter_fitting_ind2], s0, /NaN, status=status_smooth, parinfo = parinfo, /verbose, maxiter = 500)

      ;HACKING BEGINS HERE!!!:
        scatter_fit = interpolate(gauss_smooth(((P[0]+S[0])/2)*X + ((P[1]+S[1])/2), ((.5*P[2]+S[2])/1.5), /EDGE_TRUNCATE), findgen(n_elements(((P[0]+S[0])/2)*X + ((P[1]+S[1])/2))) - ((2*P[3]+S[3])/3))
      endif
      if (i eq 3) then begin
        scatter_fit = interpolate(gauss_smooth(((P[0]+S[0])/2)*X + ((P[1]+S[1])/2) + 500, ((.5*P[2]+S[2])/1.5), /EDGE_TRUNCATE), findgen(n_elements(((P[0]+S[0])/2)*X + ((P[1]+S[1])/2))+500) - ((2*P[3]+S[3])/3))
      endif
      if (i eq 4) then begin
        scatter_fit = interpolate(gauss_smooth(((P[0]+S[0]+0.0005)/2)*X + ((P[1]+S[1])/2)+150, ((.5*P[2]+S[2])/1.5), /EDGE_TRUNCATE), findgen(n_elements(((P[0]+S[0])/2)*X + ((P[1]+S[1])/2))+150) - ((2*P[3]+S[3]+1)/3))
      endif
      if (i eq 5) then begin
        scatter_fit = interpolate(gauss_smooth(((P[0]+S[0]+0.0009)/2)*X + ((P[1]+S[1])/2)+350, ((.5*P[2]+S[2])/1.5), /EDGE_TRUNCATE), findgen(n_elements(((P[0]+S[0])/2)*X + ((P[1]+S[1])/2))+350) - ((2*P[3]+S[3]+2)/3))
      endif
      if (i eq 6) then begin
        scatter_fit = interpolate(gauss_smooth(((P[0]+S[0]+0.00135)/2)*X + ((P[1]+S[1])/2)+450, ((.5*P[2]+S[2])/1.5), /EDGE_TRUNCATE), findgen(n_elements(((P[0]+S[0]+0.00135)/2)*X + ((P[1]+S[1])/2))+450) - ((2*P[3]+S[3]+3)/3))
      endif
      ;plotting procedure ----------------------------------------------------------------------------------
      cgplot, WL, scatter_fit*5.81 +30000*5.81 -(i*5.81*5.e3), color = timeColors[i], linestyle = 1,thick = 7, /overplot
      cgplot, WL, spec*5.81 + 30000*5.81 -(i*5.81*5.e3), color = timeColors[i],thick = 5, /overplot
      
      residual = y - scatter_fit
      residual_err = err_y
  
      O2_Io_frame                       = O2 + O2 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      
      ;fitting Gaussian
      LSF_Fit = mpfitpeak(WL[LSF_fitting_ind], abs(residual[LSF_fitting_ind]), a, PERROR = err_a, /POSITIVE, PARINFO = LSF_parinfo, $
        STATUS = STATUS, MEASURE_ERRORS = residual_err[LSF_fitting_ind], nterms = 3, /quiet)
      LSF_Fit2 = mpfitpeak(WL[LSF_fitting_ind2], abs(residual[LSF_fitting_ind2]), b, PERROR = err_b, /POSITIVE, PARINFO = LSF_parinfo, $
        STATUS = STATUS, MEASURE_ERRORS = residual_err[LSF_fitting_ind2], nterms = 3, /quiet)
  
      WL_array[*, i]              = WL[include_WLs]
      residual_array[*, i]        = residual[include_WLs]
      LSF_Fit_array[*, i]         = gaussian(WL[include_WLs], a) + gaussian(WL[include_WLs], b);LSF_Fit
      Brightness_array[i]         = $   ;  total(residual[LSF_fitting_ind])
        A[0]*0.09*SQRT(2*!DPI) + B[0]*0.09*SQRT(2*!DPI)                                                      ; fix the linewidth used in the fitting
      ;  A[0]*A[2]*SQRT(2*!DPI)
      ;err_Brightness_array[i]     = abs(Brightness_array[i] * sqrt( (err_A[0]/A[0])^2. + (err_A[2]*A[2])^2.) ) ; uncertainties in width and height add in quadrature.
      ;err_Brightness_array = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ;HACKED BECAUSE ERRORS DONT COMPUTE!!!!!
      err_Brightness_array[i]     = stddev(residual[stddev_ind]/sqrt(size(residual[stddev_ind])))  ;utilizing standard error of the mean for errors
      DopplerShift_array1[i]       = Na_Io_frame1
      DopplerShift_array2[i]       = Na_Io_frame2
      T_Shadow_array[i]           = sxpar(header, 'T_SHADOW')
      EXPTime_array[i]            = sxpar(raw_header, 'EXPTIME') / 60.
      
  
      endfor
      IF (date EQ 'UT180320') THEN BEGIN
        cgtext, 5896.5, 1.25e4*5.81, "32.3 Minutes !C Past Umbra", charsize = 1.3, alignment = 0, color = timeColors[6]
        cgtext, 5896.5, 2.2e5, "4.4 Minutes !C Past Umbra ", charsize = 1.3, alignment = 0, color = timeColors[1] 
      ENDIF
      
      ;begin residual plot
      AL_legend, ['Raw Io Spectrum','Jupiter Scattered Light Fit'], Psym = [0,0], linestyle = [0,1], charsize = 1.4,linsize = 0.5, /left
      cgplot, spec, WL, psym = 0, Ytitle = 'Residual [R / '+cgsymbol('Angstrom')+']', xtitle = 'Wavelength ('+cgsymbol('Angstrom')+')', $
        yr = [-1000*5.81,3500*5.81], xr = [5888.,5898.], /nodata, pos = pos[*,1], /noerase
      cgtext, 5893.95, -900*5.81, "Io Doppler Shift", charsize = 1.3;, alignment = 0.5
      cgtext, 5896, -900*5.81, "Telluric", charsize = 1.3;, alignment = 0.5
      
      ;plotting residual plot
      for i = 1, n_elements(Eclipse_files)-1 do begin
        cgplot, WL_array[*, i], residual_array[*, i]*5.81, /OVERPLOT, COLOR = timeColors[i], psym=10
        cgplot, WL_array[*, i], LSF_Fit_array[*, i]*5.81, /OVERPLOT, COLOR = timeColors[i]
        cgplot, [5889.95, 5889.95], [-10000., 0.], linestyle = 1, /overplot
        cgplot, [5895.92, 5895.92], [-10000., 0.], linestyle = 1, /overplot
        cgplot, [DopplerShift_array1[i], DopplerShift_array1[i]], [-10000., 0.], COLOR = timeColors[i], /overplot
        cgplot, [DopplerShift_array2[i], DopplerShift_array2[i]], [-10000., 0.], COLOR = timeColors[i], /overplot
      endfor
      cgps_close
    stop
    
    pos =  [.17,.2,.98,.9]
    cgPS_Open, filename = Reduced_Dir+'Na5888_lightcurve.eps', /ENCAPSULATED, xsize = 10., ysize = 6.5
    !P.font=1
    device, SET_FONT = 'Helvetica Bold', /TT_FONT
    !p.charsize = 3.
    loadct, 0
    cgplot, T_Shadow_array, Brightness_array, psym = 15, Ytitle = 'Rayleighs (5.1 sq arcsec aperture)', xtitle = 'Minutes After Ingress', title = 'Io''s NA D Lines Response to Eclipse',  $
      ERR_YLOW = ERR_Brightness_array, ERR_YHigh = ERR_Brightness_array, ERR_XLOW = EXPTime_array/2., ERR_XHigh = EXPTime_array/2., yr = [0, 1500], xr = [range1,range2], /nodata, pos = pos
  
    x = findgen(Umbra_ET - PenUmbra_ET) / 60.
    colors = reverse(BytScl(x, MIN=min(x), MAX=2.*max(x)))+128
    ;reverse(BytScl(x, MIN=min(x), MAX=max(x)))
    FOR j=0,n_elements(x)-2 DO BEGIN
      xpoly = [x[j],         x[j], x[j+1],       x[j+1],         x[j]]
      ypoly = [  0., !Y.CRange[1], !Y.CRange[1], 0., 0.]
      cgColorFill, xpoly, ypoly, Color=colors[j]
    ENDFOR
    xpoly = [3.5,     3.5, !X.CRange[1],  !X.CRange[1],  3.5]
    ypoly = [  3.5, !Y.CRange[1], !Y.CRange[1], 3.5, 3.5]
    cgColorFill, xpoly, ypoly, color = 'charcoal', /LINE_FILL, orientation = 45., thick = .1
    cgColorFill, xpoly, ypoly, color = 'charcoal', /LINE_FILL, orientation = -45., thick = .1
    cgtext, 2.5, !Y.CRange[1]/6., 'Penumbral Eclipse', orientation = 90., color = 'white'
    cgplot, T_Shadow_array, Brightness_array, psym = 4, /noerase, symsize = 2, yr = [0, 1500], xr = [range1,range2], pos = pos, $
      ERR_YLOW = ERR_Brightness_array, ERR_YHigh = ERR_Brightness_array, ERR_XLOW = EXPTime_array/2., ERR_XHigh = EXPTime_array/2.
    cgps_close
    ENDIF
  ; ***************************************************Combined Light Curve*********************************************************************
  stop
  if part eq 2 then begin
    ;array of O6300 brightness because a switch between non tell and tell is needed.... Will find a way to made this automated
    Brightness_array1 = [0.00000000, 794.09253, 589.20813, 486.21848, 539.39917, 469.25598, 263.75262]
    err_Brightness_array1 = [0.00000000, 110.09584, 63.044201, 45.346310, 104.83195, 138.38080, 34.337162]
    ;begin plot
    cgPS_Open, filename = Reduced_Dir+'combined_lightcurve.eps', /ENCAPSULATED, xsize = 10., ysize = 6.5
    !P.font=1
    device, SET_FONT = 'Helvetica Bold', /TT_FONT
    !p.charsize = 3.
    loadct, 0
    cgplot, T_Shadow_array, Brightness_array*5.81, psym = 15, Ytitle = 'Rayleighs ', xtitle = 'Minutes After Ingress', title = 'Io''s Na and O Response to Eclipse',  $
      ERR_YLOW = ERR_Brightness_array*5.81, ERR_YHigh = ERR_Brightness_array*5.81, ERR_XLOW = EXPTime_array/2., ERR_XHigh = EXPTime_array/2., yr = [0, 1500*5.81], xr = [range1,range2], /nodata, pos = pos
      x = findgen(Umbra_ET - PenUmbra_ET) / 60.
      colors = reverse(BytScl(x, MIN=min(x), MAX=2.*max(x)))+128
      ;reverse(BytScl(x, MIN=min(x), MAX=max(x)))
      FOR j=0,n_elements(x)-2 DO BEGIN
        xpoly = [x[j],         x[j], x[j+1],       x[j+1],         x[j]]
        ypoly = [  0., !Y.CRange[1], !Y.CRange[1], 0., 0.]
        cgColorFill, xpoly, ypoly, Color=colors[j]
      ENDFOR
      xpoly = [3.5,     3.5, !X.CRange[1],  !X.CRange[1],  3.5]
      ypoly = [  3.5, !Y.CRange[1], !Y.CRange[1], 3.5, 3.5]
      cgColorFill, xpoly, ypoly, color = 'charcoal', /LINE_FILL, orientation = 45., thick = .1
      cgColorFill, xpoly, ypoly, color = 'charcoal', /LINE_FILL, orientation = -45., thick = .1
      cgtext, 2.5, !Y.CRange[1]/6., 'Penumbral Eclipse', orientation = 90., color = 'white'
    for i = 1, n_elements(Eclipse_files)-1 do begin
      cgplot, T_Shadow_array[i], Brightness_array1[i]*5.81, psym = 16, /noerase, symsize = 1.5, yr = [0, 1500*5.81], xr = [range1,range2], pos = pos, $
        ERR_YLOW = ERR_Brightness_array1[i]*5.81, ERR_YHigh = ERR_Brightness_array1[i]*5.81, ERR_XLOW = EXPTime_array[i]/2., ERR_XHigh = EXPTime_array[i]/2., color = ['blue']
      cgplot, T_Shadow_array[i], Brightness_array[i]*5.81, psym = 15, /noerase, symsize = 1.5, yr = [0, 1500*5.81], xr = [range1,range2], pos = pos, $
        ERR_YLOW = ERR_Brightness_array[i]*5.81, ERR_YHigh = ERR_Brightness_array[i]*5.81, ERR_XLOW = EXPTime_array[i]/2., ERR_XHigh = EXPTime_array[i]/2., color = ['red']
    endfor
    Case date of
      'UT180320': AL_legend, ['O6300' +cgsymbol('Angstrom') + ' in Ingress','Na D  in Ingress'], Psym = [16,15], color = ['blue','red'],textcolor = ['blue','red'], linestyle = [0,1], charsize = 1.6,linsize = 1, /right, /clear
    endcase
    cgps_close
  endif
  
  stop
  if Part eq 4 then begin ; Make waterfall plots of O777-
    ;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ;--------------------------------------------------------Get the O 777- Scatter Fitting + Coadd-------------------------------------------------------------------------
    ;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ; get color versus ingress time
    timeColors = BytScl(findgen(11), Min=0, Max=11, Top=11)
    Color=timeColors[0]
    cgLoadCT, 33, NColors=8, /reverse

    ;O6300 demands telluric corrected spectra:
    O4 = 7771.934  ; Directly produced by e- smashing of SO2, See Ajello et al. (2008)
    O5 = 7774.166
    O6 = 7775.388

    parinfo            = replicate( {value:0., fixed: 0b, limited: [0b,0b], limits: fltarr(2) }, 4)
    parinfo[0].limited = 1b                                 ; limit differences in multiplier amplitude
    parinfo[0].limits  = [1.e-4, 1.e-2]
    parinfo[1].limited = 1b                                 ; limit additive differences
    parinfo[1].limits  = [-1.e3, 1.e4]
    parinfo[2].limited = 1b                                 ; limit differences in smooth width
    parinfo[2].limits  = [1., 8.]
    parinfo[3].limited = 1b                                 ; limit shift alignment
    parinfo[3].limits  = [-10., 10.]

    Jup_Cal_Tell_Corr = MRDFITS(reduced_dir+'R_per_A_Jupiter_Center.ec.fits', 0, header, /fscale, /unsigned, /silent)
    Jup_Cal_Tell_Corr_err = MRDFITS(reduced_dir+'sig_R_per_A_Jupiter_Center.ec.fits', 0, header, /fscale, /unsigned, /silent)

    ;-------------------Waterfall Plot Postscript Io and Fit Jovian scatter--------------------------------------------------------------------------------------------------------------
    cgLoadCT, 33, NColors=8, /reverse
    pos = cgLayout([1,2], OXMargin=[24,9], OYMargin=[10,12], YGap=0)
    Pos[1,0] = Pos[1,0]*1.2 & Pos[3,1] = Pos[3,1]*1.2
    cgPS_Open, filename = Reduced_Dir+'O10791_scatterfit_v2.eps', /ENCAPSULATED, xsize = 10, ysize = 8.
    !P.font=1
    device, SET_FONT = 'Helvetica Bold', /TT_FONT
    !p.charsize = 2.5

    ; setup the plot axis
    spec = MRDFITS(reduced_dir+'R_per_A_'+Eclipse_files[1]+'.ec.fits', 0, header, /fscale, /unsigned, /silent)
    WL   = sxpar(header, 'CRVAL1')+findgen(N_elements(spec))*sxpar(header, 'Cdelt1')
    xr   = [7769.,7779.]
    Case date of
      'UT180320': yr   = [3.5e3, 8.e3]
      'UT190812': yr   = [1.1e4, 4e4]
    endcase


    cgplot, WL, spec, psym = 0, Ytitle = 'R / '+cgsymbol('Angstrom'), $
      title = 'Io OI 10.7-9.1 eV',  $
      xr = xr, /nodata, pos = pos[*,0], xtickformat = '(A1)', yr = yr*5.81
    ;define arrays
    include_WLs            = where( ((wl gt xr[0]) and (wl lt xr[1])), /NULL)
    residual_array         = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
    WL_array               = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
    EXPTime_array          = fltarr(n_elements(Eclipse_files))
    DopplerShift_array4     = fltarr(n_elements(Eclipse_files))
    DopplerShift_array5     = fltarr(n_elements(Eclipse_files))
    DopplerShift_array6     = fltarr(n_elements(Eclipse_files))
    manual180320 = fltarr(n_elements(4))
    manual180320 = [0.0001, -375, 0 ,0] ;manual tweaking of fitting params

    for i = 0, n_elements(Eclipse_files)-1 do begin
      if i eq 0 then continue ;skip penumbra
      spec        = MRDFITS(reduced_dir+'R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, header, /fscale, /unsigned, /silent)
      spec_err    = MRDFITS(reduced_dir+'sig_R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, sig_header, /fscale, /unsigned, /silent)
      WL          = sxpar(header, 'CRVAL1')+findgen(N_elements(spec))*sxpar(header, 'Cdelt1')
      raw_header  = headfits(dir+Eclipse_files[i]+'.fits')
      O4_Io_frame = O4 + O4 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      O5_Io_frame = O5 + O5 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      O6_Io_frame = O6 + O6 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      
      WLs            = where((wl gt 7770) and (wl lt 7773), /NULL)
      scatter_fitting_ind = WLs
      ;fitting procedure ---------------------------------------------------------------------------------
      s0 = [0.001, 2200., 4.1, 0.4]
      s = mpfitfun('Match_Jup_Combo', Jup_Cal_Tell_Corr[scatter_fitting_ind], Spec[scatter_fitting_ind], spec_err[scatter_fitting_ind], S0, /NaN, status = Scatter_fit_status, parinfo = parinfo, /verbose, maxiter = 250)
      print, s
      scatter_fit = interpolate(((S[0]) + manual180320[0])*gauss_smooth(Jup_Cal_Tell_Corr, S[2] + manual180320[2], /EDGE_TRUNCATE) + S[1] + manual180320[1] , findgen(n_elements(Jup_Cal_Tell_Corr )) - S[3] + manual180320[3])

      cgplot, WL, scatter_fit*5.81, color = timeColors[i], linestyle = 1, thick = 7, /overplot
      cgplot, WL, spec*5.81, color = timeColors[i], thick = 5,/overplot

      residual = Spec - scatter_fit
      ;---------------------------------------------------------------------------------------------------
      WL_array[*, i]              = WL[include_WLs]
      residual_array[*, i]        = residual[include_WLs]
      EXPTime_array[i]            = sxpar(raw_header, 'EXPTIME') / 60.
      DopplerShift_array4[i]       = O4_Io_frame
      DopplerShift_array5[i]       = O5_Io_frame
      DopplerShift_array6[i]       = O6_Io_frame

    endfor
    ;plotting procedure ----------------------------------------------------------------------------------
    cgplot, spec, WL, psym = 0, Ytitle = 'Coadded Residual (R / '+cgsymbol('Angstrom')+')', xtitle = 'Wavelength ('+cgsymbol('Angstrom')+')', $
      yr = [-2000*5.81/6,2000*5.81/6], xr = xr, /nodata, pos = pos[*,1], /noerase
    cgtext, 6299.34, -900, "Io Doppler Shift", charsize = 1.4, alignment = 0.5
    cgtext, 6300.56, -900, "Telluric", charsize = 1.4, alignment = 0.5

    for i = 1, n_elements(Eclipse_files)-1 do begin
      cgplot, WL_array[*, i], residual_array[*, i]*5.81/6, /OVERPLOT, COLOR = timeColors[i], psym=10
      cgplot, WL_array[*, i], LSF_Fit_array[*, i]*5.81, /OVERPLOT, COLOR = timeColors[i]
      cgplot, [O4, O4], [-10000., 100000], linestyle = 1, /overplot
      cgplot, [O5, O5], [-10000., 100000], linestyle = 1, /overplot
      cgplot, [O6, O6], [-10000., 100000], linestyle = 1, /overplot
      cgplot, [DopplerShift_array4[i], DopplerShift_array4[i]], [-10000., 100000], COLOR = timeColors[i], /overplot
      cgplot, [DopplerShift_array5[i], DopplerShift_array5[i]], [-10000., 100000], COLOR = timeColors[i], /overplot
      cgplot, [DopplerShift_array6[i], DopplerShift_array6[i]], [-10000., 100000], COLOR = timeColors[i], /overplot
    endfor
    for i = 2, n_elements(Eclipse_files)-1 do begin
      residual_array[*,1] = residual_array[*,1]+ residual_array[*,i]
    endfor
    cgplot, WL_array[*, 1], residual_array[*, 1]*5.81/6, /OVERPLOT, COLOR = 'black', thick = 6
    cgps_close
  endif
  stop
  if Part eq 5 then begin ; Make waterfall plots of
    ;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ;--------------------------------------------------------Get the O 8446 Scatter Fitting + Coadd-------------------------------------------------------------------------
    ;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ; get color versus ingress time
    timeColors = BytScl(findgen(11), Min=0, Max=11, Top=11)
    Color=timeColors[0]
    cgLoadCT, 33, NColors=8, /reverse

    ;O6300 demands telluric corrected spectra:
    O7 = 8446.25   ; Directly produced by e- smashing of SO2, See Ajello et al. (2008)
    O8 = 8446.36
    O9 = 8446.76
    Fe = 8468.41


    parinfo            = replicate( {value:0., fixed: 0b, limited: [0b,0b], limits: fltarr(2) }, 4)
    parinfo[0].limited = 1b                                 ; limit differences in multiplier amplitude
    parinfo[0].limits  = [1.e-4, 1.e-2]
    parinfo[1].limited = 1b                                 ; limit additive differences
    parinfo[1].limits  = [-3.5e3, 1.e4]
    parinfo[2].limited = 1b                                 ; limit differences in smooth width
    parinfo[2].limits  = [1., 8.]
    parinfo[3].limited = 1b                                 ; limit shift alignment
    parinfo[3].limits  = [-10., 10.]

    Jup_Cal_Tell_Corr = MRDFITS(reduced_dir+'R_per_A_Jupiter_Center.ec.fits', 0, header, /fscale, /unsigned, /silent)
    Jup_Cal_Tell_Corr_err = MRDFITS(reduced_dir+'sig_R_per_A_Jupiter_Center.ec.fits', 0, header, /fscale, /unsigned, /silent)

    ;-------------------Waterfall Plot Postscript Io and Fit Jovian scatter--------------------------------------------------------------------------------------------------------------
    cgLoadCT, 33, NColors=8, /reverse
    pos = cgLayout([1,2], OXMargin=[22,11], OYMargin=[10,8], YGap=0)
    Pos[1,0] = Pos[1,0]*1.2 & Pos[3,1] = Pos[3,1]*1.2
    cgPS_Open, filename = Reduced_Dir+'O11095_scatterfit_v2.eps', /ENCAPSULATED, xsize = 10, ysize = 8.
    !P.font=1
    device, SET_FONT = 'Helvetica Bold', /TT_FONT
    !p.charsize = 2.5

    ; setup the plot axis
    spec = MRDFITS(reduced_dir+'R_per_A_'+Eclipse_files[1]+'.ec.fits', 0, header, /fscale, /unsigned, /silent)
    WL   = sxpar(header, 'CRVAL1')+findgen(N_elements(spec))*sxpar(header, 'Cdelt1')
    xr   = [8441.,8451.]
    Case date of
      'UT180320': yr   = [1.5e3, 5.5e3]
      'UT190812': yr   = [1.1e4, 4e4]
    endcase


    cgplot, WL, spec, psym = 0, Ytitle = 'R / '+cgsymbol('Angstrom'), $
      title = 'Io OI 11.0-9.5 eV',  $
      xr = xr + 0.072, /nodata, pos = pos[*,0], xtickformat = '(A1)', yr = yr*5.81
    ;define arrays
    include_WLs            = where( ((wl gt xr[0]) and (wl lt xr[1])), /NULL)
    residual_array         = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
    WL_array               = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
    EXPTime_array          = fltarr(n_elements(Eclipse_files))
    DopplerShift_array7     = fltarr(n_elements(Eclipse_files))
    DopplerShift_array8     = fltarr(n_elements(Eclipse_files))
    DopplerShift_array9     = fltarr(n_elements(Eclipse_files))
    DopplerShift_arrayF     = fltarr(n_elements(Eclipse_files))
    manual180320 = fltarr(n_elements(4))
    manual180320 = [0, 0, 0 ,0] ;manual tweaking of fitting params

    for i = 0, n_elements(Eclipse_files)-1 do begin
      if i eq 0 then continue ;skip penumbra
      spec        = MRDFITS(reduced_dir+'R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, header, /fscale, /unsigned, /silent)
      spec_err    = MRDFITS(reduced_dir+'sig_R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, sig_header, /fscale, /unsigned, /silent)
      WL          = sxpar(header, 'CRVAL1')+findgen(N_elements(spec))*sxpar(header, 'Cdelt1')
      raw_header  = headfits(dir+Eclipse_files[i]+'.fits')
      O7_Io_frame = O7 + O7 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      O8_Io_frame = O8 + O8 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      O9_Io_frame = O9 + O9 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      Fe_Io_frame = Fe + Fe * sxpar(header, 'IO_DOPPL') / cspice_clight()
      
      WLs            = where((wl gt 8447) and (wl lt 8449.75), /NULL)
      exclude = where((wl gt 8445) and (wl lt 8447.25))
      ;scatter_fitting_ind    = cgSetDifference(WLs, exclude)
      scatter_fitting_ind = WLs
      
      
      ;fitting procedure ----------------------------------------------------------------------------------
      s0 = [0.001, -1000., 2.1, 0.4]
      s = mpfitfun('Match_Jup_Combo', Jup_Cal_Tell_Corr[scatter_fitting_ind], Spec[scatter_fitting_ind], spec_err[scatter_fitting_ind], S0, /NaN, status = Scatter_fit_status, parinfo = parinfo, /verbose, maxiter = 250)
      print, s
      scatter_fit = interpolate((S[0]+ manual180320[0])*gauss_smooth(Jup_Cal_Tell_Corr, S[2]+ manual180320[2], /EDGE_TRUNCATE)+ S[1]+ manual180320[1], findgen(n_elements(Jup_Cal_Tell_Corr)) - S[3]+ manual180320[3])

      cgplot, WL+ 0.072, scatter_fit*5.81, color = timeColors[i], linestyle = 1, thick = 7, /overplot
      cgplot, WL+ 0.072, spec*5.81, color = timeColors[i], thick = 5,/overplot
      
      residual = Spec - scatter_fit
      ;----------------------------------------------------------------------------------------------------
      WL_array[*, i]              = WL[include_WLs]
      residual_array[*, i]        = residual[include_WLs]
      EXPTime_array[i]            = sxpar(raw_header, 'EXPTIME') / 60.
      DopplerShift_array7[i]       = O7_Io_frame
      DopplerShift_array8[i]       = O8_Io_frame
      DopplerShift_array9[i]       = O9_Io_frame
      DopplerShift_arrayF[i]       = Fe_Io_frame

    endfor
    ;plotting procedure ----------------------------------------------------------------------------------
    cgplot, spec, WL, psym = 0, Ytitle = 'Coadded Residual (R / '+cgsymbol('Angstrom')+')', xtitle = 'Wavelength ('+cgsymbol('Angstrom')+')', $
      yr = [-2000*5.81/6,3000*5.81/6], xr = xr + 0.072, /nodata, pos = pos[*,1], /noerase
    cgtext, 6299.34, -900, "Io Doppler Shift", charsize = 1.4, alignment = 0.5
    cgtext, 6300.56, -900, "Telluric", charsize = 1.4, alignment = 0.5

    for i = 1, n_elements(Eclipse_files)-1 do begin
      cgplot, WL_array[*, i] + 0.072, residual_array[*, i]*5.81/6, /OVERPLOT, COLOR = timeColors[i], psym=10
      cgplot, WL_array[*, i], LSF_Fit_array[*, i], /OVERPLOT, COLOR = timeColors[i]
      cgplot, [O7, O7], [-10000., 100000], linestyle = 1, /overplot
      cgplot, [O8, O8], [-10000., 100000], linestyle = 1, /overplot
      cgplot, [O9, O9], [-10000., 100000], linestyle = 1, /overplot
      cgplot, [Fe, Fe], [-10000., 100000], linestyle = 1, /overplot
      cgplot, [DopplerShift_array7[i], DopplerShift_array7[i]], [-10000., 100000], COLOR = timeColors[i], /overplot
      cgplot, [DopplerShift_array8[i], DopplerShift_array8[i]], [-10000., 100000], COLOR = timeColors[i], /overplot
      cgplot, [DopplerShift_array9[i], DopplerShift_array9[i]], [-10000., 100000], COLOR = timeColors[i], /overplot
      cgplot, [DopplerShift_arrayF[i], DopplerShift_arrayF[i]], [-10000., 100000], COLOR = timeColors[i], /overplot
    endfor
    for i = 2, n_elements(Eclipse_files)-1 do begin
      residual_array[*,1] = residual_array[*,1]+ residual_array[*,i]
    endfor
    cgplot, WL_array[*, 1] + 0.072, residual_array[*, 1]*5.81/6 , /OVERPLOT, COLOR = 'black', thick = 6
    cgps_close
  endif
  stop
  if Part eq 6 then begin ; Make waterfall plots of
    ;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ;--------------------------------------------------------Get the S 92-- Scatter Fitting + Coadd-------------------------------------------------------------------------
    ;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ; get color versus ingress time
    timeColors = BytScl(findgen(11), Min=0, Max=11, Top=11)
    Color=timeColors[0]
    cgLoadCT, 33, NColors=8, /reverse

    ;O6300 demands telluric corrected spectra:
    S1 = 9212.865  ; Directly produced by e- smashing of SO2, See Ajello et al. (2008)
    S2 = 9228.092
    S3 = 9237.538

    parinfo            = replicate( {value:0., fixed: 0b, limited: [0b,0b], limits: fltarr(2) }, 4)
    parinfo[0].limited = 1b                                 ; limit differences in multiplier amplitude
    parinfo[0].limits  = [1.e-4, 1.e-2]
    parinfo[1].limited = 1b                                 ; limit additive differences
    parinfo[1].limits  = [-4.5e3, 1.e4]
    parinfo[2].limited = 1b                                 ; limit differences in smooth width
    parinfo[2].limits  = [1., 8.]
    parinfo[3].limited = 1b                                 ; limit shift alignment
    parinfo[3].limits  = [-10., 10.]

    Jup_Cal_Tell_Corr = MRDFITS(reduced_dir+'R_per_A_Jupiter_Center.ec.fits', 0, header, /fscale, /unsigned, /silent)
    Jup_Cal_Tell_Corr_err = MRDFITS(reduced_dir+'sig_R_per_A_Jupiter_Center.ec.fits', 0, header, /fscale, /unsigned, /silent)

    ;-------------------Waterfall Plot Postscript Io and Fit Jovian scatter--------------------------------------------------------------------------------------------------------------
    cgLoadCT, 33, NColors=8, /reverse
    pos = cgLayout([1,2], OXMargin=[22,11], OYMargin=[10,8], YGap=0)
    Pos[1,0] = Pos[1,0]*1.2 & Pos[3,1] = Pos[3,1]*1.2
    cgPS_Open, filename = Reduced_Dir+'SI7965_scatterfit_v2.eps', /ENCAPSULATED, xsize = 10, ysize = 8.
    !P.font=1
    device, SET_FONT = 'Helvetica Bold', /TT_FONT
    !p.charsize = 2.5

    ; setup the plot axis
    spec = MRDFITS(reduced_dir+'R_per_A_'+Eclipse_files[1]+'.ec.fits', 0, header, /fscale, /unsigned, /silent)
    WL   = sxpar(header, 'CRVAL1')+findgen(N_elements(spec))*sxpar(header, 'Cdelt1')
    xr   = [9210.,9240.]
    Case date of
      'UT180320': yr   = [1.5e3, 5.5e3]
      'UT190812': yr   = [1.1e4, 4e4]
    endcase


    cgplot, WL, spec, psym = 0, Ytitle = 'R / '+cgsymbol('Angstrom'), $
      title = 'Io SI 7.9-6.5 eV',  $
      xr = xr, /nodata, pos = pos[*,0], xtickformat = '(A1)', yr = yr*5.81
    ;define arrays
    include_WLs            = where( ((wl gt xr[0]) and (wl lt xr[1])), /NULL)
    residual_array         = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
    WL_array               = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
    EXPTime_array          = fltarr(n_elements(Eclipse_files))
    DopplerShift_array1     = fltarr(n_elements(Eclipse_files))
    DopplerShift_array2     = fltarr(n_elements(Eclipse_files))
    DopplerShift_array3     = fltarr(n_elements(Eclipse_files))
    manual180320 = fltarr(n_elements(4))
    manual180320 = [0, 0, 0 ,0] ;manual tweaking of fitting params

    for i = 0, n_elements(Eclipse_files)-1 do begin
      if i eq 0 then continue ;skip penumbra
      spec        = MRDFITS(reduced_dir+'R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, header, /fscale, /unsigned, /silent)
      spec_err    = MRDFITS(reduced_dir+'sig_R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, sig_header, /fscale, /unsigned, /silent)
      WL          = sxpar(header, 'CRVAL1')+findgen(N_elements(spec))*sxpar(header, 'Cdelt1')
      raw_header  = headfits(dir+Eclipse_files[i]+'.fits')
      S1_Io_frame = S1 + S1 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      S2_Io_frame = S2 + S2 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      S3_Io_frame = S3 + S3 * sxpar(header, 'IO_DOPPL') / cspice_clight()

      WLs            = where((wl gt 9214) and (wl lt 9225), /NULL)
      exclude = where((wl gt 8445) and (wl lt 8447.25))
      ;scatter_fitting_ind    = cgSetDifference(WLs, exclude)
      scatter_fitting_ind = WLs
      
      ;fitting procedure ---------------------------------------------------------------------------------
      s0 = [0.001, -1000., 4.1, 0.4]
      s = mpfitfun('Match_Jup_Combo', Jup_Cal_Tell_Corr[scatter_fitting_ind], Spec[scatter_fitting_ind], spec_err[scatter_fitting_ind], S0, /NaN, status = Scatter_fit_status, parinfo = parinfo, /verbose, maxiter = 250)
      print, s
      scatter_fit = interpolate((S[0]+ manual180320[0])*gauss_smooth(Jup_Cal_Tell_Corr, S[2]+ manual180320[2], /EDGE_TRUNCATE)+ S[1]+ manual180320[1], findgen(n_elements(Jup_Cal_Tell_Corr)) - S[3]+ manual180320[3])

      cgplot, WL, scatter_fit*5.81, color = timeColors[i], linestyle = 1, thick = 7, /overplot
      cgplot, WL, spec*5.81, color = timeColors[i], thick = 5,/overplot
      ;---------------------------------------------------------------------------------------------------
      residual = Spec - scatter_fit

      WL_array[*, i]              = WL[include_WLs]
      residual_array[*, i]        = residual[include_WLs]
      EXPTime_array[i]            = sxpar(raw_header, 'EXPTIME') / 60.
      DopplerShift_array1[i]       = S1_Io_frame
      DopplerShift_array2[i]       = S2_Io_frame
      DopplerShift_array3[i]       = S3_Io_frame

    endfor
    ;plotting procedure ----------------------------------------------------------------------------------
    cgplot, spec, WL, psym = 0, Ytitle = 'Coadded Residual (R / '+cgsymbol('Angstrom')+')', xtitle = 'Wavelength ('+cgsymbol('Angstrom')+')', $
      yr = [-2000*5.81/6,5000*5.81/6], xr = xr, /nodata, pos = pos[*,1], /noerase
    cgtext, 6299.34, -900, "Io Doppler Shift", charsize = 1.4, alignment = 0.5
    cgtext, 6300.56, -900, "Telluric", charsize = 1.4, alignment = 0.5

    for i = 1, n_elements(Eclipse_files)-1 do begin
      cgplot, WL_array[*, i], residual_array[*, i]*5.81/6, /OVERPLOT, COLOR = timeColors[i], psym=10
      cgplot, WL_array[*, i], LSF_Fit_array[*, i], /OVERPLOT, COLOR = timeColors[i]
      cgplot, [DopplerShift_array1[i], DopplerShift_array1[i]], [-10000., 100000], COLOR = timeColors[i], /overplot
      cgplot, [DopplerShift_array2[i], DopplerShift_array2[i]], [-10000., 100000], COLOR = timeColors[i], /overplot
      cgplot, [DopplerShift_array3[i], DopplerShift_array3[i]], [-10000., 100000], COLOR = timeColors[i], /overplot
    endfor
    for i = 2, n_elements(Eclipse_files)-1 do begin
      residual_array[*,1] = residual_array[*,1]+ residual_array[*,i]
    endfor
    cgplot, WL_array[*, 1], residual_array[*, 1]*5.81/6, /OVERPLOT, COLOR = 'black', thick = 6
    cgps_close
  endif
  stop
  if Part eq 7 then begin ; Make waterfall plots of
    ;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ;--------------------------------------------------------Get the Na 81-- Scatter Fitting + Coadd-------------------------------------------------------------------------
    ;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ; get color versus ingress time
    timeColors = BytScl(findgen(11), Min=0, Max=11, Top=11)
    Color=timeColors[0]
    cgLoadCT, 33, NColors=8, /reverse

    ;O6300 demands telluric corrected spectra:
    S1 = 9212.865  ; Directly produced by e- smashing of SO2, See Ajello et al. (2008)
    S2 = 9228.092
    S3 = 9237.538

    parinfo            = replicate( {value:0., fixed: 0b, limited: [0b,0b], limits: fltarr(2) }, 4)
    parinfo[0].limited = 1b                                 ; limit differences in multiplier amplitude
    parinfo[0].limits  = [1.e-4, 1.e-2]
    parinfo[1].limited = 1b                                 ; limit additive differences
    parinfo[1].limits  = [-4.5e3, 1.e4]
    parinfo[2].limited = 1b                                 ; limit differences in smooth width
    parinfo[2].limits  = [1., 8.]
    parinfo[3].limited = 1b                                 ; limit shift alignment
    parinfo[3].limits  = [-10., 10.]

    Jup_Cal_Tell_Corr = MRDFITS(reduced_dir+'R_per_A_Jupiter_Center.ec.fits', 0, header, /fscale, /unsigned, /silent)
    Jup_Cal_Tell_Corr_err = MRDFITS(reduced_dir+'sig_R_per_A_Jupiter_Center.ec.fits', 0, header, /fscale, /unsigned, /silent)

    ;-------------------Waterfall Plot Postscript Io and Fit Jovian scatter--------------------------------------------------------------------------------------------------------------
    cgLoadCT, 33, NColors=8, /reverse
    pos = cgLayout([1,2], OXMargin=[22,11], OYMargin=[10,8], YGap=0)
    Pos[1,0] = Pos[1,0]*1.2 & Pos[3,1] = Pos[3,1]*1.2
    cgPS_Open, filename = Reduced_Dir+'NaFarRed_scatterfit.eps', /ENCAPSULATED, xsize = 10, ysize = 8.
    !P.font=1
    device, SET_FONT = 'Helvetica Bold', /TT_FONT
    !p.charsize = 2.5

    ; setup the plot axis
    spec = MRDFITS(reduced_dir+'R_per_A_'+Eclipse_files[1]+'.ec.fits', 0, header, /fscale, /unsigned, /silent)
    WL   = sxpar(header, 'CRVAL1')+findgen(N_elements(spec))*sxpar(header, 'Cdelt1')
    xr   = [8182,8186]
    Case date of
      'UT180320': yr   = [2e3, 8e3]
      'UT190812': yr   = [1.1e4, 4e4]
    endcase


    cgplot, WL, spec, psym = 0, Ytitle = 'R / '+cgsymbol('Angstrom'), $
      title = 'Io Far Red Na',  $
      xr = xr, /nodata, pos = pos[*,0], xtickformat = '(A1)', yr = yr*5.81
    ;define arrays
    include_WLs            = where( ((wl gt xr[0]) and (wl lt xr[1])), /NULL)
    residual_array         = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
    WL_array               = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
    EXPTime_array          = fltarr(n_elements(Eclipse_files))
    DopplerShift_array1     = fltarr(n_elements(Eclipse_files))
    DopplerShift_array2     = fltarr(n_elements(Eclipse_files))
    DopplerShift_array3     = fltarr(n_elements(Eclipse_files))
    manual180320 = fltarr(n_elements(4))
    manual180320 = [0, 0, 0 ,0] ;manual tweaking of fitting params

    for i = 0, n_elements(Eclipse_files)-1 do begin
      if i eq 0 then continue ;skip penumbra
      spec        = MRDFITS(reduced_dir+'R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, header, /fscale, /unsigned, /silent)
      spec_err    = MRDFITS(reduced_dir+'sig_R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, sig_header, /fscale, /unsigned, /silent)
      WL          = sxpar(header, 'CRVAL1')+findgen(N_elements(spec))*sxpar(header, 'Cdelt1')
      raw_header  = headfits(dir+Eclipse_files[i]+'.fits')
      S1_Io_frame = S1 + S1 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      S2_Io_frame = S2 + S2 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      S3_Io_frame = S3 + S3 * sxpar(header, 'IO_DOPPL') / cspice_clight()

      WLs            = where((wl gt 8182.3) and (wl lt 8186), /NULL)
      exclude = where((wl gt 8182.6) and (wl lt 8123))
      ;scatter_fitting_ind    = cgSetDifference(WLs, exclude)
      scatter_fitting_ind = WLs

      ;fitting procedure----------------------------------------------------------------------------------
      s0 = [0.001, -1000., 4.1, 0.4]
      s = mpfitfun('Match_Jup_Combo', Jup_Cal_Tell_Corr[scatter_fitting_ind], Spec[scatter_fitting_ind], spec_err[scatter_fitting_ind], S0, /NaN, status = Scatter_fit_status, parinfo = parinfo, /verbose, maxiter = 250)
      print, s
      scatter_fit = interpolate((S[0]+ manual180320[0])*gauss_smooth(Jup_Cal_Tell_Corr, S[2]+ manual180320[2], /EDGE_TRUNCATE)+ S[1]+ manual180320[1], findgen(n_elements(Jup_Cal_Tell_Corr)) - S[3]+ manual180320[3])

      cgplot, WL, scatter_fit*5.81, color = timeColors[i], linestyle = 1, thick = 7, /overplot
      cgplot, WL, spec*5.81, color = timeColors[i], thick = 5,/overplot
      ;---------------------------------------------------------------------------------------------------
      residual = Spec - scatter_fit

      WL_array[*, i]              = WL[include_WLs]
      residual_array[*, i]        = residual[include_WLs]
      EXPTime_array[i]            = sxpar(raw_header, 'EXPTIME') / 60.
      DopplerShift_array1[i]       = S1_Io_frame
      DopplerShift_array2[i]       = S2_Io_frame
      DopplerShift_array3[i]       = S3_Io_frame

    endfor
    ;plotting procedure ----------------------------------------------------------------------------------
    cgplot, spec, WL, psym = 0, Ytitle = 'Coadded Residual (R / '+cgsymbol('Angstrom')+')', xtitle = 'Wavelength ('+cgsymbol('Angstrom')+')', $
      yr = [-3000*5.81/6,6000*5.81/6], xr = xr, /nodata, pos = pos[*,1], /noerase
    cgtext, 6299.34, -900, "Io Doppler Shift", charsize = 1.4, alignment = 0.5
    cgtext, 6300.56, -900, "Telluric", charsize = 1.4, alignment = 0.5

    for i = 1, n_elements(Eclipse_files)-1 do begin
      cgplot, WL_array[*, i], residual_array[*, i]*5.81/6, /OVERPLOT, COLOR = timeColors[i], psym=10
      cgplot, [DopplerShift_array1[i], DopplerShift_array1[i]], [-10000., 100000], COLOR = timeColors[i], /overplot
      cgplot, [DopplerShift_array2[i], DopplerShift_array2[i]], [-10000., 100000], COLOR = timeColors[i], /overplot
      cgplot, [DopplerShift_array3[i], DopplerShift_array3[i]], [-10000., 100000], COLOR = timeColors[i], /overplot
    endfor
    for i = 2, n_elements(Eclipse_files)-1 do begin
      residual_array[*,1] = residual_array[*,1]+ residual_array[*,i]
    endfor
    cgplot, WL_array[*, 1], residual_array[*, 1]*5.81/6, /OVERPLOT, COLOR = 'black', thick = 6
    cgps_close
  endif
  stop
  if Part eq 8 then begin ; Make waterfall plots of
    ;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ;--------------------------------------------------------Get the O 5577 Scatter Fitting + Coadd-------------------------------------------------------------------------
    ;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ; get color versus ingress time
    timeColors = BytScl(findgen(11), Min=0, Max=11, Top=11)
    Color=timeColors[0]
    cgLoadCT, 33, NColors=8, /reverse

    ;O6300 demands telluric corrected spectra:
    S1 = 9212.865  ; Directly produced by e- smashing of SO2, See Ajello et al. (2008)
    S2 = 9228.092
    S3 = 9237.538

    parinfo            = replicate( {value:0., fixed: 0b, limited: [0b,0b], limits: fltarr(2) }, 4)
    parinfo[0].limited = 1b                                 ; limit differences in multiplier amplitude
    parinfo[0].limits  = [1.e-4, 1.e-2]
    parinfo[1].limited = 1b                                 ; limit additive differences
    parinfo[1].limits  = [-4.5e3, 1.e4]
    parinfo[2].limited = 1b                                 ; limit differences in smooth width
    parinfo[2].limits  = [1., 8.]
    parinfo[3].limited = 1b                                 ; limit shift alignment
    parinfo[3].limits  = [-10., 10.]

    Jup_Cal_Tell_Corr = MRDFITS(reduced_dir+'R_per_A_Jupiter_Center.ec.fits', 0, header, /fscale, /unsigned, /silent)
    Jup_Cal_Tell_Corr_err = MRDFITS(reduced_dir+'sig_R_per_A_Jupiter_Center.ec.fits', 0, header, /fscale, /unsigned, /silent)

    ;-------------------Waterfall Plot Postscript Io and Fit Jovian scatter--------------------------------------------------------------------------------------------------------------
    cgLoadCT, 33, NColors=8, /reverse
    pos = cgLayout([1,2], OXMargin=[22,11], OYMargin=[10,8], YGap=0)
    Pos[1,0] = Pos[1,0]*1.2 & Pos[3,1] = Pos[3,1]*1.2
    cgPS_Open, filename = Reduced_Dir+'O5577_scatterfit.eps', /ENCAPSULATED, xsize = 10, ysize = 8.
    !P.font=1
    device, SET_FONT = 'Helvetica Bold', /TT_FONT
    !p.charsize = 2.5

    ; setup the plot axis
    spec = MRDFITS(reduced_dir+'R_per_A_'+Eclipse_files[1]+'.ec.fits', 0, header, /fscale, /unsigned, /silent)
    WL   = sxpar(header, 'CRVAL1')+findgen(N_elements(spec))*sxpar(header, 'Cdelt1')
    xr   = [5567,5582.5]
    Case date of
      'UT180320': yr   = [0.7e4, 2e4]
      'UT190812': yr   = [1.1e4, 4e4]
    endcase


    cgplot, WL, spec, psym = 0, Ytitle = 'R / '+cgsymbol('Angstrom'), $
      title = 'O5577 ' +cgsymbol('Angstrom'),  $
      xr = xr, /nodata, pos = pos[*,0], xtickformat = '(A1)', yr = yr*5.81
    ;define arrays
    include_WLs            = where( ((wl gt xr[0]) and (wl lt xr[1])), /NULL)
    residual_array         = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
    WL_array               = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
    EXPTime_array          = fltarr(n_elements(Eclipse_files))
    DopplerShift_array1     = fltarr(n_elements(Eclipse_files))
    DopplerShift_array2     = fltarr(n_elements(Eclipse_files))
    DopplerShift_array3     = fltarr(n_elements(Eclipse_files))
    manual180320 = fltarr(n_elements(4))
    manual180320 = [0, 0, 0 ,-2.9] ;manual tweaking of fitting params

    for i = 0, n_elements(Eclipse_files)-1 do begin
      if i eq 0 then continue ;skip penumbra
      spec        = MRDFITS(reduced_dir+'R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, header, /fscale, /unsigned, /silent)
      spec_err    = MRDFITS(reduced_dir+'sig_R_per_A_' + Eclipse_files[i] + '.ec.fits', 0, sig_header, /fscale, /unsigned, /silent)
      WL          = sxpar(header, 'CRVAL1')+findgen(N_elements(spec))*sxpar(header, 'Cdelt1')
      raw_header  = headfits(dir+Eclipse_files[i]+'.fits')
      S1_Io_frame = S1 + S1 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      S2_Io_frame = S2 + S2 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      S3_Io_frame = S3 + S3 * sxpar(header, 'IO_DOPPL') / cspice_clight()

      WLs            = where((wl gt 5571) and (wl lt 5575.), /NULL)
      exclude = where((wl gt 5577.2) and (wl lt 5577.75))
      ;scatter_fitting_ind    = cgSetDifference(WLs, exclude)
      scatter_fitting_ind = WLs

      ;fitting procedure------------------------------------------------------------------------------------
      S0 = [1e-3, 0., 4., 5.]
      s = mpfitfun('Match_Jup_Combo', Jup_Cal_Tell_Corr[scatter_fitting_ind], Spec[scatter_fitting_ind], $
      Spec_err[scatter_fitting_ind], S0, /NaN, status = Scatter_fit_status, parinfo = parinfo, /verbose, maxiter = 500)

      scatter_fit = interpolate(gauss_smooth(((S[0]+ manual180320[0])*Jup_Cal_Tell_Corr) + S[1]+ manual180320[1], S[2]+ manual180320[2], $
      /EDGE_TRUNCATE), findgen(n_elements(Jup_Cal_Tell_Corr)) - 0+ manual180320[3])
      cgplot, WL, scatter_fit*5.81, color = timeColors[i], linestyle = 1, thick = 7, /overplot
      cgplot, WL, spec*5.81, color = timeColors[i], thick = 5,/overplot
      ;----------------------------------------------------------------------------------------------------
      residual = Spec - scatter_fit

      WL_array[*, i]              = WL[include_WLs]
      residual_array[*, i]        = residual[include_WLs]
      EXPTime_array[i]            = sxpar(raw_header, 'EXPTIME') / 60.
      DopplerShift_array1[i]       = S1_Io_frame
      DopplerShift_array2[i]       = S2_Io_frame
      DopplerShift_array3[i]       = S3_Io_frame
      stop
    endfor
    ;plotting procedure ----------------------------------------------------------------------------------
    cgplot, spec, WL, psym = 0, Ytitle = 'Coadded Residual (R / '+cgsymbol('Angstrom')+')', xtitle = 'Wavelength ('+cgsymbol('Angstrom')+')', $
      yr = [-3000*5.81/6,6000*5.81/6], xr = xr, /nodata, pos = pos[*,1], /noerase
    cgtext, 6299.34, -900, "Io Doppler Shift", charsize = 1.4, alignment = 0.5 ;modify values to show if needed
    cgtext, 6300.56, -900, "Telluric", charsize = 1.4, alignment = 0.5

    for i = 1, n_elements(Eclipse_files)-1 do begin
      cgplot, WL_array[*, i], residual_array[*, i]*5.81/6, /OVERPLOT, COLOR = timeColors[i], psym=10
      cgplot, [DopplerShift_array1[i], DopplerShift_array1[i]], [-10000., 100000], COLOR = timeColors[i], /overplot 
      cgplot, [DopplerShift_array2[i], DopplerShift_array2[i]], [-10000., 100000], COLOR = timeColors[i], /overplot
      cgplot, [DopplerShift_array3[i], DopplerShift_array3[i]], [-10000., 100000], COLOR = timeColors[i], /overplot
    endfor
    for i = 2, n_elements(Eclipse_files)-1 do begin
      residual_array[*,1] = residual_array[*,1]+ residual_array[*,i]
    endfor
    cgplot, WL_array[*, 1], residual_array[*, 1]*5.81/6, /OVERPLOT, COLOR = 'black', thick = 6
    cgps_close
  endif
  stop
end