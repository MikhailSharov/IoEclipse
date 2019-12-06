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
        return, interpolate(gauss_smooth(P[0]*X + P[1], P[2], /EDGE_TRUNCATE), findgen(n_elements(P[0]*X + P[1])) - P[3])
    end

Pro Io_Eclipses_ARCES_PEPSI_triplets_reduction, Part=Part, Date=Date

  ;================================Part 0 Generate of Telluric Absorption Spectrum==================================================================================
  if part eq 0 then begin
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

    ;----------------------Compare 2 methods to get a decent telluric absorption spectrum---------------------------

    ; Method 1: Assume absorption *is* a unity-normalized blue fast rotator stellar spectrum
    ;           Does pretty well for Oxygen 6300
    ;           Beware that interstellar sodium absorption will screw this up
      Telluric_Absorption_BFR = MRDFITS(Calibration_Dir+'\fullspecHD_159975_BFR.0001.ec.fits', 0, header, /fscale, /silent, /unsigned )
      Telluric_WL_BFR         = sxpar(header, 'CRVAL1')+findgen(N_elements(Telluric_Absorption_BFR))*sxpar(header, 'Cdelt1')
      Telluric_AIRMASS_BFR    = sxpar(header, 'AIRMASS')
      Telluric_Absorption1    = Telluric_Absorption_BFR

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
      Telluric_Absorption2 = Telluric_Absorption_A0V

    window, 0, title = 'Unity Normalized A0V Spectral Model (Blue) vs and Measured A0V (Black: HD 155379)'
    cgplot, Telluric_WL_A0V, Model_A0V_norm, color = 'blue', xr = [5570,9300]
    cgplot, Telluric_WL_A0V, A0V_norm, color = 'green', /overplot

    window, 1, title = 'Comparison of Telluric Absorption Spectra: Blue fast rotator (black) vs. A0V Stellar Model (Red)'
    cgplot, Telluric_WL_BFR, Telluric_Absorption_BFR, xr = [5000, 6270]
    cgplot, Telluric_WL_A0V, Telluric_Absorption2, /overplot, color = 'red'
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

    ; ***********Correct the Unity-Normalized Jupiter spectrum for Telluric Absorption, use the BFR for tellurics**********       
      jup_center     = MRDFITS(reduced_dir+'fullspec'+Jupiter_Center_File+'.ec.fits', 0, header, /fscale, /silent, /unsigned )     ; fullspec, normalized to one
      jup_center_err = MRDFITS(reduced_dir+'sig_fullspec'+Jupiter_Center_File+'.ec.fits', 0, header, /fscale, /silent, /unsigned ) ; 1 sig uncertainty in fullspec
      WL             = sxpar(header, 'CRVAL1')+findgen(N_elements(jup_center))*sxpar(header, 'Cdelt1')
      aligned_telluric_absorption = interpol(Telluric_Absorption_BFR, Telluric_WL_BFR, WL) ;using A0V star correction instead of BFR

      xr = [7769,7779]  ; Wavelength region for inspection
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
        jup_center_N     = jup_center*correct_IW
        jup_center_err_N = jup_center_err*correct_IW
        
        Window, 0, Title = 'Jupiter Center: Black = FullSpec, Grey = Continuum Indices, Blue = Echelle Order Interweave Normalized', xpos = -1000, ypos=100
        cgplot, WL, jup_center, xr = xr, /ynozero  
        cgplot, WL[continuum_ind], jup_center[continuum_ind], color = 'grey', /overplot
        cgplot, WL, jup_center_N, color = 'blue', /overplot
stop
      ; And now do the same for the telluric (BFR) spectrum  
        continuum_ind = []
        step = 25 
        for i = 0, N_elements(aligned_telluric_absorption) - (N_elements(jup_center) mod step) -1, step do begin
          ind = reverse(sort(aligned_telluric_absorption[i:i+step]))
          continuum_ind = [continuum_ind, i+ind[1:5]]
        endfor
        continuum_ind = continuum_ind[UNIQ(continuum_ind, SORT(continuum_ind))]
        normalize = interpol(1./aligned_telluric_absorption[continuum_ind], WL[continuum_ind], WL)
        aligned_telluric_absorption_n = aligned_telluric_absorption*normalize
        na_absorb_ind = where( ((wl gt 5889.75) and (wl lt 5890.7)) or ((wl gt 5895.5) and (wl lt 5896.65)), /NULL)
        aligned_telluric_absorption_n[na_absorb_ind] = 1.0

        Window, 1, Title = 'Blue Fast Rotator: Black = FullSpec, Grey = Continuum Indices, Blue = Echelle Order Interweave Normalized', xpos = -1000, ypos=700
        cgplot, WL, aligned_telluric_absorption, xr = xr, /ynozero  
        cgplot, WL[continuum_ind], aligned_telluric_absorption[continuum_ind], color = 'grey', /overplot
        cgplot, WL, aligned_telluric_absorption_n, color = 'blue', /overplot

      ; Match the depth of the telluric absorption to the Jupiter Disk-Center measurement. 
        telluric_fitting_ind = where( ((wl gt 6276.4) and (wl lt 6280.)) or ((wl gt 6281.3) and (wl lt 6282.7)) or $
                                      ((wl gt 6287.0) and (wl lt 6291.)) or ((wl gt 6292.0) and (wl lt 6293.0)), /NULL)
        parinfo            = replicate( {value:0., fixed: 0b, limited: [0b,0b], limits: fltarr(2) }, 3)
        parinfo[0].limited = 1b                               ; limit differences in multiplier amplitude
        parinfo[0].limits  = [0.5, 2.]                        ; multiplicative range limits
        parinfo[1].limited = 1b                               ; limit additive differences
        parinfo[1].limits  = [-0.3, 0.3]                      ; additive adjustments should be small
        parinfo[2].Fixed   = 1b                               ; Fix spectral smoothing
        parinfo[2].value   = 0.                               ; No spectral smoothing
        X     = aligned_telluric_absorption_n
        Y     = jup_center_n
        Err_Y = jup_center_err

        p0 = [1., 0., 0.] ; guess at initial Multiply P[0], Add P[1] and Smooth P[2] to
        p = mpfitfun('Match_Reference_Spectrum', x[telluric_fitting_ind], y[telluric_fitting_ind], 0.1*err_y[telluric_fitting_ind], p0, /NaN, status=status, parinfo = parinfo) ; Fit a y = mx + b function to the spectrum at every spatial bin... SLOW
        telluric_fit = P[0]*X + P[1]
        window, 2, xpos = -1000, ypos=100
        xr = [7769, 7779]
        
        cgplot, WL, jup_center_n, xr = xr, /ynozero  
        cgplot, WL, aligned_telluric_absorption_n, color = 'red', /overplot
        cgplot, WL[telluric_fitting_ind], jup_center_n[telluric_fitting_ind], /overplot, color = 'green', psym=14
        cgplot, WL, telluric_fit, color = 'blue', /overplot
stop
      expected_flux = interpol(Rayleighs_per_angstrom, WL_A, WL)   ; move expected flux to the Jovian Doppler-shift, UNITS are R / A
      smoothed_expected_flux = GAUSS_SMOOTH(expected_flux, 1.9)    ; this smoothing looks about right for APO/ARCES, UNITS are R / A

      ; We've been working with unity-normalized spectra until now. Determine the DN / S from the "Sum" files
        jup_sum        = MRDFITS(reduced_dir+'sum'+Jupiter_Center_File+'.ec.fits', 0, header, /fscale, /silent, /unsigned ) 
        raw_header     = headfits(dir+''+Jupiter_Center_File+'.fits')
        jup_sum        = jup_sum / sxpar(raw_header, 'EXPTIME')    ; Convert to DN per second
        jup_sum        = jup_sum / telluric_fit ; <-----suspect?!
        
    ; Determine the instrumental sensitivity from the expected versus measured flux at Jupiter Disk Center. This should be a smooth function    
      window, 0
      Sensitivity       = jup_sum / smoothed_expected_flux         ; Sensitivity in (DN / S) / (R / A)           
      Sensitivity_Curve = smooth(sensitivity, 5000, /edge_truncate, /nan)
      cgplot, WL, Sensitivity, xr = [3700, 10000], yr = [0.,0.0025], Ytitle = 'Measured Flux Sensitivity (DN/S) / (R/A)', Xtitle = 'Angstroms'
      fit_Sensitivity   = Sensitivity
      fit_sensitivity[where( (sensitivity gt 0.002) or (sensitivity lt 0.0), /Null)] = !values.F_NaN  ;some basic rejection criterion
      Sensitivity_Curve = smooth(fit_sensitivity, 5000, /edge_truncate, /nan)
      fit_sensitivity[where( (Sensitivity_Curve/fit_Sensitivity gt 1.5) or (Sensitivity_Curve/fit_Sensitivity lt 0.5), /Null)] = !values.F_NaN  ; further rejection criterion
      Sensitivity_Curve = smooth(fit_sensitivity, 5000, /edge_truncate, /nan)
      cgplot, WL, Sensitivity_Curve, color = 'red', /overplot  

      cspice_UTC2ET, PenUmbra_UTC, PenUmbra_ET
      cspice_UTC2ET, Umbra_UTC, Umbra_ET
      
      ;Create the jupiter calibration tellurically corrected
      tell_fix_ind = where( ((wl gt 6290.) and (wl lt 6310.)) and ((wl gt 5891.8) and (wl lt 5892.8)) and ((wl gt 5897.5) and (wl lt 5889.3)))
      
      IF (date EQ 'UT180320') THEN BEGIN
      Jup_Cal_Tell_Corr = jup_center_n
      Jup_cal_Tell_Corr[tell_fix_ind] = Jup_cal_Tell_Corr[tell_fix_ind]/telluric_fit[tell_fix_ind]
      Jup_Cal_Tell_Corr = Jup_Cal_Tell_Corr / Sensitivity_Curve
      Jup_Cal_Tell_Corr_Err = jup_center_err_n
      Jup_Cal_Tell_Corr_Err[tell_fix_ind] = Jup_Cal_Tell_Corr_Err[tell_fix_ind]/telluric_fit[tell_fix_ind]
      Jup_Cal_Tell_Corr_Err = Jup_Cal_Tell_Corr_Err / Sensitivity_Curve
      ENDIF ELSE BEGIN
      Jup_Cal_Tell_Corr = jup_center_n/telluric_fit
      Jup_Cal_Tell_Corr_Err = jup_center_err_n/telluric_fit
      Jup_Cal_Tell_Corr = Jup_Cal_Tell_Corr / Sensitivity_Curve
      Jup_Cal_Tell_Corr_Err = Jup_Cal_Tell_Corr_Err / Sensitivity_Curve
      ENDELSE
      
      
      window, 3, title = 'Io Eclipse Spectrum: Black = Interweave normalized, Blue = Telluric Abs, Red = Telluric corrected', xs = 1900
      ; Correct tellurics and write the spectra into R/A units...
      for i = 0, n_elements(Eclipse_files)-1 do begin
       
        Io_ecl     = MRDFITS(reduced_dir+'fullspec'+Eclipse_files[i]+'.ec.fits', 0, full_header, /fscale, /unsigned, /silent)
        Io_ecl_err = MRDFITS(reduced_dir+'sig_fullspec'+Eclipse_files[i]+'.ec.fits', 0, sig_header, /fscale, /unsigned, /silent )
        Y          = Io_ecl*Correct_IW
        Err_Y      = Io_ecl_err*Correct_IW
        
        Print, 'Sub-Pixel shifts are evident and should be fine-tuned here. See:' 
        print, C_CORRELATE( x[telluric_fitting_ind], y[telluric_fitting_ind], [-2,-1,0,1,2])

        p = mpfitfun('Match_Reference_Spectrum', x[telluric_fitting_ind], y[telluric_fitting_ind], err_y[telluric_fitting_ind], p0, /NaN, status=status, parinfo = parinfo, /quiet) ; Fit a y = mx + b function to the spectrum at every spatial bin... SLOW
        telluric_fit = P[0]*X + P[1]
        
      ; Inspect
        Print, 'Why does telluric correcction here depend on the input Jupiter Spectrum? Something Needs Fixing!'  
        cgplot, WL, Y, xr = [6275,6305], color = 'Black', yr = [0.65,1.15]
        cgplot, WL, telluric_fit, /overplot, color = 'Blue'
        cgplot, WL, Y/telluric_fit, /overplot, color = 'Red'
        
        Io_sum        = MRDFITS(reduced_dir+'sum'+Eclipse_files[i]+'.ec.fits', 0, header, /fscale, /silent, /unsigned )
        raw_header    = headfits(dir+Eclipse_files[i]+'.fits')
        Io_sum        = Io_sum / sxpar(raw_header, 'EXPTIME')    ; Convert to DN per second
        Io_sum        = Io_sum / Sensitivity_Curve               ; Convert to Rayleighs per Angstrom
        
        cgplot, WL, Io_Sum, ytitle = 'Rayleighs / Angstrom', xr = [3700, 10000]
        Io_Sum_Curve = smooth(Io_Sum, 2500, /edge_truncate, /nan)
        cgplot, WL, smooth(Io_Sum, 2500, /edge_truncate, /nan), /overplot, color = 'Red'
        cgplot, WL, Y*smooth(Io_Sum, 2500, /edge_truncate, /nan)/telluric_fit, /overplot, color = 'green'

        Y[tell_fix_ind] = Y[tell_fix_ind]/telluric_fit[tell_fix_ind]
        Err_Y[tell_fix_ind] = Err_Y[tell_fix_ind]/telluric_fit[tell_fix_ind]
        IF (date EQ 'UT180320') THEN BEGIN
        Spec_Cal_Tell_Corr = Y * smooth(Io_Sum, 2500, /edge_truncate, /nan);/telluric_fit ; Jupiter Spectrum Calibrated into R / A and corrected for telluric absorption
        Spec_Cal_Tell_Corr_err = Err_Y * smooth(Io_Sum, 2500, /edge_truncate, /nan);/telluric_fit ; Propogation of errors into R/A and corrected for telluric absorption
        ENDIF ELSE BEGIN
         Spec_Cal_Tell_Corr = Y * smooth(Io_Sum, 2500, /edge_truncate, /nan)/telluric_fit ; Jupiter Spectrum Calibrated into R / A and corrected for telluric absorption
         Spec_Cal_Tell_Corr_err = Err_Y * smooth(Io_Sum, 2500, /edge_truncate, /nan)/telluric_fit ; Propogation of errors into R/A and corrected for telluric absorption
        ENDELSE
        ;window, 0, xpos = -1000, ypos=100
        cgplot, wl, Spec_Cal_Tell_Corr, xr = [6285, 6305]
        ; find the instantaneous Earth-Io Doppler Shift
          cspice_UTC2ET, sxpar(raw_header, 'DATE-OBS'), ET
          ET_mid_exposure = ET + float(sxpar(raw_header, 'EXPTIME'))/2.
          cspice_spkezr, 'Io', ET_mid_exposure, 'J2000', 'LT+S', 'Earth', Io_Earth_State, ltime
          theta  = cspice_vsep(Io_Earth_state[0:2], Io_Earth_state[3:5])
          Io_wrt_Earth_Dopplershift = cos(theta) * norm(Io_Earth_State[3:5])
          SXADDPAR, Header, 'Io_DOPPL', Io_wrt_Earth_Dopplershift, 'Io-Earth V_radial in km/s (mid exposure)'
          SXADDPAR, Header, 'T_SHADOW', (ET_mid_exposure-PenUmbra_ET) / 60., 'Minutes since Penumbral ingress'

        MWRFITS, Spec_Cal_Tell_Corr, reduced_dir+'R_per_A_'+Eclipse_files[i]+'.ec.fits', header, /CREATE ;/create overwrites
        MWRFITS, Spec_Cal_Tell_Corr_err, reduced_dir+'sig_R_per_A_'+Eclipse_files[i]+'.ec.fits', sig_header, /CREATE, /silent ;/create overwrites for sig values (not fancy mate :( will be fixed later)
      endfor
    endif ; Part eq 0
stop
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
  
  ; Correct Jovian Scattered Light
    shift_array   = [2.9, 0, 0, 0, 0, 0, 0]  
    Smooth_array  = [0, 0, 0, 0, 0, 0, 0]  
    
    parinfo            = replicate( {value:0., fixed: 0b, limited: [0b,0b], limits: fltarr(2) }, 4)
    parinfo[0].limited = 1b                                 ; limit differences in multiplier amplitude
    parinfo[0].limits  = [-1.e3, 1.e7]
    parinfo[1].limited = 1b                                 ; limit additive differences
    parinfo[1].limits  = [-1.e7, 1.e7]
    parinfo[2].limited   = 1b                                 ; fix differences in smooth width
    parinfo[2].limits  = [-1.e7, 1.e7]
    parinfo[3].limited   = 1b                                 ; fix differences in smooth width
    parinfo[3].limits  = [-1.e7, 1.e7]
    parinfo[0].limits      = [0., 1.e7]
    Gaussian_parinfo       = replicate({value:0., fixed:0, limited:[0,0], limits:[0.,0.]}, 3)    
    
    include_WLs            = where( ((wl gt 6293.) and (wl lt 6308.)), /NULL)
    IF (date EQ 'UT180320') THEN BEGIN
    exclude_WLs            = where( ((wl gt 6299.6) and (wl lt 6300.2)), /NULL)
    ENDIF ELSE BEGIN
    exclude_WLs            = where( ((wl gt 6300.4) and (wl lt 6300.85)), /NULL)
    ENDELSE
    scatter_fitting_ind    = cgSetDifference(include_WLs, exclude_WLs)
    IF (date EQ 'UT180320') THEN BEGIN
    LSF_fitting_ind        =  where( ((wl gt 6299.5) and (wl lt 6300.2)), /NULL)
    ENDIF ELSE BEGIN
    LSF_fitting_ind        =  where( ((wl gt 6300.4) and (wl lt 6300.85)), /NULL)
    ENDELSE
    
    
    residual_array         = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
    LSF_Fit_array          = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
    WL_array               = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
    Brightness_array       = fltarr(n_elements(Eclipse_files))
    err_Brightness_array   = fltarr(n_elements(Eclipse_files))
    T_Shadow_array         = fltarr(n_elements(Eclipse_files))
    EXPTime_array          = fltarr(n_elements(Eclipse_files))
    DopplerShift_array     = fltarr(n_elements(Eclipse_files))
    window, 0, xs = 800, ys = 1000
    ;cgplot, WL, WL, xr = [6296,6302], yr = [0.9e4, 8.e4], ytitle = 'R/A', /nodata
    ;cgplot, WL, spec, xr = [5888,5898], yr = [1.e3, 1.e5], ytitle = 'R/A', /nodata
    
    ;begin plotting postscript scatterfit and spec--------------------------------------------------------------------------------------------------------------
    pos =  [.17,.2,.98,.9]
    cgPS_Open, filename = Reduced_Dir+'O6300_scatterfit_V1.eps', /ENCAPSULATED, xsize = 14.5, ysize = 7.

    !P.font=1
    device, SET_FONT = 'Helvetica Bold', /TT_FONT
    !p.charsize = 2.5
    loadct, 0
    cgplot, WL[include_WLs], WL[include_WLs], psym = 0, Ytitle = 'Rayleighs (5.1 sq arcsec aperture)[Time Offset]', xtitle = 'Wavelength ('+cgsymbol('Angstrom')+')', title = 'Io''s Oxygen 6300'+cgsymbol('Angstrom') + ' Response to Eclipse[Jupiter Background Subtraction]' ,  $
      yr = [0.9e4, 3.5e4], xr = [6298,6302], /nodata, pos = pos
    
    x = WL[include_WLs]
    colors = reverse(BytScl(x, MIN=min(x), MAX=2.*max(x)))+128
    
    for i = 0, n_elements(Eclipse_files)-1 do begin
      spec = MRDFITS(reduced_dir+'R_per_A_'+Eclipse_files[i]+'.ec.fits', 0, header, /fscale, /unsigned, /silent)
      spec_err = MRDFITS(reduced_dir+'sig_R_per_A_' +Eclipse_files[i]+'.ec.fits', 0, sig_header, /fscale, /unsigned, /silent)
      WL   = sxpar(header, 'CRVAL1')+findgen(N_elements(spec))*sxpar(header, 'Cdelt1')
      raw_header  = headfits(dir+Eclipse_files[i]+'.fits')
      ;cgplot, WL, spec+i*4.e3, Color=timeColors[i], ytitle = 'R/A', /overplot, thick =3
      
      X = Jup_Cal_Tell_Corr
      X_err = Jup_Cal_Tell_Corr_err
      Y = spec
      Err_Y = spec_err
      
      s0 = [10., 1200., 1.2, 1.3]
      
      p = mpfitfun('Match_Jup_Combo', X[scatter_fitting_ind], Y[scatter_fitting_ind], spec_err[scatter_fitting_ind], s0, /NaN, status=status_smooth, parinfo = parinfo, /verbose)
      scatter_fit = interpolate(gauss_smooth(P[0]*X + P[1], P[2], /EDGE_TRUNCATE), findgen(n_elements(P[0]*X + P[1])) - P[3])
  
      ;p = mpfitfun('Match_Reference_Spectrum', x[scatter_fitting_ind], y[scatter_fitting_ind], err_y[scatter_fitting_ind], p0, /NaN, status=status, parinfo = parinfo, /quiet) 
        ;scatter_fit = P[0]*X + P[1]
      ;cgplot, WL, scatter_fit + i*4.e3, /OVERPLOT, COLOR = timeColors[i], linestyle =3 ; offset each plot
      
      cgplot, WL[include_WLs], scatter_fit[include_WLs]+ 24000 -(i*4.e3),/noerase, yr = [0.9e4, 3.5e4], linestyle = 1, xr = [6298,6302], pos = pos
      cgplot, WL[include_WLs], spec[include_WLs]+ 24000 -(i*4.e3),yr = [0.9e4, 3.5e4],xr = [6298,6302],color = timeColors[i], /overplot
      
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
      
      LSF_Fit = mpfitpeak(WL[LSF_fitting_ind], residual[LSF_fitting_ind], a, PERROR = err_a, /POSITIVE, PARINFO = Gaussian_parinfo, $
                          STATUS = STATUS, ERRORS = abs(residual_err[LSF_fitting_ind]), nterms = 3, /quiet)
        WL_array[*, i]              = WL[include_WLs]
        residual_array[*, i]        = residual[include_WLs]
        LSF_Fit_array[*, i]         = gaussian(WL[include_WLs], a) ;LSF_Fit
        Brightness_array[i]         = $   ;  total(residual[LSF_fitting_ind])
                                      A[0]*A[2]*SQRT(2*!DPI)
                                      ;A[0]*0.09*SQRT(2*!DPI)                                                      ; HACK HAck HAck fix the linewidth used in the fitting
                                          ;  A[0]*A[2]*SQRT(2*!DPI)
        ;err_Brightness_array[i]     = abs(Brightness_array[i] * sqrt( (err_A[0]/A[0])^2. + (err_A[2]*A[2])^2.) ) ; uncertainties in width and height add in quadrature.
        err_Brightness_array = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ;HACKED BECAUSE ERRORS DONT COMPUTE!!!!!
        T_Shadow_array[i]           = sxpar(header, 'T_SHADOW')
        EXPTime_array[i]            = sxpar(raw_header, 'EXPTIME') / 60.
        DopplerShift_array[i]       = O2_Io_frame
      print, 'linewidth = ', a[2]
    endfor
    AL_legend, ['Raw Io Data','Jupiter Scatter Fit'], Psym = [0,0],linestyle = [0,1], /right, charsize = 1.5
    cgLoadCT, 33, NColors=8, /reverse

    cgps_close
    
    ;Plotting time series of O6300--------------------------------------------------------------------------------------------------------------------------------------------------
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
              ERR_YLOW = ERR_Brightness_array, ERR_YHigh = ERR_Brightness_array, ERR_XLOW = EXPTime_array/2., ERR_XHigh = EXPTime_array/2., yr = [0, 1700], xr = [range1,range2], /nodata, pos = pos
  
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
    IF (date EQ 'UT180320') THEN BEGIN
    AL_legend, ['O6300 in Ingress 180320'], Psym = 4, /right, charsize = 1.5
    ENDIF ELSE BEGIN
    AL_legend, ['O6300 in Umbra 190812'], Psym = 4, /right, charsize = 1.5
    ENDELSE
    cgps_close
    
stop
    ;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ;--------------------------------------------------------Get the Na Time Series and Scatter Fit-------------------------------------------------------------------------------
    ;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ; there was previously a large section here before the next section, but it did not function correctly therefore its functionality has been redone in the next section

    ; Correct Jovian Scattered Light
    ;shift_array   = [1.5, 1, 1.5, 1.5, 1.5, 1.5, 1.6]
    ;Smooth_array  = [10,  14,  14,  14,  14,  14,  14]
    
    edge1 = 5889.35
    edge2 = 5889.86
    edge3 = 5895.25
    edge4 = 5895.8
    
    parinfo[0].limits      = [0., 1.e6]
    Gaussian_parinfo       = replicate({value:0., fixed:0, limited:[0,0], limits:[0.,0.]}, 3)

    include_WLs            = where( ((wl gt 5888.) and (wl lt 5898.)), /NULL)
    exclude_WLs            = where( ((wl gt edge1) and (wl lt edge2)), /NULL)
    exclude_WLs2            = where( ((wl gt edge3) and (wl lt edge4)), /NULL)
    scatter_fitting_ind    = cgSetDifference(include_WLs, exclude_WLs)
    scatter_fitting_ind    = cgSetDifference(scatter_fitting_ind, exclude_WLs2)
    LSF_fitting_ind        =  where( ((wl gt 5889) and (wl lt 5889.7)), /NULL)

    residual_array         = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
    LSF_Fit_array          = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
    WL_array               = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
    Brightness_array       = fltarr(n_elements(Eclipse_files))
    err_Brightness_array   = fltarr(n_elements(Eclipse_files))
    T_Shadow_array         = fltarr(n_elements(Eclipse_files))
    EXPTime_array          = fltarr(n_elements(Eclipse_files))
    DopplerShift_array     = fltarr(n_elements(Eclipse_files))
    window, 0, xs = 800, ys = 1000
    ;cgplot, WL, spec, xr = [5882.,5902.], yr = [5.e3, 4.e4], ytitle = 'R/A', /nodata
    ;cgplot, WL, spec, xr = [6556.,6568.], yr = [5.e3, 4.e4], ytitle = 'R/A', /nodata
    pos =  [.17,.2,.98,.9]
    cgPS_Open, filename = Reduced_Dir+'NA5888_scatterfit_V1.eps', /ENCAPSULATED, xsize = 14.5, ysize = 7.

    !P.font=1
    device, SET_FONT = 'Helvetica Bold', /TT_FONT
    !p.charsize = 2.5
    loadct, 0
    cgplot, WL[include_WLs], WL[include_WLs], psym = 0, Ytitle = 'Rayleighs (5.1 sq arcsec aperture)[Time Offset]', xtitle = 'Wavelength ('+cgsymbol('Angstrom')+')', title = 'Io''s Na 5888'+cgsymbol('Angstrom') + ' Response to Eclipse[Jupiter Background Subtraction]' ,  $
      yr = [0.6e4, 5e4], xr = [5888.,5898.], /nodata, pos = pos

    x = WL[include_WLs]
    colors = reverse(BytScl(x, MIN=min(x), MAX=2.*max(x)))+128
    
    
    
    for i = 1, n_elements(Eclipse_files)-1 do begin
      spec = MRDFITS(reduced_dir+'R_per_A_'+Eclipse_files[i]+'.ec.fits', 0, header, /fscale, /unsigned, /silent)
      spec_err = MRDFITS(reduced_dir+'sig_R_per_A_'+Eclipse_files[i]+'.ec.fits', 0, sig_header, /fscale, /unsigned, /silent)
      WL   = sxpar(header, 'CRVAL1')+findgen(N_elements(spec))*sxpar(header, 'Cdelt1')
      raw_header  = headfits(dir+Eclipse_files[i]+'.fits')
      ;cgplot, WL, spec+i*5.e3, Color=timeColors[i], ytitle = 'R/A', /overplot, thick =3
      
      X = Jup_Cal_Tell_Corr
      Y = spec
      
      p0 = [1.1, 100, 0., 0.]
      p = mpfitfun('Match_Reference_Spectrum', x[scatter_fitting_ind], y[scatter_fitting_ind], err_y[scatter_fitting_ind], p0, /NaN, status=status_specfit, parinfo = parinfo, /quiet)
      scatter_fit = P[0]*X + P[1]
      
      s0 = [1.5, 1.0, 2.0, 2.1]
      s = mpfitfun('Jup_Shift_Smooth', scatter_fit[scatter_fitting_ind], spec[scatter_fitting_ind], spec_err[scatter_fitting_ind], s0, /NaN, status=status_smooth, parinfo = parinfo, /verbose)
      scatter_fit = interpolate(gauss_smooth(scatter_fit, S[0]), findgen(n_elements(scatter_fit)) - S[1])
      
      ;X = interpolate(smooth(Jup_Cal_Tell_Corr, smooth_array[i]), findgen(n_elements(Jup_Cal_Tell_Corr)) - shift_array[i])   using guessed scatter
      Err_Y = spec_err
      ;s0 = [5., 2.5, 1.0, 0.8]
      ;p = mpfitfun('Match_Jup_Combo', X[scatter_fitting_ind], Y[scatter_fitting_ind], spec_err[scatter_fitting_ind], s0, /NaN, status=status_smooth, parinfo = parinfo, /verbose)
      ;scatter_fit = interpolate(gauss_smooth(P[0]*X + P[1], P[2], /EDGE_TRUNCATE), findgen(n_elements(P[0]*X + P[1])) - P[3])
      ;p = mpfitfun('Match_Jup_Combo', X[include_WLs], Y[include_WLs], spec_err[include_WLs], s0, /NaN, status=status_smooth, parinfo = parinfo, /verbose)
      ;scatter_fit_2 = interpolate(gauss_smooth(P[0]*X + P[1], P[2], /EDGE_TRUNCATE), findgen(n_elements(P[0]*X + P[1])) - P[3])
      
      ;scatter_fit = (scatter_fit + i/2.5*scatter_fit_2)/(1.+i/2.5)
      
      ;s = mpfitfun('Jup_Shift_Smooth', scatter_fit[scatter_fitting_ind], Y[scatter_fitting_ind], spec_err[scatter_fitting_ind], s0, /NaN, status=status_smooth, parinfo = parinfo, /verbose)
      ;scatter_fit = interpolate(gauss_smooth(scatter_fit, S[0], /EDGE_TRUNCATE), findgen(n_elements(scatter_fit)) - S[1])
      
      ;cgplot, WL, scatter_fit + i*5.e3, /OVERPLOT, COLOR = timeColors[i], linestyle =3 ; offset each plot
      ;cgplot, [edge1, edge1] , [0, 1e6], /overplot ; plotting the edge cases of the scatter fitting indicies
      ;cgplot, [edge2, edge2] , [0, 1e6], /overplot
      ;cgplot, [edge3, edge3] , [0, 1e6], /overplot
      ;cgplot, [edge4, edge4] , [0, 1e6], /overplot
      residual = y - scatter_fit
      residual_err = err_y
      
      cgplot, WL[include_WLs], scatter_fit[include_WLs]+ 35000 -(i*5.e3),/noerase, yr = [0.6e4, 5e4], linestyle = 1, xr = [5888.,5898.], pos = pos
      cgplot, WL[include_WLs], spec[include_WLs]+ 35000 -(i*5.e3),yr = [0.6e4, 5e4],xr = [5888.,5898.],COLOR = timeColors[i], /overplot

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
      ;err_Brightness_array[i]     = abs(Brightness_array[i] * sqrt( (err_A[0]/A[0])^2. + (err_A[2]*A[2])^2.) ) ; uncertainties in width and height add in quadrature.
      err_Brightness_array = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ;HACKED BECAUSE ERRORS DONT COMPUTE!!!!!
      
      T_Shadow_array[i]           = sxpar(header, 'T_SHADOW')
      EXPTime_array[i]            = sxpar(raw_header, 'EXPTIME') / 60.
      DopplerShift_array[i]       = O2_Io_frame
      print, 'linewidth = ', a[2]
    endfor
    AL_legend, ['Raw Io Data','Jupiter Scatter Fit'], Psym = [0,0],linestyle = [0,1], /right, charsize = 1.5
    cgLoadCT, 33, NColors=8, /reverse
    cgps_close

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
      ERR_YLOW = ERR_Brightness_array, ERR_YHigh = ERR_Brightness_array, ERR_XLOW = EXPTime_array/2., ERR_XHigh = EXPTime_array/2., yr = [0, 1500], xr = [range1,range2], /nodata, pos = pos

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
    cgplot, T_Shadow_array, Brightness_array, psym = 4, /noerase, symsize = 2, yr = [0, 1500], xr = [range1,range2], pos = pos, $
      ERR_YLOW = ERR_Brightness_array, ERR_YHigh = ERR_Brightness_array, ERR_XLOW = EXPTime_array/2., ERR_XHigh = EXPTime_array/2.
    cgps_close
    ;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ;--------------------------------------------------------Get the OI triplet 10.7 - 9.1eV Scatter Fitting-------------------------------------------------------------------------
    ;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    view = [7769,7779]

    IF (date EQ 'UT180320') THEN BEGIN
      range1 = 0
      range2 = 40
    ENDIF ELSE BEGIN
      range1 = 70
      range2 = 135
    ENDELSE

    ; Correct Jovian Scattered Light
    shift_array   = [2.9, 0, 0, 0, 0, 0, 0]
    Smooth_array  = [0, 0, 0, 0, 0, 0, 0]

    parinfo            = replicate( {value:0., fixed: 0b, limited: [0b,0b], limits: fltarr(2) }, 4)
    parinfo[0].limited = 1b
    parinfo[0].limits  = [0., 1.e2]                                 ; limit differences in multiplier amplitude
    parinfo[1].limited = 1b                                 ; limit additive differences
    parinfo[1].limits  = [-1.e4, 1.e4]
    parinfo[2].limited   = 1b                                 ; fix differences in smooth width
    parinfo[2].limits  = [-1.e2, 1.e2]
    parinfo[3].limited   = 1b                                 ; fix differences in smooth width
    parinfo[3].limits  = [-1.e2, 1.e2]
    Gaussian_parinfo       = replicate({value:0., fixed:0, limited:[0,0], limits:[0.,0.]}, 3)

    include_WLs            = where( ((wl gt 7769.) and (wl lt 7779.)), /NULL)
    ;IF (date EQ 'UT180320') THEN BEGIN
    ;  exclude_WLs            = where( ((wl gt 6299.6) and (wl lt 6300.2)), /NULL))
    ;ENDIF ELSE BEGIN
    ;  exclude_WLs            = where( ((wl gt 6300.4) and (wl lt 6300.85)), /NULL))
    ;END ELSE
    scatter_fitting_ind    = include_WLs ;cgSetDifference(include_WLs, exclude_WLs)
    ;IF (date EQ 'UT180320') THEN BEGIN
    ;LSF_fitting_ind        =  where( ((wl gt 6299.5) and (wl lt 6300.2)), /NULL)
    ;ENDIF ELSE BEGIN
    ;LSF_fitting_ind        =  where( ((wl gt 6300.4) and (wl lt 6300.85)), /NULL)
    ;ENDELSE


    residual_array         = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
    LSF_Fit_array          = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
    WL_array               = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
    Brightness_array       = fltarr(n_elements(Eclipse_files))
    err_Brightness_array   = fltarr(n_elements(Eclipse_files))
    T_Shadow_array         = fltarr(n_elements(Eclipse_files))
    EXPTime_array          = fltarr(n_elements(Eclipse_files))
    DopplerShift_array4     = fltarr(n_elements(Eclipse_files))
    DopplerShift_array5     = fltarr(n_elements(Eclipse_files))
    DopplerShift_array6     = fltarr(n_elements(Eclipse_files))
    
    pos =  [.17,.2,.98,.9]
    cgPS_Open, filename = Reduced_Dir+'OI_107_91_scatterfit_V1.eps', /ENCAPSULATED, xsize = 14.5, ysize = 7.

    !P.font=1
    device, SET_FONT = 'Helvetica Bold', /TT_FONT
    !p.charsize = 2.5
    loadct, 0
    cgplot, WL[include_WLs], WL[include_WLs], psym = 0, Ytitle = 'Rayleighs (5.1 sq arcsec aperture)[Time Offset]', xtitle = 'Wavelength ('+cgsymbol('Angstrom')+')', title = 'Io''s OI 10.7-9.1 eV' ,  $
      yr = [0.6e4, 6e4], xr = view, /nodata, pos = pos

    x = WL[include_WLs]
    colors = reverse(BytScl(x, MIN=min(x), MAX=2.*max(x)))+128
    
    ;window, 0, xs = 800, ys = 1000
    ;cgplot, WL, WL, xr = [6296,6302], yr = [0.9e4, 8.e4], ytitle = 'R/A', /nodata
    ;cgplot, WL, spec, xr = [5888,5898], yr = [1.e3, 1.e5], ytitle = 'R/A', /nodata
    ;cgplot, WL, spec, xr = view, yr = [1.e4, 4.5e4], ytitle = 'R/A', /nodata
    ;timecolors = [0,1,2,3,4,5,6,7,8,9,10]
    ;begin plotting postscript scatterfit and spec--------------------------------------------------------------------------------------------------------------
    for i = 1, n_elements(Eclipse_files)-1 do begin
      spec = MRDFITS(reduced_dir+'R_per_A_'+Eclipse_files[i]+'.ec.fits', 0, header, /fscale, /unsigned, /silent)
      spec_err = MRDFITS(reduced_dir+'sig_R_per_A_' +Eclipse_files[i]+'.ec.fits', 0, sig_header, /fscale, /unsigned, /silent)
      WL   = sxpar(header, 'CRVAL1')+findgen(N_elements(spec))*sxpar(header, 'Cdelt1')
      raw_header  = headfits(dir+Eclipse_files[i]+'.fits')
      ;cgplot, WL, spec+i*4.e3, Color=timeColors[i], ytitle = 'R/A', /overplot, thick =3

      X = Jup_Cal_Tell_Corr
      X_err = Jup_Cal_Tell_Corr_err
      Y = spec
      Err_Y = spec_err

      s0 = [10., 1200., 1.2, 1.1]

      p = mpfitfun('Match_Jup_Combo', X[scatter_fitting_ind], Y[scatter_fitting_ind], spec_err[scatter_fitting_ind], s0, /NaN, status=status_smooth, parinfo = parinfo, /verbose)
      scatter_fit = interpolate(gauss_smooth(P[0]*X + P[1], P[2], /EDGE_TRUNCATE), findgen(n_elements(P[0]*X + P[1])) - P[3])

      cgplot, WL[include_WLs], scatter_fit[include_WLs]+ 35000 -(i*4.e3),/noerase, yr = [0.6e4, 6e4], linestyle = 1, xr = view, pos = pos
      cgplot, WL[include_WLs], spec[include_WLs]+ 35000 -(i*4.e3),yr = [0.6e4, 6e4],xr = view,COLOR = timeColors[i], /overplot

      residual = y - scatter_fit
      residual_err = err_y

      O4_Io_frame                       = O4 + O4 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      O5_Io_frame                       = O5 + O5 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      O6_Io_frame                       = O6 + O6 * sxpar(header, 'IO_DOPPL') / cspice_clight()

      WL_array[*, i]              = WL[include_WLs]
      residual_array[*, i]        = residual[include_WLs]
      DopplerShift_array4[i]       = O4_Io_frame
      DopplerShift_array5[i]       = O5_Io_frame
      DopplerShift_array6[i]       = O6_Io_frame
    endfor
    AL_legend, ['Raw Io Data','Jupiter Scatter Fit'], Psym = [0,0],linestyle = [0,1], /right, charsize = 1.5
    cgLoadCT, 33, NColors=8, /reverse
    cgps_close


    pos =  [.17,.2,.98,.9]
    cgPS_Open, filename = Reduced_Dir+'OI_107_91_residual_V1.eps', /ENCAPSULATED, xsize = 14.5, ysize = 7.

    !P.font=1
    device, SET_FONT = 'Helvetica Bold', /TT_FONT
    !p.charsize = 2.5
    loadct, 0
    cgplot, WL[include_WLs], WL[include_WLs], psym = 0, Ytitle = 'Rayleighs (5.1 sq arcsec aperture)[Time Offset]', xtitle = 'Wavelength ('+cgsymbol('Angstrom')+')', title = 'Io''s OI 10.7-9.1 eV' ,  $
      yr = [-2.e3, 3.e3], xr = view, /nodata, pos = pos

    x = WL[include_WLs]
    colors = reverse(BytScl(x, MIN=min(x), MAX=2.*max(x)))+128

    ;cgplot, WL, spec, xr = view, yr = [-2.e3, 7.e3], ytitle = 'Residual Jovian scattered light subtraction R/A', /nodata
    for i = 1, n_elements(Eclipse_files)-1 do begin
      cgplot, WL_array[*, i], residual_array[*, i], /OVERPLOT, COLOR = timeColors[i]
    endfor
    for i = 2, n_elements(Eclipse_files)-1 do begin
      residual_array[*,1] = residual_array[*,1]+ residual_array[*,i]
      cgplot, [DopplerShift_array4[i],DopplerShift_array4[i]], [-10000., 10000], COLOR = timeColors[i], linestyle =1, /overplot
      cgplot, [DopplerShift_array5[i],DopplerShift_array5[i]], [-10000., 10000], COLOR = timeColors[i], linestyle =1, /overplot
      cgplot, [DopplerShift_array6[i],DopplerShift_array6[i]], [-10000., 10000], COLOR = timeColors[i], linestyle =1, /overplot
    endfor
    cgplot, WL_array[*, 1], residual_array[*, 1], /OVERPLOT, COLOR = 'black', thick = 6
    
    AL_legend, ['Reduced Io','Codded Residual'], Psym = [0,0],linestyle = [0,0], thick = [0,6], /right, charsize = 1.5
    cgLoadCT, 33, NColors=8, /reverse
    cgps_close

    stop

    ;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ;--------------------------------------------------------Get the OI triplet 11.0 - 9.5eV Scatter Fitting-------------------------------------------------------------------------
    ;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    view = [8441,8451]
    
    IF (date EQ 'UT180320') THEN BEGIN
      range1 = 0
      range2 = 40
    ENDIF ELSE BEGIN
      range1 = 70
      range2 = 135
    ENDELSE

    ; Correct Jovian Scattered Light
    shift_array   = [2.9, 0, 0, 0, 0, 0, 0]
    Smooth_array  = [0, 0, 0, 0, 0, 0, 0]

    parinfo            = replicate( {value:0., fixed: 0b, limited: [0b,0b], limits: fltarr(2) }, 4)
    parinfo[0].limited = 1b
    parinfo[0].limits  = [0., 1.e2]                                 ; limit differences in multiplier amplitude
    parinfo[1].limited = 1b                                 ; limit additive differences
    parinfo[1].limits  = [-1.e4, 1.e4]
    parinfo[2].limited   = 1b                                 ; fix differences in smooth width
    parinfo[2].limits  = [-1.e2, 1.e2]
    parinfo[3].limited   = 1b                                 ; fix differences in smooth width
    parinfo[3].limits  = [-1.e2, 1.e2]
    Gaussian_parinfo       = replicate({value:0., fixed:0, limited:[0,0], limits:[0.,0.]}, 3)

    include_WLs            = where( ((wl gt 8441.) and (wl lt 8451.)), /NULL)
    ;IF (date EQ 'UT180320') THEN BEGIN
    ;  exclude_WLs            = where( ((wl gt 6299.6) and (wl lt 6300.2)), /NULL))
    ;ENDIF ELSE BEGIN
    ;  exclude_WLs            = where( ((wl gt 6300.4) and (wl lt 6300.85)), /NULL))
    ;END ELSE
    scatter_fitting_ind    = include_WLs ;cgSetDifference(include_WLs, exclude_WLs)
    ;IF (date EQ 'UT180320') THEN BEGIN
      ;LSF_fitting_ind        =  where( ((wl gt 6299.5) and (wl lt 6300.2)), /NULL)
    ;ENDIF ELSE BEGIN
      ;LSF_fitting_ind        =  where( ((wl gt 6300.4) and (wl lt 6300.85)), /NULL)
    ;ENDELSE


    residual_array         = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
    LSF_Fit_array          = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
    WL_array               = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
    Brightness_array       = fltarr(n_elements(Eclipse_files))
    err_Brightness_array   = fltarr(n_elements(Eclipse_files))
    T_Shadow_array         = fltarr(n_elements(Eclipse_files))
    EXPTime_array          = fltarr(n_elements(Eclipse_files))
    DopplerShift_array7     = fltarr(n_elements(Eclipse_files))
    DopplerShift_array8     = fltarr(n_elements(Eclipse_files))
    DopplerShift_array9     = fltarr(n_elements(Eclipse_files))
    
    pos =  [.17,.2,.98,.9]
    cgPS_Open, filename = Reduced_Dir+'OI_110_95_scatterfit_V1.eps', /ENCAPSULATED, xsize = 14.5, ysize = 7.

    !P.font=1
    device, SET_FONT = 'Helvetica Bold', /TT_FONT
    !p.charsize = 2.5
    loadct, 0
    cgplot, WL[include_WLs], WL[include_WLs], psym = 0, Ytitle = 'Rayleighs (5.1 sq arcsec aperture)[Time Offset]', xtitle = 'Wavelength ('+cgsymbol('Angstrom')+')', title = 'Io''s OI 11.0-9.5 eV' ,  $
      yr = [0.6e4, 6e4], xr = view, /nodata, pos = pos

    x = WL[include_WLs]
    colors = reverse(BytScl(x, MIN=min(x), MAX=2.*max(x)))+128
    
    ;window, 0, xs = 800, ys = 1000
    ;cgplot, WL, WL, xr = [6296,6302], yr = [0.9e4, 8.e4], ytitle = 'R/A', /nodata
    ;cgplot, WL, spec, xr = [5888,5898], yr = [1.e3, 1.e5], ytitle = 'R/A', /nodata
    ;cgplot, WL, spec, xr = view, yr = [1.e4, 4.5e4], ytitle = 'R/A', /nodata
    ;begin plotting postscript scatterfit and spec--------------------------------------------------------------------------------------------------------------
    for i = 1, n_elements(Eclipse_files)-1 do begin
      spec = MRDFITS(reduced_dir+'R_per_A_'+Eclipse_files[i]+'.ec.fits', 0, header, /fscale, /unsigned, /silent)
      spec_err = MRDFITS(reduced_dir+'sig_R_per_A_' +Eclipse_files[i]+'.ec.fits', 0, sig_header, /fscale, /unsigned, /silent)
      WL   = sxpar(header, 'CRVAL1')+findgen(N_elements(spec))*sxpar(header, 'Cdelt1')
      raw_header  = headfits(dir+Eclipse_files[i]+'.fits')
      ;cgplot, WL, spec+i*4.e3, Color=timeColors[i], ytitle = 'R/A', /overplot, thick =3

      X = Jup_Cal_Tell_Corr
      X_err = Jup_Cal_Tell_Corr_err
      Y = spec
      Err_Y = spec_err

      s0 = [10., 1200., 1.2, 1.1]

      p = mpfitfun('Match_Jup_Combo', X[scatter_fitting_ind], Y[scatter_fitting_ind], spec_err[scatter_fitting_ind], s0, /NaN, status=status_smooth, parinfo = parinfo, /verbose)
      scatter_fit = interpolate(gauss_smooth(P[0]*X + P[1], P[2], /EDGE_TRUNCATE), findgen(n_elements(P[0]*X + P[1])) - P[3])
      
      cgplot, WL[include_WLs], scatter_fit[include_WLs]+ 35000 -(i*4.e3),/noerase, yr = [0.6e4, 6e4], linestyle = 1, xr = view, pos = pos
      cgplot, WL[include_WLs], spec[include_WLs]+ 35000 -(i*4.e3),yr = [0.6e4, 6e4],xr = view,COLOR = timeColors[i], /overplot
      
      residual = y - scatter_fit
      residual_err = err_y
      
      O7_Io_frame                       = O7 + O7 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      O8_Io_frame                       = O8 + O8 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      O9_Io_frame                       = O9 + O9 * sxpar(header, 'IO_DOPPL') / cspice_clight()
      
      WL_array[*, i]              = WL[include_WLs]
      residual_array[*, i]        = residual[include_WLs]
      DopplerShift_array7[i]       = O7_Io_frame
      DopplerShift_array8[i]       = O8_Io_frame
      DopplerShift_array9[i]       = O9_Io_frame
      endfor
      
      AL_legend, ['Raw Io Data','Jupiter Scatter Fit'], Psym = [0,0],linestyle = [0,1], /right, charsize = 1.5
      cgLoadCT, 33, NColors=8, /reverse
      cgps_close
      pos =  [.17,.2,.98,.9]
      cgPS_Open, filename = Reduced_Dir+'OI_110_95_residual_V1.eps', /ENCAPSULATED, xsize = 14.5, ysize = 7.

      !P.font=1
      device, SET_FONT = 'Helvetica Bold', /TT_FONT
      !p.charsize = 2.5
      loadct, 0
      cgplot, WL[include_WLs], WL[include_WLs], psym = 0, Ytitle = 'Rayleighs (5.1 sq arcsec aperture)[Time Offset]', xtitle = 'Wavelength ('+cgsymbol('Angstrom')+')', title = 'Io''s OI 10.7-9.1 eV' ,  $
       yr = [-2.e3, 5.e3], xr = view, /nodata, pos = pos

      x = WL[include_WLs]
      colors = reverse(BytScl(x, MIN=min(x), MAX=2.*max(x)))+128
      ;cgplot, WL, spec, xr = view, yr = [-2.e3, 7.e3], ytitle = 'Residual Jovian scattered light subtraction R/A', /nodata
      for i = 1, n_elements(Eclipse_files)-1 do begin
        cgplot, WL_array[*, i], residual_array[*, i], /OVERPLOT, COLOR = timeColors[i]
      endfor
      for i = 2, n_elements(Eclipse_files)-1 do begin
        residual_array[*,1] = residual_array[*,1]+ residual_array[*,i]
        cgplot, [DopplerShift_array7[i],DopplerShift_array7[i]], [-10000., 10000], COLOR = timeColors[i], linestyle =1, /overplot
        cgplot, [DopplerShift_array8[i],DopplerShift_array8[i]], [-10000., 10000], COLOR = timeColors[i], linestyle =1, /overplot
        cgplot, [DopplerShift_array9[i],DopplerShift_array9[i]], [-10000., 10000], COLOR = timeColors[i], linestyle =1, /overplot
      endfor
        cgplot, WL_array[*, 1], residual_array[*, 1], /OVERPLOT, COLOR = 'black', thick = 6
        
        AL_legend, ['Reduced Io','Codded Residual'], Psym = [0,0],linestyle = [0,0], thick = [0,6], /right, charsize = 1.5
        cgLoadCT, 33, NColors=8, /reverse
        cgps_close
        
stop

;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
;--------------------------------------------------------Get the S I triplet 7.9 - 6.5eV Scatter Fitting-------------------------------------------------------------------------
;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

view = [9210,9240]

IF (date EQ 'UT180320') THEN BEGIN
  range1 = 0
  range2 = 40
ENDIF ELSE BEGIN
  range1 = 70
  range2 = 135
ENDELSE

; Correct Jovian Scattered Light
shift_array   = [2.9, 0, 0, 0, 0, 0, 0]
Smooth_array  = [0, 0, 0, 0, 0, 0, 0]

parinfo            = replicate( {value:0., fixed: 0b, limited: [0b,0b], limits: fltarr(2) }, 4)
parinfo[0].limited = 1b
parinfo[0].limits  = [0., 1.e2]                                 ; limit differences in multiplier amplitude
parinfo[1].limited = 1b                                 ; limit additive differences
parinfo[1].limits  = [-1.e4, 1.e4]
parinfo[2].limited   = 1b                                 ; fix differences in smooth width
parinfo[2].limits  = [-1.e2, 1.e2]
parinfo[3].limited   = 1b                                 ; fix differences in smooth width
parinfo[3].limits  = [-1.e2, 1.e2]
Gaussian_parinfo       = replicate({value:0., fixed:0, limited:[0,0], limits:[0.,0.]}, 3)

include_WLs            = where( ((wl gt 8441.) and (wl lt 8451.)), /NULL)
;IF (date EQ 'UT180320') THEN BEGIN
;  exclude_WLs            = where( ((wl gt 6299.6) and (wl lt 6300.2)), /NULL))
;ENDIF ELSE BEGIN
;  exclude_WLs            = where( ((wl gt 6300.4) and (wl lt 6300.85)), /NULL))
;END ELSE
scatter_fitting_ind    = include_WLs ;cgSetDifference(include_WLs, exclude_WLs)
;IF (date EQ 'UT180320') THEN BEGIN
;LSF_fitting_ind        =  where( ((wl gt 6299.5) and (wl lt 6300.2)), /NULL)
;ENDIF ELSE BEGIN
;LSF_fitting_ind        =  where( ((wl gt 6300.4) and (wl lt 6300.85)), /NULL)
;ENDELSE


residual_array         = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
LSF_Fit_array          = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
WL_array               = fltarr(N_elements(include_WLs), n_elements(Eclipse_files))
Brightness_array       = fltarr(n_elements(Eclipse_files))
err_Brightness_array   = fltarr(n_elements(Eclipse_files))
T_Shadow_array         = fltarr(n_elements(Eclipse_files))
EXPTime_array          = fltarr(n_elements(Eclipse_files))
DopplerShift_array1     = fltarr(n_elements(Eclipse_files))
DopplerShift_array2     = fltarr(n_elements(Eclipse_files))
DopplerShift_array3     = fltarr(n_elements(Eclipse_files))
window, 0, xs = 800, ys = 1000
;cgplot, WL, WL, xr = [6296,6302], yr = [0.9e4, 8.e4], ytitle = 'R/A', /nodata
;cgplot, WL, spec, xr = [5888,5898], yr = [1.e3, 1.e5], ytitle = 'R/A', /nodata
cgplot, WL, spec, xr = view, yr = [1.e4, 4.5e4], ytitle = 'R/A', /nodata
;begin plotting postscript scatterfit and spec--------------------------------------------------------------------------------------------------------------
for i = 1, n_elements(Eclipse_files)-1 do begin
  spec = MRDFITS(reduced_dir+'R_per_A_'+Eclipse_files[i]+'.ec.fits', 0, header, /fscale, /unsigned, /silent)
  spec_err = MRDFITS(reduced_dir+'sig_R_per_A_' +Eclipse_files[i]+'.ec.fits', 0, sig_header, /fscale, /unsigned, /silent)
  WL   = sxpar(header, 'CRVAL1')+findgen(N_elements(spec))*sxpar(header, 'Cdelt1')
  raw_header  = headfits(dir+Eclipse_files[i]+'.fits')
  cgplot, WL, spec+i*4.e3, Color=timeColors[i], ytitle = 'R/A', /overplot, thick =3

  X = Jup_Cal_Tell_Corr
  X_err = Jup_Cal_Tell_Corr_err
  Y = spec
  Err_Y = spec_err

  s0 = [10., 1200., 1.2, 1.1]

  p = mpfitfun('Match_Jup_Combo', X[scatter_fitting_ind], Y[scatter_fitting_ind], spec_err[scatter_fitting_ind], s0, /NaN, status=status_smooth, parinfo = parinfo, /verbose)
  scatter_fit = interpolate(gauss_smooth(P[0]*X + P[1], P[2], /EDGE_TRUNCATE), findgen(n_elements(P[0]*X + P[1])) - P[3])

  cgplot, WL, scatter_fit + i*4.e3, /OVERPLOT, COLOR = timeColors[i], linestyle =3 ; offset each plot

  residual = y - scatter_fit
  residual_err = err_y

  S1_Io_frame                       = S1 + S1 * sxpar(header, 'IO_DOPPL') / cspice_clight()
  S2_Io_frame                       = S2 + S2 * sxpar(header, 'IO_DOPPL') / cspice_clight()
  S3_Io_frame                       = S3 + S3 * sxpar(header, 'IO_DOPPL') / cspice_clight()

  WL_array[*, i]              = WL[include_WLs]
  residual_array[*, i]        = residual[include_WLs]
  DopplerShift_array1[i]       = S1_Io_frame
  DopplerShift_array2[i]       = S2_Io_frame
  DopplerShift_array3[i]       = S3_Io_frame
endfor

window, 1, xs = 800, ys = 1000
cgplot, WL, spec, xr = view, yr = [-2.e3, 7.e3], ytitle = 'Residual Jovian scattered light subtraction R/A', /nodata
for i = 1, n_elements(Eclipse_files)-1 do begin
  cgplot, WL_array[*, i], residual_array[*, i], /OVERPLOT, COLOR = timeColors[i]
endfor
for i = 2, n_elements(Eclipse_files)-1 do begin
  residual_array[*,1] = residual_array[*,1]+ residual_array[*,i]
  cgplot, [DopplerShift_array1[i],DopplerShift_array1[i]], [-10000., 10000], COLOR = timeColors[i], linestyle =1, /overplot
  cgplot, [DopplerShift_array2[i],DopplerShift_array2[i]], [-10000., 10000], COLOR = timeColors[i], linestyle =1, /overplot
  cgplot, [DopplerShift_array3[i],DopplerShift_array3[i]], [-10000., 10000], COLOR = timeColors[i], linestyle =1, /overplot
endfor
cgplot, WL_array[*, 1], residual_array[*, 1], /OVERPLOT, COLOR = 'black'

stop



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