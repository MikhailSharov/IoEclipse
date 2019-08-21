Pro ARCES_redux

;Written by C. Schmidt, BU Center for Space Physics, 2018, out of his personal distaste for IRAF
 
;setup P3D for use with the extraction step
defsysv,'!p3d_path','C:\IDL\Io\Apache_Point_Programs\p3d-2.2.6\'
setenv, 'p3d_path=C:\IDL\Io\Apache_Point_Programs\p3d-2.2.6\'    
Parfile = 'C:\IDL\Io\Apache_Point_Programs\ARCES_Redux\parfile.txt' ;the parameter text file that that tells P3D how to work 
Dir = 'D:\DATA\Apache Point Data\UT180320\'


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
  CSPICE_FURNSH, STRCOMPRESS('C:\SPICE\Jupiter_System\vgr1_jup230.bsp')
  cspice_ktotal, 'all', count
  Print, 'Loaded ', strtrim(string(count),2), ' new Spice kernels'


flats = ['Flat.0014.fits', 'Flat.0015.fits', 'Flat.0016.fits']
bias = ['Bias.0020.fits']

bias = MRDFITS(dir+bias[0], 0, header, /fscale, /silent, /unsigned )
s=size(bias)
sx=s(1)
sy=s(2)
window, 0, xs=1800, ys=1250

eclipsed = ['Io_penumbra.0001.fits','Io_eclipsed.0001.fits', 'Io_eclipsed.0002.fits', 'Io_eclipsed.0003.fits', 'Io_eclipsed.0004.fits', 'Io_eclipsed.0005.fits', 'Io_eclipsed.0006.fits']
bigarray = fltarr(n_elements(eclipsed[1:3]),sx,sy)
for i = 0, n_elements(bigarray[*,0,0]) - 1 do begin
  bigarray[i,*,*] = MRDFITS(dir+eclipsed[i], 0, header, /fscale, /silent, /unsigned ) - bias  
endfor  
bigarray[where(bigarray gt 2.^15)] = !values.F_NaN
ecl_mean = Median(bigarray, DIMENSION=1, /even)
;ecl_mean = Min(bigarray, DIMENSION=1, /nan)
tv, bytscl(ecl_mean, 0, 400.)
;wait, 3
Jupiter = ['Jovian_Scatter.0001.fits', 'Jovian_Scatter.0002.fits', 'Jovian_Scatter.0003.fits']
bigarray = fltarr(n_elements(Jupiter), sx, sy)
for i = 0, n_elements(bigarray[*,0,0]) - 1 do begin
  bigarray[i,*,*] = MRDFITS(dir+Jupiter[i], 0, header, /fscale, /silent, /unsigned ) - bias  
endfor  
bigarray[where(bigarray gt 2.^15)] = !values.F_NaN
Jup_mean = MEDIAN(bigarray, DIMENSION=1, /even)
;Jup_mean = Min(bigarray, DIMENSION=1, /nan)
;cghistoplot, test, min_val=00, max_val = 1000, /nan, mininput = -100, maxinput = 100

test = ecl_mean - 1.11*(smart_shift(jup_mean, 0., -0.2))
test = median(test, 3, /even, dimension=1)
test = smooth(test, [1, 5], /nan)
tv, bytscl(test, -5, 20.)


;window, 1, xs=1800, ys=900


;6300 indicies
  xind = [1028,1062]
  yind = [626 , 632]
  disp = 0.0803

ecl1 = float(MRDFITS(dir+eclipsed[1], 0, header, /fscale, /silent, /unsigned )) - float(bias)  
ecl1 = total(ecl1[ xind[0]:xind[1], yind[0]:yind[1] ], 2)
timeColors = BytScl(findgen(7), Min=0, Max=7, Top=7)
Color=timeColors[0]
cgLoadCT, 33, NColors=8
x = indgen(xind[1]-xind[0])
ind = x[where((x lt mean(x)-3) or (x gt mean(x)+3))]

Io_doppler = 6300.3 * (-17.8652767 / 3.e5)
WL = disp*(x - mean(x)) + 6300.304 + Io_Doppler


cgPS_Open, filename = strcompress('C:\IDL\Io\APO_Io_Eclipse_Gray_redux.eps', /remove_all), /ENCAPSULATED, xsize = 8.5, ysize = 6  
        !P.font=1
        device, SET_FONT = 'Helvetica Bold', /TT_FONT
        !p.charsize = 1.7    
        
      cgplot, WL, reverse(ecl1 - mean(ecl1[ind])), /ynozero, Color=timeColors[1], thick = 6, Xtitle  = 'Wavelength ' + cgsymbol('Angstrom'), ytitle = 'DN (Normalized to Zero Continuum)', $
        title = 'Apache Point ARCES Io eclipse: March 20 2018 10:53:22-11:26:16', xstyle = 1
      ecl2 = float(MRDFITS(dir+eclipsed[2], 0, header, /fscale, /silent, /unsigned )) - float(bias) 
      ecl2 = total(ecl2[ xind[0]:xind[1], yind[0]:yind[1] ], 2) 
      cgplot, WL, reverse(ecl2 - mean(ecl2[ind])), /overplot, Color=timeColors[2], thick = 6
      ecl3 = float(MRDFITS(dir+eclipsed[3], 0, header, /fscale, /silent, /unsigned )) - float(bias)  
      ecl3 = total(ecl3[ xind[0]:xind[1], yind[0]:yind[1] ], 2)
      cgplot, WL, reverse(ecl3 - mean(ecl3[ind])), /overplot, Color=timeColors[3], thick = 6
      ecl4 = float(MRDFITS(dir+eclipsed[4], 0, header, /fscale, /silent, /unsigned )) - float(bias)  
      ecl4 = total(ecl4[ xind[0]:xind[1], yind[0]:yind[1] ], 2)
      cgplot, WL, reverse(ecl4 - mean(ecl4[ind])), /overplot, Color=timeColors[4], thick = 6
      ecl5 = float(MRDFITS(dir+eclipsed[5], 0, header, /fscale, /silent, /unsigned )) - float(bias)  
      ecl5 = total(ecl5[ xind[0]:xind[1], yind[0]:yind[1] ], 2)
      cgplot, WL, reverse(ecl5 - mean(ecl5[ind])), /overplot, Color=timeColors[5], thick = 6
      ecl6 = float(MRDFITS(dir+eclipsed[6], 0, header, /fscale, /silent, /unsigned )) - float(bias)  
      ecl6 = total(ecl6[ xind[0]:xind[1], yind[0]:yind[1] ], 2)
      cgplot, WL, reverse(ecl6 - mean(ecl6[ind])), /overplot, Color=timeColors[6], thick = 6
cgPS_Close

;-------------------------------------Gray's redux-------------------------------------------
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
line = 4900.0  

dir  = 'H:\DATA\Apache Point Data\UT180320\Candace_redux\denebola-data\candaceg\Venus_data_raw\reduced\UT180320_Io\'
spec = MRDFITS(dir+'fullspecIo_eclipsed.0001.ec.fits', 0, header, /fscale, /silent, /unsigned )
WL   = sxpar(header, 'CRVAL1')+findgen(170105)*sxpar(header, 'Cdelt1')

wl_range      = 1600. ;angstroms
shift_array   = [0, 1.5, 1.75, 2., 2.25, 2.5, 2.75]
;Smooth_array  = [0,  18,  18,  18,  18,  18,  18] ; Good for H alpha 6563A
;Smooth_array  = [0,  14,  14,  14,  14,  14,  14] ; Good for Fe I    5270A
Smooth_array  = [0,  14,  14,  14,  14,  14,  14]  ; Good for 5893

window, 0, xs = 1000, ys = 800
pos = cgLayout([1,2], YGap=0)
jup = MRDFITS(dir+'fullspecJupiter_Center_Spectrum.0001.ec.fits', 0, header, /silent, /fscale, /unsigned )
junk = [min(abs(wl - (line-wl_range/2)), low_ind), min(abs(wl - (line+wl_range/2)), hi_ind)]
cgplot, wl, spec, xr = [Line-wl_range/2, line+wl_range/2], /nodata, pos = pos[*,0], yr = [-.3, 1.2], xtickformat = '(A1)'
for i = 1, 6 do begin
  spec = MRDFITS(dir+'fullspecIo_eclipsed.000'+strcompress(i, /remove_all)+'.ec.fits', 0, header, /silent, /fscale, /unsigned )
  cgplot, wl, spec-.15*i, /overplot, Color=timeColors[i], thick = 2, pos = pos[*,0] 
  cgplot, wl, smart_shift(smooth(jup, smooth_array[i]), shift_array[i], /Interp)-.15*i, /overplot, Color=timeColors[i], thick = 2, pos = pos[*,0]  
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
  residual = spec - smart_shift(smooth(jup, smooth_array[i]), shift_array[i], /Interp)
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

cgplot, wl, co_add, thick = 2, /overplot

print, tot_stat


stop
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

;--------------------------------------------------------------------------------------------

;ecl0 = float(MRDFITS(dir+eclipsed[0], 0, header, /fscale, /silent, /unsigned )) - float(bias)  
;ecl0 = total(ecl0[1038:1052,625:633], 2)
;cgplot, ecl0 - min(ecl0), /ynozero, psym = 10, /overplot, Color=timeColors[0]

stop


    
    bigarray = fltarr(n_elements(flats),sx,sy)
    for i = 0, n_elements(flats) - 1 do begin
      bigarray[i,*,*] = MRDFITS(dir+flats[i], 0, header, /fscale, /silent, /unsigned );- bias
    endfor  
    ;Take the pixel-by-pixel median of the bias images in the list
    flat = MEDIAN(bigarray, DIMENSION=1, /even) - bias
       

    ;Locate the full white-light spectra of each and every fiber in the dome flat spectral images, use in fiber extraction program
    
  img = flat
  img = img[200:1850,1:2048] 
;===================================Trace the flat and find all the gaps between the spectra=================================================

s=size(img)
sx=s(1)
sy=s(2)

window, 0, xs = sx, ys = 1000
tv, bytscl(img[*,0:1000], 0, 800)


;Specify x an positions, evenly spaced about the center of the chip, to start looking for the orders 

;Search_x_positions = [200, 500, 800, 1100, 1400]
Search_x_positions = [800, 900, 700]
Separation = 10.

dummy = img
samples = []
For n_x_pos = 0, n_elements(Search_x_positions)-1 do begin
  p3d_tracing_findspec0,  img,  Search_x_positions[n_x_pos],  10.,  .8,  ypos,  10.,  /verbose  
  s = size(samples, /dim)
      print, ypos[50]
      ;did it find all the same traces? if not replace with NaN 
      if n_x_pos eq 1 then ref = ypos[50]
      if n_x_pos ge 1 then begin
        ;if s[1] ne n_elements(ypos) then stop
        align = round((ypos[50] - ref) / Separation)
        ypos = shift(ypos, align)
        if n_elements(ypos) gt s[1] then ypos = ypos[0:s[1]-1]
;        for j = 0, s[1], do begin
;          diff = ypos[j] - samples[j]
;        endfor  
      endif
  samples = [samples, transpose(ypos)]
  for i = 0, n_elements(ypos)-1 do dummy[Search_x_positions[n_x_pos]-10:Search_x_positions[n_x_pos]+10, ypos[i]] = 0.
  tv,bytscl(dummy[*,200:1200], 100, 800)
endfor
stop


;--------------------------------------center-------------------------------------------------------
xpos = 200
p3d_tracing_findspec0,  img,  xpos,  10.,  .8,  ypos,  10.,  /verbose
dummy = img
for i = 0, n_elements(ypos)-1 do dummy[xpos-10:xpos+10, ypos[i]] = 0.
tv,bytscl(dummy[*,200:1200], 100, 800)
;oldpos=ypos
;dummy = img
;p3d_tracing_correctpos, img,  oldpos,  -1,  10.,  5.4,  8,  9,  newpos,  integral, /verbose
;for i = 0, n_elements(newpos)-1 do dummy[xpos-10:xpos+10,newpos[i]] = 0.
;tv,bytscl(dummy[*,200:1200], 100, 800)
;
;stop
;--------------------------------------center-------------------------------------------------------

xpos = 800
p3d_tracing_findspec0,  img,  xpos,  10.,  .8,  ypos,  10.,  /verbose
;dummy = img
for i = 0, n_elements(ypos)-1 do dummy[xpos-10:xpos+10, ypos[i]] = 0.
tv,bytscl(dummy[*,200:1200], 100, 800)

;oldpos=ypos
;dummy = img
;p3d_tracing_correctpos, img,  oldpos,  -1,  10.,  5.4,  8,  9,  newpos,  integral, /verbose
;for i = 0, n_elements(pos)-1 do dummy[xpos-10:xpos+10,newpos[i]] = 0.
;tv,bytscl(dummy[*,200:1200], 100, 800)

xpos = 1400
p3d_tracing_findspec0,  img,  xpos,  10.,  .8,  ypos,  10.,  /verbose
;dummy = img
for i = 0, n_elements(ypos)-1 do dummy[xpos-10:xpos+10, ypos[i]] = 0.
tv,bytscl(dummy[*,200:1200], 100, 800)

stop





p3d_tracing_trace, parfile, img, aperture_trace, daxis = 2, spec0 = 540., var_spec0 = 2, /no_display, /verbose
;p3d_tracing_trace, parfile, img, aperture_trace, daxis = 2, spec0 = 282.0, var_spec0 = 2, /no_display, /verbose
;stop

dummy = img
  for wavelength_row = 0, sx-1 do begin
    ;dummy[ aperture_trace[*, wavelength_row], wavelength_row] = 9.e99
    dummy[ wavelength_row, aperture_trace[*, wavelength_row]] = 0;9.e99
  endfor 
  window,1, xs = sx, ys = 1000, title = 'Aperature Trace of the Flat Field' 
  tv,bytscl(dummy[*,0:1000], 100, 800)
stop
;window, 1
;plot, img[900,400:450]


    p3d_tracing_calculate_lprofs, 'D:\DATA\Apache Point Data\UT180320\Flat.0016.fits', 'D:\DATA\Apache Point Data\UT180320\Bias.0020.fits', $
      Aperture_trace, lprofs, parfile = parfile, kwrdlist = parfile, dmbias = 13., gain = .7, rdnoise = 13., daxis = 1, postscriptfile = strcompress(dir+'Profile_fits.ps'), $
      monitor = 0, blockgui=1, verbose = 2


stop



end