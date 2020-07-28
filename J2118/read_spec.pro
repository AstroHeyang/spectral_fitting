;+
; NAME:
;   read_spec
;
; AUTHOR:
;   He-Yang Liu
;   My Github homepage:
;     https://github.com/AstroHeyang
;
; PURPOSE:
;   Read the spectrum and then return wave and flux after correction (both deredden and redshift) 
;
; CALLING SEQUENCE:
;   spec = read_spec([spec_file, z, ra, dec, /sdss])
;
; DESCRIPTION:
;   read the fits file 
;
; INPUT:
;   spec_file - ISIS/SDSS spectrum
;
; KEYWORD PARAMETERS:
;   z - redshifit
;   ra - right ascension
;   dec - declination
;
; OUTPUT:
;   struct, {wave, flux, flux_error, z, ra, dec} 
;-

function read_spec, spec_file, spec_z, ra, dec, sdss=sdss

  if not keyword_set(sdss) then begin
    if not keyword_set(spec_file) then begin
      spec_file = './spec/wht20180813_02705173_spew.fits'
    endif 
    spec = mrdfits(spec_file, 0, h)
    flux = spec[*, 0, 0]
    w0 = sxpar(h, 'CRVAL1')
    dw = sxpar(h, 'CD1_1')
    naxis1 = sxpar(h, 'NAXIS1')
    sigma = spec[*, 0, 3]
    wave_raw = w0 + dindgen(naxis1)*dw
    index_good = where(wave_raw gt 5600 and wave_raw lt 9000 and sigma gt 0)
    spec_wave = wave_raw[index_good]
    spec_flux = flux[index_good]/200.0
    naxis1 = n_elements(index_good)
    spec_error = sigma[index_good]/200.0
  endif

  ; sdss case
  if keyword_set(sdss) then begin
    if not keyword_set(spec_file) then begin
      spec_file = './spec/spec-0639-52146-0179.fits'
    endif
    spec = mrdfits(spec_file, 1, /silent)
    wave_raw = 10.0^(spec.loglam)
    flux = spec.flux
    index_good = where(spec.ivar gt 0 and (spec.and_mask eq 0 or  spec.and_mask eq 2l^24))
    spec_wave = wave_raw[index_good]
    spec_ivar = spec[index_good].ivar
    spec_error = spec_ivar^(-0.5)
    spec_flux = flux[index_good]
  endif
  
  if not keyword_set(spec_z) then $
    spec_z = 0.26011
  
  if not keyword_set(ra) then begin
    spec_ra = 319.720667
  endif else begin
    spec_ra = ra
  endelse
  
  if not keyword_set(dec) then begin
    spec_dec = -7.54097
  endif else begin
    spec_dec = dec
  endelse
  
  ; correct the Galactic extinction
  glactc, spec_ra, spec_dec, 2000, gal_l, gal_b, 1, /degree
  gal_ebv = 0.86*dust_getval(gal_l, gal_b, ipath='../qso_fit/SFD/maps/',/interp)
  ; Schlegel et al. 1998， Schlafly & Finkbeiner 2011
  print, 'Galactic E(B-V): ', gal_ebv, 'mag', format='(A18, D10.5, A4)'
  obs_wave = spec_wave
  fm_unred, obs_wave, spec_flux, gal_ebv, flux_unred
  fm_unred, obs_wave, spec_error, gal_ebv, error_unred
  wave = obs_wave/(1.0+spec_z)
  
  res_struct = {wave: wave, flux: flux_unred, flux_error: error_unred, z:spec_z, ra:spec_ra, dec:spec_dec}
  res = replicate(res_struct, 1)
  return, res
end

