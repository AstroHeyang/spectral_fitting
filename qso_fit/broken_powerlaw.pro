;+
; NAME:
;	broken_powerlaw
;
; AUTHOR:
;	He-Yang Liu
;	My Github homepage:
;		https://github.com/AstroHeyang
;
; PURPOSE:
;	Get a broken power-law model with broken wavelength of 5650 A
;
; CALLING SEQUENCE:
;	continuum_model = broken_powerlaw(log_wave, p)
;
; DESCRIPTION:
;	The power-law model is used to represent the AGN continuum.
;
; INPUT:
;	log_wave - an array of wavelength in log-scale
;
;	p - an array including 3 parameters of power law
;
; KEYWORD PARAMETERS:
;	NULL
;
; OUTPUT:
;	return the array of flux in shape of broken power-law 
;-

function broken_powerlaw, x, p

	if n_elements(p) lt 3 then begin
		print, "broken_powlaw: n_elements(p) lt 3. ERROR!"
		return, !VALUES.F_NAN
	endif

	restwave = exp(x)
	y = restwave[*] * 0.d
	index = where(restwave lt 5650)
	if index[0] ne -1 then y[index] = p[0] * (restwave[index] / 5000.d) ^ p[1]
	index = where(restwave ge 5650)
	if index[0] ne -1 then begin
		p2 = p[0] * (5650 / 5000.d) ^ p[1] / (5650 / 5000.d) ^ p[2]
		y[index] = p2 * (restwave[index] / 5000.d) ^ p[2]
	endif

    return,y
end
