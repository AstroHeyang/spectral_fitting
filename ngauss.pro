;+
; NAME:
;	ngauss
;
; AUTHOR:
;	He-Yang Liu
;	My Github homepage:
;		https://github.com/AstroHeyang
;
; PURPOSE:
;	Get a multiple Gaussian model
;
; CALLING SEQUENCE:
;	broadline_model = ngauss(x, p)
;
; DESCRIPTION:
;	The multiple Gaussian model is used to fit the broad line profile or narrow lines
;	with strange profiles
;
; INPUT:
;	x - an array of wavelength in log-scale
;
;	p - an array including 3*n parameters 
;
; KEYWORD PARAMETERS:
;	NULL
;
; OUTPUT:
;	return the array of flux in shape of multiple Gaussian 
;-

function ngauss, x, p

	np = n_elements(p)
	if (np mod 3) ne 0 or np lt 3 then begin
		print, "ERROR! ngauss: n_elements(p) mod 3 ne 0!"
		return, !VALUES.F_NAN
	endif

	y = x * 0
	for i = 0,np/3-1 do begin
		y += gauss1(x, [p[i*3], p[i*3 + 1], p[i*3 + 2]])
	endfor

    return, y
end
