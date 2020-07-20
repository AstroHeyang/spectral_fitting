;+
; NAME:
;	lorentz1
;
; AUTHOR:
;	He-Yang Liu
;	My Github homepage:
;		https://github.com/AstroHeyang
;
; PURPOSE:
;	Get a Lorentz model
;
; CALLING SEQUENCE:
;	model = lorentz1(x, p)
;
; DESCRIPTION:
;	The Lorentz model is used to fit the broad line profile, especially for NLS1s.
;	L(x) = 1/pi * (0.5 * Gamma) / [(x - x0)^2 + (0.5 * Gamma)^2],
;   where Gamma is the FWHM.
;
; INPUT:
;	x - an array of wavelength in log-scale
;
;	p - an array including 3 parameters, [w0, FWHM/2.35482, Flux] 
;
; KEYWORD PARAMETERS:
;	NULL
;
; OUTPUT:
;	return the array of flux in shape of Lorentz1
;
; Version:
;	v1.1 xbdong 2005
;-

function lorentz1, x, p

	np = n_elements(p)
	if np lt 3 then begin
		print, "ERROR! lorentz1: n_elements(p) lt 3!"
		return, !VALUES.F_NAN
	endif

	p[1] = p[1] * 2.35482D
	;; Prevent floating overflow, xbdong
	u = (4.d * (x - p[0])^2 + p[1]^2 ) > 1e-20  
	;; Prevent floating underflow, xbdong
    mask = u lt 1.D300  
	y = 2.D * p[2] * mask / 3.1415926D * p[1] / (u * mask )

    return, y
end
