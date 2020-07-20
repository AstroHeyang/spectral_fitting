;+
; NAME:
;	blFeIIfunc
;
; AUTHOR:
;	He-Yang Liu
;	My Github homepage:
;		https://github.com/AstroHeyang
;
; PURPOSE:
;	Get the broad-line Fe II multiplets
;
; CALLING SEQUENCE:
;	blFeII = blFeIIfunc(x, p)
;
; DESCRIPTION:
;	The Fe II multiplets are obtained from the well-known NLS1 Izw~1 (Veron-City et al. 2004) 
;	Also see Dong et al. (2011) for details.
;
; INPUT:
;	x - an array of rest-wavelength in log-scale
;
;	p - an array including 3 parameters of Fe II model; 
;		[log_w0, sigma, flux] for lorentz1 function;
;		if g4model is set, sigma is useless
;
; KEYWORD PARAMETERS:
;	extended - Besides FeII, other lines such as Cr, Ni, Ti, Si, Ca, are used
;
; OUTPUT:
;	return the array of flux for broad-line Fe II
; 
; VERSIONS:
;	v1 - xbdong 2005 
;-

function blFeIIfunc, x, p, FeII_file=FeII_file, extended=extended

	On_error, 2

	if n_elements(p) lt 3 then begin
		print, "ERROR! blFeIIfunc: n_elements(p) lt 3."
		return, !VALUES.F_NAN
	endif

	if not keyword_set(FeII_file) then begin
		FeII_file = 'VJV04_narrow_broad_FeII_data.fits.gz'
	endif

	T_type = 'default'
	if keyword_set(extended) then begin
		T_type = 'extended'
	endif 

    blFeII_VJV04_DATA = mrdfits(Feii_file, 2)

	case T_type of
		'extended': $
			this_index = where( ( strcmp(blFeII_VJV04_DATA.line, "Fe II", 5) eq 1 ) or $
								( strcmp(blFeII_VJV04_DATA.line, "Cr II", 5) eq 1 ) or $
                               	( strcmp(blFeII_VJV04_DATA.line, "Ni II", 5) eq 1 ) or $
                               	( strcmp(blFeII_VJV04_DATA.line, "Ti II", 5) eq 1 ) or $
                               	( strcmp(blFeII_VJV04_DATA.line, "Si II", 5) eq 1 ) or $
                               	( strcmp(blFeII_VJV04_DATA.line, "Ca II", 5) eq 1 ) )

         'default': $
			this_index = where(( strcmp(blFeII_VJV04_DATA.line, "Fe II", 5) eq 1 ) )
    endcase

    blFe_W0s = blFeII_VJV04_DATA[this_index].air_lambda
    airtovac, blFe_W0s
    blFe_intensity = blFeII_VJV04_DATA[this_index].intensity

    vac_W0_ref = 5019.8397
    blFe_logW0s = alog(blFe_W0s)
    dw = p[0] - alog(vac_W0_ref)

    blFe_flux = x * 0.

	for i = 0, n_elements(blFe_logW0s)-1 do begin
		pp = [ blFe_logW0s[i] + dw, p[1], blFe_intensity[i] ]
		if finite(pp[1]) eq 1 and pp[1] gt 0 then $
			blFe_flux += lorentz1(x, pp)
	endfor 

	return, p[2] * blFe_flux
end
