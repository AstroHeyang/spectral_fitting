;+
; NAME:
;	nlFeIIfunc
;
; AUTHOR:
;	He-Yang Liu
;	My Github homepage:
;		https://github.com/AstroHeyang
;
; PURPOSE:
;	Get the narrow-line Fe II multiplets
;
; CALLING SEQUENCE:
;	nlFeII = nlFeIIfunc(x, p)
;
; DESCRIPTION:
;	The Fe II multiplets are obtained from the well-known NLS1 Izw~1 (Veron-City et al. 2004) 
;	Also see Dong et al. (2011) for details.
;
; INPUT:
;	x - an array of rest-wavelength in log-scale
;
;	p - an array including 3 parameters of Fe II model; 
;		[log_w0, sigma, flux] for gauss1 function;
;
; KEYWORD PARAMETERS:
;	forbidden - only forbidden Fe II lines are used
;
;	permitted - only permitted Fe II lines are used
;
;	extended - forbidden and permitted Fe II lines and other lines 
;			   including Cr, Ni, Si, Ti, Ca, are used

; OUTPUT:
;	return the array of flux for narrow-line Fe II
; 
; VERSIONS:
;	v1 - xbdong 2005 
;-

function nlFeIIfunc, x, p, FeII_file=FeII_file, forbidden=forbidden, $
         permitted=permitted, extended=extended

	On_error, 2

	if n_elements(p) lt 3 then begin
		print, "ERROR! nlFeII_func: n_elements(p) lt 3. ERROR!"
		return, !VALUES.F_NAN
	endif

	if not keyword_set(FeII_file) then begin
		FeII_file = 'VJV04_narrow_broad_FeII_data.fits.gz'
	endif
	nlFeII_VJV04_DATA = mrdfits(FeII_file, 1)

	key_tag = keyword_set(forbidden) + keyword_set(permitted) + keyword_set(extended)
	if key_tag gt 1 then begin
		print, 'forbidden, permitted, extended: more than 1 keywords are set.'
		stop
	endif

	T_type = 'default'
	if keyword_set(extended) then begin
		T_type = 'extended'
	endif else begin
		if keyword_set(forbidden) then begin
			T_type = 'forbidden'
		endif else begin
			if keyword_set(permitted) then T_type = 'permitted'
		endelse
	endelse

	case T_type of
		'extended': $
			this_index = where(	( strcmp(nlFeII_VJV04_DATA.line, "[Fe II]", 7) eq 1 ) or $
								( strcmp(nlFeII_VJV04_DATA.line, "Fe II", 5)   eq 1 ) or $
                            	( strcmp(nlFeII_VJV04_DATA.line, "[Cr II]", 7) eq 1 ) or $
                            	( strcmp(nlFeII_VJV04_DATA.line, "[Ni II]", 7) eq 1 ) or $
                            	( strcmp(nlFeII_VJV04_DATA.line, "[Ti II]", 7) eq 1 ) or $
                            	( strcmp(nlFeII_VJV04_DATA.line, "Cr II", 5) eq 1 ) or $
                            	( strcmp(nlFeII_VJV04_DATA.line, "Ni II", 5) eq 1 ) or $
                            	( strcmp(nlFeII_VJV04_DATA.line, "Ti II", 5) eq 1 ) or $
                            	( strcmp(nlFeII_VJV04_DATA.line, "Si II", 5) eq 1 ) or $
                            	( strcmp(nlFeII_VJV04_DATA.line, "Ca II", 5) eq 1 ) or $
                            	( strcmp(nlFeII_VJV04_DATA.line, "[Si II]", 7) eq 1 ) or $
                            	( strcmp(nlFeII_VJV04_DATA.line, "[Ca II]", 7) eq 1 ) )
		'forbidden':$
			this_index = where( strcmp(nlFeII_VJV04_DATA.line, "[Fe II]", 7) eq 1)
		'permitted':$
			this_index = where( strcmp(nlFeII_VJV04_DATA.line, "Fe II", 5)   eq 1)
		'default':$
			this_index = where( ( strcmp(nlFeII_VJV04_DATA.line, "[Fe II]", 7) eq 1 ) or $
                            ( strcmp(nlFeII_VJV04_DATA.line, "Fe II", 5)   eq 1 ) )
		else: begin
			print, 'wrong T_type:', T_type
			stop
			end
	endcase

	nlFe_W0s = nlFeII_VJV04_DATA[this_index].air_lambda
	airtovac, nlFe_W0s
	nlFe_intensity = nlFeII_VJV04_DATA[this_index].intensity

	vac_w0_ref = 5019.8397 
	;; FeII42 \lambda 5018. My Convention.
    ;; So the parameter structure of FeII looks like a single line
	nlFe_logW0s = alog(nlFe_W0s)
	dw = p[0] - alog(vac_w0_ref)

	nlFe_flux = x * 0.
	for i = 0, n_elements(nlFe_logW0s)-1 do begin
		pp = [ nlFe_logW0s[i] + dw, p[1], nlFe_intensity[i] ]
		if finite(pp[1]) eq 1 and pp[1] gt 0 then $
			nlFe_flux += gauss1(x, pp)
	endfor

	return, p[2] * nlFe_flux
end 
