;+
; NAME:
;   read_par
;
; AUTHOR:
;   He-Yang Liu
;   My Github homepage:
;       https://github.com/AstroHeyang
;
; PURPOSE:
;   Read the default parameter file
;
; CALLING SEQUENCE:
;   par = read_par(parfile)
;
; DESCRIPTION:
;	NULL
;
; INPUT:
;   parfile - an default parameter file
;
; KEYWORD PARAMETERS:
;   NULL
;
; OUTPUT:
;   return a struct including the default parameters
;-
function read_par, parfile

	if not keyword_set(parfile) then begin
		parfile = 'default.par'
	endif

	par = {	par_base,			$
		OIII_4959:		0,		$
		OIII_5007:		0,		$
		bl_Hbeta_4861:	0,		$
		bl_Halpha_6563:	0,		$
		lorentz:		0		$
		}

	par_name = tag_names(par)
	npar = n_elements(par_name)

	lines = djs_readlines(parfile)
	nlines = n_elements(lines)

	for i=0L, nlines-1 do begin
		this_line = strtrim(lines[i], 2)
		first = strmid(this_line, 0, 1)
		if first NE '%' and first NE '' then begin
			this_par = strsplit(this_line, /extra)
			this_par_name = strupcase(strtrim(this_par[0], 2))
			for j=0L, npar-1 do begin
				if this_par_name EQ par_name[j] then $
					par.(j) = this_par[1]
			endfor
		endif
	endfor

	par.OIII_4959 = par.OIII_5007
	return, par

end
