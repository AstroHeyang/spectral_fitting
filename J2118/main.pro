;+
; NAME:
;   main
;
; AUTHOR:
;   He-Yang Liu
;   My Github homepage:
;     https://github.com/AstroHeyang
;
; PURPOSE:
;	The main program of this project.
;
; CALLING SEQUENCE:
;   main
;
; DESCRIPTION:
;   NULL
;
; INPUT:
;   NULL
;
; KEYWORD PARAMETERS:
;   NULL
;
; OUTPUT:
;	A series of images showing the fitting results 
;-

PRO main
  ;---------------- SDSS spectrum -----------------
  ; TODO: generate the filename automatically
  ; 1G1G for Hbeta region
  index_fit = [4755, 5045]
  spec_fit, index_fit, hb_func='2g', /sdss
  plot_result, index_fit=index_fit, hb_func='2g', /SDSS, filename='J2118_SDSS_4755_5045_hb_1g1g.ps'
  
  ; 1G1L for Hbeta region
  spec_fit, index_fit, hb_func='1g1l', /sdss
  plot_result, index_fit=index_fit, hb_func='1g1l', /SDSS, filename='J2118_SDSS_4755_5045_hb_1g1l.ps'
  
  ; 1G2G for Hbeta region
  spec_fit, index_fit, hb_func='3g', /sdss
  plot_result, index_fit=index_fit, hb_func='3g', /SDSS, filename='J2118_SDSS_4755_5045_hb_1g2g.ps'
  
  ;---------------- ISIS spectrum -----------------
  ; 1G1G for Hbeta region
  index_fit = [4755, 5045]
  spec_fit, index_fit, hb_func='2g'
  plot_result, index_fit=index_fit, hb_func='2g', filename='J2118_ISIS_4755_5045_hb_1g1g.ps'
  
  ; 1G1L for Hbeta region
  spec_fit, index_fit, hb_func='1g1l'
  plot_result, index_fit=index_fit, hb_func='1g1l', filename='J2118_ISIS_4755_5045_hb_1g1l.ps'
  
  ; 1G2G for Hbeta region
  spec_fit, index_fit, hb_func='3g'
  plot_result, index_fit=index_fit, hb_func='3g', filename='J2118_ISIS_4755_5045_hb_1g2g.ps'

  ; 1G2G for Hbeta region, fixed Narrow line width
  nlfixed_array = lindgen(7)*50 + 400
  for i=0l,n_elements(nlfixed_array)-1 do begin
    nlfixed = nlfixed_array[i]
	filename = 'J2118_ISIS_' + strtrim(string(index_fit[0]), 2) + '_' + strtrim(string(index_fit[1]), 2) + $
	  '_hb_1g2g_' + 'nlfixedto' + strtrim(string(nlfixed), 2) + '.ps'
    spec_fit, index_fit, hb_func='3g', nlfixed=nlfixed
    plot_result, index_fit=index_fit, hb_func='3g', filename=filename, /nlfixed
  endfor

  ; 1G2G for Hbeta region, fixed one Broad line width
  blfixed_array = lindgen(13)*100 + 1300
  for i=0l,n_elements(blfixed_array)-1 do begin
    blfixed = blfixed_array[i]
	filename = 'J2118_ISIS_' + strtrim(string(index_fit[0]), 2) + '_' + strtrim(string(index_fit[1]), 2) + $
	  '_hb_1g2g_' + 'blfixedto' + strtrim(string(blfixed), 2) + '.ps'
    spec_fit, index_fit, hb_func='3g', blfixed=blfixed
    plot_result, index_fit=index_fit, hb_func='3g', filename=filename, /blfixed
  endfor
END
