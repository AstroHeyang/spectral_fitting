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
  ; 1G1G for Hbeta region
  index_fit = [4755, 5045]
  spec_fit, index_fit, hb_func='2g', /sdss
  plot_result, index_fit=index_fit, hb_func='2g', /SDSS, filename='J2118_SDSS_hb_1g1g.ps'
  
  ; 1G1L for Hbeta region
  spec_fit, index_fit, hb_func='1g1l', /sdss
  plot_result, index_fit=index_fit, hb_func='1g1l', /SDSS, filename='J2118_SDSS_hb_1g1l.ps'
  
  ; 1G2G for Hbeta region
  spec_fit, index_fit, hb_func='3g', /sdss
  plot_result, index_fit=index_fit, hb_func='3g', /SDSS, filename='J2118_SDSS_hb_1g2g.ps'
  
  ;---------------- ISIS spectrum -----------------
  ; 1G1G for Hbeta region
  index_fit = [4755, 5045]
  spec_fit, index_fit, hb_func='2g'
  plot_result, index_fit=index_fit, hb_func='2g', filename='J2118_ISIS_hb_1g1g.ps'
  
  ; 1G1L for Hbeta region
  spec_fit, index_fit, hb_func='1g1l'
  plot_result, index_fit=index_fit, hb_func='1g1l', filename='J2118_ISIS_hb_1g1l.ps'
  
  ; 1G2G for Hbeta region
  spec_fit, index_fit, hb_func='3g'
  plot_result, index_fit=index_fit, hb_func='3g', filename='J2118_ISIS_hb_1g2g.ps'
END
