; This code reproduces the historical energy budget constraint
; on radiative forcing and Figure 1 as presented in Andrews and
; Forster (2019; Nature Climate Change).

PRO idl_af2019_ncc

 binwidth=0.02 ; for histogram
 max_F=10 ; for histogram
 min_F=-10 ; for histogram
 scale=1.645 ; converts between sigma=1 and 5-95% range
 nsamples=1000000

 ; Define normal distributions with a mean of zero
 ; and a standard deviation of one.
 norm1=RANDOMN(SEED1,nsamples)
 norm2=RANDOMN(SEED2,nsamples)
 norm3=RANDOMN(SEED3,nsamples)
 norm4=RANDOMN(SEED4,nsamples)

 ; Define dT, dN, lambda and uncertainty
 dT_mean=0.98 
 dT_5095=0.20
 dN_mean=0.61
 dN_5095=0.35
 lambda_mean=1.74
 lambda_5095=0.48

 ; Define non-aerosol forcing and uncertainty
 dF_nonAER_mean=3.06
 dF_nonAER_5095=0.45

 ; Generate dT, dN, lambda and F_nonAER distributions
 dT=dT_mean+norm1*dT_5095/scale
 dN=dN_mean+norm2*dN_5095/scale
 lambda=lambda_mean+norm3*lambda_5095/scale
 dF_nonAER=dF_nonAER_mean+norm4*dF_nonAER_5095/scale

 ; Calculate dF and dF_AER
 dF=dN+lambda*dT
 dF_AER=dF-dF_nonAER

 ; Calculate forcing PDF and percentiles
 pdf=HISTOGRAM(dF,MIN=min_F,MAX=max_F,BINSIZE=binwidth,LOCATIONS=xbin)
 pdf=pdf/TOTAL(pdf,/DOUBLE)/binwidth
 cdf=TOTAL(pdf,/CUMULATIVE,/DOUBLE)*binwidth
 k=WHERE(cdf GE 0.05)
 l=WHERE(cdf GE 0.95)
 m=WHERE(cdf GE 0.50)
 pdf_dF=pdf
 p50_dF=xbin(m[0])+binwidth/2
 p05_dF=xbin(k[0])+binwidth/2.
 p95_dF=xbin(l[0])+binwidth/2.

 ; Calculate aerosol forcing PDF and percentiles
 pdf=HISTOGRAM(dF_AER,MIN=min_F,MAX=max_F,BINSIZE=binwidth,LOCATIONS=xbin)
 pdf=pdf/TOTAL(pdf,/DOUBLE)/binwidth
 cdf=TOTAL(pdf,/CUMULATIVE,/DOUBLE)*binwidth
 k=WHERE(cdf GE 0.05)
 l=WHERE(cdf GE 0.95)
 m=WHERE(cdf GE 0.50)
 pdf_dF_AER=pdf
 p50_dF_AER=xbin(m[0])+binwidth/2
 p05_dF_AER=xbin(k[0])+binwidth/2.
 p95_dF_AER=xbin(l[0])+binwidth/2.

 ; Define IPCC AR5 forcing and uncertainty
 dF_AR5_med=2.30-0.09
 dF_AR5_p95=3.30-0.09
 dF_AR5_p05=1.10-0.09
 dF_AR5_AER_med=-0.90+0.21
 dF_AR5_AER_p95=-0.10+0.21
 dF_AR5_AER_p05=-1.90+0.21

 ; Define CMIP5 historical forcing and uncertainty 
 dF_CMIP5_med=1.7+0.2
 dF_CMIP5_p95=1.7+0.9+0.2
 dF_CMIP5_p05=1.7-0.9+0.2
 
 ; Define CMIP5 aerosol forcings
 dF_AER_CMIP5_med=-1.17
 dF_AER_CMIP5_p95=-1.17+0.49
 dF_AER_CMIP5_p05=-1.17-0.49

 PRINT,'IPCC AR5 Median, 5-95%     ',dF_AR5_med,  dF_AR5_p05,  dF_AR5_p95,  ' AEROSOL ',dF_AR5_AER_med,  dF_AR5_AER_p05,  dF_AR5_AER_p95
 PRINT,'CMIP5 Median, 5-95%        ',dF_CMIP5_med,dF_CMIP5_p05,dF_CMIP5_p95,' AEROSOL ',dF_AER_CMIP5_med,dF_AER_CMIP5_p05,dF_AER_CMIP5_p95
 PRINT,'dF Constraint Median, 5-95%',p50_dF,      p05_dF,      p95_dF,      ' AEROSOL ',p50_dF_AER,      p05_dF_AER,      p95_dF_AER


 ; ********** PLOT FIGURE 1 *******************
 A=FIndGen(16)*(!PI*2/16.)
 UserSym,cos(A),sin(A),/fill
 !P.MULTI = [0,2,2]
 TEK_COLOR
 !P.FONT = 0
 !P.THICK = 5
 !X.THICK = 5
 !Y.THICK = 5
 !P.CHARTHICK=2
 PR, file='Figure1.ps',/CPS,/LANDSCAPE
   PLOT,xbin,pdf_dF,xrange=[-3.5,4.5],yrange=[0,1.4],xstyle=9,ystyle=9,THICK=5,$
        xtitle='Effective Radiative Forcing (Wm!E-2!N)',$
        ytitle='Probability Distribution Function',/NODATA
 
    ; ********** TOTAL FORCING *******************
   
    ; Plot forcing PDF
    OPLOT,xbin,pdf_dF,THICK=5,COLOR=0
   
    ; Plot 5-95% bar for AR5 forcing
    col = 2
    ypos1 = +1.3
    ypos=ypos1
    ysize=0.04
    p95=dF_AR5_p95
    p05=dF_AR5_p05
    med=dF_AR5_med
    OPLOT, [med, med],[ypos, ypos+ysize],color=col, thick=4
    OPLOT, [p05, p95],[ypos+ysize/2., ypos+ysize/2.],color=col, thick=4  
    OPLOT, [p95,p95],[ypos, ypos+ysize],color=col, thick=4
    OPLOT, [p05,p05],[ypos, ypos+ysize],color=col, thick=4

    ; Plot 5-95% bar for CMIP5 forcing & individual models
    col = 15
    ypos2 = ypos1-0.08
    ypos=ypos2
    ysize=0.06
    p95=dF_CMIP5_p95
    p05=dF_CMIP5_p05
    med=dF_CMIP5_med
    OPLOT, [med, med],[ypos, ypos+ysize],color=col, thick=4
    OPLOT, [p05, p95],[ypos+ysize/2., ypos+ysize/2.],color=col, thick=2
    OPLOT, [p95,p95],[ypos, ypos+ysize],color=col, thick=4
    OPLOT, [p05,p05],[ypos, ypos+ysize],color=col, thick=4
    ; Define individual CMIP5 models
    dF_CMIP5_models=[1.1,2.2,2.2,2.0,2.5,1.5,0.9,2.3,1.1,$
                  2.0,2.0,2.3,2.5,0.8,1.7,1.9,1.0,1.6,$
                  1.1,2.1,2.3,1.2,1.4]+0.2
    FOR j=0,N_ELEMENTS(dF_CMIP5_models)-1 DO BEGIN
      OPLOT,[dF_CMIP5_models[j],dF_CMIP5_models[j]],[ypos+ysize/2.,ypos+ysize/2.],PSYM=8,COLOR=col,SYMSIZE=0.4
    ENDFOR
      ; Offset individual CMIP5 models that sit on top of each other
      OPLOT,[1.3,1.3],[ypos,ypos],PSYM=8,COLOR=col,SYMSIZE=0.4
      OPLOT,[1.3,1.3],[ypos+0.06,ypos+0.06],PSYM=8,COLOR=col,SYMSIZE=0.4
      OPLOT,[2.2,2.2],[ypos,ypos],PSYM=8,COLOR=col,SYMSIZE=0.4
      OPLOT,[2.2,2.2],[ypos+0.06,ypos+0.06],PSYM=8,COLOR=col,SYMSIZE=0.4
      OPLOT,[2.4,2.4],[ypos,ypos],PSYM=8,COLOR=col,SYMSIZE=0.4
      OPLOT,[2.5,2.5],[ypos,ypos],PSYM=8,COLOR=col,SYMSIZE=0.4
      OPLOT,[2.5,2.5],[ypos+0.06,ypos+0.06],PSYM=8,COLOR=col,SYMSIZE=0.4
      OPLOT,[2.7,2.7],[ypos,ypos],PSYM=8,COLOR=col,SYMSIZE=0.4

    ; Plot 5-95% bar for forcing constraint
    col = 0
    ypos3 = ypos2-0.08
    ypos=ypos3
    ysize=0.04
    p95=p95_dF
    p05=p05_dF
    med=p50_dF
    OPLOT, [med, med],[ypos, ypos+ysize],color=col, thick=4
    OPLOT, [p05, p95],[ypos+ysize/2., ypos+ysize/2.],color=col, thick=4  
    OPLOT, [p95,p95],[ypos, ypos+ysize],color=col, thick=4
    OPLOT, [p05,p05],[ypos, ypos+ysize],color=col, thick=4

    ; ********** AEROSOL FORCING *******************

    ; Plot aerosol forcing PDF
    OPLOT,xbin,pdf_dF_AER,THICK=5,COLOR=11

    ; Plot 5-95% bar for AR5 Aerosol forcing
    col = 2
    ypos1 = +1.3
    ypos=ypos1
    ysize=0.04
    p95=dF_AR5_AER_p95
    p05=dF_AR5_AER_p05
    med=dF_AR5_AER_med
    OPLOT, [med, med],[ypos, ypos+ysize],color=col, thick=4
    OPLOT, [p05, p95],[ypos+ysize/2., ypos+ysize/2.],color=col, thick=4  
    OPLOT, [p95,p95],[ypos, ypos+ysize],color=col, thick=4
    OPLOT, [p05,p05],[ypos, ypos+ysize],color=col, thick=4

    ; Plot 5-95% bar for CMIP5 aerosol forcing & individual models
    col = 15
    ypos2 = ypos1-0.08
    ypos=ypos2
    ysize=0.06
    p95=dF_AER_CMIP5_p95
    p05=dF_AER_CMIP5_p05
    med=dF_AER_CMIP5_med
    OPLOT, [med, med],[ypos, ypos+ysize],color=col, thick=4
    OPLOT, [p05, p95],[ypos+ysize/2., ypos+ysize/2.],color=col, thick=2
    OPLOT, [p95,p95],[ypos, ypos+ysize],color=col, thick=4
    OPLOT, [p05,p05],[ypos, ypos+ysize],color=col, thick=4
    zzzz=[0.68,0.84,0.97,1.13,1.37]*(-1.) ; CMIP5 aerosol forcings that are well seperated
    FOR j=0,N_ELEMENTS(zzzz)-1 DO BEGIN
      OPLOT,[zzzz[j],zzzz[j]],[ypos+ysize/2.,ypos+ysize/2.],PSYM=8,COLOR=col,SYMSIZE=0.4
    ENDFOR
      ; Offset individual CMIP5 models that sit on top of each other
      OPLOT,[-1.24,-1.24],[ypos+0.01,ypos+0.01],PSYM=8,COLOR=col,SYMSIZE=0.4
      OPLOT,[-1.24,-1.24],[ypos+0.05,ypos+0.05],PSYM=8,COLOR=col,SYMSIZE=0.4
      OPLOT,[-1.53,-1.53],[ypos+0.01,ypos+0.01],PSYM=8,COLOR=col,SYMSIZE=0.4
      OPLOT,[-1.55,-1.55],[ypos+0.05,ypos+0.05],PSYM=8,COLOR=col,SYMSIZE=0.4
    
    ; Plot 5-95% bar for aerosol forcing constraint
    col = 11
    ypos3 = ypos2-0.08
    ypos=ypos3
    ysize=0.04
    p95=p95_dF_AER
    p05=p05_dF_AER
    med=p50_dF_AER
    OPLOT, [med, med],[ypos, ypos+ysize],color=col, thick=4
    OPLOT, [p05, p95],[ypos+ysize/2., ypos+ysize/2.],color=col, thick=4  
    OPLOT, [p95,p95],[ypos, ypos+ysize],color=col, thick=4
    OPLOT, [p05,p05],[ypos, ypos+ysize],color=col, thick=4

    ; Add labels
    XYOUTS,3.5,ypos1-0.005,'IPCC AR5',COL=2,CHARSIZE=0.8
    XYOUTS,3.5,ypos2-0.005,'CMIP5',COL=9,CHARSIZE=0.8
    XYOUTS,-1,0.85,'Aerosols',COL=11,CHARSIZE=0.8,ALIGNMENT=0.5
    XYOUTS,1.15,0.95,'Total forcing',COL=0,CHARSIZE=0.8,ALIGNMENT=0.5

    OPLOT,[-20,20],[0,0],COL=0,THICK=3

  
 PREND,/keep,/noprint

END
