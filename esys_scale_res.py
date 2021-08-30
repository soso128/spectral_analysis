from numpy import sqrt,exp,poly1d

def esys_rel(e, sknum, scale_poly, res_poly, ntag = False):
    sesys = 1
    resys = 1
    if sknum < 4:
        if sknum != 2:
            sesys=0.9997+0.0008082*e+9.573E-6*e**2+2.128E-8*e**3
            resys=0.9922+0.000509*e-8.8878E-6*e**2+3.08E-8*e**3
            sesys = ((sesys-1)*2)+1
        else:
            sesys=0.9951+0.00145*e+1.9834E-5*e**2+1.9652E-7*e**3
            resys=0.9916+0.0005729*e-1.4034E-5*e**2+7.8073E-8*e**3
            sesys = ((sesys-1)*1.5)+1
    else:
        sesys = poly1d(scale_poly[f'relic{int(ntag)}'])(e)
        resys = poly1d(res_poly[f'relic{int(ntag)}'])(e)
    return sqrt((1-sesys)**2 + (1-resys)**2)

def esys_nue(e, sknum, scale_poly, res_poly, ntag = False):
    sesys = 1
    resys = 1
    if sknum < 4:
        if sknum != 2:
            sesys=0.8864+0.01049*e-4.887E-4*e**2+1.206E-5*e**3- 1.632E-7*e**4+1.1485E-9*e**5-3.2798E-12*e**6
            resys=1.011-0.0003702*e+5.145E-10*e**2+7.501E-9*e**3
            sesys = ((sesys-1)*2)+1
        else:
            sesys=0.8911+0.006248*e-2.143E-4*e**2+3.705E-6*e**3  -3.1287E-8*e**4+1.0361E-10*e**5
            resys=1.012-0.0007307*e+6.0877E-6*e**2-2.1226E-8*e**3
            sesys = ((sesys-1)*1.5)+1
    else:
        sesys = poly1d(scale_poly[f'ccnuemed{int(ntag)}'])(e)
        resys = poly1d(res_poly[f'ccnuemed{int(ntag)}'])(e)
    return sqrt((1-sesys)**2 + (1-resys)**2)

def esys_ncel_low(e,sknum, scale_poly, res_poly, ntag = False):
    sesys = 1
    resys = 1
    if sknum < 4:
        if sknum != 2:
            sesys=1.015+0.0003052*e+2.1056E-6*e**2 +  0.02042*exp(-0.5*((e-26.78)/11.16)**2)
            resys=0.9909+0.0009907*e-2.411E-5*e**2+1.3837E-7*e**3
            sesys = ((sesys-1)*2)+1
        else:
            sesys=1.036+0.000073*e+3.9564E-6*e**2+   0.05188*exp(-0.5*((e-33.54)/12.57)**2)
            resys=0.959+0.00427*e-1.3474E-4*e**2+ 1.5786E-6*e**3-6.4322E-9*e**4
            sesys = ((sesys-1)*1.5)+1
    else:
        sesys = poly1d(scale_poly[f'nclow{int(ntag)}'])(e)
        resys = poly1d(res_poly[f'nclow{int(ntag)}'])(e)
    return sqrt((1-sesys)**2 + (1-resys)**2)

def esys_ncel_med(e,sknum, scale_poly, res_poly, ntag = False):
    sesys = 1
    resys = 1
    if sknum < 4:
        if sknum != 2:
            sesys=0.9684+0.00177*e-2.411E-5*e**2+1.227E-7*e**3  + 0.08*exp(-0.5*((e-20.21)/6.838)**2)
            resys=0.9717+0.001483*e-2.7838E-5*e**2+1.4329E-7*e**3  + 0.020663*exp(-0.5*((e-20.9)/8.3981)**2)
            sesys = ((sesys-1)*2)+1
        else:
            sesys=1.011-0.0002791*e+2.776E-6*e**2+   0.1487*exp(-0.5*((e-22.13)/6.278)**2)
            resys=1.01-0.00054*e+2.4056E-6*e**2+  0.02402*exp(-0.5*((e-24.6)/6.413)**2)
            sesys = ((sesys-1)*1.5)+1
    else:
        sesys = poly1d(scale_poly[f'ncmed{int(ntag)}'])(e)
        resys = poly1d(res_poly[f'ncmed{int(ntag)}'])(e)
    return sqrt((1-sesys)**2 + (1-resys)**2)

def esys_ncel_high(e,sknum, scale_poly, res_poly, ntag = False):
    sesys = 1
    resys = 1
    if sknum < 4:
        if sknum != 2:
            sesys=1.015+0.0003052*e+2.1056E-6*e**2 +  0.02042*exp(-0.5*((e-26.78)/11.16)**2)
            resys=0.9909+0.0009907*e-2.411E-5*e**2+1.3837E-7*e**3
            sesys = ((sesys-1)*2)+1
        else:
            sesys=1.036+0.000073*e+3.9564E-6*e**2+   0.05188*exp(-0.5*((e-33.54)/12.57)**2)
            resys=0.959+0.00427*e-1.3474E-4*e**2+ 1.5786E-6*e**3-6.4322E-9*e**4
            sesys = ((sesys-1)*1.5)+1
    else:
        sesys = poly1d(scale_poly[f'nchigh{int(ntag)}'])(e)
        resys = poly1d(res_poly[f'nchigh{int(ntag)}'])(e)
    return sqrt((1-sesys)**2 + (1-resys)**2)

def esys_mupi_low(e,sknum, scale_poly, res_poly, ntag = False):
    sesys = 1
    resys = 1
    if sknum < 4:
        if sknum != 2:
            sesys=0.9974+0.0002159*e-5.58985E-6*e**2+6.546E-8*e**3
            resys=0.9991+0.0001937*e-8.5064E-6*e**2+4.9135E-8*e**3
            sesys = ((sesys-1)*2)+1
        else:
            sesys=0.9888+0.00016*e-4.0465E-6*e**2+7.36E-8*e**3
            resys=1.0038-0.0002594*e-1.8287E-6*e**2+    2.0917E-8*e**3
            sesys = ((sesys-1)*1.5)+1
    else:
        sesys = poly1d(scale_poly[f'mupilow{int(ntag)}'])(e)
        resys = poly1d(res_poly[f'mupilow{int(ntag)}'])(e)
    return sqrt((1-sesys)**2 + (1-resys)**2)

def esys_mupi_med(e,sknum, scale_poly, res_poly, ntag = False):
    sesys = 1
    resys = 1
    if sknum < 4:
        if sknum != 2:
            sesys=0.9760+0.004208*e-2.32E-4*e**2+5.293E-6*e**3-   5.5545E-8*e**4+2.2825E-10*e**5
            resys=0.9985+0.0002817*e-1.0873E-5*e**2+6.558E-8*e**3
            sesys = ((sesys-1)*2)+1
        else:
            sesys=0.9966+0.0001365*e-3.4537E-6*e**2+ 5.9911E-8*e**3
            resys=1.002-0.0001752*e-3.3445E-6*e**2+   2.9225E-8*e**3
            sesys = ((sesys-1)*1.5)+1
    else:
        sesys = poly1d(scale_poly[f'mupimed{int(ntag)}'])(e)
        resys = poly1d(res_poly[f'mupimed{int(ntag)}'])(e)
    return sqrt((1-sesys)**2 + (1-resys)**2)

def esys_mupi_high(e,sknum, scale_poly, res_poly, ntag = False):
    sesys = 1
    resys = 1
    if sknum < 4:
        if sknum != 2:
            sesys=0.6771+0.059734*e-4.3985E-3*e**2+1.7147E-4*e**3  -3.8088E-6*e**4+4.8773E-8*e**5-3.3778E-10*e**6  +9.8481E-13*e**7
            resys=0.9994+6.4007E-5*e-2.9394E-6*e**2
            sesys = ((sesys-1)*2)+1
        else:
            sesys=1.0583-0.0065657*e+3.2169E-4*e**2-5.9384E-6*e**3+4.2056E-8*e**4-8.8362E-11*e**5
            resys=0.9978+0.0001586*e-9.3226E-6*e**2+   6.0812E-8*e**3
            sesys = ((sesys-1)*1.5)+1
    else:
        sesys = poly1d(scale_poly[f'mupihigh{int(ntag)}'])(e)
        resys = poly1d(res_poly[f'mupihigh{int(ntag)}'])(e)
    return sqrt((1-sesys)**2 + (1-resys)**2)

def esys(e,sknum,region,sigbkg, scale_poly = None, res_poly = None, ntag = False):
    if region == 0:
        if sigbkg == 0: return esys_nue(e,sknum, scale_poly, res_poly, ntag)
        if sigbkg == 2: return esys_ncel_low(e,sknum, scale_poly, res_poly, ntag)
        if sigbkg == 3: return esys_mupi_low(e,sknum, scale_poly, res_poly, ntag)
        if sigbkg == 5: return esys_rel(e,sknum, scale_poly, res_poly, ntag)
    elif region == 1:
        if sigbkg == 0: return esys_nue(e,sknum, scale_poly, res_poly, ntag)
        elif sigbkg == 2: return esys_ncel_med(e,sknum, scale_poly, res_poly, ntag)
        elif sigbkg == 3: return esys_mupi_med(e,sknum, scale_poly, res_poly, ntag)
        if sigbkg == 5: return esys_rel(e,sknum, scale_poly, res_poly, ntag)
    elif region == 2:
        if sigbkg == 0: return esys_nue(e,sknum, scale_poly, res_poly, ntag)
        if sigbkg == 2: return esys_ncel_high(e,sknum, scale_poly, res_poly, ntag)
        elif sigbkg == 3: return esys_mupi_high(e,sknum, scale_poly, res_poly, ntag)
        if sigbkg == 5: return esys_rel(e,sknum, scale_poly, res_poly, ntag)
    return 0
