""" bpz_tools.py: Contains useful functions for I/O and math.
    TO DO:
      Include higher order interpolations
"""

from numpy import *
import numpy as N
from string import *
import useful as U
from useful import *
import numpy
import os,sys
from glob import *
from scipy import optimize,interpolate
import re
import matplotlib.pyplot as plt
# from pylab import *

try:
    import multiprocessing as MP
    mp=1
except:
    print "multiprocessing not available"
    mp=0

global z_ab

h_cgs=6.6260755e-27 # ergs s
clight_AHz=2.99792458e18
Vega='Vega_reference'

#Smallest number accepted by python
tiny=N.finfo(N.double).tiny
etiny=log(tiny)

#This quantities are used by the AB files
zmax_ab=12.
dz_ab=0.001
ab_clip=1e-6
z_ab=N.arange(0.,zmax_ab,dz_ab) 

#Initialize path info
bpz_dir=os.getenv('BPZPATH')
fil_dir=bpz_dir+'/FILTER/'
sed_dir=bpz_dir+'/SED/'
ab_dir=bpz_dir+'/AB/'

#Fit_dm
dmB=0.0005
mB=12.

def get_lib(lib="B10.list"):
    if lib[-5:]<>".list":
        lib+=".list"
    try:
        templates=get_str(lib,0)
    except:
        templates=get_str(sed_dir+lib,0)
    return templates
    
def seesed(cosa="LRG.sed",x1=2000.,x2=12000.,norm=1.,tipo=None,sed_dir=sed_dir):
    if type(cosa)==type(""):
        if cosa[-4:]=="list":
            try:
                templates=get_str(sed_dir+cosa,0)
            except:
                templates=get_str(cosa,0)            
        else:
            templates=[cosa]
    else:
        templates=cosa

    for i in range(len(templates)):
        x,y=get_sed(templates[i],sed_dir=sed_dir)
        g=greater_equal(x,x1)*less_equal(x,x2)
        print templates[i],max(y[g]),median(y[g])
        y=y*norm/N.trapz(y[g],x[g])*(x2-x1)
        ymin=min(y[g])*0.95
        ymax=max(y[g])*1.05
        if tipo==None:
            plot(x,y)
        else:
            plot(x,y,tipo)
        axis((x1,x2,ymin,ymax))

def seefilters(filters="JPAS.filters",gr=1):
    try:
        filters=U.get_str(filters,0)
    except:
        filters=U.get_str(fil_dir+filters,0)
    xmin=11000.
    xmax=0.
    ymax=0.
    for i in range(len(filters)):
        x,y=get_filter(filters[i])
        if x[0]<xmin: xmin=x[0]
        if x[-1]>xmax: xmax=x[-1]
        if y.max()>ymax: ymax=y.max()
        plt.plot(x,y)
    if gr: plt.grid()
    plt.axis((xmin-100.,xmax+100,0.,ymax*1.06))
    plt.xlabel("Wavelength")
    plt.ylabel("Throughput")
    plt.show()

class load_mags:
    def __init__(self,cat="COSMOSIB_DR2.cat",
                 cols="COSMOSIB_DR2_B10.columns"):
        filters=get_filter_list(cols)
        self.m={}
        self.dm={}
        dict=get_filter_dict(cols,cap=0)
        ms=get_2Darray(cat)
        for i in range(len(filters)):
            i1,i2=int(dict[filters[i]][0])-1,int(dict[filters[i]][1])-1
            self.m[filters[i]],self.dm[filters[i]]=ms[:,i1],ms[:,i2]
            self.m[filters[i]]+=float(dict[filters[i]][4]) # Correct for zp
            print "%s loaded as .m['%s'],.dm['%s']" % (filters[i],filters[i],filters[i])
        self.filters=filters

def load_m0(cat="COSMOSIB_DR2.cat",cols="COSMOSIB_DR2_B10.columns"):
    im0=int(params_file(cols)["M_0"])-1
    return get_data(cat,im0)

class convert_fluxes2mags:
    #Given a catalog defined in fluxes, it will convert it to BPZ-like mags
    #(BPZ can run on fluxes, but sometimes this is an interesting exercise)
    #It doesn't deal with Z_S, etc., those will have to be changed by hand
    def __init__(self,fluxcat="cosmos-1.spec.cat",fluxcols="cosmos-1.spec.columns",
                 magcat="cosmos-1.spec_mag.cat",zp=25.):
        a=load_mags(fluxcat,fluxcols)
        fils=get_filter_list(fluxcols)
        d=params_file(fluxcols)
        buffer=""
        c=0
        mags=[]
        for k in fils:
            c+=2
            m,dm=sex2bpzmags(a.m[k],a.dm[k],zp=zp)
            buffer+="%s %i,%i AB 0 0\n" % (k,c-1,c)
            mags.append(m)
            mags.append(dm)
        if "Z_S" in d.keys():
            print d["Z_S"]
            z_s=get_data(fluxcat,int(d["Z_S"])-1)
            mags.append(z_s)
            buffer+="%s %i\n" % ("Z_S",c+1)


        put_data(magcat,tuple(mags))
        open(os.path.splitext(magcat)[0]+".columns","w").write(buffer)
            
#Auxiliary synthetic photometry functions
def flux(xsr,ys,yr,ccd='yes',units='nu'):
    """
    Flux of spectrum ys observed through response yr,
    both defined on xsr 
    Both f_nu and f_lambda have to be defined over lambda
    If units=nu, it gives f_nu as the output
    """
    if ccd=='yes': yr=yr*xsr
    norm=N.trapz(yr,xsr)
    f_l=N.trapz(ys*yr,xsr)/norm
    if units=='nu':
	lp=sqrt(norm/N.trapz(yr/xsr/xsr,xsr))      #Pivotal wavelenght	
	return f_l*lp**2/clight_AHz
    else: return f_l

def sed2photons(sed,filtro,seddir=sed_dir,fildir=fil_dir,
                ccd="yes",cal="AB",
                mirror=4.41,
                area_object=1., #This is the area of the sky which contains mag (sq. arcsec)
                area_det=1., #This is the area where we count photons (sq. arcsec)
                t_exp=250,
                x_thru=None,
                y_thru=None):
    
    """If SED properly calibrated in ergs/s/cm^2/A
    returns photons through filter & telescope """
    x,y,r=get_sednfilter(sed,filtro,seddir,fildir)
    ph=fl2photons(x,y)
    if x_thru<>None: r=match_resol(x_thru,y_thru,x)        
    return N.trapz(ph*r,x)*mirror*t_exp*area_det/area_object
    
def mag2photons(mag,
                filtro,
                ccd='yes',
                cal="AB",
                mirror=4.41,
                area_object=1., #This is the area which contains mag (sq. arcsec)
                area_det=1., #This is the area where we count photons (sq. arcsec)
                t_exp=250.):
    """Converts AB magnitudes to total photons s^-1 m^-2"""
    if cal<>"AB": mag=VegatoAB(mag,"filter")
    lp=filter_center(filtro,ccd)
    return fl2photons(lp,ABtofl(mag,filtro,ccd))*effective_width(filtro,ccd)*mirror*t_exp*area_det/area_object

def fl2mag(sed,filtro,seddir=sed_dir,fildir=fil_dir,ccd="yes",cal="AB"):
    """If SED properly calibrated in ergs/s/cm^2/A
       returns AB magnitude through filter """
    if type(sed)==type("filtro"):
        x,y,r=get_sednfilter(sed,filtro,seddir,fildir)
    else:
        x_sed=sed[0]
        y_sed=sed[1]
        nsed=len(x_sed)
        x_res,y_res=get_filter(filtro)
        nres=len(x_res)
        i1=searchsorted(x_sed,x_res[0])-1
        i1=maximum(i1,0)
        i2=searchsorted(x_sed,x_res[nres-1])+1
        i2=minimum(i2,nsed-1)
        r=match_resol(x_res,y_res,x_sed[i1:i2])
        r=where(less(r,0.),0.,r) #Transmission must be >=0
        x=x_sed[i1:i2]
        y=y_sed[i1:i2]
        
    fnu=flux(x,y,r,ccd=ccd,units="nu")
    if cal=="AB":
        return AB(fnu)
    else:
        return ABtoVega(AB(fnu),filtro)


def ABtofl(ABmag,filtro,ccd='yes'):
    """Converts AB magnitudes to flux in ergs s^-1 cm^-2 AA^-1
    Accepts normal filter names or a wavelenght array"""
    try:
        lp=pivotal_wl(filtro,ccd)
        f=AB2Jy(ABmag)
        return f/lp**2*clight_AHz*1e-23
    except:
        f=AB2Jy(ABmag)
        return f/filtro**2*clight_AHz*1e-23

def fnu2fl(l,flux):
    return flux*clight_AHz/l**2

def fl2fnu(l,flux):
    return flux/clight_AHz*l**2

def fl2photons(l,flux):
    """input in ergs s^-1 cm^-2 A^-1
       output photons/s/m^2/A"""
    return flux/(h_cgs*clight_AHz/l)*1e4

def photons2fl(l,flux):
    """input photons/s/m^2/A
    output in ergs s^-1 cm^-2 A^-1"""
    return flux*((h_cgs*clight_AHz/l)/1e4)

def pivotal_wl(filtro,ccd='yes'):
    xr,yr=get_filter(filtro)
    if ccd=='yes': yr=yr*xr
    norm=N.trapz(yr,xr)
    return sqrt(norm/N.trapz(yr/xr/xr,xr))  

def effective_width(filtro,ccd="yes"):
    xr,yr=get_filter(filtro)
    if ccd=='yes': yr=yr*xr/mean(xr)
    return N.trapz(yr,xr)

def filter_center(filtro,ccd='yes',fil_dir=fil_dir):
    """Estimates the central wavelenght of the filter"""
    if type(filtro)==type(""):
        xr,yr=get_filter(filtro)
    else:
        xr=filtro[0]
        yr=filtro[1]
    if ccd=='yes': yr=yr*xr
    return N.trapz(yr*xr,xr)/N.trapz(yr,xr)

def effective_wavelength(filtro,ccd='yes'):
    """Estimates the effective wavelenght of the filter"""
    if type(filtro)==type(""):
        xr,yr=get_filter(filtro)
    else:
        xr=filtro[0]
        yr=filtro[1]
    if ccd=='yes': yr=yr*xr
    return N.trapz(yr*xr*xr,xr)/N.trapz(yr*xr,xr)

def filter_fwhm(filtro,ccd='yes'):
    if type(filtro)==type(""):
        xr,yr=get_filter(filtro)
    else:
        xr=filtro[0]
        yr=filtro[1]
    np=len(xr)
    if ccd=='yes': yr=yr*xr/mean(xr)
    imax=argmax(yr)
    ymax=yr[imax]
    xmax=xr[imax]
    fhr=match_resol(yr[::-1][:(np-imax)]/ymax,xr[::-1][:(np-imax)],0.5)
    fhl=match_resol(yr[:imax]/ymax,xr[:imax],0.5)
    return fhr-fhl

def filter_fwhm_i(filtro,ccd='yes'):
    xr,yr=get_filter(filtro)
    if ccd=='yes': yr=yr*xr/mean(xr)
    imax=argmax(yr)
    ymax=yr[imax]
    xmax=xr[imax]
    ih_1=argmin(abs(yr[:imax]-ymax/2.))
    ih_2=argmin(abs(yr[imax:]-ymax/2.))+imax
    return xr[ih_2]-xr[ih_1]
   
def AB(flux):
    """AB magnitude from f_nu"""
    return -2.5*log10(flux)-48.60

def fnu(AB):
    """f_nu from AB magnitude"""
    return 10.**(-0.4*(AB+48.60))

def flux2mag(flux):
    """Convert arbitrary flux to magnitude"""
    return -2.5*log10(flux) 

def Jy2AB(flux):
    """Convert flux in Jy to AB magnitudes"""
    return -2.5*log10(flux*1e-23)-48.60

def AB2Jy(ABmag):
    """Convert AB magnitudes to Jansky"""
    return 10.**(-0.4*(ABmag+48.60))/1e-23    
    
def mag2flux(mag):
    """Convert flux to arbitrary flux units"""
    return 10.**(-.4*mag)

def e_frac2mag(fracerr):
    """Convert fractionary flux error to mag error"""
    return 2.5*log10(1.+fracerr)
            
def e_mag2frac(errmag):
    """Convert mag error to fractionary flux error"""
    return 10.**(.4*errmag)-1.

def flux_det(aperture,pixelnoise,s2noise=1):
    """Given an aperture, the noise per pixel and the 
       signal to noise, it estimates the detection flux limit"""
    npixels=pi*(aperture/2.)**2
    totalnoise=sqrt(npixels)*pixelnoise
    return s2noise*totalnoise

def get_limitingmagnitude(m,dm,n_sigma=1.,dm_int=0.2,tag=None,plots=0):
    """Given a list of magnitudes and magnitude errors,
    calculate by extrapolation the n_sigma error limit"""
    np=len(m)
    g=less(m,99.)*greater(m,00.)
    y,x=autobin_stats(compress(g,dm),compress(g,m),n_points=15,stat="median")

    # In case the error dm_int is not contained in the data set
    if dm_int >= y[-2] or dm_int < y[0]: 
        dm_int=y[-3] # Take third point from the end to avoid limit effects
    mlim=match_resol(y,x,dm_int)-flux2mag(1./n_sigma/e_mag2frac(dm_int))

    if plots:
        #BPZsim approximation
        sn=n_sigma/mag2flux(mlim)*mag2flux(x)
        plot(dm[g],m[g],".")
        plot(y,x)
        plot(e_frac2mag(1./sn),x)
        show()

    return mlim

def fit_dm(m,dm,k=0.045,alpha=1.85,dmB=dmB,mB=mB,
                       free=(1,0,0,0),plots=1):
    # It return the parameters of a fit to the dm vs m relationship of the shape
    # dm = dm0*exp(k*(m-mo)**alpha)
    k=Parameter(k)
    alpha=Parameter(alpha)
    dmB=Parameter(dmB)
    mB=Parameter(mB)
    
    def f(x): return dmB()*exp(k()*(x-mB())**alpha())
    
    pars0=[k,alpha,dmB,mB]
    pars=[]
    if free[0]: pars.append(k)
    if free[1]: pars.append(alpha)
    if free[2]: pars.append(dmB)
    if free[3]: pars.append(mB)

    g=greater(m,0)*less(m,30.)*less(dm,0.2)
    fit(f,pars,y=dm[g],x=m[g])
    
    xm=arange(min(m[g])-.1,max(m[g])+.1,0.01)
    if plots:
        plot(m[g],dm[g],".")
        plot(xm,f(xm))
        show()

    chi2=sum((dm[g]-f(m[g]))/dm[g])**2/g.sum()
    return map(lambda x: x(),pars0)

def generate_dm(m,k=0.045,alpha=1.85,dmB=dmB,mB=mB,minsn=5.,
                plots=1):
    # The error goes as dmB*exp(k*(m-mB)**alpha) up to 5 sigma
    # At fainter magnitudes the absolute flux error is kept constant 
    # Calculate magnitude at which we fix the minimum flux error (at 5 sigma)
    dmsn=e_frac2mag(1./minsn)
    msn=mB+(1./k*log(dmsn/dmB))**(1./alpha)
    dfomin=mag2flux(msn)/minsn #Minimum error
    #Now calculate magnitude errors
    try:
        m.shape
    except:
        m=array([m])
    dm=where(m<=msn,
             dmB*exp(k*(m-mB)**alpha),             
             e_frac2mag(dfomin/mag2flux(m)))
    return dm

def k_fromdm(m=24.,dm=0.2,alpha=1.85,dmB=dmB,mB=mB):
    return log(dm/dmB)/(m-mB)**alpha

def m_sigma(sn=1.,k=0.045,alpha=1.85,dmB=dmB,mB=mB): 
    dm=e_frac2mag(1./sn)
    def f(x): return (generate_dm(x,k=k,alpha=alpha,dmB=dmB,mB=mB)[0]-dm)**2
    return float(optimize.fmin_powell(f,mB+1.))

def limmag(cat,columns,n_sigma=1.,dm_int=0.2,plots=0):
    """Give a photometric catalog and a colums file, 
       calculate the limiting magnitudes in each of the filters"""
    #Get filters
    filters=get_filter_list(columns)
    mags=zeros(len(filters),'float')
    dict=get_filter_dict(columns)
    for i in range(len(filters)):
        i1,i2=int(dict[filters[i]][0])-1,int(dict[filters[i]][1])-1
        m,dm=get_data(cat,(i1,i2))
        mags[i]=get_limitingmagnitude(m,dm,n_sigma,dm_int,tag=filters[i],plots=plots)
        print filters[i],mags[i]
    return filters,mags

#Synthetic photometry 

#class vetau_meiksin:
#    def __init__(self):
#        if not os.path.exists(ab_dir+"tau_meiksin.npy"):
#            meiksin=get_2Darray(sed_dir+"transmission_table.dat")
#            l=meiksin[:,0]
            

class vetau_madau:
    """
    Madau 1995 extinction for several galaxy spectra at redshift z 
    defined on a wavelenght grid wl (for all z)
    """
    def __init__(self):
        l=arange(912.,12000.,1.)  
        if not os.path.exists(sed_dir+"tau_madau.npy"):
            taumadau=ones((len(l),len(z_ab)),"float")
            for i in range(len(z_ab)):
                taumadau[:,i]=etau_madau(l,z_ab[i])
            numpy.save(sed_dir+"tau_madau.npy",taumadau)
        else:
            taumadau=numpy.load(sed_dir+"tau_madau.npy")
        self.taumadau=taumadau
        self.z_ab=z_ab
        self.l=l
   
    def e(self,wl,z):
        iz=searchsorted(self.z_ab,z)
        il=searchsorted(self.l,wl)
        il=clip(il,0,self.taumadau.shape[0]-1)
        return self.taumadau[il,:][:,iz]

def etau_madau_old(wl,z):
    """
    Madau 1995 extinction for a galaxy spectrum at redshift z 
    defined on a wavelenght grid wl
    """
    ll=912.
    c=array([3.6e-3,1.7e-3,1.2e-3,9.3e-4])
    l=array([1216.,1026.,973.,950.])
    tau=zeros_like(wl)
    xe=1.+z

    #Lyman series
    for i in range(len(l)):
        tau=where(wl<=l[i]*xe,tau+c[i]*np.power(wl/l[i],3.46),tau)

    #Photoelectric absorption
    xc=wl/ll
    xc3=np.power(xc,3)
    tau=where(wl<=ll*xe,
              tau+0.25*xc3*(np.power(xe,.46)-np.power(xc,0.46))\
                  +9.4*np.power(xc,1.5)*(np.power(xe,0.18)-np.power(xc,0.18))\
                  -0.7*xc3*(np.power(xc,-1.32)-np.power(xe,-1.32))\
                  -0.023*(np.power(xe,1.68)-np.power(xc,1.68)),
              tau)

    return where(tau > 700., 0., exp(-tau))


def etau_madau(wl,z):
    """
    Madau 1995 extinction for a galaxy spectrum at redshift z 
    defined on a wavelenght grid wl
    """
    ll=912.
    c=array([3.6e-3,1.7e-3,1.2e-3,9.3e-4])
    l=array([1216.,1026.,973.,950.])
    tau=zeros_like(wl)
    xe=1.+z

    # Lyman series
    for i in range(len(l)):
        i1=searchsorted(wl,l[i]*xe)
        tau[:i1]+=c[i]*(wl[:i1]/l[i])**3.46

    #Photoelectric absorption
    il=searchsorted(wl,ll*xe)
    xc=(wl/ll)[:il]
    xc3=xc**3
    tau[:il]+=(0.25*xc3*(xe**.46-xc**0.46)\
                   +9.4*xc**1.5*(xe**0.18-xc**0.18)\
                   -0.7*xc3*(xc**-1.32-xe**-1.32)\
                   -0.023*(xe**1.68-xc**1.68))

    return exp(-tau)


def etau(wl,z):
    """
    Madau 1995 and Scott 2000 extinction for a galaxy spectrum
    at redshift z observed on a wavelenght grid wl
    """

    n=len(wl)
    l=array([1216.,1026.,973.,950.])
    xe=1.+z

    #If all the spectrum is redder than (1+z)*wl_lyman_alfa 
    if wl[0]> l[0]*xe: return zeros(n)+1.

    #Extinction coefficients

    c=array([1.,0.47,0.33,0.26])
    if z>4.:
        #Numbers from Madau paper
        coeff=0.0036
        gamma=2.46
    elif z<3:
        #Numbers from Scott et al. 2000 paper
        coeff=0.00759
        gamma=1.35
    else:
        #Interpolate between two numbers
        coeff=.00759+(0.0036-0.00759)*(z-3.)
        gamma=1.35+(2.46-1.35)*(z-3.)
    c=coeff*c

    ll=912.
    tau=wl*0.
    i1=searchsorted(wl,ll)
    i2=n-1
    #Lyman series absorption
    for i in range(len(l)):
	i2=searchsorted(wl[i1:i2],l[i]*xe)
	tau[i1:i2]=tau[i1:i2]+c[i]*(wl[i1:i2]/l[i])**(1.+gamma)
        
    if ll*xe < wl[0]: return exp(-tau)

    #Photoelectric absorption
    xe=1.+z
    i2=searchsorted(wl,ll*xe)
    xc=wl[i1:i2]/ll
    xc3=xc**3
    tau[i1:i2]=tau[i1:i2]+\
		(0.25*xc3*(xe**.46-xc**0.46)\
		+9.4*xc**1.5*(xe**0.18-xc**0.18)\
		-0.7*xc3*(xc**(-1.32)-xe**(-1.32))\
		-0.023*(xe**1.68-xc**1.68))
    return exp(-tau)

def get_sednfilter(sed,filtro,seddir=sed_dir,fildir=fil_dir):
    #Gets a pair of SED and filter from the database
    #And matches the filter resolution to that of the spectrum
    #where they overlap
    """Usage:
    xs,ys,yr=get_sednfilter(sed,filtro)
    """
    #Figure out the correct names
    #if seddir==None:
    #    if sed[-4:]<>'.sed':sed=sed+'.sed'
    #    sed=sed_dir+sed
    #else:
    #    sed=seddir+"/"+sed
    #if fildir==None:
    #    if filtro[-4:]<>'.res':filtro=filtro+'.res'
    #    filter=fil_dir+filtro
    #Get the data
    x_sed,y_sed=get_sed(sed,seddir)
    nsed=len(x_sed)
    x_res,y_res=get_filter(filtro,fildir)
    nres=len(x_res)
    if not ascend(x_sed):
        print
        print 'Warning!!!'
        print 'The wavelenghts in %s are not properly ordered' % sed
        print 'They should start with the shortest lambda and end with the longest'        
    if not ascend(x_res):
        print
        print 'Warning!!!'
        print 'The wavelenghts in %s are not properly ordered' % filtro
        print 'They should start with the shortest lambda and end with the longest'
        
    #Define the limits of interest in wavelenght
    i1=searchsorted(x_sed,x_res[0])-1
    i1=maximum(i1,0)
    i2=searchsorted(x_sed,x_res[nres-1])+1
    i2=minimum(i2,nsed-1)
    r=match_resol(x_res,y_res,x_sed[i1:i2])
    r=where(less(r,0.),0.,r) #Transmission must be >=0
    return x_sed[i1:i2],y_sed[i1:i2],r
        
def get_sed(sed,sed_dir=sed_dir):
    #Get x_sed,y_sed from a database spectrum
    """Usage:
    xs,ys=get_sed(sed)
    """
    if sed[-4:]<>'.sed':sed=sed+'.sed'
    sed=sed_dir+sed
    #Get the data
    x,y=get_data(sed,range(2))
    if not ascend(x):
        print
        print 'Warning!!!'
        print 'The wavelenghts in %s are not properly ordered' % sed
        print 'They should start with the shortest lambda and end with the longest'
    return x,y


def get_filter_list(columns):
    #Get a filter list from a .columns file
    f=get_str(columns,0)
    for item in ["Z_S","M_0","OTHER","ID","Z_S","X","Y"]:
        if item in f: f.remove(item)
    return f

def get_filter_dict(columns,cap=0):
    #Get a filter dictionary from a .columns file
    d=params_file(columns,cap=cap)
    for item in ["Z_S","M_0","OTHER","ID","Z_S","X","Y"]:
        try:
            del(d[item])
        except:
            continue
    return d

def get_filter(filtro,fil_dir=fil_dir):
    #Get x_res,y_res from a database spectrum
    """Usage:
    xres,yres=get_filter(filtro)
    """
    #Figure out the correct names
    if filtro[-4:]<>'.res':filtro=filtro+'.res'
    filtro=fil_dir+filtro
    #Get the data
    x,y= get_data(filtro,range(2))
    if not ascend(x):
        print
        print 'Warning!!!'
        print 'The wavelenghts in %s are not properly ordered' % filtro
        print 'They should start with the shortest lambda and end with the longest'
    return x,y


def sed_corrcoef(sed1,sed2,dir=None):
    if dir==None: dir=sed_dir
    #Calculate the correlation coefficient between two seds
    xs1,ys1=get_sed(sed1,sed_dir=dir)
    xs2,ys2=get_sed(sed2,sed_dir=dir)
    dx1=(xs1[-1]-xs1[0])/float(len(xs1))
    dx2=(xs2[-1]-xs2[0])/float(len(xs2))
    #Define the limits of overlap in wavelenght
    x1=max(xs1[0],xs2[0])
    x2=min(xs1[-1],xs2[-1])

    if dx1<=dx2: 
        g=less(xs1,x2)*greater(xs1,x1)
        x=xs1[g]
        y1=ys1[g]
        y2=match_resol(xs2,ys2,x)
    else:
        g=less(xs2,x2)*greater(xs2,x1)
        x=xs2[g]
        y2=ys2[g]
        y1=match_resol(xs1,ys1,x)

    return numpy.corrcoef(y1,y2)[0][1]

def redshift(wl,flux,z,abs="yes"):
    """ Redshift spectrum y defined on axis x 
      to redshift z
      Usage:
	 y_z=redshift(wl,flux,z) 
    """
    if z==0.: return flux
    else: 
	f=match_resol(wl,flux,wl/(1.+z))
        if abs=="yes":
            return where(less(f,0.),0.,f)
        else:
            return f

def normalize(x,y,m,filtro="i_SDSS",units="nu",np=100):
    x0,r0=get_filter(filtro)
    xs=linspace(x0[0],x0[-1],np)
    r=match_resol(x0,r0,xs)
    ys=match_resol(x,y,xs)
    fl=y*ABtofl(m,filtro)/(N.trapz(ys*r*xs,xs)/N.trapz(r*xs,xs))
    if units=="nu":
        return fl2fnu(x,fl)
    else:
        return fl

def observe(sed,z=0.01,ABmag=20.,filtro="i_SDSS",units="lambda",madau=1):
    xs,ys=get_sed(sed)
    yz=redshift(xs,ys,z)
    if madau: yz*=etau_madau(xs,z)
    # Units ergs/s/cm^2/A (inputmag has to be AB)
    return xs,normalize(xs,yz,ABmag,filtro,units)    
    
#class Normalize:
#    def __init__(self,x_sed,y_sed,m,filter='F814W_WFPC2',units='nu'):
#        """Normalizes a spectrum (defined on lambda) to 
#        a broad band (AB) magnitude and transforms the 
#        spectrum to nu units""
#        Usage:
#        normflux=normalize(wl,spectrum,m,filter='F814W_WFPC2')
#        """
#        if filter[-4:]<>'.res':filter=filter+'.res'
#        filter=fil_dir+filter
#        x_res,y_res=get_data(filter,range(2))
#        nres=len(x_res)
#        nsed=len(x_sed)
#        i1=searchsorted(x_sed,x_res[0])-1
#        i1=maximum(i1,0)
#        i2=searchsorted(x_sed,x_res[nres-1])+1
#        i2=minimum(i2,nsed-1)
#        r=match_resol(x_res,y_res,x_sed[i1:i2])
#        r=where(less(r,0.),0.,r) #Transmission must be >=0
#        flujo=flux(x_sed[i1:i2],y_sed[i1:i2],r,ccd='yes',units='nu')
#        self.norm=flujo/mag2flux(m)
#        if units=='nu': self.flux_norm = y_sed*x_sed*x_sed/clight_AHz/self.norm
#        else:           self.flux_norm = y_sed/self.norm

class NormalizeSED:
    """Given a spectrum fl defined in ergs/s/AA/cm^2,
      a set of filters and a set of magnitudes in those filters
      It renormalizes the spectrum flux to match the photometry
      using a polynomial fit or a direct interpolation"""
    def __init__(self,
                 sed="sky_dark_Puxley.sed",
                 filters=["U_Johnson","B_Harris","V_Harris","R_Harris","I_Cousins","z_SupCam"],
                 mags=[22.0,22.7,21.9,21.0,20.0,18.8],
                 cal="Vega",
                 seddir=sed_dir,
                 fildir=fil_dir,
                 poly=1,
                 ):
        nf=arange(len(mags))
        fluxes=zeros(len(mags),'float')
        f_sed=zeros(len(mags),'float')
        mab=zeros(len(mags),'float')
        eff_l=zeros(len(mags),'float')

        for i in nf:
            if cal=="Vega": mab[i]=VegatoAB(mags[i],filters[i])
            else: mab[i]=mags[i]
           #Calculate fluxes of photometric data
            fluxes[i]=ABtofl(mab[i],filters[i])
            xs,ys,r=get_sednfilter(sed,filters[i],seddir,fildir)
            eff_l[i]=N.trapz(ys*r*xs*xs,xs)/N.trapz(ys*r*xs,xs) #xs extra factor to account for CCD
            #Now calculate integrals of spectrum
            f_sed[i]=N.trapz(ys*r*xs,xs)/trapz(r*xs,xs)
            
        xs,ys=get_sed(sed)
        eff_l=array([xs[0]]+list(eff_l)+[xs[-1]])
        rr=fluxes/f_sed
        rr=array([rr[0]]+list(rr)+[rr[-1]])
        if poly:
            polyc=polyfit(eff_l,rr,len(rr)-1)
            r=polyval(polyc,xs)
        else:
            r=match_resol(eff_l,rr,xs)
        r=where(less(xs,eff_l[1]),rr[1]+xs*0.,r)
        r=where(greater(xs,eff_l[-2]),rr[-2]+xs*0.,r)
        #plot(xs,r,eff_l,rr,"d")
        #show()
        #plot(xs,ys*1e16,xs,ys*r*1e16,eff_l[1:-1],fluxes*1e16,"d")
        #plot(xs,ys*1e16)
        #axis((3500.,11000.,0.,1.))
        #show()
        self.xs=xs
        self.ys=ys*r
  

def obs_spectrum(sed,z,madau=1,x_new=None):
    """Generate a redshifted and madau extincted spectrum"""
    #Figure out the correct names
    if sed[-4:]<>'.sed':sed=sed+'.sed'
    sed=sed_dir+sed
    #Get the data
    x_sed,y_sed=get_data(sed,range(2))

    ##ys_z will be the redshifted and corrected spectrum
    #yr=redshift(x_sed,y_sed,z)

    if x_new==None:
        ys_z=match_resol(x_sed,y_sed,x_sed/(1.+z))
        if madau: ys_z=etau_madau(x_sed,z)*ys_z
        return x_sed,ys_z
    else:
        ys_z=match_resol(x_sed,y_sed,x_new/(1.+z))
        if madau: ys_z=etau_madau(x_new,z)*ys_z
        return x_new,ys_z


#def observe(sed,z,w1=None,w2=None):
#    x,y=get_sed(sed)
#    if w1<>None:
#        g=greater(x,x1)
#    if w2<>None:
#        g*=less(x,x2)
#    else:
#        g=x*1
#    redshift(x,y,z)[g]*etau_madau(


#def order(lib,dir=None,ref="El_CWW.sed"):
#    "Order library by similarity among templates"
#    if dir==None: dir=sed_dir
#    templates=get_str(dir+lib,0)
#    corre=zeros(len(templates))
#    newt=[]
#    corr=[]
#    i=0.
#    while len(newt)<len(templates):
#        if i==0: 
#            t0=ref
#        else: 
#            t0=newt[-1]
#        for j in range(len(templates)):
#            if templates[j] in newt or templates[j]==t0: 
#                corre[j]=-1.
#            else:
#                corre[j]=sed_corrcoef(templates[j],t0,dir=dir)
#            #print t0,templates[j],corre[j]
#        newt.append(templates[argmax(corre)])
#        corr.append(max(corre))
#        print newt[-1],corr[-1]
#        i+=1
#    return newt

def order2(lib,dir=None,ref1="El_CWW.sed",ref2="flat.sed"):
    """Order library by similarity among templates"""
    if type(lib)==type([]) or type(lib)==type((2,)):
        templates=lib
    else:
        if dir==None: dir=sed_dir
        templates=get_str(dir+lib,0)
    score=zeros(len(templates))*1.
    for i in range(len(templates)):
        corr1=sed_corrcoef(templates[i],ref1,dir=dir)
        corr2=sed_corrcoef(templates[i],ref2,dir=dir)
        score[i]=corr1-corr2
        print templates[i],score[i]
    return take(templates,argsort(-score))

#def SFR_halpha(sed):
#    # Do not use for anything serious
#    """Assumes that SED is in ergs s^-1
#       Taken from Calzetti 0707.0467
#       Units M_sun/yr"""
#    x,y=get_sed(sed)
#    lc=add.reduce(match_resol(x,y,(6500,6800.)))*0.5
#    l_halpha=fl2fnu(6563.,maximum(0.,match_resol(x,y,6563)-lc))
#    l24um=fl2fnu(240000,match_resol(x,y,240000.))
#    print sed,l_halpha,l24um
#    return 5.3e-42*(l_halpha+0.031*l24um)

def D4000(sed,units="lambda",red=None):
    """return D4000 defined as the ratio between average fluxes in
       the 4050-4250 AA region and the 3750-3950 AA region"""

    try:
        x,y=get_sed(sed)
    except:
        x=sed[0]
        y=sed[1]

    gblue=less_equal(x,3950)*greater_equal(x,3750.)
    gred=less_equal(x,4250)*greater_equal(x,4050.)

    #Hamilton 1985 & Spinrad
    red=fl2fnu(mean(x[gred]),mean(y[gred]))
    blue=fl2fnu(mean(x[gblue]),mean(y[gblue]))
    D4000=red/blue
    return D4000
        

def nf_z_sed(sed,filtro,z=array([0.]),ccd='yes',units='lambda',madau='yes'):
    """Returns array f with f_lambda(z) or f_nu(z) through a given filtro 
       Takes into account intergalactic extinction. 
       Flux normalization at each redshift is arbitrary 
    """
    if type(z)==type(0.): z=array([z])

    #Figure out the correct names
    if sed[-4:]<>'.sed':sed=sed+'.sed'
    sed=sed_dir+sed
    if filtro[-4:]<>'.res':filtro=filtro+'.res'
    filtro=fil_dir+filtro

    #Get the data
    x_sed,y_sed=get_data(sed,range(2))
    nsed=len(x_sed)
    x_res,y_res=get_data(filtro,range(2))
    nres=len(x_res)

    #Wavelenght range of interest as a function of z
    wl_1=x_res[0]/(1.+z)
    wl_2=x_res[-1]/(1.+z)
    n1=clip(searchsorted(x_sed,wl_1)-1,0,1000000)
    n2=clip(searchsorted(x_sed,wl_2)+1,0,nsed-1)
    
    #Change resolution of filtro
    x_r=x_sed[n1[0]:n2[0]]
    r=match_resol(x_res,y_res,x_r)
    r=where(less(r,0.),0.,r) #Transmission must be >=0

    #Operations necessary for normalization and ccd effects
    if ccd=='yes': r=r*x_r
    norm_r=N.trapz(r,x_r)
    if units=='nu': const=norm_r/N.trapz(r/x_r/x_r,x_r)/clight_AHz
    else: const=1.
    const=const/norm_r
    
    nz=len(z)
    f=zeros(nz)*1.
    for i in range(nz):
        i1,i2=n1[i],n2[i]
        ys_z=match_resol(x_sed[i1:i2],y_sed[i1:i2],x_r/(1.+z[i]))
        if madau<>'no': ys_z=etau_madau(x_r,z[i])*ys_z
        f[i]=N.trapz(ys_z*r,x_r)*const        
    if nz==1: return f[0]
    else: return f

def efflam_z_sed(sed,filtro,z=array([0.]),ccd='yes',madau='yes',sed_dir=sed_dir,fil_dir=fil_dir):
    """
    Returns the effective lambda of a filtro + SED combination
    Takes into account intergalactic extinction. 
    """

    try:
        len(z)
    except:
        z=array([z])

    #Figure out the correct names
    if "/" in sed: sed=sed.split("/")[-1]
    if sed[-4:]<>'.sed':sed=sed+'.sed'

    if "/" in filtro: filtro=filtro.split("/")[-1]
    if filtro[-4:]<>'.res':filtro=filtro+'.res'

    # Get the data
    x_sed,y_sed=get_sed(sed,sed_dir)
    nsed=len(x_sed)
    x_res,y_res=get_filter(filtro,fil_dir)
    nres=len(x_res)
    
    if not ascend(x_sed):
        print
        print 'Warning!!!'
        print 'The wavelenghts in %s are not properly ordered' % sed
        print 'They should start with the shortest lambda and end with the longest'        
        print 'This will probably crash the program'

    if not ascend(x_res):
        print
        print 'Warning!!!'
        print 'The wavelenghts in %s are not properly ordered' % filtro
        print 'They should start with the shortest lambda and end with the longest'
        print 'This will probably crash the program'

    if x_sed[-1]<x_res[-1]: #The SED does not cover the whole filter interval
        print 'Extrapolating the spectrum'
        #Linear extrapolation of the flux using the last 4 points
        #slope=mean(y_sed[-4:]/x_sed[-4:])
        d_extrap=(x_sed[-1]-x_sed[0])/len(x_sed)
        x_extrap=arange(x_sed[-1]+d_extrap,x_res[-1]+d_extrap,d_extrap)
        extrap=lsq(x_sed[-5:],y_sed[-5:])
        y_extrap=extrap.fit(x_extrap)
        y_extrap=clip(y_extrap,0.,max(y_sed[-5:]))
        x_sed=concatenate((x_sed,x_extrap))
        y_sed=concatenate((y_sed,y_extrap))
        #connect(x_sed,y_sed)
        #connect(x_res,y_res)

    #Wavelenght range of interest as a function of z
    wl_1=x_res[0]/(1.+z)
    wl_2=x_res[-1]/(1.+z)
    n1=clip(searchsorted(x_sed,wl_1)-1,0,100000)
    n2=clip(searchsorted(x_sed,wl_2)+1,0,nsed-1)
    
    #Typical delta lambda
    delta_sed=(x_sed[-1]-x_sed[0])/len(x_sed)
    delta_res=(x_res[-1]-x_res[0])/len(x_res)

    #Change resolution of filter
    if delta_res>delta_sed:
        x_r=arange(x_res[0],x_res[-1]+delta_sed,delta_sed)
        print 'Changing filter resolution from %.2f AA to %.2f AA' % (delta_res,delta_sed)
        r=match_resol(x_res,y_res,x_r)
        r=where(less(r,0.),0.,r) #Transmission must be >=0
    else:
        x_r,r=x_res,y_res

    #Operations necessary for normalization and ccd effects
    if ccd=='yes': r=r*x_r

    nz=len(z)
    efflam=zeros(nz)*1.
    for i in range(nz):
        i1,i2=n1[i],n2[i]
        ys_z=match_resol(x_sed[i1:i2],y_sed[i1:i2],x_r/(1.+z[i]))
        if madau<>'no': ys_z=etau_madau(x_r,z[i])*ys_z
        efflam[i]=N.trapz(ys_z*r*x_r*x_r,x_r)/N.trapz(ys_z*r*x_r,x_r)
    if nz==1: return efflam[0]
    else: return efflam


def lf_z_sed(sed,filtro,z=array([0.]),ccd='yes',units='lambda',madau='yes'):
    """
    Returns array f with f_lambda(z) or f_nu(z) through a given filter 
    Takes into account intergalactic extinction. 
    Flux normalization at each redshift is arbitrary 
    """

    try:
        len(z)
    except:
        z=array([z])

    # Figure out the correct names
    if sed[-4:]<>'.sed':sed=sed+'.sed'
    sed=sed_dir+sed
    if filtro[-4:]<>'.res':filtro=filtro+'.res'
    filtro=fil_dir+filtro

    # Get the data
    x_sed,y_sed=get_data(sed,range(2))
    nsed=len(x_sed)
    x_res,y_res=get_data(filtro,range(2))
    nres=len(x_res)
    
    if not ascend(x_sed):
        print
        print 'Warning!!!'
        print 'The wavelenghts in %s are not properly ordered' % sed
        print 'They should start with the shortest lambda and end with the longest'        
        print 'This will probably crash the program'

    if not ascend(x_res):
        print
        print 'Warning!!!'
        print 'The wavelenghts in %s are not properly ordered' % filtro
        print 'They should start with the shortest lambda and end with the longest'
        print 'This will probably crash the program'

    if x_sed[-1]<x_res[-1]: #The SED does not cover the whole filter interval
        print 'Extrapolating the spectrum'
        #Linear extrapolation of the flux using the last 4 points
        #slope=mean(y_sed[-4:]/x_sed[-4:])
        d_extrap=(x_sed[-1]-x_sed[0])/len(x_sed)
        x_extrap=arange(x_sed[-1]+d_extrap,x_res[-1]+d_extrap,d_extrap)
        extrap=lsq(x_sed[-5:],y_sed[-5:])
        y_extrap=extrap.fit(x_extrap)
        y_extrap=clip(y_extrap,0.,max(y_sed[-5:]))
        x_sed=concatenate((x_sed,x_extrap))
        y_sed=concatenate((y_sed,y_extrap))
        #connect(x_sed,y_sed)
        #connect(x_res,y_res)

    #Wavelenght range of interest as a function of z
    wl_1=x_res[0]/(1.+z)
    wl_2=x_res[-1]/(1.+z)
    n1=clip(searchsorted(x_sed,wl_1)-1,0,100000)
    n2=clip(searchsorted(x_sed,wl_2)+1,0,nsed-1)
    
    #Typical delta lambda
    delta_sed=(x_sed[-1]-x_sed[0])/len(x_sed)
    delta_res=(x_res[-1]-x_res[0])/len(x_res)

    #Change resolution of filter
    if delta_res>delta_sed:
        x_r=arange(x_res[0],x_res[-1]+delta_sed,delta_sed)
        print 'Changing filter resolution from %.2f AA to %.2f AA' % (delta_res,delta_sed)
        r=match_resol(x_res,y_res,x_r)
        r=where(less(r,0.),0.,r) #Transmission must be >=0
    else:
        x_r,r=x_res,y_res

    #Operations necessary for normalization and ccd effects
    if ccd=='yes': r=r*x_r
    norm_r=N.trapz(r,x_r)
    if units=='nu': const=norm_r/N.trapz(r/x_r/x_r,x_r)/clight_AHz
    else: const=1.

    const=const/norm_r
    
    nz=len(z)
    f=zeros(nz)*1.
    for i in range(nz):
        i1,i2=n1[i],n2[i]
        ys_z=match_resol(x_sed[i1:i2],y_sed[i1:i2],x_r/(1.+z[i]))
        if madau<>'no': ys_z=etau_madau(x_r,z[i])*ys_z
        f[i]=N.trapz(ys_z*r,x_r)*const        
    if nz==1: return f[0]
    else: return f


def of_z_sed(sed,filtro,z=array([0.]),ccd='yes',units='lambda',madau='yes'):
    """Returns array f with f_lambda(z) or f_nu(z) through a given filter 
       Takes into account intergalactic extinction. 
       Flux normalization at each redshift is arbitrary 
    """
    if type(z)==type(0.): z=array([z])

    #Figure out the correct names
    if sed[-4:]<>'.sed':sed=sed+'.sed'
    sed=sed_dir+sed
    if filtro[-4:]<>'.res':filtro=filtro+'.res'
    filtro=fil_dir+filtro
    
    #Get the data
    x_sed,y_sed=get_data(sed,range(2))
    nsed=len(x_sed)
    x_res,y_res=get_data(filtro,range(2))
    nres=len(x_res)

    #Define the limits of interest in wl
    i1=searchsorted(x_sed,x_res[0])-1
    i1=maximum(i1,0)
    i2=searchsorted(x_sed,x_res[-1])+1
    i2=minimum(i2,nsed-1)
    if x_sed[-1]<x_res[-1]: #The SED does not cover the whole filter interval
        #Linear extrapolation of the flux using the last 4 points
        #slope=mean(y_sed[-4:]/x_sed[-4:])
        d_extrap=(x_sed[-1]-x_sed[0])/len(x_sed)
        x_extrap=arange(x_sed[-1]+d_extrap,x_res[-1]+d_extrap,d_extrap)
        extrap=lsq(x_sed[-5:],y_sed[-5:])
        y_extrap=extrap.fit(x_extrap)
        y_extrap=clip(y_extrap,0.,max(y_sed[-5:]))
        x_sed=concatenate((x_sed,x_extrap))
        y_sed=concatenate((y_sed,y_extrap))
        i2=len(y_sed)-1
    r=match_resol(x_res,y_res,x_sed[i1:i2])
    r=where(less(r,0.),0.,r) #Transmission must be >=0
    nz=len(z)
    f=zeros(nz)*1.
    for i in range(nz):
	ys_z=match_resol(x_sed,y_sed,x_sed/(1.+z[i]))
	if madau<>'no': ys_z[i1:i2]=etau_madau(x_sed[i1:i2],z[i])*ys_z[i1:i2]
	f[i]=flux(x_sed[i1:i2],ys_z[i1:i2],r,ccd,units)
    if nz==1: return f[0]
    else: return f

f_z_sed=lf_z_sed
#f_z_sed=nf_z_sed
#f_z_sed=of_z_sed

def f_z_sed_AB(sed,filtro,zz=array([0.]),units='lambda',dz_ab=dz_ab,newAB=0):
    # It assumes ccd=yes,madau=yes by default 
    #AB (checks that the AB file is all right and generates a new one if not)
    ab_file=ABflux(sed,filtro,dz_ab=dz_ab,newAB=newAB)    
    lp=pivotal_wl(filtro)
    f_ab=get_data(ab_file,1)
    fnu=U.match_resol(z_ab,f_ab,zz)
    if units=='nu':      return fnu
    elif units=='lambda': return fnu/lp**2*clight_AHz
    else:
        print 'Units not valid'

def vABflux(sed,filtro,madau=1,plots=0,dz_ab=dz_ab,sed_dir=sed_dir,ab_dir=ab_dir,ccd="yes",units="nu",max_sampling=500):
    """
    Calculates a AB file like the ones used by bpz
    It will set to zero all fluxes
    which are ab_clip times smaller than the maximum flux.
    This eliminates residual flux which gives absurd
    colors at very high-z
    """

    #Figure out the correct names
    if "/" in sed: sed=sed.split("/")[-1]
    if sed[-4:]<>'.sed':sed=sed+'.sed'
    sedroot=sed[:-4]

    if "/" in filtro: filtro=filtro.split("/")[-1]
    if filtro[-4:]<>'.res':filtro=filtro+'.res'
    filtroroot=filtro[:-4]

    # Get the data
    xf,yf=get_filter(filtro,fil_dir=fil_dir)
    xs,ys=get_sed(sed,sed_dir=sed_dir)
    xs,ys=multicompress(less_equal(xs,xf[-1])*greater_equal(xs,xf[0]/(1.+zmax_ab)),(xs,ys))
    
    # Resample to have enough points
    nspectral=sum(less_equal(xs,xf[-1])*greater_equal(xs,xf[0]))
    if nspectral< max_sampling:
        sampling=nspectral
    else:
        sampling=max_sampling

    yf=match_resol(xf,yf,linspace(xf[0],xf[-1],sampling))
    xf=linspace(xf[0],xf[-1],sampling)

    # Operations necessary for normalization and ccd effects
    if ccd=='yes': yf*=xf
    if units=='nu': 
        const=1./N.trapz(yf/xf/xf,xf)/clight_AHz
    else: 
        const=1./N.trapz(yf,xf)

    yf*=const

    xfz=outer(xf,1./(1.+z_ab))
    yfz=outer(yf,0.*z_ab+1.)

    maxy=yfz.max()
    yfz=clip(yfz,0.,maxy)
    if madau:
        a=vetau_madau()
        yfz*=a.e(xf,z_ab)

    yfz=resize(match_resol(xs,ys,ravel(xfz))*ravel(yfz),xfz.shape)

    fAB=N.trapz(yfz,xf,axis=0)
    ABoutput=ab_dir+sedroot+"."+filtroroot+".AB"
    print 'Writing AB file ',ABoutput
    put_data(ABoutput,(z_ab,fAB))
    #return z_ab,fAB


def get_AB(ab_file,madau=1,dz_ab=dz_ab,sed_dir=sed_dir,ab_dir=ab_dir,ccd="yes",units="nu",newAB=0):
    #Returns z,f_AB, creates file if necessary
    sed=ab_file.split(".")[0]+".sed"
    filtro=ab_file.split(".")[1]+".res"
    if newAB:
        ABflux(sed,filtro,madau=madau,dz_ab=dz_ab,sed_dir=sed_dir,
               ab_dir=ab_dir,ccd=ccd,units=units)
    else:
        try:
            d=get_data(ab_dir+ab_file,(0,1))
        except:
            ABflux(sed,filtro,madau=madau,dz_ab=dz_ab,sed_dir=sed_dir,
                   ab_dir=ab_dir,ccd=ccd,units=units)
    d=get_data(ab_dir+ab_file,(0,1))
    return d


def ABflux_M(seds,filters,madau=1,plots=0,dz_ab=dz_ab,
             sed_dir=sed_dir,ab_dir=ab_dir,ccd="yes",units="nu",newAB=0):
    # Calculate AB fluxes for a set of filters and SEDs using multiprocessing
    if not mp:
        print "Multiprocessing not available"
        print "Using sequential processing"
        for filtro in filters:
            print filtro
            for s in seds:
                print s
                ABflux(s,filtro,madau=madau,plots=plots,dz_ab=dz_ab,sed_dir=sed_dir,ab_dir=ab_dir,ccd="yes",units="nu",newAB=newAB)
    else:
        print filters
        p=MP.Pool()
        if not newAB:
            p.map_async(ABtuple,[(s,filtro) for s in seds for filtro in filters])
        else:
            print "newAB",newAB
            p.map_async(ABtuple,[(s,filtro,newAB) for s in seds for filtro in filters])
        p.close()
        p.join()

def ABflux(sed,filtro,madau=1,plots=0,dz_ab=dz_ab,sed_dir=sed_dir,
           ab_dir=ab_dir,ccd="yes",units="nu",newAB=0):
    """
    Calculates a AB file like the ones used by bpz
    It will set to zero all fluxes which are ab_clip times smaller than the maximum flux.
    This eliminates residual flux which gives absurd colors at very high-z
    """
    #Figure out the correct names
    if sed[-4:]<>'.sed':sed=sed+'.sed'
    sedroot=sed[:-4]
    sed=sed_dir+sed
    if filtro[-4:]<>'.res':filtro=filtro+'.res'
    filtroroot=filtro[:-4]
    filtro=fil_dir+filtro

    ABoutput=ab_dir+sedroot+"."+filtroroot+".AB"
    if os.path.exists(ABoutput) and newAB==0:
        try:
            d1,d2=get_data(ABoutput,(0,1))
            if len(d1)==len(z_ab):
                return ABoutput
        except:
            pass
        
        print "%s already exists, but it seems to be corrupted" % ABoutput
        print "Generating file again"
                    
    #Get the data
    x_sed,y_sed=get_data(sed,range(2))
    nsed=len(x_sed)
    x_res,y_res=get_data(filtro,range(2))
    nres=len(x_res)

    if not ascend(x_sed):
        print
        print 'Warning!!!'
        print 'The wavelenghts in %s are not properly ordered' % sed
        print 'They should start with the shortest lambda and end with the longest'        
        print 'This will probably crash the program'

    if not ascend(x_res):
        print
        print 'Warning!!!'
        print 'The wavelenghts in %s are not properly ordered' % filtro
        print 'They should start with the shortest lambda and end with the longest'
        print 'This will probably crash the program'

    if x_sed[-1]<x_res[-1]: #The SED does not cover the whole filter interval
        print 'Extrapolating the spectrum'
        #Linear extrapolation of the flux using the last 4 points
        #slope=mean(y_sed[-4:]/x_sed[-4:])
        d_extrap=(x_sed[-1]-x_sed[0])/len(x_sed)
        x_extrap=arange(x_sed[-1]+d_extrap,x_res[-1]+d_extrap,d_extrap)
        extrap=lsq(x_sed[-5:],y_sed[-5:])
        y_extrap=extrap.fit(x_extrap)
        y_extrap=clip(y_extrap,0.,max(y_sed[-5:]))
        x_sed=concatenate((x_sed,x_extrap))
        y_sed=concatenate((y_sed,y_extrap))
        #plot(x_sed,y_sed,
        #     x_res,y_res);show()
        #connect(x_sed,y_sed)
        #connect(x_res,y_res)

    # Wavelenght range of interest as a function of z_ab
    wl_1=x_res[0]/(1.+z_ab)
    wl_2=x_res[-1]/(1.+z_ab)

    n1=clip(searchsorted(x_sed,wl_1)-1,0,100000000)
    n2=clip(searchsorted(x_sed,wl_2)+1,0,nsed-1)

    # Change filter resolution to SED resolution
    delta_sed=median(x_sed[1:]-x_sed[:-1])
    g=less_equal(x_sed,x_res[-1]+delta_sed)*greater_equal(x_sed,x_res[0]-delta_sed)
    x_r=x_sed[g]
    r=clip(match_resol(x_res,y_res,x_r),0.,1e10)
    # plot(x_res,y_res,x_r,r,"s")
    # show()

    #Operations necessary for normalization and ccd effects
    if ccd=='yes': r=r*x_r
    norm_r=N.trapz(r,x_r)
    if units=='nu': 
        const=norm_r/N.trapz(r/x_r/x_r,x_r)/clight_AHz
    else: const=1.

    const=const/norm_r
    
    nz_ab=len(z_ab)
    f=zeros(nz_ab)*1.
    w=0
    for i in range(nz_ab):
        i1,i2=n1[i],n2[i]
        try:
            ys_z=match_resol(x_sed[i1:i2],y_sed[i1:i2],x_r/(1.+z_ab[i]))
            if madau: ys_z=etau_madau(x_r,z_ab[i])*ys_z
            f[i]=N.trapz(ys_z*r,x_r)*const        
        except:
            if plots:
                plot(x_sed[i1:i2],y_sed[i1:i2],
                     x_r/(1.+z_ab[i]),y_sed[0]*1.+x_r*0.)
                show()
            if w==0:
                print "Warning, setting f[i]==0 for z_ab="
                w=1
            print z_ab[i],
            f[i]=0.
            
    print 'Writing AB file ',ABoutput
    put_data(ABoutput,(z_ab,f))
    return ABoutput


#Useful for multiprocessing
def ABtuple(sf,madau=1,plots=0,dz_ab=dz_ab,sed_dir=sed_dir,ab_dir=ab_dir,ccd="yes",units="nu"):
    if len(sf)==2:
        ABflux(sf[0],sf[1],madau=madau,plots=plots,dz_ab=dz_ab,sed_dir=sed_dir,ab_dir=ab_dir,ccd=ccd,units=units)
    else:
        print "newAB",sf[2]
        ABflux(sf[0],sf[1],madau=madau,plots=plots,dz_ab=dz_ab,sed_dir=sed_dir,ab_dir=ab_dir,ccd=ccd,units=units,newAB=sf[2])        

def VegatoAB(m_vega,filtro,Vega=Vega):
    cons=AB(f_z_sed(Vega,filtro,z=0.,units='nu',ccd='yes'))
    return m_vega+cons
    
def ABtoVega(m_ab,filtro,Vega=Vega):
    cons=AB(f_z_sed(Vega,filtro,z=0.,units='nu',ccd='yes'))
    return m_ab-cons

#Photometric redshift functions

def vlikelihood(f,ef,ft_z):
    """
    Usage: ps[:nz,:nt,:ng]=vlikelihood(f[:nf,:ng],ef[:nf,:ng],ft_z[:nz,:nt,:nf])
    """
    
    global minchi2
    nz=ft_z.shape[0]
    nt=ft_z.shape[1]
    ng=f.shape[-1]

    #Define auxiliary arrays
    chi2=zeros((nz,nt,ng),'float')
    ftt=zeros((nz,nt,ng),'float')
    fot=zeros((nz,nt,ng),'float')
   
    ief2=1./(ef*ef) #ief2=ief2[:nf,:ng]
    foo=add.reduce(f*f*ief2,0) # foo=foo[:ng]
    factor=ft_z[:nz,:nt,:nf,newaxis]*ief2[newaxis,newaxis,:nf,:ng]
    #factor=factor[:nz,:nt,:nf,:ng]

    ftt[:nz,:nt,:ng]=add.reduce(ft_z[:nz,:nt,:nf,newaxis]*factor,-2)
    fot[:nz,:nt,:ng]=add.reduce(f[newaxis,newaxis,:nf,:ng]*factor,-2)
    chi2[:nz,:nt,:ng]=foo[newaxis,newaxis,:ng]-power(fot[:nz,:nt,newaxis],2)/ftt[:nz,:nt,newaxis]
   
    min_chi2=min(chi2)
    minchi2=min(min_chi2)
#   chi2=chi2-minchi2
    chi2=clip(chi2,0.,-2.*etiny)
   
    p=where(greater_equal(chi2,-2.*etiny),0.,exp(-chi2/2.))
    
    norm=add.reduce(add.reduce(add.reduce(p)))
    return p/norm


def likelihood(f,ef,ft_z):
   """ 
   Usage: ps[:nz,:nt]=likelihood(f[:nf],ef[:nf],ft_z[:nz,:nt,:nf])
   """
   global minchi2
   axis=ft_z.shape
   nz=axis[0]
   nt=axis[1]
   
   chi2=zeros((nz,nt),'float')
   ftt=zeros((nz,nt),'float')
   fgt=zeros((nz,nt),'float')
   
   ief2=1./(ef*ef)
   fgg=add.reduce(f*f*ief2)
   factor=ft_z[:nz,:nt,:]*ief2
   
   ftt[:nz,:nt]=add.reduce(ft_z[:nz,:nt,:]*factor,-1)
   fgt[:nz,:nt]=add.reduce(f[:]*factor,-1)
   chi2[:nz,:nt]=fgg-power(fgt[:nz,:nt],2)/ftt[:nz,:nt]
   
   minchi2=min(min(chi2))
#   chi2=chi2-minchi2
   #chi2=clip(chi2,0.,-2.*etiny)
   
   p=where(greater_equal(chi2,-2.*etiny),0.,exp(-chi2/2.))

   norm=add.reduce(add.reduce(p))
   return p/norm

def new_likelihood(f,ef,ft_z):
    """ 
    Usage: ps[:nz,:nt]=likelihood(f[:nf],ef[:nf],ft_z[:nz,:nt,:nf])
    """
    global minchi2
    rolex=reloj()
    rolex.set()
    nz,nt,nf=ft_z.shape

    foo=add.reduce((f/ef)**2)

    fgt=add.reduce(
	f[newaxis,newaxis,:nf]*ft_z[:nz,:nt,:nf]/ef[newaxis,newaxis,:nf]**2
	,-1)

    ftt=add.reduce(
	ft_z[:nz,:nt,:nf]*ft_z[:nz,:nt,:nf]/ef[newaxis,newaxis,:nf]**2
	,-1)

    ao=fgt/ftt
#    print mean(ao),std(ao)

    chi2=foo-fgt**2/ftt+(1.-ao)**2*ftt

    minchi2=min(min(chi2))
    chi2=chi2-minchi2
    chi2=clip(chi2,0.,-2.*etiny)
    p=exp(-chi2/2.)
    norm=add.reduce(add.reduce(p))
    return p/norm

#class p_c_z_t:
#    def __init__(self,f,ef,ft_z):
#	self.nz,self.nt,self.nf=ft_z.shape
#	self.foo=add.reduce((f/ef)**2)
#	self.fgt=add.reduce(
#	    f[NewAxis,NewAxis,:]*ft_z[:,:,:]/ef[NewAxis,NewAxis,:]**2
#	    ,-1)
#	self.ftt=add.reduce(
#	    ft_z[:,:,:]*ft_z[:,:,:]/ef[NewAxis,NewAxis,:]**2
#	    ,-1)
#        #When all the model fluxes are equal to zero
#        self.chi2=self.foo-(self.fgt**2+1e-100)/(self.ftt+1e-100)
#	self.chi2_minima=loc2d(self.chi2[:self.nz,:self.nt],'min')
#	self.i_z_ml=self.chi2_minima[0]
#	self.i_t_ml=self.chi2_minima[1]
#	self.min_chi2=self.chi2[self.i_z_ml,self.i_t_ml]
#	self.likelihood=exp(-0.5*clip((self.chi2-self.min_chi2),0.,1400.))
#        self.likelihood=where(equal(self.chi2,1400.),0.,self.likelihood)
#        #Add the f_tt^-1/2 multiplicative factor in the exponential
#        self.chi2+=-0.5*log(self.ftt+1e-100)
#        min_chi2=min(min(self.chi2))
#	self.Bayes_likelihood=exp(-0.5*clip((self.chi2-min_chi2),0.,1400.))
#        self.Bayes_likelihood=where(equal(self.chi2,1400.),0.,self.Bayes_likelihood)
#        
#        #plo=FramedPlot()
#        #for i in range(self.ftt.shape[1]):
#        #    norm=sqrt(max(self.ftt[:,i]))
#        #    # plo.add(Curve(arange(self.ftt.shape[0]),self.ftt[:,i]**(0.5)))
#        #    plo.add(Curve(arange(self.ftt.shape[0]),self.likelihood[:,i],color='red'))
#        #    plo.add(Curve(arange(self.ftt.shape[0]),self.likelihood[:,i]*sqrt(self.ftt[:,i])/norm))
#        #plo.show()#
#
#    def bayes_likelihood(self):
#        return self.Bayes_likelihood


class p_c_z_t_old:
    def __init__(self,f,ef,ft_z):
        self.nz,self.nt,self.nf=ft_z.shape
        
        # Define likelihood quantities taking into account non-observed objects

        self.foo=add.reduce((f/ef)**2)
	self.fot=add.reduce(f[newaxis,newaxis,:]*ft_z[:,:,:]/ef[newaxis,newaxis,:]**2,-1)
        self.ftt=add.reduce(ft_z[:,:,:]*ft_z[:,:,:]/ef[newaxis,newaxis,:]**2,-1)

        #If ftt==0, set value of chi2 to maximum to avoid overflows
        self.chi2=where(equal(self.ftt,0.),
                        1./tiny,
                        self.foo-self.fot**2/self.ftt)
	self.chi2_minima=loc2d(self.chi2[:self.nz,:self.nt],'min')
 	self.i_z_ml=self.chi2_minima[0]
	self.i_t_ml=self.chi2_minima[1]
	self.min_chi2=self.chi2[self.i_z_ml,self.i_t_ml]
        self.likelihood=exp(-0.5*(self.chi2-self.min_chi2))
        
    def plots(self):        
        # Normalize and collapse likelihoods (without prior)
        for it in range(nt):
            plot(arange(self.nz),self.likelihood[:,it])
        show()
        ask('More?')


class p_c_z_t:
    def __init__(self,f,ef,ft_z):
        self.nz,self.nt,self.nf=ft_z.shape
        
        # Define likelihood quantities taking into account non-observed objects

        self.foo=(f/ef).sum()**2
        #self.foo=add.reduce((f/ef)**2)
	self.fot=add.reduce(f[newaxis,newaxis,:]*ft_z[:,:,:]/ef[newaxis,newaxis,:]**2,-1)
        self.ftt=add.reduce(ft_z[:,:,:]*ft_z[:,:,:]/ef[newaxis,newaxis,:]**2,-1)
        self.am=self.fot/self.ftt


        #If ftt==0, set value of chi2 to maximum to avoid overflows
        self.chi2=where(equal(self.ftt,0.),
                        1./tiny,
                        self.foo-self.fot**2/self.ftt)

	self.chi2_minima=loc2d(self.chi2[:self.nz,:self.nt],'min')
 	self.i_z_ml=self.chi2_minima[0]
	self.i_t_ml=self.chi2_minima[1]
	self.min_chi2=self.chi2[self.i_z_ml,self.i_t_ml]
        self.likelihood=exp(-0.5*(self.chi2-self.min_chi2))

        # plot(1./sqrt(self.ftt)*(1.+U.erf(self.fot/sqrt(2.*self.ftt)))*self.likelihood)
        # plot(self.likelihood,"--")
        # show()

        
    def plots(self):        
        # Normalize and collapse likelihoods (without prior)
        for it in range(nt):
            plot(arange(self.nz),self.likelihood[:,it])
        show()
        ask('More?')



class p_c_z_t_F:
    def __init__(self,f,ef,ft_z,dz,z):
        self.nz,self.nt,self.nf=ft_z.shape
        step=int(0.01/dz)
        #Define likelihood quantities taking into account non-observed objects
        self.foo=add.reduce(where(less(f/ef,1e-4),0.,(f/ef)**2))
        nonobs=less(f[newaxis,newaxis,:]/ef[newaxis,newaxis,:]+ft_z[::step,:,:]*0.,1e-4)
	self.fot=add.reduce(
            where(nonobs,0.,f[newaxis,newaxis,:]*ft_z[::step,:,:]/ef[newaxis,newaxis,:]**2)
            ,-1)
        self.ftt=add.reduce(
            where(nonobs,0.,ft_z[::step,:,:]*ft_z[::step,:,:]/ef[newaxis,newaxis,:]**2)
            ,-1)
    
        self.chi2=where(equal(self.ftt,0.),
                        1./tiny,
                        self.foo-self.fot**2/self.ftt)
	self.chi2_minima=loc2d(self.chi2[:self.nz,:self.nt],'min')
	self.i_z_ml=self.chi2_minima[0]
	self.i_t_ml=self.chi2_minima[1]
	self.min_chi2=self.chi2[self.i_z_ml,self.i_t_ml]
	self.i_z_ml*=step
	self.i_t_ml*=step
        self.likelihood=self.chi2*0.
        for jt in range(self.nt):
            self.likelihood[:,jt]=match_resol(z[::step],exp(-0.5*(self.chi2-self.min_chi2))[:,jt],z)

        
    def plots(self):        
        # Normalize and collapse likelihoods (without prior)
        for it in range(nt):
            plot(arange(self.nz),self.likelihood[:,it])
        show()
        ask('More?')


class p_c_z_t_2:
    def __init__(self,f,ef,ft_z):
        self.nz,self.nt,self.nf=ft_z.shape
        #Define likelihood quantities taking into account non-observed objects
        #nonobs=less(f[newaxis,newaxis,:]/ef[newaxis,newaxis,:]+ft_z[:,:,:]*0.,1e-4)
        #self.foo=add.reduce(where(less(f/ef,1e-4),0.,(f/ef)**2))
        #nonobs=less(f[newaxis,newaxis,:]/ef[newaxis,newaxis,:]+ft_z[:,:,:]*0.,1e-4)
	self.fot=add.reduce(f[newaxis,newaxis,:]*ft_z[:,:,:]/ef[newaxis,newaxis,:]**2,-1)
        self.ftt=add.reduce((ft_z[:,:,:]/ef[newaxis,newaxis,:])**2,-1)
        self.foo=add.reduce((f/ef)**2)            
        #Define chi2 adding eps to the ftt denominator to avoid overflows
        self.chi2=where(equal(self.ftt,0.),
                        self.foo,
                        self.foo-(self.fot**2)/(self.ftt+tiny))
	self.chi2_minima=loc2d(self.chi2[:self.nz,:self.nt],'min')
	self.i_z_ml=self.chi2_minima[0]
	self.i_t_ml=self.chi2_minima[1]
	self.min_chi2=self.chi2[self.i_z_ml,self.i_t_ml]
        self.likelihood=exp(-0.5*(self.chi2-self.min_chi2))

        #Now we add the Bayesian f_tt^-1/2 multiplicative factor to the exponential
        #(we don't multiply it by 0.5 since it is done below together with the chi^2
        #To deal with zero values of ftt we again add an epsilon value.
        #self.expo=where(
        #    equal(self.ftt,0.),
        #    self.chi2,
        #    self.chi2+log(self.ftt+eps)
        #    )
        #Renormalize the exponent to preserve dynamical range
        #self.expo_minima=loc2d(self.expo,'min')
        #self.min_expo=self.expo[self.expo_minima[0],self.expo_minima[1]]
        #self.expo-=self.min_expo
        #self.expo=clip(self.expo,0.,-2.*etiny)
        #Clip very low values of the probability
        #self.Bayes_likelihood=where(
        #    equal(self.expo,-2.*etiny),
        #    0.,
        #    exp(-0.5*self.expo))

    #def bayes_likelihood(self):
    #    return self.Bayes_likelihood
        
    #def various_plots(self):        
    #    #Normalize and collapse likelihoods (without prior)
    #    norm=add.reduce(add.reduce(self.Bayes_likelihood))
    #    bl=add.reduce(self.Bayes_likelihood/norm,-1)
    #    norm=add.reduce(add.reduce(self.likelihood))
    #    l=add.reduce(self.likelihood/norm,-1)
    #    plo=FramedPlot()
    #    plo.add(Curve(arange(self.nz),bl,color='blue'))
    #    plo.add(Curve(arange(self.nz),l,color='red'))
    #    plo.show()


        #plo2=FramedPlot()
        #for i in range(self.ftt.shape[1]):
        #for i in range(2):
        #    #plo2.add(Curve(arange(self.nz),log(self.fot[:,i]*self.fot[:,i])))
        #    plo2.add(Curve(arange(self.nz),log(self.ftt[:,i])))
        #plo2.show()


        #for i in range(self.ftt.shape[1]):
        #for i in range(2):
        #    plo2.add(Curve(arange(self.nz),-0.5*(self.fot[:,i]*self.fot[:,i]/self.ftt[:,i]+log(self.ftt[:,i]))))
        #    plo2.add(Curve(arange(self.nz),-0.5*(self.fot[:,i]*self.fot[:,i]/self.ftt[:,i]),color='red'))
        #plo2.show()

        #plo3=FramedPlot()
        #for i in range(self.ftt.shape[1]):
        #    norm=sqrt(max(self.ftt[:,i]))
        #    plo3.add(Curve(arange(self.nz),self.fot[:,i]*self.fot[:,i]/self.ftt[:,i]))
        #plo3.show()

class p_c_z_t_color:
    def __init__(self,f,ef,ft_z):
	self.nz,self.nt,self.nf=ft_z.shape
	self.chi2=add.reduce(
	((f[newAxis,newAxis,:]-ft_z[:,:,:])/ef[newAxis,newAxis,:])**2
	,-1)
	self.chi2_minima=loc2d(self.chi2[:self.nz,:self.nt],'min')
	self.i_z_ml=self.chi2_minima[0]
	self.i_t_ml=self.chi2_minima[1]
	self.min_chi2=self.chi2[self.i_z_ml,self.i_t_ml]
	self.likelihood=exp(-0.5*clip((self.chi2-self.min_chi2),0.,1400.))
    def bayes_likelihood(self):
	return self.likelihood

#def gr_likelihood(f,ef,ft_z):
#    #Color-redshift Likelihood a la Rychards et al. (SDSS QSOs)
#    global minchi2
#    nf=f.shape[0]
#    nz=ft_z.shape[0]
#    nt=ft_z.shape[1]
#    print f,ef,ft_z[:10,0,:]
#    chi2=add.reduce(
#	((f[NewAxis,NewAxis,:nf]-ft_z[:nz,:nt,:nf])/ef[NewAxis,NewAxis,:nf])**2
#	,-1)
#    minchi2=min(min(chi2))
#    chi2=chi2-minchi2
#    chi2=clip(chi2,0.,1400.)
#    p=exp(-chi2/2.)
#    norm=add.reduce(add.reduce(p))
#    return p/norm

#def p_and_minchi2(f,ef,ft_z):
#    p=gr_likelihood(f,ef,ft_z)
#    return p,minchi2

#def new_p_and_minchi2(f,ef,ct):
#    p=color_likelihood(f,ef,ct)
#    return p,minchi2


def prior3D(z,m,name="UgCd_eB11",nt=11,n_interp=0,m_step=0.1,new=0):
    """
    Given the magnitude m, produces the prior  p(z|T,m)
    Usage: pi[:nz,:nt]=prior(z[:nz],m,name='eB11',nt)
    """    
    global prior_d
    try:
        prior_d.keys()
        hay_prior=1
    except:
        hay_prior=0

    if not hay_prior or new:
	prior_dict={}
        npzfile=bpz_dir+"/"+name+"_prior.npz"
        prior_d=N.load(npzfile)

    accuracy=str(len(str(int(1./m_step)))-1)#number of decimals kept
    forma='%.'+accuracy+'f'        
    m_dict=forma %m    
    
    if not prior_dict.has_key(m_dict):
        it=linspace(1.,nt,n_interp*(nt-1)+nt)
        im=searchsorted(prior_d["xm"],m)
        if im==0 or im<0:
            data=prior_d["pmzt"][0,:,:]
        if im>=len(prior_d["xm"])-1:
            data=prior_d["pmzt"][-1,:,:]
        else:
            #Interpolate across closest magnitudes 
            w1=(float(m_dict)-prior_d["xm"][im])
            w2=(prior_d["xm"][im+1]-float(m_dict))
            data=(w1*prior_d["pmzt"][im,:,:]+w2*prior_d["pmzt"][im+1,:,:])/(w1+w2)
            #Apply final interpolation
        prior_dict[m_dict]=interpolate.RectBivariateSpline(prior_d["xz"],prior_d["xt"],data,kx=2,ky=2,s=0)(z,it)
    return prior_dict[m_dict]

def prior(z,m,name='hdfn',nt=6,n_interp=0,x=None,y=None,filter_m="HST_ACS_WFC_F814W",lib="eB11.list",m_step=0.1,new=0):
    """
    Given the magnitude m, produces the prior  p(z|T,m)
    Usage: pi[:nz,:nt]=prior(z[:nz],m,name='hdfn',nt)
    """
    global prior_dict
    try:
        len(prior_dict.keys())
    except:
        prior_dict={}

    if "eB11" in name:
        return prior3D(z,m,name,nt=nt,n_interp=n_interp,m_step=m_step,new=new)   
    elif name in ['none','flat','NONE','None','FLAT']:
        return
    elif name[:2]=="SM":
        import cosmology as C
    else:
        exec('import prior_%s as prior_%s' % (name,name))
            
    #We estimate the priors at m_step intervals
    #and keep them in a dictionary.
    m_dict=round(m/m_step)*m_step
    print 'm_dict', m_dict

    if not prior_dict.has_key(m_dict) or new or name=='lensing': #if lensing, the magnitude alone is not enough
        if name[:2]=="SM":
            import cosmology as C
            if name[:2]=="SM":
                letras=list(name)
                use_alpha,use_phi=0,0    
            if "A" in letras:
                use_alpha=1
            elif "F" in letras:
                use_phi=1
            #Use a very faint magnitude for those cases in which there is no reasonable magnitude measurement
            if m_dict<14.:
                m=32.
            elif m_dict>32:
                m=32
            else:
                m=m_dict
                #[::10]
            prior0=C.nzt_mass(z,arange(nt),m,filter_m=filter_m,lib=lib.split("/")[-1],
                              area=1.,use_phi=use_phi,use_alpha=use_alpha)
            print 'prior0', np.shape(prior0)
        elif name=='lensing':
            exec("prior0=prior_%s.function(z,m_dict,nt,x,y)" % name)
        else:
            exec("prior0=prior_%s.function(z,m_dict,nt)" % name)
        if n_interp:
            if name[:2]=="SM":
                #[::10]
                prior_dict[m_dict]=interpolate.RectBivariateSpline(z,arange(nt),prior0[:,:],kx=1,ky=1)(z,linspace(0,nt-1,nt+(nt-1)*n_interp))
            else:
                prior_dict[m_dict]=interpolate.interp1d(arange(nt),prior0[:,:],kind="slinear",axis=1
                                                        )(linspace(0,nt-1,nt+(nt-1)*n_interp))
        else:
            if name[:2]=="SM":
                prior_dict[m_dict]=prior0
            else:
                prior_dict[m_dict]=prior0
            
        # if n_interp:
	#    pp_i=prior_dict[m_dict]
	#    nz=pp_i.shape[0]
	#    nt=pp_i.shape[1]
	#    nti=nt+(nt-1)*int(n_interp)
	#    tipos=arange(nt)*1.
	#    itipos=arange(nti)*1./(1.+float(n_interp))
	#    buffer=zeros((nz,nti))*1.
	#    for iz in range(nz):
        #		buffer[iz,:]=match_resol(tipos,pp_i[iz,:],itipos)
	#    prior_dict[m_dict]=buffer

    return prior_dict[m_dict]

def interval(p,x,ci=.99):
    """Gives the limits of the confidence interval
       enclosing ci of the total probability
       i1,i2=limits(p,0.99)
    """
    q1=(1.-ci)/2.
    q2=1.-q1
    cp=add.accumulate(p)
    if cp[-1]<>1.: cp=cp/cp[-1]
    i1=searchsorted(cp,q1)-1
    i2=searchsorted(cp,q2)
    i2=minimum(i2,len(p)-1)
    i1=maximum(i1,0)
    return x[i1],x[i2] 
 
def odds1(p,x,x1,x2):
    """Estimate the fraction of the total probability p(x) enclosed by the interval x1,x2"""
    cp=add.accumulate(p)
    i1=searchsorted(x,x1)-1
    i2=searchsorted(x,x2)
    if i1<0:
        return cp[i2]/cp[-1]
    if i2>len(x)-1:
        return 1.-cp[i1]/cp[-1]
    return (cp[i2]-cp[i1])/cp[-1]

def odds2(p,x,x1,x2):
    """Estimate the fraction of the total probability p(x)
    enclosed by the interval x1,x2"""
    # New way of calculating the odds
    cp=add.accumulate(p)
    #cp/=cp[-1]
    x1=max(x1,x[0])
    x2=min(x2,x[-1])
    o=match_resol(x,cp,(x1,x2))
    return (o[1]-o[0])/(cp[-1]-cp[0])


def odds3(p,x,x1,x2):
    """Estimate the fraction of the total probability p(x)
    enclosed by the interval x1,x2"""
    # New way of calculating the odds
    #nx=arange(x[0],x[-1],x[1]-x[0])
    #p=clip(match_resol(x,p,nx),0.,1e10)
    #x=nx
    cp=add.accumulate(p)
    # cp/=cp[-1]
    x1=max(x1,x[0])
    x2=min(x2,x[-1])
    o=match_resol(x,cp,(x1,x2))
    return (o[1]-o[0])/(cp[-1]-cp[0])

odds=odds1

class p_bayes:
    #This class reads the information contained in the files produced by BPZ
    #when the option -PROBS_LITE is on
    def __init__(self,file):
        self.file=file
        dummy=get_2Darray(file)
        self.id_list=map(int,list(dummy[:,0]))
        self.p=dummy[:,1:]
        del(dummy)
        header=get_header(file)
        #header=split(header,'(')[2]
        #header=split(header,')')[0]
        #zmin,zmax,dz=map(float,tuple(split(header,',')))
        header=header.split('(')[2]
        header=header.split(')')[0]
        zmin,zmax,dz=map(float,tuple(header.split(',')))
        self.z=arange(zmin,zmax,dz)
        print zmin,zmax,dz

    def p_j(self,id):
        j=self.id_list.index(int(id))
        print j
        return self.p[j,:]        

    def plot_p(self,id,limits=None):
        if type(id)<>type((1,)):
            try:
                j=self.id_list.index(int(id))
                p_j=self.p[j,:]
                plot(self.z,p_j)
                if limits<>None:
                    axis((limits[0],limits[1],limits[2],limits[3]))
                show()
            except:
                print 'Object %i not in the file %s' % (id,self.file)
            self.prob=p_j/max(p_j)
        else:
            #p=FramedPlot()
            #p.frame1.draw_grid=1
            pall=self.p[0,:]*0.+1.
            pmax=0.
            for i in id:
                try:
                    j=self.id_list.index(int(i))
                    p_j=self.p[j,:]
                    if max(p_j)>pmax: pmax=max(p_j)
                    pall*=p_j
                    plot(self.z,p_j)
                except:
                    print 'Object %i not in the file %s' % (id,self.file)
            plot(self.z,pall/max(pall)*pmax,'r')
            if limits<>None:
                axis((limits[0],limits[1],limits[2],limits[3]))
            show()
            self.prob=pall/max(pall)

    def maxima(self,limits=(0.,6.5)):
        g=greater_equal(self.z,limits[0])*less_equal(self.z,limits[1])
        z,p=multicompress(g,(self.z,self.prob))
        imax=argmax(p)
        xp=add.accumulate(p)
        xp/=xp[-1]
        self.q66=match_resol(xp,z,array([0.17,0.83]))
        self.q90=match_resol(xp,z,array([0.05,0.95]))
        #print self.q66
        #print self.q90
        return z[imax]

                
    #def hist_p(self,dz=0.25):
    #    self.pt=sum(self.p)
    #    self.xz=arange(self.z[0],self.z[-1]+dz,dz)
    #    self.hb=bin_stats(self.z,self.pt,self.xz,'sum')
        
#Misc stuff

def runbpz(inputs="test.cat -ZMAX 10."):
    command="python "+os.path.expanduser("~/Dropbox/bpz/bpz.py ")+inputs
    print command
    os.system(command)


def checkpars(input="-PRIOR B11",benchmark="-PRIOR eB11v2",
              cats=["fireworks_spec.cat","global_spzcal.cat","COSMOSIB_DR2.cat","UDF_zspec.cat","deep2_dr4.cat"],
              cols=["fireworks_spec_B10.columns","global_spzcal.columns","COSMOSIB_DR2_B10.columns","UDF_zspec.columns","deep2_dr4.columns"],
              pars=[' -EXCLUDE "*WFI","IRAC_*","*ISAAC",HST_ACS_WFC_435W ',
                    ' -EXCLUDE F_J,F_KS,F_H,"F_3*","F_4*" ',
                    ' -EXCLUDE "*SDSS","*CFHTLS","IA*","NB*","*UV","*K*",B_SupCam,V_SupCam,gprime_SupCam ',
                    ' -EXCLUDE "*NIC3" ',
                    " "],
              new=1,
              full=0,
              cut=0.5,
              plots=0):
    out1=map(lambda x: x.split(".")[0]+"_"+"_".join(input[1:].split())+".bpz",cats)
    out2=map(lambda x: x.split(".")[0]+"_"+"_".join(benchmark[1:].split())+".bpz",cats)
    if full:
        pars=[" "," "," "," "," "]
        out1=map(lambda x: x.split(".")[0]+"_full.bpz",out1)
        out2=map(lambda x: x.split(".")[0]+"_full.bpz",out2)
    for i in range(len(cats)):
        if not os.path.exists(out1[i]) or new:
            map(lambda x: os.remove(x),glob(cats[i].split(".")[0]+"*.npy"))
            print out1[i]
            runbpz_m("%s %s -COLUMNS %s -OUTPUT %s %s" % (
                    cats[i],input,cols[i],out1[i],pars[i]))
        a=d_stats(out1[i],cut); print out1[i]; print a.nice(),a.n
        xz=linspace(min(a.zs),max(a.zs),50)
        if plots: 
            hist(a.zb,xz,color="b",histtype="step")
        if not os.path.exists(out2[i]) or new:
            print out2[i]
            runbpz_m("%s %s -COLUMNS %s -OUTPUT %s %s" % (
                    cats[i],benchmark,cols[i],out2[i],pars[i]))
        a=d_stats(out2[i],cut); print out2[i]; print a.nice(),a.n
        if plots: 
            hist(a.zb,xz,color="r",histtype="step")
            hist(a.zs,xz,color="k",histtype="step",linewidth=4)
            show()

def nm(m,Dm=1.):
    #Numbers counts in F814W for 1sq-deg
    #From Annis et al. (Stripe 82) and Benitez et al. 2004
    #Good match with COSMOS, excellent match with CFHTLS (taking into account the difference in filters)
    nmb=10.**(0.60*(m-0.19)-8.45)*Dm#F814W mags are fainter than iSDSS
    nmf=10.**(0.33*m-3.29)*Dm*1.25
    return where(nmb>nmf,nmf,nmb)

class check_prior:
    def __init__(self,prior="prior_B13",
                 #cats=["fireworks_spec.cat","global_spzcal.cat","COSMOSIB_DR2.cat","UDF_zspec.cat","deep2_dr4.cat"],
                 #cols=["fireworks_spec_B10.columns","global_spzcal.columns","COSMOSIB_DR2_B10.columns","UDF_zspec.columns","deep2_dr4.columns"],
                 cats=["fireworks_spec.cat","COSMOSIB_DR2.cat","UDF_zspec.cat"],
                 cols=["fireworks_spec_B10.columns","COSMOSIB_DR2_B10.columns","UDF_zspec.columns"],
                 new=0,
                 nt=11):
        P=__import__(prior)
        M=[]
        Z=[]
        T=[]
        out=map(lambda x: x.split(".")[0]+"_OT.bpz",cats)
        for i in range(len(cats)):
            if not os.path.exists(out[i]) or new:
                runbpz_m("%s -PRIOR %s -COLUMNS %s -OUTPUT %s -CACHE yes -ONLY_TYPE yes" % (
                    cats[i],prior.split("_")[1],cols[i],out[i]))
            m,z,t=get_data(out[i],(10,9,4))
            M+=list(m)
            Z+=list(z)
            T+=list(t)
        M=array(M)
        Z=array(Z)
        T=array(map(lambda x: int(round(x-1.)),T))
        xz=arange(0.01,6.,.1)
        nz=len(xz)
        M,Z,T=multicompress(greater(M,0)*less(M,99),
                            (M,Z,T))
        Nz=zeros((nz,nt),"float")
        Pz=zeros((nz,nt),"float")
        for j in range(11):
            for k in range(len(M[T==j])):
                Pz[:,j]+=P.function(xz,M[T==j][k])[:,j]
            Nz[:,j]=histogram_old(Z[T==j],xz)
        Nz/=Nz.mean()
        Pz/=Pz.mean()
        #for j in range(11):
        #    print Pz[:,j].sum(),Nz[:,j].sum()
        #    plot(xz,Pz[:,j],xz,Nz[:,j])
        #    show()
        plot(xz,Pz[:,:5].sum(-1),"r",xz,Nz[:,:5].sum(-1),"r--")
        show()
        plot(xz,Pz[:,5:7].sum(-1),"g",xz,Nz[:,5:7].sum(-1),"g--")
        show()
        plot(xz,Pz[:,7:].sum(-1),"b",xz,Nz[:,7:].sum(-1),"b--")
        show()
        plot(xz,Pz[:,:].sum(-1),"k",xz,Nz[:,:].sum(-1),"k--")
        show()


class runbpz_m:
    def __init__(self,inputs="test.cat -ZMAX 10.",cleanup=1,n=10):
        bits=inputs.split()
        try:
            output=bits[bits.index("-OUTPUT")+1]
        except:
            output=bits[0][:-4]+".bpz"            

        l=array(map(lambda x: x[0]!="#", open(inputs.split()[0],"r").readlines())).sum()
        ll=map(int,linspace(0,l,n+1))
        outfiles=[]
        arguments=[]
        for i in range(n):            
            arguments.append(
                inputs+ " -LINES %i-%i -OUTPUT out%i%i.bpz -CACHE yes" % (ll[i],ll[i+1],ll[i],ll[i+1]))
            if i>0: arguments[-1]+=" -VERBOSE no "
            outfiles.append("out%i%i.bpz" % (ll[i],ll[i+1]))

        runbpz(arguments[0]+" -GET_Z no")
            
        p=MP.Pool()
        p.map_async(runbpz,arguments)
        p.close()
        p.join()

        #Put the output together again
        for i in range(n):
            if i==0:
                bag=loadtxt(outfiles[0])
                header=get_header(outfiles[0])
            else:
                bag=concatenate((bag,loadtxt(outfiles[i])))
        nc=bag.shape[1]
        savetxt(output,bag,fmt="%s "+(nc-1)*"%.5f ")
        #Clean up header
        bits=header.split("\n")
        for i in range(len(bits)):
            if "LINES" in bits[i]: j=i
            if "File" in bits[i]:
                l=bits[i].split()
                l[2]=output
                bits[i]=" ".join(l)
                if "OUTPUT" in bits[i]: bits[i]="##OUTPUT=%s" % output
                
        bits.remove(bits[j])
        put_header(output,"\n".join(bits))
        if cleanup:
            for i in range(n):
                os.remove(outfiles[i])
                
def get_datasex(file,cols,purge=1,mag=(2,99.),emag=(4,.44),flag=(24,4),detcal='none'):
    """
      Usage:
      x,y,mag,emag=get_datasex('file.cat',(0,1,24,12))
      If purge=1, the function returns the corresponding columns 
      of a SExtractor output file, excluding those objects with 
      magnitude <mag[1], magnitude error <=emag[1] and flag <=flag[1]
      mag[0],emag[0] and flag[0] indicate the columns listing 
      these quantities in the file
    """
    if type(cols)==type(0): nvar=1
    else:nvar=len(cols)

    if purge:
	if nvar>1:datos=get_2Darray(file,cols)
	else: datos=get_data(file,cols)
	if detcal=='none': detcal=file
	m,em,f=get_data(detcal,(mag[0],emag[0],flag[0]))
	good=less_equal(f,flag[1])*less_equal(em,emag[1])*less(m,mag[1])
	datos=compress(good,datos,0)
	lista=[]
	if nvar>1:
	    for i in range(datos.shape[1]):lista.append(datos[:,i])
	    return tuple(lista)
	else: return datos
    else:
	return get_data(file,cols)


def cleanforbpz(datos,zp=0.,sn_min=1.,m_lim=None,plots=1):
    #datos = (m1,em1,m2,em2,...)
    nc=len(datos)
    salida=[]

    for i in range(0,nc,2):
        m,em=datos[i:i+2]
        print "len(m),len(em)",len(m),len(em)
        if plots:
            plot(m,em,".")
            axis((18.,30.,0.,1.))
            ask()

        # Define reasonable range for BPZ
        seen=greater(m,6.)*less(m,99.)*greater(em,0.)*less(em,0.75257)
        # Get limiting magnitude
        if m_lim==None:
            ml=get_limitingmagnitude(compress(seen,m),compress(seen,em))        
        
        print "Filter",i/2.
        print "m_lim",ml
        nondetected=seen*greater(m,ml)
        detected=seen*less_equal(m,ml)
        notseen=logical_not(seen+nondetected)

        print "Not detected",nondetected.sum()
        print "Detected",detected.sum()
        print "Not observed",notseen.sum()

        mc=m*0.
        emc=em*0.

        mc=where(nondetected,99.,mc)
        mc=where(notseen,-99.,mc)
        mc=where(detected,m+zp,mc)
        
        emc=where(nondetected,ml+zp,emc)
        emc=where(notseen,0.,emc)
        emc=where(detected,em,emc)

        salida.append(mc)
        salida.append(emc)
        if plots:
            plot(mc,emc,".")
            axis((18.,30.,0.,1.))
            show()
        
    return salida

def sex2bpzmags(f,ef,zp=0.,sn_min=1.,m_lim=None):
    """
    This function converts a pair of flux, error flux measurements from SExtractor
    into a pair of magnitude, magnitude error which conform to BPZ input standards:
    - Nondetections are characterized as mag=99, errormag=m_1sigma
    - Objects with absurd flux/flux error combinations or very large errors are
      characterized as mag=-99 errormag=0.
    """
    nondetected=less_equal(f,0.)*greater(ef,0) #Flux <=0, meaningful phot. error
    nonobserved=less_equal(ef,0.) #Negative errors
    #Clip the flux values to avoid overflows
    f=numpy.clip(f,1e-100,1e10)
    ef=numpy.clip(ef,1e-100,1e10)
    nonobserved+=equal(ef,1e10)
    nondetected+=less_equal(f/ef,sn_min) #Less than sn_min sigma detections: consider non-detections
    
    detected=logical_not(nondetected+nonobserved)
    
    m=zeros(len(f))*1.
    em=zeros(len(ef))*1.
    
    m = where(detected,-2.5*log10(f)+zp,m)
    m = where(nondetected,99.,m)
    m = where(nonobserved,-99.,m)

    em = where(detected,2.5*log10(1.+ef/f),em)
    if not m_lim:
        em = where(nondetected,-2.5*log10(ef)+zp,em)
    else:
        em = where(nondetected,m_lim,em)        
    em = where(nonobserved,0.,em)
    #plot(compress(detected,m),compress(detected,em),".")
    #show()
    return m,em



class d_stats:
    def __init__(self,bpz="bobo.bpz",cut=0.,
                 plots=None,zmin=None,zmax=None,mmin=None,mmax=None,
                 omin=None,omax=None,dmin=None,dmax=None,
                 Dmin=None,Dmax=None,tmin=None,tmax=None,chi2max=None,
                 cols="bpz"):

        if cols=="bpz":
            try: zb,zb1,zb2,tb,zs,o,chi2,m,zml=get_data(bpz,(1,2,3,4,10,5,9,11,7))
            except: zb,zb1,zb2,tb,zs,o,chi2,m,zml=get_data(bpz,(1,2,3,4,9,5,8,10,6))
        else:
            try: zb,zb1,zb2,tb,zs,o,chi2,m,zml=get_data(bpz,colss)
            except: zb,zb1,zb2,tb,zs,o,chi2,m,zml=get_data(bpz,cols)            
        self.ng=len(zb)
        self.nt=int(round(max(tb)))
        so=sort(o)
        no=arange(len(o))*1.
        no*=float(max(no))**-1.
        o_thr=match_resol(no,so,cut)
        if cut==0.: o_thr=0.
        g=o>=o_thr
        d=(zb-zs)/(1.+zs)
        if zmin<>None: g=g*greater(zs,zmin)
        if zmax<>None: g=g*less_equal(zs,zmax)
        if mmax<>None: g=g*less_equal(m,mmax)
        if mmin<>None: g=g*greater(m,mmin)
        if tmax<>None: g=g*less_equal(tb,tmax)
        if tmin<>None: g=g*greater(tb,tmin)
        if omin<>None: g=g*greater_equal(o,omin)
        if omax<>None: g=g*less(o,omax)
        if dmax<>None: g=g*less_equal(abs(d),dmax)
        if dmin<>None: g=g*greater(abs(d),dmin)
        if Dmax<>None: g=g*less_equal(abs(zb-zs),Dmax)
        if Dmin<>None: g=g*greater(abs(d),Dmin)
        if chi2max<>None: g=g*less_equal(chi2,chi2max)
        #d=(zml[g]-zs[g])/(1.+zs[g])
        d=d[g]
        self.tb=tb[g]
        self.o=o[g]
        self.m=m[g]
        self.zs=zs[g]
        self.zb=zb[g]
        self.chi2=chi2[g]
        ds=sort(d)
        ns=arange(len(d))/float(len(d))
        dm1s=match_resol(ns,ds,0.17)
        dp1s=match_resol(ns,ds,0.83)
        d1s=0.5*(dp1s-dm1s)
        self.med=median(d)
        self.sigma=d1s
        self.std_mad=std_mad(d)
        self.std_robust=std_robust(d,3.,5.)
        self.std_phat=std(compress(less_equal(abs(d),0.15),d))
        self.std=std(d)
        self.std_ci=std_ci(d,n_sigma=1.)
        self.dzb=median(concatenate((zb-zb1,zb2-zb)))
        self.d=d    
        self.no=no
        self.so=so
        self.n=len(d)*1.
        self.ni=len(g)*1.
        # Defined in the same sense as BPZ
        self.out3=greater(abs(d),3.*std_mad(d)).sum()/float(len(d))
        self.out5=greater(abs(d),5.*std_mad(d)).sum()/float(len(d))
        self.out15=greater(abs(d),.15).sum()/float(len(d))
        if plots:
            plot(zs[g],zb[g],"o")
            plot(zb[g],zb[g])
            title(bpz)
            xlabel("z_s")
            ylabel("z_B")
            #plot(zs[g],zml[g],"+")
            show()


    def out(self,x,nsigma=5.):
        return greater(abs(self.d),nsigma*x).sum()/float(len(self.d))

    def nice(self,x=0.03,nsigma=5.):
        output="  med    std_mad  std_ci std_phat  std  n>5sigma  n>5.*%.4f num. \n" % x
        output+= 8*" %.4f " % (self.med,self.std_mad,self.std_ci,self.std_phat,self.std,self.out5,self.out(x),self.n)
        return output

    def types(self,lib=None):
        if lib:
            try:
                temps=get_lib(lib)
            except:
                temps=get_str(lib,0)
        ind=[]
        dzp=[]
        rms=[]
        out15=[]
        ht=[]

        output="type   med    std_mad  f>0.15 N_t\n" 
        format="%s "+ 3*" %.4f"+ " %i \n" 
        for it in range(1,self.nt+1):
            ind.append(it-1)
            if lib: ind.append(temps[it-1])
            g=less_equal(self.tb,it+0.5)*greater(self.tb,it-0.5)
            ht.append(g.sum())
            tag=str(ind[-1])
            if g.sum()>0:
                dzp.append(median(self.d[g]))            
                rms.append(std_mad(self.d[g]))
                out15.append(greater(abs(self.d[g]),0.15).sum()/float(g.sum()))
            else:
                dzp.append(0.)            
                rms.append(0.)
                out15.append(0.)
            output+=format % (tag,dzp[-1],rms[-1],out15[-1],ht[-1])                 

        self.dzp=array(dzp)
        self.rms=array(rms)
        self.out15=array(out15)
        self.ht=array(ht)
        return output

class bpz_diagnosis:
    def __init__(self,bpz_file='/home/txitxo/bpz/TEST/hdfn.bpz',
                 columns=(1,4,5,6,9,10)):
        # columns correspond to the positions of the variables z_b,odds,z_ml, z_s and m_0
        # in the bpz file
        self.zb,self.tb,self.odds,self.zm,self.zs,self.mo=get_data(bpz_file,columns)
    def stats(self,type='rms',
              odds_min=0.99,
              mo_min=0.,mo_max=99.,
              zs_min=0.,zs_max=6.5,
              t_min=0,t_max=100,
              plots='yes',
              thr=.2):
        good=greater_equal(self.mo,mo_min)
        good*=less_equal(self.mo,mo_max)
        good*=greater_equal(self.zs,zs_min)
        good*=less_equal(self.zs,zs_max)
        good*=greater_equal(self.tb,t_min)
        good*=less_equal(self.tb,t_max)
        self.n_total=len(good)
        self.good=good*greater_equal(self.odds,odds_min)
        self.n_selected=sum(self.good)
        self.d=compress(self.good,(self.zb-self.zs)/(1.+self.zs))
        #try: self.std_thr=std_thr(self.d,thr)
        #except: self.std_thr=1e10
        #try: self.med_thr=med_thr(self.d,thr)
        #except: self.med_thr=1
        #try: self.n_thr=self.n_selected-out_thr(self.d,thr)
        #except: self.n_thr=0
        b=stat_robust(self.d,3,5)
        b.run()
        self.n_remaining=b.n_remaining
        self.n_outliers=b.n_outliers
        self.rms=b.rms
        self.med=b.median
        self.std_log=std_log(self.d)
        if plots=='yes':
            #points(compress(self.good,self.zs),compress(self.good,self.zb),(0.,zs_max,0.,zs_max))
            p=FramedPlot()
            xmin=min(compress(self.good,self.zs))
            xmax=max(compress(self.good,self.zs))
            #print xmin,xmax
            x=arange(xmin,xmax,.01)
            p.add(Curve(x,x,width=3))
            p.add(Curve(x,x+3.*self.rms*(1.+x),width=1))
            p.add(Curve(x,x-3.*self.rms*(1.+x),width=1))
            p.add(Points(compress(self.good,self.zs),compress(self.good,self.zb)))
            p.xlabel=r"$z_{spec}$"
            p.ylabel=r"$z_b$"
            p.show()
            if ask("save plot?"):
                name=raw_input("name?")
                p.write_eps(name)


class bpz_diagnosis_old:
    #This class characterized the quality of a bpz run by comparing
    #the output with the input spectroscopic redshifts
    
    def __init__(self,bpz_file='/home/txitxo/bpz/TEST/hdfn.bpz',
                 columns=(1,5,6,9,10)):
        #columns correspond to the positions of the variables z_b,odds,z_ml, z_s and m_0
        #in the bpz file
        self.zb,self.odds,self.zm,self.zs,self.mo=get_data(bpz_file,columns)
#        print self.zb[:10],self.odds[:10],self.zm[:10],self.zs[:10],self.mo[:10]

    def stats(self,odds_thr=(0.95,0.99),z_lim=(0.,10.)):
        #z_lim selects in the *estimated* quantities
        #parameters controlling the 'purging' of outliers
        d_thr=3. #Should remove only 1% of all points
        n=5
        #Produce stats characterizing the quality of the results
        #rms and fraction of outliers for zm: rms_zm, fo_zm
        if z_lim[0]<>0. or z_lim[1]<>8.:
            good_b=greater_equal(self.zb,z_lim[0])*less_equal(self.zb,z_lim[1])
            good_m=greater_equal(self.zm,z_lim[0])*less_equal(self.zm,z_lim[1])
            good_s=greater_equal(self.zs,z_lim[0])*less_equal(self.zs,z_lim[1])
            zb,zsb,oddsb=multicompress(good_b,(self.zb,self.zs,self.odds))
            zm,zsm=multicompress(good_m,(self.zm,self.zs))
            nzs=sum(good_s)
        else:
            zb=self.zb
            zm=self.zm
            zsb=self.zs
            oddsb=self.odds
            zsm=zsb
            nzs=len(self.zs)
        dzb=(zb-zsb)/(1.+zsb)
        dzm=(zm-zsm)/(1.+zsm)
        nzb=len(dzb)
        nzm=len(dzm)

        greater_equal(self.zm,z_lim[0])*less_equal(self.zm,z_lim[1])

        print "Number of galaxies selected using zs=",nzs
        print "Number of galaxies selected using zb=",nzb
        print "Number of galaxies selected using zm=",nzm
        
        #ZB
#        zb_stat=stat_robust(dzb,d_thr,n)
#        zb_stat.run()
#        mean_zb,rms_zb,n_out_zb,frac_zb=\
#          zb_stat.mean,zb_stat.rms,zb_stat.n_outliers,zb_stat.fraction
#        print "Z_B vs Z_S"
#        print "<z_b-z_s>/(1+zs)=%.4f, rms=%.4f, n_outliers=%i, fraction outliers=%.2f" %\
#              (mean_zb,rms_zb,n_out_zb,frac_zb)

        #ZM
#        zm_stat=stat_robust(dzm,d_thr,n)
#        zm_stat.run()
#        mean_zm,rms_zm,n_out_zm,frac_zm=\
#          zm_stat.mean,zm_stat.rms,zm_stat.n_outliers,zm_stat.fraction
#        print "Z_M vs Z_S"
#        print "<z_m-z_s>/(1+zs)=%.4f, rms=%.4f, n_outliers=%i, fraction outliers=%.2f" %\
#              (mean_zm,rms_zm,n_out_zm,frac_zm)

        #Total Fraction of zm with dz larger than rms_zb
#        f_zm_rms_zb=sum(greater(abs(self.dzm),3*rms_zb))/ float(self.nz)
#        print "Fraction of zm with |zm-zs|/(1+zs) > 3*rms_dzb= %.2f" % f_zm_rms_zb
                
        #Total Fraction of zb with dz larger than rms_zm
#        f_zb_rms_zm=sum(greater(abs(self.dzb),3.*rms_zm))/ float(self.nz)
#        print "Fraction of zb with |zb-zs|/(1+zs) > 3*rms_dzm= %.2f" % f_zb_rms_zm

        #Total Fraction of zb with dz larger than 0.06(1+zs)
        f_zb_0p06=sum(greater(abs(dzb),3*0.06))/ float(nzb)
        print "Fraction of zb with <zb-zs> > 3*0.06(1+z)= %.2f" % f_zb_0p06

        #Total Fraction of zm with dz larger than 0.06(1+zs)
        f_zm_0p06=sum(greater(abs(dzm),3*0.06))/ float(nzm)
        print "Fraction of zm with <zm-zs> > 3*0.06(1+z)= %.2f" % f_zm_0p06

        print "\nSelect objects using odds thresholds\n"
        for i in range(len(odds_thr)):            
            goodo=greater_equal(oddsb,odds_thr[i])
            print "# of objects with odds > %.2f = %.2f " % (odds_thr[i],sum(goodo))
            zbo,zso=multicompress(goodo,(zb,zsb))
            dzbo=(zbo-zso)/(1.+zso)
            zbo_stat=stat_robust(dzbo,d_thr,n)
            zbo_stat.run()
            mean_zbo,rms_zbo,n_out_zbo,frac_zbo=\
               zbo_stat.mean,zbo_stat.rms,zbo_stat.n_outliers,zbo_stat.fraction
            print "     Z_BO vs Z_S"
            print "     <z_bo-z_s>=%.4f, rms=%.4f, n_outliers=%i, fraction outliers=%.2f" %\
                  (mean_zbo,rms_zbo,n_out_zbo,frac_zbo)

            #Total Fraction of zb with dz larger than 0.06(1+zs)
            f_zbo_0p06=sum(greater(abs(dzbo),3*0.06))/ float(len(dzbo))
            print "Fraction of zbo with <zbo-zso> > 3*0.06(1+z)= %.2f" % f_zbo_0p06

            #Plot
            p=FramedPlot()
            p.add(Points(zbo,zso,type='circle'))
            p.xlabel=r"$z_s$"
            p.ylabel=r"$z_b$"
            p.add(Slope(1.,type='dotted'))
            p.show()
            p.write_eps('plot_'+str(odds_thr[i])+'.eps')

            #Otroplot
            p=FramedPlot()
            xz=arange(-1.,1.,0.05)
            hd=hist(dzbo,xz)
            p.add(Histogram(hd,xz[0],0.05))
            p.show()
            
            #Completeness fractions as a function of redshift
            #Odds fractions as a function of magnitude
        
    def plots(self):
        pass

    def webpage(self):
        pass
        #Produce a webpage with a summary of the numeric estimates and plots



def extinction(filtro="u_SDSS",ebv=0.008):
    leff=array(map(float,array((3372.,4404,5428,6509,8090,3683,4393,5519,6602,8046,12660,16732,22152,38079,
                5244,6707,7985,9055,6993,4690,3502,4676,4127,4861,5479,3546,4925,6335,7799,
                9294,3047,4711,5498,6042,7068,8066,4814,6571,8183))))
    a_over_eb=array((5.434,4.315,3.315,2.673,1.940,4.968,4.325,3.240,2.634,1.962,0.902,0.576,
                     0.367,0.153,3.476,2.590,1.991,1.540,2.467,4.035,5.231,4.049,4.552,
                     3.858,3.277,5.155,3.793,2.751,2.086,1.479,5.849,4.015,3.252,2.889,
                     2.435,1.948,3.907,2.649,1.893))
    return match_resol(leff,a_over_eb,filter_center(filtro))*ebv                     

class newfilter:
    # Create and properly store a new filter transmission curve
    # It can take three types of inputs
    # a) A raw filter name
    # b) The raw filter curve
    # c) The filter edges
    # If norm=0, it will not apply QE,atm, etc. corrections (it will then call the file RAW)
    
    def __init__(self,
                 rawfilter,
                 plots=0,
                 Deltawing=20.,
                 wing="BARR",
            	 atmosphere="atm_apo2lapalma.dat",
                 airmass=1.2,
                 aluminum="aluminum2.dat",
                 filter_eff="filter_trans_guess.dat",
                 apply_filter_eff=0,
		 #l=array([3600.,4300.,5500.,6500.,8200.,9500.]),
                 #eff_inst=array([0.45,0.71,0.89,0.89,0.89,0.89]),
                 ccd_type="E2V_MB_DeepDepleted",
                 #mirror=0.85,
                 norm=1,
                 clobber=0,
                 newfilter=None,
                 prefix="JSCH",verbose=1):
        try:
            #If the input is a tuple with two scalars, it will generate the filter shape
            float(rawfilter[0]),float(rawfilter[1])
            internal=1
            apply_filter_eff=1
        except:
            # If not, it will assume that we provide the raw filter shape 
            internal=0.
            if type(rawfilter)==str : 
                xs,yr=get_filter(rawfilter)
            else:
                xs,yr=rawfilter[0],rawfilter[1]
                
        if ccd_type=="LBNL":
            xccd=array([3000.,3250.,3500.,4000.,4500.,5000.,5500.,6000.,6500.,7000.,
                        8000.,9000.,9500.,10000.,10500.])
            ccd=array([0.05,0.20,0.33,0.63,0.71,0.75,0.77,0.80,0.85,0.90,
                       0.91,0.93,0.87,0.58,0.11])
        elif ccd_type=="Marconi":
            xccd=array([3600.,4300.,5500.,6500.,8200.,9500.])
            ccd=array([0.48,0.67,0.90,0.90,0.70,0.30]) #Marconi
        elif ccd_type=="OSIRIS": #http://www.iac.es/project/OSIRIS/
            xccd=array([  3000.,   3500.,   4050.,   4360.,   5000.,   5460.,   6000.,
                          6500.,   7000.,   8000.,   8500.,   9000.,   9500.,   9900.,
                          10500.])
            ccd=array([ 0.     ,  0.239  ,  0.47203,  0.56459,  0.67213,  0.73886,
                        0.80809,  0.86585,  0.89288,  0.88287,  0.79217,  0.60334,
                        0.36898,  0.17646,  0.     ])
        elif ccd_type=="STA_red":
            xccd=array([3000.,3200.,3400.,3600.,3800.,4000.,
                        4500.,5000.,5500.,6000.,6500.,7000.,
                        7500.,8000.,8500.,9000.,9500.,10000.,
                        10500.,11000.])
            ccd=array([0.45,0.43,0.44,0.42,0.61,0.69,
                       0.73,0.72,0.77,0.80,0.83,0.88,
                       0.92,0.93,0.96,0.83,0.62,0.31,
                       0.09,0.05])
        elif ccd_type=="E2V_MB_DeepDepleted":
            #http://www.e2v.com/products-and-services/imaging/space---scientific-imaging/qe-curves/
            xccd=array([2.5e3,3500.,4500.,5000.,5500.,6500.,7500.,8500.,9500.,10500.])
            ccd=array( [0.06, 0.20, 0.79, 0.95, 0.97, 0.94, 0.89, 0.78, 0.45,0.07])

        elif ccd_type=="JPAS_QE":
            #Sent by Toni Marin
            xccd,ccd=get_data("JPAS_QE.dat",(0,1))

        #elif ccd_type=="E2V_pepsi":
        #    ##http://www.e2v.com/products-and-services/imaging/space---scientific-imaging/qe-curves/
        #    xccd=array([2.5e3,3000.,3500.,4000.,4500.,5000.,5500.,6000.,6500.,7000.,7500.,8000.,8500.,9000.,9500.,10000,10500.,11000.])
        #    ccd=array( [0.06, 0.60, 0.76, 0.94, 0.97, 0.92, 0.92, 0.89, 0.86, 0.84, 0.80, 0.74, 0.62, 0.46, 0.30,  0.15, 0.07,     0.])

        # Get "raw" filter curves
        if internal:
            # Generate raw transmission curve
            x1=float(rawfilter[0])
            x2=float(rawfilter[1])
            if wing=="BARR":
                xrmin=x1-Deltawing*1./0.8 #Wing width def. betw. 0.9-0.1 trans.
                xrmax=x2+Deltawing*1./0.8 #Wing width def. betw. 0.9-0.1 trans.
                #Make sure grid contains x1 and x2
                xs=fixed_grid(xrmin-2.,xrmax+2.,1.,(x1,x2))
                y1=0.8/Deltawing*(xs-x1)+1.
                y2=1.-0.8/Deltawing*(xs-x2)                
            elif wing=="Gauss":
                xrmin=x1-Deltawing*5.
                xrmax=x2+Deltawing*5.
                #Make sure grid contains x1 and x2
                xs=fixed_grid(xrmin,xrmax,1.,(x1,x2))
                y1=exp(-(xs-x1)**2/(2.*Deltawing**2))
                y2=exp(-(xs-x2)**2/(2.*Deltawing**2))                 
            yr=xs*0.+1.
            yr=where(less(xs,x1),y1,yr)
            yr=where(greater(xs,x2),y2,yr)
            yr=clip(yr,0.,1.)
            # Search for duplicated data
            xs,yr=multisort(xs,(xs,yr))
            diff=xs-xs[-1:len(xs)-1]
            if not N.all(diff):
                xs,yr=multicompress(diff,(xs,yr))
                print "%i Duplicates found in the filter wavelenghts" % logical_not(diff).sum()
            ask()

        else:
            if max(yr)>1.:
                if verbose:
                    print ""
                    print "Input filter does not appear to be normalized!"
                    print ""
                    print "Setting maximum transmission equal to 1..."
                    print "Will apply effective filter correction"
                yr/=max(yr)            
                apply_filter_eff=1

        if norm:
            la,at=U.get_data(atmosphere,(0,1))
            at=at**airmass
            lal,al=U.get_data(aluminum,(0,1))
            lf,fe=U.get_data(filter_eff,(0,1))
            # Calculate effective filter transmission
            yc=clip(U.match_resol(xccd,ccd,xs),0.,1.)
            ya=clip(U.match_resol(la,at,xs),0.,1.)
            yal=clip(U.match_resol(lal,al,xs),0.,1.)
            if apply_filter_eff: #This correction tries to take into accout typical transp.
                yfe=clip(U.match_resol(lf,fe,xs),0.,1.)
            else:
                yfe=yal*0.+1.
            if plots:
                plot(xs,ya,xs,yal,xs,yfe,xs,yc)
                leyenda=["atm","alum","gen_filter_eff","ccd"]
            y=yr*ya*yal*yfe*yc
            y=clip(y,0.,1.)
            self.ya=ya
            self.yc=yc
            self.yfe=yfe
            self.yal=yal
        else:
            y=yr


        if newfilter==None:
            if not type(rawfilter)==str:
                #lmean=N.trapz(xs*yr,xs)/trapz(yr,xs)
                #lwidth=N.trapz(yr,xs)
                lmean=N.trapz(xs*yr,xs)/trapz(yr,xs)
                lwidth=(match_resol(yr[::-1]*xs[::-1]/max(yr)/mean(xs),xs[::-1],0.5)-
                        match_resol(yr*xs/max(yr)/mean(xs),xs,0.5))
                newfilter=prefix+str(int(round(lmean)))+"_"+str(int(round(lwidth)))
                if not norm: newfilter+="_raw"
            else:
                if norm: newfilter=rawfilter+"_thru"
                else: newfilter=rawfilter+"_raw"

        if plots:
            plot(xs,yr,xs,y)
            try:
                leyenda=leyenda+["filter","final"]
            except:
                leyenda=["filter","final"]
            legend(leyenda)
            show()

        if os.path.exists(fil_dir+newfilter+".res") and not clobber:
            if verbose:
                print "Filter %s already exists " % fil_dir+newfilter+".res"
                print "If you want to rewrite the filter set clobber=1"
                print "Exiting without doing anything"
        else:
            put_data(fil_dir+newfilter+".res",(xs,y))
            if verbose:
                print "Filter transmission"+fil_dir+newfilter+".res written"
            for fileAB in glob(ab_dir+"*."+newfilter+".AB"):
                if verbose:
                    print "Deleting "+fileAB
                os.remove(fileAB)

            
        self.xs=xs
        self.yr=yr
        self.y=y
        print newfilter,alltrue(xs[1:]-xs[:-1])
        self.newfilter=newfilter


def genfilterfwhm(wlc,fwhm,
                  filter0="schott_filter_A.res",
                  ccd_type="JPAS_QE",norm=0,prefix="JT250",plots=0,clobber=1,verbose=1):
    dx=1.
    xr,yr=get_filter(filter0)
    xc=filter_center(filter0)
    yr/=max(yr) # Normalize filter to max=1
    fwhm0=filter_fwhm(filter0)
    padding=(fwhm-fwhm0)
    print fwhm0,padding
    #Set filter center at 0
    xr-=xc
    xr=where(xr>=0,xr+padding/2.+wlc,xr-padding/2.+wlc)
    x=arange(min(xr),max(xr)+dx,dx)
    y=match_resol(xr,yr,x)
    x,y=multicompress(greater(y,.001),(x,y))
    x =concatenate(
        (x[0]-array((2.,1.))*dx,
         x,
            x[-1]+array((1.,2.))))
    y= concatenate(
        (array((0.,0.)),y,array((0.,0.)))
        )
    a=newfilter((x,y),clobber=clobber,ccd_type=ccd_type,
                apply_filter_eff=1.,prefix=prefix,norm=norm,plots=plots,verbose=verbose)
    return a.newfilter

def genfilterfwhm_old(wlc,fwhm,
                      filter0="schott_filter_A.res",
                      ccd_type="JPAS_QE",norm=1,prefix="JT250",plots=0,clobber=1,verbose=1):
    xr,yr=get_filter(filter0)
    dx=median(xr[1:]-xr[:-1])/2.
    fwhm0=filter_fwhm(filter0)
    im=argmax(yr)
    yr/=yr[im] # Normalize filter to max=1
    xr+=-xr[im] # Center filter at 0
    # Split filter into two halves, create wings (r<99%)
    l=less(xr,xr[im])*less_equal(yr,0.99)
    r=greater(xr,xr[im])*less_equal(yr,0.99)
    yrl=yr[l]
    xrl=xr[l]
    yrr=yr[r]
    xrr=xr[r]
    deltaw=fwhm0-(xrr[0]-xrl[-1])
    xrr+=-xrr[0]
    xrr+=dx
    xrl+=-xrl[-1]
    yrr/=max(yrr)
    yrl/=max(yrl)
    xrl,yrl=multicompress(greater_equal(yrl,0.002),(xrl,yrl))
    xrr,yrr=multicompress(greater_equal(yrr,0.002),(xrr,yrr))
    xrl,yrl=pad_data(xrl,yrl,2,0.,dx,r=0)
    xrr,yrr=pad_data(xrr,yrr,2,0.,dx,l=0)
    xw=concatenate((xrl,xrr))
    yw=concatenate((yrl,yrr))
    xc=filter_center((xw+1000.,yw))-1000.
    fc=filter_fwhm((xw+1000.,yw))
    #Constant central part of length delta
    delta=fwhm-fc-2.
    x=concatenate((xrl-delta/2-xc,      arange(-delta/2.-xc,delta/2.-xc,dx),xrr+delta/2.-xc))+wlc+1.
    y=concatenate((yrl        ,1.+0.*arange(-delta/2.,delta/2.,dx),yrr         ))
    g=concatenate(([1],diff(x)))
    x,y=multicompress(g,(x,y))
    a=newfilter((x,y),clobber=clobber,ccd_type=ccd_type,
                          apply_filter_eff=1.,prefix=prefix,norm=norm,plots=plots,verbose=verbose)
    return a.newfilter

def PSF_corr(seeing,area_phot):
    """For a gaussian PSF with FWHM= seeing it calculates the difference
    between the total magnitude and that enclosed by a aperture of area_phot"""
    sigma=seeing/2.35
    r=arange(0.,seeing*6.,sigma/20.)
    r_phot=sqrt(area_phot/pi)
    int_profile=add.accumulate(exp(-r*r/2./sigma**2)*r)
    return flux2mag(match_resol(r,int_profile,r_phot)/int_profile[-1])

#class optimum_phot:
#    """ 
#    Calculates S/N for a given profile assuming that optimal
#    photometric weighting is performed, i.e. if f(r) is the
#    radial flux of an object and Fo -+is its total flux, it calculates
#    
#    S/N =F_w / sigma_pixel
#    
#    where F_w = int w(r) * f(r) * 2 pi r dr
#    
#    and w(r)= f(r)/Fo
#    """
#    def __init__(self,
#                 seeing=0.8,
#                 model="gauss",
#                 re=0.37,
#                 pixel=0.4):
#        npixel=1./pixel**2
#        noise_pixel=sqrt(npixel)
#        sigma=seeing/2.354
#        r=arange(0.,5.*seeing,0.001)
#        if model=="gauss":
#            f=1./2./pi/sigma**2*exp(-r*r/2./sigma**2)
#            norm=add.reduce(f*2.*pi*r)
#            w=f/norm
#        else:
#            f=exp(-7.67*(r/re)**0.25)
#            norm=add.reduce(f*2.*pi*r)
#            w=f/norm
#        F_w=add.reduce(w*f*2.*pi*r)/add.reduce(w*2.*pi*r)
#        noise=add.reduce(w*2.*pi*r*
                         
        
class signal:
    def __init__(self,
                 mag=24.,
                 cal="AB",
                 filtername="V_Johnson",
                 t_exp=1000.,
                 pixel_size=0.4534,
                 n_exp=4.,
                 readout=6.,
                 area_object=7.1,# arcsec2
                 area_det=7.1,
                 sky="sky_lapalma",
                 fudge="INT",
                 mirror=3.89,
                 realistic=0.):

        # If realistic, increase background noise by 40% to reflect that the real noise is about 20% higher
        # due to the presence of unresolved objects for typical apertures
        self.realistic=realistic
        
        area_sky=1. #All sky spectra we use are defined on 1 sq. arcsec^2

        # Calculate photon numbers quantities per area
        if cal=="Vega": mag=VegatoAB(mag,filtername)

        if fudge=="INT":
            #This fudge factor makes our results normalized to those of SIGNAL at the INT 
            #It assumes a loss of transmission which is not accounted for
            xf=array((3353.,4347.2,5491.7,6414.1,8027.40))
            #ff=array((1.,0.751,0.745,0.89,1.))
            ff=array((0.75,0.75,0.75,0.9,.9))
            fudge=clip(match_resol(xf,ff,filter_center(filtername)),0.,1.)
        else:
            fudge=1.

        c_sky=sed2photons(sky,filtername,mirror=mirror,area_det=area_det,area_object=area_sky,t_exp=t_exp)*fudge
        c_obj=mag2photons(mag,filtername,mirror=mirror,area_det=area_det,area_object=area_object,t_exp=t_exp)*fudge
        n_pixel=area_det/pixel_size**2
        
        self.n_pixel=n_pixel
        self.readout=readout
        self.n_exp=n_exp
        self.mirror=mirror
        self.area_det=area_det
        self.area_sky=area_sky
        self.area_object=area_object
        self.t_exp=t_exp
        self.fudge=fudge
        self.filtername=filtername
        self.c_obj=c_obj
        self.c_sky=c_sky
        self.readnoise=sqrt(readout**2*n_exp*n_pixel)
        if self.realistic:
            self.ratio_sky_sn=sqrt(c_sky*1.4)/self.readnoise
            self.sn=c_obj/sqrt(c_sky*1.4+c_obj+n_exp*n_pixel*readout**2)
        else:
            self.ratio_sky_sn=sqrt(c_sky)/self.readnoise
            self.sn=c_obj/sqrt(c_sky+c_obj+n_exp*n_pixel*readout**2)
        self.mag=mag
        self.t_exp=t_exp
        self.sky=sky

    def sn_mag(self,mag,t):
        cobj=mag2photons(mag,self.filtername,
                         mirror=self.mirror,
                         area_det=self.area_det,
                         area_object=self.area_object,
                         t_exp=t)*self.fudge
        csky=sed2photons(self.sky,self.filtername,
                         mirror=self.mirror,
                         area_det=self.area_det,
                         area_object=self.area_sky,
                         t_exp=t)*self.fudge
        if self.realistic:
            sn_mag=cobj/sqrt(csky*1.4+cobj+self.n_exp*self.n_pixel*self.readout**2)
            self.ratio_sky_sn=sqrt(csky*1.4+cobj)/self.readnoise
        else:
            sn_mag=cobj/sqrt(csky+cobj+self.n_exp*self.n_pixel*self.readout**2)
            self.ratio_sky_sn=sqrt(csky+cobj)/self.readnoise            
        return sn_mag

    def mag_sn(self,sn=5.):
        def f(x): 
            return abs(self.sn_mag(mag=x,t=self.t_exp)-sn)
        m=optimize.golden(f,brack=(0.,40.),tol=0.001)        
        return m

    def t_sn(self,sn=5.):
        def f(x): return abs(self.sn_mag(mag=self.mag,t=x)-sn)
        t=optimize.golden(f,tol=0.01)        
        return t

    def mag5(self):
        return self.mag_sn(5.)

    def t5(self):
        return self.t_sn(5.)


def filtersfromcols(cols="hdfn_z.columns"):
    filters=[]
    lineas=open(cols,"r").readlines()
    for l in lineas:
        if l[0]=="#": continue
        if "," in l: 
            filtro=l.split()[0]
            if "." in filtro: 
                filtro=filtro.split(".")[0]
            filters.append(filtro)

    return filters


class fluxcomparison:
    def __init__(self,
                 bpz="hdfn_z.bpz",
                 comp="hdfn_z.flux_comparison",
                 cols="hdfn_z.columns",
                 lib=None,
                 exclude=None,
                 new=0,
                 sed_dir=sed_dir
                 ):
        #exclude is a tuple containing filter names which have been excluded when running BPZ
        # If the library is included, it will calculate the effective wavels. (takes a while)

        # Read the filters
        filtros=filtersfromcols(cols)
        nfil0=len(filtros)
        dzp=get_data(cols,3,nrows=nfil0)
        dmin=float(re.split('MIN_DM=',
                            open(bpz,"r").read())[1][:5].split("#")[0])
        
        d0=sqrt(e_mag2frac(dzp)**2+e_mag2frac(dmin)**2)
        
        if exclude<>None:
            for f in exclude:
                print f
                filtros.remove(f)

        nf=len(filtros)
            
        # Read the redshifts and the types
        t,zs,m,o=get_data(bpz,(4,9,10,5))
        ng=len(zs)

        t=ravel(outer(t,ones(nf)))
        it=map(lambda x: round(x-1.),t)
        m=ravel(outer(m,ones(nf)))
        o=ravel(outer(o,ones(nf)))
        zs=ravel(outer(zs,ones(nf)))
        i_filtros=ravel(outer(ones(ng),arange(nf)))

        pivot=array(map(lambda x: pivotal_wl(x),filtros))
        lp=ravel(outer(ones(ng),pivot))

        centers=array(map(lambda x: filter_center(x), filtros))
        c=ravel(outer(ones(ng),centers))
        lc=c/(1.+zs)
        d0=ravel(outer(ones(ng),d0))


        # Now determine the effective wavelengths of the filters, taking 
        # into account the SED we are integrating
        if lib<>None:
            efflam_file=comp.split(".")[0]+".efflam"
            if niet(efflam_file) or new:
                print "Calculating effective wawelenghts for %s" % comp
                sed=get_lib(lib)
                nsed=len(sed)
                l0=ravel(outer(ones(ng)*1.,ones(nf)))
                xz=arange(0.,10.,0.01)
                for i in range(nsed):
                    print sed[i],
                    for j in range(nf):
                        g=equal(it,i)*equal(j,i_filtros)
                        l0[g]=efflam_z_sed(sed[i],filtros[j],zs[g],sed_dir=sed_dir)
                    print g.sum()
                put_data(efflam_file,(l0,))
                print "Saving effective wawelenghts to %s" % efflam_file
            else:
                print "Reading effective wawelenghts from %s" % efflam_file
                l0=get_data(efflam_file,0)
            l=l0/(1.+zs)
        else:
            l0=c
            l=lc
            

        # Read the fluxes
        a=get_data(comp,4)
        a=ravel(outer(a,centers*0.+1.))
        ft=ravel(get_2Darray(comp,range(5,5+nf)))
        fo=ravel(get_2Darray(comp,range(5+nf,5+2*nf)))
        dfo=ravel(get_2Darray(comp,range(5+2*nf,5+3*nf)))
        fo=fnu2fl(lp,fo)
        dfo=fnu2fl(lp,dfo)
        ft=fnu2fl(lp,ft)
        self.l0=l0
        self.l=l
        self.centers=centers
        self.lc=lc
        self.fo=fo
        self.ft=ft
        self.dfo=dfo
        self.t=t
        self.m=m
        self.o=o
        self.zs=zs
        self.c=c
        self.a=a
        self.i_filtros=i_filtros
        self.filtro=map(lambda x: filtros[int(x)],i_filtros)
        self.filtros=filtros
        self.d0=d0

class fluxbasis:
    # Given a set of basis functions (e.g. Chebyshev polynomials)
    # defined as spectra and kept in the SPECTRA folder
    # and a set of redshifts and filters, this class calculates
    # the AB fluxes (in the BPZ sense) through those filters
    # (they can be later used for fits, etc.)
    def __init__(self,zs,filtros,basis_list):
        basis=get_str(basis_list,0)
        Result=zeros(len(zs),'float',len(basis))
        for i in range(len(basis)):
            for j in range(len(filtros)):  
                g=filtros==filtros[j]
                print i,j,basis[i],filtros[j]
                Result[g,i]=f_z_sed_AB(basis[i],filtros[j],zs[g])

        self.Result=Result



class testsignal:
    def __init__(self,ccd="Marconi",sky="sky_lapalma",airmass=1.,c=0):

        a=newfilter("U_RGO",newfilter="U_RGO_%s" % ccd,apply_filter_eff=0,clobber=c,
                    airmass=airmass,plots=1,ccd_type=ccd)
        a=newfilter("B_Harris",newfilter="B_Harris_%s" % ccd,apply_filter_eff=0,clobber=c,
                    airmass=airmass,plots=1,ccd_type=ccd)
        a=newfilter("V_Harris",newfilter="V_Harris_%s" % ccd,apply_filter_eff=0,clobber=c,
                    airmass=airmass,plots=1,ccd_type=ccd)
        a=newfilter("R_Harris",newfilter="R_Harris_%s" % ccd,apply_filter_eff=0,clobber=c,
                    airmass=airmass,plots=1,ccd_type=ccd)
        a=newfilter("I_Cousins",newfilter="I_Cousins_%s" % ccd,apply_filter_eff=0,clobber=c,
                    airmass=airmass,plots=1,ccd_type=ccd)

        print "filtername, SN, SN_sig, C_sky, C_sky_sig,C_obj,C_obj_sig"
        a=signal(24.,"Vega","U_RGO_Marconi",250.,
                 pixel_size=0.33,n_exp=1.,readout=4.,
                 area_det=0.33**2,sky=sky,#1.*500./542.65,
                 mirror=4.41)

        print a.filtername,a.sn,1.07,a.c_sky,118.34,a.c_obj,13.

        a=signal(24.,"Vega","B_Harris_Marconi",250.,
                 pixel_size=0.33,n_exp=1.,readout=4.,area_det=0.33**2,
                 sky=sky,#0.9*720./695.,
                 mirror=4.41)
        print a.filtername,a.sn,4.02,a.c_sky,460.,a.c_obj,96.

        a=signal(24.,"Vega","V_Harris_Marconi",250.,
                    pixel_size=0.33,n_exp=1.,readout=4.,area_det=0.33**2,
                 sky=sky,#0.7*860./826.,
                 mirror=4.41)

        print a.filtername,a.sn,3.14,a.c_sky,1102.,a.c_obj,110.

        a=signal(24.,"Vega","R_Harris_Marconi",250.,
                    pixel_size=0.33,n_exp=1.,readout=4.,area_det=0.33**2,
                 sky=sky,#0.7*1330./1127.,
                 mirror=4.41)
        print a.filtername,a.sn,2.32,a.c_sky,2954.,a.c_obj,129.

        a=signal(24.,"Vega","I_Cousins_Marconi",250.,
                    pixel_size=0.33,n_exp=1.,readout=4.,area_det=0.33**2,
                 sky=sky,#0.9*1400./1487.,
                 mirror=4.41)
        print a.filtername,a.sn,1.49,a.c_sky,5271.,a.c_obj,92.

        #NB tests

        a=signal(24.,"Vega","F4331_118",250.,
                 pixel_size=0.33,n_exp=1.,readout=4.,area_det=0.33**2,
                 sky=sky,#0.9*1400./1487.,
                 mirror=4.41)
        print a.filtername,a.sn,1.52,a.c_sky,75.46,a.c_obj,15.77


        a=signal(24.,"Vega","F5446_118",250.,
                    pixel_size=0.33,n_exp=1.,readout=4.,area_det=0.33**2,
                 sky=sky,#0.9*1400./1487.,
                 mirror=4.41)
        print a.filtername,a.sn,1.12,a.c_sky,151.24,a.c_obj,15.12

        a=signal(24.,"Vega","F6470_118",250.,
                    pixel_size=0.33,n_exp=1.,readout=4.,area_det=0.33**2,
                 sky=sky,#0.9*1400./1487.,
                 mirror=4.41)
        print a.filtername,a.sn,0.67,a.c_sky,262.16,a.c_obj,11.44

        a=signal(24.,"Vega","F8143_118",250.,
                    pixel_size=0.33,n_exp=1.,readout=4.,area_det=0.33**2,
                 sky=sky,#0.9*1400./1487.,
                 mirror=4.41)
        print a.filtername,a.sn,0.36,a.c_sky,444.28,a.c_obj,7.72


def test():
    """ Tests some functions defined in this module"""

    test='flux'
    Testing(test)

    x=arange(912.,10001.,.1)
    r=exp(-(x-3500.)**2/2./200.**2)
    f=1.+sin(x/100.)
    
    e_ccd=add.reduce(f*r*x)/add.reduce(r*x)
    e_noccd=add.reduce(f*r)/add.reduce(r)

    r_ccd=flux(x,f,r,ccd='yes',units='lambda')
    r_noccd=flux(x,f,r,ccd='no',units='lambda')

    if abs(1.-e_ccd/r_ccd)>1e-6 or abs(1.-e_noccd/r_noccd)>1e-6: raise test

    #print '        f_lambda          '
    #print 'Results                  Expected'
    #print 'CCD ',r_ccd,e_ccd
    #print 'No CCD ',r_noccd,e_noccd

    nu=arange(1./x[-1],1./x[0],1./x[0]/1e2)*clight_AHz
    fn=(1.+sin(clight_AHz/100./nu))*clight_AHz/nu/nu
    xn=clight_AHz/nu
    rn=match_resol(x,r,xn)
    e_ccd=add.reduce(fn*rn/nu)/add.reduce(rn/nu)
    e_noccd=add.reduce(fn*rn)/add.reduce(rn)
    r_ccd=flux(x,f,r,ccd='yes',units='nu')
    r_noccd=flux(x,f,r,ccd='no',units='nu')

    #print '           f_nu           '
    #print 'Results                  Expected'
    #print 'CCD',r_ccd,e_ccd
    #print 'no CCD',r_noccd,e_noccd

    if abs(1.-e_ccd/r_ccd)>1e-6 or abs(1.-e_noccd/r_noccd)>1e-6: raise test

    test='AB'
    Testing(test)
    if AB(10.**(-.4*48.60))<>0.: raise test
    
    test='flux2mag and mag2flux'
    Testing(test)
    m,f=20.,1e-8
    if mag2flux(m)<>f: raise test
    if flux2mag(f)<>m: raise test

    test='e_frac2mag and e_mag2frac'
    Testing(test)
    f=1e8
    df=1e7/f
    m=flux2mag(f)
    dm=m-flux2mag(f*(1.+df))
    if abs(e_frac2mag(df)-dm)>1e-12: 
	print abs(e_frac2mag(df)-dm)
	raise test
    if abs(e_mag2frac(dm)-df)>1e-12: 
	print e_mag2frac(dm),df
	raise test

    test='etau_madau'
    #Un posible test es generar un plot de la absorpcion a distintos redshifts
    #igual que el que viene en el paper de Madau.

    test='f_z_sed'
    Testing(test)
    #Estimate fluxes at different redshift for a galaxy with a f_nu\propto \nu spectrum
    # (No K correction) and check that their colors are constant
    x=arange(1.,10001.,10.)
    f=1./x
    put_data(sed_dir+'test.sed',(x,f))
    z=arange(0.,10.,.25)
    b=f_z_sed('test','B_Johnson.res',z,ccd='no',units='nu',madau='no')
    v=f_z_sed('test','V_Johnson.res',z,ccd='no',units='nu',madau='no')
    c=array(map(flux2mag,b/v))
    if(sometrue(greater(abs(c-c[0]),1e-4))): 
	print c-c[0]
	raise test

    test='VegatoAB' # To be done
    test='ABtoVega'
    test='likelihood'
    
    #Test: generar un catalogo de galaxias con colores, e intentar recuperar 
    #sus redshifts de nuevo utilizando solo la likelihood

    test='p_and_minchi2' # To be done
    test='prior'

    test='interval'
    test='odds'

    test=' the accuracy of our Johnson-Cousins-Landolt Vega-based zero-points'
    Testing(test)
    #filters=['U_Johnson.res','B_Johnson.res','V_Johnson.res','R_Cousins.res',
    #'I_Cousins.res']

    filters=[
        'HST_ACS_WFC_F435W',
        'HST_ACS_WFC_F475W',
        'HST_ACS_WFC_F555W',
        'HST_ACS_WFC_F606W',
        'HST_ACS_WFC_F625W',
        'HST_ACS_WFC_F775W',
        'HST_ACS_WFC_F814W',
        'HST_ACS_WFC_F850LP'
        ]
    

    ab_synphot=array([
        -0.10719,
        -0.10038,
        8.743e-4,
        0.095004,
        0.174949,
        0.40119,
        0.44478,
        0.568605
        ])
       

    f_l_vega=array([
        6.462e-9,
        5.297e-9,
        3.780e-9,
        2.850e-9,
        2.330e-9,
        1.270e-9,
        1.111e-9,
        7.78e-10])


    print '     f_l for Vega'
    sufix='cgs A^-1'
    print '                               f_lambda(Vega)     synphot(IRAF)   difference %'
    for i in range(len(filters)):
        f_vega=f_z_sed(Vega,filters[i],ccd='yes')
	tupla=(ljust(filters[i],16),f_vega,f_l_vega[i],f_vega/f_l_vega[i]*100.-100.)
	print '     %s         %.6e       %.6e      %.4f'%tupla +"%"


    print '    '
    print '    AB zeropoints for Vega '
    sufix='cgs Hz'
    tipo='nu'
    print "                                AB zero point     synphot(IRAF)   difference"
    for i in range(len(filters)):
	f_vega=f_z_sed(Vega,filters[i],units=tipo,ccd='yes')
	tupla=(ljust(filters[i],16),AB(f_vega),ab_synphot[i],AB(f_vega)-ab_synphot[i])
	print '     %s         %.6f       %.6f      %.6f' % tupla
   

    print '    '
    print '    AB zeropoints for a c/lambda^2 spectrum (flat in nu)'
    sufix='cgs Hz'
    tipo='nu'
    print "                                 Result             Expected  "
    for i in range(len(filters)):
	f_flat=f_z_sed('flat',filters[i],units=tipo,ccd='yes')
	tupla=(ljust(filters[i],16),AB(f_flat),0.)
	print '     %s         %.6e       %.6f' % tupla
   

    print ''
    print '         Everything OK    in   bpz_tools '
    print ''


if __name__ == '__main__':
    test()
else:
    pass
#    print 'bpz_tools loaded as module'
