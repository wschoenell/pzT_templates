#This module contains several functions which calculate 
#observational quantities affected by cosmology

from bpz_tools import *
from useful import get_data
import glob
import numpy as N
try:
    from scipy import integrate
    from scipy import unique
    from scipy import interpolate
except:
    print "Could not import scipy properly"

cho=2.99e3       #c/H_0 in Mpc
ht=9.7776e9      #hubble time in h^-1 yr
# Planck http://arxiv.org/abs/1303.5076
Omega_matter=0.315
Omega_lambda=1.-Omega_matter
H_over_H0=0.673
cosmo=(Omega_matter,Omega_lambda,H_over_H0)

# Comment out if memory limited
global dict_reobs
dict_reobs={}

# K-corrections and the like
def kcor(z,sed,filtro):
    """K-correction in a giver filter for the spectrum SED at redshift z
    ( m=M+5log(D_L/10pc)+K(z)   )"""
    fo=f_z_sed_AB(sed,filtro)
    try:
        len(z)
    except:
        z=N.array([z])    
    k=2.5*N.log10((1.+z)*fo/f_z_sed_AB(sed,filtro,z))
    if len(k)==1: 
        return k[0]
    else: 
        return k

def ab_name(seds,filters): 
    if filters[-4:]==".res": filters=filters[:-4]
    if seds[-4:]==".sed": seds=seds[:-4]
    return seds+'.'+filters+'.AB'

def genmzt(cat="fireworks_spec.cat",cols="fireworks_spec.columns",output=None,
           seds="f0419102.list",interp=2,onlytype="yes",mzt=None,tail=None):
    if not output:  output=cat[:-4]+"_"+seds[:-5]+".bpz"
    command="python /Users/txitxo/bpz/bpz.py %s -COLUMNS %s -SPECTRA %s -INTERP %i -ONLY_TYPE %s -VERBOSE no" % (cat,cols,seds,interp,onlytype)
    if tail: command+=" %s" % tail    
    ejecuta(command)
    if "Z_S" in open(output,"r").read():
        m,z,t=get_data(output,(9+1,1,4))    
    else:
        m,z,t=get_data(output,(9,1,4))            
    if not mzt: mzt=output[:-4]+".mzt"
    put_data(mzt,(m,z,t))
    print "Data saved in ",mzt


class genmzt_fromprior:
    def __init__(self,prior="prior_B13",out="test.mzt",
                 area=1.73,
                 mmin=18.,
                 mmax=23.5,
                 zmin=0.001,
                 zmax=5.,
                 lib="eB11.list",nt=11,
                 filter0="HST_ACS_WFC_F814W",
                 filter_output=None):
                 #Everything is assumed to be AB
        nx=1000.
        P=__import__(prior)
        mmax0=mmax
        xz=linspace(zmin,zmax,nx)
        nz=len(xz)
        if filter_output<>None and filter0<>filter_output:
            seds=get_lib(lib)
            #Determine the maximum color (exclude last template, too noisy)
            c=zeros((nz,nt-1),"float")
            M=c*0.
            for i in range(len(seds)-1):
                M[:,i]=m_abs(mmax,xz,seds[i],filter_m=filter_output,filter_M="B_Johnson",cal_M="Vega")
                c[:,i]=color_z(seds[i],filter_output,filter0,xz,calibration="AB")
                #plot(M[:,i],c[:,i])
            #show()
            Dc=c[M>-26].min()
            print "Dc=",Dc
            mmax=mmax0-Dc #Make the sample deeper in the I band to avoid losing bluer objects
        else:
            mmax=mmax0
        Dm=mmax-mmin
        dm=Dm/nx
        #Determine total number of galaxies based on number counts
        xm=linspace(mmin,mmax,nx)
        #Generate the prior and calculate its maximum
        pztm=P.function(xz,xm,nt=nt)
        M=[]
        Z=[]
        T=[]
        for im in range(len(xm)):
            pmax=pztm[:,:,im].max()
            expand=pztm[:,:,im].size/(pztm[:,:,im].sum()/pmax)
            n=int(nm(xm[im],dm)*area*expand)
            #print "m,n,nm=",xm[im],n,nm(xm[im],dm)*area
            mr=N.random.random(size=n)*dm+(xm[im]-dm/2.)
            itr=N.random.randint(nt,size=n)
            izr=N.random.randint(nz,size=n)
            l=less(N.random.random(size=n)*pmax,pztm[izr,itr,im].ravel())
            m,it,iz=multicompress(l,(mr,itr,izr))
            M.append(m)
            Z.append(xz[iz])
            T.append(it+1)
        m=concatenate(M)
        t=concatenate(T)
        z=concatenate(Z)
        t=clip(t+N.random.normal(0.,0.5,len(t)),1.,nt)
        if filter_output<>None and filter0<>filter_output:
            a=getcolors(lib,filter_output,filter0,(m,z,t),99.)
            #Use 99 to avoid excluding any objects, we cut later
            g=less_equal(a.mags,mmax0)*greater_equal(a.mags,mmin)
            m,z,t=multicompress(g,(a.mags,z,t))
        self.m,self.z,self.t=m,z,t
        put_data(out,(self.m,self.z,self.t),format="%5.3f %5.3f %4.2f")
                    

class getcolors:
    def __init__(self,
                 library="pegase.list",
                 filters="f05a_225.filters", # A .columns file will do too
                 filter0="HST_ACS_WFC_F775W.res",
                 mzt="FIREWORKS_mock_f05a225.mzt",
                 momax=None,
                 dz_ab=dz_ab,
                 colorfile=None,
                 new=1,
                 newAB=0,verbose=0):
        if type(mzt)<>type(""): 
            m,z,t=mzt
            mzt="mzt.mzt"
        else:
            m,z,t=get_data(mzt,(0,1,2))
        
        if colorfile==None:
            colorfile=mzt[:-4]+"_"+filters.split(".")[0]+".col"
        if os.path.exists(colorfile) and not new:
            if verbose: print colorfile," exists"
            else: pass
        else:
            if verbose: print "Generating ",colorfile
            seds=get_str(sed_dir+library,0)
            try:
                filtros=get_filter_list(filters)
            except:
                if type(filters)==tuple or type(filters)==type([1,2]):                    
                    filtros=filters
                else:
                    filtros=[filters]
            if momax<>None:
                good=less_equal(m,momax)
                m,z,t=multicompress(good,(m,z,t))
            ng=len(m)
            if verbose: print "ng=",ng
            t+=-1
            nf=len(filtros)
            mags=zeros((ng,nf),'float')
            tipos=unique(t)

            d={}

            if verbose: print "types=",tipos

            suma=len(t)
            checksum=0.
            for i in tipos:
                g=equal(t,i)
                checksum+=g.sum()
                if g.sum()==0: continue
                if floor(i)==i: # For integer types
                    j=int(i)
                    #print i,"j=",j
                    # Load m0
                    if verbose: print j,seds[j],g.sum()
                    ab_file=ab_name(seds[j],filter0)
                    try: d[ab_file]
                    except: d[ab_file]=get_AB(ab_file,newAB=newAB)
                    fold=match_resol(d[ab_file][0],d[ab_file][1],z[g])
                    # Load the rest of the filters
                    for k in range(nf):
                        ab_file=ab_name(seds[j],filtros[k])
                        try: d[ab_file]
                        except: d[ab_file]=get_AB(ab_file,newAB=newAB)
                        fnew=clip(match_resol(d[ab_file][0],d[ab_file][1],z[g]),tiny,1./tiny)
                        mags[g,k]=clip(2.5*(log10(fold/fnew))+m[g],-99.,99.)                            
                else:
                    l=int(floor(i))
                    r=int(ceil(i))
                    wr=i-l #Checked with bpzsim
                    wl=r-i
                    #print i,"l,r",l,r,wr,wl
                    if verbose: print i,seds[l],seds[r],g.sum()
                    ab_file_l=ab_name(seds[l],filter0)
                    try: d[ab_file_l]
                    except: d[ab_file_l]=get_AB(ab_file_l)
                    fold_l=match_resol(d[ab_file_l][0],d[ab_file_l][1],z[g])

                    ab_file_r=ab_name(seds[r],filter0)
                    try: d[ab_file_r]
                    except: d[ab_file_r]=get_AB(ab_file_r)
                    fold_r=match_resol(d[ab_file_r][0],d[ab_file_r][1],z[g])

                    fold=wl*fold_l+wr*fold_r
                    # Load the rest of the filters
                    for k in range(nf):
                        ab_file_l=ab_name(seds[l],filtros[k])
                        try: d[ab_file_l]
                        except: d[ab_file_l]=get_AB(ab_file_l,newAB=newAB)
                        fnew_l=clip(match_resol(d[ab_file_l][0],d[ab_file_l][1],z[g]),tiny,1./tiny)

                        ab_file_r=ab_name(seds[r],filtros[k])
                        try: d[ab_file_r]
                        except: d[ab_file_r]=get_AB(ab_file_r,newAB=newAB)
                        try:
                            fnew_r=match_resol(d[ab_file_r][0],d[ab_file_r][1],z[g])
                        except:
                            print ab_file_r
                            #plot(d[ab_file_r][0],d[ab_file_r][1])
                            #show()

                        fnew=wl*fnew_l+wr*fnew_r
                        mags[g,k]=clip(2.5*log10(fold/(fnew+tiny))+m[g],-99.,99.)
                

            put_2Darray(colorfile,mags)

        try:
            self.mags=get_2Darray(colorfile)
        except:
            self.mags=get_data(colorfile,0)
            

def testgetcolors(ng=10):
    rolex=watch()
    filter0="HST_ACS_WFC_F775W"
    library="UDF.list"
    seds=get_str(sed_dir+library,0)
    nt=len(seds)
    put_str("test.sim",(["HST_ACS_WFC_F435W",
                         "HST_ACS_WFC_F475W",
                         "HST_ACS_WFC_F555W",
                         "HST_ACS_WFC_F606W",
                         "HST_ACS_WFC_F775W",
                         "HST_ACS_WFC_F850LP"],6*["AB"],ones(6)*24.,ones(6)*5,ones(6)*1.))
    m=linspace(18.,25.,ng)
    z=linspace(0.,6.,ng)
    t=asarray((sin(m)+1.)/2.*nt+1,'int')
    mzt=put_data("test%i.mzt" % ng,(m,z,t))
    # Calculate in the getcolors way
    rolex.set()
    a=getcolors(library,"test.sim",filter0,"test.mzt")
    colors=a.mags
    rolex.check()
    # Calculate the old way
    m,z,t=get_data("test.mzt",(0,1,2))
    ocolors=colors*0.
    filtros=get_str("test.sim",0)
    rolex.set()
    for i in range(len(m)):
        for j in range(len(filtros)):
            print i,j
            ocolors[i,j]=reobs(seds[int(t[i]-1)],m[i],z[i],filter0,z[i],filtros[j])
    rolex.check()

    diff=colors-ocolors
    print diff.sum()
        

def reobs(sed,m=0.,z_0=0.,oldfilter='I_LRIS',
          z_new=0.,newfilter='V_LRIS',cosmology=cosmo,madau='yes',dz_ab=dz_ab):
    """Arguments: sed,m,z_0,oldfilter,z_new,newfilter,cosmology
    Takes a galaxy with m at redshift z_0 in oldfilter, 
    SED=sed and produces its new magnitude in newfilter at z_new. 
    Takes into account cosmological dimming and intergalactic Madau absorption
    The tuple cosmology=(omega,lambda,hubble_constant)
    """
    if sed[-4:]=='.sed': sed=sed[:-4]
    single_z=type(z_new)==type(z_0)
    if single_z:
        if z_0==z_new and oldfilter==newfilter: return m
        z_new=array([z_new])
        
    #Calculate fnew
    model=join([sed,newfilter,'AB'],'.')
    model_path=os.path.join(ab_dir,model)
    #Check whether there are already AB files
    if madau=='yes':
        if os.path.exists(model_path):#model[:-3] in ab_db:
            try:
                zo,f_mod_0=dict_reobs[model]
            except:
                dict_reobs[model]=get_data(model_path,(0,1))
                zo,f_mod_0=dict_reobs[model]
            fnew=match_resol(zo,f_mod_0,z_new)
        else:
            fnew=f_z_sed_AB(sed,newfilter,z_new,'nu',dz_ab=dz_ab)
    else:
        fnew=f_z_sed(sed,newfilter,z_new,units='nu',madau=madau)

    #fnew=where(less(fnew,fnew*ab_clip),99.,fnew) #if the new flux is smaller
    #than ab_clip times the maximum flux, returns 99. (code non-detection) 
    
    #Calculate f_old
    model=join([sed,oldfilter,'AB'],'.')
    model_path=os.path.join(ab_dir,model)
    
    #Check whether there are already AB files
    if madau=='yes':
        if os.path.exists(model_path):#model[:-3] in ab_db:
            try:
                zo,f_mod_0=dict_reobs[model]
            except:
                #print "Read %s" % model 
                dict_reobs[model]=get_data(model_path,(0,1))
                zo,f_mod_0=dict_reobs[model]
            f_old=match_resol(zo,f_mod_0,z_0)
        else:
            f_old=f_z_sed_AB(sed,oldfilter,array([z_0]),units='nu',dz_ab=dz_ab)
    else:
        f_old=f_z_sed(sed,oldfilter,array([z_0]),units='nu')

    # k=2.5*log10(((1.+z_new)/(fnew+1e-100))*((f_old+1e-100)/(1.+z_0)))

    #k=2.5*log10(((f_old+1e-100)/(fnew+1e-100))*(1.+z_new)/(1.+z_0))
    if fnew==0: return 99.
    k=2.5*(log10(f_old)-log10(fnew)+log10(1.+z_new)-log10(1.+z_0))
    
    if single_z:
        if z_0==z_new[0]:
            m_obs=m+k
            return m_obs[0]

    #Distance modulus
    dist=dist_mod(z_new,cosmology)-dist_mod(z_0,cosmology)
    m_obs=m+dist+k
    if single_z:
        return m_obs[0]
    else:
        return m_obs

class getmags:
    def __init__(self,sed,m0,phot_sys0,filter0,filter1,phot_sys1):
        # - sed has to be the name of a SED file contained in $BPZPATH/SED,
        # - phot_sys0 has to be either 'Vega' or 'AB',
        # - filter0 has to be the name of a filter transmission file contained
        #    in $BPZPATH/FILTERS
        # - filter1 is a list containing the names of filter transmission files contained
        #   in $BPZPATH/FILTERS
        # - phot_sys1 is a list containing either 'Vega' or 'AB' as its elements
        nf=len(filter1)
        # If input magnitude is Vega, transform it to AB
        if phot_sys0=="Vega": m0=VegatoAB(m0,filter0)
        self.m=ones(nf)*1.
        for i in range(nf):
            self.m[i]=reobs(sed,m=m0,z_0=0.,oldfilter=filter0,z_new=0.,newfilter=filter1[i])
            if phot_sys1[i]=="Vega": self.m[i]=ABtoVega(m[i],filter1[i])

def color_z(sed,filter_new,filter_old,z=arange(0.,1.5,0.5),calibration='AB',file=None):
    """
    Calculates the color filter_new-filter_old at the redshift vector z.
    It can return the color in Vega or AB calibrations
    Usage:
    gr=color_z('El_cww','g_WFC','r_WFC',z=arange(0.,2.,0.1),'Vega')
    It also works with scalars, e.g.
    gr=color_z('El_cww','g_WFC','r_WFC',1.2)
    """

    try:
        n=len(z)
    except:
        z=array([z])
        n=1
    color=z*0.
    for i in range(len(z)):
        color[i]=reobs(sed,0.,z[i],filter_old,z[i],filter_new)
        
    if calibration=='Vega': color+=ABtoVega(0.,filter_new)-ABtoVega(0.,filter_old)
    
    if file==None:
        if n==1: color=color[0]
        return color
    else: put_data(file,(z,color),header='z     %s-%s(%s) ' % (filter_new,filter_old,calibration))
    
def m_abs(m,z,sed,filter_m,cosmology=cosmo,filter_M="B_Johnson",cal_m="AB",cal_M="Vega"):
    """Arguments: m,z,sed,filter_m,cosmology,filter_M
    If filter2 is used, returns the absolute magnitude 
    in a filter different from the input apparent magnitude"""
    #print z,sed,filter
    mabs=m-dist_mod(z,cosmology)-kcor(z,sed,filter_m)
    #Correct rest frame color
    if filter_M and filter_M<>filter_m:
	mabs=reobs(sed,m=mabs,z_0=0.,oldfilter=filter_m,z_new=0.,newfilter=filter_M)
    if cal_m<>cal_M:
        if cal_m=="AB":
            mabs=ABtoVega(mabs,filter_M)
        elif cal_m=="Vega":
            mabs=VegatoAB(mabs,filter_M)        
    return mabs

M_abs=m_abs

# def mzt_m_abs(m,z,t,lib="eB11.list",
#               filter_m="i_SDSS",filter_M="B_Johnson",cal_m="AB",cal_M="Vega",
#               cosmology=cosmo):
def mzt_m_abs(m,z,t,lib,filter_m,filter_M,cal_m,cal_M,cosmology=cosmo):
    try:
        z*t*m
    except:
        print "LENGTHS OF VECTORS m,z,t INCOMPATIBLE"
    dz=0.001 # step in distance modulus approximation
    xz=z_ab+dz_ab/2.
    nz=len(xz)
    seds=get_lib(lib)
    nt=len(seds)
    # print 'filter_M',filter_M
    mabsdat=lib.split(".")[0]+filter_m+cal_m+filter_M+cal_M+".mabs"
    # print 'mabsdat',mabsdat
    # pausa = raw_input('pepe')
    data=zeros((nz,nt))
    if not os.path.exists(mabsdat):
        for i in range(nt):
            data[:,i]=-kcor(xz,seds[i],filter_m)
            if filter_M<>filter_m:
                data[:,i]+=reobs(seds[i],0.,0.,filter_m,0.,filter_M)
        if cal_M<>cal_m: 
            if cal_M=="Vega" and cal_m=="AB":
                data+=ABtoVega(0.,filter_M)
            elif cal_M=="AB" and cal_m=="Vega":
                data+=VegatoAB(0.,filter_M)

        data=clip(data,-99.,99.)
        put_2Darray(mabsdat,data)
    else:
        data=get_2Darray(mabsdat)
    data+=-dist_mod(xz,cosmology)[:,newaxis]
    d=interpolate.RectBivariateSpline(xz,arange(nt),data,kx=1,ky=1,s=0)
    M=m+d.ev(z,t-1.)
    return M

class sm_grid:
    def __init__(self,lib,filter_m,zmax,cal_m="AB",m0=20.):    
        #cached_sm20="SM20_"+pars.d["SPECTRA"][:-5]+"_"+str(n_interp)+"_"+str(pars.d["DZ"])+"_"+str(pars.d["ZMAX"])+"_"+str(pars.d["FILTER_0"])+".npy" 
        cached_sm20="SM20_%s_%s_%s.npy" % (lib,filter_m,zmax)
        z=N.linspace(1e-4,zmax,1000)
        nz=len(z)
        nt=len(get_lib(lib))
        t=arange(nt)+1.
        if not os.path.isfile(cached_sm20):
            ixz,ixt=N.mgrid[:nz,:nt] 
            ixz=N.ravel(ixz)
            ixt=N.ravel(ixt)
            sm20=mzt_mass_star_T(m0+z[ixz]*0.,z[ixz],t[ixt],lib=lib,
                                 filter_m=filter_m,cal_m=cal_m)
            sm20=N.resize(sm20,(nz,nt))
            N.save(cached_sm20,sm20)
            print "Writing stellar mass information to cache"
        else:
            sm20=N.load(cached_sm20)
            #print "Reading stellar mass information from cache"
        self.sm20=sm20
        self.z=z
        self.t=t

def mzt_m_abs_bak(m,z,t,lib="eB11.list",
                  filter_m="i_SDSS",filter_M="B_Johnson",cal_m="AB",cal_M="Vega",
                  cosmology=cosmo):
    try:
        z*t*m
    except:
        print "LENGTHS OF VECTORS m,z,t INCOMPATIBLE"
    dz=0.001 # step in distance modulus approximation
    #NOTE THAT IT ASSUMES THAT TYPES START AT i=1
    tb0=t-1
    tl=floor(tb0).astype(int)
    tr=ceil(tb0).astype(int)
    dt=tb0-tl
    seds=get_lib(lib)
    ns=len(seds)
    # Subtract distance modulus
    xz=arange(min(z),max(z)+dz,dz)
    M=m-match_resol(xz,dist_mod(xz,cosmology),z)
    for j in range(ns):
        g=equal(tl,j)
        if g.sum()>0.: 
            if j < ns-1:
                # Interpolate between types if required
                M[g]+=-(
                    (1.-dt[g])*kcor(z[g],seds[j],filter_m)+dt[g]*kcor(z[g],seds[j+1],filter_m)
                    )
            else: # Last template
                M[g]+=-kcor(z[g],seds[j],filter_m)            
            
            if filter_M<>filter_m:
                M[g]+=reobs(seds[j],0.,0.,filter_m,0.,filter_M)
    
    if cal_M<>cal_m: 
        if cal_M=="Vega" and cal_m=="AB":
            M+=ABtoVega(0.,filter_M)
        elif cal_M=="AB" and cal_m=="Vega":
            M+=VegatoAB(0.,filter_M)
                        
    return M
    
def v_m_abs(bpzfile,filter_M,
            cal="Vega",
            lib="eB10.list",colsmzt=(9,1,4),
            filter_m="HST_ACS_WFC_F850LP",cosmology=cosmo):
    if type(bpzfile)==type(""):
        m0,zb,tb=get_data(bpzfile,colsmzt)
    else:
        m0,zb,tb=bpzfile
    tb0=tb-1
    tl=floor(tb0).astype(int)
    tr=ceil(tb0).astype(int)
    dt=tb0-tl
    seds=get_lib(lib)
    ns=len(seds)
    # Subtract distance modulus
    xz=arange(min(zb),max(zb)+0.01,0.01)
    M=m0-match_resol(xz,dist_mod(xz,cosmology),zb)

    for j in range(ns):
        g=equal(tl,j)
        print j,g.sum()
        if g.sum()>0.: 
            if j < ns-1:
                # Interpolate between types if required
                M[g]+=-(
                    (1.-dt[g])*kcor(zb[g],seds[j],filter_m)+dt[g]*kcor(zb[g],seds[j+1],filter_m)
                    )
            else: # Last template
                M[g]+=-kcor(zb[g],seds[j],filter_m)            
            
            if filter_M<>filter_m:
                M[g]+=reobs(seds[j],0.,0.,filter_m,0.,filter_M)
    
    if cal=="Vega": 
        M+=ABtoVega(0.,filter_M)
    return M

def test_v_m_abs():
    m=linspace(18,25,100)
    z=linspace(0.1,5.,100)
    t=randint(6,11,size=100)
    put_data("testvmabs.mzt",(m,z,t))
    seds=get_lib("eB10.list")
    s=map(lambda x: seds[x],around(t-1).astype('int'))
    t0=time()
    M=array(map(lambda x: m_abs(m[x],z[x],s[x],filter_m="i_SDSS",filter_M="g_SDSS"),arange(len(m))))
    t1=time()
    print t1-t0
    Mv=v_m_abs("testvmabs.mzt","g_SDSS",cal="AB",lib="eB10.list",colsmzt=(0,1,2),filter_m="i_SDSS")
    print time()-t1
    print "Diff=",mean(Mv-M),"+-",std_mad(Mv-M)

def v_rest_color(bpzfile,filter_left,filter_right,
                 cal="Vega",
                 lib="t200s.list",cols=(1,4)):
    zb,tb=get_data(bpzfile,cols)
    seds=get_lib(lib)
    xt=arange(len(seds))+1.
    ct=xt*0.
    for j in range(len(seds)):
        ct[j]=color_z(seds[j],filter_left,filter_right,0.,calibration=cal)
    return match_resol(xt,ct,tb)

def rest_color(sed,filter_left="U_Johnson",filter_right="B_Johnson",cal="Vega"):
    return color_z(sed,filter_left,filter_right,0.,calibration=cal)

def mzt_rest_color(t,lib="eB11.list",filter_left="g_SDSS",filter_right="i_SDSS",cal="AB"):
    seds=get_lib(lib)
    t0=arange(0.,len(seds))+1.
    c0=zeros(len(seds),"float")
    for i in range(len(seds)):
        c0[i]=color_z(seds[i],filter_left,filter_right,0.,calibration=cal)
    return match_resol(t0,c0,t)

def sun_abs_mag(filter_M,cal="Vega"):
    MB=4.74 #AB, http://www.ucolick.org/~cnaw/sun.html
    MS=MB+color_z("sun_reference",filter_M,"HST_ACS_WFC_F606W",0.,calibration="AB")
    if cal=="Vega":
        return ABtoVega(MS,filter_M)
    else:
        return MS

def solar_luminosity(M,filter_M="B_Johnson",cal="Vega"):
    #Solar luminosity in filter_M
    return mag2flux(M-sun_abs_mag(filter_M,cal=cal))

Lsolar=solar_luminosity
def MfromLsolar(Lsolar,filter_L="B_Johnson",cal="Vega"):
    return flux2mag(Lsolar)+sun_abs_mag(filter_L,cal=cal)

def solar_luminosity_m(m=20.,z=0.1,sed="LRG",filter_m="i_SDSS",cal_m="AB",filter_L="B_Johnson",cosmology=cosmo):
    M=m_abs(m,z=z,sed=sed,filter_m=filter_m,filter_M=filter_L,cal_m=cal_m,cal_M="Vega",cosmology=cosmology)
    return Lsolar(M,filter_M=filter_L,cal="Vega")

Lsolarm=solar_luminosity_m

def mfromLsolar(Lsolar,z,sed="LRG",filter_L="B_Johnson",filter_m="i_SDSS",cal_m="AB"):
    M=MfromLsolar(Lsolar,filter_L=filter_L,cal="Vega")
    return M-m_abs(0.,z,sed,filter_m=filter_m,filter_M=filter_L,cal_M="Vega",cal_m=cal_m)

def halpha(sed,dha=10.,dcomp=100.):
    l,f=get_sed(sed)
    gha=greater(l,6563.-10)*less(l,6563.+10.)
    goha=greater(l,6563.-100)*less(l,6563.+100.)*logical_not(gha)
    return mean(f[gha])/mean(f[goha])

#def mass_star_KG(ABmag=24.,filter_m="i_SDSS",z=1.,sed="Sbc_cww",cosmology=cosmo): 
#    # Very rough estimate of the stellar mass for a template
#    # Based Modified Bell et al. 2003 stellar mass estimate
#    # from Kannappan & Gawiser 2007
#    # log(M/L_K)=....
#    # This is matched to the BC03 results, should be very similar
#    M_k=ABtoVega(m_abs(ABmag,z,sed,filter_m,filter_M="KS_2MASS",cosmology=cosmo),"KS_2MASS")
#    L_k=solar_luminosity(M_k,"KS_2MASS","Vega")
#    BR=rest_color(sed,"B_Johnson","R_Cousins","Vega")
#    if BR> 1.2: logML_k=-0.616+0.34*BR
#    else: logML_k=-0.808+0.5*BR
#    ML_k=10.**(logML_k)
#    Mstar=ML_k*L_k
#    #Additional factors
#    # If dust present (the presence of Halpha is used as proxy), 
#    # decrease M* by 15% (do not include)
#    #ha=halpha(sed)>1.15
#    #if ha or BR<1.3: Mstar*=0.85
#    #Now reduce mass by additional 10% if galaxy is not a massive E/SO
#    #(To account for the presence of starburts)
#    if BR>1.3 and Mstar>1e11: Mstar*=0.9
#    # This factor matches the results to a BC calibration
#    # See Figure 1h from KG07 offset ~1.5-1.8
#    factor=1.6
#    return Mstar/factor

def mass_star_SL(ABmag=24.,
                 filter_m="i_SDSS",z=1.,sed="Sbc_cww",cosmology=cosmo,IMF="Salpeter"): 
    # Rough estimate of the stellar mass for a template
    # Stellar mass estimate from Starlight fits
    if sed[-4:]==".sed": sed=sed[:-4]
    d={}
    d["Ell7_A_0"]=2.354  
    d["Ell6_A_0"]=1.933  
    d["Ell5_A_0"]=1.507  
    d["Ell4_A_0"]=1.233  
    d["ES0_B10"]=2.544 
    d["Sbc_B10"]=1.847
    d["Scd_B10"]=1.900
    d["SB1_B10"]=1.279
    d["SB2_B10"]=0.5502
    d["SB3_B10"]=0.5983
    d["SB11_A_0_l"]=0.09684
    Mi=m_abs(ABmag,z,sed,filter_m,filter_M="i_SDSS",cosmology=cosmology)
    #Li=solar_luminosity_m(m=ABmag,z=z,sed=sed,filter_m=filter_m,cal_m="AB",filter_L="i_SDSS",cosmology=cosmology)
    Li=10.**(-.4*Mi) #We use an AB definition for the luminosity
    gamma=d[sed]
    Mstar=Li*10.**gamma
    if IMF=="Chabrier": Mstar*=10.**(-.2)
    print "SL",sed,gamma,log10(Mstar)-log10(Li)
    return Mstar

def mass_star_T(ABmag=24.,
                filter_m="i_SDSS",z=1.,sed="Sbc_cww",cosmology=cosmo,IMF="Chabrier"): 
    # Rough estimate of the stellar mass for a template
    # Stellar mass estimate from Taylor et al. 2011
    M_i=m_abs(ABmag,z,sed,filter_m,filter_M="i_SDSS",cal_M="AB",cosmology=cosmology)
    #L_i=solar_luminosity(M_i,"i_SDSS","AB")
    L_i=10.**(-.4*M_i) #They use an AB definition for the luminosity
    gi=rest_color(sed,"g_SDSS","i_SDSS","AB")
    gamma=(1.15+0.7*gi)
    Mstar=L_i*10.**gamma
    if IMF=="Salpeter": #Apply typical correction for Salpeter IMF
        Mstar*=10.**(.2)
    print "T",gamma,log10(Mstar)-log10(L_i)
    return Mstar

def mass_star_W2013(ABmag=24.,filter_m="i_SDSS",z=1.,sed="Sbc_cww",
                    cosmology=cosmo,IMF="Salpeter"): 
    # Rough estimate of the stellar mass for a template
    # Wilkins et al. 2013, astro-ph/1304.4421
    # We use the z=0 relationship to be able to work with the rest frame color
    filter1="HST_ACS_WFC_F435W"
    filter2="HST_ACS_WFC_F606W"
    filter_L="HST_ACS_WFC_F606W"
    cal="AB"
    M_V=m_abs(ABmag,z,sed,filter_m,filter_M=filter_L,
              cal_M=cal,cosmology=cosmology)
    L_V=solar_luminosity(M_V,filter_L,cal)    
    color=rest_color(sed,filter1,filter2,cal=cal)
    gamma=10.**(-0.5+1.1*color)
    Mstar=gamma*L_V
    if IMF=="Chabrier": Mstar*=10.**-.2
    return Mstar

def mzt_mass_star_T(m,z,t,lib="eB11.list",
                    filter_m="i_SDSS",cal_m="AB",cosmology=cosmo,IMF="Chabrier"):
    # Assumes that the BPZ mag is AB
    # We assume t is BPZ-output like, i.e. t=1,2,3,...
    try:
        z*t*m
    except:
        print "LENGTHS OF VECTORS m,z,t NOT COMPATIBLE"
    # print "Calculating absolute magnitude"
    filter_M="i_SDSS" # "HST_ACS_WFC_F814W"
    cal_M="AB"
    M_i=mzt_m_abs(m,z,t,lib=lib,
                  filter_m=filter_m,cal_m=cal_m,
        filter_M=filter_M,cal_M=cal_M,cosmology=cosmology)
    M_i_sun=sun_abs_mag(filter_M=filter_M,cal=cal_M)
    L_i=10.**(-.4*M_i) #They use an AB definition for the luminosity
    #L_i=mag2flux(M_i-M_i_sun)
    # print "Calculating rest frame color"
    gi=mzt_rest_color(t,lib=lib,filter_left="g_SDSS",filter_right="i_SDSS",cal="AB")
    # Mstar=10.**(1.15+0.7*gi)*L_i #Original from Txitxo
    Mstar_red =10.**(1.160+0.694*gi)*L_i ## For red galaxies. #Alberto
    Mstar_blue=10.**(1.336+0.447*gi)*L_i ## for blue galaxies.#Alberto
    Mstar = N.where(t<5.6,Mstar_red,Mstar_blue)
    if IMF=="Salpeter": Mstar*=10.**.2
    return Mstar

def mass_star_T_cz(mL=25.,mR=24.,z=0.187,
                   filter_L="HST_ACS_WFC_F475W",
                   filter_R="HST_ACS_WFC_F814W",
                   lib="B14.list",cal="AB",mmin=0.):
    # Avoid absurd values
    gR=logical_not(isnan(mR)+isinf(mR))*greater(mR,mmin)
    gL=logical_not(isnan(mL)+isinf(mL))*greater(mL,mmin)
    g=gR*gL
    sm=mR*0.
    t=mR*0.
    # Get the type
    seds=get_lib(lib)
    ts=arange(len(seds))
    C0=ts*0.
    for i in ts:
        C0[i]=color_z(seds[int(i)],filter_L,filter_R,z,calibration=cal)
    t[g]=clip(match_resol(ts+1.,C0,mL[g]-mR[g]),1.,len(ts)+1.)
    sm[g]=mzt_mass_star_T(mR[g],t[g]*0.+z,t[g],lib=lib,filter_m=filter_R,cal_m=cal)
    sm[isnan(sm)+isinf(sm)]=0.
    return sm

def Phi_Ms_M2013(Ms,z,t,lib="eB11",use_phi=1,use_alpha=1):
    #Stellar mass function from the Muzzin et al. paper
    #Assumes we are using the eB11.list
    
    if lib.split(".")[0]=="eB11":
        q=less_equal(t,6.25)
        sf=greater(t,6.25)
    elif lib.split(".")[0]=="B14":
        q=less_equal(t,3.25)
        sf=greater(t,3.25)        

    phi=zeros(len(Ms),"float")
    xz=array((0.38,0.78,1.25,1.75,2.25,2.75,3.5))
    #z=clip(z,min(xz),10.)
    #Quiescent galaxies 
    msq=interpolate.interp1d(xz,
                             array((10.75,10.84,10.83,10.80,10.79,10.81,11.00)),"slinear")(clip(z[q],min(xz),max(xz)))
    fsq=interpolate.interp1d(xz,
                             array((30.65,14.38, 7.48, 3.61, 1.14, 0.66, 0.05))*1e-4,"slinear")(clip(z[q],min(xz),max(xz)))
    aq=-.4

    #SF galaxies 
    mssf=interpolate.interp1d(xz,
                              array((10.75,10.82,10.82,10.91,11.03,11.14,11.47)),"slinear")(clip(z[sf],min(xz),max(xz)))
    fssf=interpolate.interp1d(xz,
                              array((13.58,10.95, 7.20, 4.49, 2.01, 1.09, 0.08))*1e-4,"slinear")(clip(z[sf],min(xz),max(xz)))
    asf=-1.3    

    #hist(concatenate((msq,mssf)))
    #show()
    
    phi[q]= log(10.)*exp(-10.**(Ms[q] -msq))
    phi[sf]=log(10.)*exp(-10.**(Ms[sf]-mssf))

    #You can use only the exponential part, but then the overall normalization will be wrong
    if use_alpha:
        phi[q]*= 10.**((Ms[q]- msq )*(1.+aq ))
        phi[sf]*=10.**((Ms[sf]-mssf)*(1.+asf))

    if use_phi:
        phi[q]*=fsq
        phi[sf]*=fssf
        
    return phi

#def nzt_mass(xz,xt,m,filter_m="HST_ACS_WFC_F814W",
#             lib="eB11.list",area=1.,use_phi=1,use_alpha=1):
#    xtg,xzg=meshgrid(xt,xz)
#    #t=xtg.ravel(order="F")
#    #z=xzg.ravel(order="F")
#    t=xtg.ravel()
#    z=xzg.ravel()
#    ms=log10(mzt_mass_star_T(m*ones(len(z),"float"),z,t,lib=lib,filter_m=filter_m))
#    phi=resize(Phi_Ms_M2013(ms,z,t,lib=lib,use_alpha=use_alpha,use_phi=use_phi),
#               (len(xz),len(xt)))*dVc(xz)[:,newaxis]
#    return phi/phi.sum()

def nzt_mass(xz,xt,m,filter_m="HST_ACS_WFC_F814",
                  lib="eB11.list",area=1.,use_phi=1,use_alpha=1):    
    dz=median(xz[1:]-xz[:-1])
    zmin=min(xz)
    zmax=max(xz)
    sm=sm_grid(lib,filter_m,zmax)
    # print 'sm.z, sm.t', sm.z, sm.t
    ms=interpolate.RectBivariateSpline(sm.z,sm.t,log10(sm.sm20*mag2flux(m-20.)),kx=1,ky=1,s=0)(xz,xt)
    nt=len(xt)
    nz=len(xz)
    
    phi=resize(Phi_Ms_M2013(ravel(ms),
                            ravel(outer(xz,xt*0+1.)),
                            ravel(outer(xz*0+1.,xt)),                     
                            lib=lib,use_alpha=use_alpha,use_phi=use_phi),ms.shape)*dVc(xz)[:,newaxis]
    if not use_phi*use_alpha:
        phi=phi/phi.sum(0)

    # print 'phi', phi
    #plot(xz,phi)
    #show()

    #phi=resize(Phi_Ms_M2013(ms,xz,xt,lib=lib,use_alpha=use_alpha,use_phi=use_phi),
    #           (len(xz),len(xt)))*dVc(xz)[:,newaxis]
    return phi/phi.sum()
    
# Stellar Mass Function from the Marchesini paper
def Phi_Ms_2D(z,Ms):
    z=clip(z,0.1,3.5)
    xz=array((0.1,1.65,2.5,3.5))
    a0=array((-1.18,-0.99,-1.01,-1.39))
    ms0=array((10.96,10.91,10.96,11.38))
    fs0=array((30.87,10.17,3.95,0.53))*1e-4
    #alfa=match_resol(xz,a0,z)
    alfa=-1.
    ms=match_resol(xz,ms0,z)
    fs=match_resol(xz,fs0,z)
    return (log(10.)*fs[:,newaxis]*
            10**((Ms[newaxis,:]-ms[:,newaxis])*(1.+alfa))*
        exp(-10**(Ms[newaxis,:]-ms[:,newaxis])))

def Phi_Ms(z,Ms):
    if z.shape > Ms.shape:
        print "Shapes are different!!!"
        return 0
    z=clip(z,0.1,3.5)
    xz=array((0.1,1.65,2.5,3.5))
    a0=array((-1.18,-0.99,-1.01,-1.39))
    ms0=array((10.96,10.91,10.96,11.38))
    fs0=array((30.87,10.17,3.95,0.53))*1e-4
    #alfa=match_resol(xz,a0,z)
    alfa=-1.
    ms=match_resol(xz,ms0,z)
    fs=match_resol(xz,fs0,z)
    return log(10.)*fs*10**((Ms-ms)*(1.+alfa))*exp(-10**(Ms-ms))
        


def mass_star_B(ABmag=24.,
                filter_m="i_SDSS",z=1.,sed="Sbc_cww",cosmology=cosmo): 
    # Very rough estimate of the stellar mass for a template
    # Stellar mass estimate from Bell et al. 2003
    #Do not use it; DEPRECATED  
    M_i=m_abs(ABmag,z,sed,filter_m,filter_M="i_SDSS",cosmology=cosmology)
    gi=rest_color(sed,"g_SDSS","i_SDSS","AB")
    Mstar=10.**(1.68+0.518*gi-0.4*M_i)
    return Mstar

def mzt_mass_star_B(m,z,t,lib="eB11.list",
                    filter_m="i_SDSS",cal_m="AB",cosmology=cosmo):
    # Assumes that the BPZ mag is AB
    # We assume t is BPZ-output like, i.e. t=1,2,3,...

    try:
        z*t*m
    except:
        print "LENGTHS OF VECTORS m,z,t NOT COMPATIBLE"
    M_i=mzt_m_abs(m,z,t,lib=lib,filter_m=filter_m,cal_m=cal_m,filter_M="i_SDSS",cal_M="AB",cosmology=cosmology)
    gi=mzt_rest_color(t,lib=lib,filter_left="g_SDSS",filter_right="i_SDSS",cal="AB")
    Mstar=10.**(1.68+0.518*gi-0.4*M_i)
    return Mstar


def mzt_mass_star_KG(m,z,t,lib="eB10.list",
                     filter_m="i_SDSS",cal_m="AB",cosmology=cosmo):
    # Very rough estimate of the stellar mass for a template
    # Based Modified Bell et al. 2003 stellar mass estimate
    # from Kannappan & Gawiser 2007
    # log(M/L_K)=....
    # Assumes that the BPZ mag is AB
    # Matches Lin et al. 2007 within 0.066 dex for the COSMOS catalog
    # We assume t is BPZ-output like, i.e. t=1,2,3,...

    try:
        z*t*m
    except:
        print "LENGTHS OF VECTORS m,z,t NOT COMPATIBLE"

    print m,z,t
    print "t=",t
    M_k=mzt_m_abs(m,z,t,lib=lib,filter_m=filter_m,cal_m=cal_m,filter_M="KS_2MASS",cal_M="Vega",cosmology=cosmology)
    print "t=",t
    seds=get_lib(lib)
    print "M_k",M_k
    L_k=mag2flux(M_k)*solar_luminosity(0.,"KS_2MASS","Vega")
    print "L_k=",L_k
    # Template data
    logML=zeros(len(seds),"float")
    ha=zeros(len(seds),"float")
    BR=zeros(len(seds),"float")
    sb=zeros(len(seds),"float")

    for j in range(len(seds)):
        ha[j]=halpha(seds[j])>1.15
        BR[j]=rest_color(seds[j],"B_Johnson","R_Cousins","Vega")
        print BR[j]
        if BR[j]> 1.2: logML[j]=-0.616+0.34*BR[j]
        else: logML[j]=-0.808+0.5*BR[j]
        #sb[j]=ha[j] or BR[j] < 1.3
        # If dust present (the presence of Halpha is used as proxy), 
        # decrease M* by 15%
        # Do not apply this correction; results are worse with it
        #if sb[j]: logML[j]+=-0.15
        #if ha[j]: logML[j]+=-0.15
        #if sb[j]: logML[j]+=0.
    #Now take into account the interpolation between templates

    it=array(map(lambda x: int(round(x-1)),t))
    tl=floor(t-1.).astype(int)
    tr=ceil(t-1.).astype(int)
    dt=t-1.-tl
    logML_k=((1.-dt)*take(logML,tl)+dt*take(logML,tr))

    Mstar=10.**(logML_k)*L_k

    # Now reduce mass by additional 10% if galaxy is not a massive E/SO
    # (To account for the presence of starburts)
    smallorblue=take(sb,it)+less(Mstar,1e11)
    Mstar[smallorblue.astype(int)]*=0.9
    # This factor matches the results to a BC calibration, also to Lin et al. 2007
    # See Figure 1h from KG07 offset ~1.5-1.8
    factor=1.6
    return Mstar/factor


def v_mass_star_KG(bpzfile,lib="CWWSB_B2004a.list",cols=(1,4,9),
                  filter="i_SDSS",cosmology=cosmo):
    # Very rough estimate of the stellar mass for a template
    # Based Modified Bell et al. 2003 stellar mass estimate
    # from Kannappan & Gawiser 2007
    # log(M/L_K)=....
    # Assumes that the BPZ mag is AB
    # Matches Lin et al. 2007 within 0.066 dex for the COSMOS catalog

    print "Careful with the cols parameter"
    print "They correspond to z,t,m"
    
    print cols
    M_k=v_m_abs(bpzfile,"KS_2MASS","Vega",lib,cols,filter,cosmology)
    z,t,m=get_data(bpzfile,cols)
    seds=get_lib(lib)

    L_k=mag2flux(M_k)*solar_luminosity(0.,"KS_2MASS","Vega")

    # Template data
    logML=zeros(len(seds),"float")
    ha=zeros(len(seds),"float")
    BR=zeros(len(seds),"float")
    sb=zeros(len(seds),"float")
    for j in range(len(seds)):
        ha[j]=halpha(seds[j])>1.15
        BR[j]=rest_color(seds[j],"B_Johnson","R_Cousins","Vega")
        if BR[j]> 1.2: logML[j]=-0.616+0.34*BR[j]
        else: logML[j]=-0.808+0.5*BR[j]
        #sb[j]=ha[j] or BR[j] < 1.3
        # If dust present (the presence of Halpha is used as proxy), 
        # decrease M* by 15%
        # Do not apply this correction; results are worse with it
        #if sb[j]: logML[j]+=-0.15
        #if ha[j]: logML[j]+=-0.15
        #if sb[j]: logML[j]+=0.
    #Now take into account the interpolation between templates
    it=array(map(lambda x: int(round(x-1.)),t))
    tl=floor(t-1.).astype(int)
    tr=ceil(t-1.).astype(int)
    dt=t-1.-tl
    logML_k=((1.-dt)*take(logML,tl)+dt*take(logML,tr))
    Mstar=10.**(logML_k)*L_k

    # Now reduce mass by additional 10% if galaxy is not a massive E/SO
    # (To account for the presence of starburts)
    smallorblue=take(sb,it)+less(Mstar,1e11)
    Mstar[smallorblue.astype(int)]*=0.9
    # This factor matches the results to a BC calibration, also to Lin et al. 2007
    # See Figure 1h from KG07 offset ~1.5-1.8
    factor=1.6
    return Mstar/factor

def p_mass_star_KG(p_i,mo,filter="i_SDSS",lib="CWWSB_B2004a.list",cosmology=cosmo):
    # Very rough estimate of the stellar mass for a template
    # Based Modified Bell et al. 2003 stellar mass estimate
    # from Kannappan & Gawiser 2007
    # log(M/L_K)=....
    # Assumes that the BPZ mag is AB
    # Matches Lin et al. 2007 within 0.066 dex for the COSMOS catalog

    print "Careful with the cols parameter"
    print "They correspond to z,t,m"
    
    print cols
    M_k=m_abs

    M_k=v_m_abs(bpzfile,"KS_2MASS","Vega",lib,cols,filter,cosmology)
    z,t,m=get_data(bpzfile,cols)
    seds=get_lib(lib)

    L_k=mag2flux(M_k)*solar_luminosity(0.,"KS_2MASS","Vega")

    # Template data
    logML=zeros(len(seds),"float")
    ha=zeros(len(seds),"float")
    BR=zeros(len(seds),"float")
    sb=zeros(len(seds),"float")
    for j in range(len(seds)):
        ha[j]=halpha(seds[j])>1.15
        BR[j]=rest_color(seds[j],"B_Johnson","R_Cousins","Vega")
        if BR[j]> 1.2: logML[j]=-0.616+0.34*BR[j]
        else: logML[j]=-0.808+0.5*BR[j]
        #sb[j]=ha[j] or BR[j] < 1.3
        # If dust present (the presence of Halpha is used as proxy), 
        # decrease M* by 15%
        # Do not apply this correction; results are worse with it
        #if sb[j]: logML[j]+=-0.15
        #if ha[j]: logML[j]+=-0.15
        #if sb[j]: logML[j]+=0.
    #Now take into account the interpolation between templates
    it=array(map(lambda x: int(round(x-1.)),t))
    tl=floor(t-1.).astype(int)
    tr=ceil(t-1.).astype(int)
    dt=t-1.-tl
    logML_k=((1.-dt)*take(logML,tl)+dt*take(logML,tr))
    Mstar=10.**(logML_k)*L_k

    # Now reduce mass by additional 10% if galaxy is not a massive E/SO
    # (To account for the presence of starburts)
    smallorblue=take(sb,it)+less(Mstar,1e11)
    Mstar[smallorblue.astype(int)]*=0.9
    # This factor matches the results to a BC calibration, also to Lin et al. 2007
    # See Figure 1h from KG07 offset ~1.5-1.8
    factor=1.6
    return Mstar/factor

def v_mass_star_lin(bpzfile,lib="CWWSB_B2004a.list",cols=(1,4,9),
                  filter="i_SDSS",cosmology=cosmo):
    # Very rough estimate of the stellar mass for a template
    # Lin et al. 2007
    # http://adsabs.harvard.edu/abs/2007ApJ...660L..51L
    # Agrees within 0.25 dex with Bundy et al. 2006

    M_B=v_m_abs(bpzfile,"B_Johnson","Vega",lib,cols,filter,cosmology)
    z,t,m=get_data(bpzfile,cols)
    seds=get_lib(lib)

    # Template data
    BV_t=zeros(len(seds),"float")
    UB_t=zeros(len(seds),"float")
    BR_t=zeros(len(seds),"float")
    logMt=zeros(len(seds),"float")
    for j in range(len(seds)):
        BV_t[j]=rest_color(seds[j],"B_Johnson","V_Johnson","Vega")
        UB_t[j]=rest_color(seds[j],"U_Johnson","B_Johnson","Vega")
        BR_t[j]=rest_color(seds[j],"B_Johnson","R_Cousins","Vega")
        logMt[j]=-0.4*-5.48+1.737*BV_t[j]+0.098*UB_t[j]-0.130*UB_t[j]**2-1.003

    it=array(map(lambda x: int(round(x-1.)),t))
    tl=floor(t-1.).astype(int)
    tr=ceil(t-1.).astype(int)
    dt=t-1.-tl
    #BV=((1.-dt)*take(BV_t,tl)+dt*take(BV_t,tr))
    #UB=((1.-dt)*take(UB_t,tl)+dt*take(UB_t,tr))
    #BR=((1.-dt)*take(BR_t,tl)+dt*take(BR_t,tr))
    
    logMstar=((1.-dt)*take(logMt,tl)+dt*take(logMt,tr))-0.4*M_B-0.268*z
    
    return 10.**logMstar

#v_mass_star=v_mass_star_lin
v_mass_star=v_mass_star_KG

def test_v_mass_star(bpz="test.bpz",lib="frwks0_140.list",cols=(1,4,10)):
    # Various galaxy data
    M_K=v_m_abs(bpz,"K_Johnson","Vega",lib,cols,"i_SDSS")
    M_B=v_m_abs(bpz,"B_Johnson","Vega",lib,cols,"i_SDSS")
    seds=get_lib(lib)
    z,t,m=get_data(bpz,cols)
    g=greater(M_B,-30.)
    m,z,t,M_B,M_K=multicompress(g,(m,z,t,M_B,M_K))    
    BV_t=zeros(len(seds),"float")
    UB_t=zeros(len(seds),"float")
    for j in range(len(seds)):
        BV_t[j]=rest_color(seds[j],"B_Johnson","V_Johnson","Vega")
        UB_t[j]=rest_color(seds[j],"U_Johnson","B_Johnson","Vega")

    it=array(map(lambda x: int(round(x-1.)),t))
    tl=floor(t-1.).astype(int)
    tr=ceil(t-1.).astype(int)
    dt=t-1.-tl
    BV=((1.-dt)*take(BV_t,tl)+dt*take(BV_t,tr))
    UB=((1.-dt)*take(UB_t,tl)+dt*take(UB_t,tr))

    print "lms_KG"
    lms_KG=log10(v_mass_star_KG(bpz,lib,cols))
    print "lms_lin"
    lms_lin=log10(v_mass_star_lin(bpz,lib,cols))
    show()
    d=lms_lin-lms_KG
    lms_lin,lms_KG,d=multicompress(g,(lms_lin,lms_KG,d))

    print "median d",median(d)
    print "std mad d",std_mad(d)
    print "std d",std(d)
    print "std_mad(d+0.268*z)",std_mad(d+0.268*z)
    print "median(d+0.268*z)",median(d+0.268*z)
    plot(t,d,"s")
    plot(t,d+0.268*z- 0.11,"s")
    show()

    
def L(M):
    """Calculates the luminosity in ergs/s/Hz
      for a galaxy of known absolute magnitude"""
    return fnu(M)*4.*pi*(10.*3.08568025*1e18)**2


def M2300(sed,z=1.,ABmag=20.,filter="i_SDSS"):
    """Calculates the luminosity at 2300AA 
    for a galaxy of known redshift and AB magnitude"""
    x,y=get_sed(sed)
    yz=redshift(x,y,z)
    f=normalize(x,yz,m=ABmag,filter=filter)
    f2300=match_resol(x,f,2300.)
    m2300=AB(f2300)
    return m2300-dist_mod(z)

def SFR_Muv(sed,z,ABmag,filter="i_SDSS"):
    """From Kenicutt 1998, ARAAA"""
    return 1.4e-28*L(M2300(sed,z,ABmag))

def f_OII(M_uv,z):
    """From Kennicutt 1998, flux in OII (ergs/s/cm2 from SFR (through M_2300AA"""
    return 10.**(-0.4*M_uv+10.65-dist_mod(z)/2.5)*1e-17

def addlines(sed,dry=1,factor=1.,plots=0,erase=1):
    """ Add emission lines to a template following a procedure similar to that of Ilbert et al. 2009 """
    #wavs=array([1216.,3727.30,4861.33,5006.,6563.])
    #ratios=array((2.,1.,0.61,0.36,1.77))


    #http://www.sdss.org/dr7/algorithms/linestable.html    
    wavs=array([2799.117,3727.092,4102.89,4341.68,4862.68,4960.295,5008.240,6549.86,6564.61,6585.27,6718.29,6732.67])
    ratios=array([1.0,5.0,0.5,1.0,2.0,2.0,3.0,3.0,8.0,3.0,3.0,3.0])
    ratios*=1./5.

    #Let's build the lines as gaussians with dl=250./3e5*l (250 is the velocity dispersion of a M=-21 galaxy)
    sigma=250./3e5*wavs
    #res=ratios*0.
    x,y=get_sed(sed)
    # First normalize template to M=-21 in the B band
    x0=x
    fnorm=normalize(x,y,-21.,"B_Johnson",units="lambda")
    y=fnorm
    f2300=match_resol(x,fl2fnu(x,fnorm),2300.)
    m2300=AB(f2300)
    foii=10.**(-0.4*m2300+10.65)*1e-17
    #print foii
    f=foii*ratios
    for i in range(len(wavs)):
        xx=arange(-10.*sigma[i]+wavs[i],10.*sigma[i]+wavs[i]+1.,1.)
        yy=factor*f[i]/sqrt(2.*pi)/sigma[i]*exp(-(xx-wavs[i])**2/2./sigma[i]**2)
        i1=searchsorted(x,xx[0])
        i2=searchsorted(x,xx[-1])
        nyy=match_resol(x,y,xx)
        if erase:
            if mean(nyy)>(y[i1]+y[i2])*0.5: #i.e. if there is an emission line
                nyy=nyy*0.+(y[i1]+y[i2])*0.5
        x=concatenate([x[:i1],xx,x[i2:]])
        y=concatenate((y[:i1],nyy+yy,y[i2:]))

    
    if plots:
        plot(x0,fnorm)
        plot(x,y,"r")
        axis((1000.,7000.,0.,10.))
        show()


    if sed[-4:]==".sed":
        name=sed[:-4]+"_l.sed"
    else:
        name=sed+"_l.sed"

    print "Saving data to %s " % name
    if not dry:  put_data(name,(x,y))


def mag_corr(filter="g_SDSS",
             refsed="LRG.sed",
             corrsed="LRG_svd1.sed",
             z=0.16,
             a=1.):
    """
    Calculate the magnitude correction to a given spectrum refsed at redshift z,
    observed in a given filter if we add to it a correction spectrum multiplied by a factor a
    """
               
    fo=f_z_sed_AB(refsed,filter,z,units="nu")
    fn=f_z_sed_AB(corrsed,filter,z,units="nu")*a+fo
    return flux2mag(fn/fo)

class mag_corr:
    """
    Calculate the magnitude correction to a given spectrum refsed at redshift z,
    observed in a given filter if we add to it a correction spectrum multiplied by a factor a
    """               
    def __init__(self,filtro="g_SDSS",
                 refsed="LRG.sed",
                 corrsed="LRG_svd1.sed",
                 z=0.16,
                 c=1.):        
        self.fo=f_z_sed_AB(refsed,filtro,z,units="nu")
        self.df=f_z_sed_AB(corrsed,filtro,z,units="nu")
        self.mag_corr=flux2mag(1.+c*self.df/self.fo)

class Rectify:
    """This class takes an observed spectrum defined by
    x_spec,y_spec, with observed magnitudes""" 
    def __init__(self,x_spec,y_spec,z,
                 mags=[20.,20.,20.,20.,20.],
                 errmags=[0.1,0.1,0.1,0.1,0.1],
                 filters=['HST_ACW_WFC_F435W',
                          'HST_ACW_WFC_F475W',
                          'HST_ACW_WFC_F625W',
                          'HST_ACW_WFC_F775W',
                          'HST_ACW_WFC_F850LP'],
                 calibration=['AB','AB','AB','AB','AB'],
                 spec_library='CWWSB_B2004a.list'):
        self.z=z
        self.x_spec=x_spec
        self.y_spec=y_spec
        self.mags=mags
        self.errmags=errmags
        self.filters=filters
        self.calibration=calibration
        self.spec_library=spec_library
    def run(self):
        # 1. Calcular el template que mejor aproxima los colores observados
        # 2., 
        pass
    
def normpz_hdfn(m):
    x=arange(20.,28.25,0.25)
    y=array([0.086534672316250477, 0.098785977916272596, 0.11215243920861845, 0.1268597519035902, 0.14344744528968084, 0.16249029793737405, 0.18449216283357278, 0.2099790708309997, 0.23955764300686597, 0.27393028403326586, 0.31389628183448337, 0.36035059976587069, 0.4142831278917275, 0.47677864948657256, 0.54901731087036942, 0.63227538264274508, 0.72792614776526487, 0.83744070201018816, 0.96238816055967302, 1.1044340972182007, 1.265335068984103, 1.4469262177643025, 1.6510988076404765, 1.8797655567349656, 2.1348136443629042, 2.4180477202928898, 2.7311273523843784, 3.0755045453151841, 3.4523670387023242, 3.8625921905605574, 4.3067147255294973, 4.7849098803394785, 5.2969918482574201])
    return match_resol(x,y,m)

def pz_hdfn(z,m=24.):
    #For a number of galaxies, return their probality of being at certain
    #Redshift given their magnitude in the AB i band
    momin=20.
    momax=28.
    m=clip(m,momin,momax)
    a=array((2.465,1.806,1.806,0.906,0.906,0.906))
    zo=array((0.431,0.390,0.390,0.0626,0.0626,0.0626))
    km=array((0.0913,0.0636,0.0636,0.123,0.123,0.123))
    fo_t=array((0.35,0.25,0.25))
    k_t=array((0.450,0.147,0.147))
    zt_at_a=numpy.power.outer(z,a)
    f_t=zeros((len(a),),'float')
    dm=clip(m-momin,0.,8.)
    f_t[:3]=fo_t*exp(-k_t*dm)
    f_t[3:]=(1.-add.reduce(f_t[:3]))/3.
    zmt=clip(zo+km*dm,0.01,15.)
    return add.reduce(z**a*exp(-z**a/zmt**2),-1)/normpz_hdfn(m)

def vpz_hdfn(z,m):
    try:
        len(m)
    except:
        m=zeros(len(z))+m
    salida=zeros(len(z))
    for i in range(len(z)):
        salida[i]=pz_hdfn(z[i],m[i])
    return salida


#Incluir algo que pase de magnitudes absolutas a luminosidades solares
#def luminosity(m,filter):
#    Basicamente hallar la magnitud absoluta del sol en el filtro que corresponda
#    m_sun y normalizar
#    return m-dist_mod(z,cosmo[0],cosmo[1],cosmo[2])

#COSMOLOGICAL DISTANCES

#def dl_lambda(z,omega=.3,h=1.):
#    """Aproximation for the luminosity distance 
#    for flat cosmologies with cosmological constant
#    ApJSS, Ue Li Pen 120:4950, 1999"""
#    
#    if omega<0.2: raise """omega less than 0.2: outside
#    parameter range for the aproximation"""
#    
#    if h>1. or h<.4 :
#	print "Wrong value for h",h
#	sys.exit()
#	
#    def eta(a,om):
#	s=((1.-om)/om)**(1./3.)
#	return 2.*sqrt(s**3+1.)*\
#	       (a**(-4)-0.1540*s/a**3+0.4304*s**2/a**2+
#		0.19097*s**3/a+0.066941*s**4)**(-1./8.)
#	
#    return 2.9979*1e5/(h*100.)*(1.+z)*\
#	       (eta(1.,omega)-eta(1./(1.+z),omega))
#
#def dl_nolambda(z,omega=.3,h=.7):
#    """Luminosity distance for a lambda=0
#    universe"""
#    cosa=sqrt(1.+omega*z)
#    return 2.9979*1e5/(h*100.)*z*\
#	   (1.+cosa+z)/(1.+cosa+omega*z/2.)

#def dl(z,cosmology=(.3,.7,.7)):
#    omega,l,h=cosmology
#    if l>0.: 
#	if l+omega<>1.: raise 'lambda>0 but no flat cosmology!'
#	return dl_lambda(z,omega,h)
#    if l==0: return dl_nolambda(z,omega,h)
#    if l<0: raise 'lambda<0!!'
#    
#def da(z,cosmology=(.3,.7,.7)):
#    return dl(z,cosmology)/(1.+z)**2

#######################################################################
######New distance definitions. Hogg 1999, astro-ph/9905116############
###### Tested with Ned Wright's calculator ############################
#######################################################################

def dh(cosmology=cosmo,units="Mpc"):
    #Hubble distance (units can be m or Mpc)
    h=cosmology[2]
    if units=="Mpc":
        return 3000./h
    else:
        return 9.26e25/h

def th(cosmology=cosmo,tunits="yr"):
    #Hubble time (units can be s or yr)
    h=cosmology[2]
    if tunits=="yr":
        return 9.78e9/h
    else:
        return 3.09e17/h

def omega_k(cosmology=cosmo):
    return 1.-cosmology[0]-cosmology[1]

def e(z,cosmology=cosmo):
    o_k=omega_k(cosmology)
    o_m=cosmology[0]
    o_l=cosmology[1]
    return sqrt(o_m*(1.+z)**3+o_k*(1.+z)**2+o_l)

def dc(z,cosmology=cosmo,units="Mpc"):
    #Comoving distance
    def f(x): return 1./e(x,cosmology)
    try:
        len(z)
        uno=0
    except:
        uno=1
    
    if uno:
        return dh(cosmology,units)*integrate.quad(f,0.,z)[0]
    else:
        xz=concatenate(([0.],z))
        res=xz*0.
        for i in range(1,len(xz)):
            res[i]=dh(cosmology,units)*integrate.quad(f,xz[i-1],xz[i])[0]
        return add.accumulate(res)[1:]

def dm(z,cosmology=cosmo,units="Mpc"):
    #Transverse comoving distance
    o_k=omega_k(cosmology)
    if o_k>0:
        res= (dh(cosmology,units)/sqrt(o_k)*
              sinh(sqrt(o_k)*dc(z,cosmology,units)/dh(cosmology,units)))
    elif o_k==0:
        res= dc(z,cosmology,units)
    else:
	res= (dh(cosmology,units)/sqrt(-o_k)*
              sin(sqrt(-o_k)*dc(z,cosmology,units)/dh(cosmology,units)))
    return res

def da(z,cosmology=cosmo,units="Mpc"):
    return dm(z,cosmology,units)/(1.+z)

def dl(z,cosmology=cosmo,units="Mpc"):
    return dm(z,cosmology,units)*(1.+z)

def dVc(z,cosmology=cosmo,units="Mpc",fov=41252.9612):
    return dh(cosmology,units)*(1.+z)**2*da(z,cosmology,units)**2/e(z,cosmology)*fov/(4.*pi*(180./pi)**2)

def Vc(z,cosmology=cosmo,units="Mpc",fov=1.):
    uno=0
    try: len(z)
    except: z=array([z]);uno=1
    
    f=fov/(4.*pi*(180./pi)**2) #Convert to fraction of the full sphere
    #Comoving volume
    k=omega_k(cosmology)
    Dh=dh(cosmology,units)
    sk=sqrt(abs(k))
    y=z*0.
    for i in range(len(z)):
        Dm=dm(z[i],cosmology,units)
        if k>0:
            vc=4.*pi*Dh**3/(2.*k)*(Dm/Dh*sqrt(1.+k*Dm**2/Dh**2)-
                                   1./sk*arcsinh(sk*Dm/Dh))
        elif k==0:
            vc=4.*pi/3.*Dm**3
        else:
            vc=4.*pi*Dh**3/(2.*k)*(Dm/Dh*sqrt(1.+k*Dm**2/Dh**2)-
                                   1./sk*arcsin(sk*Dm/Dh))
        y[i]=vc*f
    if uno: return y[0]
    else: return y

def tl(z,cosmology=cosmo,tunits="yr"):
    # Lookback time
    def f(x): return 1./((1.+x)*e(x,cosmology))
    return th(cosmology,tunits)*integrate.quad(f,0.,z)[0]

def testHogg():
    x=arange(0.,5.,0.01)
    # Fig 1
    def f1(x):
        y=x*0.
        for i in range(len(x)): y[i]=dm(x[i],(1.,0.,1.))/dh((1.,0.,1.))
        return y
    def f2(x):
        y=x*0.
        for i in range(len(x)): y[i]=dm(x[i],(0.05,0.,1.))/dh((0.05,0.,1.))
        return y
    def f3(x):
        y=x*0.
        for i in range(len(x)): y[i]=dm(x[i],(0.2,0.8,1.))/dh((0.2,0.8,1.))
        return y

    plot(x,f1(x),x,f2(x),"-",x,f3(x),".")
    xticks(arange(0.,5.2,0.2))
    yticks(arange(0.,3.2,0.2))
    show()

    # Fig 2
    def f1(x):
        y=x*0.
        for i in range(len(x)): y[i]=da(x[i],(1.,0.,1.))/dh((1.,0.,1.))
        return y
    def f2(x):
        y=x*0.
        for i in range(len(x)): y[i]=da(x[i],(0.05,0.,1.))/dh((0.05,0.,1.))
        return y
    def f3(x):
        y=x*0.
        for i in range(len(x)): y[i]=da(x[i],(0.2,0.8,1.))/dh((0.2,0.8,1.))
        return y

    plot(x,f1(x),x,f2(x),"-",x,f3(x),".")
    xticks(arange(0.,5.2,0.2))
    yticks(arange(0.,0.52,0.02))
    show()



######################## Hogg 1999 ########################################

class vc:
    def __init__(self,z=0.57,sed='Sbc_cww',m=20.,em=0.02,filter='B_Johnson',
		 cosmo=(0.3,0.7,0.7),vc_filter='B_Johnson'):
	
	"""Generates velocity dispersion and error for a galaxy using TF or Faber--Jackson 
	Inputs: redshift, magnitude, error, filter, spectral type, cosmo and filter for TF in the rest frame
	(FB Jackson always uses BJ as the rest frame filter)
	Usage:
	cosa=vc(0.55,'El_cww',20.,0.02,'I_Cousins',(0.3,0.7,0.7))
	cosa=vc(0.55,'Scd_cww',20.,0.02,'I_Cousins',(0.3,0.7,0.7),'H_Johnson')
	Uses the closest rest frame filter by default
	Assumes that the input magnitudes are AB
	"""
	self.sed=sed

	#Everything has to be properly transformed from AB to Vega!!
	 
	#Info about TF
	#Pierce and Tully 1992
	#vc=158.1*10.**(-(mabs+constant_TF)/slope_TF)
	#It actually only works for Sbc galaxies 
	#For bluer stuff it is better to use the reddest filter, I_Cousins
	#H_Johnson gives weird results, it may be due to template problems
	
	self.filters_TF=['B_Johnson','R_Cousins','I_Cousins','H_Johnson']
	self.centers_TF=[4477.8,6648.33,8086.4,16509.64]
	self.slope_TF=[7.48,8.23,8.72,9.50]
	self.constant_TF=[19.55,20.46,20.94,21.67]
	self.error_TF=[0.14,0.10,0.10,0.08]

	#Info about FB	
	# Kochanek 1994, ApJ
	# sigma_*=(225+-22.5)*(L/L_*)**(.24+-0.03)
	# the error is approximate Kochanek 1996, magnitude is BJ
	# with M_B (BJ) = -19.9+5*log10(h)
	# Using L/L*=10.**[-0.4[M-M_*]] 
	# sigma_*=225.*10.(-(mabs+19.9)/10.42) 

	if sed=='El_cww':
	    self.m_abs=ABtoVega(m_abs(m,z,sed,filter,cosmo,'BJ'),'BJ')
	    self.v_c=225.*10.**(.4*(-self.m_abs-(19.9-5.*log10(cosmo[2])))*0.24)
	    self.e_v_c=(25./225.)*self.v_c
	    self.filter_v_c='BJ'
	else:
	    if sed=='Sbc_cww' or sed=='Scd_cww': 
		fc=filter_center(filter)
	        #Look for the closest filter
		k=argmin(abs(array(self.centers_TF)-fc/(1.+z)))
	    elif sed=='Im_cww' or sed=='SB2_kin' or sed=='SB3_kin': 
	        k=2 #Use I_Cousins		
	    self.m_abs=ABtoVega(m_abs(m,z,sed,filter,cosmo,self.filters_TF[k]),self.filters_TF[k])
	    self.v_c=158.1*10.**(-(self.m_abs+self.constant_TF[k])/self.slope_TF[k])
	    self.e_v_c=self.v_c*2.3/self.slope_TF[k]*sqrt(em**2+self.error_TF[k]**2)
	    self.filter_v_c=self.filters_TF[k]
	    
def dist_mod(z,cosmology=cosmo):
    """Usage: dist_mod(z,cosmology)"""
    return 25.+ 5.*log10(dl(z,cosmology))    
       
def angular_size(length,z,cosmology=cosmo):
    """Usage: angular_size(length,z,cosmology)
       Input: Mpc, Output: arcsec"""
    return length/da(z,cosmology)/3.141592654*180.*3600.

def physical_size(angle,z,cosmology=cosmo):
    """Usage: physical_size(angle,z,cosmology)
       Units: arcseconds, Mpc"""
    return angle/360./60./60.*2.*3.141592654*da(z,cosmology)

#def lookback_time_open(z,omega=0.3,h=0.7):
#    """Usage: lookback_time_open(z,omega,h)
#       Units: Myr Approximation from Peacock"""
#    omega_z=omega*(1.+z)/(1.+omega*z)
#    h_z=h*(1.+z)*sqrt(1.+omega*z)
#    t=h_z*(1.+omega_z**.6/2.)
#    return ht/t/1e9

def obs_E_star(z,filter,cosmology=cosmo,sed="LRG.sed"):
    """Calculate observed flux of a L* red galaxy at redshift z, through a given filter"""
    h=cosmo[2]
    zstar=array((0.1,0.3,0.5,0.7,0.9,1.1))#first point from 2dfGRS, 2nd from DEEP2  (Brown et al. 2007)             
    mstar=array((-19.43,-19.78,-19.92,-20.12,-20.26,-20.67))+5.*log10(h)
    mst=match_resol(zstar,mstar,z)
    return mst-m_abs(0.,z,sed="LRG",filter_m=filter,cosmology=cosmology,filter_M="B_Johnson")


def test():
    test='reobs'
    Testing(test)
    z1,z2,f1,f2,t,c=0.2,0.8,'V_LRIS','I_LRIS','El_cww',(0.3,0.7,0.7)
    dr=reobs(t,0.,z1,f1,z2,f2,c)
    ds=(reobs(t,0.,0.,f1,0.,f2,c)+5.*log10(dl(z2,c)/dl(z1,c))+(kcor(z2,t,f2)-kcor(z1,t,f1)))
    print 'dr,ds'
    print dr,ds

    ask('More?')

    #The reobs function has been tested indirectly by using
    #bpz to estimate redshifts of objects whose colors were generated
    #by reobs, the agreement is perfect. 
    #The rest of the cosmological functions can be tested 
    #by comparing them with plots in the original references
    pass

    test='Distance modulus'
    Testing(test)
    print 'Compare with Peebles Physical Cosmology, page 329'
    print 'Values of Omega are 0.2,0.5,1.'
    z=arange(0.0001,10.,.01)
    omega=[0.2,0.5,1.]
    d=[]
    dlambda=[]
    p1=FramedPlot()
    p1.title='Lambda = 0'
    p1.xrange=-0.2,10.
    p1.yrange=41.,52.
    p2=FramedPlot()
    p2.title='Flat universes'
    p2.xrange=-0.2,10.
    p2.yrange=41.,52.
    for i in range(len(omega)):
	d.append(dist_mod(z,(omega[i],0.,1.)))
	dlambda.append(dist_mod(z,(omega[i],1.-omega[i],1.)))
	p1.add(Curve(z,d[i]))
	p2.add(Curve(z,dlambda[i]))
    p1.show()
    p2.show()

    print 
    print 
    print 

    ask('More?')

    test='Cosmological distances'
    Testing(test)
    z=arange(0.,4.,.01)
    da1=da(z,(1.,0.,1.))/cho
    da2=da(z,(.3,0.,1.))/cho
    da3=da(z,(.3,0.7,1.))/cho
    p=FramedPlot()
    p.add(Curve(z,da1))
    p.add(Curve(z,da3))
    p.add(Curve(z,da2,style='dashed'))
    p.yrange=0.,1.
    p.show()
    print "Compare with Cosmological Physics, page 93"


    print 
    print 
    print 

    ask('More?')

    test='K-corrections'
    Testing(test)
    print 'Compare with Physical cosmology, page 331'
    
    z=arange(0.,1.5,.01)
    p=FramedPlot()
    p.xrange=0.,1.5
    p.yrange=-1.,5.
    for tipo in ['El_cww','Sbc_cww','Scd_cww']:
	k=kcor(z,tipo,'B_Johnson')
	p.add(Curve(z,k))
    p.show()

    print 
    print 
    print 


if __name__ == '__main__':
    test()
else:
    pass



