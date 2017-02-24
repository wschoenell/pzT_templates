"""
bpz: Bayesian Photo-Z estimation
Reference: Benitez 2000, ApJ, 536, p.571
z=Usage:
python bpz.py catalog.cat 
Needs a catalog.columns file which describes the contents of catalog.cat
"""

#import pylab as P
import useful as U
import bpz_tools as B
import cosmology as C
import numpy as N
import os,string,time,sys
# import pylab as P
import re
from scipy import interpolate
import numexpr as ne

# Initialization and definitions

tiny=B.tiny

#if pytables installed use hdf5 files
try:
    import tables as tb
    tb_yes=1.
except:
    tb_yes=0.

tb_yes=1.    

#If pylab installed show plots

try: 
    import pylab as P
    plots=1
except: 
    plots=0

plots=0
#Current directory
homedir=os.getcwd()

#Parameter definition 
pars=U.params()

pars.d={
    "ABSOLUTE_MAGNITUDE":"no",  # Set to yes if you want to calculate absolute magnitudes
    "ABSOLUTE_MAGNITUDE_FILTER": "B_JOHNSON", # Filter in which absolute magnitudes will be calculated
    "ABSOLUTE_MAGNITUDE_CAL": "Vega", # Filter in which absolute magnitudes will be calculated
    "ADD_CONTINUOUS_PROB":None,  #This can be used to add external priors or likelihoods to the estimation
    "ADD_SPEC_PROB":None,        
    "AB_DIR": B.ab_dir,
    "CACHE": "no",          # Use caches
    "CHECK": 'no',          # Perform some checks, compare observed colors with templates, etc.
    "COLOR":'no',           # Use colors instead of fluxes
    "CONVOLVE_P":'no',      # Smooth p(z)
    "DELTA_M0": 0.,         # Add Dm to magnitude (useful for lensed objects)
    "DZ":      0.001,       # redshift resolution
    "DZ_AB": None,          # redshift resolution for the AB files (dont' touch unless you know w.y.a.d.)
    "DZ_FIRST": 0.05,       # Resolution at which we run the first pass
    #    "DZ_LIN": "yes",         # Use linear sampling for the redshift
    'EXCLUDE': 'none',      # Filters to be excluded from the estimation, use "*IRAC*" to exclude all files with that name
    "FASTER": "no",         # Very coarse resolution, runs faster but with more mistakes
    "FC":None,              # Fractions of galaxies in clusters in the field
    "FILTER_0": "HST_ACS_WFC_F814W", #Default filter for M_0, tries to update it with the .colums information (to be used for m_abs and stellar mass calculations)
    'FILTER_DIR': B.fil_dir,
    'GET_Z': 'yes',         # Actually obtain photo-z
    "HDF5": 'no',          #Save probabilities and other data in a HDF5 file 
    'INTERACTIVE':'no',     # Don't query the user
    'INTERP': 3,            # Number of interpolated templates between each of the original ones
    'LINES': None,          # Process only objects belonging from lines l1 to l2 in the python array sense e.g. 0-100 will provide the first 100 number
    'MADAU':'yes',          # Apply Madau correction to spectra
    'MAG':     'yes',       # Data in magnitudes?
    "M0_MIN": None,         # Apply a magnitude cat in the input data
    "M0_MAX": None,         # Apply a mangitude cat in the input data  
    'MERGE_PEAKS':'no',
    'MIN_DM': 0.01,         # Minimum photometric error; mimics intrinsic color variations in galaxies 
    'MP': "no",             # Use average as redshift estimate
    #"MEDP": "yes",         # Use median as redshift estimate (we take this as default)
    'NEW_AB': 'no',         # If yes, generate new AB files even if they already exist
    'N_PEAKS':1,            # If it is desired to obtain information about multiple peaks
    'NOISE_FACTOR': 1.,     # To expand the photometric input errors
    'ODDS':0.6827,            # Odds thresold; affects confidence limits definition
    'ONLY_TYPE':'no',       # Use spectroscopic redshifts instead of photo-z; basically performs spectral classification
    "OUTPUT":None,          # Output .bpz file. The default is the same root as the .cat file
    'PLOTS':'no',           # Don't produce plots 
    'PRIOR':   'SM', # Prior name
    'PROBS_LITE': 'no',     # Save only the final probability distribution
    'P_CAT':1e-5,           # Probability of catastrophic error
    'P_MIN':1e-8,           # Set to 0 all probabilities smaller than this x pmax when saving data
    'SED_DIR': B.sed_dir,
    'SIGMA_CONVOLVE':0.01,  # Smooth p(z) convolving with a gaussian with this width
    'SIGMA_EXPECTED': 0.03, # This is the expected DZ precision; it is used to calculate the odds (integrates within +-1 sigma)
    'SLOW':"no",            # Run more slowly, explore full resolution 
    'SPECTRA':'eB11.list',  # Template list 
    "STELLAR_MASS": "no",   # Calculate stellar masses
    "USE_Z_S": "yes",       # If no, it will disregard information about spectroscopic redshifts in the .columns file
    'VERBOSE': 'yes',       # Print estimated redshifts to the standard output
    'ZC': None,             #To add a cluster spike-like prior
    'ZMAX':    7.0,         # maximum redshift
    'ZMIN':    1e-4,        # minimum redshift
    'ZP_ERRORS':0.,
    'ZP_OFFSETS':0.,
    'Z_THR':0,              #Integrate probability for z>z_thr
}               

# Define the default values of the parameters 
pars.d['INPUT']=sys.argv[1]       # catalog with the photometry
obs_file=pars.d['INPUT']          #
root=os.path.splitext(pars.d['INPUT'])[0]
pars.d['COLUMNS']=root+'.columns' # column information for the input catalog
if not pars.d["OUTPUT"]:
    pars.d['OUTPUT']= root+'.bpz'     # output 

nargs=len(sys.argv)
ipar=2

if nargs>2: #Check for parameter file and update parameters
    if  sys.argv[2]=='-P': 
	pars.fromfile(sys.argv[3])
	ipar=4
# Update the parameters using command line additions
pars.fromcommandline(sys.argv[ipar:]) 

rolex=U.watch()


if pars.d["HDF5"].lower()=="no":
    hdf5=0
else:
    hdf5=1
    if pars.d["HDF5"].lower()=="yes":
        pars.d["HDF5"]=pars.d["OUTPUT"].split(".")[0]+".hdf5"

if hdf5:
    print ""
    print "Saving probabilities and other data to %s" % pars.d["HDF5"]
    print ""
    
#This allows to change the auxiliary directories used by BPZ
if pars.d['SED_DIR']<>B.sed_dir:
    print "Changing sed_dir to ",pars.d['SED_DIR'] 
    sed_dir=pars.d['SED_DIR']
    if sed_dir[-1]<>'/': B.sed_dir+='/'
else:
    sed_dir=B.sed_dir

if pars.d['AB_DIR']<>B.ab_dir:
    print "Changing ab_dir to ",pars.d['AB_DIR'] 
    ab_dir=pars.d['AB_DIR']
    if ab_dir[-1]<>'/': ab_dir+='/'
else:
    ab_dir=B.ab_dir

if pars.d['FILTER_DIR']<>B.fil_dir:
    print "Changing fil_dir to ",pars.d['FILTER_DIR'] 
    fil_dir=pars.d['FILTER_DIR']
    if fil_dir[-1]<>'/': fil_dir+='/'
else:
    fil_dir=B.fil_dir

if pars.d["DZ_AB"]<>None:
    dz_ab=float(pars.d["DZ_AB"])
else:
    dz_ab=B.dz_ab

#Better safe than sorry
if pars.d['OUTPUT']==obs_file or pars.d['HDF5']==obs_file or pars.d['PROBS_LITE']==obs_file:
    print "This would delete the input file!"
    sys.exit()
if pars.d['OUTPUT']==pars.d['COLUMNS'] or pars.d['PROBS_LITE']==pars.d['COLUMNS'] or pars.d['HDF5']==pars.d['COLUMNS']:
    print "This would delete the .columns file!"
    sys.exit()    

dmin_dm=float(pars.d['MIN_DM'])
dz_exp=float(pars.d["SIGMA_EXPECTED"])
#Set floor on dmin_dm to avoid overflows
if dmin_dm<1e-4: dmin_dm=1e-4

stmass=0

if pars.d["STELLAR_MASS"][0].lower()=="y":
    stmass=1
    
absmag=0
if pars.d["ABSOLUTE_MAGNITUDE"][0].lower()=="y":
    absmag=1

if pars.d['PLOTS'][0].lower()=='n': plots=0
if pars.d['INTERACTIVE'][0].lower()=='n': inter=0
else: inter=1
if pars.d['VERBOSE'][0].lower()=='y': 
    print "Current parameters"
    U.view_keys(pars.d)
pars.d['N_PEAKS']=int(pars.d['N_PEAKS'])
if pars.d["ADD_SPEC_PROB"]<>None:
    specprob=1
    specfile=pars.d["ADD_SPEC_PROB"]
    spec=U.get_2Darray(specfile)
    ns=spec.shape[1]
    if ns/2<>(ns/2.):
        print "Number of columns in SPEC_PROB is odd"
        sys.exit()
    z_spec=spec[:,:ns/2]
    p_spec=spec[:,ns/2:]
    # Write output file header
    header="#ID "
    header+=ns/2*" z_spec%i"
    header+=ns/2*" p_spec%i"
    header+="\n"
    header=header % tuple(list(range(ns/2))+list(range(ns/2)))
    specout=open(specfile.split()[0]+".p_spec","w")
    specout.write(header)
else:
    specprob=0
pars.d['DELTA_M0']=float(pars.d['DELTA_M0'])    

#Some misc. initialization info useful for the .columns file
nofilters=['M_0','OTHER','ID','Z_S','X','Y']

#Numerical codes for nondetection, etc. in the photometric catalog
unobs=-99. #Objects not observed
undet= 99.  #Objects not detected


#Define the z-grid
zmin=float(pars.d['ZMIN'])
zmax=float(pars.d['ZMAX'])
if zmin > zmax : raise 'zmin < zmax !'
dz=float(pars.d['DZ'])

#if pars.d["DZ_LIN"][0].lower()=="y":
#    linear=1
#else:
#    linear=0

#if linear:
#   z=N.arange(zmin,zmax+dz,dz)
#else:
#    x1=N.log10(1.+zmin)
#    x2=N.log10(1.+zmax)
#    dx=N.log10(1.+zmin+dz)-x1
#    x=N.arange(x1,x2,dx)
#    z=10.**x-1.

z=N.arange(zmin,zmax+dz,dz)

dzf=float(pars.d["DZ_FIRST"])
nsz=int(dzf/dz)
zf=z[::nsz]

#Get a list with the filter names and check whether they are in stock
col_file=pars.d['COLUMNS']
filters=U.get_str(col_file,0)

for cosa in nofilters: 
    if filters.count(cosa):filters.remove(cosa)
    
for i in range(len(filters)):
    if filters[i][-4:]=='.res': filters[i]=filters[i][:-4]
    
    if not os.path.exists(fil_dir+filters[i]+".res"):
	print 'filter ', filters[i], 'not in filter directory',fil_dir, ':'
        sys.exit()

if pars.d['EXCLUDE']<>'none':
    if type(pars.d['EXCLUDE'])==type(' '):
	pars.d['EXCLUDE']=[pars.d['EXCLUDE']]
    else:
        pars.d['EXCLUDE']=list(pars.d['EXCLUDE'])

    for cosa in pars.d['EXCLUDE']:
        if "*" in cosa:
            for j in N.arange(len(filters)):
                if re.match(cosa.replace("*",".*"),filters[j]):
                    pars.d["EXCLUDE"].append(filters[j])

    for filtro in pars.d['EXCLUDE']:
        if filters.count(filtro): filters.remove(filtro)
    print "Using filter system: ", filters
    #

# Get a list with the spectrum names and check whether they're in stock
# Look for the list in the home directory first, 
# if it's not there, look in the SED directory
spectra_file=os.path.join(homedir,pars.d['SPECTRA'])
if not os.path.exists(spectra_file):
    spectra_file=os.path.join(sed_dir,pars.d['SPECTRA'])

spectra=B.get_lib(spectra_file)
for i in range(len(spectra)):
    if spectra[i][-4:]=='.sed': spectra[i]=spectra[i][:-4]

nf=len(filters)
nt=len(spectra)
nz=len(z)

#Get the model fluxes
f_mod=N.zeros((nz,nt,nf),'float')
abfiles=[]

#Parameters to define interpolation between the colors
nt0=nt
n_interp=int(pars.d['INTERP'])
xt= N.linspace(0,nt0-1,(nt0-1)*(n_interp+1)+1)
nt=len(xt)

cached_fmod=os.path.splitext(pars.d["COLUMNS"])[0]+"_"+pars.d["SPECTRA"][:-5]+"_"+str(n_interp)+"_"+str(pars.d["DZ"])+"_"+str(pars.d["ZMAX"])+".npy"
readfmod=0
if os.path.isfile(cached_fmod):
    readfmod=  pars.d["NEW_AB"][0].lower()=="n"
    readfmod*= os.path.getmtime(cached_fmod) > os.path.getmtime(spectra_file)
    readfmod*= os.path.getmtime(cached_fmod) > os.path.getmtime(pars.d["COLUMNS"])

if readfmod and pars.d["CACHE"][0].lower()=="y":
    f_mod=N.load(cached_fmod)
    print "Reading f_mod from cache",cached_fmod
else:
    B.ABflux_M(spectra,filters)
    for it in range(nt0):
        for jf in range(nf):
            if filters[jf][-4:]=='.res': filtro=filters[jf][:-4]
            else: filtro=filters[jf]
            model=".".join([spectra[it],filtro,'AB'])
            model_path=os.path.join(ab_dir,model)
            abfiles.append(model)
            # Generate new ABflux files if not present
            # or if new_ab flag on
            newab=0.
            if pars.d["NEW_AB"][0].lower()=="y": newab=1
            else:
                # This catches corrupt AB files
                try:
                    # zo,f_mod_0=U.get_data(model_path,(0,1))
                    zo,f_mod_0=N.loadtxt(model_path,unpack=True)
                    f_mod[:,it,jf]=U.match_resol(zo,f_mod_0,z)
                except:
                    newab=1                
            if newab:
                if not os.path.exists(sed_dir+spectra[it]+".sed"):
                    print 'SED ', spectra[it], 'not in ',sed_dir
                    sys.exit()
                print '     Generating ',model,'....'
                B.ABflux(spectra[it],filtro,dz_ab=dz_ab,sed_dir=sed_dir)
                zo,f_mod_0=U.get_data(model_path,(0,1))
                #zo,f_mod_0=numpy.loadtxt(model_path,unpack=True)
                f_mod[:,it,jf]=U.match_resol(zo,f_mod_0,z)
                
            if N.less(f_mod[:,it,jf],0.).any():
                print 'Warning: some values of the model AB fluxes are <0'
                print 'due to the interpolation '
                print 'Clipping them to f>=0 values'
                # To avoid rounding errors in the calculation of the likelihood
                f_mod[:,it,jf]=N.clip(f_mod[:,it,jf],tiny,1./tiny)

                # We forbid f_mod to take values in the (0,1e-100) interval
                # f_mod[:,it,jf]=where(N.less(f_mod[:,it,jf],1e-100)*N.greater(f_mod[:,it,jf],0.),0.,f_mod[:,it,jf])
    if n_interp:
        xt0=N.linspace(0,nt0-1,nt0)
        buffer=N.zeros((nz,nt,nf),'float')
        for iz in N.arange(nz):
            for jf in range(nf):
                buffer[iz,:,jf]=U.match_resol(xt0,f_mod[iz,:,jf],xt)
        f_mod=buffer
    N.save(cached_fmod,f_mod)
    print "Writing f_mod to cache",cached_fmod
            
# Load all the parameters in the columns file to a dictionary   
col_pars=U.params(cap=0)
col_pars.fromfile(col_file)
for filt in col_pars.d.keys():
    if "." in filt:
        col_pars.d[filt.split(".")[0]]=col_pars.d[filt]
        del(col_pars.d[filt])
        
# Read which filters are in which columns
flux_cols=[]
eflux_cols=[]
cals=[]
zp_errors=[]
zp_offsets=[]

for filt in filters:
    datos=col_pars.d[filt]
    flux_cols.append(int(datos[0])-1)
    eflux_cols.append(int(datos[1])-1)
    cals.append(datos[2])
    zp_errors.append(datos[3])
    zp_offsets.append(datos[4])

zp_errors=N.array(map(float,zp_errors))
zp_offsets=N.array(map(float,zp_offsets))

if pars.d['ZP_OFFSETS']:
    try:
        zpoff=float(pars.d['ZP_OFFSETS'])
    except:
        zpoff=N.array(map(float,pars.d['ZP_OFFSETS']))
    #Add zp_offsets to the already existing in the .columns file
    print "Adding command line zp offsets to those in the .columns file"
    zp_offsets+=zpoff

    pars.d['ZP_OFFSETS']=",".join(map(lambda x: str(x),zp_offsets))

if pars.d['ZP_ERRORS']:
    try:
        zperr=float(pars.d['ZP_ERRORS'])
    except:
        zperr=N.array(map(float,pars.d['ZP_ERRORS']))        
    #Substitute input zp_errors in the .cols file for those in the command line
    print "Substituting offsets in the .columns file for those in the command line"
    zp_errors=zp_offsets*0.+zperr # Make sure both have same dimension

#########
flux_cols=tuple(flux_cols)
ncols=len(flux_cols)
iflux=N.arange(ncols)
cols=list(flux_cols)
eflux_cols=tuple(eflux_cols)
ieflux=N.arange(ncols,2*ncols)
ncols=2*ncols
cols=cols+list(eflux_cols)
if col_pars.d.has_key('M_0'):
    m_0_col=int(col_pars.d['M_0'])-1
    cols=cols+[m_0_col]
    ncols=len(cols)
    im0=ncols-1
spec=pars.d["USE_Z_S"][0].lower()=="y" and col_pars.d.has_key("Z_S")
if spec:
    z_s_col=int(col_pars.d['Z_S'])-1
    cols=cols+[z_s_col]
    ncols=len(cols)
    izs=ncols-1
if col_pars.d.has_key('X'):
    x_col=int(col_pars.d['X'])-1
    cols=cols+[x_col]
    ncols=len(cols)
    ix=ncols-1
if col_pars.d.has_key('Y'):
    y_col=int(col_pars.d['Y'])-1
    cols=cols+[y_col]
    ncols=len(cols)
    iy=ncols-1

cols=tuple(cols)

#READ the flux and errors from obs_file
cached_file=obs_file[:-3]+"npy"
if os.path.isfile(cached_file) and os.path.getmtime(cached_file)> os.path.getmtime(obs_file) and pars.d["CACHE"][0].lower()=="y":
    obs=N.load(cached_file)
    print "Reading photometric data from cache"
else:
    obs=N.loadtxt(obs_file,usecols=cols)
    N.save(cached_file,obs)
    print "Writing data to cache"

ng0=obs.shape[0]
if pars.d["LINES"]:
    try:
        l1,l2=map(int,pars.d["LINES"].split("-"))
    except:
        l1,l2=0,int(pars.d["LINES"])
else:
    l1,l2=0,ng0

gg=N.greater_equal(range(ng0),l1)*N.less(range(ng0),l2)
if pars.d["M0_MIN"]:
    gg*=N.greater_equal(obs[:,im0],float(pars.d["M0_MIN"]))
if pars.d["M0_MAX"]:
    gg*=N.less(obs[:,im0],float(pars.d["M0_MAX"]))
obs=N.compress(gg,obs,0)
    
f_obs=obs[:,iflux]
ef_obs=obs[:,ieflux]*float(pars.d["NOISE_FACTOR"])
ng=f_obs.shape[0]

if col_pars.d.has_key('M_0'):
    try:
        pars.d["FILTER_0"]=filters[list(flux_cols).index(m_0_col)]
        print "FILTER_0 set to", pars.d["FILTER_0"]
    except:
        print "Could not find value of FILTER_0"
        print "FILTER_0 set to default", pars.d["FILTER_0"]
        print "Calculations of stellar mass and absolute magnitudes may be erroneous"
    #m_0=obs[:,(im0,)]
    m_0=obs[:,im0]
    if pars.d["MAG"]=="no":
        print "Taking log10 of flux to calculate magnitude for prior"
        m_0=B.flux2mag(m_0)
    if pars.d["DELTA_M0"]:
        print "Adding DELTA_M0",pars.d['DELTA_M0']
        m_0+=pars.d['DELTA_M0']
    m_0=N.clip(m_0,0.,99.)

if stmass:
    # Now calculate the tabulated stellar masses for all combinations of types and redshifts
    sm=C.sm_grid(pars.d["SPECTRA"],pars.d["FILTER_0"],zmax)
    #sm20=interpolate.RectBivariateSpline(sm.z,sm.t,sm.sm20,kx=1,ky=1,s=0)(z,xt)
    sm20=interpolate.RectBivariateSpline(sm.z,sm.t-1.,sm.sm20,kx=1,ky=1,s=0)(z,xt)

if absmag:
    cached_mabs="MABS_"+pars.d["SPECTRA"][:-5]+"_"+str(n_interp)+"_"+str(pars.d["DZ"])+"_"+str(pars.d["ZMAX"])+"_"+str(pars.d["FILTER_0"])+".npy"
    # Now calculate the tabulated stellar masses for all combinations of types and redshifts
    if not os.path.isfile(cached_mabs):
        ixz,ixt=N.mgrid[:nz,:nt]
        ixz=N.ravel(ixz)
        ixt=N.ravel(ixt)
        
        ## Original line from Tx. commented out by Alberto
        ## mabs20=C.mzt_m_abs(20.+z[ixz]*0.,z[ixz],xt[ixt]+1.,lib=spectra_file,filter_m=pars.d["FILTER_0"],filter_M=pars.d["ABSOLUTE_MAGNITUDE_FILTER"],cal_m="AB",cal_M=pars.d["ABSOLUTE_MAGNITUDE_CAL"])
        # Using this line or the former one yields different results...!!!
        mabs20=C.mzt_m_abs(20.+z[ixz]*0.,z[ixz],xt[ixt]+1.,'eB11.list','HST_ACS_WFC_F814W','B_JOHNSON','AB','AB')
        
        # mabs20=C.mzt_m_abs(20.+z[ixz]*0.,z[ixz],xt[ixt]+1.,lib=spectra_file,
        #                    filter_m=pars.d["FILTER_0"],filter_M=pars.d["ABSOLUTE_MAGNITUDE_FILTER"],
        #     cal_m="AB",cal_M=pars.d["ABSOLUTE_MAGNITUDE_CAL"])
        mabs20=N.resize(mabs20,(nz,nt))
        N.save(cached_mabs,mabs20)
        print "Writing Absolute Magnitude information to cache"
    else:
        mabs20=N.load(cached_mabs)
        print "Reading Absolute Magnitude information from cache"
        
# Get spectroscopic redshifts (if present and USE_Z_S=yes)
if spec:
    z_s=obs[:,izs]

#Get the X,Y coordinates
if col_pars.d.has_key('X'):
    x=obs[:,ix]

if col_pars.d.has_key('Y'):
    y=obs[:,iy]

#Get the objects ID (as a string)
if col_pars.d.has_key('ID'):
    #    print col_pars.d['ID']
    id_col=int(col_pars.d['ID'])-1
    iden=U.get_str(obs_file,id_col)
    #id=map(str,obs[:,(id_col,)])
else:
    iden=map(str,range(1,len(f_obs[:,0])+1))
    
if N.rank(f_obs)==1:
    print "One line input file!"
    print "BPZ doesn't like such files"
    print "Duplicating the size of the input to avoid problems"
    f_obs=N.resize(f_obs,(2,f_obs.shape[0]))
    ef_obs=N.resize(ef_obs,(2,ef_obs.shape[0]))
    iden=N.resize(iden,2)
    try: y=N.resize(y,2)
    except: pass
    try: x=N.resize(y,2)
    except: pass
    try: z_s=N.resize(z_s,2)
    except: pass
    try: m_0=N.resize(m_0,2)
    except: pass
    ng=f_obs.shape[0]


#Convert them to arbitrary fluxes if they are in magnitudes
if pars.d['MAG']=='yes':
    seen=N.asarray(N.greater(f_obs,0.)*N.less(f_obs,undet),'int')
    no_seen=N.asarray(N.equal(f_obs,undet),'int')
    no_observed=N.asarray(N.equal(f_obs,unobs),'int')
    todo=(seen+no_seen+no_observed)*1.
    #The minimum photometric error is dmin_dm
    #ef_obs=ef_obs+seen*equal(ef_obs,0.)*dmin_dm
    ef_obs=N.where(N.greater_equal(ef_obs,0.),N.clip(ef_obs,dmin_dm,1e10),ef_obs)
    # Check that all inputs are kosher
        
    if not N.alltrue(todo):
	print 'Objects with unexpected magnitudes!'
	print """Allowed values for magnitudes are 
	0<m<"""+`undet`+" m="+`undet`+"(non detection), m="+`unobs`+"(not observed)" 
        formatillo="% i "+ len(f_obs[i,:])*"%.2f "+ len(f_obs[i,:])*"%.2f "

        for j in N.arange(len(filters)):
            if not N.alltrue(todo[:,j]):
                print filters[j]+" contains wrong values"
                if U.ask("Do you want to see them?"):
                    for i in range(len(todo[:,j])):
                        if not todo[i,j]: print i,f_obs[i,j],ef_obs[i,j]
                    print "N= %i wrong values" % N.equal(todo[:,j],0.).sum()

                U.ask()
                        
	for i in range(len(todo)):
	    if not N.alltrue(todo[i,:]):
		print formatillo % tuple([i+1]+list(f_obs[i,:])+list(ef_obs[i,:]))


	sys.exit()
 
    #Detected objects

    try: 
        f_obs=N.where(seen,10.**(-.4*f_obs),f_obs)
    except OverflowError:
        print 'Some of the input magnitudes have values which are >700 or <-700'
        print 'Purge the input photometric catalog'
        print 'Minimum value',f_obs.min()
        print 'Maximum value',f_obs.max()
        print 'Indexes for minimum values',N.argmin(f_obs,0.)
        print 'Indexes for maximum values',N.argmax(f_obs,0.)
        print 'Bye.'
        sys.exit()

    try:
        ef_obs=N.where(seen,(10.**(.4*ef_obs)-1.)*f_obs,ef_obs)
    except OverflowError:
        print 'Some of the input magnitude errors have values which are >700 or <-700'
        print 'Purge the input photometric catalog'
        print 'Minimum value',ef_obs.min()
        print 'Maximum value',ef_obs.max()
        print 'Indexes for minimum values',N.argmin(ef_obs,0.)
        print 'Indexes for maximum values',N.argmax(ef_obs,0.)
        print 'Bye.'
        sys.exit()

    # Looked at, but not detected objects (mag=99.)
    # We take the flux equal to zero, and the error in the flux equal to the 1-sigma detection error.
    # If m=99, the corresponding error magnitude column in supposed to be dm=m_1sigma, to avoid errors
    # with the sign we take the absolute value of dm 
    f_obs=N.where(no_seen,0.,f_obs)
    ef_obs=N.where(no_seen,10.**(-.4*abs(ef_obs)),ef_obs)

    # Objects not looked at (mag=-99.)
    f_obs=N.where(no_observed,0.,f_obs)
    ef_obs=N.where(no_observed,1e10,ef_obs) #Mimic non-observation by using a huge observation error
    # print f_obs
    # print ef_obs

#Flux codes:
# If f>0 and ef>0 : normal objects
# If f==0 and ef>0 :object not detected
# If f==0 and ef==0: object not observed
#Everything else will crash the program

#Check that the observed error fluxes are reasonable 
#if sometrue(N.less(ef_obs,0.)): raise 'Negative input flux errors'
#if sum(sum(N.asarray(N.less(ef_obs,0.),'int'))) > 0:
#    raise 'Negative input flux errors'

if (N.less(ef_obs,0.)*N.greater(ef_obs,0.)).sum():
    print 'Negative input flux errors'
    raise 

f_obs=N.where(N.less(f_obs,0.),0.,f_obs) #Put non-detections to 0
ef_obs=N.where(N.less(f_obs,0.),N.maximum(1e-100,f_obs+ef_obs),ef_obs) # Error equivalent to 1 sigma upper limit

#if sometrue(N.less(f_obs,0.)) : raise 'Negative input fluxes'

seen=N.greater(f_obs,0.)*N.greater(ef_obs,0.)*1.
no_seen=N.equal(f_obs,0.)*N.greater(ef_obs,0.)*1.
no_observed=N.equal(f_obs,0.)*N.equal(ef_obs,0.)*1.

todo=N.asarray(seen+no_seen+no_observed,'int')

#Convert (internally) objects with zero flux and zero error(non observed)
#to objects with almost infinite (~1e108) error and still zero flux
#This will yield reasonable likelihoods (flat ones) for these objects
ef_obs=N.where(no_observed,1e108,ef_obs)

#Include the zero point errors
zp_frac=B.e_mag2frac(zp_errors)
#zp_frac=10.**(.4*zp_errors)-1.
ef_obs=N.where(seen,N.sqrt(ef_obs*ef_obs+(zp_frac*f_obs)**2),ef_obs)
ef_obs=N.where(no_seen,N.sqrt(ef_obs*ef_obs+(zp_frac*(ef_obs/2.))**2),ef_obs)

#Add the zero-points offset
#The offsets are defined as m_new-m_old
zp_offsets=N.array(map(float,zp_offsets))
zp_offsets=N.where(N.not_equal(zp_offsets,0.),10.**(-.4*zp_offsets),1.)
f_obs=f_obs*zp_offsets
ef_obs=ef_obs*zp_offsets

#Convert fluxes to AB if needed
for i in range(f_obs.shape[1]):
    if cals[i]=='Vega':
	const=B.mag2flux(B.VegatoAB(0.,filters[i]))
	f_obs[:,i]=f_obs[:,i]*const
	ef_obs[:,i]=ef_obs[:,i]*const
    elif cals[i]=='AB':continue
    else:
	print 'AB or Vega?. Check '+col_file+' file'
	sys.exit()
		
#If 'check' on, initialize some variables
check=pars.d['CHECK']

if check=='yes':
    r=N.zeros((ng,nf),'float')+1.
    dm=N.zeros((ng,nf),'float')+1.
    w=r*0.
    wg=r*0.

#Visualize the colors of the galaxies and the templates 

#When there are spectroscopic redshifts available, you can plot the color/redshift comparison for all templates
if inter and spec and plots and check=='yes' and not U.ask('Pass without plotting colors vs spectroscopic redshifts?'):
    color_m=N.zeros((nz,nt,nf-1),'float')
    P.figure(1)
    nrows=2
    ncols=(nf-1)/nrows
    if (nf-1)%nrows: ncols+=1
    for i in range(nf-1):
	# Check for overflows
	fmu=f_obs[:,i+1]
	fml=f_obs[:,i]
	good=N.greater(fml,1e-100)*N.greater(fmu,1e-100)
	zz,fmu,fml=U.multicompress(good,(z_s,fmu,fml))
	colour=fmu/fml
	colour=N.clip(colour,1e-5,1e5)
	colour=2.5*N.log10(colour)
        P.subplot(nrows,ncols,i+1)
        P.plot(zz,colour,"bo")
	for it in range(nt):
	    # Prevent overflows
	    fmu=f_mod[:,it,i+1]
	    fml=f_mod[:,it,i]
	    good=N.greater(fml,1e-100)
	    zz,fmu,fml=U.multicompress(good,(z,fmu,fml))
	    colour=fmu/fml
	    colour=N.clip(colour,1e-5,1e5)
	    colour=2.5*N.log10(colour)
            P.plot(zz,colour,"r")
        P.xlabel(r'$z$')
        P.ylabel('%s - %s' %(filters[i],filters[i+1]))

    P.show()

# Get other information which will go in the output file (as strings)
if col_pars.d.has_key('OTHER'):
    if col_pars.d['OTHER']<>'all':
	other_cols=col_pars.d['OTHER']
        if type(other_cols)==type((2,)):
            other_cols=tuple(map(int,other_cols))
        else:
            other_cols=(int(other_cols),)
	other_cols=map(lambda x: x-1,other_cols)
	n_other=len(other_cols)
    else:
	n_other=get_2Darray(obs_file,cols='all',nrows=1).shape[1]
	other_cols=range(n_other)

    others=U.get_str(obs_file,other_cols)

    if len(other_cols)>1:
	other=[]
	for j in range(len(others[0])):
	    lista=[]
	    for i in range(len(others)):
		lista.append(others[i][j])
	    other.append(" ".join(lista))
    else:
	other=others

if pars.d['GET_Z']=='no': get_z=0
else: get_z=1

# Prepare the output file
out_name=pars.d['OUTPUT']
resfile_name=os.path.splitext(out_name)[0]+".txt"
resfile=open(resfile_name,"w")

if get_z:
    if os.path.exists(out_name):
        if pars.d["VERBOSE"]=="yes": print "File %s exists. Copying it to %s.bak" % (out_name,out_name)
        os.system('cp %s %s.bak' % (out_name,out_name))
    output=open(out_name,'w')

# This generates a file with m,z,T and observed/expected colors
if check=='yes':
    pars.d['FLUX_COMPARISON']=os.path.splitext(pars.d["OUTPUT"])[0]+'.flux_comparison'
    
if pars.d['PROBS_LITE']=='no': 
    save_probs=0
else: 
    save_probs=1

# Include some header information

#   File name and the date...
time_stamp=time.ctime(time.time())
if get_z: output.write('## File '+out_name+'  '+time_stamp+'\n')

# and also the parameters used to run bpz...
if get_z:output.write("""##
##Parameters used to run BPZ:
##
""")
claves=pars.d.keys()
claves.sort()
for key in claves:
    if type(pars.d[key])==type((1,)):
	cosa=",".join(list(pars.d[key]))
    else:
	cosa=str(pars.d[key])
    if get_z: output.write('##'+key.upper()+'='+cosa+'\n')

if hdf5:
    filtros=tb.Filters(complevel=5,complib="lzo") #lz0 is much faster than zlib
    fp_file=tb.openFile(pars.d["HDF5"],mode="w",title="BPZ outputs")
    zh=             fp_file.createArray(fp_file.root, "redshift",z)
    th=             fp_file.createArray(fp_file.root, "type",   xt)
    m0h=             fp_file.createArray(fp_file.root, "m_0",   m_0)
    full_probs=     fp_file.createCArray(fp_file.root,"FullProbability",
                                         tb.Float32Atom(),shape=(ng,nz,nt),
                                         chunkshape=(1,nz,nt),
                                         filters=filtros)
    full_likelihood=fp_file.createCArray(fp_file.root,"Likelihood",
                                         tb.Float32Atom(),
                                         shape=(ng,nz,nt),
                                         chunkshape=(1,nz,nt),filters=filtros)
    if stmass:
        sm20h=             fp_file.createArray(fp_file.root, "Stellar_Mass_zT_for_m0eq20",   sm20)
    if absmag:
        mabs20h=           fp_file.createArray(fp_file.root, "Absolute_Magnitude_zT_for_m0eq20",   mabs20)
        print "Creating %s" % pars.d["HDF5"]

if save_probs:
    probs=open(pars.d['PROBS_LITE'],'w')
    probs.write('# ID  p_bayes(z)  where z=N.arange(%.4f,%.4f,%.4f) \n' % (zmin,zmax+dz,dz))
    
# Use a empirical prior?
tipo_prior=pars.d['PRIOR']
useprior=0
if col_pars.d.has_key('M_0'): has_mags=1
else: has_mags=0
if has_mags and tipo_prior<>'none' and tipo_prior<>'flat': useprior=1

# Add cluster 'spikes' to the prior?
cluster_prior=0.
if pars.d['ZC'] : 
    cluster_prior=1
    if type(pars.d['ZC'])==type(""): zc=N.array([float(pars.d['ZC'])])
    else:    zc=N.array(map(float,pars.d['ZC']))
    if type(pars.d['FC'])==type(""): fc=N.array([float(pars.d['FC'])])
    else:    fc=N.array(map(float,pars.d['FC']))    

    fcc=N.add.reduce(fc)
    if fcc>1. : 
	print ftc
	raise 'Too many galaxies in clusters!'
    pi_c=N.zeros((nz,nt),'float')    
    #Go over the different cluster spikes
    for i in range(len(zc)):
	#We define the cluster within dz=0.015 limits
	cluster_range=N.less_equal(N.abs(z-zc[i]),.01)*1.
	#Clip values to avoid overflow
	exponente=N.clip(-(z-zc[i])**2/2./(0.003)**2,-700.,0.)
	#Outside the cluster range g is 0
	g=N.exp(exponente)*cluster_range
	norm=g.sum()
	pi_c[:,0]=pi_c[:,0]+g/norm*fc[i]

    #Go over the different types
    if pars.d['SPECTRA']=="eB11.list":
        print 'We only apply the cluster prior to the LRGs'
        for i in range(1,5+2*n_interp):
            pi_c[:,i]=pi_c[:,i]+pi_c[:,0]
    else:
        for i in range(1,2+2*n_interp):
            pi_c[:,i]=pi_c[:,i]+pi_c[:,0]
        

#Output format
format='%'+`N.maximum(5,len(iden[0]))`+'s' #ID format
peakformat=' %.4f %.3f %.3f %5.2f %10.6f'
if stmass: peakformat+=" %.4f"
if absmag: peakformat+=" %.4f"
format=format+pars.d['N_PEAKS']*peakformat+' %6.3f %5.2f %10.4f'

#Add header with variable names to the output file
sxhdr="""##
##Column information
##
# 1 ID"""
k=1

if pars.d['N_PEAKS']>1:
    for j in range(pars.d['N_PEAKS']):
        sxhdr+="""
# %i Z_B_%i
# %i Z_B_MIN_%i
# %i Z_B_MAX_%i
# %i T_B_%i
# %i ODDS_%i""" % (k+1,j+1,k+2,j+1,k+3,j+1,k+4,j+1,k+5,j+1)
        k+=5
        if stmass:
            k+=1
            sxhdr=sxhdr+'#\n %i STELLAR MASS_%i(log10(M_sun))' % (k+1,j+1)
        if absmag:
            k+=1
            sxhdr=sxhdr+'#\n %i M_ABS_%i' % (k+1,j+1)
else:
    sxhdr+="""
# %i Z_B
# %i Z_B_MIN
# %i Z_B_MAX
# %i T_B
# %i ODDS""" % (k+1,k+2,k+3,k+4,k+5)
    k+=5
    if stmass:
        k+=1
        sxhdr=sxhdr+'\n# %i STELLAR MASS(log10(M_sun))' % k
    if absmag:
        k+=1
        sxhdr=sxhdr+'\n# %i M_ABS' % k
    
sxhdr+="""    
# %i Z_ML
# %i T_ML
# %i CHI-SQUARED\n""" % (k+1,k+2,k+3)

nh=k+4
if spec:
    sxhdr=sxhdr+'# %i Z_S\n' % nh
    format=format+'  %.5f'
    nh+=1
if has_mags: 
    format=format+'  %.3f'
    sxhdr=sxhdr+'# %i M_0\n' % nh
    nh+=1
if col_pars.d.has_key('OTHER'):
    sxhdr=sxhdr+'# %i OTHER\n' % nh
    format=format+' %s'
    nh+=n_other

if get_z: output.write(sxhdr+'##\n')

#Proceed to redshift estimation
if check=='yes': buffer_flux_comparison=""

if pars.d['CONVOLVE_P']=='yes':
    # Convolve with gaussian to make probabilities smoother
    # This is recommended; if not there may be too many close peaks
    sigma_g=float(pars.d["SIGMA_CONVOLVE"])
    x=N.arange(-3.*sigma_g,3.*sigma_g+dz/100.,dz)
    gaus=N.exp(-(x/sigma_g)**2)

# Catastrophic error probability (it is an effective 1-2%)
if pars.d['P_CAT']:
    pC=N.ones((nz,nt),"float")*float(pars.d['P_CAT'])/nt/(zmax-zmin)
    normpC=pC.sum()

if pars.d["VERBOSE"][0].lower()=="y":
    print "Time spent in initialization"
    tinit=rolex.t()
    print rolex.check()

eff=0.
tm0=0.
tm1=0.
for ig in range(ng):    
    if not get_z: continue
    if useprior:
        if pars.d['PRIOR']=='lensing':
            p_i=B.prior(z,m_0[ig],tipo_prior,nt0,n_interp,x[ig],y[ig])
        elif pars.d['PRIOR'][:2]=="SM":
            #Use the rest of the letter to signal whether we use alpha and phi
            p_i=B.prior(z,m_0[ig],tipo_prior,nt0,n_interp,filter_m=pars.d["FILTER_0"],lib=spectra_file)
        else:
            p_i=B.prior(z,m_0[ig],tipo_prior,nt0,n_interp)
    else:
        p_i=N.ones((nz,nt))/float(nz*nt)

    #print m_0[ig],(p_i[:,:]*z[:,N.newaxis]).mean()

    if cluster_prior:p_i=(1.-fcc)*p_i+pi_c

    if pars.d['COLOR']=='yes':
        likelihood=p_c_z_t_color(f_obs[ig,:nf],ef_obs[ig,:nf],f_mod[:nz,:nt,:nf])
    else:
        f=f_obs[ig,:nf]
        ef=ef_obs[ig,:nf]
        fef=f/ef
        fefef=fef/ef
        foo=N.add.reduce((fef)**2)
        #ft_z=f_mod[:,:nt,:nf]
        if pars.d["SLOW"][0].lower()=="y":
            likelihood=B.p_c_z_t(f,ef,f_mod)
            iz_ml=likelihood.i_z_ml
            t_ml=likelihood.i_t_ml
            red_chi2=likelihood.min_chi2/float(nf-1.)
            norm=likelihood.likelihood.max()
            p=likelihood.likelihood/norm*N.exp(-0.5*red_chi2) #Normalize peak probability to reduced chi2
        else:
            #Quick sample of the probability (introduce randint to avoid systematics)
            iz0=N.random.randint(nsz)
            ftz=f_mod[iz0::nsz,::n_interp+1,:]
            Gchi2=(foo-
                   N.add.reduce(ftz*fefef[N.newaxis,N.newaxis,:]  ,-1)**2
                /N.add.reduce((ftz/ef[N.newaxis,N.newaxis,:])**2,-1)
                )
            GPz=N.add.reduce(N.exp(-0.5*(Gchi2-Gchi2.min())),-1)
            GPz/=GPz.max()            
            #Interpolate calculate at points in which they != 0
            g=U.match_resol(z[iz0::nsz],GPz,z) > (float(pars.d["P_MIN"])/nt)
            eff+=g.sum()/float(nz)

            #Now run full probability analysis
            likelihood=N.zeros((nz,nt),"float")
            ft=N.compress(g,f_mod,axis=0)
            likelihood[g,:]=N.exp(-.5*(
                foo-
                N.add.reduce(fefef[N.newaxis,N.newaxis,:]*ft,-1)**2/
                N.add.reduce((ft/ef[N.newaxis,N.newaxis,:])**2,-1)
                ))
            lmax=likelihood.max()

            #P.plot(z,likelihood)
            #P.plot(z,p_i)
            #P.show()
            if lmax>0:
                chi2min=-2.*N.log(lmax)
            else:
                chi2min=99.
                likelihood=likelihood*0.+1.
                lmax=1.
            red_chi2=chi2min/float(nf-1.)
            likelihood/=lmax
            p=likelihood*N.exp(-0.5*red_chi2) #Normalize to red_chi2 level
            iz_ml,t_ml=U.loc2d(p,"max")
            pmax=p.max()
            if N.isnan(pmax) or N.isinf(pmax): 
                p=N.ones((nz,nt),"float")
                iz_ml,t_ml,red_chi2=0.,0.,0.
            p*=N.greater(p,pars.d["P_MIN"]*pmax)
        
    #if pb.sum()==0.:
    #    pb=N.exp(N.log(p_i+tiny)+N.log(p+tiny))
    #    #pb=N.exp(N.log(p_i)+N.log(p))
    #    P.contourf(pb,500)
    #    P.show()

    if pars.d['P_CAT']: 
        pb=p_i*(p+pC) 
    else:
        pb=p_i*p

    pbmax=pb.max()
    # Clip low prob values
    pb*=N.greater(pb,pars.d["P_MIN"]*pbmax)
    norm=pb.sum()
    pb/=norm 


    if hdf5:
        #full_probs[ig,:,:]=     N.where(pb[:nz,:nt]>pbmax*float(pars.d["P_MIN"]),pb,0.)
        #full_likelihood[ig,:,:]=N.where( p[:nz,:nt]>pmax *float(pars.d["P_MIN"]),p, 0.)
        full_probs[ig,:,:]=      pb
        full_likelihood[ig,:,:]=  p

    #Collapse redshift probability along nt
    p_bayes=N.add.reduce(pb[:nz,:nt],-1)

    # Convolve with a gaussian
    if pars.d['CONVOLVE_P']=='yes': 
        if pars.d["ONLY_TYPE"]=="yes":
            sigma_g=pars.d["SIGMA_CONVOLVE"]*(1.+z_s[ig])
            x=N.arange(-3.*sigma_g,3.*sigma_g+dz/100.,dz)
            gaus=N.exp(-(x/sigma_g)**2)
        p_bayes=N.convolve(p_bayes,gaus,1)    

    #if pars.d['CONVOLVE_P']=='yes' and pars.d['ONLY_TYPE']=='no': p_bayes=convolve(p_bayes,gaus,1)    
        
    # Eliminate all low level features in the prob. distribution
    #if pars.d["P_MIN"]>0.:
    #    pmax=max(p_bayes)
    #    p_bayes=where(greater(p_bayes,pmax*float(pars.d['P_MIN'])),p_bayes,0.)
    
    p_bayes=p_bayes/p_bayes.sum()

    if specprob:
        p_spec[ig,:]=U.match_resol(z,p_bayes,z_spec[ig,:])*p_spec[ig,:]
        norma=p_spec[ig,:].sum()
        if norma==0.: norma=1.
        p_spec[ig,:]/=norma
        vyjod=tuple([iden[ig]]+list(z_spec[ig,:])+list(p_spec[ig,:]))
        formato="%s "+5*" %.4f"
        formato+=5*" %.3f"
        formato+="\n"
        print formato % vyjod
        specout.write(formato % vyjod)

    if pars.d['N_PEAKS']>1:
        # Identify  maxima and minima in the final probability
        g_max=N.less(p_bayes[2:],p_bayes[1:-1])*N.less(p_bayes[:-2],p_bayes[1:-1])
        g_min=N.greater(p_bayes[2:],p_bayes[1:-1])*N.greater(p_bayes[:-2],p_bayes[1:-1])
    
        g_min+=N.equal(p_bayes[1:-1],0.)*N.greater(p_bayes[2:],0.)
        g_min+=N.equal(p_bayes[1:-1],0.)*N.greater(p_bayes[:-2],0.)
    
        i_max=N.compress(g_max,N.arange(nz-2))+1
        i_min=N.compress(g_min,N.arange(nz-2))+1                      

        # Check that the first point and the last one are not minima or maxima,
        # if they are, add them to the index arrays

        if p_bayes[0]>p_bayes[1]:
            i_max=N.concatenate([[0],i_max])
            i_min=N.concatenate([[0],i_min])
        if p_bayes[-1]>p_bayes[-2]:
            i_max=N.concatenate([i_max,[nz-1]])
            i_min=N.concatenate([i_min,[nz-1]])
        if p_bayes[0]<p_bayes[1]:
            i_min=N.concatenate([[0],i_min])
        if p_bayes[-1]<p_bayes[-2]:
            i_min=N.concatenate([i_min,[nz-1]])

        jm=N.searchsorted(i_min,i_max)
        p_tot=[]
        for i in range(len(i_max)):
            if jm[i]>0.:
                p_tot.append(N.sum(p_bayes[i_min[jm[i]-1]:i_min[jm[i]]]))
            else:
                p_tot.append(N.sum(p_bayes[i_min[jm[i]]:i_min[jm[i]]]))
                
        p_tot,i_max=U.multisort(-N.array(p_tot),(p_tot,i_max))
        z_min=[]
        z_max=[]
        z_peak=[]
        t_peak=[]
        i_peak=[]
        jm=N.searchsorted(i_min,i_max)

        for i in range(len(i_max)):
            z_min.append(z[i_min[jm[i]-1]])
            z_max.append(z[i_min[jm[i]]])
            z_aver=N.sum(p_bayes[i_min[jm[i]-1]:i_min[jm[i]]]*z[i_min[jm[i]-1]:i_min[jm[i]]])/p_tot[i]
            if jm[i]>0.:
                z_med=U.match_resol(N.add.accumulate(p_bayes[i_min[jm[i]-1]:i_min[jm[i]]])/p_tot[i],z[i_min[jm[i]-1]:i_min[jm[i]]],0.5)
            else:
                z_med=z[i_min[jm[i]]]
                
            if pars.d["MP"][0].lower()=="y":
                z_peak.append(z_aver)
            else:
                z_peak.append(z_med)
            i_peak.append(N.searchsorted(z,z_peak[-1]))
            # WE TAKE THE TYPE WHICH CORRESPONDS TO THE REDSHIFT ESTIMATE (NO MARGINALIZING HERE)
            t_peak.append(pb[i_peak[-1],:nt].argmax())

        if n_interp:
            t_peak=list(N.array(t_peak)/(1.+n_interp))
            
        if pars.d['MERGE_PEAKS']=='yes':
            # Merge peaks which are very close 
            merged=N.zeros(len(z_peak))
            for k in range(len(z_peak)):
                #Sanity check
                if z_min[k]==zmax or z_max[k]==zmin: merged[k]=1.
                for j in range(len(z_peak)):
                    if j>k and merged[k]==0 and merged[j]==0:
                        if abs(z_peak[k]-z_peak[j])<float(pars.d["SIGMA_EXPECTED"])*(1.+z_peak[j]):
                            # Modify the element which receives the accretion
                            z_peak[k]=(z_peak[k]*p_tot[k]+z_peak[j]*p_tot[j])/(p_tot[k]+p_tot[j])
                            i_peak[k]=N.searchsorted(z,z_peak[k])
                            z_min[k]=N.minimum(z_min[k],z_min[j])
                            z_max[k]=N.maximum(z_max[k],z_max[j])
                            p_tot[k]+=p_tot[j]
                            #print z_min[k],z_max[k],z_min[j],z_max[j],p_tot[k]
                            #WE DON'T CHANGE THE TYPE!!
                            # Put the merged element in the list
                            merged[j]=1
                            
 
            #print z_min
            #print z_max

            rem=N.logical_not(merged)
            p_tot=p_tot[rem]
            z_min=N.array(z_min)[rem]
            z_max=N.array(z_max)[rem]
            z_peak=N.array(z_peak)[rem]
            t_peak=N.array(t_peak)[rem]
            i_peak=N.array(i_peak)[rem]
            

        # Now sort them properly
        p_tot,z_min,z_max,z_peak,t_peak,i_peak=U.multisort(-N.array(p_tot),(p_tot,z_min,z_max,z_peak,t_peak,i_peak))

        #Cap them to N=N_PEAKS
        npe=pars.d["N_PEAKS"]
        p_tot=p_tot[:npe]
        z_peak=z_peak[:npe]
        z_min=z_min[:npe]
        z_max=z_max[:npe]
        t_peak=t_peak[:npe]
        i_peak=i_peak[:npe]

        #Renormalize probabilities (SUM MAY BE DIFFERENT FROM 0)
        p_tot=N.array(p_tot)/N.array(p_tot).sum()


        
        #Now find proper confidence intervals for each peak
        #print z_min
        #print z_max
        for ip in range(len(z_peak)):
            if z_min[ip]<>z_max[ip] and z_min[ip]<z_max[ip]:
                imin=N.searchsorted(z,z_min[ip])
                imax=N.searchsorted(z,z_max[ip])
                pacc=N.add.accumulate(p_bayes[imin:imax])/p_bayes[imin:imax].sum()
                zx=z[imin:imax]
                # 66% confidence interval
                z_min[ip]=U.match_resol(pacc,zx,0.17)
                z_max[ip]=U.match_resol(pacc,zx,0.83)
                    

        # Collapse stellar mass probability along nt, but fixing the redshift
        if stmass:
            smzb_peak=[]
            for ip in range(len(z_peak)):
                # Scale stellar mass by magnitude
                smzb_peak.append(N.log10(N.sum(
                            sm20[i_peak[ip],:]*pb[i_peak[ip],:],-1)/N.sum(pb[i_peak[ip],:],-1
                                                                        )*B.mag2flux(m_0[ig])/B.mag2flux(20.)))
            
        if absmag:
            mabszb_peak=[]
            for ip in range(len(z_peak)):
                # This is mabs(|z_B)
                mabszb_peak.append(B.flux2mag(N.sum(B.mag2flux(B.mabs20[i_peak[ip],:])*pb[i_peak[ip],:],-1)/N.sum(pb[i_peak[ip],:],-1))-20.+m_0[ig])

                
    # Define the peak
    if pars.d['ONLY_TYPE']=='yes': # Use only the redshift information, no priors
        zb=z_s[ig]
        iz_b=abs(z-z_s[ig]).argmin()
    else:
        #Weighted redshift
        zmp=N.add.reduce(z*p_bayes) # Weighted redshift
        #Modal redshift
        zmod=z[p_bayes.argmax()]
        #Median redshift
        z_med=U.match_resol(N.add.accumulate(p_bayes),z,0.5)
        if pars.d["MP"][0].lower()=="y": 
            zb=zmp
        else: 
            zb=z_med
        # If the average redshift is too different from the modal one, use the modal
        if abs(zmod-zb)< dz_exp*(1.+zb): zb=zmod            
        iz_b=abs(z-zb).argmin()
        # print 'zb:',zb
        # print 'iz_b:',iz_b
        # print 'z[iz_b]:',z[iz_b]
        # pausa = raw_input('pausa')

    #Now define main quantities in an environment of +-dz_exp of the redshift estimate
    # Calculate point estimate of the odds if dz_exp=0 
    if dz_exp==0.:
        # At the spectroscopic redshift is available and less than 9
        if spec and z_s[ig]<9.9:
            o=U.match_resol(z,p_bayes,z_s[ig])/dz
        else:
            o=U.match_resol(z,p_bayes,zb)/dz            
    else:
        #Integrate probability within expected limits around the redshift estimate
        zo1=zb-dz_exp*(1.+zb)
        zo2=zb+dz_exp*(1.+zb)
        if pars.d['Z_THR']>0:
            zo1=float(pars.d['Z_THR'])
            zo2=float(pars.d['ZMAX'])
        o=B.odds(p_bayes[:nz],z,zo1,zo2)
        
    ## The type is calculated integrating in a SIGMA_EXPECTED environment of the best type        
    #izo1=N.searchsorted(z,zo1)
    #izo2=N.searchsorted(z,zo2)
    #t_b=N.sum(pb[izo1:izo2,:nt]*xt)/N.sum(pb[izo1:izo2,:nt])

    # The type is calculated at the best redshift, collapsing the probability so it is p(t|z_b)    
    # NEW from August 2013    
    t_b=(pb[iz_b,:]*xt).mean()/pb[iz_b,:].mean()
    it_b=N.searchsorted(xt,t_b)
    # print 't_b,it_b',t_b,it_b
    # pausa = raw_input('paused')
    if N.any(N.isnan(t_b)) or pars.d["ONLY_TYPE"][0].lower()=="y":
        it_b=pb[iz_b,:nt].argmax()
        t_b=xt[it_b]
    
    # Redshift confidence limits (defined by "ODDS")
    if pars.d['ONLY_TYPE']=="yes":
        z1=zb
        z2=zb
    else:
        z1=U.match_resol(N.add.accumulate(p_bayes),z,0.5-float(pars.d['ODDS'])/2.)
        z2=U.match_resol(N.add.accumulate(p_bayes),z,0.5+float(pars.d['ODDS'])/2.)
        #z1,z2=interval(p_bayes[:nz],z,0.6827)

    tt_b=xt[t_b]
    tt_ml=xt[t_ml]

    #Stellar mass probability
    if stmass:
        # Scale stellar mass by magnitude
        #smzb=N.log10(sum(sm20[izo1:izo2,:]*pb[izo1:izo2,:])/sum(pb[izo1:izo2,:])*B.mag2flux(m_0[ig])/B.mag2flux(20.))
        if m_0[ig]<=0:
            smzb=0.
        else:
            if pars.d['ONLY_TYPE']=="yes":
                smzb=sm20[iz_b,it_b]
            else:
                #smzb=(sm20[izo1:izo2,:]*pb[izo1:izo2,:]).mean()/pb[izo1:izo2,:].mean()
                #We use the photo-z point estimate to calculate the stellar mass
                smzb=sm20[iz_b,it_b]
                
            smzb=N.log10(smzb*B.mag2flux(m_0[ig])/B.mag2flux(20.))
            
    if absmag:
        #This is mabs(|z_B)
        if m_0[ig]<=0:
            mabszb=0.
        else:
            if pars.d['ONLY_TYPE']=="yes":
                mabszb=mabs20[iz_b,it_b]-20.+m_0[ig]
            else:
                #mabszb=B.flux2mag(
                #    (B.mag2flux(mabs20[izo1:izo2,:])*pb[izo1:izo2,:]).mean()/pb[izo1:izo2,:].mean()-20.+m_0[ig])
                #We use the photo-z point estimate
                mabszb=mabs20[iz_b,it_b]-20.+m_0[ig]
                # print 'mabs20[iz_b,it_b]',mabs20[iz_b,it_b]
                # print 'iz_b,it_b',iz_b,it_b
                # print 'm_0[ig]:',m_0[ig]
                # print 'mabszb=mabs20[iz_b,it_b]-20.+m_0[ig]: ',mabszb
                # pausa = raw_input('paused')
            
    if pars.d['N_PEAKS']==1:
        salida=[iden[ig],zb,z1,z2,t_b+1,o]
        if stmass: salida=salida+[smzb]
        if absmag: salida=salida+[mabszb]
        salida+=[z[iz_ml],tt_ml+1,red_chi2]
    else:
        salida=[iden[ig]]
        for k in range(pars.d['N_PEAKS']):
            if k<= len(p_tot)-1:
                salida=salida+[z_peak[k],z_min[k],z_max[k],t_peak[k]+1,p_tot[k]]
                if stmass: 
                    salida=salida+[smzb_peak[k]]
                if absmag: 
                    salida=salida+[mabszb_peak[k]]
            else:
                salida+=[-1.,-1.,-1.,-1.,-1.]
                if stmass:
                    salida=salida+[-1.]
                if absmag: 
                    salida=salida+[-1.]

        salida+=[z[iz_ml],tt_ml+1,red_chi2]
        
    if spec:salida.append(z_s[ig])
    if has_mags: salida.append(m_0[ig])
    if col_pars.d.has_key('OTHER'):salida.append(other[ig])

    if get_z: 
        output.write(format % tuple(salida)+'\n')
    if pars.d['VERBOSE']=='yes': 
        print format % tuple(salida)
    else:
        if ig/100==float(ig/100.): print ".",

    if check=='yes':
        ft=f_mod[iz_b,it_b,:]
        fo=f_obs[ig,:]
	efo=ef_obs[ig,:]	
	factor=ft/efo/efo
	ftt=N.add.reduce(ft*factor)
	fot=N.add.reduce(fo*factor)
	am=fot/ftt
	ft*=am   

        flux_comparison=[iden[ig],m_0[ig],z[iz_b],t_b+1,am]+list(N.concatenate([ft,fo,efo]))
	nfc=len(flux_comparison)

	format_fc='%s  %.2f  %.2f   %.2f'+(nfc-4)*'   %.3e'+'\n'
	buffer_flux_comparison=buffer_flux_comparison+ format_fc % tuple(flux_comparison)

    if save_probs:
        texto='%s ' % str(iden[ig])
        texto+= len(p_bayes)*'%.3e '+'\n'
        probs.write(texto % tuple(p_bayes))
        
if check=='yes' and get_z:
    header_flux_comparison="#  1 ID\n#  2 M_0\n#  3 Z_B\n#  4 T_B\n#  5 A_norm\n"
    #nfil=len(filters)
    for ik in range(6,nf+6):
        header_flux_comparison+="# %2i F_T(%s)\n" % (ik,filters[ik-6])
    for ik in range(6+nf,6+nf*2):
        header_flux_comparison+="# %2i F_obs(%s)\n" % (ik,filters[ik-nf-6])
    for ik in range(6+nf*2,6+nf*3):
        header_flux_comparison+="# %2i dF_obs(%s)\n" % (ik,filters[ik-2*nf-6])
    open(pars.d['FLUX_COMPARISON'],'w').write(header_flux_comparison)
    open(pars.d['FLUX_COMPARISON'],'a').write(buffer_flux_comparison)

if get_z:    output.close()

if save_probs: probs.close()
    

if check=='yes':
    zb,zm,zb1,zb2,o,tb=U.get_data(out_name,(1,6,2,3,5,4))
    #Plot the comparison between z_spec and z_B

    if spec and inter:
        if U.ask('Compare z_B vs z_spec?'):
            good=N.less(abs(z_s),9.99)
            print >> resfile, 'Total initial number of objects with spectroscopic redshifts= ',N.asarray(good,'int').sum()
            print 'Total initial number of objects with spectroscopic redshifts= ',N.asarray(good,'int').sum()
            od_th=float(pars.d["ODDS"]) 
            if not U.ask("By default we exclude objects with odds > %.2f\n OK?" % od_th):
                od_th=input('Odds threshold? ')
                if has_mags:
                    mg_min=input('Bright magnitude limit?  ')
                    mg_max=input('Faint magnitude limit?  ')
                zt_min=input('Minimum redshift limit?  ')
                zt_max=input('Maximum redshift limit?  ')
                t_min=input('Minimum type?  ')
                t_max=input('Maximum type?  ')
                good*=less(m_0,mg_max)*N.greater_equal(m_0,mg_min)*N.less(z_s,zt_max)*N.greater_equal(z_s,zt_min)*N.less(tb,t_max)*N.greater_equal(tb,t_min)
            print good.shape
            print >> resfile, "Total number of galaxies,odds", good.shape[0], od_th    
            good=good*N.greater_equal(o,od_th)
            #print "nt0=",nt0
            for t in N.arange(.5,nt0+1.5):
                if t< nt0+0.5: 
                    print >> resfile, "Type %i" % (int(t)+1.,),
                    g=good*N.less_equal(tb,t+1.)*N.greater_equal(tb,t)
                    if g.sum()==0:
                        print >> resfile,  "0 galaxies of this type"
                        continue
                else:
                    print >> resfile, "All types",
                    g=good
                    if g.sum()==0.:
                        print >> resfile, "Odds too restrictive"
                        print >> resfile, "Using odds=0"
                        g=o*0.+1
                
                zmo,zso,zbo,zb1o,zb2o,tbo=U.multicompress(g,(zm,z_s,zb,zb1,zb2,tb))
                #print 'Number of objects with odds > %.2f= %i '% (od_th,len(zbo))
                print >> resfile, len(zbo),
                deltaz=(zbo-zso)/(1.+zso)
                rms=U.std_mad(deltaz)
                midpt=N.median(deltaz)
                #sz=stat_robust(deltaz,3.,3)
                #sz.run()
                #outliers=asarray(abs(deltaz-midpt) > 3.*0.04,'int')
                outliers=N.asarray(abs(deltaz-midpt) > 5.*rms,'int')
                #print 'Number of outliers [dz >%.2f*(1+z)]=%i' % (3.*sz.rms,add.reduce(outliers))
                print >> resfile,  outliers.sum(),
                catastrophic=N.asarray(N.greater_equal(deltaz*(1.+zso),1.),'int')
                n_catast=catastrophic.sum()
                #print 'Number of catastrophic outliers [dz >1]=',n_catast            
                print >> resfile, n_catast,            
                print >> resfile, '   Delta z/(1+z) = %.4f +- %.4f' % (midpt,rms)
                if inter and plots:
                    P.figure(2)
                    P.subplot(211)
                    P.plot(N.arange(min(zso),max(zso)+0.01,0.01),
                         N.arange(min(zso),max(zso)+0.01,0.01),
                         "r")
                P.errorbar(zso,zbo,[abs(zbo-zb1o),abs(zb2o-zbo)],fmt="bo")
                P.xlabel(r'$z_{spec}$')
                P.ylabel(r'$z_{bpz}$')
                P.subplot(212)
                P.plot(zso,zmo,"go",zso,zso,"r")
                P.xlabel(r'$z_{spec}$')
                P.ylabel(r'$z_{ML}$')
            if plots: P.show()
            
resfile.close()
os.system("cat %s" % resfile_name)

if hdf5: fp_file.close()

if pars.d["VERBOSE"][0].lower()=="y":
    rolex.check()
    print "Elapsed time after initialization",rolex.t()-tinit,"s"
print "Average eff",eff/float(ng)
