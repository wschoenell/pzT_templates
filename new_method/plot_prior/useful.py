#Useful functions and definitions

import os,sys
from types import *
from time import *
import scipy
from scipy import optimize
from scipy import signal
#Some functions already use N, others do not, clean this up
from numpy import *
import numpy as np
import string
from string import *
# import pylab as pl

pi=3.14159265358979323846

def ejecuta(command=None,verbose=1,dry=0):
    import os
    if verbose: 
        print
        print command
        print
    if not dry:
        r=watch()
        os.system(command)
        if verbose: print r.t(),"s"

def niet(file):
    return not os.path.exists(file)    

def ask(what="?"):
    """
    Usage:
    ans=ask(pregunta)
    This function prints the string wjat, 
    (usually a question) and asks for input 
    from the user. It returns the value 0 if the 
    answer starts by 'n' and 1 otherwise, even 
    if the input is just hitting 'enter'
    """
    if what[-1]<>'\n': what=what+'\n'
    ans=raw_input(what)
    try:
	if ans[0]=='n': return 0
    except:
	pass
    return 1
    
#Input/Output subroutines

#Read/write headers

def get_header(file):
    """ Returns a string containing all the lines 
    at the top of a file which start by '#'"""
    buffer=''
    for line in open(file).readlines():
	if line[0]=='#': buffer=buffer+line
	else: break
    return buffer

def put_header(file,text,comment=1):
    """Adds text (starting by '#' and ending by '\n')
    to the top of a file."""
    if len(text)==0: return
    if text[0]<>'#' and comment: text='#'+text
    if text[-1]<>'\n':text=text+'\n'
    buffer=text+open(file).read()
    open(file,'w').write(buffer)

#Files containing strings

def get_str(file,cols="all",nrows='all',separator=None):
    """ 
        Reads strings from a file
        Usage: 
	     x,y,z=get_str('myfile.cat',(0,1,2))
             x,y,z are returned as string lists
    """
    if type(cols)==type(0):
	cols=(cols,)
	nvar=1
    else: nvar=len(cols)
    lista=[]
    for i in range(nvar): lista.append([])
    buffer=open(file).readlines() 
    if nrows=='all': 
        nrows=ones(len(buffer))
    elif type(nrows)==type(100):
        nrows=ones(nrows)        
    ln=len(nrows)

    counter=0
    for lines in buffer:
        if lines[0]=='#': continue # Do not count lines with headers
	if not nrows[counter]: 
            counter+=1
            continue 
        if separator:
            pieces=string.split(lines,separator)
        else:
            pieces=string.split(lines)
 	if len(pieces)==0: continue
   	for j in range(nvar):lista[j].append(pieces[cols[j]])
	counter+=1
        if counter>=ln: break

    if nvar==1: return lista[0]
    else: return tuple(lista) 

def get_str_test(file,cols="all",nrows='all',separator=None):
    """ 
        Reads strings from a file
        Usage: 
	     x,y,z=get_str('myfile.cat',(0,1,2))
             x,y,z are returned as string lists
    """
    if type(cols)==type(0):
	cols=(cols,)
	nvar=1
    else: nvar=len(cols)
    lista=[]
    for i in range(nvar): lista.append([])
    buffer=open(file).readlines() 
    if nrows=='all': 
        nrows=ones(len(buffer))
    elif type(nrows)==type(100):
        nrows=ones(nrows)
    ln=len(nrows)

    counter=0
    for lines in buffer:
        if counter>ln: break
        if lines[0]=='#': continue # Do not count lines with headers
	if nrows[counter]: 
            if separator: 
                pieces=string.split(lines,separator)
            else:
                pieces=string.split(lines)
            if len(pieces)==0: continue
            for j in range(nvar):lista[j].append(pieces[cols[j]])
            counter+=1

    if nvar==1: return lista[0]
    else: return tuple(lista) 

def put_str(file,tupla):
    """ Writes tuple of string lists to a file
        Usage:
	  put_str(file,(x,y,z))
    """    
    if type(tupla)<>type((2,)):
	raise 'Need a tuple of variables'
    f=open(file,'w')    
    for i in range(1,len(tupla)):
	if len(tupla[i])<>len(tupla[0]):
	    raise 'Variable lists have different lenght'
    for i in range(len(tupla[0])):
	cosas=[]
	for j in range(len(tupla)):cosas.append(str(tupla[j][i]))
	f.write(join(cosas)+'\n')
    f.close()

#Files containing data

def get_data(file,cols=0,nrows='all',separator=None):
    """ Returns data in the columns defined by the tuple
    (or single integer) cols as a tuple of float arrays 
    (or a single float array)"""

    if cols=='all':
        #Get the number of columns in the file
	for line in open(file).readlines():
            if separator:
                pieces=line.split(separator)
            else:
                pieces=line.split()
	    if len(pieces)==0: continue
	    if line[0]=='#':continue
	    cols=range(len(pieces))
	    break
    elif type(cols)==type(0):
	cols=(cols,)
    data=get_str(file,cols,nrows,separator=separator)
    nvar=len(cols)
    if nvar==1: return array(map(float,data))
    else:
	data=list(data)
	for j in range(nvar): 
            try:
                data[j]=array(map(float,data[j]))
            except:
                print j,data[j]
	return tuple(data) 


#def get_data(file,cols=False,nrows="all",separator=None):
#    #Improved version suggested by Martin Borstad Eriksen
#    cols=(cols,) if isinstance(cols,int) else cols
#    nvar=len(cols)    
#    print cols
#    if nrows=="all": 
#        if cols:
#            return loadtxt(file,usecols=cols,unpack=True)
#        else:
#            return loadtxt(file,unpack=True)
#    else:
#        data=get_str(file,cols,nrows,separator=separator)
#        if nvar==1: return array(map(float,data))
#        else:
#            data=list(data)
#            for j in range(nvar): 
#                try:
#                    data[j]=array(map(float,data[j]))
#                except:
#                    print j,data[j]
#            return tuple(data) 
        
        
def put_data(file,variables,header='',format='',append='no'):
    """ Writes tuple of list/arrays to a file 
        Usage:
	  put_data(file,(x,y,z),header,format)
	where header is any string  
        and format is a string of the type:
           '%f %f %i ' 
	The default format is all strings
    """    
    if type(variables)<>type((2,)):
	raise 'Need a tuple of variables'
    if format=='': format='%s  '*len(variables)
    if append=='yes': f=open(file,'a')
    else: f=open(file,'w')
    if header<>"":
	if header[0]<>'#': header='#'+header
	if header[-1]<>'\n': header=header+'\n'
	f.write(header)
    for i in range(len(variables[0])):
	cosas=[]
	for j in range(len(variables)):
	    cosas.append(variables[j][i])
	line=format % tuple(cosas)             
	f.write(line+'\n')
    f.close()


#Read/write 2D arrays

def get_2Darray(file,cols='all',nrows='all',verbose='no'):
    """Read the data on the defined columns of a file to an 2 array
    Usage:
    x=get_2Darray(file)
    x=get_2Darray(file,range(len(p))
    x=get_2Darray(file,range(0,10,2),nrows=5000)
    Returns x(nrows,ncols)
    """
    if cols=='all':
        #Get the number of columns in the file
	for line in open(file).readlines():
	    pieces=line.split()
	    if len(pieces)==0: continue
	    if line[0]=='#':continue
	    nc=len(pieces)
	    cols=range(nc)
            if verbose=='yes': print 'cols=',cols
	    break
    else:
	nc=len(cols)
    
    lista=get_data(file,cols,nrows)
    nl=len(lista[0])
    x=zeros((nl,nc))*1.
    for i in range(nc):
        x[:,i]=array(lista[i])
    return x

def put_2Darray(file,array,header='',format='',append='no'):
    """ Writes a 2D array to a file, where the first 
    index changes along the lines and the second along
    the columns
    Usage: put_2Darray(file,a,header,format)
	where header is any string  
        and format is a string of the type:
           '%f %f %i ' 
    """
    lista=[]
    for i in range(array.shape[1]):lista.append(array[:,i])
    lista=tuple(lista)
    put_data(file,lista,header,format,append)
       

class watch:
    def __init__(self):
	self.time0=time()        
    def t(self):
        return time()-self.time0
    def check(self):
        print
        print "Elapsed time", strftime('%H:%M:%S',gmtime(self.t()))
        print

def params_file(file,cap=1,numbers=1):
    """ 
	Read a input file containing the name of several parameters 
        and their values with the following format:
        
	KEY1   value1,value2,value3   # comment
        KEY2   value

        Returns the dictionary
        dict['KEY1']=(value1,value2,value3)
        dict['KEY2']=value
    """
    dict={}
    for line in open(file,'r').readlines():
        if line[0]==' ' or line[0]=='#': continue 
	halves=line.split('#')
	#replace commas in case they're present
	halves[0]=replace(halves[0],',',' ') 	
	pieces=halves[0].split()
	if len(pieces)==0: continue
        if cap:
            key=pieces[0].upper()
        else:
            key=pieces[0]            
#	if type(key)<>type(''):
#	    raise 'Keyword not string!'
	if len(pieces)<2:
	    mensaje='No value(s) for parameter  '+key
	    raise mensaje
        dict[key]=tuple(pieces[1:]) 
	if len(dict[key])==1: dict[key]=dict[key][0]

    if numbers:
        for k in dict.keys():
            try: dict[k]=float(dict[k])
            except: pass

    return dict

def params_commandline(lista,dict=None,cap=1):
    """ Read an input list (e.g. command line) 
	containing the name of several parameters 
        and their values with the following format:
        
        ['-KEY1','value1,value2,value3','-KEY2','value',etc.] 
          
         Returns a dictionary containing 
        dict['KEY1']=(value1,value2,value3)
        dict['KEY2']=value 
	etc.
    """
    if len(lista)%2<>0:
        print 'Error: The number of parameter names and values does not match'
        sys.exit()
    if dict==None: dict={}
    for i in range(0,len(lista),2):
        key=lista[i]
	if type(key)<>type(''):
	    raise 'Keyword not string!'
	#replace commas in case they're present
	if key[0]=='-':
            if cap:
                key=key[1:].upper()
            else:
                key=key[1:]                
        lista[i+1]=replace(lista[i+1],',',' ')
	values=tuple(lista[i+1].split())
        if len(values)<1:
	    mensaje='No value(s) for parameter  '+key
	    raise mensaje
        dict[key]=values
	if len(dict[key])==1: dict[key]=dict[key][0]
    return dict

def view_keys(dict):
    """Prints sorted dictionary keys"""
    claves=dict.keys()
    claves.sort()
    for line in claves:
	print upper(line),'  =  ',dict[line]

class params:
    """This class defines and manages a parameter dictionary"""
    def __init__(self,d=None,cap=1):
        if d==None:self.d={}
        else: self.d=d
        self.cap=cap

    # Define a few useful methods:

    def fromfile(self,file):
	"""Update the parameter dictionary with a file"""
	self.d.update(params_file(file,cap=self.cap)) 
    
    def fromcommandline(self,command_line):
	"""Update the parameter dictionary with command line options (sys.argv[i:])"""
	self.d.update(params_commandline(command_line,cap=self.cap))

    def update(self,dict):
	"""Update the parameter information with a dictionary"""
        for key in dict.keys():
            if self.cap:
                self.d[key]=dict[key.upper()]
            else:
                self.d[key]=dict[key]                
                   
    def check(self):
        """Interactively check the values of the parameters"""
        view_keys(self.d)
        paso1=raw_input('Do you want to change any parameter?(y/n)\n')
        while paso1[0] == 'y':
            key=raw_input('Which one?\n')
            if not self.d.has_key(key):
                paso2=raw_input("This parameter is not in the dictionary.\
Do you want to include it?(y/n)\n")
                if paso2[0]=='y':
                    value=raw_input('value(s) of '+key+'?= ')
                    self.d[key]=tuple(split(replace(value,',',' ')))
                else:continue
            else:
                value=raw_input('New value(s) of '+key+'?= ')
                self.d[key]=tuple(split(replace(value,',',' ')))
	    view_keys(self.d)
            paso1=raw_input('Anything else?(y/n)\n')            

    def write(self,file):
        claves=self.d.keys()
        claves.sort()
        buffer=''
        for key in claves:
            if type(self.d[key])==type((2,)):
                values=map(str,self.d[key])
                line=key+' '+string.join(values,',')
            else:
                line=key+' '+str(self.d[key])
            buffer=buffer+line+'\n'
        print line
        open(file,'w').write(buffer)

    def put_header(self,file):
        claves=self.d.keys()
        claves.sort()
        buffer=''
        for key in claves:
            if type(self.d[key])==type((2,)):
                values=map(str,self.d[key])
                line=key+' '+string.join(values,',')
            else:
                line=key+' '+str(self.d[key])
            buffer=buffer+"#"+line+'\n'
        put_header(file,buffer)

#List of colors from biggles:
def biggles_colors():
    try: import biggles
    except: pass
    return get_str('/home/txitxo/Python/biggles_colors.txt',0)

#Some miscellaneous numerical functions

def ascend(x):
    """True if vector x is monotonically ascendent, false otherwise 
       Recommended usage: 
       if not ascend(x): sort(x) 
    """
    return alltrue(greater_equal(x[1:],x[0:-1]))

#def fixed_grid(xmin,xmax,dx,x0):
#    # Define an equally spaced grid of points which contains xmin,xmax and th#e points in x0
#    x=arange(xmin,xmax+dx,dx)
#    l=list(x)
#    if l[-1]>xmax: l=l[:-1]
#    x=array(l)
#    indies=searchsorted(x,x0)
#    for j in range(len(x0)):
#        if x[indies[j]]==x0[j]:
#            continue
#        else:
#            l=l[:indies[j]]+[x0[j]]+l[indies[j]:]
#    return array(l)

def fixed_grid(xmin,xmax,dx,x0):
    n=int(round((xmax-xmin)/dx))+1
    x=sort(list(set(list(linspace(xmin,xmax,n))+list(x0))))    
    return array(x)

    # Define an equally spaced grid of points which contains xmin,xmax and the points in x0
    x=arange(xmin,xmax+dx,dx)
    l=list(x)
    if l[-1]>xmax: l=l[:-1]
    x=array(l)
    indies=searchsorted(x,x0)
    for j in range(len(x0)):
        if x[indies[j]]==x0[j]:
            continue
        else:
            l=l[:indies[j]]+[x0[j]]+l[indies[j]:]
    return array(l)


def match_resol_old(xg,yg,xf):
    """ 
    Interpolates and/or extrapolate yg, defined on xg, onto the xf coordinate set.
    Usage:
    ygn=match_resol(xg,yg,xf)
    """
    #If only one point available
    if len(xg)==1 and len(yg)==1: return xf*0.+yg
    
    ng=len(xg)
    d=(yg[1:]-yg[0:-1])/(xg[1:]-xg[0:-1])
    #Get positions of the new x coordinates
    ind=clip(searchsorted(xg,xf)-1,0,ng-2)
    try:
        len(ind)
        one=0
    except:
        one=1
        ind=array([ind])
    ygn=take(yg,ind)+take(d,ind)*(xf-take(xg,ind))
    if one: ygn=ygn[0]
    return ygn

def match_resol(xg,yg,xf):
    """ 
    Interpolates and/or extrapolate yg, defined on xg, onto the xf coordinate set.
    Usage:
    ygn=match_resol(xg,yg,xf)
    """
    if alltrue(xf>xg[0]) and alltrue(xg[-1]>xf): 
        return np.interp(xf,xg,yg)
    else:
        # If only one point available
        if len(xg)==1 and len(yg)==1: return xf*0.+yg
    
        ng=len(xg)
        d=(yg[1:]-yg[0:-1])/(xg[1:]-xg[0:-1])
        # Get positions of the new x coordinates
        ind=clip(searchsorted(xg,xf)-1,0,ng-2)
        try:
            len(ind)
            one=0
        except:
            one=1
            ind=array([ind])
        ygn=take(yg,ind)+take(d,ind)*(xf-take(xg,ind))
        if one: ygn=ygn[0]
        return ygn


def overlap(x,y):
    """Returns 1 if vectors x and y overlap, 0 otherwise"""
    if (x[0]<=y[-1] and x[-1]>y[0]) or (y[0]<=x[-1] and y[-1]>x[0]):
	return 1
    else: return 0

def match_objects(coords1,coords2,tail1=(),tail2=(),accuracy=1.):
    """
    where coords1 and coords2 are tuples containing 1-D arrays,
    and tail1 and tail2 are tuples containing sequences of 
    arbitrary types
    Usage:
    results=match_objects((x1,y1),(x2,y2),(a1,b1,c1),(d2,e2),accuracy=.5)
    It returns the sequence x1,y1,a1,b1,c1,d2,e2 for those objects 
    which have dist(x1,y1-x2,y2)< accuracy
    """
    acc2=accuracy**2
    nc=len(coords1)
    np1=len(coords1[0])
    np2=len(coords2[0])
    a1=array(coords1)
    a2=array(coords2)
    nt1=len(tail1)
    for i in range(nt1): 
	if len(tail1[i])<> np1: raise 'Not the same lenght as coordinates 1'
    nt2=len(tail2)
    for i in range(nt2): 
	if len(tail2[i])<> np2: raise 'Not the same lenght as coordinates 2'
    match=zeros(np1,'int')-1
    for j in range(np1):
	dist=add.reduce((a1[:,j,newaxis]-a2[:,:])**2)
	i_min=dist.argmin()
	if dist[i_min]<acc2:match[j]=i_min
    good=greater_equal(match,0)
    n1=compress(good,range(np1))    
    match=compress(good,match)
    a1=compress(good,a1)
    salida=list(a1)
    for i in range(nt1):
	if type(tail1[i][0])==type('si'):
	    t=[]
	    for j in n1: t.append(tail1[i][j])
	else:
	    t=take(tail1[i],n1)
	salida.append(t)
    for i in range(nt2):
	if type(tail2[i][0])==type('si'):
	    t=[]
	    for j in match: t.append(tail2[i][j])
	else:
            print tail2[i]
            print match
	    t=take(tail2[i],match)
	salida.append(t)
    return salida


def match_min(coords1,coords2,tail1=(),tail2=()):
    """
    where coords1 and coords2 are tuples containing 1-D arrays,
    and tail1 and tail2 are tuples containing sequences of 
    arbitrary types

    Usage:

    results=match_min((x1,y1),(x2,y2),(a1,b1,c1),(d2,e2))
    It returns the sequence x1,y1,a1,b1,c1,d2,e2, dist_min 
    where dist_min is the minimal value of dist(x1,y1-x2,y2)
    The match looks for the objects with minimal distance
    """
    nc=len(coords1)
    np1=len(coords1[0])
    np2=len(coords2[0])
    a1=array(coords1)
    a2=array(coords2)
    nt1=len(tail1)
    for i in range(nt1): 
	if len(tail1[i])<> np1: 
            print "i=",i
            raise 'Not the same lenght as coordinates 1'
    nt2=len(tail2)
    for i in range(nt2): 
	if len(tail2[i])<> np2:
            print "i=",i 
            raise 'Not the same lenght as coordinates 2'

    match=zeros(np1)-1

    dist_min=zeros(np1)*1.

    for j in range(np1):
	dist=np.sqrt(np.add.reduce((a1[:,j,newaxis]-a2[:,:])**2))
	i_min=dist.argmin()
        dist_min[j]=dist[i_min]
	match[j]=int(i_min)

    match=asarray(match,"int")

    salida=list(a1)

    for i in range(nt1):salida.append(tail1[i])
    
    for i in range(nt2):
        print i
	if type(tail2[i][0])==type('si'):
	    t=[]
	    for j in match: t.append(tail2[i][j])
	else:
	    t=tail2[i][match]
	salida.append(t)

    salida.append(dist_min)
    return tuple(salida)

def match_min2(coords1,coords2,tail1=(),tail2=()):
    """
    where coords1 and coords2 are tuples containing 1-D arrays,
    and tail1 and tail2 are tuples containing sequences of 
    arbitrary types

    Usage:

    results=match_min((x1,y1),(x2,y2),(a1,b1,c1),(d2,e2))
    It returns the sequence x1,y1,x2,y2,a1,b1,c1,d2,e2, dist_min 
    where dist_min is the minimal value of dist(x1,y1-x2,y2)
    The match looks for the objects with minimal distance
    """
    nc=len(coords1)
    np1=len(coords1[0])
    np2=len(coords2[0])
    a1=array(coords1)
    a2=array(coords2)
    nt1=len(tail1)
    for i in range(nt1): 
	if len(tail1[i])<> np1: raise 'Not the same lenght as coordinates 1'
    nt2=len(tail2)
    for i in range(nt2): 
	if len(tail2[i])<> np2: raise 'Not the same lenght as coordinates 2'
    match=zeros(np1)-1
    dist_min=zeros(np1)*1.
    x2=zeros(np1)*1.
    y2=zeros(np1)*1.
    for j in range(np1):
	dist=sqrt(add.reduce((a1[:,j,newaxis]-a2[:,:])**2))
	i_min=dist.argmin()
        dist_min[j]=dist[i_min]
        x2[j],y2[j]=a2[0,i_min],a2[1,i_min]
	match[j]=i_min
        
    salida=list(a1)
    salida.append(x2)
    salida.append(y2)

    for i in range(nt1):salida.append(tail1[i])
    
    for i in range(nt2):
	if type(tail2[i][0])==type('si'):
	    t=[]
	    for j in match: t.append(tail2[i][j])
	else:
            print match
	    t=tail2[match]
	salida.append(t)

    salida.append(dist_min)
    return tuple(salida)

def match_ID(ID1,ID2):
    #Return a logical vector (1 or 0) which tells 
    #whether objects in ID1 are contained in ID2
    return map(lambda x: x in ID2,ID1)


def dist(x,y,xc=0.,yc=0.):
    """Distance between point (x,y) and a center (xc,yc)"""
    return sqrt((x-xc)**2+(y-yc)**2)

def loc2d(a,extremum='max'):
    """ Locates the maximum of an 2D array
        Usage:
	max_vec=max_loc2d(a)
    """
    if extremum=="max":
        return unravel_index(a.argmax(),a.shape)
    elif extremum=="min":
        return unravel_index(a.argmin(),a.shape)
    else:
	raise 'Which extremum are you looking for?'
    
def histogram_old(a,bins):
    """
    Histogram of 'a' defined on the bin grid 'bins'
       Usage: h=hist(p,xp)
    """
    n=searchsorted(sort(a),bins)
    n=concatenate([n,[len(a)]])
    n=array(map(float,n))
#    n=array(n)
    return n[1:]-n[:-1]

#def hist2D(a,xbins,ybins):
#    """
#    Histogram of 'a' defined on the grid xbins X ybins
#       Usage: h=hist2D(p,xp,yp)
#       Points larger than xbins[-1],ybins[-1] are asigned to
#       the 'last' bin
#    """   
#    nx=len(xbins)
#    ny=len(ybins)
#    #We use searchsorted differenty from the 1-D case
#    hx=searchsorted(xbins,a)
#    hy=searchsorted(ybins,a)        
#    h=zeros((nx,ny))
#    for i in range(len(hx)):
#        for j in range(len(hy)):
#            h[hx[i],hy[i]]=+1
#    for k in range(len(a)):
#        for i in range(len(xbins)):
#            for j in range(len(ybins)):
#                if a[k]>xbins[i] and a[k]<xbins[i+1] \
#                   and a[k]>ybins[i] and a[k]< ybins[i+1]:
#                    h[i,j]=h[i,j]+1
#                    break
#                else:
                                        
def bin2D(x,y,
          xmin=None,xmax=None,
          ymin=None,ymax=None,
          nx=10,ny=10,
          p=None,norm=0):
    if xmin==None:
        xmin=x.min()
    if xmax==None:
        xmax=x.max()
    if ymin==None:
        ymin=y.min()
    if ymax==None:
        ymax=y.max()
    ix=array(map(round,(x-xmin)/(xmax-xmin)*nx),"int")
    iy=array(map(round,(y-ymin)/(ymax-ymin)*ny),"int")
    N=zeros((nx,ny),"float")
    for i in range(nx):
        for j in range(ny):
            g=equal(i,ix)*equal(j,iy)
            try:
                N[i,j]=p[g].sum()
            except:
                N[i,j]=g.sum()*1.
    if norm:
        return N/N.sum()
    else:
        return N

def minfil(x,y,xbins,f=5):
    xm=[]
    ym=[]
    nbins=len(xbins)
    for i in range(nbins):
	if i<nbins-1:
	    good=greater_equal(x,xbins[i])*less(x,xbins[i+1])
	else:
            good=greater_equal(x,xbins[-1])
        if sum(good)>14:
            nnp=sum(good)
            yp=sort(compress(good,y))[nnp/f]
            good=good*less_equal(y,yp)
            xm.append(median(compress(good,x)))
            ym.append(median(compress(good,y)))
        else:
            print 'Bin starting at xbins[%i] has %i points' % (i,sum(good))
            xm.append(median(compress(good,x)))
            ym.append(0.)
    return match_resol(array(xm),array(ym),x)


def mad(x):
    return median(abs(x-median(x)))


def std_mad_neg(x):
    #This calculates the equivalent to the rms of the negative fluctuations
    xm=median(x)
    return -1.4826*median(x[x<0])
    
def std_mad(x):
    return 1.4826*mad(x)

def negmirror(x):
    return concatenate((x[x<0],-x[x<0]))

def posmirror(x):
    return concatenate((x[x>=0],-x[x>=0]))

def std_ci(x,n_sigma=1.):
    xm=median(x)
    p=confidence_interval(n_sigma)
    return match_resol(linspace(0.,1.,len(x)+1)[1:],sort(abs(x-xm)),p)/float(n_sigma)

def sigma_mixture(s,w=None):
    s=array(s)
    if w==None: 
        w=s*0.+1./len(s)
    else: 
        w=array(w,'float')
        w*=(1./sum(w))
    return sqrt(sum(w*s**2))

#This is much more robust than the weighted mean
def weighted_median(x,w):
    x,w=multisort(x,(x,w))
    s=add.accumulate(w)/add.reduce(w)
    return match_resol(s,x,0.5)

def weighted_mean(x,w):
    return add.reduce(x*w)/add.reduce(w)

def weighted_rms(x,w):
    center=weighted_mean(x,w)
    return sqrt(add.reduce((x-center)**2*w**2)/add.reduce(w**2))

def pad_data(x,y,n_points=2,value=0.,dl=20.,l=1,r=1):
    n=arange(n_points)
    xleft=x.min()-arange(n_points,0,-1)*dl
    yleft=xleft*0.+value
    xright=x.max()+arange(1,n_points+1)*dl
    yright=xright*0.+value        
    if l:
        x=concatenate((xleft,x))
        y=concatenate((yleft,y))
    if r:
        x=concatenate((x,xright))
        y=concatenate((y,yright))
    return x,y

def median_filter(x,y,dx_kernel=None,nx_kernel=11,pad=1,ypad=0.):
    #For two vectors x,y, generate a vector ym which 
    #is the median of either nx_kernel points or 
    # the median of the points within dx_kernel
    dk=int((nx_kernel-1.)/2)
    xs,ys=multisort(x,(x,y))
    if pad:
        xs,ys=pad_data(xs,ys,n_points=dk,value=ypad)
          
    ym=ys*0.
    for i in range(len(xs)):
        i1=max(i-dk,0)
        i2=min(i+dk+1,len(xs)-1)
        ym[i]=median(ys[i1:i2])
    return xs,ym
                    

def mean_robust_filter(x,y,dx_kernel=None,nx_kernel=11,pad=1,ypad=0.):
    #For two vectors x,y, generate a vector ym which 
    #is the median of either nx_kernel points or 
    # the median of the points within dx_kernel
    # and which takes into account the position of the points

    dk=int((nx_kernel-1.)/2)
    xs,ys=multisort(x,(x,y))
    if pad:
        xs,ys=pad_data(xs,ys,n_points=dk,value=ypad)          
    ym=ys*0.
    xm=ys*0.
    for i in range(len(xs)):
        i1=max(i-dk,0)
        i2=min(i+dk,len(xs)-1)
        ym[i]=mean_robust(ys[i1:i2])
        xm[i]=mean_robust(xs[i1:i2])
    return xm,ym
                    

def weighted_median_filter(x,y,w,dx_kernel=None,nx_kernel=11,pad=1,ypad=0.):
    #For two vectors x,y, generate a vector ym which 
    #is the weighted mean or median of either nx_kernel points or 
    # the the points within dx_kernel
    # and which takes into account properly the position of the points
    # It won't have the same dimensions

    dk=int((nx_kernel-1.)/2)
    x0,y0,w0=multisort(x,(x,y,w))
    if pad:
        xs,ys=pad_data(x0,y0,n_points=dk,value=ypad)          
        xxs,ws=pad_data(x0,w0,n_points=dk,value=ypad)          
    ym=ys*0.
    xm=ys*0.
    xmm=xm*0.
    ymm=ym*0.
    for i in range(len(xs)):
        i1=max(i-dk,0)
        i2=min(i+dk,len(xs)-1)
        y=ys[i1:i2]
        x=xs[i1:i2]
        w=ws[i1:i2]
        s=std_mad(y)
        g=less(abs(y-median(y)),3.*s)
        xmm[i]=mean(x)
        ymm[i]=median(y)
        if g.sum()>10.: 
            xm[i]=weighted_mean(x[g],w[g])
            ym[i]=weighted_median(y[g],w[g])
        else:
            xm[i]=xmm[i]
            ym[i]=ymm[i]            
        if abs(ym[i])>1.:
            xm[i]=xmm[i]
            ym[i]=ymm[i]

    xm,ym=multisort(xm,(xm,ym))
    return xm[1:-1],ym[1:-1]
                    

def bin_stats(x,y,xbins,stat='average',out_n_points=0,verbose=0):
    """Given the variable y=f(x), and 
    the bins limits xbins, return the 
    corresponding statistics, e.g. <y(xbins)>
    Options are rms, median y average
    """
    nbins=len(xbins)
    if   stat=='average' or stat=='mean':          func=mean
    elif stat=='median':                           func=median
    elif stat=='rms' or stat=='std'              : func=std
    elif stat=='std_mad'                         : func=std_mad
    elif stat=='std_ci'                         : func=std_ci
    elif stat=='std_robust' or stat=='rms_robust': func=std_robust
    elif stat=='mean_robust':                      func=mean_robust
    elif stat=='median_robust':                    func=median_robust
    elif stat=='sum':                              func=sum
    elif stat=='product':                          func=product
    elif stat=='sigma_mixture':                    func=sigma_mixture
    elif stat=='n_outliers':                       func=n_outliers
    results=[]
    nnp=[]
    for i in range(nbins):
	if i<nbins-1:
	    good=asarray(greater_equal(x,xbins[i])
                         *less(x,xbins[i+1]),int)
	else: good=asarray(greater_equal(x,xbins[-1]),int)
        if sum(good)>1.:
            results.append(func(compress(good,y)))
        else:
            results.append(0.)
            if verbose:
                print 'Bin starting at xbins[%i] has %i points' % (i,sum(good))
        nnp.append(good.sum())
    if not out_n_points:
        return array(results)
    else:
        return array(nnp),array(results)


def bin_aver(x,y,xbins):
    return bin_stats(x,y,xbins,stat='average')

def autobin_stats(x,y,n_bins=8,stat='average',n_points=None,xmed=0):
    """
    Given the variable y=f(x), form n_bins, distributing the 
    points equally among them. Return the average x position 
    of the points in each bin, and the corresponding statistic stat(y).
    n_points supersedes the value of n_bins and makes the bins 
    have exactly n_points each
    Usage:
      xb,yb=autobin_stats(x,y,n_bins=8,'median')
      xb,yb=autobin_stats(x,y,n_points=5)
    """
    
    if not ascend(x):
	ix=np.argsort(x)
	x=np.take(x,ix)
	y=np.take(y,ix)
    n=len(x)
    if n_points==None: 
	#This throws out some points
	n_points=n/n_bins
    else: 
	n_bins=n/n_points
	#if there are more that 2 points in the last bin, add another bin
	if n%n_points>2: n_bins=n_bins+1

    if n_points<=1:
	print 'Only 1 or less points per bin, output will be sorted input vector with rms==y'
	return x,y
    xb,yb=[],[]

    if   stat=='average' or stat=='mean':          func=mean
    elif stat=='median':                           func=median
    elif stat=='rms' or stat=='std'              : func=std
    elif stat=='std_robust' or stat=='rms_robust': func=std_robust
    elif stat=='std_mad'                         : func=std_mad
    elif stat=='mean_robust':                      func=mean_robust
    elif stat=='median_robust':                    func=median_robust
    elif stat=='product':                          func=product
    elif stat=='sigma_mixture':                    func=sigma_mixture
    elif stat=='n_outliers':                       func=n_outliers
    elif stat=='gt0p02':                           func=gt0p02
    elif stat=='sum':                              func=sum
    
    for i in range(n_bins):
        if xmed:
            newx=median(x[i*n_points:(i+1)*n_points])
        else:
            newx=mean(x[i*n_points:(i+1)*n_points])
        if not isfinite(newx): continue
        xb.append(newx)
	if func==std and n_points==2:
	    print 'n_points==2; too few points to determine rms'
	    print 'Returning abs(y1-y2)/2. in bin as rms'
	    yb.append(abs(y[i*n_points]-y[i*n_points+1])/2.)
	else:
	    yb.append(func(y[i*n_points:(i+1)*n_points]))
	if i>2 and xb[-1]==xb[-2]: 
	    yb[-2]=(yb[-2]+yb[-1])/2.
	    xb=xb[:-1]
	    yb=yb[:-1]
    return array(xb),array(yb)
	
def gt0p02(x): return greater(abs(x),0.02).sum()


class stat_robust:
    #Generates robust statistics using a sigma clipping
    #algorithm. It is controlled by the parameters n_sigma
    #and n, the number of iterations

    def __init__(self,x,n_sigma=3,n=5,reject_fraction=None):
        self.x=x
        self.n_sigma=n_sigma
        self.n=n
        self.reject_fraction=reject_fraction
    def run(self):
        good=ones(len(self.x))
        nx=sum(good)
        if self.reject_fraction==None:
            for i in range(self.n):
                if i>0: xs=compress(good,self.x)
                else: xs=self.x
                #            aver=mean(xs)
                aver=median(xs)
                rms=std(xs)
                good=good*less_equal(abs(self.x-aver),self.n_sigma*rms)
                nnx=sum(good)
                if nnx==nx: break
                else: nx=nnx
        else:
            np=float(len(self.x))
            nmin=int((0.5*self.reject_fraction)*np)
            nmax=int((1.-0.5*self.reject_fraction)*np)
            orden=argsort(self.x)
            connect(arange(len(self.x)),sort(self.x))
            good=greater(orden,nmin)*less(orden,nmax)
        self.good=array(good,'int')
        self.remaining=compress(good,self.x)
        self.max=max(self.remaining)
        self.min=min(self.remaining)
        self.mean=mean(self.remaining)
        self.rms=std(self.remaining)
        self.median=median(self.remaining)
        self.outliers=compress(logical_not(good),self.x)
        self.n_remaining=len(self.remaining)
        self.n_outliers=len(self.outliers)
        self.fraction=1.-(float(self.n_remaining)/float(len(self.x)))   


class stat_mad:
    #Generates robust statistics using MAD
    def __init__(self,x,nsigma=3.):
        self.x=x
        self.nsigma=3.
        self.rms=std_mad(self.x)
        self.midpt=median(self.x)
        good=less(abs(self.x-self.midpt),nsigma*self.rms)
        self.good=array(good,'int')
        self.remaining=compress(good,self.x)
        self.max=max(self.remaining)
        self.min=min(self.remaining)
        self.mean=mean(self.remaining)
        self.rms=std_mad(self.remaining)
        self.median=median(self.remaining)
        self.outliers=compress(logical_not(good),self.x)
        self.n_remaining=len(self.remaining)
        self.n_outliers=len(self.outliers)
        self.fraction=1.-(float(self.n_remaining)/float(len(self.x)))
    
#def std_robust(x,n_sigma=3.,n=5):
#    x=purge_outliers(x,n_sigma,n)
#    return std(x-mean(x))

#def mean_robust(x,n_sigma=3.,n=5):
#    x=purge_outliers(x,n_sigma,n)
#    return mean(x)

#def median_robust(x,n_sigma=3.,n=5):
#    x=purge_outliers(x,n_sigma,n)
#    return median(x)



def purge_outliers(x,n_sigma=5.,n=5):
    #Experimental yet. Only 1 dimension
    for i in arange(n):
	med=median(x)
	rms=std_mad(x)
	x=compress(less_equal(abs(x-med),n_sigma*rms),x)
    return x        

def std_robust(x,n_sigma=3.,n=5):
    x=purge_outliers(x,n_sigma,n)
    return std(x-median(x))

def n_outliers(x,n_sigma=5.,n=5):
    n0=len(x)
    for i in arange(n):
	med=median(x)
	rms=std_mad(x)
	x=compress(less_equal(abs(x-med),n_sigma*rms),x)
    return n0-len(x)

def mean_robust(x,n_sigma=3.,n=5):
    x=purge_outliers(x,n_sigma,n)
    return mean(x)

def median_robust(x,n_sigma=3.,n=5):
    x=purge_outliers(x,n_sigma,n)
    return median(x)


def std_log(x,fa=sqrt(20.)):
    #This function does not seem to do anything really useful
    dx=std(x)
    print "std(x)",dx,
    #if abs(dx)<1e-100:dx=mean(abs(x))
    a=fa*dx
    print sqrt(average(a*a*log(1.+x*x/(a*a)))),
    print std_robust(x,3,3),
    print len(x)
    return sqrt(average(a*a*log(1.+x*x/(a*a))))


#def std_log(x,fa=20.):
#    dx=median(abs(x))
#    if abs(dx)<1e-100:dx=mean(abs(x))
#    a=fa*dx
#    return sqrt(average(a*a*log10(1.+x*x/(a*a))))

def med_thr(x,thr=0.2,max_it=10):
    xm=median(x)
    xm0=xm+thr
    for i in range(max_it):
        good=less_equal(x-xm,thr)*greater_equal(x-xm,-thr)
        xm=median(compress(good,x))
        if abs(xm-xm0)<thr/1000.: break
        xm0=xm
        # print xm
    return xm

def std_thr(x,thr=0.2,max_it=10):
    xm=med_thr(x,thr,max_it)
    good=less_equal(x-xm,thr)*greater_equal(x-xm,-thr)
    return std(compress(good,x))

def out_thr(x,thr=0.2,max_it=10):
    xm=med_thr(x,thr,max_it)
    good=less_equal(x-xm,thr)*greater_equal(x-xm,-thr)
    return len(x)-sum(good)


#def bin_aver(x,y,xbins):
#    """Given the variable y=f(x), and 
#    the bins limits xbins, return the 
#    average <y(xbins)>"""
#    a=argsort(x)
#    nbins=len(xbins)
#    y=take(y,a)
#    x=take(x,a)
#    n=searchsorted(x,xbins)
#    results=xbins*0.
#    num=hist(x,xbins)
#    for i in range(nbins):
#	if i< nbins-1:
#	    results[i]=add.reduce(y[n[i]:n[i+1]])
#	else:
#	    results[i]=add.reduce(y[n[i]:])
#	if num[i]>0:
#	    results[i]=results[i]/num[i]
#    return results

def multicompress(condition,variables):
    lista=list(variables)
    n=len(lista)
    for i in range(n):	lista[i]=compress(condition,lista[i])
    return tuple(lista)

def multisort(first,followers):
    #sorts the vector first and matches the ordering 
    # of followers to it
    #Usage:
    # new_followers=multi_sort(first,followers)
    order=argsort(first)
    if type(followers)<> type((1,)):
	return take(followers,order)
    else:
	nvectors=len(followers)
	lista=[]
	for i in range(nvectors):
	    lista.append(take(followers[i],order))
	return tuple(lista)
	
def erfc(x):
    """
    Returns the complementary error function erfc(x)
    erfc(x)=1-erf(x)=2/sqrt(pi)*\int_x^\inf e^-t^2 dt   
    """
    try: x.shape
    except: x=array([x])
    z=abs(x)
    t=1./(1.+0.5*z)
    erfcc=t*exp(-z*z-
        1.26551223+t*(
        1.00002368+t*(
        0.37409196+t*(
        0.09678418+t*(
       -0.18628806+t*(
        0.27886807+t*(
       -1.13520398+t*(
        1.48851587+t*(
       -0.82215223+t*0.17087277)
        ))))))))
    erfcc=where(less(x,0.),2.-erfcc,erfcc)
    return erfcc

def erf(x):
    """
    Returns the error function erf(x)
    erf(x)=2/sqrt(pi)\int_0^x \int e^-t^2 dt
    """
    return 1.-erfc(x)

def cumulative_normal(x):
    return 0.5*(1.+erf(x/sqrt(2.)))

def confidence_interval(x):
    #Confidence interval for x sigma
    return 2.*cumulative_normal(x)-1.

#def erf_brute(x):
#    step=0.00001
#    t=arange(0.,x+step,step)
#    f=2./sqrt(pi)*exp(-t*t)
#    return sum(f)*step

#def erfc_brute(x):
#    return 1.-erf_brute(x)


#def gauss_int_brute(x=arange(0.,3.,.01),average=0.,sigma=1.):
#    step=x[1]-x[0]
#    gn=1./sqrt(2.*pi)/sigma*exp(-(x-average)**2/2./sigma**2)
#    return add.accumulate(gn)*step

def gauss_int_erf(x=(0.,1.),average=0.,sigma=1.):
    """
    Returns integral (x) of p=int_{-x1}^{+x} 1/sqrt(2 pi)/sigma exp(-(t-a)/2sigma^2) dt
    """
    x=(x-average)/sqrt(2.)/sigma
    return (erf(x)-erf(x[0]))*.5

gauss_int=gauss_int_erf

def inv_gauss_int(p):
    #Brute force approach. Limited accuracy for >3sigma
    #find something better 
    #DO NOT USE IN LOOPS (very slow)
    """
    Calculates the x sigma value corresponding to p
    p=int_{-x}^{+x} g(x) dx
    """
    if p<0. or p>1.:
	print 'Wrong value for p(',p,')!'
	sys.exit()
    step=.00001
    xn=arange(0.,4.+step,step)
    gn=1./sqrt(2.*pi)*exp(-xn**2/2.)
    cgn=add.accumulate(gn)*step
    p=p/2.
    ind=searchsorted(cgn,p)
    return xn[ind]


class NumberCounts:
    #Define number counts and produce some plots
    def __init__(self,m,dm=1.,mmin=10.,mmax=35.,area=1.,xcor=None,ycor=None,type_cor='negative'):
        #xcor and ycor are corrections to the total number counts, e.g. area or incompleteness
        if mmin==10. and mmax==35.:
            xm=arange(10.,35.,dm)
            imin=searchsorted(xm,min(m))
            imax=searchsorted(xm,max(m))
            self.xm=xm[imin-1:imax]
            print 'min(m),max(m)',min(m),max(m)
            print 'self.xm[0],self.xm[-1]',self.xm[0],self.xm[-1]
        else:
            self.xm=arange(mmin,mmax+dm,dm)
            
        self.dnc=hist(m,self.xm)
        self.xm=self.xm+dm*0.5
        self.dnc=self.dnc/area/dm
        if xcor<>None and ycor<>None:
            if type_cor=='negative':
                self.dnc-=match_resol(xcor,ycor,self.xm)
                self.dnc=clip(self.dnc,0.,1e50)
            elif type_cor=='positive': self.dnc+=match_resol(xcor,ycor,self.xm)
            elif type_cor=='multiplicative': self.dnc*=match_resol(xcor,ycor,self.xm)

        self.cnc=add.accumulate(self.dnc)
        try:
            self.ldnc=log10(self.dnc)
        except:
            print 'Differential numbers counts contains bins with zero galaxies'
            print 'We set those values to 1e-1'
            dnc=where(equal(self.dnc,0.),1e-2,self.dnc)
            self.ldnc=log10(dnc)
                        
        try:
            self.lcnc=log10(self.cnc)
        except:
            print 'Could not calculate log of cumulative numbers counts'
            
class lsq:
    #Defines a least squares minimum estimator given two 
    #vectors x and y
    def __init__(self,x,y,dy=0.):
	try: dy.shape
	except: dy=x*0.+1.
	dy2=dy**2
	s=add.reduce(1./dy2)
	sx=add.reduce(x/dy2)
	sy=add.reduce(y/dy2)
	sxx=add.reduce(x*x/dy2)
	sxy=add.reduce(x*y/dy2)
	delta=s*sxx-sx*sx
	self.a=(sxx*sy-sx*sxy)/delta
	self.b=(s*sxy-sx*sy)/delta
	self.da=sqrt(sxx/delta)
	self.db=sqrt(s/delta)
    def fit(self,x):
	return self.b*x+self.a

def rotation(x,y,angle):
    xp=x*cos(angle)+y*sin(angle)
    yp=-x*sin(angle)+y*cos(angle)
    return xp,yp


class Parameter:
#http://www.scipy.org/Cookbook/FittingData?highlight=(fitting)
    def __init__(self, value):
        self.value = value
        
    def set(self, value):
        self.value = value
   
    def __call__(self):
        return self.value
   
def fit(function, parameters, y, x = None):
# http://www.scipy.org/Cookbook/FittingData?highlight=(fitting)
# mu = Parameter(7)
# sigma = Parameter(3)
# height = Parameter(5)
# def f(x): return height() * exp(-((x-mu())/sigma())**2)
# fit(f, [mu, sigma, height], data)
    def f(params):
        i = 0
        for p in parameters:
            p.set(params[i])
            i += 1
        return y - function(x)
   
    if x is None: x = arange(y.shape[0])
    p = [param() for param in parameters]
    return optimize.leastsq(f, p)


def savitzky_golay(data, kernel = 11, order = 4):
    #From scipy/Cookbook
    """
        applies a Savitzky-Golay filter
        input parameters:
        - data => data as a 1D numpy array
        - kernel => a positiv integer > 2*order giving the kernel size
        - order => order of the polynomal
        returns smoothed data as a numpy array

        invoke like:
        smoothed = savitzky_golay(<rough>, [kernel = value], [order = value]
    """
    try:
            kernel = abs(int(kernel))
            order = abs(int(order))
    except ValueError, msg:
        raise ValueError("kernel and order have to be of type int (floats will be converted).")
    if kernel % 2 != 1 or kernel < 1:
        raise TypeError("kernel size must be a positive odd number, was: %d" % kernel)
    if kernel < order + 2:
        raise TypeError("kernel is to small for the polynomals\nshould be > order + 2")

    # a second order polynomal has 3 coefficients
    order_range = range(order+1)
    half_window = (kernel -1) // 2
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    # since we don't want the derivative, else choose [1] or [2], respectively
    m = np.linalg.pinv(b).A[0]
    window_size = len(m)
    half_window = (window_size-1) // 2

    # precompute the offset values for better performance
    offsets = range(-half_window, half_window+1)
    offset_data = zip(offsets, m)

    smooth_data = list()

    # temporary data, with padded zeros (since we want the same length after smoothing)
    data = np.concatenate((np.zeros(half_window), data, np.zeros(half_window)))
    for i in range(half_window, len(data) - half_window):
            value = 0.0
            for offset, weight in offset_data:
               value += weight * data[i + offset]
            smooth_data.append(value)
    return np.array(smooth_data)

def sgolay2d( z, window_size, order, derivative=None):
    """
    """
    # number of terms in the polynomial expression
    n_terms = ( order + 1 ) * ( order + 2)  / 2.0

    if  window_size % 2 == 0:
        raise ValueError('window_size must be odd')

    if window_size**2 < n_terms:
        raise ValueError('order is too high for the window size')

    half_size = window_size // 2

    # exponents of the polynomial. 
    # p(x,y) = a0 + a1*x + a2*y + a3*x^2 + a4*y^2 + a5*x*y + ... 
    # this line gives a list of two item tuple. Each tuple contains 
    # the exponents of the k-th term. First element of tuple is for x
    # second element for y.
    # Ex. exps = [(0,0), (1,0), (0,1), (2,0), (1,1), (0,2), ...]
    exps = [ (k-n, n) for k in range(order+1) for n in range(k+1) ]

    # coordinates of points
    ind = np.arange(-half_size, half_size+1, dtype=np.float64)
    dx = np.repeat( ind, window_size )
    dy = np.tile( ind, [window_size, 1]).reshape(window_size**2, )

    # build matrix of system of equation
    A = np.empty( (window_size**2, len(exps)) )
    for i, exp in enumerate( exps ):
        A[:,i] = (dx**exp[0]) * (dy**exp[1])

    # pad input array with appropriate values at the four borders
    new_shape = z.shape[0] + 2*half_size, z.shape[1] + 2*half_size
    Z = np.zeros( (new_shape) )
    # top band
    band = z[0, :]
    Z[:half_size, half_size:-half_size] =  band -  np.abs( np.flipud( z[1:half_size+1, :] ) - band )
    # bottom band
    band = z[-1, :]
    Z[-half_size:, half_size:-half_size] = band  + np.abs( np.flipud( z[-half_size-1:-1, :] )  -band )
    # left band
    band = np.tile( z[:,0].reshape(-1,1), [1,half_size])
    Z[half_size:-half_size, :half_size] = band - np.abs( np.fliplr( z[:, 1:half_size+1] ) - band )
    # right band
    band = np.tile( z[:,-1].reshape(-1,1), [1,half_size] )
    Z[half_size:-half_size, -half_size:] =  band + np.abs( np.fliplr( z[:, -half_size-1:-1] ) - band )
    # central band
    Z[half_size:-half_size, half_size:-half_size] = z

    # top left corner
    band = z[0,0]
    Z[:half_size,:half_size] = band - np.abs( np.flipud(np.fliplr(z[1:half_size+1,1:half_size+1]) ) - band )
    # bottom right corner
    band = z[-1,-1]
    Z[-half_size:,-half_size:] = band + np.abs( np.flipud(np.fliplr(z[-half_size-1:-1,-half_size-1:-1]) ) - band )

    # top right corner
    band = Z[half_size,-half_size:]
    Z[:half_size,-half_size:] = band - np.abs( np.flipud(Z[half_size+1:2*half_size+1,-half_size:]) - band )
    # bottom left corner
    band = Z[-half_size:,half_size].reshape(-1,1)
    Z[-half_size:,:half_size] = band - np.abs( np.fliplr(Z[-half_size:, half_size+1:2*half_size+1]) - band )

    # solve system and convolve
    if derivative == None:
        m = np.linalg.pinv(A)[0].reshape((window_size, -1))
        return scipy.signal.fftconvolve(Z, m, mode='valid')
    elif derivative == 'col':
        c = np.linalg.pinv(A)[1].reshape((window_size, -1))
        return scipy.signal.fftconvolve(Z, -c, mode='valid')
    elif derivative == 'row':
        r = np.linalg.pinv(A)[2].reshape((window_size, -1))
        return scipy.signal.fftconvolve(Z, -r, mode='valid')
    elif derivative == 'both':
        c = np.linalg.pinv(A)[1].reshape((window_size, -1))
        r = np.linalg.pinv(A)[2].reshape((window_size, -1))
        return scipy.signal.fftconvolve(Z, -r, mode='valid'), scipy.signal.fftconvolve(Z, -c, mode='valid')


#Tests 

def Testing(test):
    print 'Testing ',test,'...'

def test():
    """ Tests the functions defined in this module"""

    Testing("I/O FUNCTIONS")
    
    test="put_str and get_str"
    x=arange(100.)
    y=map(str,sin(x))
    x=map(str,x)
    put_str('test.txt',(x,y))
    xn,yn=get_str('test.txt',range(2))
    if xn<>x:raise test
    if yn<>y:raise test

    test='put_data and get_data'
    Testing(test)
    x=arange(100.)
    y=sin(x)
    put_data('test.dat',(x,y),'header of test.dat','%.18f %.18f')
    xn,yn=get_data('test.dat',range(2))
    if sometrue(not_equal(xn,x)): raise test
    if sometrue(not_equal(yn,y)): raise test

    test="put_header and get_header"
    Testing(test)
    f=open('test.head','w')
    f.close()

    lines=["This is a test header line",
	   "This is another test header line",
	   "This is yet another test header line"]

    put_header('test.head',lines[2])
    put_header('test.head',lines[1])
    put_header('test.head',lines[0])

    nlines=get_header("test.head")
    lines=string.join(map(lambda x: '#'+x+'\n',lines),'')

    if lines<>nlines: 
	print `lines`
	print `nlines`
	raise test

    test='put_2Darray and get_2Darray'
    Testing(test)
    x=arange(200.)
    y=reshape(x,(40,-1))
    put_2Darray('test.dat',y,'header of test.dat')
    yn=get_2Darray('test.dat')
    comp=not_equal(yn,y)
    for i in range(yn.shape[1]):
	if sometrue(comp[:,i]): raise test
    
    #Testing("MISC NUMERICAL FUNCTIONS")

    test='ascend'
    Testing(test)
    y=sin(arange(100))
    z=sort(y)
    if not ascend(z): raise test
    z=z[::-1]
    if ascend(z): raise test

#    test="hist"
#    Testing(test)
#    x=arange(0.,100.,.1)
#    y=arange(100.)
#    h=hist(x,y)
#    points(h,ones(100)*10)
#    if sometrue(not_equal(h,ones(100)*10)): raise test

    test="bin_aver"
    Testing(test)
    x=arange(0.,10.1,.1)
    xb=arange(10.)+.01
    y=x
    yb=bin_aver(x,y,xb)+.45
    yr=arange(10)+1.
    if sometrue(greater_equal(yb-yr,1e-10)): raise test
    
    test='dist'
    Testing(test)
    a=arange(0,2.*pi,pi/6.)
    x=10.*sin(a)+3.
    y=5*cos(a)-3.
    d=sqrt((((x-3.)**2+(y+3.)**2)))
    nd=map(dist,x,y,ones(len(x))*3.,ones(len(x))*-3.)
    if sometrue(not_equal(d,nd)):
	print d
	print nd
	raise test


    test="loc2d"
    Testing(test)
    m=fromfunction(dist,(10,10))
    if loc2d(m)<>(9,9): raise test
    if loc2d(m,'min')<>(0,0): raise test


    test="match_objects"
    Testing(test)
    x=arange(10.)
    y=arange(10,20.)
    t1=(map(str,x),map(str,y))
    x=x+np.random(x.shape)*.707/2.
    y=y+RandomArray.random(y.shape)*.707/2.
    x0=arange(10000.)
    y0=arange(10.,10010.)
    t2=(map(str,x0),map(str,y0))
    cosas1=match_objects((x,y),(x0,y0),t1,t2,accuracy=.5)
    if not (cosas1[2]==cosas1[4] and cosas1[3]==cosas1[5]):
	raise test

    test="match_min"
    Testing(test)
    x=arange(10.)
    y=arange(10,20.)
    t1=(map(str,x),map(str,y))
    x=x+RandomArray.random(x.shape)*.707/2.
    y=y+RandomArray.random(y.shape)*.707/2.
    x0=arange(10000.)
    y0=arange(10.,10010.)
    t2=(map(str,x0),map(str,y0))
    cosas1=match_min((x,y),(x0,y0),t1,t2)
    #put_data('bobo',cosas1)
    #os.system('more bobo')
    if not (cosas1[2]==cosas1[4] and cosas1[3]==cosas1[5]):
	raise test


    test="match_resol"
    Testing(test)
    xobs=arange(0.,10.,.33)
    yobs=cos(xobs)*exp(-xobs)
    xt=arange(0.,10.,.33/2.)
    yt=match_resol(xobs,yobs,xt)
    ytobs=cos(xt)*exp(-xt)
    if plots:
	plot=FramedPlot()
	plot.add(Points(xobs,yobs,color="blue",type='cross'))
	plot.add(Curve(xt,yt,color="red"))
	plot.add(Points(xt,yt,color="red",type='square'))
	plot.show()
	print "The crosses are the original data"
	print "The continuous line/red squares represents the interpolation"
    else:
	print '   X     Y_Interp  Y_expected'
	for i in range(len(x)):
            print 3*'%7.4f  '% (xt[i],yt[i],ytobs[i])

    test="gauss_int"
    Testing(test)
    x=arange(0.,3,.5)
    p=gauss_int(x)
    pt=array(
	[0.,.19146,.34134,.43319,.47725,.49379])
    diff=abs(p-pt)
    print diff
    if sometrue(greater(diff,2e-5)): raise test

    test="inv_gauss_int"
    Testing(test)
    for i in range(len(pt)):
	z=inv_gauss_int(2.*pt[i])
	if abs(x[i]-z) > 1e-3 : raise test

    print 'Everything tested seems to be OK in useful.py'


if __name__ == '__main__':test()
else: pass


