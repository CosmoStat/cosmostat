'''
Created on Mar 30, 2015

@author: Ming  Jiang and Jean-Luc Starck

                    CLASS FUNCTION  class starlet2d()
            
    Allow to perform a starlet transform, manipulate it (visualisation, thresholding, statistics, etc),
    and to reconstruct.
    If pysap is installed, then the pysparse module should be available and 
    the code will used C++ binding for fast calculation.
    Otherwise full python code is used.
    
    Details of the starlet transform can be found in 
    J.L. Starck, F. Murtagh, and J. Fadili, 
    Sparse Image and Signal Processing: Wavelets and 
    Related Geometric Multiscale Analysis, 
    Cambridge University Press, Cambridge (GB),  2016.
    
    or 
    
    J.-L. Starck, J. Fadili and F. Murtagh,
    "The Undecimated Wavelet Decomposition and its Reconstruction", 
    IEEE Transaction on Signal Processing , 16, 2, pp 297--309, 2007.
        
    Example how to use the Class:
        CW = starlet2d()      # Create the class
        CW.transform(Image)   # Starlet transform of a 2D np array
        CW.stat()             # print statistics of all scales
        r = CW.recons()       #  reconstruct an image from its coefficients
    more examples are given at the end of this file.
'''
import numpy as np
import scipy.signal as psg
# import pcosmostat.sparsity.sparse2d.param as pm
from pycs.tools.cosmostat_init import *
from pycs.tools.stats import *

import sys
 
import imp
PYSAP_CXX = True
try:
	import pysparse
#    imp.find_module('pysparse')
except ImportError:
    PYSAP_CXX=False

#if 'pysparse' in sys.modules:
#    import pysparse
#    PYSAP_CXX = True

if PYSAP_CXX is False:
    print("Warning in starlet.py: do not find pysad bindings ==> use slow python code. ")

#print("PYSAP_CXX = ", PYSAP_CXX)


def test_ind(ind,N):
    """
    function to handle the border using a mirror effect.
    If the index is < 0 or >= N, where N is the size of image in one direction,
    it returns the correct index in [0,N-1], using mirror effect.
    Parameters
    ----------
    ind : TYPE
        DESCRIPTION.
    N : TYPE
        DESCRIPTION.

    Returns
    -------
    res : TYPE
        DESCRIPTION.

    """
    res = ind
    if ind < 0 :
        res = -ind
        if res >= N:
            res = 2*N - 2 - ind
    if ind >= N :
        res = 2*N - 2 - ind
        if res < 0:
            res = -ind
    return res


def b3splineTrans(im_in,step):
    """
    Apply a 2d B-spline smmothing to an image, using holes in the smoothing 
    kernel (a-trous algorithm)
    Parameters
    ----------
    im_in : np.ndarray
            input image.
    step : int
        the hole size.



    Returns
    -------
    im_out : 2D np.ndarray
        smoothed image.
    """
    (nx,ny) = np.shape(im_in)
    im_out = np.zeros((nx,ny))
    c1 = 1./16
    c2 = 1./4
    c3 = 3./8

    buff = np.zeros((nx,ny))

    for i in np.arange(nx):
        for j in np.arange(ny):
            jl = test_ind(j-step,ny)
            jr = test_ind(j+step,ny)
            jl2 = test_ind(j-2*step,ny)
            jr2 = test_ind(j+2*step,ny)
            buff[i,j] = c3 * im_in[i,j] + c2 * (im_in[i,jl] + im_in[i,jr]) + c1 * (im_in[i,jl2] + im_in[i,jr2])

    for j in np.arange(ny):
        for i in np.arange(nx):
            il = test_ind(i-step,nx)
            ir = test_ind(i+step,nx)
            il2 = test_ind(i-2*step,nx)
            ir2 = test_ind(i+2*step,nx)
            im_out[i,j] = c3 * buff[i,j] + c2 * (buff[il,j] + buff[ir,j]) + c1 * (buff[il2,j] + buff[ir2,j])

    return im_out

def b3spline_fast(step):
    """
    Kernel computation for fast smoothing using convolve2d function
    Parameters
    ----------
    step : TYPE
        the hole size.
    Returns
    -------
    kernel2d : 2D np.ndarray
        calculated kernel.
    """
    step_hole = int(step)
    c1 = 1./16.
    c2 = 1./4.
    c3 = 3./8.
    length = int(4*step_hole+1)
    kernel1d =  np.zeros((1,length))
    kernel1d[0,0] = c1
    kernel1d[0,-1] = c1
    kernel1d[0,step_hole] = c2
    kernel1d[0,-1-step_hole] = c2
    kernel1d[0,2*step_hole] = c3
    kernel2d = np.dot(kernel1d.T,kernel1d)
    return kernel2d

def star2d(im, scale, gen2=False, bord=1, nb_procs=0, fast=True, verb=0):
    """
    Routie to calculate the 1st and 2nd generation starlet transform.
    if the global variable PYSAP_CXX is True, a C++ code will be used through
    binding for this calculation.
    Parameters
    ----------
    im : 2D np.ndarray
            input image.
     scale : int. 
        number of scales.
    gen2 : bool, optional
        if True, performs the second generation starlet transform
    bord : int, optional
        Type of border used to handle the border effect. The default is 1.
        this parameter is only used if the C++ pysap code is available.
    nb_procs : int, optional
        Numper of preocessor to use. Only used if the C++ pysap code is available 
        and if openmp is available. The default is 0.
    fast : bool, optional
        for python implementation only. If true, the convolve2d routine is used, 
        which is faster. The default is True.
    verb : bool, optional
        Verbose mode. The default is 0.
    Returns
    -------
    3D np.ndarray
            output wavelet transform  [0:scale,0:nx,0:ny]
    """
#    print ('IN STAR2D 2')
    (nx,ny) = np.shape(im)
    nz = scale
    # Normalized transfromation
    if PYSAP_CXX is True:
        # print("BINDING: ", head, ", norm = ", l2norm)
        # verb=1
        ima = np.zeros((nx,ny))
        ima[:,:]=im
        psWT = pysparse.MRStarlet(bord, gen2, nb_procs,verb)
        wl = psWT.transform(ima.astype(np.float),nz)
        wt = (np.stack(wl)).astype(np.double)
    else:
        wt = np.zeros((nz,nx,ny))
        step_hole = int(1)
        im_in = np.copy(im)
        for i in np.arange(nz-1):
            if fast:
                kernel2d = b3spline_fast(step_hole)
                im_out = psg.convolve2d(im_in, kernel2d, boundary='symm',mode='same')
            else:
                im_out = b3splineTrans(im_in,step_hole)
            if gen2:
                if fast:
                    im_aux = psg.convolve2d(im_out, kernel2d, boundary='symm',mode='same')
                else:
                    im_aux = b3splineTrans(im_out,step_hole)
                wt[i,:,:] = im_in - im_aux
            else:
                wt[i,:,:] = im_in - im_out
            im_in = np.copy(im_out)
            step_hole *= 2
        wt[nz-1,:,:] = np.copy(im_out)
    return wt


def istar2d(wt, gen2=True, bord=0, nb_procs=0, fast=True, verb=0):
    """
    Routie to calculate the 1st and 2nd generation incerse starlet transform.
    if the global variable PYSAP_CXX is True, a C++ code will be used through
    binding for this calculation.
    Parameters
    ----------
    wt : 3D np.ndarray
            input wavelet transform.
    gen2 : bool, optional
        if True. assume the second generation starlet reconstruction.
    bord : int, optional
        Type of border used to handle the border effect. The default is 1.
        this parameter is only used if the C++ pysap code is available.
    nb_procs : int, optional
        Numper of preocessor to use. Only used if the C++ pysap code is available 
        and if openmp is available. The default is 0.
    fast : bool, optional
        for python implementation only. If true, the convolve2d routine is used, 
        which is faster. The default is True.
    verb : bool, optional
        Verbose mode. The default is 0.

    Returns
    -------
    2D np.ndarray:
        Reconstructed image
    """
    (nz,nx,ny) = np.shape(wt)

    # PYSAP_CXX=0
    if PYSAP_CXX is True:
        # print("RECBINDING: ", head, ", norm = ", l2norm)
        dat_list = []
        for s in range(nz):    
            dat_list.append( wt[s,:,:].astype(np.float))
        psWT = pysparse.MRStarlet(bord, gen2, nb_procs,verb)
        imRec = (psWT.recons(dat_list)).astype(np.double)
    else:
        # trans = 1 if gen2 else 2
        if gen2:
            '''
            h' = h, g' = Dirac
            '''
            step_hole = int(pow(2,nz-2))
            imRec = np.copy(wt[nz-1,:,:])
            for k in np.arange(nz-2,-1,-1):
                if fast:
                    kernel2d = b3spline_fast(step_hole)
                    im_out = psg.convolve2d(imRec, kernel2d, boundary='symm',mode='same')
                else:
                    im_out = b3splineTrans(imRec,step_hole)
                imRec = im_out + wt[k,:,:]
                step_hole /= 2
        else:
            '''
            h' = Dirac, g' = Dirac
            '''
    #         imRec = np.sum(wt,axis=0)
            '''
            h' = h, g' = Dirac + h
            '''
            imRec = np.copy(wt[nz-1,:,:])
            step_hole = int(pow(2,nz-2))
            for k in np.arange(nz-2,-1,-1):
                if fast:
                    kernel2d = b3spline_fast(step_hole)
                    imRec = psg.convolve2d(imRec, kernel2d, boundary='symm',mode='same')
                    im_out = psg.convolve2d(wt[k,:,:], kernel2d, boundary='symm',mode='same')
                else:
                    imRec = b3splineTrans(imRec,step_hole)
                    im_out = b3splineTrans(wt[k,:,:],step_hole)
                imRec += wt[k,:,:]+im_out
                step_hole /= 2
 
    return imRec    
        

def adstar2d(wtOri, gen2=True, bord=0, nb_procs=0, fast=True, verb=0):
    """
    Routine to calculate the 1st and 2nd generation adjoint starlet operator.
    if the global variable PYSAP_CXX is True, a C++ code will be used through
    binding for this calculation. This routine is generally used when the gradient 
    of a functional involving a starlet transform operator is required.
    Parameters
    ----------
    wtOri : 3D np.ndarray
            input wavelet transform.
    gen2 : bool, optional
        if True. assume the second generation starlet reconstruction.
    bord : int, optional
        Type of border used to handle the border effect. The default is 1.
        this parameter is only used if the C++ pysap code is available.
    nb_procs : int, optional
        Numper of preocessor to use. Only used if the C++ pysap code is available 
        and if openmp is available. The default is 0.
    fast : bool, optional
        for python implementation only. If true, the convolve2d routine is used, 
        which is faster. The default is True.
    verb : bool, optional
        Verbose mode. The default is 0.
    Returns
    -------
    2D np.ndarray:
        Reconstructed image
    """
    (nz,nx,ny) = np.shape(wtOri)
    wt = np.copy(wtOri)
    if PYSAP_CXX is True:
        # print("BINDING")
        dat_list = []
        for s in range(nz):
            dat_list.append((wt[s,:,:]).astype(float))
        psWT = pysparse.MRStarlet(bord, gen2, nb_procs, verb)
        imRec = (psWT.recons(dat_list,True)).astype(double)   
    else:
        # print("NO BINDING")
        # Unnormalization step
        # !Attention: wt is not the original wt after unnormalization
        imRec = np.copy(wt[nz-1,:,:])
        step_hole = pow(2,nz-2)
        for k in np.arange(nz-2,-1,-1):
            if fast:
                kernel2d = b3spline_fast(step_hole)
                imRec = psg.convolve2d(imRec, kernel2d, boundary='symm',mode='same')
                im_out = psg.convolve2d(wt[k,:,:], kernel2d, boundary='symm',mode='same')
                if gen2:
                    im_out2 = psg.convolve2d(im_out, kernel2d, boundary='symm',mode='same')
                    imRec += wt[k,:,:] -im_out2
                else: imRec += wt[k,:,:] -im_out
            else:
                imRec = b3splineTrans(imRec,step_hole)
                im_out = b3splineTrans(wt[k,:,:],step_hole)
                if gen2:
                    im_out2 = b3splineTrans(im_out,step_hole)
                    imRec += wt[k,:,:] -im_out2
                else: imRec += wt[k,:,:]-im_out
            step_hole /= 2
    return imRec


#==========================================================================
#======================= Beginning of the STARLET CLASS ==================
#==========================================================================

class starlet2d():
    """
    Class for the starlet decomposition and reconstruction
    """
    name = "wt"     # name of the class
    gen2 = True     # if true, it will the second genereal starlet transform    
    l2norm=False    # if true, consider a l2 normalisation
    nx=0            # image size first axis
    ny=0            # image size second axis
    ns=0            # number of scales
    coef=0.         # Starlet coefficients
    TabNorm=0.      # Coefficient normalixation table
    SigmaNoise = 1. # noise standard deviation
    TabNsigma = 0   # detection level per scale
    Starlet_Gen1TabNorm =0 # Normalization table for the first generation starlet transform

    # __init__ is the constructor
    def __init__(self, name='wt', gen2=True,l2norm=True, bord=1, verb=False, nb_procs=0): 
        """
        Constructor

        Parameters
        ----------
        name : string, optional
            name of transform. Used when information is printed. The default is 'wt'.
        gen2 : bool, optional
            if True. assume the second generation starlet reconstruction.
        l2norm : bool, optional
            if True, assume a l2 normalisation of the wavelet coefficients.
        bord : int, optional
            Type of border used to handle the border effect. The default is 1.
            this parameter is only used if the C++ pysap code is available.
                    # case 0:  bord = I_ZERO;  
                    # case 1:  bord = I_CONT;  
                    # case 2:  bord = I_MIRROR;  
                    # case 3:  bord = I_PERIOD;
        nb_procs : int, optional
            Numper of preocessor to use. Only used if the C++ pysap code is available 
            and if openmp is available. The default is 0.
        verb : bool, optional
        Returns
        -------
        None.

        """
        self.name = name      # self.name is an object variable
        self.gen2=gen2
        self.l2norm=l2norm
        self.verb=verb
        self.nb_procs=nb_procs
        self.bord=bord

    def get_gen1_starlet_tabnorm(self):
        """
        Compute the normalisation coefficients at each scale of the firast generation 
        starlet transform.    
        Returns
        -------
        tabNs : TYPE
            DESCRIPTION.
        """
        im = np.zeros((self.nx,self.ny))
        im = im.astype('float64')
        im[int(self.nx/2),int(self.ny/2)] = np.float64(1.)
        wt = star2d(im,self.ns,gen2=False)
        tmp = wt**2
        tabNs = np.sqrt(np.sum(np.sum(tmp,1),1))
        return tabNs
  
    def init_starlet(self, nx, ny, nscale=0):
        """
        Initialize the scale for a given image size and a number of scales.
        Parameters
        ----------
        nx, ny : int
            image size.
        nscale : int, optional
            Number of wavelet scales. The default is 0.
            If it is 0, the numnber of scales is fixed to
                  log( MIN([nx,ny]))  
        Returns
        -------
        None.
        """
        self.nx = int(nx)
        self.ny = int(ny)
        if nscale == 0:
            mins = np.min( [nx,ny])
            nscale = int(np.log(mins)  // 1)
        self.ns = int(nscale)
        self.Starlet_Gen1TabNorm = self.get_gen1_starlet_tabnorm()
        if self.l2norm:
            self.TabNorm = np.ones(self.ns, dtype=float)
        else:
            self.TabNorm = self.get_gen1_starlet_tabnorm()
        # for pysparse 
        self.nb_procs=0

    def info(self):          # sound is a method (a method is a function of an object)
        """
        Print information relative to the intialisation.
        """
        print(self.name, ": Nx = ", self.nx, ", Ny = ", self.ny, ", Ns = ", self.ns)
        if self.gen2:
            print("starlet 2nd  generation")
        else:
            print("starlet 1st  generation")
        if self.l2norm:
            print("l2 normalisation")
        else:
            print("l1 normalisation")
        # print("Coef TabSize = ", np.shape(self.coef))

    def stat(self):
        """
        Print Min, Max, Mean and standard deviation of all scales.
        Parameters
        ----------
        None.

        Returns
        -------
        None.

        """
        print(self.name, ": Nx = ", self.nx, ", Ny = ", self.ny, ", Ns = ", self.ns)
        for j in range(self.ns):
            s = (self.coef)[j]
            print("%s Scale %2d: Min = %f, Max = %f, Mean = %f, std = %f" % (self.name, j+1,s.min(), s.max(), s.mean(), s.std()))

#    def transform(im,nscale,gen2=self.gen2,normalization=self.l2norm):
    def transform(self, im, WTname=None):
        """
        Apply the starlet transform to image. Coeffients are stored in 
        self.coef[:,:,:].  self.coef[s,:,:] is the wavelet scale at scale s.
        See class routines get_scale, get_ptr_scale, put_scale to manipulate
        the coefficients.
        Parameters
        ----------
        im : 2D np.ndarray
             input image..
        WTname : string, optional
            Name given to the decomposition. The default is None.
        Returns
        -------
        None. 
        """
        (Nx,Ny) = im.shape
        if self.ns <=1 or self.nx != Nx or self.ny != Ny :
            self.init_starlet(Nx, Ny, nscale=0)
        if WTname is not None: 
            self.name = WTname
        self.coef = star2d(im, self.ns, self.gen2, self.bord, self.nb_procs, True, self.verb)
        if self.l2norm:
            for i in np.arange(self.ns):
                self.coef[i,:,:] /= self.Starlet_Gen1TabNorm[i]


    def recons(self, adjoint=False):
        """
        Reconstruct an image from its calculated starlet coefficients.
        Parameters
        ----------
        adjoint : bool, optional
            If true, used the adjoint operator instead of the exact reconstruction one. 
            The default is False.

        Returns
        -------
        rec : 2D np.ndarray
            Reconstructed image.
        """
        wt = np.copy(self.coef)
        if self.l2norm:
            for i in np.arange(self.ns):
                wt[i,:,:] *= self.Starlet_Gen1TabNorm[i]
        if adjoint:
            rec =  adstar2d(wt,gen2=self.gen2,bord=self.bord, nb_procs=self.nb_procs, fast=True, verb=self.verb)
        else:
            rec =  istar2d(wt, gen2=self.gen2,bord=self.bord, nb_procs=self.nb_procs, fast=True, verb=self.verb)
        return rec

    def denoising(self, Image, SigmaNoise=0, Nsigma=3,ThresCoarse=False, hard=True):
        """
        Do a denoising of the input image, by taking the wavelet decomposition, 
        thresholding it, and reconstructing the denoised image.
        Parameters
        ----------
        Image : 2D np.ndarray
            DESCRIPTION.
        SigmaNoise : float, optional
            Standard deviation of the noise. Default is 0.
        Nsigma: float, optional
            Detection level. Defautl is 3  (.e. 3 SigmaNoise).
        ThresCoarse : bool, optional
            IF true the coarsest scale is removed. The default is False.
        hard : bool, optional
            Type of threshold, true for hard thresholding and false
            for soft thresholding. The default is True.
        Returns
        -------
        2D np.ndarray
            Denoised image.
        """
        if SigmaNoise == 0:
            SigmaNoise = get_noise(Image)
        self.SigmaNoise = SigmaNoise
        self.transform(Image)
        self.threshold(SigmaNoise=SigmaNoise, Nsigma=Nsigma, ThresCoarse=ThresCoarse, hard=hard)
        return self.recons()


    def pos_transform(self,im, nscale=0, Niter=100,fast=True,hard=False,den=False, KillCoarse=False, pos=True, SigmaNoise=0, Nsigma=3.,verb=False):
        """
        Iterative method to make a decomposition on positive coefficients.
        Coeffients are stored in self.coef[:,:,:].  
        See class routines get_scale, get_ptr_scale, put_scale to manipulate
        the coefficients.
        
        Parameters
        ----------
            ----------
        im : 2D np.ndarray
            input image.
        hard : bool, optional
            if True, use hard thresholding, and soft thresholding otherwise.
            Default is False
        den : bool, optional
            if true, denoise also the coefficeints. Default is False
        KillCoarse : bool, optional
            IF true the coarsest scale is removed. The default is False.
        fast : bool, optional
            for python implementation only. If true, the convolve2d routine is used, 
            which is faster. The default is True.
        pos: bool, optional
            it true, keep only positive wavelet coefficients. Default is True. 
        SigmaNoise: float, optional
            Standard deviation of the noise. Default is 0.
        Nsigma: float, optional
            Detection level. Defautl is 3  (.e. 3 SigmaNoise).  
        verb : bool, optional
            Verbose mode. The default is 0.
        Raises
        ------
        ValueError
            Can only be used if the number of scales > 1

        Returns
        -------
        None.

        """
        self.l2norm=True
        (Nx,Ny) = im.shape
        self.init_starlet(Nx, Ny, nscale=self.ns)
        if self.ns <= 1:
            raise ValueError('Number of scales must be > 1 ! '
                             'Input value = {} and is of type {}.'.format(nscale, type(nscale)))

        rsd = np.copy(im)
        self.transform(im)
        mwt = self.coef.max()
        # wt = np.copy(self.coef)  
        wt = np.zeros((self.ns,self.nx,self.ny))
        for it in np.arange(Niter):
            ld = mwt * (1. - (it+1.)/Niter)
            if ld < 0:
                ld = 0
            if verb:
                print ("Iter ", it, ": lambda="+str(ld), ", Resi = ", np.std(rsd))
            self.transform(rsd)
            wt += self.coef
            if den:
                if SigmaNoise != 0:
                    noise = mad(wt[0])
                else:
                    noise = SigmaNoise
    #            print(noise)
                hard_thresholding(wt,Nsigma*noise)
            if hard:
                hard_thresholding(wt,ld)
            else:
                soft_thresholding(wt,ld)
            if pos is True:
                wt[wt<0] = 0
            if KillCoarse is True:
                    wt[self.ns-1,:,:] = 0
            self.coef = np.copy(wt) 
            rec = self.recons()  
    #        print (rec>=0).all()
            rsd = im - rec
    #         fits.writeto('pstar2d'+str(it)+'.fits',rsd,clobber=True)
    #        print ((np.abs(rsd)).sum())


    def get_scale(self, j):
        """
        Return a copy of a given scale of the decomposition. 
        Parameters
        ----------
        j : int
            Scale number. It must be in [0:self.ns]
        Returns
        -------
        Scale : 2D np.ndarray
            jth wavelet scale of the decomposition.
        """
        Scale = np.zeros((self.nx,self.ny))
        Scale[:,:]=(self.coef)[j,:,:]
        return Scale

    def get_ptr_scale(self, j):
        """
        Return a pointer to the jth scale. Modifying the return array will 
        impact the coefficients self.coef of the class.

        Parameters
        ----------
        j : int
            Scale number. It must be in [0:self.ns]
        Returns
        -------
        Scale : 2D np.ndarray
            jth wavelet scale of the decomposition.
        """
        return (self.coef)[j]
        
    def put_scale(self, ScaleCoef, j):
        """
        Replace the scale j in self.coef by the 2D array ScaleCoef.
        Parameters
        ----------
        ScaleCoef : 2D np.ndarray
            New coefficients at scale j to be inserted in the class.
        j : int
             Scale number. It must be in [0:self.ns].
        Returns
        -------
        None.

        """
        self.coef[j,:,:] = ScaleCoef

    def tvs(self, j):
        """
        Display the scale j
        Parameters
        ----------
        j : int
             Scale number. It must be in [0:self.ns].
        Returns
        -------
        Window appearing showing scale j.
        """
        s = self.get_ptr_scale(j)
        tvilut(s)

    def dump(self):
        """
        Print all variable and function names of the class
        Returns
        -------
        None.

        """
        print(self.__dict__)

    def get_noise(self):
        """
        Estimate the noise in the data from the first wavelet scale
        Returns
        -------
        SigmaNoise : float
            estimated noise standard deviation.
        """
        s = (self.coef)[0]
        SigmaNoise = mad(s)
        return SigmaNoise
        
    def tvsl(self, j, SigmaNoise=0, Levels=[5]):
        """
        Display the scale j, with contours corresponding to the noise detect
        levels given in Levels. Several contour levels can be given.
        Parameters
        ----------
        j : int
             Scale number. It must be in [0:self.ns].
        Returns
        -------
        Window appearing showing scale j, with contours around structures 
        detected a specified levels.
        """
        if SigmaNoise == 0:
            SigmaNoise = self.get_noise()
#        print(self.TabNorm)
#        print("noise ", SigmaNoise)
        Norm = self.TabNorm[j]
#        print("norm = ", Norm)
        TabLevels = np.array(Levels)
#        print("TabLevels = ", TabLevels)
        TabLevels = self.TabNorm[j] * SigmaNoise * np.array(Levels)
#        print(TabLevels)
        TabLevels[ TabLevels > (self.coef)[j].max()] = (self.coef)[j].max()
        tvimacont((self.coef)[j], TabLevels, vmin=0, vmax=0, gamma=0.5, cmap='gist_stern')
        
    def tvall(self, scales=None, multiview=False):
        """
        Display a window with all scales.
        Parameters
        ----------
        scales : list, optional
            selection of scales. The default is None.
        multiview : int, optional
            multiview. The default is False.

        Returns
        -------
        None.
        """
        tv_frames(self.coef, scales=None, multiview=False)

    def get_tabsigma(self, nscale, Nsigma=3):
        """
        Create the detection table TabNsigma[0:nsale-1], for diffent type
        of calling. 
        By default, it is 4 at the finest scale at 3 at the others.
        If Nsigma is an array small than the number of scales, the last value
        of Nsigma is repeated.
        exemple of call:
            print(CLASS.get_tabsigma(4))              => array([4., 3., 3., 3.])
            print(CLASS.get_tabsigma(4, Nsigma=[3,4]) => array([3, 4., 4., 4.])

        Parameters
        ----------
        nscale : int
            number of scales.
        Nsigma : int or 1D np.ndarray, optional
            Detect level [per scale]. The default is [4,3,..,3]

        Returns
        -------
        TabNsigma : 1D np.ndarray
            Detection level per scale.
        """
        TabNsigma = np.zeros(nscale)
        for j in np.arange(nscale):
            vssig = vsize(Nsigma)
            if vssig[0] == 0:
                  TabNsigma[j] = Nsigma
                  if j==0:
                      TabNsigma[j] += 1
            else:
                if vssig[1] > j:
                    TabNsigma[j] = Nsigma[j]
                else:
                    TabNsigma[j] = Nsigma[vssig[1]-1]
        return TabNsigma
                    
    def threshold(self, SigmaNoise=0, Nsigma=3, ThresCoarse=False, hard=True,  FirstDetectScale=0, KillCoarse=False, Verbose=False):
        """
        Apply a hard or a soft thresholding on the coefficients self.coef
        Parameters
        ----------
        SigmaNoise : float, optional
            Noise standard deviation. The default is 0.
            If it is 0, it will be automatically estimated from the first scale.
        Nsigma : 1D np.ndarray, optional
            Detect level [per scale]. The default is [4,3,..,3]  
        ThresCoarse : bool, optional
            If true the coarsest scale is also thresholded. The default is False.
        hard : bool, optional
            IF true, apply hard thresholding, and soft-thresholding otherwise. 
            The default is True.
        FirstDetectScale : int, optional
            Remove the first FirstDetectScale scales. The default is 0.
        KillCoarse :  bool, optional
            IF true the coarsest scale is removed. The default is False.
        Verbose : TYPE, optional
            DESCRIPTION. The default is False.
        Returns
        -------
        None.

        """
        if ThresCoarse:
            Last = self.ns
        else:
            Last = self.ns - 1
        vs=vsize(SigmaNoise)
        dim = vs[0]
        if dim == 0:
            if SigmaNoise == 0:
                SigmaNoise = self.get_noise()
        self.SigmaNoise = SigmaNoise
        if Verbose:
            print("SigmaNoise = ", SigmaNoise, ", vsize(SigmaNoise) = ", vs)
        TabNsigma = self.get_tabsigma(self.ns, Nsigma=Nsigma)
        if Verbose:
            print("TabNsigma = ", TabNsigma)
        for j in np.arange(Last):
            s = self.get_ptr_scale(j)
            if dim == 0:
                Thres = SigmaNoise * TabNsigma[j] * self.TabNorm[j]
            elif dim == 1:
                Thres = SigmaNoise[j] * TabNsigma[j]
            elif dim == 2:
                Thres = SigmaNoise * TabNsigma[j]  * self.TabNorm[j]
            else:
                # print(SigmaNoise.shape)
                Nsig = TabNsigma[j]
                Thres = SigmaNoise[j,:,:] * Nsig
            self.TabNsigma = TabNsigma
            if hard:
                hard_thresholding(s,Thres)
            else:
                soft_thresholding(s,Thres)
            if Verbose:
                print("     scale ",j+1, ", % of non zeros = ", np.count_nonzero(s) * 100. / float(self.nx * self.ny))
        if FirstDetectScale > 0:
            self.coef[0:FirstDetectScale,:,:] = 0.
        if KillCoarse:
            self.coef[self.ns - 1,:,:] = 0.

    def copy(self, name="wt"):
        """
        Duplicate the class, making copy of the coefficients.
        Parameters
        ----------
        name : TYPE, optional
            DESCRIPTION. The default is "wt".

        Returns
        -------
        NewClass : starlet2d
            Copy of the class.
        """
        x = self
        x.name = name
        x.coef=np.zeros((x.ns,x.nx,x.ny))
        x.TabNorm= np.copy(self.TabNorm)
        return x

################################  END CLASS ######################


if __name__ == '__main__':
    print ( "Main :)")
        
    i = readfits("/Users/starck/Main/python/data/ngc2997.fits")
    #[1]
    # PYSAP_CXX=True
# In[1]:

    ns=5  
    testbinding=1
    if testbinding:
        print("TEST BINDING FUNCTION")
        gen2=1
        WT = pysparse.MRStarlet()
        wl = WT.transform(i,ns)
        w = np.stack(wl)
        # WT.info()
        dat_list = []
        for s in range(5):
            dat_list.append(w[s,:,:])
        r = WT.recons(dat_list)
        info(r-i, name="   => resi blinding")
        if (r-i).std() < 1e-5:
            print ("OK  TEST BINDING FUNCTION")
        else:
            print  ("Error in  TEST BINDING FUNCTION")
        print  (" ")
    
# In[2]:
        
    testroutines=1
    if testroutines:
        print("TEST routines starlets")
        bord=2
        gen2=True
        verb=0
        w = star2d(i, ns, gen2=gen2, bord=bord, verb=verb)
        r = istar2d(w, gen2=gen2, bord=bord, verb=verb)
        info(i-r, name="  ==> resi")
        if (r-i).std() < 1e-5:
            print ("OK  TEST 1 routines starlets")
        else:
            print  ("Error in TEST 1 routines starlets")
        
        gen2=False
        w = star2d(i, ns, gen2=gen2)
        r = istar2d(w, gen2=gen2)
        info(i-r, name="  ==> resi")
        if (r-i).std() < 1e-5:
            print ("OK  TEST 2 routines starlets")
        else:
            print  ("Error in TEST 2 routines starlets")
        print  (" ")
            
    
    testclass=1
    if testclass: 
        gen2=False
        l2norm=False
        CW= starlet2d(gen2=gen2, l2norm=l2norm, name="wt C")
        CW.transform(i)
        r = CW.recons()
        info(i-r, name="  ==> resi")
        if (r-i).std() < 1e-5:
            print ("OK  TEST 1 Class starlet(gen1,l1norm)")
        else:
            print  ("Error in TEST 1 Class starlet")
                
        n = np.random.normal(loc=0.0, scale=1., size=(256,256))
        gen2=True
        l2norm=True
        CW= starlet2d(gen2=gen2, l2norm=l2norm, name="wt C2")
        CW.transform(n, WTname='noise')
        CW.stat()
        r = CW.recons()
        info(n-r, name="  ==> resi")
        if (r-n).std() < 1e-5:
            print ("OK  TEST 1 Class starlet (l2norm,gen2)")
        else:
            print  ("Error in TEST 1 Class starlet")
        print  (" ")    
    
    testdenoise=1
    if testdenoise: 
        CW= starlet2d()
        CW.transform(i)
        r = CW.denoising(i)
        info(i-r, name="  ==> resi")
        s=CW.SigmaNoise
        print( s)
        if (r-i).std() < 1.5 * CW.SigmaNoise:
            print ("OK  TEST 1 Class denoise")
        else:
            print  ("Error in TEST 1 Class denoise")
        print  (" ")    
    
    testpos=1
    if testpos:
        CW= starlet2d()
        CW.pos_transform(i, verb=False, pos=True)
        CW.info()
        CW.stat()
        r= CW.recons()
        info(r,name='REC')
        info(i-r, name="  ==> resi")
        ra = (i-r).max()
        if ra.max() < 1.:
            print ("OK  TEST Pos starlet")
        else:
            print  ("Error in TEST Pos starlet")
        print  (" ")    
        #for s in range(CW.ns):
        #    CW.tvs(s)

    testttv=0
    if testttv:
        print ("OK  TEST TV 1")
        CW= starlet2d()
        CW.transform(i)
        for s in range(CW.ns):
            CW.tvs(s)
        #  CW.tvall()
        print ("OK  TEST TV 2 ")
        # CW.tvall(multiview=True)
     


