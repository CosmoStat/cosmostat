##########################################################################
# XXX - Copyright (C) XXX, 2017
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# System import
import numpy
import scipy.signal


def nsNorm(nx,ny,nz,trans=1):
    im = numpy.zeros((nx,ny))
    im[nx/2,ny/2] = 1
    if trans == 1:                       # starlet transform 2nd generation
        wt = star2d(im,nz,fast=True,gen2=True,normalization=False)
        tmp = wt**2
    elif trans == 2:                      # starlet transform 1st generation
        wt = star2d(im,nz,fast=True,gen2=False,normalization=False)
        tmp = wt**2
    tabNs = numpy.sqrt(numpy.sum(numpy.sum(tmp,1),1))
    return tabNs

def wavOrth2d(im, nz, wname='haar', wtype=1, mode="symmetric"):
    sx,sy = numpy.shape(im)
    scale = nz
    if scale > numpy.ceil(numpy.log2(sx))+1 or scale > numpy.ceil(numpy.log2(sy))+1:
        print("Too many decomposition scales! The decomposition scale will be set to default value: 1!")
        scale = 1
    if scale < 1:
        print("Decomposition scales should not be smaller than 1! The decomposition scale will be set to default value: 1!")
        scale = 1

    band = numpy.zeros((scale+1,len(numpy.shape(im))))
    band[-1] = numpy.shape(im)

    if wname =='haar' or wname == 'db1' or wname == 'db2' or wname == 'db3' or wname == 'db4' or wname == 'db5':
        wtype = 1
    else:
        wtype = 2

    (h0,g0) = wavFilters(wname,wtype,'d')
    lf = numpy.size(h0)
    wt = numpy.array([])
    start = numpy.array([1,1])

    for sc in numpy.arange(scale-1):
        end = numpy.array([sx + lf - 1, sy + lf - 1])
        lenExt = lf - 1
        imExtCol = numpy.lib.pad(im, ((0,0),(lenExt,lenExt)), mode)      # Extension of columns
        tmp = scipy.signal.convolve2d(imExtCol,h0[numpy.newaxis,:],'valid')
        im = convdown(tmp,h0[numpy.newaxis,:],lenExt,start,end)                             # Approximation
        hor = convdown(tmp,g0[numpy.newaxis,:],lenExt,start,end)                             # Horizontal details
        tmp = scipy.signal.convolve2d(imExtCol,g0[numpy.newaxis,:],'valid')
        vet = convdown(tmp,h0[numpy.newaxis,:],lenExt,start,end)                             # Vertical details
        dig = convdown(tmp,g0[numpy.newaxis,:],lenExt,start,end)                             # Diagonal details
        wt = numpy.hstack([hor.flatten(),vet.flatten(),dig.flatten(),wt])
        sx,sy = numpy.shape(im)
        band[-2-sc] = numpy.array([sx,sy])
    wt = numpy.hstack([im.flatten(),wt])
    band[0] = numpy.shape(im)
    return wt,band


def convdown(x,F,lenExt,start,end):
    im = numpy.copy(x[:,start[1]:end[1]:2])                  # Downsampling
    y = numpy.lib.pad(im, ((lenExt,lenExt),(0,0)), 'symmetric')      # Extension of rows
    y = scipy.signal.convolve2d(y.T,F,'valid')
    y = y.T
    y = y[start[0]:end[0]:2,:]
    return y


def wavFilters(wname,wtype,mode):
    if wtype == 1:
        F = scaleFilter(wname,1)
        (Lo_D,Hi_D,Lo_R,Hi_R) = orthWavFilter(F)
    elif wtype == 2:
        (Rf,Df) = scaleFilter(wname,2)
        [Lo_D,Hi_D1,Lo_R1,Hi_R,Lo_D2,Hi_D,Lo_R,Hi_R2] = biorWavFilter(Rf,Df)
    if mode =='d':
        return (Lo_D,Hi_D)
    elif mode =='r':
        return (Lo_R,Hi_R)
    elif mode == 'l':
        return (Lo_D,Lo_R)
    elif mode == 'h':
        return (Hi_D,Hi_R)
    elif mode == "all":
        return Hi_D, Lo_D, Hi_R, Lo_R, len(Hi_D), len(Hi_R)


def haar():
    return wavFilters(wname="haar", wtype=1, mode="all")


def bio97():
    return wavFilters(wname="9/7", wtype=2, mode="all")


def scaleFilter(wname,wtype):
    if wtype == 1:
        if wname =='haar' or wname == 'db1':
            F = numpy.array([0.5,0.5])

        elif wname == 'db2':
            F = numpy.array([0.34150635094622,0.59150635094587,0.15849364905378,-0.09150635094587])

        elif wname == 'db3':
            F = numpy.array([0.23523360389270,0.57055845791731,0.32518250026371,-0.09546720778426,\
                          -0.06041610415535,0.02490874986589])

        elif wname == 'db4':
            F = numpy.array([0.16290171402562,0.50547285754565,0.44610006912319,-0.01978751311791,\
                          -0.13225358368437,0.02180815023739,0.02325180053556,-0.00749349466513])

        elif wname == 'db5':
            F = numpy.array([0.11320949129173,0.42697177135271,0.51216347213016,0.09788348067375,\
                           -0.17132835769133,-0.02280056594205,0.05485132932108,-0.00441340005433,\
                           -0.00889593505093,0.00235871396920])
        return F
    elif wtype == 2:
        if wname == '9/7' or wname == 'bior4.4':
            Df = numpy.array([0.0267487574110000,-0.0168641184430000,-0.0782232665290000,0.266864118443000,\
                           0.602949018236000,0.266864118443000,-0.0782232665290000,-0.0168641184430000,\
                           0.0267487574110000])
            Rf = numpy.array([-0.0456358815570000,-0.0287717631140000,0.295635881557000,0.557543526229000,\
                           0.295635881557000,-0.0287717631140000,-0.0456358815570000])
        elif wname == 'bior1.1':
            Df = numpy.array([0.5])
            Rf = numpy.array([0.5])
            Df = numpy.hstack((Df,Df[::-1]))
            Rf = numpy.hstack((Rf,Rf[::-1]))
        elif wname == 'bior1.3':
            Df = numpy.array([-1./16,1./16,1./2])
            Rf = numpy.array([0.5])
            Df = numpy.hstack((Df,Df[::-1]))
            Rf = numpy.hstack((Rf,Rf[::-1]))
        elif wname == 'bior1.5':
            Df = numpy.array([3./256,-3./256,-11./128,11./128,1./2])
            Rf = numpy.array([0.5])
            Df = numpy.hstack((Df,Df[::-1]))
            Rf = numpy.hstack((Rf,Rf[::-1]))
        elif wname == 'bior2.2':
            Df = numpy.array([-1./8,1./4,3./4,1./4,-1./8])
            Rf = numpy.array([1./4,1./2,1./4])
        elif wname == 'bior2.4':
            Df = numpy.array([3./128,-3./64,-1./8,19./64,45./64,19./64,-1./8,-3./64,3./128])
            Rf = numpy.array([1./4,1./2,1./4])
        elif wname == 'bior2.6':
            Df = numpy.array([-5./1024,5./512,17./512,-39./512,-123./1024,81./256,175./256,81./256,-123./1024,-39./512,17./512,5./512,-5./1024])
            Rf = numpy.array([1./4,1./2,1./4])
        elif wname == 'bior2.8':
            Df = 1./(2**15)*numpy.array([35,-70,-300,670,1228,-3126,-3796,10718,22050,10718,-3796,-3126,1228,670,-300,-70,35])
            Rf = numpy.array([1./4,1./2,1./4])
        elif wname == 'bior3.1':
            Df = 1./4*numpy.array([-1,3])
            Rf = 1./8*numpy.array([1,3])
            Df = numpy.hstack((Df,Df[::-1]))
            Rf = numpy.hstack((Rf,Rf[::-1]))
        elif wname == 'bior3.3':
            Df = 1./64*numpy.array([3,-9,-7,45])
            Rf = 1./8*numpy.array([1,3])
            Df = numpy.hstack((Df,Df[::-1]))
            Rf = numpy.hstack((Rf,Rf[::-1]))
        elif wname == 'bior3.5':
            Df = 1./512*numpy.array([-5,15,19,-97,-26,350])
            Rf = 1./8*numpy.array([1,3])
            Df = numpy.hstack((Df,Df[::-1]))
            Rf = numpy.hstack((Rf,Rf[::-1]))
        elif wname == 'bior3.7':
            Df = 1./(2**14)*numpy.array([35,-105,-195,865,363,-3489,-307,11025])
            Rf = 1./8*numpy.array([1,3])
            Df = numpy.hstack((Df,Df[::-1]))
            Rf = numpy.hstack((Rf,Rf[::-1]))
        elif wname == 'bior3.9':
            Df = 1./(2**17)*numpy.array([-63,189,469,-1911,-1308,9188,1140,-29676,190,87318])
            Rf = 1./8*numpy.array([1,3])
            Df = numpy.hstack((Df,Df[::-1]))
            Rf = numpy.hstack((Rf,Rf[::-1]))
        return (Rf,Df)

def orthWavFilter(F):
    p = 1
#     h1 = numpy.copy(F)
    Lo_R = numpy.sqrt(2)*F/numpy.sum(F)
#     Lo_R = F/numpy.sqrt(numpy.sum(F**2))
    Hi_R = numpy.copy(Lo_R[::-1])
    first = 2-p%2
#     print first
#     print tmp
    Hi_R[first::2] = -Hi_R[first::2]
    Hi_D=numpy.copy(Hi_R[::-1])
    Lo_D=numpy.copy(Lo_R[::-1])
    return (Lo_D,Hi_D,Lo_R,Hi_R)


def biorWavFilter(Rf,Df):
    lr = len(Rf)
    ld = len(Df)
    lmax = max(lr,ld)
    if lmax%2:
        lmax += 1
    Rf = numpy.hstack([numpy.zeros((lmax-lr)/2),Rf,numpy.zeros((lmax-lr+1)/2)])
    Df = numpy.hstack([numpy.zeros((lmax-ld)/2),Df,numpy.zeros((lmax-ld+1)/2)])

    [Lo_D1,Hi_D1,Lo_R1,Hi_R1] = orthWavFilter(Df)
    [Lo_D2,Hi_D2,Lo_R2,Hi_R2] = orthWavFilter(Rf)

    return (Lo_D1,Hi_D1,Lo_R1,Hi_R1,Lo_D2,Hi_D2,Lo_R2,Hi_R2)
