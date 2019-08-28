# This file provides functions needed to measure signal and profiles
# from Healpix maps at given source catalogs.
# The most important inputs are positions and a healpix map.

import healpy as hp
import numpy as np


def generate_inmaskcat(mask, pos_src, out_edge, tol=1):

    '''
    Return unmasked sources.
    Input: mask: mask healpix map;
           pos_src: position of sources (with shape (Nsource, 3));
           out_edge: out edge of the disk centred at each source(in arcmin);
           tol: tolerance between 0 and 1 indicating a threshold of fraction of
                unmasked pixels in a disk
    Output: indices of unmasked sources.
    '''

    Ns = hp.npix2nside(mask.size)
    Nsource = pos_src.shape[0]
    inmask_ind = np.zeros(Nsource)
    for i in range(Nsource):
        if i % 10000 == 0:
            print i
        list_out = hp.query_disc(Ns, pos_src[i], np.radians(out_edge / 60.))
        npix_out = list_out.size
        neff_out = np.sum(mask[list_out])

        if(neff_out < tol * npix_out):
            continue
        else:
            inmask_ind[i] = 1

    print str(inmask_ind.sum()) + ' sources in mask.'
    return inmask_ind


def get_disk(mask, pos_src, out_edge, in_edge, tol=1):

    '''
    Return unmasked pixel indices of disks centred at a source.
    For usage purpose the input includes two disks.
    Input: mask: mask healpix map;
           pos_src: position of source (with shape=3);
           out_edge: outer disk centred at each source(in arcmin);
           in_edge: inner disk centred at each source(in arcmin);
           tol: tolerance between 0 and 1 indicating a threshold of fraction of
                unmasked pixels in a disk
    Output: pixel indices of both disks.
    '''

    Ns = hp.npix2nside(mask.size)

    list_out = hp.query_disc(Ns, pos_src, np.radians(out_edge / 60.))
    npix_out = len(list_out)

    list_out_unmasked = list_out[mask[list_out] > 0]

    if(list_out_unmasked.size < tol * npix_out):

        return [0, 0]
    lon, lat = hp.vec2dir(pos_src, lonlat=True)
    list_in = hp.query_disc(Ns, pos_src, np.radians(in_edge / 60.))
    list_in_unmasked = list_in[mask[list_in] > 0]
    return list_out_unmasked, list_in_unmasked


def get_bkg(skymap, mask, list_out_unmasked, list_in_unmasked):

    '''
    Return mean background signal for a source which is calculated by 
    taking mean of signals in a ring around a source.

    Input: skymap: healpix skymap;
           mask: mask healpix map, must have the same Nside as skymap;
           list_out_unmasked: pixel indices of outer disks of the ring;
           list_in_unmasked: pixel indices of inner disks of the ring.
    Output: mean background signal
    '''

    if np.asarray(list_out_unmasked).sum() == 0:
        return 0
    neff_in = sum(mask[list_in_unmasked])
    neff_out = sum(mask[list_out_unmasked])
    bkg = (sum(skymap[list_out_unmasked] * mask[list_out_unmasked]) -
           sum(skymap[list_in_unmasked] * mask[list_in_unmasked])) / \
           (neff_out-neff_in)  # Background
    return bkg


def get_totalsignal(skymap, mask, pos_src_all, out_edge, in_edge, signal_edge=15,
                  tol=1):

    '''
    Return integrated signal signal in a disk around a collection of sources

    Input: skymap: healpix skymap;
           mask: mask healpix map;
           pos_src_all: position of sources (with shape (Nsource, 3));
           out_edge: outer disk of the ring to calculate background(in arcmin);
           in_edge: inner disk of the ring to calculate background(in arcmin);
           signal_edge: radius of the disk to calculate total signal;
           tol: tolerance between 0 and 1 indicating a threshold of fraction of
                unmasked pixels in a disk
    Output: total signal corresponding to each source
    '''

    Nsource = pos_src_all.shape[0]
    skymap_masked = skymap * mask
    Ns = hp.npix2nside(mask.size)
    tot_signal = np.zeros(Nsource)
    fraction = 0
    print('Total source number: '+str(Nsource))
    print('Calculating...')
    for i in range(Nsource):
        fraction_next = np.int(i * 10 / np.float(Nsource))
        if fraction_next > fraction:
            print(int(i*100 / np.float(Nsource))), '%'
            fraction = fraction_next
        list_out_unmasked, list_in_unmasked = get_disk(mask, pos_src_all[i], out_edge,
                                                       in_edge, tol)
        bkg = get_bkg(skymap, mask, list_out_unmasked, list_in_unmasked)
        list_signal = hp.query_disc(Ns, pos_src_all[i], np.radians(signal_edge / 60.))
        tot_signal[i] = ((skymap_masked)[list_signal]-bkg).sum()
    return tot_signal


def get_centralsignal(skymap, mask, pos_src_all, out_edge, in_edge, fwhm,
                    tol=1):

    '''
    Return central signal of a collection of sources which is calculated by
    taking the signal of the pixel which countains the source then multiply by 
    beam normalization factor.

    Input: skymap: healpix skymap;
           mask: mask healpix map;
           pos_src_all: position of sources (with shape (Nsource, 3));
           out_edge: outer disk of the ring to calculate background(in arcmin);
           in_edge: inner disk of the ring to calculate background(in arcmin);
           fwhm: FWHM of the gaussian beam of the skymap
           tol: tolerance between 0 and 1 indicating a threshold of fraction of
                unmasked pixels in a disk
    Output: central signal corresponding to each source
    '''

    sigma = fwhm / 2.35482004503
    Nsource = pos_src_all.shape[0]
    skymap_masked = skymap * mask
    Ns = hp.npix2nside(mask.size)
    tot_signal = np.zeros(Nsource)
    fraction = 0
    print('Total source number: '+str(Nsource))
    print('Calculating...')
    for i in range(Nsource):
        fraction_next = np.int(i * 10 / np.float(Nsource))
        if fraction_next > fraction:
            print(int(i*100 / np.float(Nsource))), '%'
            fraction = fraction_next
        list_out_unmasked, list_in_unmasked = get_disk(mask, pos_src_all[i], out_edge,
                                                       in_edge, tol)
        bkg = get_bkg(skymap, mask, list_out_unmasked, list_in_unmasked)
        pos_src_ind = hp.vec2pix(Ns, pos_src_all[i][0], pos_src_all[i][1],
                                 pos_src_all[i][2])
        tot_signal[i] = ((skymap_masked)[pos_src_ind]-bkg)
    tot_signal *= (sigma * np.sqrt(2*np.pi))
    return tot_signal


def get_npix(mask, rbins, pos_src, out_edge, tol=1):

    '''
    Return the number of sources within each r bins around a source.

    Input: mask: mask healpix map;
           rbins: r bins specified by the user. Need to be an int or an array
           pos_src: position of the source (with size=3);
           out_edge: radius of the disk that is taken into account. in arcmin;
           tol: tolerance between 0 and 1 indicating a threshold of fraction of
                unmasked pixels in a disk
    Output: number of sources within each r bins around a source.
    '''

    list_out_unmasked, list_in_unmasked = get_disk(mask, pos_src, out_edge, 0,
                                                   tol)
    if np.asarray(list_out_unmasked).sum() == 0:
        return np.zeros(rbins.size-1)
    Ns = hp.npix2nside(mask.size)
    vector_disc = np.array(hp.pix2vec(Ns, list_out_unmasked))
    # Get the unit vector of each pixels within the out disc
    innprod = np.round(np.dot(pos_src, vector_disc), 15)
    innprod[innprod > 1] = 1
    rtheta = np.degrees(np.arccos(innprod)) * 60.

    normal = np.histogram(rtheta, bins=rbins)[0]
    return normal


def get_1d_profile(skymap, mask, pos_src, rbins, out_edge, in_edge, tol=1):
    list_out_unmasked, list_in_unmasked = get_disk(mask, pos_src, out_edge,
                                                   in_edge, tol)

    '''
    Calculate the 1d signal profile of a single source

    Input: skymap: healpix skymap;
           mask: mask healpix map;
           pos_src: position of the source (with size=3);
           rbins: r bins specified by the user. Need to be an int or an array
           out_edge: outer disk of the ring to calculate background(in arcmin);
           in_edge: inner disk of the ring to calculate background(in arcmin);
           tol: tolerance between 0 and 1 indicating a threshold of fraction of
                unmasked pixels in a disk
    Output: center of each r bin, mean signal in each r bin; standard
            deviation in each r bin; number of pixels in each r bin
    '''

    r_center = 0.5 * (rbins[1:] + rbins[:-1])
    if np.asarray(list_out_unmasked).sum() == 0:
        return r_center, np.zeros(rbins.size-1),\
            np.zeros(rbins.size-1), np.zeros(rbins.size-1)

    bkg = get_bkg(skymap, mask, list_out_unmasked, list_in_unmasked)

    Ns = hp.npix2nside(mask.size)
    vector_disc = np.array(hp.pix2vec(Ns, list_out_unmasked))
    innprod = np.round(np.dot(pos_src, vector_disc), 15)
    innprod[innprod > 1] = 1
    rtheta = np.degrees(np.arccos(innprod)) * 60.

    pureT_1d = skymap[list_out_unmasked] - bkg
    signal, xedges = np.histogram(rtheta, bins=rbins, weights=pureT_1d)
    npix, xedges = np.histogram(rtheta, bins=rbins)
    signal2, xedges = np.histogram(rtheta, bins=rbins, weights=pureT_1d**2)
    std = (signal2 * npix - signal**2) ** 0.5 / npix
    signal = signal / npix
    return r_center, signal, std, npix


def get_2d_profile(skymap, mask, pos_src, xbins, ybins, out_edge, in_edge,
                   tol=1):

    '''
    Calculate the 2d signal profile(map) of a single source

    Input: skymap: healpix skymap;
           mask: mask healpix map;
           pos_src: position of the source (with size=3);
           xbins: array, bins in x direction specified by the user.
           ybins: array, bins in y direction specified by the user.
           out_edge: outer disk of the ring to calculate background(in arcmin);
           in_edge: inner disk of the ring to calculate background(in arcmin);
           tol: tolerance between 0 and 1 indicating a threshold of fraction of
                unmasked pixels in a disk
    Output: 2d signal profile(map) of a single source
    '''

    list_out_unmasked, list_in_unmasked = get_disk(mask, pos_src, out_edge,
                                                   in_edge, tol)
    bkg = get_bkg(skymap, mask, list_out_unmasked, list_in_unmasked)

    xybin = np.transpose([np.tile(xbins, len(ybins)),
                          np.repeat(ybins, len(xbins))])
    Ns = hp.npix2nside(mask.size)
    nbins = xbins.size
    x, y = xybin.T
    lon, lat = hp.vec2dir(pos_src, lonlat=True)

    lat_pix = y / 60. + lat  # in degrees
    lon_pix = x / 60. / np.cos(np.radians(lat)) + lon
    try:
        x_pix, y_pix, z_pix = hp.ang2vec(lon_pix, lat_pix, lonlat=True).T
    except ValueError, e:
        print(e)
        return np.zeros((nbins, nbins))
    list_eff = hp.vec2pix(Ns, x_pix, y_pix, z_pix)
    pureT_2d = skymap[list_eff] - bkg  # Effective sz signal
    hit_2d = mask[list_eff]
    if hit_2d.sum() < tol * hit_2d.size:
        return np.zeros((nbins, nbins))
    return pureT_2d.reshape(nbins, nbins) / hit_2d.reshape(nbins, nbins)


def stack_1d_profile(skymap, mask, pos_src_all, rbins, out_edge, in_edge,
                     tol=1):

    '''
    Stack the 1d signal profile of a collection of sources.

    Input: skymap: healpix skymap;
           mask: mask healpix map;
           pos_src_all: position of the source (with shape=(Nsource, 3));
           rbins: r bins specified by the user. Need to be an int or an array
           out_edge: outer disk of the ring to calculate background(in arcmin);
           in_edge: inner disk of the ring to calculate background(in arcmin);
           tol: tolerance between 0 and 1 indicating a threshold of fraction of
                unmasked pixels in a disk
    Output: center of each r bin, stacked signal in each r bin; standard
            deviation in each r bin; covariance matrix; correlation
            coefficients
    '''

    Nsource = pos_src_all.shape[0]
    N_inmask = pos_src_all.shape[0]
    nbins = rbins.size - 1
    profile_all = np.zeros((Nsource, nbins))
    std_all = np.zeros((Nsource, nbins))
    std_weight = np.zeros((Nsource, nbins))
    npix_all = np.zeros((Nsource, nbins))
    hit_all = np.zeros((Nsource, nbins))
    fraction = 0
    print('Total source number: '+str(Nsource))
    print('Calculating...')

    for i in range(Nsource):
        fraction_next = np.int(i * 10 / np.float(Nsource))
        if fraction_next > fraction:
            print(int(i*100 / np.float(Nsource))), '%'
            fraction = fraction_next
        r_center, profile_all[i], std_all[i], \
            npix_all[i] = get_1d_profile(skymap, mask,
                                         pos_src_all[i],
                                         rbins, out_edge,
                                         in_edge, tol=1)
        if npix_all[i].sum() == 0:
            N_inmask = N_inmask - 1
        hit_all[i] = (npix_all[i] > 0)
        std_weight[i] = 1 / std_all[i]
        std_weight[i][np.where(npix_all[i] == 0)[0]] = 0

    npix_all_stacked = np.sum(npix_all, axis=0)
    if np.where(npix_all_stacked == 0)[0].size > 0:
        raise ValueError('Some rbins get no pixel at all!')

    hit_all_stacked = np.sum(hit_all, axis=0)

    profile_masked = np.ma.masked_values(profile_all, 0.0)
    profile_stacked = np.asarray(np.ma.mean(profile_masked, axis=0))
    profile_cov = np.asarray(np.ma.cov(profile_masked, rowvar=False)) / \
        np.outer(hit_all_stacked, hit_all_stacked) ** 0.5
    profile_corr = np.asarray(np.ma.corrcoef(profile_masked, rowvar=False))
    profile_std = np.diag(profile_cov) ** 0.5
    print("Number of Objects (in mask): ", N_inmask)
    print("Used " + str(N_inmask/float(Nsource) * 100) + '%')

    return r_center, profile_stacked, profile_std, profile_cov, profile_corr


def stack_2d_profile(skymap, mask, pos_src_all, xbins, ybins, out_edge,
                     in_edge, tol=1):

    '''
    Stack the 2d signal profiles(maps) of a collection of sources

    Input: skymap: healpix skymap;
           mask: mask healpix map;
           pos_src_all: position of the source (with shape=(Nsource, 3));
           xbins: array, bins in x direction specified by the user.
           ybins: array, bins in y direction specified by the user.
           out_edge: outer disk of the ring to calculate background(in arcmin);
           in_edge: inner disk of the ring to calculate background(in arcmin);
           tol: tolerance between 0 and 1 indicating a threshold of fraction of
                unmasked pixels in a disk
    Output: stacked 2d signal profile(map)
    '''

    Nsource = pos_src_all.shape[0]
    N_inmask = pos_src_all.shape[0]
    nbins = xbins.size
    profile_all = np.zeros((Nsource, nbins, nbins))
    fraction = 0
    print('Total source number: '+str(Nsource))
    print('Calculating...')
    for i in range(Nsource):
        fraction_next = np.int(i * 10 / np.float(Nsource))
        if fraction_next > fraction:
            print(int(i*100 / np.float(Nsource))), '%'
            fraction = fraction_next

        profile_all[i] = get_2d_profile(skymap, mask, pos_src_all[i], xbins,
                                        ybins, out_edge, in_edge)
        if profile_all[i].sum() == 0:
            N_inmask = N_inmask - 1
    profile_masked = np.ma.masked_values(profile_all, 0.0)
    profile_stacked = np.asarray(np.ma.mean(profile_masked, axis=0))

    print("Number of Objects (in mask) ", N_inmask)
    print("Used " + str(N_inmask/float(Nsource) * 100) + '%')

    return profile_stacked
