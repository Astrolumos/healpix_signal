# This file defines classes needed to measure signal and profiles
# from Healpix maps at given source catalogs.
# The most important inputs are positions and a healpix map.

import healpy as hp
import numpy as np


class signal_profile:

    def __init__(self, skymap, mask, pos_src_all, fwhm):
        self.skymap = skymap
        self.mask = mask
        self.pos_src_all = pos_src_all
        self.Nside = hp.npix2nside(mask.size)
        self.Nsource = pos_src_all.shape[0]
        self.fwhm = fwhm
        self.sigma = fwhm / 2.35482004503

    def generate_inmaskcat(self, out_edge, tol=1):

        '''
        Return unmasked sources.
        Input: out_edge: out edge of the disk centred at each source(in arcmin);
               tol: tolerance between 0 and 1 indicating a threshold of fraction of
                unmasked pixels in a disk
        Output: indices of unmasked sources.
        '''
        mask = self.mask
        Ns = self.Nside
        Nsource = self.Nsource
        pos_src_all = self.pos_src_all
        inmask_ind = np.zeros(Nsource)
        for i in range(Nsource):
            if i % 10000 == 0:
                print i
            list_out = hp.query_disc(Ns, pos_src_all[i], np.radians(out_edge / 60.))
            npix_out = list_out.size
            neff_out = np.sum(mask[list_out])

            if(neff_out < tol * npix_out):
                continue
            else:
                inmask_ind[i] = 1

        print str(inmask_ind.sum()) + ' sources in mask.'
        return inmask_ind

    def get_disk(self, src_ind, out_edge, in_edge, tol=1):

        '''
        Return unmasked pixel indices of disks centred at a source.
        For usage purpose the input includes two disks.
        Input: src_ind: index of source;
               out_edge: outer disk centred at each source(in arcmin);
               in_edge: inner disk centred at each source(in arcmin);
               tol: tolerance between 0 and 1 indicating a threshold of 
               fraction of unmasked pixels in a disk
        Output: pixel indices of both disks.
        '''
        mask = self.mask
        Ns = self.Nside
        pos_src = self.pos_src_all[src_ind]

        list_out = hp.query_disc(Ns, pos_src, np.radians(out_edge / 60.))
        npix_out = len(list_out)

        list_out_unmasked = list_out[mask[list_out] > 0]

        if(list_out_unmasked.size < tol * npix_out):
            return [0, 0]
        lon, lat = hp.vec2dir(pos_src, lonlat=True)
        list_in = hp.query_disc(Ns, pos_src, np.radians(in_edge / 60.))
        list_in_unmasked = list_in[mask[list_in] > 0]
        return list_out_unmasked, list_in_unmasked

    def get_bkg(self, src_ind, out_edge, in_edge, tol=1):

        '''
        Return mean background signal for a source which is calculated by 
        taking mean of signals in a ring around a source.

        Input: src_ind: index of source;
               list_out_unmasked: pixel indices of outer disks of the ring;
               list_in_unmasked: pixel indices of inner disks of the ring.
        Output: mean background signal
        '''

        skymap = self.skymap
        mask = self.mask
        list_out_unmasked, list_in_unmasked = self.get_disk(src_ind,
                                                            out_edge,
                                                            in_edge, tol)
        if np.asarray(list_out_unmasked).sum() == 0:
            return 0
        neff_in = sum(mask[list_in_unmasked])
        neff_out = sum(mask[list_out_unmasked])
        bkg = (sum(skymap[list_out_unmasked] * mask[list_out_unmasked]) -
               sum(skymap[list_in_unmasked] * mask[list_in_unmasked])) / \
            (neff_out-neff_in)  # Background
        return bkg

    def get_totalsignal(self, out_edge, in_edge, signal_edge=15,
                        tol=1):

        '''
        Return integrated signal signal in a disk around a collection of
        sources

        Input: out_edge: outer disk of the ring to calculate background
                         (in arcmin);
               in_edge: inner disk of the ring to calculate background
                        (in arcmin);
               signal_edge: radius of the disk to calculate total signal;
               tol: tolerance between 0 and 1 indicating a threshold of
                    fraction of unmasked pixels in a disk
        Output: total signal corresponding to each source
        '''

        Nsource = self.Nsource
        skymap = self.skymap
        mask = self.mask
        skymap_masked = skymap * mask
        Ns = self.Nside
        pos_src_all = self.pos_src_all
        tot_signal = np.zeros(Nsource)
        fraction = 0
        print('Total source number: '+str(Nsource))
        print('Calculating...')
        for i in range(Nsource):
            fraction_next = np.int(i * 10 / np.float(Nsource))
            if fraction_next > fraction:
                print(int(i*100 / np.float(Nsource))), '%'
                fraction = fraction_next
            list_out_unmasked, list_in_unmasked = self.get_disk(mask,
                                                                pos_src_all[i],
                                                                out_edge,
                                                                in_edge, tol)
            bkg = self.get_bkg(skymap, mask, list_out_unmasked,
                               list_in_unmasked)
            list_signal = hp.query_disc(Ns, pos_src_all[i],
                                        np.radians(signal_edge / 60.))
            tot_signal[i] = ((skymap_masked)[list_signal]-bkg).sum()
        return tot_signal

    def get_centralsignal(self, out_edge, in_edge, tol=1):

        '''
        Return central signal of a collection of sources which is calculated by
        taking the signal of the pixel which countains the source then
        multiply by beam normalization factor.

        Input: out_edge: outer disk of the ring to calculate background
                         (in arcmin);
               in_edge: inner disk of the ring to calculate background
                        (in arcmin);
               tol: tolerance between 0 and 1 indicating a threshold of
                    fraction of unmasked pixels in a disk
        Output: central signal corresponding to each source
        '''

        Nsource = self.Nsource
        skymap = self.skymap
        mask = self.mask
        skymap_masked = skymap * mask
        Ns = self.Nside
        pos_src_all = self.pos_src_all
        sigma = self.sigma
        tot_signal = np.zeros(Nsource)
        fraction = 0
        print('Total source number: '+str(Nsource))
        print('Calculating...')
        for i in range(Nsource):
            fraction_next = np.int(i * 10 / np.float(Nsource))
            if fraction_next > fraction:
                print(int(i*100 / np.float(Nsource))), '%'
                fraction = fraction_next
            list_out_unmasked, list_in_unmasked = self.get_disk(mask,
                                                                pos_src_all[i],
                                                                out_edge,
                                                                in_edge, tol)
            bkg = self.get_bkg(skymap, mask, list_out_unmasked,
                               list_in_unmasked)
            pos_src_ind = hp.vec2pix(Ns, pos_src_all[i][0], pos_src_all[i][1],
                                     pos_src_all[i][2])
            tot_signal[i] = ((skymap_masked)[pos_src_ind]-bkg)
        tot_signal *= (sigma * np.sqrt(2*np.pi))
        return tot_signal

    def get_npix(self, rbins, src_ind, out_edge, tol=1):

        '''
        Return the number of sources within each r bins around a source.

        Input: rbins: r bins specified by the user. Need to be an int or an array
               src_ind: index of source;
               out_edge: radius of the disk that is taken into account. in arcmin;
               tol: tolerance between 0 and 1 indicating a threshold of fraction of
                    unmasked pixels in a disk
        Output: number of sources within each r bins around a source.
        '''
        mask = self.mask
        pos_src = self.pos_src_all[src_ind]
        list_out_unmasked, list_in_unmasked = self.get_disk(mask, pos_src,
                                                            out_edge, 0,
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

    def get_1d_profile(self, src_ind, rbins, out_edge, in_edge, tol=1):

        '''
        Calculate the 1d signal profile of a single source

        Input: src_ind: index of source;
               rbins: r bins specified by the user. Need to be an int or an array
               out_edge: outer disk of the ring to calculate background(in arcmin);
               in_edge: inner disk of the ring to calculate background(in arcmin);
               tol: tolerance between 0 and 1 indicating a threshold of fraction of
                    unmasked pixels in a disk
        Output: center of each r bin, mean signal in each r bin; standard
                deviation in each r bin; number of pixels in each r bin
        '''

        mask = self.mask
        skymap = self.skymap
        pos_src = self.pos_src_all[src_ind]
        list_out_unmasked, list_in_unmasked = self.get_disk(src_ind,
                                                            out_edge,
                                                            in_edge, tol)

        r_center = 0.5 * (rbins[1:] + rbins[:-1])
        if np.asarray(list_out_unmasked).sum() == 0:
            return r_center, np.zeros(rbins.size-1),\
                np.zeros(rbins.size-1), np.zeros(rbins.size-1)

        bkg = self.get_bkg(src_ind, out_edge, in_edge, tol)

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
        signal[npix==0] = 0.0
        return r_center, signal, std, npix

    def stack_1d_profile(self, rbins, out_edge, in_edge,
                         tol=1):

        '''
        Stack the 1d signal profile of a collection of sources.

        Input: rbins: r bins specified by the user. Need to be an int or array
               out_edge: outer disk of the ring to calculate background(in arcmin);
               in_edge: inner disk of the ring to calculate background(in arcmin);
               tol: tolerance between 0 and 1 indicating a threshold of fraction of
                    unmasked pixels in a disk
        Output: center of each r bin, stacked signal in each r bin; standard
                deviation in each r bin; covariance matrix; correlation
                coefficients
        '''

        Nsource = self.Nsource
        N_inmask = self.generate_inmaskcat(out_edge, tol).sum()
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
                npix_all[i] = self.get_1d_profile(i, rbins, out_edge,
                                                  in_edge, tol)
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

    def get_2d_profile(self, src_ind, xbins, ybins, out_edge, in_edge,
                       tol=1):

        '''
        Calculate the 2d signal profile(map) of a single source

        Input: src_ind: index of source;
               xbins: array, bins in x direction specified by the user.
               ybins: array, bins in y direction specified by the user.
               out_edge: outer disk of the ring to calculate background(in arcmin);
               in_edge: inner disk of the ring to calculate background(in arcmin);
               tol: tolerance between 0 and 1 indicating a threshold of fraction of
                    unmasked pixels in a disk
        Output: 2d signal profile(map) of a single source
        '''

        mask = self.mask
        skymap = self.skymap
        pos_src = self.pos_src_all[src_ind]

        list_out_unmasked, list_in_unmasked = self.get_disk(src_ind,
                                                            out_edge,
                                                            in_edge, tol)
        bkg = self.get_bkg(src_ind, out_edge, in_edge, tol)

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

    def stack_2d_profile(self, xbins, ybins, out_edge,
                         in_edge, tol=1):

        '''
        Stack the 2d signal profiles(maps) of a collection of sources

        Input: xbins: array, bins in x direction specified by the user.
               ybins: array, bins in y direction specified by the user.
               out_edge: outer disk of the ring to calculate background(in arcmin);
               in_edge: inner disk of the ring to calculate background(in arcmin);
               tol: tolerance between 0 and 1 indicating a threshold of fraction of
                    unmasked pixels in a disk
        Output: stacked 2d signal profile(map)
        '''

        Nsource = self.Nsource
        N_inmask = self.generate_inmaskcat(out_edge, tol).sum()
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

            profile_all[i] = self.get_2d_profile(i, xbins,
                                                 ybins, out_edge, in_edge)
        profile_masked = np.ma.masked_values(profile_all, 0.0)
        profile_stacked = np.asarray(np.ma.mean(profile_masked, axis=0))

        print("Number of Objects (in mask) ", N_inmask)
        print("Used " + str(N_inmask/float(Nsource) * 100) + '%')

        return profile_stacked
