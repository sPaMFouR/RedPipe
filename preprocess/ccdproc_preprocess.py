import numpy as np
import ccdproc as ccdp
from astropy.io import fits
from astropy import units as u
from astropy.stats import mad_std
from astropy.nddata import CCDData
from photutils import detect_sources
from astropy.nddata.blocks import block_replicate


READ_NOISE = 4.87  #Units in electrons/photon
GAIN = 1.22 * u.electron / u.adu
SATURATION = 800000 #Units in electrons/photon

def inv_median(a):
    """Used as scaling function in flat combine task."""
    return 1 / np.median(a)

def bias_combine(textlist_bias, master_bias='mbias.fits'):
    """Combines bias files to make a master bias.

    Keyword arguments:
    textlist_bias -- A python list object with paths/names to the
    individual files.
    master_bias -- Output master bias file (default 'mbias.fits')
    """

    new_list = [CCDData.read(filename, unit=u.adu) for filename in textlist_bias]

    combined_bias = ccdp.combine(img_list=new_list, method='average',
                                sigma_clip=True,sigma_clip_low_thresh=5,
                                sigma_clip_high_thresh=5, 
                                sigma_clip_func=np.ma.median,
                                sigma_clip_dev_func=mad_std)

    combined_bias.meta['combined'] = True
    combined_bias.data = combined_bias.data.astype('float32')
    combined_bias.write(master_bias, hdu_mask=None, hdu_uncertainty=None)

def subtract_bias(textlist_files, master_bias='mbias.fits', prefix_str='bs_'):
    """Subtract bias from the given files list
    
    Keyword arguments:
    textlist_files -- A python list object with paths/names to the
    individual files.
    master_bias -- Master bias used for subtraction (default 'mbias.fits')
    prefix_str -- String to be prefixed for subtracted files
    """
    master = CCDData.read(master_bias, unit=u.adu)

    for filename in textlist_files:
        ccd = CCDData.read(filename, unit=u.adu)
        bias_subtracted = ccdp.subtract_bias(ccd=ccd, master=master)
        bias_subtracted.meta['biassub'] = True
        bias_subtracted.data = bias_subtracted.data.astype('float32')
        bias_subtracted.write(prefix_str+filename, hdu_mask=None, hdu_uncertainty=None)

def flat_combine(textlist_files, band, outfile='mflat'):
    """Combines multiple flats for a given input files list

    Keyword arguments:
    textlist_files -- A python list object with paths/names to the
    individual files.
    band -- Filter name associated with the flat files
    outfile -- Master flat name (default outfile + band + '.fits')
    """

    new_list = [CCDData.read(filename, unit=u.adu) for filename in textlist_files]

    combined_flat = ccdp.combine(img_list=new_list, method='average', scale=inv_median,
                                sigma_clip=True, sigma_clip_low_thresh=5, 
                                sigma_clip_high_thresh=5, sigma_clip_func=np.ma.median,
                                sigma_clip_dev_func=mad_std)

    combined_flat.meta['combined'] = True
    file_name = outfile + band + '.fits'
    combined_flat.data = combined_flat.data.astype('float32')
    combined_flat.write(file_name, hdu_mask=None, hdu_uncertainty=None)

def flat_correction(textlist_files, master_flat, prefix_str='f'):
    """To flat field the science images using master flat

    Keyword arguments:
    textlist_files -- A python list object with paths/names to the
    individual files.
    master_flat -- Master flat used to flat field the science image
    prefix_str -- String to be added to newly created sceince image
    """

    master = CCDData.read(master_flat, unit=u.adu)
    for filename in textlist_files:
        ccd = CCDData.read(filename, unit=u.adu)
        flat_corrected = ccdp.flat_correct(ccd=ccd, flat=master, min_value=0.9)
        flat_corrected.meta['flatcorr'] = True
        flat_corrected.data = flat_corrected.data.astype('float32')

        flat_corrected.write(prefix_str + filename, hdu_mask=None, hdu_uncertainty=None)

def cosmic_ray_corr(textlist_files, prefix_str='c'):
    """Gain correction and Cosmic ray correction using LA Cosmic method.

    Keyword arguments:
    textlist_files -- A python list object with paths/names to the
    individual files.
    prefix_str -- String appended to the name of newly created file
    """

    for filename in textlist_files:

        file_corr = CCDData.read(filename, unit=u.adu)
        file_corr = ccdp.gain_correct(file_corr, gain=GAIN)
        new_ccd = ccdp.cosmicray_lacosmic(file_corr, readnoise=READ_NOISE, sigclip=7,
                                             satlevel=SATURATION, niter=4,
                                             gain_apply=False, verbose=True)
        new_ccd.meta['crcorr'] = True
        new_ccd.data = new_ccd.data.astype('float32')
        new_ccd.write(prefix_str + filename, hdu_mask=None, hdu_uncertainty=None)