
import sys
import numpy as N

sys.path.append('')

def plotim():
    """ Plots the image and overlays the island borders with the island number. Also draws the detected gaussians
        at their fwhm radius, with each source being a colour (and line style). """
    bdsm.analysis.plotresults(img)

def getisl(c):
    """ Plots the image and overlays the island borders with the island number. Also draws the detected gaussians
        at their fwhm radius, with each source being a colour (and line style). """
    islid = bdsm.analysis.getisland(img, c)
    return img.islands[islid]

def plot_morph_isl(img, isls=None, thr=None):
    bdsm.analysis.plot_morph_isl(img, isls, thr)

def call_pybdsm(version, parameters):

    if version not in ['stable', 'david', 'test']: raise RuntimeError(version+" Version unknown")
    if version == 'stable': import bdsm_stable as bdsm
    if version == 'david': import bdsm_david as bdsm
    if version == 'test': import bdsm_test as bdsm
    img = bdsm.execute(bdsm.fits_chain, parameters)

    return img, bdsm
 

#img, bdsm = call_pybdsm('test', {'fits_name': "subim.fits", 'beam' : (0.0015, 0.0015, 0.0), 'thresh':"hard", 'atrous_do' : False})

#img, bdsm = call_pybdsm('test', {'fits_name': "concatenated-003-002.restored.fits", 'thresh':"hard", 'atrous_do' : False, 'stop_at' : 'isl'})

#img = bdsm.execute(bdsm.fits_chain,{'fits_name': "HydraA_74MHz_image.fits", 'thresh':"hard", 'atrous_do' : True, 'atrous_bdsm_do' : False, 'atrous_jmax' : 6, 'solnname' : 'del-norms_nobeam_deeper', 'ini_gausfit' : 'nobeam', 'opdir_overwrite' : 'append', 'mean_map' : False, 'rms_map' : False, 'thresh_pix' : 60, 'thresh_isl' : 45})

#img = bdsm.execute(bdsm.fits_chain,{'fits_name': "HydraA_74MHz_image.pybdsm.atrous.w6_norms_deep.fits", 'thresh':"hard", 'atrous_do' : False, 'solnname' : 'nobeam', 'opdir_overwrite' : 'append', 'mean_map' : False, 'rms_map' : False, 'mean_map' : False, 'ini_gausfit' : 'nobeam', 'flag_smallsrc' : False, 'flag_minsnr' : 0.2, 'flag_maxsnr' : 3.0, 'flag_maxsize_isl' : 5.0, 'flag_maxsize_bm' : 45.0, 'flag_minsize_bm' : 0.2})

#img, bdsm = call_pybdsm('test', {'fits_name': "A2255_85CM_BEAM_cut.fits", 'beam' : (0.0167, 0.0167, 0.0), 'thresh':"hard", 'atrous_do' : True, 'atrous_bdsm_do' : True, 'solnname' : 'del', 'ini_gausfit' : 'fbdsm', 'opdir_overwrite' : 'append'})

#img = bdsm.execute(bdsm.fits_chain,{'fits_name': "matteo_mfs.im.fits", 'beam' : (0.002, 0.002, 0.0), 'thresh':"hard"})

#img = bdsm.execute(bdsm.fits_chain,{'fits_name': "WN65341H.fits", 'beam': (.0165515, .01500, 0.0), 'thresh':"hard", 'atrous_do' : False})

#img, bdsm = call_pybdsm('test', {'fits_name': "WN35078H.fits", 'beam': (.0261, .01500, 0.0), 'thresh':"hard", 'atrous_do' : True, 'shapelet_do' : False, 'ini_gausfit' : 'default' })

#img, bdsm = call_pybdsm('test', {'fits_name': "3C274-P.FITS", 'beam': (.00654, .00654, -45.0), 'thresh':"hard", 'atrous_do' : True, 'atrous_jmax' : 5, 'bbs_patches' : 'single', 'solnname' : 'new', 'ini_gausfit' : 'default', 'opdir_overwrite' : 'append', 'atrous_bdsm_do' : True, 'rms_map' : False, 'mean_map' : False, 'thresh_pix' : 100, 'thresh_isl' : 60})

#img = bdsm.execute(bdsm.fits_chain,{'fits_name': "Cas_A-P.models.FITS", 'thresh':"hard", 'atrous_do' : False, 'rms_map' : False, 'mean_map' : False })

#img = bdsm.execute(bdsm.fits_chain,{'fits_name': "VIRA-4.MOD.FITS", 'thresh':"hard", 'atrous_do' : True })

#img = bdsm.execute(bdsm.fits_chain,{'fits_name': "VIRA-4.MOD.pybdsm.atrous.w6.fits", 'thresh':"hard", 'rms_box' : (63, 21)})

#img = bdsm.execute(bdsm.fits_chain,{'fits_name': "Cyg_A-P_mod.FITS", 'thresh':"hard", 'atrous_do' : False , 'rms_map' : False })

#img = bdsm.execute(bdsm.fits_chain,{'fits_name': "Cyg_A-4.model.FITS", 'thresh':"hard", 'atrous_do' : False, 'rms_map' : False , 'thresh_pix' : 6})

#img = bdsm.execute(bdsm.fits_chain,{'fits_name': "Cyg_A-P-cut.fits", 'thresh':"hard", 'atrous_do' : True , 'rms_map' : False, 'mean_map' : 'const', 'thresh_pix' : 1000, 'thresh_isl' : 800, 'ini_gausfit' : 'default', 'solnname' : 'del', 'atrous_bdsm_do' : False})

img, bdsm = call_pybdsm('test' ,{'fits_name': "Cyg_A-P-cut.pybdsm.atrous.w12.fits", 'thresh':"hard", 'atrous_do' : False , 'rms_map' : False, 'mean_map' : 'const', 'ini_gausfit' : 'fbdsm', 'solnname' : 'del', 'opdir_overwrite' : 'append', 'stop_at' : 'isl'})

#img, bdsm = call_pybdsm('test' ,{'fits_name': "Cyg_A-P-cut.pybdsm.atrous.w12.fits", 'thresh':"hard", 'atrous_do' : False , 'rms_map' : False, 'mean_map' : 'const', 'thresh_pix' : 30, 'thresh_isl' : 20, 'ini_gausfit' : 'fbdsm'})

#img = bdsm.execute(bdsm.fits_chain,{'fits_name': "SB128_138-002-002.fits", 'thresh':"hard", 'solnname' : 'try' })

#img = bdsm.execute(bdsm.fits_chain,{'fits_name': "concatenated-000-004.restored.fits", 'rms_box' : (130, 40), 'thresh':"hard", 'atrous_do' : False, 'shapelet_do' : False })

#img = bdsm.execute(bdsm.fits_chain,{'fits_name': "mi_spam.fits", 'beam': (.0222, .0222, 0.0), 'thresh':"hard", 'atrous_do' : False })

#img = bdsm.execute(bdsm.fits_chain,{'fits_name': "marijke.fits", 'beam': (.004, .004, 0.0), 'thresh':"hard", 'atrous_do' : True, 'thresh_isl' : 20, 'thresh_pix' : 100})

#img = bdsm.execute(bdsm.fits_chain,{'fits_name': "SST1cont.image.restored.fits", 'beam': (.008333, .008333, 0.0), 'thresh':"hard", 'atrous_do' : False})

#img = bdsm.execute(bdsm.fits_chain,{'fits_name': "bootbig.FITS", 'beam': (.00154, .00154, 0.0), 'thresh':"hard", \
#                'atrous_do' : True, 'atrous_bdsm_do' : False})

#img = bdsm.execute(bdsm.fits_chain,{'fits_name': "WN35060H", 'beam': (.0165515, .01500, 0.0), 'thresh':"hard"})

#img = bdsm.execute(bdsm.fits_chain,{'fits_name': "lock_cube1.fits", 'beam': (.0054, .0044, 0.0), \
#                'collapse_mode' : 'average', 'collapse_wt' : 'unity', 'beam_sp_derive' : \
#                True, 'atrous_do' : True, 'debug_figs' : True})

#img = bdsm.execute(bdsm.fits_chain,{'fits_name': "newcube1.fits", 'beam': (.00389, .003056, 0.0), \
#                'collapse_mode' : 'average', 'collapse_wt' : 'rms', 'thresh' : 'hard', 'atrous_do' : True})

#img = bdsm.execute(bdsm.fits_chain,{'fits_name': "sim1.1.FITS", 'beam': (.00143, .00143, 0.0),\
#                'collapse_mode' : 'average', 'collapse_wt' : 'rms', 'thresh' : 'hard', 'thresh_pix' : '30'})

#img = bdsm.execute(bdsm.fits_chain,{'fits_name': "A2255_I.fits", 'beam': (.018, .014, 5.0), 'collapse_mode'
#        : 'average', 'collapse_wt' : 'rms', 'thresh' : 'hard', 'thresh_isl' : 20.0, 'thresh_pix' : 50.0,
#        'polarisation_do': True, 'atrous_do' : True})

#img = bdsm.execute(bdsm.fits_chain,{'fits_name': "try.fits", 'beam': (.056, .028, 160.0), 'thresh':"hard", 'thresh_pix':20.})

