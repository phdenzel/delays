import sys
import os
import pprint
import numpy as np
from scipy import ndimage
import matplotlib as mpl
import matplotlib.pyplot as plt
root = "/Users/phdenzel/gleam"
sys.path.append(root)
import gleam
from gleam.lensobject import LensObject
from gleam.utils.lensing import LensModel
from gleam.reconsrc import ReconSrc, run_model
from gleam.utils.plotting import kappa_map_plot, kappa_profiles_plot
from gleam.utils.plotting import arrival_time_surface_plot
from gleam.utils.plotting import plot_scalebar, plot_labelbox
from gleam.utils.rgb_map import radial_mask
import gleam.utils.colors as gcl
gcl.GLEAMcmaps.register_all()

defaultrc = mpl.rcParams.copy()
mpl.rcParams = defaultrc

# mpl.rcParams.update(mpl.rcParamsDefault)
# mpl.rc('text', usetex=True)
# mpl.rc('font',**{'family':'serif','serif':['Palatino']})


VERBOSE = False
VIEWSTATE_PLOTS = False
PROFILES_PLOT = False
SYNTH_PLOTS = True

objs = ['B1608+656', 'DESJ0408-5354', 'HE0435-1223', 'PG1115+080',
        'RXJ0911+0551', 'RXJ1131-1231', 'SDSSJ1004+4112', 'WFIJ2033-4723']
sigf = [80, 100, 60, 600, 140, 4000, 80, 80]
# sigf = [20, 100, 12, 160, 35, 2000, 2, 140]
squarify = [0, 0, 0, 0, 0, 0, 0, 0]
lmaxf = [0.2500, 0.0075, 0.0075, 0.1000, 0.2000, 0.0100, 0.0100, 0.1000]
smaxf = [1., 10., 10., 0.5, 0.9, .75, 1., 1.]

fitsdir = 'data/delay_qsos/'
jsondir = 'jsons/'
statedir = 'states/'
statefiles = ['11doubles_dg45.state', '11doubles_dg60.state',
              '11doubles_CMB_dg60.state', '11doubles_SNeIa_dg60.state',
              '7quads_dg45.state', '7quads_CMB_dg45.state',
              '7quads_SNeIa_dg45.state', '7quads_dg60.state',
              '7quads_CMB_dg60.state', '7quads_SNeIa_dg60.state', 
              'all_dg60.state', 'all_SNeIa_dg60.state']
statefile = statefiles[7]  # 7  8/9 5/6
print(statefile)

lm = LensModel(statedir+statefile)

lens_idcs = [int(i) for i in sys.argv[1:]] or range(0, 8)

# for objidx in range(0, 8):
# for objidx in [0, 1, 2, 3, 4, 5, 7]:
# for objidx in [6,]:
# for objidx in [4, 5]:
for objidx in lens_idcs:
    lens = objs[objidx]
    print(lens)

    fitsfile = fitsdir + '{}.fits'.format(lens)
    if VERBOSE:
        print(fitsfile)
    jsonfile = jsondir + '{}.json'.format(lens)
    if VERBOSE:
        print(jsonfile)

    # LensObject
    with open(jsonfile) as f:
        lo = LensObject.from_json(f)
    lo.squarify(squarify[objidx])
    if VERBOSE:
        print(lo.__v__)

    # LensModel
    lm.obj_idx = objidx
    if VERBOSE:
        print(lm.__v__)

################################################################################
    if VIEWSTATE_PLOTS:
        mpl.rcParams['figure.figsize'] = 4, 4
        # Kappa
        if 0:  # for 4: RXJ0911
            flip_extent = list(lm.extent[0:2][::-1]) + list(lm.extent[2:])
            extnt = flip_extent
        else:
            extnt = lm.extent
        kappa_map_plot(lm, factor=lm.dlsds, obj_index=objidx, extent=extnt,
                       contours=True, levels=7, delta=0.1, label=lm.obj_name)
        plot_scalebar(R=lm.maprad, length=max(int(lm.maprad/2), 1))
        plt.axis('off')
        plt.gcf().axes[0].get_xaxis().set_visible(False)
        plt.gcf().axes[0].get_yaxis().set_visible(False)
        plt.savefig('results/{}_{}_kappa.pdf'.format(lm.obj_name, lm.filename.replace('.state', '')),
                    transparent=True, bbox_inches='tight', pad_inches=0)
        plt.close()

        # Kappa profile
        kappa_profiles_plot(lm, obj_index=objidx, ensemble_average=True,
                            refined=True, interpolate=150, levels=10,
                            as_range=True, maprad=lm.maprad, pixrad=lm.pixrad,
                            adjust_limits=True, annotation_color='white',
                            label_axes=True, fontsize=20)
        plt.tight_layout()
        plt.savefig('results/{}_{}_profiles.pdf'.format(lm.obj_name, lm.filename.replace('.state', '')),
                    transparent=True, bbox_inches='tight', pad_inches=0)
        plt.close()

        # Arrival-time surface
        if 0:  # for 4: RXJ0911
            lm.minima = lm.minima * np.array([-1, 1])
            saddles = lm.saddle_points*1
            lm.saddle_points = lm.saddle_points * np.array([-1, 1])
        else:
            saddles = None
        arrival_time_surface_plot(lm, geofactor=1./lm.dlsds, psifactor=1., obj_index=objidx, draw_images=True,
                                  contours=True, levels=60, sad_contour_saddles=saddles, extent=extnt,
                                  sad_contour_shift=None, scalebar=False,
                                  label=lm.obj_name, color='black')
        # plt.plot(0, 0, color=gcl.blue, marker='+', markersize=16)
        plot_scalebar(R=lm.maprad, length=max(int(lm.maprad/2), 1), color='black')
        plt.axis('off')
        plt.gcf().axes[0].get_xaxis().set_visible(False)
        plt.gcf().axes[0].get_yaxis().set_visible(False)
        plt.tight_layout()
        plt.savefig('results/{}_{}_arriv.pdf'.format(lm.obj_name, lm.filename.replace('.state', '')),
                    transparent=True, bbox_inches='tight', pad_inches=0)
        plt.close()

################################################################################
    if not SYNTH_PLOTS:
        continue
    # ReconSrc
    if objidx in [4, 5]:
        lo.data = np.fliplr(lo.data)
    Mo = 80
    if objidx in [5]:
        Mo = 250
    elif objidx in [6]:
        Mo = 600
    reconsrc = ReconSrc(lo, lm, M=Mo, M_fullres=400, mask_keys=['circle'],
                        zfactor=lm.dlsds)
    reconsrc.chmdl(-1)
    sig2 = sigf[objidx] * np.abs(reconsrc.lensobject.data)
    sig2[sig2 == 0] = sig2[sig2 != 0].min()

    if objidx in [5,]:
        # reconsrc.rotation = 115.2  # 1*reconsrc.lensobject.hdr['ORIENTAT']
        reconsrc.rotation = -98.2
    if objidx in [9,]:
        reconsrc.rotation = 1*reconsrc.lensobject.hdr['ORIENTAT']

    # dij = reconsrc.lens_map(mask=True)
    # plt.imshow(dij, extent=reconsrc.extent, origin='lower',
    #            vmin=0, vmax=lmaxf[objidx]*np.max(dij), interpolation='bicubic',
    #            cmap='gravic')
    # phcmap = plt.get_cmap('phoenix')
    # min_clr, sad_clr = phcmap(.2), phcmap(.5)
    # plt.plot(lm.minima.T[0], lm.minima.T[1], color=min_clr, marker='o', lw=0)
    # plt.plot(lm.saddle_points.T[0], lm.saddle_points.T[1], color=sad_clr, marker='o', lw=0)
    # plt.show()
    # exit(1)

    # PSF
    if objidx in [0, 1, 2, 7]:
        reconsrc.calc_psf("psf/tinytim_ACS.fits", normalize=True,
                          window_size=8, verbose=VERBOSE)
    elif objidx in [5,]:
        reconsrc.calc_psf("psf/tinytim_ACSF555W.fits", normalize=True,
                          window_size=8, verbose=VERBOSE)
    elif objidx in [3, 4, 6]:
        reconsrc.calc_psf("psf/tinytim_WFC3.fits", normalize=True,
                          window_size=8, verbose=VERBOSE)
    # reconsrc.calc_psf("psf/tinytim_SBC.fits", normalize=True,
    #                   window_size=8, verbose=True)


    # Calculating lens and source plane
    kw = dict(method='minres', use_psf=True, use_mask=True,
              use_filter=False, sigma2=sig2.copy(), cached=False)
    wrad = 0.8
    if objidx in [3, 4, 7]:
        wrad = 0.5
    reduced = False
    if objidx in [6]:
        reduced = True

    reconsrc.chmdl(-1)
    if objidx in [6]:
        reconsrc.inv_proj_matrix(use_mask=False, r_max=10)
    elif objidx in [7]:
        reconsrc.inv_proj_matrix(use_mask=False, r_max=2)
    else:
        reconsrc.inv_proj_matrix(use_mask=False)
    dij = reconsrc.lens_map(mask=True)
    src = reconsrc.plane_map(**kw)
    synth = reconsrc.reproj_map(from_cache=False, save_to_cache=False, **kw)
    residmap = reconsrc.residual_map(nonzero_only=True, within_radius=wrad,
                                     from_cache=False, save_to_cache=False,
                                     **kw)
    chi2 = reconsrc.reproj_chi2(reduced=reduced, nonzero_only=True,
                                within_radius=wrad, from_cache=False,
                                save_to_cache=False, **kw)
    print("Chi2: {}".format(chi2))
    if objidx in [4, 5]:
        dij = np.fliplr(dij)
        synth = np.fliplr(synth)
        residmap = np.fliplr(residmap)
        reconsrc.rotation = 115.2 - 98.2
        dij = ndimage.rotate(dij, reconsrc.rotation, reshape=False)
        synth = ndimage.rotate(synth, reconsrc.rotation, reshape=False)
        residmap = ndimage.rotate(residmap, reconsrc.rotation, reshape=False)

    # Plotting
    lmax = lmaxf[objidx] * np.max(dij)
    sbf = smaxf[objidx]*(reconsrc.lensobject.px2arcsec[0]/reconsrc.src_pxscale)**2

    # Lens data
    if objidx in [3]:
        dij = np.log10(dij+1)
        synth = np.log10(synth+1)
        lmax = np.log10(0.2*lmax+1)
        
    plt.imshow(dij, extent=reconsrc.extent, cmap='gravic', origin='lower',
               interpolation='bicubic', vmin=0, vmax=lmax)
    c = reconsrc.lensobject.roi.buffer.center
    c = (c - reconsrc.lensobject.center.xy) * reconsrc.lensobject.px2arcsec[0]
    plt.plot(*c, marker='+', markersize=10, color='grey')
    if objidx in [3, 4, 5]:
        plt.xlim(-3, 3)
        plt.ylim(-3, 3)
        pdng = (0.55, 0.55)
        if objidx == 5: pdng = (0.45, 0.45)
        plot_scalebar(R=reconsrc.lensobject.mapr, length=1,
                      padding=pdng, barheight=0.03*0.5)
    else:
        plot_scalebar(R=reconsrc.lensobject.mapr, length=1)
    plot_labelbox(lens, position='top left')
    plt.axis('off')
    plt.gcf().axes[0].get_xaxis().set_visible(False)
    plt.gcf().axes[0].get_yaxis().set_visible(False)
    plt.tight_layout()
    plt.savefig("results/{}_data.pdf".format(lens), transparent=True,
                bbox_inches='tight', pad_inches=0)
    plt.close()

    # Source plane
    plt.imshow(src, extent=reconsrc.src_extent, cmap='gravic', origin='lower',
               interpolation='bicubic', vmin=0, vmax=lmax*sbf/20)
    plot_scalebar(R=reconsrc.r_max, length=1)
    plot_labelbox(lens, position='top left')
    plt.axis('off')
    plt.gcf().axes[0].get_xaxis().set_visible(False)
    plt.gcf().axes[0].get_yaxis().set_visible(False)
    plt.tight_layout()
    # plt.savefig(
    #     "results/{}_{}_src.pdf".format(lens,
    #                                    lm.filename.replace('.state', '')),
    #     transparent=True, bbox_inches='tight', pad_inches=0)
    plt.close()

    # Synthetic image
    if objidx == 6:
        lmax = lmax * 0.5
    plt.imshow(synth, extent=reconsrc.extent, cmap='gravic', origin='lower',
               interpolation='bicubic', vmin=0, vmax=lmax)
    c = reconsrc.lensobject.roi.buffer.center
    c = (c - reconsrc.lensobject.center.xy) * reconsrc.lensobject.px2arcsec[0]
    plt.plot(*c, marker='+', markersize=10, color='grey')
    if objidx in [3, 4, 5]:
        plt.xlim(-3, 3)
        plt.ylim(-3, 3)
        pdng = (0.55, 0.55)
        if objidx == 5: pdng = (0.45, 0.45)
        plot_scalebar(R=reconsrc.lensobject.mapr, length=1,
                      padding=pdng, barheight=0.03*0.5)
    else:
        plot_scalebar(R=reconsrc.lensobject.mapr, length=1)
    plot_labelbox(lens, position='top left')
    plt.axis('off')
    plt.gcf().axes[0].get_xaxis().set_visible(False)
    plt.gcf().axes[0].get_yaxis().set_visible(False)
    plt.tight_layout()
    plt.savefig(
        "results/{}_{}_synth.pdf".format(lens,
                                         lm.filename.replace('.state', '')),
        transparent=True, bbox_inches='tight', pad_inches=0)
    plt.close()

    # Residual maps
    plt.imshow(residmap, extent=reconsrc.extent, cmap='vilux', origin='lower')
    c = reconsrc.lensobject.roi.buffer.center
    c = (c - reconsrc.lensobject.center.xy) * reconsrc.lensobject.px2arcsec[0]
    plt.plot(*c, marker='+', markersize=10, color='grey')
    if objidx in [3, 4, 5]:
        plt.xlim(-3, 3)
        plt.ylim(-3, 3)
        pdng = (0.55, 0.55)
        if objidx == 5: pdng = (0.45, 0.45)
        plot_scalebar(R=reconsrc.lensobject.mapr, length=1,
                      padding=pdng, barheight=0.03*0.5)
    else:
        plot_scalebar(R=reconsrc.lensobject.mapr, length=1)
    plot_labelbox(lens, position='top left')
    plt.axis('off')
    plt.gcf().axes[0].get_xaxis().set_visible(False)
    plt.gcf().axes[0].get_yaxis().set_visible(False)
    plt.tight_layout()
    # plt.savefig(
    #     "results/{}_{}_resid.pdf".format(lens,
    #                                      lm.filename.replace('.state', '')),
    #     transparent=True, bbox_inches='tight', pad_inches=0)
    plt.close()


if PROFILES_PLOT:
    fig, axes = plt.subplots(4, 2, sharex=False, sharey=True, figsize=(7, 11))
    for i in range(len(objs)):
        lm.obj_idx = i
        print(i)
        yextent = [0.65, 3.00]
        plt.sca(axes[i // 2][i % 2])
        plots, _, _ = kappa_profiles_plot(lm, obj_index=i, ensemble_average=True, refined=True,
                                          interpolate=150, levels=10, as_range=True, maprad=lm.maprad,
                                          pixrad=lm.pixrad, adjust_limits=True, annotation_color='white',
                                          einstein_radius_indicator=True, kappa1_line=True,
                                          label_axes=False, fontsize=22)
        xlim = list(plt.gca().axes.get_xlim())
        plt.contourf(np.ones((4, 4))*plots[0].levels[0], extent=xlim+yextent,
                     cmap='agaveglitch', levels=plots[0].levels, zorder=-99)
        plot_labelbox(objs[i], position='top right', padding=(0.03, 0.04), color='white',
                      fontsize=8, fontname='Computer Modern', family='sans-serif')
        plt.ylim(*yextent)
        # plt.rcParams['mathtext.fontset'] = 'stixsans'
        if (i % 2) == 0:
            if (i // 2) == 2:
                axes[i // 2][i % 2].set_ylabel(r'$\mathsf{\kappa}_{<\mathsf{R}}$', fontsize=22)
        if (i // 2) > 2:
            if (i % 2) == 0:
                axes[i // 2][i % 2].set_xlabel(r'R [arcsec]', fontsize=16)
        # else:
        #     axes[i // 2][i % 2].set_xticklabels([""]*len(axes[i // 2][i % 2].axes.get_xticklabels()))
        axes[i // 2][i % 2].set_frame_on(True)
    axes[3][0].xaxis.set_label_coords(1, -0.175)
    axes[2][0].yaxis.set_label_coords(-0.1, 1)
    fig.subplots_adjust(wspace=0.025)
    plt.savefig("results/profiles.pdf", transparent=True, bbox_inches='tight', pad_inches=0)
    plt.close()
    # plt.show()
