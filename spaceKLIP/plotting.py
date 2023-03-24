import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import astropy
import astropy.units as u
import jwst.datamodels

def annotate_compass(ax, image, wcs, xf=0.9, yf=0.1, length_fraction=0.07, color='white', fontsize=12, ):
    """ Plot a compass annotation onto an image, to indicate North and East

    Makes use of the methods from jdaviz, but positions the compass differently:
    jdaviz defaults to putting it in the center, and filling most of the image.
    Here we want to make a small compass in the corner.

    Parameters
    ----------
    ax : matplotlib.Axes
        axis to draw into
    image : ndarray
        2D image to be annotated (used just to get the image dimensions)
    wcs : astropy.wcs.WCS
        World Coordinate System information
    xf, yf : floats
        X and Y fractions of the image, to specify the location where the compass should be displayed.
        Values should be between 0 and 1.
    length_fraction : float
        Length of the compass, as a fraction of the size of the entire image
    color : str
        Color
    fontsize : float
        Font size

    """
    import jdaviz

    # Use jdaviz to compute arrow positions
    x, y, xn, yn, xe, ye, degn, dege, xflip = jdaviz.configs.imviz.wcs_utils.get_compass_info(wcs, image.shape,
                                                                                              r_fac=length_fraction)

    # but then apply offsets to recenter the
    xo = image.shape[1] * xf - x
    yo = image.shape[0] * yf - y
    x += xo
    xn += xo
    xe += xo
    y += yo
    yn += yo
    ye += yo

    # plot like in jdaviz:
    ax.plot(x, y, marker='o', color=color, markersize=4)
    ax.annotate('N', xy=(x, y), xytext=(xn, yn),
                arrowprops={'arrowstyle': '<-', 'color': color, 'lw': 1.5},
                color=color, fontsize=fontsize, va='center', ha='center')
    ax.annotate('E', xy=(x, y), xytext=(xe, ye),
                arrowprops={'arrowstyle': '<-', 'color': color, 'lw': 1.5},
                color=color, fontsize=fontsize, va='center', ha='center')


def annotate_scale_bar(ax, image, wcs, length=1 * u.arcsec, xf=0.1, yf=0.1, color='white', lw=3, fontsize=10):
    """Plot a scale bar on an image

    Parameters
    ----------
    ax : matplotlib.Axes
        axis to draw into
    image : ndarray
        2D image to be annotated (used just to get the image dimensions)
    wcs : astropy.wcs.WCS
        World Coordinate System information
    length : astropy.Quantity
        Length of the scale bar, in arcsec or equivalent unit
    xf, yf : floats
        X and Y fractions of the image, to specify the location where the compass should be displayed.
        Values should be between 0 and 1.
    color : str
        Color
    fontsize : float
        Font size
    lw : float
        line width

    """

    pixelscale = astropy.wcs.utils.proj_plane_pixel_scales(wcs).mean() * u.deg
    sb_length = (length / pixelscale).decompose()

    xo = image.shape[1] * xf
    yo = image.shape[0] * yf

    ax.plot([xo, xo + sb_length], [yo, yo], color=color, lw=lw)
    ax.text(xo + sb_length / 2, yo + 0.02 * image.shape[0], length, color=color,
            horizontalalignment='center', fontsize=fontsize)


def annotate_secondary_axes_arcsec(ax, image, wcs):
    """ Update an image display to add secondary axes labels in an arcsec

    Parameters
    ----------
    ax : matplotlib.Axes
        axis to draw into
    image : ndarray
        2D image to be annotated (used just to get the image dimensions)
    wcs : astropy.wcs.WCS
        World Coordinate System information
 
    """
    # define the forward and inverse transforms needed for secondary_axes.
    # see https://matplotlib.org/3.1.0/gallery/subplots_axes_and_figures/secondary_axis.html

    pixelscale = astropy.wcs.utils.proj_plane_pixel_scales(wcs)
    pix2as_x = lambda x: (x - wcs.wcs.crpix[0]) * pixelscale[0] * 3600
    pix2as_y = lambda y: (y - wcs.wcs.crpix[1]) * pixelscale[1] * 3600
    as2pix_x = lambda x: x / pixelscale[0] / 3600 + wcs.wcs.crpix[0]
    as2pix_y = lambda y: y / pixelscale[1] / 3600 + wcs.wcs.crpix[1]

    secax = ax.secondary_xaxis('top', functions=(pix2as_x, as2pix_x))
    secax.set_xlabel('Offset [arcsec]', fontsize='small')
    secay = ax.secondary_yaxis('right', functions=(pix2as_y, as2pix_y))
    secay.set_ylabel('Offset [arcsec]', fontsize='small')
    secax.tick_params(labelsize='small')
    secay.tick_params(labelsize='small')


def display_coron_image(filename):
    """Display and annotate a coronagraphic image

    Shows image on asinh scale along with some basic metadata, scale bar, and compass.

    Parameters
    ----------
    filename : str
        Filename
    """

    if ('uncal' in filename):
        raise RuntimeError("Display code does not support showing stage 0 uncal files. Reduce the data further before trying to display it.")
    elif ('rateints' in filename) or ('calints' in filename):
        modeltype = jwst.datamodels.CubeModel
        cube = True
    else:
        modeltype = jwst.datamodels.ImageModel
        cube = False

    model = modeltype(filename)  # cubemodel needed for rateints

    if cube:
        image = np.nanmean(model.data, axis=0)
        dq = model.dq[0]
    else:
        image = model.data
        dq = model.dq

    bpmask = np.zeros_like(image) + np.nan
    bpmask[(model.dq[0] & 1) == True] = 1

    # Set up image stretch
    #  including reasonable min/max for asinh stretch
    stats = astropy.stats.sigma_clipped_stats(image)
    low = stats[0] - stats[2]  # 1 sigma below image mean.
    high = np.nanmax(image)

    interval = astropy.visualization.ManualInterval(low, high)
    stretch = astropy.visualization.AsinhStretch(a=.0001)

    norm = astropy.visualization.ImageNormalize(image,
                                                interval=interval,
                                                stretch=stretch)

    # Display image. Overplot DQ
    fig, ax = plt.subplots(figsize=(16, 9))
    im = ax.imshow(image, norm=norm)

    imdq = ax.imshow(bpmask, vmin=0, vmax=1.5, cmap=matplotlib.cm.inferno)

    # Colorbar
    cb = fig.colorbar(im, pad=0.1, aspect=30, label=model.meta.bunit_data)
    cb.ax.set_yscale('asinh')

    # Annotations
    ax.text(0.01, 0.99,
            f"{model.meta.target.proposer_name}\n{model.meta.instrument.filter}, {model.meta.exposure.readpatt}:{model.meta.exposure.ngroups}:{model.meta.exposure.nints}\n{model.meta.exposure.effective_exposure_time:.2f} s",
            transform=ax.transAxes, color='white', verticalalignment='top', fontsize=10)
    ax.set_title(os.path.basename(filename) + "\n", fontsize=14)

    ax.set_xlabel("Pixels", fontsize='small')
    ax.set_ylabel("Pixels", fontsize='small')
    ax.tick_params(labelsize='small')

    if cube:
        ax.text(0.99, 0.99, f"Showing average of {model.meta.exposure.nints} ints",
                style='italic', fontsize=10, color='white',
                horizontalalignment='right', verticalalignment='top', transform=ax.transAxes)

    try:
        wcs = model.meta.wcs
        # I don't know how to deal with the slightly different API of the GWCS class
        # so, this is crude, just cast it to a regular WCS and drop the high order distortion stuff
        # This suffices for our purposes in plotting compass annotations etc.
        # (There is almost certainly a better way to do this...)
        wcs = astropy.wcs.WCS(model.meta.wcs.to_fits()[0])
    except:
        wcs = model.get_fits_wcs()
        if cube:
            wcs = wcs.dropaxis(2)  # drop the nints axis
    annotate_compass(ax, image, wcs, yf=0.07)
    annotate_scale_bar(ax, image, wcs, yf=0.07)

    # Annotate secondary axes in arcsec relative to coron mask center (really, relative to V2/V3Ref)
    annotate_secondary_axes_arcsec(ax, image, wcs)

    # TODO:
    #   add second panel with zoom in on center


def display_coron_dataset(database, restrict_to=None, save_filename=None):
    """ Display multiple files in a coronagraphic dataset

    Parameters
    -----------
    database : spaceklip.Database
        database of files to plot
    restrict_to : str, optional
        Optional query string to only display some data. Only datasets whose
        database concatenation (file group) name includes this string will be shown.
        Most simply, set this to a filter name to only plot images with that filter.
    save_filename : str
        If provided, the plots will be saved to a PDF file with this name.


    # TODO potentially provide other ways of filtering the data, e.g. to show
    only the PSF stars or only references, etc. 
    """
    if save_filename:
        from matplotlib.backends.backend_pdf import PdfPages
        pdf = PdfPages(save_filename)

    for key in database.obs:
        if (restrict_to is None) or (restrict_to in key):
            obstable = database.obs[key]
            for typestr in ['SCI', 'REF']:
                filenames = obstable[obstable['TYPE'] == typestr]['FITSFILE']

                for fn in filenames:
                    display_coron_image(fn)
                    if save_filename:
                        pdf.savefig(plt.gcf())
    if save_filename:
        pdf.close()
