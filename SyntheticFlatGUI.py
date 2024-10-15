import tkinter as tk
import tkinter.font
import tkinter.messagebox # for pyinstaller
from tkinter.filedialog import askopenfilename
from configparser import ConfigParser
import glob, os, sys, copy, numpy as np
import cv2
from rawpy import imread
from scipy.stats import sigmaclip
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter
import threading
import datetime as dt
from astropy.io import fits

GUINAME = "SyntheticFlatGUI"
VERSION = '2.0'

# STATIC SETTINGS =============================================================
IGNORE_EDGE = 10          # ignore pixels close to the image edges
IGNORE_RADII = 5          # ignore pixels extreme radii (close to 0 and maximum)
RADIAL_RESOLUTION = 1000000  # should be larger than 16 bit (65536)
DEBUG_MODE = False

RAWTYPES = ['arw', 'crw', 'cr2', 'cr3', 'nef', 'raf', 'rw2']
TIFTYPES = ['tif', 'tiff']
FITSTYPES = ['fit', 'fits']
OTHERTYPES = ['jpeg', 'jpg', 'png']
FILETYPES = RAWTYPES + TIFTYPES + FITSTYPES + OTHERTYPES

# STATIC MAJOR FUNCTIONS =====================================================
def load_image(file):

    # measure time
    start = dt.datetime.now()

    print("")
    print(file)
    if not os.path.isfile(file):
        raise ValueError("File does not exist!\n\n")

    # load image depending on type
    imagetype = os.path.basename(file).split(".")[1].lower()
    header = ""
    if imagetype in RAWTYPES:
        print("read raw image (" + imagetype.upper() + ") ...")
        im_load = imread(file).raw_image_visible
    elif imagetype in TIFTYPES:
        print("read tif type image (" + imagetype.upper() + ") ...")
        im_load = cv2.imread(file, flags=cv2.IMREAD_UNCHANGED)  # cv2.IMREAD_ANYDEPTH, cv2.IMREAD_UNCHANGED
    elif imagetype in FITSTYPES:
        print("read fits type image (" + imagetype.upper() + ") ...")
        im_load, header = fits.getdata(file, ext=0, header=True)
        print(" ".join(str(header).split()))
    elif imagetype in OTHERTYPES:
        print("read other type image (" + imagetype.upper() + ") ...")
        im_load = cv2.imread(file, flags=cv2.IMREAD_UNCHANGED)  # cv2.IMREAD_ANYDEPTH, cv2.IMREAD_UNCHANGED
    else:
        raise ValueError("image type not supported!!")

    # remember original shape
    origshape = im_load.shape
    print_image_info(im_load)

    # order axes
    im_deb = order_axes(im_load, type='HWC')

    # check colors and debayer, if no colors
    if len(im_deb.shape) == 3:
        print("no need to debayer ...")
    elif len(im_deb.shape) == 2:
        print("debayer ...")
        im_deb = debayer(im_load) # returns 4 colors!
    else:
        raise ValueError('Bad input shape')

    print_image_info(im_deb)

    # print execution time
    stop = dt.datetime.now()
    print("execution time:", int((stop - start).total_seconds() + 0.5), "seconds")

    return im_deb, origshape, header


def corr_gradient(image, resolution_factor=4, sigma_clip=3):

    # measure time
    start = dt.datetime.now()

    # measure slopes
    height, width, colors = image.shape
    slopes_x = []
    slopes_y = []
    for c in range(colors):
        if colors > 3 and c == 3:
            slopes_x.append(slopes_x[-1])
            slopes_y.append(slopes_y[-1])
            continue

        rowmeans = []
        rowselects = []
        for row in range(height):
            if row < IGNORE_EDGE or row > height - IGNORE_EDGE:
                continue
            if row % resolution_factor == 0:
                this_rowmean = sigma_clip_mean(image[row,IGNORE_EDGE:-IGNORE_EDGE, c], sigma_clip=sigma_clip)
                if this_rowmean:
                    rowmeans.append(this_rowmean)
                    rowselects.append(row)
        slope_y = np.mean(np.diff(rowmeans) / np.diff(rowselects)) * height
        slopes_y.append(slope_y)

        colmeans = []
        colselects = []
        for col in range(width):
            if col < IGNORE_EDGE or col > width - IGNORE_EDGE:
                continue
            if col % resolution_factor == 0:
                this_colmean = sigma_clip_mean(image[IGNORE_EDGE:-IGNORE_EDGE, col, c], sigma_clip=sigma_clip)
                if this_colmean:
                    colmeans.append(this_colmean)
                    colselects.append(col)
        slope_x = np.mean(np.diff(colmeans) / np.diff(colselects)) * width
        slopes_x.append(slope_x)

    # print
    print("gradient slopes x: ", slopes_x)
    print("gradient slopes y: ", slopes_y)

    for c in range(colors):
        x = np.linspace(0, 1, width)
        y = np.linspace(0, 1, height)
        X, Y = np.meshgrid(x, y)
        gradient = X * slopes_x[c] + Y * slopes_y[c] - (slopes_x[c] + slopes_y[c]) / 2
        image[:, :, c] = image[:, :, c] - gradient

    # print execution time
    stop = dt.datetime.now()
    print("execution time:", int((stop-start).total_seconds() + 0.5), "seconds")

    print_image_info(image)
    return image


def calc_histogram(image, circular=False):

    image = copy.deepcopy(image)

    # shape
    height, width, colors = image.shape
    centerdist_map = calc_centerdist_map(image.shape)

    for c in range(colors):
        if circular:
            mask = (centerdist_map > height / 2) | (centerdist_map > width / 2)
            image[mask, c] = np.nan

        # caLculate histogram
        counts, bins = np.histogram(
            image[:, :, c].flatten(),
            np.linspace(0, getmax(np.nanmax(image[:, :, c])), 2 ** 8)
        )
        bins = bins[1:]
        if c == 0:
            data = bins
        data = np.column_stack((data, counts))
    return data


def calc_rad_profile(image, statistics=2, extrapolate_max=True, resolution_factor=4):

    # measure time
    start = dt.datetime.now()

    # shape
    height, width, colors = image.shape
    centerdist_map = calc_centerdist_map(image.shape)

    # acquire pixel values
    radii = []
    rad_counts = {}
    for i in range(height):
        for j in range(width):
            if not (i % resolution_factor == 0 and j % resolution_factor == 0):
                continue
            if i < IGNORE_EDGE or j < IGNORE_EDGE or i > height-IGNORE_EDGE or j > width-IGNORE_EDGE:
                continue
            if not np.min(image[i, j, :] > 0):
                continue
            rad = int(centerdist_map[i, j])
            if not rad in radii:
                radii.append(rad)
                rad_counts[rad] = []
                for c in range(colors):
                    rad_counts[rad].append([])
            for c in range(colors):
                rad_counts[rad][c].append(image[i, j, c])
    maxrad = centerdist_map[0, 0]
    rad_profile = np.zeros((len(radii), colors + 1))
    rad_profile_raw_mean = np.zeros((len(radii), colors + 1))
    index = 0
    for rad in sorted(radii):
        rad_profile[index, 0] = rad / maxrad
        for c in range(colors):
            rad_profile[index, c + 1] = apply_statistics(rad_counts[rad][c], statistics)
        rad_profile_raw_mean[index, 0] = rad / maxrad
        for c in range(colors):
            rad_profile_raw_mean[index, c + 1] = np.mean(rad_counts[rad][c])
        index += 1
    rad_profile = rad_profile[~np.isnan(rad_profile).any(axis=1), :]
    print("rad_profile shape: ", rad_profile.shape)

    # safety
    if not rad_profile.shape[0] > IGNORE_RADII * 2:
        print("\nSomething ist wrong!!\nInterrupt!\n\n")

    # cut edges
    if IGNORE_RADII:
        mask = rad_profile[:, 0] > IGNORE_RADII / maxrad
        rad_profile_cut = rad_profile[mask, :]
        mask = rad_profile_cut[:, 0] < 1 - IGNORE_RADII / maxrad
        rad_profile_cut = rad_profile_cut[mask, :]
    else:
        rad_profile_cut = copy.deepcopy(rad_profile)

    # cut data inside maximum
    if extrapolate_max:
        maxind = []
        for c in range(colors):
            this_filtered = rad_profile_cut[:, c + 1]
            this_filtered = savgol_filter(this_filtered, window_length=odd_int(rad_profile_cut.shape[0] / 10),
                                                       polyorder=2, mode='interp')
            this_filtered = savgol_filter(this_filtered, window_length=odd_int(rad_profile_cut.shape[0] / 5),
                                                       polyorder=2, mode='interp')
            mi = np.argmax(this_filtered)
            maxind.append(mi)
        print("maximum index: ", maxind)
        if np.max(maxind) > int(rad_profile_cut.shape[0]/2):
            print("maximum index too large -> skip maximum cut")
        else:
            rad_profile_cut = rad_profile_cut[max(maxind):, :]

    # smooth data and interpolate
    radii = np.linspace(0, 1, RADIAL_RESOLUTION)
    rad_profile_smoothed = radii
    slopes_inner = []
    slopes_outer = []
    for c in range(colors):
        y_new = savgol_filter(rad_profile_cut[:, c + 1], window_length=odd_int(rad_profile_cut.shape[0] / 10), polyorder=2, mode='interp')
        y_new = savgol_filter(y_new, window_length=odd_int(rad_profile_cut.shape[0] / 5), polyorder=2, mode='interp')

        # extrapolate inner
        x_new = list(rad_profile_cut[:, 0].flatten())
        y_new = list(y_new.flatten())
        xdist = x_new[1] - x_new[0]
        ylast = y_new[0]
        xlast = x_new[0]
        slope = (y_new[1] - y_new[0]) / xdist
        if slope > 0:
            slope = 0
        slopes_inner.append(int(slope))
        xadd = np.arange(0, xlast - xdist, xdist)
        # quadratic function with given value and slope at xlast and zero derivative at x=0
        yadd = ylast + slope / 2 * (xadd ** 2 / xlast - xlast)
        x_new = list(xadd) + x_new
        y_new = list(yadd) + y_new

        # extrapolate outer
        x_new = x_new[:-4]  # whatever reason, otherwise step in profile
        y_new = y_new[:-4]  # whatever reason, otherwise step in profile
        xdist = x_new[-1] - x_new[-2]
        ylast = y_new[-1]
        xlast = x_new[-1]
        slope = (y_new[-1] - y_new[-2]) / xdist
        slopes_outer.append(int(slope))
        xadd = np.arange(xlast + xdist, 1, xdist)
        # linear function
        yadd = ylast + slope * (xadd - xlast)
        x_new = x_new + list(xadd)
        y_new = y_new + list(yadd)

        # interpolate onto continuous values
        f = interp1d(x_new, y_new, kind='quadratic', fill_value='extrapolate')
        y_new = f(radii)
        rad_profile_smoothed = np.column_stack((rad_profile_smoothed, y_new))

    # print extrapolation info
    print("extrapolated with inner slopes: ", slopes_inner)
    print("extrapolated with outer slopes: ", slopes_outer)

    # normalize by smoothed curve
    for c in range(colors):
        rad_profile_raw_mean[:, c + 1] = rad_profile_raw_mean[:, c + 1] / np.max(rad_profile_smoothed[:, c + 1])
        rad_profile[:, c + 1] = rad_profile[:, c + 1] / np.max(rad_profile_smoothed[:, c + 1])
        rad_profile_cut[:, c + 1] = rad_profile_cut[:, c + 1] / np.max(rad_profile_smoothed[:, c + 1])
        rad_profile_smoothed[:, c + 1] = rad_profile_smoothed[:, c + 1] / np.max(rad_profile_smoothed[:, c + 1])

    # print minimum value
    print("minimum:", np.min(rad_profile_smoothed[:, 1:]))

    # print execution time
    stop = dt.datetime.now()
    print("execution time:", int((stop-start).total_seconds() + 0.5), "seconds")

    return rad_profile_raw_mean, rad_profile, rad_profile_cut, rad_profile_smoothed


def calc_synthetic_flat(rad_profile, grey_flat=False, out_size=(4024, 6024)):

    # measure time
    start = dt.datetime.now()

    # shape, prepare
    height, width = out_size[:2]
    centerdist_map = calc_centerdist_map(out_size)
    im_syn = np.zeros(out_size)

    # normalize (again)
    colors_in = rad_profile.shape[1] - 1
    for c in range(colors_in):
        rad_profile[:, c + 1] = rad_profile[:, c + 1] / np.max(rad_profile[:, c + 1])

    # convert centerdist_map into indices of radial profile
    r_indices = (RADIAL_RESOLUTION - 1) * centerdist_map / centerdist_map[0, 0]
    r_indices = r_indices.astype(int)

    # 3D output
    if len(out_size) >= 3:
        for c in range(out_size[2]):
            if grey_flat:
                # write green flat in each color
                im_syn[:, :, c] = rad_profile[r_indices, 2]
            else:
                # write separate flat in each color
                im_syn[:, :, c] = rad_profile[r_indices, c + 1]

    # 2D output
    else:
        if grey_flat:
            # write green flat once
            im_syn = rad_profile[r_indices, 2]
        else:
            # write separate flats in each color and bayer afterwards
            im_syn = np.zeros([out_size[0], out_size[1], colors_in])
            for c in range(colors_in):
                im_syn[:, :, c] = rad_profile[r_indices, c + 1]
            im_syn = bayer(im_syn, keep_size=True)

    # print flat info
    if np.count_nonzero(im_syn == 0) > 0:
        print("zeros (%):", "{:.10f}".format(100.0 - 100.0 * np.count_nonzero(im_syn) / np.prod(im_syn.shape)))
    if np.count_nonzero(np.isnan(im_syn)) > 0:
        print("nan (%):", "{:.10f}".format(100.0 * np.count_nonzero(np.isnan(im_syn)) / np.prod(im_syn.shape)))

    # print execution time
    stop = dt.datetime.now()
    print("execution time:", int((stop-start).total_seconds() + 0.5), "seconds")

    print_image_info(im_syn)
    return im_syn



# STATIC MINOR FUNCTIONS =====================================================

def print_image_info(image):
    minmedmax = str(np.min(image)) + " > " + str(np.median(image)) + " > " + str(np.max(image))
    uniques = "d" + str(len(np.unique(image))) + "/" + str(int(np.max(image) / np.min(np.diff(sorted(np.unique(image))))))
    print(image.shape, uniques, minmedmax, sep=' | ')

def get_axes_order(shape):
    first_axis = -1
    type = ""
    for c in range(len(shape)):
        if shape[c] < 5:
            type += "C"
        elif first_axis > 0:
            if shape[c] < first_axis:
                type += "H"
                type = type.replace("X", "W")
            else:
                type += "W"
                type = type.replace("X", "H")
        else:
            type += "X"
            first_axis = shape[c]
    return type

def order_axes(image, type='HWC'):
    image_shape = image.shape
    axes_sort = np.argsort(image_shape)[::-1]
    image = np.transpose(image, axes=axes_sort)
    if type == 'HWC':
        image = np.swapaxes(image, 0, 1)
    if len(image.shape) == 3 and type == 'CHW':
        image = np.swapaxes(image, 0, 2)
    return image

def get_depth(image):
    if np.max(image) > 2 ** 16 - 1:
        depth = 32
    elif np.max(image) > 2 ** 8 - 1:
        depth = 16
    elif np.max(image) > 1:
        depth = 8
    else:
        depth = 0
    print("image depth:", depth)
    return depth

def set_depth(image, depth):
    depth_image = copy.deepcopy(image)

    # check and print
    if not depth:
        print("depth ist None...")
        return depth_image
    elif depth > 0:
        print("set depth uint" + str(depth) + " (max " + str(2 ** depth - 1) + ")")
    else:
        print("set depth float 32" + " (max " + str(2 ** 32 - 1) + ")")

    # if input image is 0-1 and depth is integer, multiply
    if np.max(depth_image) <= 1 and depth > 0:
        depth_image = (2 ** depth -1) * depth_image
    else:
        depth_image = depth_image

    # clip
    if depth > 0:
        clip_value = 2 ** depth - 1
    else:
        clip_value = 1
    depth_image = np.clip(depth_image, 0, clip_value)

    # set image depth
    if depth == 32:
        depth_image = depth_image.astype(np.uint32)
    elif depth == 16:
        depth_image = depth_image.astype(np.uint16)
    elif depth == 8:
        depth_image = depth_image.astype(np.uint8)
    else:
        depth_image = depth_image.astype(np.float32)

    # check if some nonzeros have been set to zero
    set_zero = image.flatten()[(image.flatten() > 0) & (depth_image.flatten() == 0)]
    if len(set_zero) > 0:
        print(">>> Some values were set to zero!!")
        print(set_zero)
    return depth_image

def getmax(maxval):
    if maxval <= 1:
        return 1
    else:
        vallog = int(np.log2(maxval)) + 1
        if vallog % 2 == 1:
            vallog += 1
    return 2 ** vallog - 1

def debayer(image):
    rows = image.shape[0]
    cols = image.shape[1]
    db_image = np.zeros((int(rows / 2), int(cols / 2), 4))
    db_image[:,:,0] = image[0::2, 0::2]
    db_image[:,:,1] = image[0::2, 1::2]
    db_image[:,:,2] = image[1::2, 0::2]
    db_image[:,:,3] = image[1::2, 1::2]
    return db_image


def bayer(image, keep_size=False):
    if len(image.shape) < 3:
        return image
    rows, cols, colors = image.shape
    print("bayer " + str(colors) + " colors image" + " (keep size)" * keep_size)
    if colors == 3:
        color_order = [0, 1, 1, 2]
    else:
        color_order = [0, 1, 2, 3]
    if keep_size:
        b_image = np.zeros((rows, cols))
        b_image[0::2, 0::2] = image[0::2, 0::2, color_order[0]]
        b_image[0::2, 1::2] = image[0::2, 1::2, color_order[1]]
        b_image[1::2, 0::2] = image[1::2, 0::2, color_order[2]]
        b_image[1::2, 1::2] = image[1::2, 1::2, color_order[3]]
    else:
        b_image = np.zeros((rows * 2, cols * 2))
        b_image[0::2, 0::2] = image[:,:,color_order[0]]
        b_image[0::2, 1::2] = image[:,:,color_order[1]]
        b_image[1::2, 0::2] = image[:,:,color_order[2]]
        b_image[1::2, 1::2] = image[:,:,color_order[3]]

    return b_image


def calc_centerdist_map(shape):
    height, width = shape[:2]
    xax = np.arange(width)
    yax = np.arange(height)
    X, Y = np.meshgrid(xax, yax)
    X = X - width  / 2 + 0.5
    Y = Y - height / 2 + 0.5
    distances = np.sqrt(np.power(X, 2) + np.power(Y, 2))
    return distances


def create_folder(folder_path):
    try:
        os.mkdir(folder_path)
    except:
        pass
    return folder_path


def odd_int(number):
    if int(number) % 2 == 1:
        return int(number)
    else:
        return int(number) + 1


def contains(string, substrings):
    if not isinstance(substrings, list):
        substrings = [substrings]
    for sub in substrings:
        if string.lower().find(sub.lower()) >= 0:
            return True
    return False


def rgbtohex(r, g, b):
    return f'#{r:02x}{g:02x}{b:02x}'


def apply_statistics(array, statistics):
    if statistics[:10] == 'sigma clip':
        sigma_clip = float(statistics[10:].strip())
        result = sigma_clip_mean(array, sigma_clip=sigma_clip)
    elif statistics == 'median':
        result = np.median(array)
    elif statistics == 'min':
        result = np.min(array)
    elif statistics == 'max':
        result = np.max(array)
    else:
        result = np.mean(array)
    return result


def sigma_clip_mean(array, sigma_clip=2.0):
    reduced = sigmaclip(array, low=sigma_clip, high=sigma_clip)[0]
    if not reduced.any():
        result = None
    else:
        result = np.mean(reduced)
    return result


def write_csv(data, savepath, original_file, suffix):
    columns = data.shape[1]
    fmt = ['%.10f'] * columns
    # fmt = '%.10f', '%.10f', '%.10f', '%.10f'
    np.savetxt(savepath + os.sep + os.path.basename(original_file).split('.')[0] + suffix + ".csv",
               data, delimiter=",", fmt=fmt)


# IMAGE CLASS =====================================================
class Image():
    def __init__(self, file, debug=False):

        # init attributes
        self.file = file
        self.debug = debug
        self.image = None
        self.header = ''
        self.origpath = os.path.dirname(file)
        self.origtype = os.path.basename(self.file).split(".")[1].lower()
        self.origshape = None
        self.origdepth = None
        self.outpath = ''
        self.outpath_csv = ''
        self.outtype = ''
        self.radprof = None
        self.image_flat = None

        # set paths
        self.origpath = os.path.dirname(self.file)
        self.outpath = create_folder(self.origpath + os.sep + GUINAME)
        self.outpath_csv = create_folder(self.outpath + os.sep + "csv")

        self.set_debug()

    def set_debug(self):
        if not self.debug:
            return
        self.file = self.file.split('.')[0] + '.tif' # fake some ending
        self.image = np.ones((400, 600, 3))
        self.origshape = (400, 600, 3)
        radprof_x = np.linspace(0, 1, RADIAL_RESOLUTION)
        radprof_y = np.linspace(1, 0.7, RADIAL_RESOLUTION)
        self.radprof = np.column_stack((radprof_x, radprof_y, radprof_y, radprof_y))
        self.set_paths_types()

    def load(self):
        if self.debug:
            return
        self.image, self.origshape, self.header = load_image(self.file)

        # origtype, outtype
        self.origtype = os.path.basename(self.file).split(".")[1].lower()
        if self.origtype in RAWTYPES:
            self.outtype = 'tif'
        else:
            self.outtype = self.origtype

        # origdepth
        self.origdepth = get_depth(self.image)

        # float conversion necessary for unclipped calculations
        # handle clipping and depth upon writing images
        self.image = self.image.astype(np.float64)

    def write_image(self, suffix, flat=False):
        print("write image \"" + suffix + "\"")
        newfile = self.outpath + os.sep + os.path.basename(self.file).split('.')[0] + suffix + "." + self.outtype
        print(newfile)

        # choose
        if flat:
            image_write = self.image_flat
        else:
            image_write = self.image

        # the image is certainly debayered now
        # -> if original was bayered, bayer it again before writing
        if len(self.origshape) == 2:
            image_write = bayer(image_write)
        else:
            image_write = order_axes(image_write, get_axes_order(self.origshape))

        # write
        #image_write = image_write[18:22, 28:32, :]
        #print(image_write[:,:,0])
        image_write = set_depth(image_write, self.origdepth)
        print_image_info(image_write)
        if self.outtype in FITSTYPES:
            fits.PrimaryHDU(image_write, header=self.header).writeto(newfile, overwrite=True)
        else:
            cv2.imwrite(newfile, image_write)

    def gradcorr(self, resolution_factor):
        if self.debug:
            return
        self.image = corr_gradient(self.image, resolution_factor=resolution_factor)

    def subtract_bias(self, bias_value):
        if self.debug:
            return
        self.image = self.image - bias_value

    def calc_histogram(self, circular=False):
        if self.debug:
            return
        data = calc_histogram(
            self.image,
            circular=circular,
        )
        write_csv(data, self.outpath_csv, self.file, "_histogram")

    def calc_rad_profile(self, statistics=2, extrapolate_max=True, resolution_factor=4):
        if self.debug:
            return
        radprof1, radprof2, radprof3, radprof4 = calc_rad_profile(
            self.image,
            statistics=statistics,
            extrapolate_max=extrapolate_max,
            resolution_factor=resolution_factor,
        )
        write_csv(radprof1, self.outpath_csv, self.file, "_radprof_0_raw_mean")
        write_csv(radprof2, self.outpath_csv, self.file, "_radprof_1_clipped")
        write_csv(radprof3, self.outpath_csv, self.file, "_radprof_2_cut")
        write_csv(radprof4, self.outpath_csv, self.file, "_radprof_3_smooth")
        self.radprof = radprof4

    def calc_synthetic_flat(self, grey_flat=False):
        self.image_flat = calc_synthetic_flat(
            self.radprof,
            grey_flat=grey_flat,
            out_size=self.image.shape,
        )

    def apply_synthetic_flat(self, clip=True):
        if clip:
            clip_val = np.max(self.image)
        self.image = self.image / self.image_flat
        self.image[self.image > clip_val] = clip_val


# TKINTER CLASS =====================================================

class NewGUI():
    def __init__(self):
        self.root = tk.Tk()
        self.root.title(GUINAME + " v" + VERSION)
        self.lastpath = ""
        self.loaded_files = []
        self.bias_value = 0
        self.running = False
        self.asked_stop = False

        padding = 5

        self.root.protocol("WM_DELETE_WINDOW",  self.on_close)

        # icon and DPI
        try:
            self.root.iconbitmap(GUINAME + ".ico")
            self.root.update() # important: recalculate the window dimensions
        except:
            print("Found no icon.")

        # menu bar
        menubar = tk.Menu(self.root)
        self.root.config(menu=menubar)

        # mainoptions
        options = tk.Menu(menubar, tearoff=0)
        self.opt_gradient  = tk.BooleanVar()
        self.opt_histogram = tk.BooleanVar()
        self.opt_radprof   = tk.BooleanVar()
        self.opt_synthflat = tk.BooleanVar()
        options.add_checkbutton(label="Correct gradient", onvalue=1, offvalue=0, variable=self.opt_gradient)
        options.add_checkbutton(label="Calculate histogram", onvalue=1, offvalue=0, variable=self.opt_histogram)
        options.add_checkbutton(label="Calculate radial profile", onvalue=1, offvalue=0, variable=self.opt_radprof, command=self.toggle_radprof)
        options.add_checkbutton(label="Export synthetic flat", onvalue=1, offvalue=0, variable=self.opt_synthflat, command=self.toggle_synthflat)
        menubar.add_cascade(label="Options", menu=options)

        # settings
        settings = tk.Menu(menubar, tearoff=0)
        self.set_export_corr_input = tk.BooleanVar()
        self.set_circular_hist   = tk.BooleanVar()
        self.set_grey_flat       = tk.BooleanVar()
        self.set_extrapolate_max = tk.BooleanVar()
        settings.add_checkbutton(label="Histogram of largest circle", onvalue=1, offvalue=0, variable=self.set_circular_hist)
        settings.add_checkbutton(label="Extrapolate inside max", onvalue=1, offvalue=0, variable=self.set_extrapolate_max)
        settings.add_checkbutton(label="Export corrected input images", onvalue=1, offvalue=0, variable=self.set_export_corr_input)
        settings.add_checkbutton(label="Grey synthetic flat", onvalue=1, offvalue=0, variable=self.set_grey_flat)
        menubar.add_cascade(label="Settings", menu=settings)

        # statistics
        statistics = tk.Menu(menubar, tearoff=0)
        self.radio_statistics = tk.StringVar(self.root)
        self.sigmas = ["mean", "median", "min", "max", "sigma clip 0.5", "sigma clip 1.0",
                       "sigma clip 2.0", "sigma clip 3.0", "sigma clip 4.0", "sigma clip 8.0"]
        for opt in self.sigmas:
            statistics.add_radiobutton(label=opt, value=opt, variable=self.radio_statistics)
        menubar.add_cascade(label="Statistics", menu=statistics)

        # resolution
        resolution = tk.Menu(menubar, tearoff=0)
        self.radio_resolution = tk.StringVar(self.root)
        self.resolution_factor = ['full', '1/2', '1/4', '1/8', '1/16']
        for opt in self.resolution_factor:
            resolution.add_radiobutton(label=opt, value=opt, variable=self.radio_resolution)
        menubar.add_cascade(label="Resolution", menu=resolution)

        # Reset
        menubar.add_command(label="Reset", command=self.reset_config)

        # buttons
        self.button_load = tk.Button(text="Load files", command=self.load_files)
        self.button_load.grid(row=0, column=0, sticky='NWSE', padx=padding, pady=padding, columnspan=2)
        self.button_load = tk.Button(text="Set bias value", command=self.ask_bias)
        self.button_load.grid(row=1, column=0, sticky='NWSE', padx=padding, pady=padding)
        self.button_load = tk.Button(text="Bias from file", command=self.ask_bias_file)
        self.button_load.grid(row=1, column=1, sticky='NWSE', padx=padding, pady=padding)
        self.button_start = tk.Button(text="Start", command=lambda: threading.Thread(target=self.process).start())
        self.button_start.grid(row=2, column=0, sticky='NWSE', padx=padding, pady=padding)
        self.button_stop = tk.Button(text="Stop", command=self.stop)
        self.button_stop.grid(row=2, column=1, sticky='NWSE', padx=padding, pady=padding)

        # labels
        self.label_files_var = tk.StringVar()
        self.label_files_var.set("0 files")
        self.label_files = tk.Label(textvariable=self.label_files_var, justify='left')
        self.label_files.grid(row=0, column=2, sticky='NWSE', padx=padding, pady=padding)
        self.label_bias_var = tk.StringVar()
        self.label_bias_var.set(0)
        self.label_bias = tk.Label(textvariable=self.label_bias_var, justify='left')
        self.label_bias.grid(row=1, column=2, sticky='NWSE', padx=padding, pady=padding)
        self.label_status_var = tk.StringVar()
        self.label_status_var.set("ready")
        self.label_status = tk.Label(textvariable=self.label_status_var, justify='left')
        self.label_status.grid(row=2, column=2, sticky='NWSE', padx=padding, pady=padding)
        self.update_labels()

        # configure
        self.root.grid_columnconfigure(0, weight=1)
        self.root.grid_columnconfigure(1, weight=1)
        self.root.grid_columnconfigure(2, weight=3)
        for i in range(3):
            self.root.grid_rowconfigure(i, weight=1)

        # default configs
        self.load_config_file()

        # mainloop
        self.root.mainloop()

    def stop(self):
        print("asked_stop:", self.asked_stop, end=' ')
        self.asked_stop = True
        print("to", self.asked_stop)
        self.update_labels(status="stopping...")

    def check_stop(self):
        if self.asked_stop:
            self.running = False
            self.update_labels(status="interrupted")
            print("\nInterrupted!\n")
            self.asked_stop = False
            return True
        else:
            return False

    def on_close(self):
        print("... save config file")
        config_object = ConfigParser()

        config_object["BASICS"] = {}
        config_object["BASICS"]["window size"]      = self.root.winfo_geometry()
        config_object["BASICS"]["lastpath"]         = self.lastpath
        config_object["BASICS"]["radio_statistics"] = self.radio_statistics.get()
        config_object["BASICS"]["radio_resolution"] = self.radio_resolution.get()
        config_object["BASICS"]["bias_value"]   = str(self.bias_value)

        config_object["OPTIONS"] = {}
        config_object["OPTIONS"]["opt_gradient"]  = str(self.opt_gradient.get())
        config_object["OPTIONS"]["opt_histogram"] = str(self.opt_histogram.get())
        config_object["OPTIONS"]["opt_radprof"]   = str(self.opt_radprof.get())
        config_object["OPTIONS"]["opt_synthflat"] = str(self.opt_synthflat.get())

        config_object["SETTINGS"] = {}
        config_object["SETTINGS"]["set_export_corr_input"]  = str(self.set_export_corr_input.get())
        config_object["SETTINGS"]["set_circular_hist"]    = str(self.set_circular_hist.get())
        config_object["SETTINGS"]["set_grey_flat"]        = str(self.set_grey_flat.get())
        config_object["SETTINGS"]["set_extrapolate_max"]  = str(self.set_extrapolate_max.get())

        with open(GUINAME + ".conf", 'w') as conf:
            config_object.write(conf)

        self.root.destroy()

    def reset_config(self, reset_window=False):
        config_object = ConfigParser()
        config_object["BASICS"] = {}
        if reset_window:
            config_object["BASICS"]["window size"] = '318x128+313+94'
        else:
            config_object["BASICS"]["window size"] = self.root.winfo_geometry()
        config_object["BASICS"]["lastpath"] = './'
        config_object["BASICS"]["radio_statistics"] = 'sigma clip 2.0'
        config_object["BASICS"]["radio_resolution"] = '1/4'
        config_object["BASICS"]["bias_value"] = '0'

        config_object["OPTIONS"] = {}
        config_object["OPTIONS"]["opt_gradient"] = 'True'
        config_object["OPTIONS"]["opt_histogram"] = 'True'
        config_object["OPTIONS"]["opt_radprof"] = 'True'
        config_object["OPTIONS"]["opt_synthflat"] = 'True'

        config_object["SETTINGS"] = {}
        config_object["SETTINGS"]["set_export_corr_input"] = 'True'
        config_object["SETTINGS"]["set_circular_hist"] = 'True'
        config_object["SETTINGS"]["set_grey_flat"] = 'True'
        config_object["SETTINGS"]["set_extrapolate_max"] = 'True'

        self.apply_config(config_object)
        self.update_labels(status="ready")

    def apply_config(self, config_object):
        self.root.geometry(config_object["BASICS"]["window size"])
        self.lastpath = config_object["BASICS"]["lastpath"]
        self.radio_statistics.set(config_object["BASICS"]["radio_statistics"])
        self.radio_resolution.set(config_object["BASICS"]["radio_resolution"])
        self.bias_value = int(config_object["BASICS"]["bias_value"])

        self.opt_gradient.set(config_object["OPTIONS"]["opt_gradient"]   == 'True')
        self.opt_histogram.set(config_object["OPTIONS"]["opt_histogram"] == 'True')
        self.opt_radprof.set(config_object["OPTIONS"]["opt_radprof"]     == 'True')
        self.opt_synthflat.set(config_object["OPTIONS"]["opt_synthflat"] == 'True')

        self.set_export_corr_input.set(config_object["SETTINGS"]["set_export_corr_input"]   == 'True')
        self.set_circular_hist.set(config_object["SETTINGS"]["set_circular_hist"] == 'True')
        self.set_grey_flat.set(config_object["SETTINGS"]["set_grey_flat"]         == 'True')
        self.set_extrapolate_max.set(config_object["SETTINGS"]["set_extrapolate_max"] == 'True')

        self.update_labels()

    def toggle_radprof(self):
        if not self.opt_radprof.get():
            self.opt_synthflat.set(0)

    def toggle_synthflat(self):
        if self.opt_synthflat.get():
            self.opt_radprof.set(1)

    def load_config_file(self):
        # read
        if os.path.exists(GUINAME + ".conf"):
            config_object = ConfigParser()
            config_object.read(GUINAME + ".conf")
            self.apply_config(config_object)

        # default
        else:
            self.reset_config(reset_window=True)

    def load_files(self):
        if self.running:
            return
        user_input = askopenfilename(initialdir=self.lastpath, multiple=True,
              filetypes=[('Images', ";".join(["*." + x for x in FILETYPES])), ('all', '.*')])
        if user_input:
            self.loaded_files = user_input
            self.update_labels(file=str(len(self.loaded_files)) + " files")
            self.update_labels(status="ready")
            self.lastpath = os.path.dirname(self.loaded_files[0])
        return

    def ask_bias(self):
        if self.running:
            return
        user_input = tk.simpledialog.askinteger(title="Bias value", prompt="Which bias value should be subtracted from the image?")
        self.bias_value = int(user_input)
        self.label_bias_var.set(self.bias_value)
        self.update_labels()
        return

    def ask_bias_file(self):
        if self.running:
            return
        self.update_labels(status="calc bias...")
        user_input_file = askopenfilename(initialdir=self.lastpath, multiple=True,
              filetypes=[('RAW format (supported)', RAW_TYPES), ('Image format (not supported)', IMAGE_TYPES), ('all', '.*')])
        im_raw = imread(user_input_file[0]).raw_image_visible
        self.bias_value = int(sigma_clip_mean(im_raw))
        self.label_bias_var.set(self.bias_value)
        self.update_labels(status="ready")
        return

    def update_labels(self, file="", status=""):
        if file:
            self.label_files_var.set(file)
        if status:
            self.label_status_var.set(status)
            print("\n>> ", status)
        if contains(self.label_status_var.get().lower(), ["ready", "finish"]):
            self.label_status.configure(background=rgbtohex(180, 230, 180))
        elif contains(self.label_status_var.get().lower(), ["error", "interr", "stop", "no file"]):
            self.label_status.configure(background=rgbtohex(250, 180, 180))
        elif contains(self.label_status_var.get().lower(), ["write"]):
            self.label_status.configure(background=rgbtohex(240, 230, 180))
        else:
            self.label_status.configure(background=rgbtohex(240, 210, 180))
        if len(self.loaded_files) > 0:
            self.label_files.configure(background=rgbtohex(210, 230, 255))
        else:
            self.label_files.configure(background=rgbtohex(250, 180, 180))
        if self.bias_value > 0:
            self.label_bias.configure(background=rgbtohex(210, 230, 255))
        self.label_bias_var.set(self.bias_value)
        self.root.update()
        return

    def process(self):
        if self.running:
            return
        else:
            self.running = True
            self.update_labels(status="running...")
        self.asked_stop = False
        try:
            counter = 0

            # safety check
            if not self.loaded_files:
                self.running = False
                self.update_labels(status="no file chosen.")
                return

            # resolution factor
            if self.radio_resolution.get()[2:].isnumeric():
                resolution_factor = int(self.radio_resolution.get()[2:])
            else:
                resolution_factor = 1

            for file in self.loaded_files:

                # set and display
                counter += 1
                self.update_labels(file=os.path.basename(file) + " (" + str(counter) + "/" + str(len(self.loaded_files)) + ")")

                # load
                self.update_labels(status="load...")
                if DEBUG_MODE:
                    imobj = Image(file, debug=True)
                else:
                    imobj = Image(file)
                imobj.load()
                if self.check_stop(): return

                # write original image
                if self.set_export_corr_input.get():
                    self.update_labels(status="write original...")
                    imobj.write_image("_0_input")
                    if self.check_stop(): return

                # subtract bias
                if abs(self.bias_value) > 0:
                    self.update_labels(status="subtract bias...")
                    imobj.subtract_bias(self.bias_value)
                    if self.check_stop(): return

                # gradient
                if self.opt_gradient.get():
                    self.update_labels(status="calc gradient...")
                    imobj.gradcorr(resolution_factor)
                    if self.check_stop(): return
                    imobj.gradcorr(resolution_factor)
                    if self.check_stop(): return

                    # write gradient-corrected image
                    if self.set_export_corr_input.get():
                        self.update_labels(status="write gradcorr...")
                        imobj.write_image("_1_gradcorr")
                        if self.check_stop(): return

                # histogram
                if self.opt_histogram.get():
                    self.update_labels(status="calc histogram...")
                    imobj.calc_histogram(circular=self.set_circular_hist.get())
                    if self.check_stop(): return

                # radial profile
                if self.opt_radprof.get():
                    self.update_labels(status="calc radial profile...")
                    imobj.calc_rad_profile(statistics=self.radio_statistics.get(), extrapolate_max=self.set_extrapolate_max.get(), resolution_factor=resolution_factor)
                    if self.check_stop(): return

                # synthetic flat
                if self.opt_synthflat.get():

                    # calculate synthetic flat
                    self.update_labels(status="calc synthetic flat...")
                    imobj.calc_synthetic_flat(grey_flat=self.set_grey_flat.get())
                    if self.check_stop(): return

                    # write synthetic flat tif
                    self.update_labels(status="write flat...")
                    imobj.write_image("_2_synthflat", flat=True)
                    if self.check_stop(): return

                    # correct input
                    if self.set_export_corr_input.get():
                        self.update_labels(status="apply synthetic flat...")
                        imobj.apply_synthetic_flat()
                        if self.check_stop(): return

                        self.update_labels(status="write flatcorr...")
                        imobj.write_image("_3_flatcorr")
                        if self.check_stop(): return

                self.update_labels(status="finished.")
                print("Finished file.")

        except Exception as e:
            print("\nERROR!!")
            print("during status:", self.label_status_var.get())
            print("message:", e)
            self.update_labels(status="unknown error...")
            raise e
            return
        finally:
            self.running = False
            self.update_labels(file=str(len(self.loaded_files)) + " files")
        return


if __name__ == '__main__':
    new = NewGUI()

