#! /Users/sehuang/opt/anaconda3/bin/python3.8

import numpy as np
from numpy import ma
from astropy.io import fits
from scipy.optimize import curve_fit
import argparse
import sys

desc = 'aperture photometry'
parser = argparse.ArgumentParser(description = desc)

parser.add_argument('-m1', '--mag1', type = float, default = 1.0, help  = 'magnitude of star 1, default = 1.0')
parser.add_argument('-x1', '--x1', type = int, default = 1, help = 'rought x coordinate of star 1, default = 1')
parser.add_argument('-y1', '--y1', type = int, default = 1, help = 'rought y coordinate of star 1, default = 1')
parser.add_argument('-x2', '--x2', type = int, default = 2, help = 'rought x coordinate of star 2, default = 2')
parser.add_argument('-y2', '--y2', type = int, default = 2, help = 'rought y coordinate of star 2, default = 2')
parser.add_argument('file', nargs = 1, help  = 'input fits file')

args = parser.parse_args()

mag1 = args.mag1
x1, y1 = args.x1, args.y1
x2, y2 = args.x2, args.y2
file = args.file[0]

if file[-5:] != '.fits':
	print('error: must input fits file')
	sys.exit()

fits_file = fits.open(file)
data = fits_file[0].data
header = fits_file[0].header
fits_file.close()
data_x = header['NAXIS1']
data_y = header['NAXIS2']

err_msg = 0
if x1>=data_x:
	err_msg += 1
	print(f'x1 out of image range: {data_x}')
if x2>=data_x:
	err_msg += 1
	print(f'x2 out of image range: {data_x}')
if y1>=data_y:
	err_msg += 1
	print(f'y1 out of image range: {data_y}')
if y2>=data_y:
	err_msg += 1
	print(f'y2 out of image range: {data_y}')
if err_msg>0:
	sys.exit() 

def moff(x, a, b, c):
    return a*(1+(x/c)**2)**-b

w1_err, w2_err = 0, 0
for i in range(10):
	w = 60+i*10
	image1_xl, image1_xh = int(x1-w/2), int(x1+w/2)+1
	image1_yl, image1_yh = int(y1-w/2), int(y1+w/2)+1
	image2_xl, image2_xh = int(x2-w/2), int(x2+w/2)+1
	image2_yl, image2_yh = int(y2-w/2), int(y2+w/2)+1
	if image1_xl<1 or image1_xh>data_x or image1_yl<1 or image1_yh>data_y:
		w1_err += 1
	if image2_xl<1 or image2_xh>data_x or image2_yl<1 or image2_yh>data_y:
		w2_err += 1
	if w1_err>0:
		print('star 1 is near the edge')
		break
	if w2_err>0:
		print('star 2 is near the edge')
		break
	
	image1 = data[image1_yl:image1_yh,image1_xl:image1_xh]
	image2 = data[image2_yl:image2_yh,image2_xl:image2_xh]
	image_x, image_y = np.meshgrid(np.arange(0, w+1, 1), np.arange(0, w+1, 1))
	
	mask1 = image1<=2*image1.std()
	mask2 = image2<=2*image2.std()
	maimage1 = ma.array(image1, mask = mask1)
	maimage2 = ma.array(image2, mask = mask2)
	
	com_x1, com_y1 = ma.sum(maimage1*image_x)/ma.sum(maimage1), ma.sum(maimage1*image_y)/ma.sum(maimage1)
	com_x2, com_y2 = ma.sum(maimage2*image_x)/ma.sum(maimage2), ma.sum(maimage2*image_y)/ma.sum(maimage2)

	image_r1 = np.sqrt((image_x-com_x1)**2+(image_y-com_y1)**2)
	image_r2 = np.sqrt((image_x-com_x2)**2+(image_y-com_y2)**2)

	a1, b1 = curve_fit(moff, image_r1.flatten(), image1.flatten(), bounds = [0, [np.inf, np.inf, np.inf]])
	amp1, beta1, R1 = a1
	a2, b2 = curve_fit(moff, image_r2.flatten(), image2.flatten(), bounds = [0, [np.inf, np.inf, np.inf]])
	amp2, beta2, R2 = a2

	hwhm1 = R1*np.sqrt(2**(1/beta1)-1)
	hwhm2 = R2*np.sqrt(2**(1/beta2)-1)
	ap = 3
	anu_r1 = 6
	anu_r2 = 10
	
	if hwhm1*anu_r2>=w/2 or hwhm2*anu_r2>=w/2:
		continue
	else:
		mask_r1 = (image_r1>ap*hwhm1)
		mask_r2 = (image_r2>ap*hwhm2)
		mask_anu1 = (image_r1<anu_r1*hwhm1) | (image_r1>anu_r2*hwhm1)
		mask_anu2 = (image_r2<anu_r1*hwhm2) | (image_r2>anu_r2*hwhm2)

		image_ap1 = ma.array(image1, mask = mask_r1)
		image_ap2 = ma.array(image2, mask = mask_r2)
		image_anu1 = ma.array(image1, mask = mask_anu1)
		image_anu2 = ma.array(image2, mask = mask_anu2)
		
		backg1 = ma.median(image_anu1)
		backg2 = ma.median(image_anu2)
		
		flux1 = ma.sum(image_ap1)
		flux2 = ma.sum(image_ap2)
		
		net_flux1 = flux1-backg1*(mask_r1==0).sum()
		net_flux2 = flux2-backg2*(mask_r2==0).sum()
		
		error1 = np.sqrt(flux1+backg1*(mask_r1==0).sum())
		error2 = np.sqrt(flux2+backg2*(mask_r2==0).sum())
		
		print(f"net flux1: {net_flux1:.0f} ± {error1:.0f}")
		print(f"net flux2: {net_flux2:.0f} ± {error2:.0f}")
		
		mag2 = mag1-2.5*np.log10(net_flux2/net_flux1)
		
		mag_error1 = 2.5*np.log10(1+error1/net_flux1)
		mag_error2 = 2.5*np.log10(1+error2/net_flux2)
		mag_err_tot = np.sqrt(mag_error1**2+mag_error2**2)
		
		print(f'star1: {mag1:.3f}')
		print(f'star2: {mag2:.3f} ± {mag_err_tot:.5f}')
		
		break


