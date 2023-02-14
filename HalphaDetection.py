#! /Users/sehuang/opt/anaconda3/bin/python3.8

from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from webdriver_manager.chrome import ChromeDriverManager
import requests
import argparse as arg
import matplotlib.pyplot as plt
from astropy.io import fits
from tqdm import tqdm
import sys
import pandas as pd
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from astroquery.simbad import Simbad

desc = 'find the H-alpha counter part'
parser = arg.ArgumentParser(description = desc)

parser.add_argument('-o', '--object', default = ' ', help = 'Input object name, default: " "')
parser.add_argument('-c', '--coord', default = ' ', help = 'Input decimal object coordinates in RA and Dec, e.g., "241.1743,-39.2209". Default: "0.0,0.0"')

args = parser.parse_args()

object = args.object
coord = args.coord

print('Verifing the required info...')
if object == ' ':
	print('    Please input object name')
	sys.exit()
else:
	print(f'    Target: {object}')

if (coord == ' ') and (object != ' '):
	coord = Simbad.query_object(object)
	starcoord = SkyCoord(coord['RA'], coord['DEC'], unit = (u.hourangle, u.deg))
	ra0, dec0 = starcoord.ra.degree[0], starcoord.dec.degree[0]
	print(f'    Coordinates: {ra0:.5f}, {dec0:.5f}')
elif coord == ' ':
	print('    Please input coordinates of the object')
	sys.exit()
else:
	coord = coord.split(',')
	if len(coord) != 2:
		print(f'    Your input coordinates are not correct: {coord}')
		print('    Please refer to "241.1743,-39.2209"')
		sys.exit()
	print(f'    Coordinates: {coord[0]:.5f}, {coord[1]:.5f}')
	starcoord = SkyCoord(float(coord[0]), float(coord[1]), frame = 'icrs', unit = 'deg')
	ra0, dec0 = starcoord.ra.degree, starcoord.dec.degree

print('')
print('Downloading the data...')

options = Options()
options.add_argument("--disable-notifications")
options.add_argument("--headless")

print('    Connecting to website...')
driver = webdriver.Chrome(ChromeDriverManager().install(), chrome_options = options)
driver.get("http://www-wfau.roe.ac.uk/sss/halpha/hapixel.html")

coord = driver.find_element_by_name('coords')
sizex = driver.find_element_by_name('size')
sizey = driver.find_element_by_name('sizey')
coord.send_keys(f'{starcoord.to_string("hmsdms", sep = " ")}')
sizex.send_keys('15')
sizey.send_keys('15')

coord.submit()
print('    Saving the data...')

#next page
download = driver.find_element_by_xpath('//table[1]/tbody/tr[2]/td[1]/a')
link = download.get_attribute('href')
response = requests.get(link, timeout = 1)

if response.status_code != 200:
	response.raise_for_status()
	driver.close()
	sys.exit()

with open(f'{object}.fits', 'wb') as file:
	file.write(response.content)
	file.close()
driver.close()
print(f'        Save the file as {object}.fits')

print('Processing the data...')
fits_file = fits.open(f'{object}.fits')
fits_header = fits_file[0].header
fits_image = fits_file[0].data
fits_data = fits_file[1].data
fits_file.close()

dataframe = pd.DataFrame(columns = range(len(fits_data[0])), index = range(3))

for i in tqdm(range(len(fits_data))):
	dataframe.loc[i] = fits_data[i]

pic_xlen = fits_header['NAXIS1']
pic_xcen = int(pic_xlen/2)
digitperdegree = pic_xlen/(15/60)
r = 1/60    

fig = plt.figure(dpi = 576)
plt.imshow(fits_image, cmap = 'binary', origin = 'lower', vmax = fits_image.mean()+7*fits_image.std())
plt.axis('square')
fig.set_facecolor('w')
plt.xlim(pic_xcen-r*digitperdegree, pic_xcen+r*digitperdegree)
plt.ylim(pic_xcen-r*digitperdegree, pic_xcen+r*digitperdegree)
plt.colorbar()
plt.xlabel('pixel')
plt.ylabel('pixel')
plt.savefig(f'{object}_Im.png')
print(f'    Image ----> {object}_Im.png')

fig, ax = plt.subplots(dpi = 576)
plt.scatter(dataframe.iloc[:,0], dataframe.iloc[:,1], marker = 'o', color = 'k', s = 10)
plt.scatter(ra0-r/20, dec0, marker = 1, color = 'r', s = 30)
plt.scatter(ra0, dec0+r/20, marker = 2, color = 'r', s = 30)
ax.set_xlabel('RA offset [arcsec]')
ax.set_ylabel('DEC offset [arcsec]')
plt.axis('square')
plt.xticks(np.arange(ra0-r, ra0+2*r, 0.25/60))
plt.yticks(np.arange(dec0-r, dec0+2*r, 0.25/60))
def xtickl(x, pos):
	#"""The two arguments are the value and tick position."""
	return f'{(x-ra0)*3600:.0f}'
def ytickl(y, pos):
	return f'{(y-dec0)*3600:.0f}'
ax.xaxis.set_major_formatter(xtickl)
ax.yaxis.set_major_formatter(ytickl)
plt.xlim(ra0+r, ra0-r)
plt.ylim(dec0-r, dec0+r)
plt.title(f'{object} ({ra0:.5f}, {dec0:.5f})')
fig.set_facecolor('w')
plt.savefig(f'{object}_FC.png')
print(f'    Finding Chart ----> {object}_FC.png')

esstar_n = ((ra0-dataframe[0])**2+(dec0-dataframe[1])**2)<(2/3600)**2
esstar = dataframe[esstar_n]
if len(esstar) == 0:
	print('    No counter part matched')
	sys.exit()
elif len(esstar) > 1:
	print('    More than one counter part matched')
	sys.exit()

#find the mean and std value within 0.5 mag by the target
esstar_sr = float(esstar[5])
esstar_ha = float(esstar[4])
es_srmag_n = np.abs(dataframe[5]-esstar_sr)<=0.25
es_srmag = dataframe[es_srmag_n]

x, y = es_srmag[5], es_srmag[5]-es_srmag[4]

def sigclip(x, y, sig, n):
	for i in range(n):
		y_mean, y_std = y.mean(), y.std()
		n_sig_y = np.abs(y-y_mean)<sig*y_std
		y = y[n_sig_y]
		x = x[n_sig_y]
		if len(y[n_sig_y]==0) == 0:
 			break
	return x, y
x, y = sigclip(x, y, 3, 20)
y_mean, y_std = y.mean(), y.std()
sig_esstar = ((esstar_sr-y_mean)-esstar_ha)/y_std

fig = plt.figure(dpi = 576)
plt.scatter(dataframe[5], dataframe[4], marker = 'o', color = 'k', s = 5)
plt.plot([0,21],[0,21], linestyle = '--', color = 'k', linewidth = 0.3)
plt.scatter(x, x-y, marker = 'o', color = 'g', s = 5)
plt.scatter(esstar_sr, esstar_sr-y_mean, marker = '_', color = 'b', s = 50)
plt.scatter(esstar_sr, esstar_sr-y_mean-3*y_std, marker = '_', color = 'b', s = 50)
plt.scatter(esstar_sr, esstar_ha, marker = '+', color = 'r', s = 100)
plt.axis('square')
plt.xticks(np.arange(0, 21, 2))
plt.yticks(np.arange(0, 21, 2))
if (esstar_sr > 17.5):
	plt.xlim(20.5, 14.5)
	plt.ylim(20.5, 14.5)
elif (esstar_sr < 8.5):
	plt.xlim(11.5, 5.5)
	plt.ylim(11.5, 5.5)
else:
	plt.xlim(esstar_sr+3, esstar_sr-3)
	plt.ylim(esstar_sr+3, esstar_sr-3)
plt.xlabel('SR [mag]')
plt.ylabel(r'H$_\alpha$ [mag]')
plt.title(f'{object} H$_\\alpha$: {sig_esstar:.1f} $\sigma$')
fig.set_facecolor('w')
plt.savefig(f'{object}_ha2sr.png')
print(f'    H-alpha to Short Red ----> {object}_ha2sr.png')

