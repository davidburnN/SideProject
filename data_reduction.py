#! /Users/sehuang/opt/anaconda3/bin/python3.8

import argparse
import numpy as np
from numpy import ma
from astropy.io import fits
from astropy import stats
import pandas as pd
import sys

desc = 'data reduction'
parser = argparse.ArgumentParser(description = desc)

rejection_l = ['none','sigclip']
cenfunc_l = ['mean','median']
parser.add_argument('-r', '--rejection', choices = rejection_l, default = 'sigclip', help = 'rejection algorithm, default: sigclip')
parser.add_argument('-t', '--threshold', type = float, default = 3.0, help = 'sigclip sigma threshold, default: 3.0')
parser.add_argument('-i', '--iteration', type = int, default = 10, help = 'maximum number of iteration, default: 10')
parser.add_argument('-c', '--cenfunc', choices = cenfunc_l, default = 'median', help = 'cenfunc of sigclip, default: median')
parser.add_argument('files', nargs = '+', help = 'input files')

args = parser.parse_args()

reject = args.rejection
thresh = args.threshold
iter = args.iteration
cenf = args.cenfunc
files = args.files

columns = ['fileadd','filename','xpix','ypix','object','imagetyp','filter','exptime','number']
file_table = pd.DataFrame(columns = columns, index = range(3))
file_d_table = pd.DataFrame(columns = columns, index = range(3))
file_df_table = pd.DataFrame(columns = columns, index = range(3))

for i in range(len(files)):
	fileadd = files[i]
	if fileadd[-5:] != '.fits':
		continue
	file = fits.open(fileadd)
	filename = fileadd.split('/')[-1]
	header = file[0].header
	file.close()
	ypix = header['NAXIS1']
	xpix = header['NAXIS2']
	imagetyp = header['IMAGETYP']
	exptime = header['EXPTIME']
	if 'OBJECT' in header:
		object = header['OBJECT']
	else:
		object = '__NONE__'
	if 'filter' in header:
		filter = header['FILTER']
	else:
		filter = '__NONE__'
	file_table.loc[i,'fileadd'] = fileadd
	file_table.loc[i,'filename'] = filename
	file_table.loc[i,'xpix'] = xpix
	file_table.loc[i,'ypix'] = ypix
	file_table.loc[i,'object'] = object
	file_table.loc[i,'imagetyp'] = imagetyp
	file_table.loc[i,'filter'] = filter
	file_table.loc[i,'exptime'] = exptime
	file_table.loc[i,'number'] = 1

#checking files completeness
print('#')
print('checking files self-consistant\n')

#
print('#')
print('check which kind of filter do the files use\n')
light_file_table = file_table[file_table['imagetyp']=='LIGHT']

if len(light_file_table) == 0:
	print('there\'s no light frame')
	sys.exit()
else:
	light_filter = light_file_table.groupby('filter').sum().loc[:,'number']
print(light_filter,'\n')

filter_list_pack = light_filter.index
filter_list = []

for i in range(len(filter_list_pack)):
	filter_list.append(filter_list_pack[i])

print('done\n')

#
print('#')
print('check if each filter has the flat field\n')
flat_file_table = file_table[file_table['imagetyp']=='FLAT']

if len(flat_file_table) == 0:
	print('thers\'s no flat field')
	sys.exit()
else:
	flat_filter = flat_file_table.groupby('filter').sum().loc[:,'number']

flat_filter_list_pack = flat_filter.index
flat_filter_list = []
flat_filter_num = []

for i in range(len(flat_filter_list_pack)):
	flat_filter_list.append(flat_filter_list_pack[i])
	flat_filter_num.append(flat_filter.iloc[i])
for i in filter_list:
	print(f'for {i} band\n')
	k = 0
	for j in range(len(flat_filter_list)):
		if i != flat_filter_list[j]:
			continue
		else:
			k += 1
			num = flat_filter_num[j]
	if k == 0:
		print('    do not have enough flat field\n')
		sys.exit()
	else:
		print(f'    this band has {num} flat field\n')

print('done\n')

#
print('#')
print('check if each exptime has the dark field\n')

light_flat_cri1 = file_table['imagetyp']=='LIGHT'
light_flat_cri2 = file_table['imagetyp']=='FLAT'
light_flat_file_table = file_table[light_flat_cri1 | light_flat_cri2]
dark_file_table = file_table[file_table['imagetyp']=='DARK']

light_flat_exptime = light_flat_file_table.groupby('exptime').sum().loc[:,'number']
dark_exptime = dark_file_table.groupby('exptime').sum().loc[:,'number']

light_flat_exptime_list_pack = light_flat_exptime.index
light_flat_exptime_list = []
dark_exptime_list_pack = dark_exptime.index
dark_exptime_list = []
dark_exptime_num = []

for i in range(len(light_flat_exptime_list_pack)):
	light_flat_exptime_list.append(light_flat_exptime_list_pack[i])
for i in range(len(dark_exptime_list_pack)):
	dark_exptime_list.append(dark_exptime_list_pack[i])
	dark_exptime_num.append(dark_exptime.iloc[i])

for i in light_flat_exptime_list:
	k = 0
	print(f'for {i} exposure time')
	for j in range(len(dark_exptime_list)):
		if i != dark_exptime_list[j]:
			continue
		else:
			k += 1
			num = dark_exptime_num[j]
	if k == 0:
		print('    there is no dark frame')
		sys.exit()
	else:
		print(f'    there are {num} dark frames\n')

print('files are self-consistant\n')

#dark combination
print('#')
print('now processing the dark combination\n')

for i in light_flat_exptime_list:
	print(f'for {i} exposure time')
	l = 0
	dark_file_exptime_table = dark_file_table[dark_file_table['exptime']==i]
	for j in range(len(dark_file_exptime_table)):
		fileadd = dark_file_exptime_table.iloc[j,0]
		file = fits.open(fileadd)
		header = file[0].header
		if reject == 'none':
			data = file[0].data
			mask = np.zeros((dark_file_exptime_table.iloc[j,3], dark_file_exptime_table.iloc[j,2]))
		else:
			data = stats.sigma_clip(file[0].data, sigma = thresh, maxiters = iter, cenfunc = cenf)
			mask = data.mask
		if l == 0:
			datacube = data
			maskcube = mask
		elif l == 1:
			datacube = np.concatenate(([datacube], [data]), axis = 0)
			maskcube = np.concatenate(([maskcube], [mask]), axis = 0)
		else:
			datacube = np.concatenate((datacube, [data]), axis = 0)
			maskcube = np.concatenate((maskcube, [mask]), axis = 0)
		l += 1	
		file.close()
		sys.stdout.write('\r')
		sys.stdout.write(f'    [{"="*int(10*j/len(dark_file_exptime_table)):10}]')
		sys.stdout.flush()
	cube = ma.array(datacube, mask = maskcube)
	if reject == 'none':
		dark_image = ma.average(cube, axis = 0)
	else:
		cube = stats.sigma_clip(cube, sigma = thresh, maxiters = iter, cenfunc = cenf, axis = 0)
		dark_image = ma.average(cube, axis = 0)
	header['comment'] = 'combined dark frame'
	header['comment'] = f'  rejection: {reject}'
	if reject == 'sigclip':
		header['comment'] = f'  sigma:   {thresh}'
		header['comment'] = f'  maxiter: {iter}'
		header['comment'] = f'  cenfunc: {cenf}'
	header['comment'] = f"    mean:   {ma.mean(dark_image):.2f}"
	header['comment'] = f'    median: {ma.median(dark_image):.2f}'
	header['comment'] = f'    std:    {ma.std(dark_image):.2f}'
	header['comment'] = f'    max:    {ma.max(dark_image):.2f}'
	header['comment'] = f'    min:    {ma.min(dark_image):.2f}'
	outputfile = f'dark_{int(i):04d}.fits'
	dark_image = ma.filled(dark_image, fill_value = ma.median(dark_image))
	fits.writeto(outputfile, dark_image, header = header)

	sys.stdout.write('\r')
	sys.stdout.write(f'    [{"="*10:10}]')
	sys.stdout.flush()

	print(f'    output: {outputfile}')
print('')
print('done\n')


#dark subtraction
print('#')
print('now processing dark subtraction\n')

for i in range(len(light_flat_file_table)):
	fileadd = light_flat_file_table.iloc[i,0]
	filename = light_flat_file_table.iloc[i,1]
	exptime = light_flat_file_table.iloc[i,7]
	file = fits.open(fileadd)
	filter = light_flat_file_table.iloc[i,6]
	object = light_flat_file_table.iloc[i,4]
	imagetyp = light_flat_file_table.iloc[i,5]
	data = file[0].data
	header = file[0].header
	file.close()
	darkname = f'dark_{int(exptime):04d}.fits'
	darkfile = fits.open(darkname)
	dark = darkfile[0].data
	data_subtrac = data-dark
	outputfile = filename.split('.')[0]+'_d.fits'
	header['comment'] = 'frame has been dark subtracted'
	header['comment'] = f'    subtracted file: {darkname}'
	fits.writeto(outputfile, data_subtrac, header = header)
	sys.stdout.write('\r')
	sys.stdout.write(f'[{"="*int(10*i/len(light_flat_file_table)):10}]')
	sys.stdout.flush()
	file_d_table.loc[i,'fileadd'] = outputfile
	file_d_table.loc[i,'filename'] = outputfile
	file_d_table.loc[i,'exptime'] = exptime
	file_d_table.loc[i,'object'] = object
	file_d_table.loc[i,'filter'] = filter
	file_d_table.loc[i,'imagetyp'] = imagetyp
	
sys.stdout.write('\r')
sys.stdout.write(f'[{"="*10:10}]')
sys.stdout.flush()
print('\n')
print('done\n')


#flat combination
print('#')
print('now processing flat combination\n')

flat_d_file_table = file_d_table[file_d_table['imagetyp']=='FLAT']

for i in flat_filter_list:
	print(f'for {i} band')
	flat_d_filter_file_table = flat_d_file_table[flat_d_file_table['filter']==i]
	l = 0
	s = 0
	for j in range(len(flat_d_filter_file_table)):
		fileadd = flat_d_filter_file_table.iloc[j,0]
		filename = flat_d_filter_file_table.iloc[j,1]
		file = fits.open(fileadd)
		filedata = file[0].data
		header = file[0].header
		file.close()
		if s == 0:
			scale = np.median(filedata)
		data = filedata*scale/np.median(filedata)
		if l == 0:
			datacube = data
		elif l == 1:
			datacube = ma.concatenate(([datacube], [data]), axis = 0)
		else:
			datacube = ma.concatenate((datacube, [data]), axis = 0)
		l += 1
		s += 1
		sys.stdout.write('\r')
		sys.stdout.write(f'    [{"="*int(10*j/len(dark_file_exptime_table)):10}]')
		sys.stdout.flush()
	header['comment'] = 'combined flat field'
	if reject == 'sigclip':
		cube = stats.sigma_clip(datacube, sigma = thresh, maxiters = iter, cenfunc = cenf, axis = 0)
		header['comment'] = f'    rejection = {reject}'
		header['comment'] = f'    threshold = {thresh}'
		header['comment'] = f'    maxiters = {iter}'
		header['comment'] = f'    cenfunc = {cenf}'
	elif reject == 'none':
		cube = datacube
		header['comment'] = f'    rejection = {reject}'
	flatimage = ma.average(cube, axis = 0)
	flatimage = ma.filled(flatimage, fill_value = ma.median(flatimage))
	flatimage = flatimage/ma.mean(flatimage)
	output = f'flat_{i}.fits'
	fits.writeto(output, flatimage, header = header)
	sys.stdout.write('\r')
	sys.stdout.write(f'    [{"="*10:10}]')
	sys.stdout.flush()
	print(f'    output: {output}\n')
print('')
print('done\n')

#flatfielding
print('#')
print('now processing flatfielding\n')

light_d_file_table = file_d_table[file_d_table['imagetyp']=='LIGHT']

for i in range(len(light_d_file_table)):
	fileadd = light_d_file_table.iloc[i,0]
	filename = light_d_file_table.iloc[i,1]
	filter = light_d_file_table.iloc[i,6]
	object = light_d_file_table.iloc[i,4]
	file = fits.open(fileadd)
	header = file[0].header
	data = file[0].data
	file.close()
	flatname = f'flat_{filter}.fits'
	flatfile = fits.open(flatname)
	flat = flatfile[0].data
	flatfile.close()
	data = data/flat
	header['comment'] = 'file has been flatfielded'
	header['comment'] = f'    flatfield = {flatname}'
	output = f'{filename.split(".")[0]}f_{object}.fits'
	fits.writeto(output, data, header = header)
	sys.stdout.write('\r')
	sys.stdout.write(f'[{"="*int(10*i/len(light_d_file_table)):10}]')
	sys.stdout.flush()
	
sys.stdout.write('\r')
sys.stdout.write(f'[{"="*10:10}]')
sys.stdout.flush()
print('\n')
print('done\n')


