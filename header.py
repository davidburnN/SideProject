#! /Users/sehuang/opt/anaconda3/bin/python3.8

import argparse as arg
from astropy.io import fits
import sys

desc = 'get information from header'
parser = arg.ArgumentParser(description = desc)

modekey = ['comp','log']
parser.add_argument('-m', '--mode', choices = modekey, help = 'printing mode, can be complete = \'comp\' or log = \'log\'')
parser.add_argument('-k', '--keywords', default = '', help = 'keyword in header, type in \'1,2,3\'')
parser.add_argument('files', nargs = '+', help = 'fits files')

args = parser.parse_args()

mode = args.mode
keywords = args.keywords.split(',')
files = args.files

Keywords = []
for i in keywords:
	Keywords.append(i.upper())

if mode == 'comp':
	if len(files) > 1:
		print('can only input single file in complete mode')
	else:
		files = files[0]
		if files[-5:] != '.fits':
			print('file must be fits')
			sys.exit()
		else:
			file = fits.open(files)
			header = file[0].header
			print(f'file name: {i.split("/")[-1]}')
			print(repr(header))
elif mode == 'log':
	vlen_max_list = []
	for i in Keywords:
		vlen_max = 0
		for j in files:
			if j[-5:] != '.fits':
				print(f'{j.split("/")[-1]} is not a fits')
				sys.exit()
			file = fits.open(j)
			header = file[0].header
			if i in header:
				vlen = len(str(header[i]))
			else:
				vlen = 8
			if vlen>vlen_max:
				vlen_max = vlen
			file.close()
		vlen_max_list.append(vlen_max)
	title = f'{"file name":<25}'
	for i in range(len(Keywords)):
		if vlen_max_list[i]>len(Keywords[i]):
			splen = vlen_max_list[i]-len(Keywords[i])
			fronts = int(splen/2)
			backs = splen-fronts
			title += f'|{" "*fronts}{Keywords[i]}{" "*backs}'
		else:
			title += f'| {Keywords[i]} '
	print(title)
	for i in files:
		file = fits.open(i)
		header = file[0].header
		file_name = i.split('/')[-1]
		value = ''
		for j in range(len(Keywords)):
			if Keywords[j] in header:
				Keylen = len(Keywords[j])
				v = header[Keywords[j]]
				vlen = len(str(v))
				if vlen_max_list[j]>Keylen:
					splen = vlen_max_list[j]-vlen
					fronts = int(splen/2)
					backs = splen-fronts
					value += f'|{" "*fronts}{v}{" "*backs}'
				else:
					splen = Keylen-vlen+2
					fronts = int(splen/2)
					backs = splen-fronts
					value += f'|{" "*fronts}{v}{" "*backs}'
			else:
				Keylen = len(Keywords[j])
				if vlen_max_list[j]>Keylen:
					splen = vlen_max_list[j]-8
					fronts = int(splen/2)
					backs = splen-fronts
					value += f'|{" "*fronts}__NONE__{" "*backs}'
				else:
					splen = Keylen-8
					fronts = int(splen/2)
					backs = splen-fronts
					value += f'|{" "*fronts}__NONE__{" "*backs}'
		print(f'{file_name:<25}{value}')


