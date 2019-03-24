## This program will run int.py on the stars in the list.

import os, sys
from numpy import loadtxt
from Interpolator import interpolator
from time import time

## Specify the input and output files.
infile = 'data_files/MLtesterr3bclump.txt'
outfile = infile.replace('data_files','new_outfiles').replace('.txt','.out')

## Read in the stars.
stars = loadtxt(infile,usecols=(0,1,2,3,4,5))
ids, mass, logg, feoh, alfe, teff = zip(*stars)

## Now read in the lines for outputting later. Cut out the header.
stars = open(infile).readlines()[1:]

## Start the output file if it hasn't been opened. Otherwise,
## figure out which stars we've interpolated and skip them
if not os.path.exists(outfile): 
	open(outfile,'a').write('## KICID      Mass         Logg       Fe/H          Al/Fe         Teff       C/N           ')
	open(outfile,'a').write('M_err       Logg_err     Fe/H_err     Al/Fe_err      Teff_err      IntTeff      IntML      ')
	open(outfile,'a').write('     IntAge          IntC12          IntC13          IntN14          IntXsurf\n')

DoneStars = loadtxt(outfile,usecols=(0,)).tolist()
try: len(DoneStars)
except TypeError: DoneStars = [DoneStars,]

## For the full Tayar sample.
OnesToSkip = (5371676,6792689,8285712,8459156,8609704,5113061,3222519,3543433,8569991,9726045)

t1 = time()
## Now run the loop.
for i in xrange(len(ids)):

	#if i > 20: break

	I, M, F, A, L, T = ids[i], [mass[i],], [feoh[i],], [alfe[i],], [logg[i],], [teff[i],]

	print '\n#', i, '\t', I, '\n'
	## Check that this star hasn't been interpolated before.
	if I in DoneStars: continue

	## Check if its one we want to skip.
	if I in OnesToSkip:
		print 'Skipping this one!'
		lin = '#'+'%.0d'%I+'\t-99.000000\n'
		open(outfile,'a').write(lin)
		continue

	if True: ## Include Age and C/N
		thetef, theml, theage, thec12, thec13, then14, thehyd = interpolator(M, F, A, L, T)[0]

		lin  = stars[i].replace('\n','')
		lin += '       %.1f'%thetef+'\t'
		lin += '   %.6f'%theml+'\t'
		lin += '   %.4e'%theage+'\t'
		lin += '   %.4e'%thec12+'\t'
		lin += '   %.4e'%thec13+'\t'
		lin += '   %.4e'%then14+'\t'
		lin += '   %.4e'%thehyd+'\n'
		open(outfile,'a').write(lin)

	else:
		thetef, theml = interpolator(M, F, A, L, T)[0]

		lin  = stars[i].replace('\n','')
		lin += '       %.1f'%thetef+'\t'
		lin += '   %.6f'%theml+'\n'
		open(outfile,'a').write(lin)


print '\n\nThis whole run took', time()-t1, 'seconds!'



















