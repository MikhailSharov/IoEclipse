##########################################
###
###  Echelle Reduction script for use with APO 3.5-m telescope.
###
###   Directory structure required:
###      WD/:
###        echelleReduction.py
###        filter_redflat.fits
###        badcols
###        dcr
###        dcr.par
###        database/:
###          apechtrace130522
###          ecarcnewref.ec
###        raw/:
###          DATE/
###            unprocessed data
###        reduced/: 
###          (stuff will go here)
###
#########
###
###  Edited by Jean McKeever 11/22/13
###
#########

### Import PyRAF packages
import sys,subprocess,os,datetime as dt
from subprocess import call
startdir = os.getcwd()

from pyraf import iraf
from pyraf.iraf import noao,imred,ccdred,proto,twodspec,apextract,echelle,images,imgeom,generic,astutil,onedspec,crutil
os.chdir(startdir)

### A few simple definitions
yes='yes'
no='no'
INDEF='INDEF'
STDOUT=sys.stdout
STDERR=sys.stderr

### Set up necessary files/directories
date=str(sys.argv[1])
call(['mkdir','reduced/'+date])
logfile=open('reduced/'+date+'/reductionlog.txt','w')
sys.stdout=logfile
sys.stderr=logfile

print 'Copying data files and other necessary files'
    #call(['mkdir','reduced/'+date+'/database'])
params = 'cp -r database/ reduced/'+date+'/database/'
call(params,shell=True)
    #sample='3350:3561,3574:3578,3590:3731,3740:3746,3752:3756,3762:3818,3822:3837,3843:3856,3862:3877,3882:3926,3940:3961,3975:4000,4007:4027,4035:4140,4150:4170,4174:4219,4230:4265,4273:4279,4282.5:4297.5,4305:4321,4325:4379,4382.8:4425,4435:4477.5,4482:4492,4503:4645,4650:4685,4690:5138,5140:5165,5185:5323,5330:5368,5372:6555,6570:7593,7615:10900'


#call(['cp','badcols','reduced/'+date])
call(['cp','dcr.par','reduced/'+date])
call(['cp','echmask.pl','reduced/'+date])
call(['cp','filter_redflat.fit','reduced/'+date])
params='cp raw/'+date+'/*.fits reduced/'+date+'/'
call(params,shell=True)
print 'Changing to reduced data directory'
os.chdir('reduced/'+date)

# SET DATA TYPE
# Set the type of data being extracted to 'Echelle'. Forced review of the 
# parameters being set is suppressed.
iraf.setinstrument('echelle',review=no)

# MAKE IMAGE LISTS
# note there are blue flats and 'no filter' (red) flats for ARCES.
# this is because the instrument is undersensitive in the blue.
# lists include: allflats, flat_blue, flat_red, biases, arcs (Th-Ar), targets
# also makes the list objflat (contains both targets and flats)
# (darks are also an option, but are omitted by default)
print 'Generating file lists...'
iraf.hselect('*.fits','$I','IMAGETYP == "flat"',Stdout='allflats')
#iraf.hselect('*.fits','$I','IMAGETYP == "Flat"',Stdout='allflats')
iraf.hselect('@allflats','$I','FILTER == "Blue"',Stdout='flat_blue')
iraf.hselect('@allflats','$I','FILTER == "Open"',Stdout='flat_red')
iraf.hselect('*.fits','$I','IMAGETYP == "zero"',Stdout='biases')
#iraf.hselect('*.fits','$I','IMAGETYP == "bias"',Stdout='biases')
#iraf.hselect('*.fits','$I','IMAGETYP == "Bias"',Stdout='biases')
#iraf.hselect('*.fits','$I','IMAGETYP == "dark"',Stdout='darks')
iraf.hselect('*.fits','$I','IMAGETYP == "comp"',Stdout='arcs')
#iraf.hselect('*.fits','$I','IMAGETYP == "Comp"',Stdout='arcs')
iraf.hselect('*.fits','$I','IMAGETYP == "object"',Stdout='targets')
#iraf.hselect('*.fits','$I','IMAGETYP == "Object"',Stdout='targets')
outfile=open('objflat','w')
iraf.hselect('*.fits','$I','IMAGETYP == "flat"',Stdout=outfile)
#iraf.hselect('*.fits','$I','IMAGETYP == "Flat"',Stdout=outfile)
iraf.hselect('*.fits','$I','IMAGETYP == "object"',Stdout=outfile)
#raf.hselect('*.fits','$I','IMAGETYP == "Object"',Stdout=outfile)
outfile.close()

# SET DISPAXIS (dispersion axis)
# you can set Dispaxis = 1 (horizontal) or 2 (vertical)
print 'Setting dispaxis...'
iraf.hedit('*.fits','dispaxis','1',add=yes,verify=no,show=yes,update=yes)

# CREATE ECHELLE BAD PIXEL MASK
# Use the task TEXT2MASK to create the mask from the ASCII 'badcols' file
#print 'Creating bad pixel mask from file "badcols"...'
#iraf.text2mask('badcols','echmask',ncols=2128,nlines=2068,linterp=1,cinterp=1,
#               square=1,pixel=1)
#  Note: This ^ doesn't work currently on 64-bit systems. Does not produce 
#  a working mask. Use 32-bit IRAF to produce mask.

#  Forces interpolation of bad pixel masks along the
iraf.fixpix('@objflat','echmask.pl',linterp=1,verbose=yes)
iraf.fixpix('@arcs','echmask.pl',linterp=1,verbose=yes)


# MAKE AN AVERAGE (FIDUCIAL) BIAS AND DARK (if applicable)
# by default, we skip the darks, but this can be changed if desired
iraf.imdelete('bias_fid',verify=no)
iraf.zerocombine('@biases',output='bias_fid',combine='median',
                 reject='none',ccdtype='',process=no,delete=no,
                 scale='none')
#iraf.imdelete('dark_fid',verify=no)
#iraf.darkcombine('@darks',output='dark_fid',combine='average',
#                 reject='avsigclip',ccdtype='',process=no,delete=no,
#                 scale='none')


# CALIBRATE ARCS, OBJECT AND FLAT FIELD SPECTRA
# Calibrate object images and flat field images with the task ccdproc.  The
# calibrations applied at this step are the bias level subtraction
# and trimming.  The trim section is hardcoded here.
trimsec='[200:1850,1:2048]'
print 'Calibrating object images and flats (bias, and trimming)...'
iraf.ccdproc('@objflat',output='',ccdtype='object',trimsec=trimsec,
             zerocor=yes,darkcor=no,flatcor=no,trim=yes,
             fixpix=no,overscan=no,zero='bias_fid',interactive=no)
print 'Calibrating object images and flats (bias, and trimming)...'
iraf.ccdproc('@flat_red',output='',ccdtype='flat',trimsec=trimsec,
             zerocor=yes,darkcor=no,flatcor=no,trim=yes,
             fixpix=no,overscan=no,zero='bias_fid',interactive=no)
print 'Calibrating object images and flats (bias, and trimming)...'
iraf.ccdproc('@flat_blue',output='',ccdtype='flat',trimsec=trimsec,
             zerocor=yes,darkcor=no,flatcor=no,trim=yes,
             fixpix=no,overscan=no,zero='bias_fid',interactive=no)
print 'Calibrating ThAr arcs...'
iraf.ccdproc('@arcs',output='',ccdtype='comp',trimsec=trimsec,
             zerocor=yes,darkcor=no,flatcor=no,trim=yes,
             fixpix=no,overscan=no,zero='bias_fid',interactive=no)


### Cosmic ray removal using DCR. Histogram method. astro-ph/0311290
flist=open('objflat')
for ifile in flist:
    ifile=ifile.strip()
    params='../../dcr '+ifile+' '+ifile+' crmask >> reductionlog.txt'
    call(params,shell=True)

flist.close()


# COMBINE 'RED' AND 'BLUE' FLATS SEPARATELY
# Separately combine 'red' and 'blue' flat field frames to produce mean flats
print 'Combining red and blue flats into separate fiducial flats...'
iraf.flatcombine('@flat_blue',output='flat_blue_mean.fits',combine='median',
                 ccdtype='',reject='avsigclip',process=no,delete=no,
                 scale='none',mclip=yes,lsigma=3.0,hsigma=3.0)
iraf.flatcombine('@flat_red',output='flat_red_mean.fits',combine='median',
                 ccdtype='',reject='avsigclip',process=no,delete=no,
                 scale='none',mclip=yes,lsigma=3.0,hsigma=3.0)

# Combine the 'red' and 'blue' mean flats to result in a 'superflat' or 
# fiducial flat by which object images will eventually be divided
print 'Combining average red and blue flats into a superflat...'
iraf.imarith('flat_red_mean.fits','*','filter_redflat.fit',
             'flat_red_mean_filt.fits',verbose=yes)
iraf.imarith('flat_blue_mean.fits','+','flat_red_mean_filt.fits','junk',
             verbose=yes)
iraf.imarith('junk','/',2.0,'flat_fid',verbose=yes)
iraf.imdelete('junk',verify=no)


# EXTRACT AND NORMALIZE SUPERFLAT
# Magnify superflat by 4 in the cross-dispersion direction.
# Then apply model apertures, recenter, resize, trace, and extract.
# Scattered light subtraction isn't necessary because superflat will be
# normalized anyway (we only want pixel-to-pixel variation).
 
print 'Resample the superflat in the y direction...'
iraf.magnify(input='flat_fid.fits',output='flat_fid_mag.fits',xmag=1,ymag=4)
iraf.hedit('flat_fid_mag.fits',fields='CCDSEC',value='[200:1850,1:8189]',
           add=no,addonly=no,delete=no,verify=no,show=yes,update=yes)

print 'Modeling and extracting the superflat...'
iraf.apall(input='flat_fid_mag.fits',ref='echtrace130522',format='echelle',
           interactive=no,find=no,recenter=yes,resize=yes,edit=no,trace=yes,
           fittrace=no,extract=yes,extras=no,review=no,line=825,nsum=10,
           lower=-14.0,upper=14.0,b_function='chebyshev',b_order=2,
           b_niterate=3,b_naverage=-3,b_sample='-22:-15,15:22',width=18.0,
           radius=18.0,npeaks=INDEF,shift=no,ylevel=.05,t_nsum=5,t_step=1,
           t_function='legendre',t_order=5,t_naverage=3,t_niterate=3,
           t_low_reject=2.5,t_high_reject=2.5,t_nlost=3,t_sample='*',
           background='fit',weights='variance',clean=no,saturation=INDEF,
           gain='GAIN',readnoise='RDNOISE',pfit='fit2d',apertures='6-96')


print 'Normalizing the superflat...'
iraf.imcopy('flat_fid_mag.ec.fits','nonorm_flat_fid_mag.ec.fits')
iraf.imdelete('flat_fid_mag.ec.fits',verify=no)

# The extraction of the flat ends up with some negative values, which 
# skews the fit in the following steps, so I added a small bias level
# equal to the minimum value in the

minval=iraf.imstat('nonorm_flat_fid_mag.ec.fits',fields='min',Stdout=1)[-1]
minval=float(minval)
if minval < 0:
    iraf.imcopy('nonorm_flat_fid_mag.ec.fits','nonorm_flat_fid_mag_bias.ec.fits')
    iraf.imdelete('nonorm_flat_fid_mag.ec.fits',verify=no)
    print 'adding an offset to avoid negative'
    iraf.sarith('nonorm_flat_fid_mag_bias.ec.fits','-',minval-.02,
                'nonorm_flat_fid_mag.ec.fits',verbose=yes)

iraf.sfit(input='nonorm_flat_fid_mag.ec.fits',output='flat_fid_mag.ec.fits',
          type='ratio',replace=no,wavesca=no,logscal=no,override=yes,
          listonly=no,interac=no,sample='*',naverage=-5,function='spline3',
          order=12,low_rej=2,high_rej=3.5,niterate=10,grow=0)



# APPLY MODEL APERTURES TO WAVECALS
# Apply model apertures. No resizing, recentering or tracing.
# Also extract spectra here if doing 1D flatfielding, else extract
# after apflatten
call(['rm','all.list'])
call(['cp','arcs','all.list'])
flist=open('all.list')
for ifile in flist:
    ifile=ifile.strip()
    print 'Resampling the arc by a factor of 4 in the y direction...'
    iraf.magnify(input=ifile,output=ifile,xmag=1,ymag=4)
    iraf.hedit(images=ifile,fields='CCDSEC',value='[200:1850,1:8189]',add=no,
               addonly=no,delete=no,verify=no,show=yes,update=yes)
    print 'Applying model apertures to the arc and extracting spectra...'
    iraf.apall(input=ifile,reference='flat_fid_mag',format='echelle',
               interactive=no,find=no,recenter=no,resize=no,edit=no,
               trace=no,fittrace=no,extract=yes,extras=no,review=no,
               shift=no,background='none',weights='none',apertures='6-96')

flist.close()
iraf.hselect('*.ec.fits','$I','IMAGETYP == "comp"',Stdout='arcs_extracted')

# APPLY MODEL APERTURES TO OBJECTS
#
# To properly account for and deal with the spatial aliasing of ARCES spectra,
# the spectra must be resampled in the order direction.  (Uses the
# images and imgeom packages to do this with the task magnify)
#
# Magnify objects by 4 in the cross-dispersion direction.
# Then apply model apertures, recenter, resize, and trace.
# Manual interactive scattered light subtraction because interactive apscat
#   in automated script doesn't work. If doing automatic apscat, turn off
#   interaction and check apscat1 and apscat2 for params.
# Then recenter and resize apertures again for better fit.
# Then extract the spectra (for 1D flatfielding). 

#Values from old script. Probably should refit interactively to confirm
iraf.apscat1.function='spline3'
iraf.apscat1.order=25
iraf.apscat1.low_reject=10
iraf.apscat1.high_reject=0.5
iraf.apscat1.niterate=5
iraf.apscat2.function='spline3'
iraf.apscat2.order=6
iraf.apscat2.low_reject=3
iraf.apscat2.high_reject=0.5
iraf.apscat2.niterate=5


call(['rm','all.list'])
call(['cp','targets', 'all.list'])
flist=open('all.list')
for ifile in flist:
    print 'Resampling the object by a factor of 4 in the y direction...'
    s1=ifile.strip()[:-5]
    iraf.magnify(input=s1,output=s1,xmag=1,ymag=4)
    iraf.hedit(images=s1,fields='CCDSEC',value='[200:1850,1:8189]',add=no,
               addonly=no,delete=no,verify=no,show=yes,update=yes)
    print 'Applying model apertures to the object...'
    iraf.apall(input=s1,ref='flat_fid_mag',format='echelle',interactive=no,
               find=no,recenter=yes,resize=yes,edit=no,trace=no,fittrace=no,
               extract=no,extras=yes,shift=no,width=18.0,radius=18.0,
               ylevel=0.05,background='fit',b_function='chebyshev',b_order=2,
               b_sample='-22:-15,15:22',b_naverage=-3,weights='variance',
               pfit='fit2d',readnoise='RDNOISE',gain='GAIN',apertures='6-96')
    print "Removing scattered inter-order light from the spectrum..."
    s2='noscat'+s1
    iraf.imcopy(s1,s2)
    iraf.imdel(s1,verify=no)
    print "Removing scattered light with apscatter..."
    iraf.apscatter(input=s2,output=s1,ref=s1,interac=no,find=no,recent=no,
                   resize=no,edit=no,trace=no,fittrace=no,subtrac=yes,
                   smooth=yes,fitscat=no,fitsmoo=no,nsum=-10)
    s3=s1+'.ec.fits'
    print "Fine-tuning apertures and extracting spectra..."
    iraf.apall(input=s1,output=s3,ref=s1,format="echelle",interactive=no,
               find=no,recenter=no,resize=no,edit=no,trace=no,fittrace=no,
               extract=yes,extras=yes,review=no,shift=no,background='fit',
               b_function='chebyshev',b_order=1,b_sample='-22:-15,15:22',
               b_naverage=-3,weights='variance',clean=yes,pfit='fit2d',
               readnoise='RDNOISE',gain='GAIN',apertures='6-96')
             

    
flist.close()
iraf.hselect('*.ec.fits','$I','IMAGETYP == "object"',Stdout='targets_extracted')


# ADD 'REFSPEC1' PARAMETER TO OBJECT SPECTRA HEADERS
# Use UNIX awk command to add image names to comparison spectra under the
# "REFSPEC1" header parameter
flist=open('arcs_extracted')
for ifile in flist:
    ifile=ifile.strip()
    iraf.hedit(ifile,'REFSPEC1',ifile,add=yes,verify=no)

    
flist.close()


# DIVIDE OBJECT SPECTRA AND WAVECALS BY THE SUPERFLAT (1D FLATFIELDING)
#
# Divide by the normalized, extracted fiducial flat prior to dispersion 
# correcting object images.  The reasoning is as follows:
#
#	"Flatfielding should be done prior to dispersion correction, i.e.
#	strictly in pixel space. Dispcor resamples the spectra when it 
#	linearizes the dispersion and thus loses some FPN information." 
#	(J. Thorburn)
#
# Added the provision for flat-fielding the wavecals in hopes of recovering 
# response particularly in the blue, leading to the ID'ing of more lines.
#
# PL: 2D flatfielding is bad in our case due to narrow apertures.
#     So we use this 1D flatfielding, assuming response is same in the
#     cross-dispersion direction.

print "Flatfielding extracted arcs and object spectra..."
flist=open('arcs_extracted')
for ifile in flist:
    ifile=ifile.strip()
    iraf.sarith(ifile,'/','flat_fid_mag.ec.fits',ifile,clobber=yes)

flist.close()

flist=open('targets_extracted')
for ifile in flist:
    ifile=ifile.strip()
    iraf.sarith(ifile+'[*,*,4]','/','flat_fid_mag.ec.fits','sig_'+ifile,clobber=yes)
    iraf.sarith(ifile+'[*,*,1]','/','flat_fid_mag.ec.fits',ifile,clobber=yes)

flist.close()

iraf.hedit('sig*.fits','bandid1','',delete=yes,verify=no,show=yes,update=yes)


##
## Set UTMIDDLE Header keyword. = Dateobs +exptime/2, and jd keywords
## Set airmass sets utmiddle too!!!
call(['rm','all.list'])
params='ls sig_* > all.list; cat arcs_extracted targets_extracted >> all.list'
call(params,shell=True)


flist=open('all.list')
for ifile in flist:
    ifile=ifile.strip()
    iraf.hedit(ifile,'UTMIDDLE',value='0',addonly=yes,delete=no,verify=no)
    datobs=iraf.hsel(ifile,'date-obs','yes',Stdout=1)[0]
    exptime=iraf.hsel(ifile,'exptime','yes',Stdout=1)[0]
    dateobs=dt.datetime.strptime(datobs,'%Y-%m-%dT%H:%M:%S.%f')
    utmiddle=dateobs+dt.timedelta(0,float(exptime)/2)
    iraf.hedit(ifile,'UTMIDDLE',value=utmiddle.strftime('%H:%M:%S.%f'),
               add=yes,delete=no,addonly=no,verify=no)


flist.close()


print 'Setting Julian Date Keywords'
iraf.setjd('@targets_extracted',date='date-obs',exposure='exptime',ra='ra',
           dec='dec',epoch='equinox',utdate=yes,uttime=yes,listonly=no)
iraf.setairmass('@targets_extracted',date='date-obs',exposure='exptime',ra='ra',
                dec='dec',equinox='equinox',utmiddle='utmiddle',st='lst',
                airmass='airmass',show=yes,update=yes,override=yes,
                intype='beginning',outtype='middle',ut='date-obs')


iraf.setjd('sig_//@targets_extracted',date='date-obs',exposure='exptime',ra='ra',
           dec='dec',epoch='equinox',utdate=yes,uttime=yes,listonly=no)
iraf.setairmass('sig_//@targets_extracted',date='date-obs',exposure='exptime',ra='ra',
                dec='dec',equinox='equinox',utmiddle='utmiddle',st='lst',
                airmass='airmass',show=yes,update=yes,override=yes,
                intype='beginning',outtype='middle',ut='date-obs')


# DISPERSION CALIBRATE WAVECALS
#print "Median combining arcs and dispersion calibrating the fiducial arc..."
#call(['rm','arc_fid.ec.fits'])
#iraf.scombine('@arcs_extracted',output='arc_fid.ec.fits',group="apertures",
#              combine="median",first=no,w1=INDEF,w2=INDEF,dw=INDEF,nw=INDEF,
#              scale="none",log=no)
#iraf.ecreidentify(images='arc_fid.ec.fits',reference='arcnewref.ec',cradius=2,
#                  threshold=50,shift=INDEF)
iraf.ecreidentify(images='@arcs_extracted',reference='arcnewref.ec',
                  shift=INDEF,cradius=2,threshold=45,refit=yes)


# ASSIGN WAVECALS TO OBJECT SPECTRA
# Assign comparison spectra to object images based on times the images were 
# made, with weights for later dispersion solution assignment
iraf.refspectra('@targets_extracted',references='@arcs_extracted',
                select='interp',sort='utmiddle',group='observat',time=yes,
                confirm=no,override=yes,assign=yes)

iraf.refspectra('sig_//@targets_extracted',references='@arcs_extracted',
                select='interp',sort='utmiddle',group='observat',time=yes,
                confirm=no,override=yes,assign=yes)

iraf.refspectra('@arcs_extracted',references='@arcs_extracted',
                select='match',confirm=no,override=yes,assign=yes)


# APPLY DISPERSION SOLUTION TO OBJECT SPECTRA
# Apply the dispersion solution to object spectra based on the previously-
# assigned comparison spectra
iraf.dispcor(input='@targets_extracted',output='@targets_extracted',flux=no,
             linearize=no,log=no,verbose=yes,samedisp=no)

iraf.dispcor(input='sig_//@targets_extracted',output='sig_//@targets_extracted',flux=no,
             linearize=no,log=no,verbose=yes,samedisp=no)

# Just as a check to the dispersion solution.
iraf.dispcor(input='@arcs_extracted',output='z//@arcs_extracted',flux=no,
             linearize=no,log=no,verbose=yes,samedisp=no)



###
###  Continuum fitting and removal different orders for different parts...
###  Trial and error to get it working.
###
print 'Continuum fitting the spectra...'

iraf.continuum(input='@targets_extracted',output='fitcont//@targets_extracted',
               lines='*',type='fit',replace=no,wavescale=yes,logscale=no,
               listonly=no,interactive=no,sample='*',naverage=1,
               function='spline3',order=4,niterate=20,low_reject=1.15,
               high_reject=2.5,markrej=no,grow=0,override=no)


# Good fits for old 90-92 Ca H and K
#if grating == 'old':
#    iraf.continuum(input='@targets_extracted',
#                   output='fitcont//@targets_extracted',lines='90,91,92',
#                   type='fit',replace=no,wavescale=yes,logscale=no,listonly=no,
#                   interactive=no,sample='3850:3920,3945:3960,3980:4050',
#                   naverage=-3,function='spline3',order=5,niterate=10,
#                   low_reject=2.5,high_reject=5.0,markrej=no,grow=0,
#                   override=yes)

# Order 92 of the new grating was the only one that gave issues. CA H or K.
#if grating == 'new':
#    iraf.continuum(input='@targets_extracted',
#                   output='fitcont//@targets_extracted',lines='92',type='fit',
#                   replace=no,wavescale=yes,logscale=no,listonly=no,
#                   interactive=no,sample='3860:3920,3945:3960',naverage=-3,
#                   function='spline3',order=5,niterate=10,low_reject=1.5,
#                   high_reject=5.0,markrej=no,grow=0,override=yes)
print 'Combining spectra...'
iraf.scombine('@targets_extracted','actual_flux//@targets_extracted',group='images',
              combine='average',reject='none',first=no,w1=INDEF,w2=INDEF,dw=INDEF,
              nw=INDEF,log=no,scale='none',weight='none',zero='none')


print 'Combining spectra...'
iraf.scombine('@targets_extracted','sum//@targets_extracted',group='images',
              combine='sum',reject='none',first=no,w1=INDEF,w2=INDEF,dw=INDEF,
              nw=INDEF,log=no,scale='none',weight='none',zero='none')

iraf.scombine('fitcont//@targets_extracted','sumfitcont//@targets_extracted',
              group='images',combine='sum',reject='none',first=no,w1=INDEF,
              w2=INDEF,dw=INDEF,nw=INDEF,log=no,scale='none',weight='none',
              zero='none')

iraf.sarith('sum//@targets_extracted','/','sumfitcont//@targets_extracted','fullspec//@targets_extracted')

iraf.sarith('@targets_extracted','/','fitcont//@targets_extracted','nocontinuum//@targets_extracted')


iraf.sarith('sig_//@targets_extracted','/','fitcont//@targets_extracted','sig_flat//@targets_extracted')


iraf.scombine('sig_flat//@targets_extracted','sig_fullspec//@targets_extracted',
              group='images',combine='average',reject='minmax',first=no,w1=INDEF,
              w2=INDEF,dw=INDEF,nw=INDEF,log=no,scale='none',weight='none',zero='none',
              nlow='0',nhigh='90',nkeep='1')

### Shouldnt be necessary, but incase there are huge spikes cut them off.
iraf.imreplace(images='fullspec*',value=50,lower=50,upper=INDEF)
iraf.imreplace(images='fullspec*',value=-10,lower=INDEF,upper=-10)

### Fix exposure times. They are off by the number of orders in the trace
### from the summing of orders part

norders=91.0   #115.0 updated to reflect ignoring edge orders
    
flist=open('targets_extracted')
for ifile in flist:
    ifile='fullspec'+ifile.strip()
    exptime=iraf.hsel(ifile,'exptime','yes',Stdout=1)[0]
    print ifile,exptime
    iraf.hedit(ifile,'exptime',value=float(exptime)/norders,add=no,
               verify=no,show=yes,update=yes,addonly=no,delete=no)


flist.close()

#Clean up some middle steps....
#params='rm sum* fitcont* noscat* sig_flat*'
#call(params,shell=True)


print 'ALL DONE! Yay! (=^.^=)  '
os.chdir('../../')
logfile.close()
sys.stdout=STDOUT
sys.stderr=STDERR
