procedure lacos_spec(input,output,outmask)

# cosmic ray rejection script for long slit spectra
# by Pieter van Dokkum, April 2001
#
# info at http://www.astro.yale.edu/dokkum/lacosmic/


char input{"",prompt="input spectrum"}
char output{"",prompt="cosmic ray cleaned output spectrum"}
char outmask{"",prompt="output bad pixel map (.pl)"}
real gain{2.,prompt="gain (electrons/ADU)"}
real readn{6.,prompt="read noise (electrons)"}
int xorder{9,prompt="order of object fit (0=no fit)"}
int yorder{3,prompt="order of sky line fit (0=no fit)"}
real sigclip{4.5,prompt="detection limit for cosmic rays (sigma)"}
real sigfrac{0.5,prompt="fractional detection limit for neighbouring pixels"}
real objlim{1.,prompt="contrast limit between CR and underlying object"}
int niter{4,prompt="maximum number of iterations"}
bool verbose{yes}

begin

 char blk,lapla,deriv2,med5,sub5,sub5abs,med3,med7,imstatout,inputmask
 char noise,sigmap,firstsel,starreject,gfirstsel,finalsel
 char oldoutput,sigima,galaxy,skymod
 char kernel,gkernel
 real midpt,skylev,sig,sigcliplow
 char l1,l2
 int i
 bool stop
 int previous,npix

 if (!deftask("imcalc")) error(2,"Load package stsdas")
 if (gain<=0) error(2,"Gain is required")

 if (verbose) {
  print("")
  print("_______________________________________________________________________________")
  print("")
  print("                 L.A.Cosmic: Laplacian cosmic ray removal")
  print("")
  print("                   P. van Dokkum, 2001, PASP 113, 1420")
  print("")
  print("                  Spectroscopy version 1.0  (April 2001)")
  print("_______________________________________________________________________________")
  print("")
  }

 # make temporary files

 blk = mktemp("lacos")
 lapla = mktemp("lacos")
 deriv2 = mktemp("lacos")
 kernel = mktemp("lacos")
 gkernel=mktemp("lacos")
 med3 = mktemp("lacos")
 med5 = mktemp("lacos")
 med7 = mktemp("lacos")
 sub5 = mktemp("lacos")
 sub5abs = mktemp("lacos")
 imstatout = mktemp("lacos")
 noise = mktemp("lacos")
 sigmap = mktemp("lacos")
 firstsel = mktemp("lacos")
 starreject=mktemp("lacos")
 gfirstsel = mktemp("lacos")
 finalsel = mktemp("lacos")
 inputmask = mktemp("lacos")
 oldoutput = mktemp("lacos")
 sigima = mktemp("lacos")
 galaxy=mktemp("lacos")
 skymod = mktemp("lacos")

 # set some parameters in standard IRAF tasks

 convolve.bilinear=no
 convolve.radsym=no

 # create Laplacian kernel

 print("0 -1 0;",>> kernel)
 print("-1 4 -1;",>> kernel)
 print("0 -1 0;",>> kernel)

 # create growth kernel

 print("1 1 1;",>> gkernel)
 print("1 1 1;",>> gkernel)
 print("1 1 1;",>> gkernel)

 # initialize loop

 i=1
 stop=no
 previous=0

 imcopy(input,oldoutput,verb-)
 imcopy(input,outmask,verb-)
 imreplace(outmask,0,upper=INDEF,lower=INDEF)

  # subtract object spectra if desired

  if (xorder>0) {
   if (verbose) {
    print("Subtracting object spectra")
    print("  fit order = "//xorder)
    print("")
    }
   fit1d(oldoutput,galaxy,"fit",axis=1,order=xorder,func="leg",low=4.,
        high=4.,nav=1,inter-,sample="*",niter=3,grow=0,cursor="")
   imarith(oldoutput,"-",galaxy,oldoutput)
   }

  else {
   imcopy(oldoutput,galaxy,verb-)
   imreplace(galaxy,0,lower=INDEF,upper=INDEF)
   }

  if (yorder>0) {
   if (verbose) {
    print("Subtracting sky lines")
    print("  fit order = "//yorder)
    print("")
    }
   fit1d(oldoutput,skymod,"fit",axis=2,order=yorder,func="leg",low=4.,high=4.,
        inter-,sample="*",nav=1,niter=3,grow=0,cursor="")
   imarith(oldoutput,"-",skymod,oldoutput)
   }

  else {
   imcopy(oldoutput,skymod,verb-)
   imreplace(skymod,0,lower=INDEF,upper=INDEF)
   }

  # add object spectra to sky model

  imarith(skymod,"+",galaxy,skymod)

 # start iterations

 while(!stop) {

  if (verbose) {
   print("")
   if (i<10) print("_______________________________ Iteration "//i//" ___________________________________")
   if (i>9) print("_______________________________ Iteration "//i//"___________________________________")
   print("")
   print("")
   }

  # add median of residuals to sky model

  median(oldoutput,med5,5,5,zlor=INDEF,zhir=INDEF,verb-)
  imarith(skymod,"+",med5,med5)

  # take second-order derivative (Laplacian) of input image
  # kernel is convolved with subsampled image, in order to remove negative
  # pattern around high pixels

  if (verbose) {
   print("Convolving with Laplacian kernel")
   print("")
   }
  blkrep(oldoutput,blk,2,2)
  convolve(blk,lapla,kernel)
  imreplace(lapla,0,upper=0,lower=INDEF)
  blkavg(lapla,deriv2,2,2,option="average")

  if (verbose) {
   print("Creating noise model using:")
   print("  gain = "//gain//" electrons/ADU")
   print("  readnoise = "//readn//" electrons")
   print("")
   }

  # create noise model

  imcalc(med5,noise,"sqrt(im1*"//gain//" + "//readn//"**2)/"//gain,verb-)
  imreplace(med5,0.00001,upper=0,lower=INDEF)

  # divide Laplacian by noise model

  imarith(deriv2,"/",noise,sigmap)

  # Laplacian of blkreplicated image counts edges twice:

  imarith(sigmap,"/",2.,sigmap)

  # removal of large structure (bright, extended objects)

  imdel(med5)
  median(sigmap,med5,5,5,zlo=INDEF,zhi=INDEF,verb-)
  imarith(sigmap,"-",med5,sigmap)

  # find all candidate cosmic rays
  # this selection includes sharp features such as stars and HII regions

  if (verbose) {
   print("Selecting candidate cosmic rays")
   print("  sigma limit = "//sigclip)
   print("")
   }
  imcopy(sigmap,firstsel,verb-)
  imreplace(firstsel,0,upper=sigclip,lower=INDEF)
  imreplace(firstsel,1,lower=0.1,upper=INDEF)

  # compare candidate CRs to median filtered image
  # this step rejects bright, compact sources from the initial CR list

  if (verbose) {
   print("Removing suspected emission lines")
   print("  selecting cosmic rays > "//objlim//" times object flux")
   print("")
   }

  # subtract background and smooth component of objects

  median(oldoutput,med3,3,3,zlo=INDEF,zhi=INDEF,verb-)
  median(med3,med7,7,7,zlo=INDEF,zhi=INDEF,verb-)
  imarith(med3,"-",med7,med3)
  imarith(med3,"/",noise,med3)
  imreplace(med3,0.01,upper=0.01,lower=INDEF)

  # compare CR flux to object flux

  imcalc(firstsel//","//sigmap//","//med3,starreject,"(im1*im2)/im3",verb-)

  # discard if CR flux <= objlim * object flux

  imreplace(starreject,0,upper=objlim,lower=INDEF)
  imreplace(starreject,1,lower=0.5,upper=INDEF)
  imarith(firstsel,"*",starreject,firstsel)

  # grow CRs by one pixel and check in original sigma map

  convolve(firstsel,gfirstsel,gkernel)
  imreplace(gfirstsel,1,lower=0.5,upper=INDEF)
  imarith(gfirstsel,"*",sigmap,gfirstsel)
  imreplace(gfirstsel,0,upper=sigclip,lower=INDEF)
  imreplace(gfirstsel,1,lower=0.1,upper=INDEF)

  # grow CRs by one pixel and lower detection limit

  sigcliplow = sigfrac * sigclip

  if (verbose) {
   print("Finding neighbouring pixels affected by cosmic rays")
   print("  sigma limit = "//sigcliplow)
   print("")
   }

  convolve(gfirstsel,finalsel,gkernel)
  imreplace(finalsel,1,lower=0.5,upper=INDEF)
  imarith(finalsel,"*",sigmap,finalsel)
  imreplace(finalsel,0,upper=sigcliplow,lower=INDEF)
  imreplace(finalsel,1,lower=0.1,upper=INDEF)

  # determine number of CRs found in this iteration

  imdel(gfirstsel)
  imcalc(finalsel//","//outmask,gfirstsel,"(1-im2)*im1",verb-)
  imstat(gfirstsel,fields="npix",lower=0.5,upper=INDEF,for-) | scan(npix)

  # create cleaned output image; use 3x3 median with CRs excluded

  if (verbose) {
   print("Creating output:")
   print("  bad pixel mask: "//outmask)
   print("  cleaned image: "//output)
   print("")
   }

  imdel(med5)
  imarith(outmask,"+",finalsel,outmask)
  imreplace(outmask,1,lower=1,upper=INDEF)
  imcalc(outmask//","//oldoutput,inputmask,"(1-10000*im1)*im2",verb-)
  median(inputmask,med5,5,5,zloreject=-9999,zhi=INDEF,verb-)
  imarith(outmask,"*",med5,med5)
  if (i>1) imdel(output)
  imcalc(outmask//","//oldoutput//","//med5,output,"(1-im1)*im2+im3",verb-)

  # add sky and object spectra back in

  imdel(oldoutput)
  imcopy(output,oldoutput,verb-)

  imarith(output,"+",skymod,output)

  # cleanup and get ready for next iteration

  if (verbose) {
   print("Cleaning up")
   print("")
   }

   if (verbose) {
   print("")
   print(npix//" cosmic rays found in iteration "//i)
   print("")
   }

  if (npix==0) stop=yes

  i=i+1
  if (i>niter) stop=yes

  # delete temp files

  imdel(blk//","//lapla//","//deriv2//","//med5)
  imdel(med3//","//med7//","//noise//","//sigmap)
  imdel(firstsel//","//starreject//","//gfirstsel)
  imdel(finalsel//","//inputmask)

  }

 imdel(oldoutput//","//skymod//","//galaxy)
 delete(kernel//","//gkernel)

end



