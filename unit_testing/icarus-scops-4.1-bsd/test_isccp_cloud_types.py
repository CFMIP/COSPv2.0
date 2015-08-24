#############################################################################
#
# The purpose of this program is to compare the ISCCP simulator test code 
# output to a reference dataset. Specifically, only the test-based 1/+/- cloud 
# representation is compared.
#
# Coded by Dustin Swales 2014 CIRES/NOAA-ESRL-PSD
#
#############################################################################
import os,sys,numpy

wrkDir, filename = os.path.split(os.path.abspath(__file__))
fileNEW  = wrkDir+'/data/output/'+'stdout'
fileREF  = wrkDir+'/data/output/'+'stdout.expected'
print '###################################################################################'
print 'Comparing: '
print '    stdout'
print 'to the reference file: '
print '    stdout.expected'
print '####################################################################################'

#############################################################################
# Determine which parts of the files to compare
#############################################################################

# Number of lines in reference file
fIDref    = open(fileREF,'r')
flinesRef = fIDref.readlines()
nlinesRef = len(flinesRef)

# Number of lines in new file
fIDnew    = open(fileNEW,'r')
flinesNew = fIDnew.readlines()
nlinesNew = len(flinesNew)

# If files are of different length quit.
if (nlinesRef != nlinesNew):
    print "ERROR:  Reference file and comparison file are of different length"
    print "EXITING!!"
    quit()

# Pull out number of points, levels and columns
fIDref   = open(fileREF,'r')
linesREF = fIDref.readlines()
npoints = int(linesREF[20])
nlevs   = int(linesREF[22])
ncols   = int(linesREF[24])
nconfig = 18 # Hardcoded in for # configurations used (3x3x2)
errorThresh = 1e-5

#############################################################################
# Loop over all lines of file and read in tau histogram, grid-mean statistics and 
# subcolumn results.
#############################################################################
hist_ref       = numpy.zeros(shape=(7,7,npoints,nconfig))
cldArea_ref    = numpy.zeros(shape=(npoints,nconfig))
meanPtop_ref   = numpy.zeros(shape=(npoints,nconfig))
meanTauCld_ref = numpy.zeros(shape=(npoints,nconfig))
meanAlbCld_ref = numpy.zeros(shape=(npoints,nconfig))
meanTb_ref     = numpy.zeros(shape=(npoints,nconfig))
meanTbClr_ref  = numpy.zeros(shape=(npoints,nconfig))
hist_new       = numpy.zeros(shape=(7,7,npoints,nconfig))
cldArea_new    = numpy.zeros(shape=(npoints,nconfig))
meanPtop_new   = numpy.zeros(shape=(npoints,nconfig))
meanTauCld_new = numpy.zeros(shape=(npoints,nconfig))
meanAlbCld_new = numpy.zeros(shape=(npoints,nconfig))
meanTb_new     = numpy.zeros(shape=(npoints,nconfig))
meanTbClr_new  = numpy.zeros(shape=(npoints,nconfig))
readLines1     = 'FALSE'
readLines2     = 'FALSE'
nc             = 0 # Configuration counter
np             = 0 # Point counter
errorcount     = 0
for ij in range(0,nlinesRef):

    # Pull out cloud-top height direction type
    if("adjusted top IR only" in flinesRef[ij]):
        top=int(flinesRef[ij+1].split()[0])

    # Determine when and where to stop reading
    if ("START" in flinesRef[ij] or readLines1 == "TRUE" ):
        readLines1 = "TRUE"
    if ("RESUME" in flinesRef[ij]):
        readLines1 = "FALSE"
        readLines2 = "TRUE"
    if ("STOP" in flinesRef[ij]):
        readLines2 = "FALSE"

    # Read in joint-histogram of tau amd gridbox mean statistics.
    if ("taumin" in flinesRef[ij] and readLines1 == "TRUE"):
        for i1 in range(1,8):
            hist_ref[i1-1,:,np,nc]=flinesRef[ij+i1*2+1].split()
            hist_new[i1-1,:,np,nc]=flinesNew[ij+i1*2+1].split()

        hist_test1 = numpy.max(abs(hist_ref[:,:,np,nc]-hist_new[:,:,np,nc]))
        if (hist_test1 != 0.0):
            print '####################################################################################'
            print "          Differences in joint histograms for configuration #",nc+1
            print "REFERENCE (@line",ij,"): "
            print hist_ref[:,:,np,nc]
            print "NEW (@line",ij,"):       "
            print hist_new[:,:,np,nc]
            errorcount=errorcount+1


        if (top != 2):
            cldArea_ref[np,nc]    = float(flinesRef[ij+22].split()[2])
            meanPtop_ref[np,nc]   = float(flinesRef[ij+23].split()[2])
            meanTauCld_ref[np,nc] = float(flinesRef[ij+24].split()[2])
            meanAlbCld_ref[np,nc] = float(flinesRef[ij+25].split()[2])
            meanTb_ref[np,nc]     = float(flinesRef[ij+26].split()[2])
            meanTbClr_ref[np,nc]  = float(flinesRef[ij+27].split()[2])
            cldArea_new[np,nc]    = float(flinesNew[ij+22].split()[2])
            meanPtop_new[np,nc]   = float(flinesNew[ij+23].split()[2])
            meanTauCld_new[np,nc] = float(flinesNew[ij+24].split()[2])
            meanAlbCld_new[np,nc] = float(flinesNew[ij+25].split()[2])
            meanTb_new[np,nc]     = float(flinesNew[ij+26].split()[2])
            meanTbClr_new[np,nc]  = float(flinesNew[ij+27].split()[2])
        if (top == 2):
            cldArea_ref[np,nc]    = float(flinesRef[ij+20].split()[2])
            meanPtop_ref[np,nc]   = float(flinesRef[ij+21].split()[2])
            meanTauCld_ref[np,nc] = float(flinesRef[ij+22].split()[2])
            meanAlbCld_ref[np,nc] = float(flinesRef[ij+23].split()[2])
            meanTb_ref[np,nc]     = -999.00
            meanTbClr_ref[np,nc]  = -999.00
            cldArea_new[np,nc]    = float(flinesNew[ij+20].split()[2])
            meanPtop_new[np,nc]   = float(flinesNew[ij+21].split()[2])
            meanTauCld_new[np,nc] = float(flinesNew[ij+22].split()[2])
            meanAlbCld_new[np,nc] = float(flinesNew[ij+23].split()[2])
            meanTb_new[np,nc]     = -999.00
            meanTbClr_new[np,nc]  = -999.00
 
        # Compare grid-mean statistics
        stat_test1 = numpy.max(abs(cldArea_ref[np,nc]   - cldArea_new[np,nc]))/abs(cldArea_ref[np,nc])
        stat_test2 = numpy.max(abs(meanPtop_ref[np,nc]  - meanPtop_new[np,nc]))/abs(meanPtop_ref[np,nc])
        stat_test3 = numpy.max(abs(meanTauCld_ref[np,nc]- meanTauCld_new[np,nc]))/abs(meanTauCld_ref[np,nc])
        stat_test4 = numpy.max(abs(meanAlbCld_ref[np,nc]- meanAlbCld_new[np,nc]))/abs(meanAlbCld_ref[np,nc])
        stat_test5 = numpy.max(abs(meanTb_ref[np,nc]    - meanTb_new[np,nc]))/abs(meanTb_ref[np,nc])
        stat_test6 = numpy.max(abs(meanTbClr_ref[np,nc] - meanTbClr_new[np,nc]))/abs(meanTbClr_ref[np,nc])
        if (stat_test1 > errorThresh):
            print '####################################################################################'
            print "          Differences in grid-mean statistics (totalcldarea) for configuration #",nc+1
            print "REFERENCE (@line",ij,"): ",cldArea_ref[np,nc]
            print "NEW (@line",ij,"):       ",cldArea_new[np,nc]
            errorcount=errorcount+1
        if (stat_test2 > errorThresh):
            print '####################################################################################'
            print "          Differences in grid-mean statistics (meanptop) for configuration #",nc+1
            print "REFERENCE (@line",ij,"): ",meanPtop_ref[np,nc]
            print "NEW (@line",ij,"):       ",meanPtop_new[np,nc]
            errorcount=errorcount+1
        if (stat_test3 > errorThresh):
            print '####################################################################################'
            print "          Differences in grid-mean statistics (meantaucld) for configuration #",nc+1
            print "REFERENCE (@line",ij,"): ",meanTauCld_ref[np,nc]
            print "NEW (@line",ij,"):       ",meanTauCld_new[np,nc]
            errorcount=errorcount+1
        if (stat_test4 > errorThresh):
            print '####################################################################################'
            print "          Differences in grid-mean statistics (meanalbedocld) for configuration #",nc+1
            print "REFERENCE (@line",ij,"): ",meanAlbCld_ref[np,nc]
            print "NEW (@line",ij,"):       ",meanAlbCld_new[np,nc]
            errorcount=errorcount+1
        if (stat_test5 > errorThresh):
            print '####################################################################################'
            print "          Differences in grid-mean statistics (meantb) for configuration #",nc+1
            print "REFERENCE (@line",ij,"): ",meanTb_ref[np,nc]
            print "NEW (@line",ij,"):       ",meanTb_new[np,nc]
            errorcount=errorcount+1
        if (stat_test6 > errorThresh):
            print '####################################################################################'
            print "          Differences in grid-mean statistics (meantbclr) for configuration #",nc+1
            print "REFERENCE (@line",ij,"): ",meanTbClr_ref[np,nc]
            print "NEW (@line",ij,"):       ",meanTbClr_new[np,nc]
            errorcount=errorcount+1

        np = np+1
        # Reset point counter and iterate configuration counter
        if (np == 24):
            nc = nc+1
            np = 0

    # Compare cloud subcolumn representation 
    if (readLines2 == "TRUE"):
        if (flinesRef[ij] != flinesNew[ij]):
            print '####################################################################################'
            print "          Differences in subcolumn results for configuration #",nc 
            print "REFERENCE (@line",ij,"): ",flinesRef[ij],
            print "NEW (@line",ij,"):       ",flinesNew[ij],
            errorcount = errorcount+1

#############################################################################
# Print summary information
#############################################################################
if (errorcount == 0):
    print 'Using and error-threshold of ',errorThresh
    print "NO DIFFERENCES WERE FOUND BETWEEN THE REFERENCE FILE AND THE NEW FILE"
if (errorcount != 0):
    print '####################################################################################'
    print 'Using and error-threshold of ',errorThresh
    print errorcount," DIFFERENCES WERE FOUND"



 


            
