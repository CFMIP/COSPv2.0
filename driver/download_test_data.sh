DRIVER_DIR=${PWD}
GDFILE='https://docs.google.com/uc?export=download&id=17eK4_DVEvFOE9Uf6siXJDpWZJKT1aqkU'
OUTPATH=data/inputs/UKMO/cosp_input.um_global.nc.gz
wget --no-check-certificate $GDFILE -O $OUTPATH
gunzip ${OUTPATH}
cd data/inputs/UKMO
md5sum -c cosp_input.um_global.nc.md5
cd ${DRIVER_DIR}
GDFILE='https://docs.google.com/uc?export=download&id=10fjcxnmHpt8Go6ipHWUEnN_Siwdtdaqb'
OUTPATH=data/outputs/UKMO/cosp2_output.um_global.gfortran.kgo.nc.gz
wget --no-check-certificate $GDFILE -O $OUTPATH
gunzip ${OUTPATH}
cd data/outputs/UKMO
md5sum -c cosp2_output.um_global.gfortran.kgo.nc.md5


