KGO_VERSION=v001
# Input: global UM
DRIVER_DIR=${PWD}
GDFILE='https://docs.google.com/uc?export=download&id=17eK4_DVEvFOE9Uf6siXJDpWZJKT1aqkU'
OUTPATH=data/inputs/UKMO/cosp_input.um_global.nc.gz
wget --no-check-certificate $GDFILE -O $OUTPATH
gunzip ${OUTPATH}
cd data/inputs/UKMO
md5sum -c cosp_input.um_global.nc.md5
# KGO: global UM
cd ${DRIVER_DIR}
GDFILE='https://docs.google.com/uc?export=download&id=1uQBPUEXlniQWEp2nU3iC6d8CO3GM9cvP'
OUTPATH=data/outputs/UKMO/cosp2_output.um_global.gfortran.kgo.$KGO_VERSION.nc.gz
wget --no-check-certificate $GDFILE -O $OUTPATH
gunzip ${OUTPATH}
cd data/outputs/UKMO
md5sum -c cosp2_output.um_global.gfortran.kgo.$KGO_VERSION.nc.md5
# KGO:  UM
cd ${DRIVER_DIR}
GDFILE='https://docs.google.com/uc?export=download&id=1gSEdJJpqhfElsFNcTF_r4A_0vIMEWGla'
OUTPATH=data/outputs/UKMO/cosp2_output_um.gfortran.kgo.$KGO_VERSION.nc.gz
wget --no-check-certificate $GDFILE -O $OUTPATH
gunzip ${OUTPATH}
cd data/outputs/UKMO
md5sum -c cosp2_output_um.gfortran.kgo.$KGO_VERSION.nc.md5
