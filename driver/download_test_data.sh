KGO_VERSION=v002
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
GDFILE='https://docs.google.com/uc?export=download&id=1s5Ha6Hqnv_hWbRUs8KQpJ4Lxy8uvJDar'
OUTPATH=data/outputs/UKMO/cosp2_output.um_global.gfortran.kgo.$KGO_VERSION.nc.gz
wget --no-check-certificate $GDFILE -O $OUTPATH
gunzip ${OUTPATH}
cd data/outputs/UKMO
md5sum -c cosp2_output.um_global.gfortran.kgo.$KGO_VERSION.nc.md5
# KGO:  UM
cd ${DRIVER_DIR}
GDFILE='https://docs.google.com/uc?export=download&id=11dKcIL3EQr7s6jbo4f9GsoW0SufesGbq'
OUTPATH=data/outputs/UKMO/cosp2_output_um.gfortran.kgo.$KGO_VERSION.nc.gz
wget --no-check-certificate $GDFILE -O $OUTPATH
gunzip ${OUTPATH}
cd data/outputs/UKMO
md5sum -c cosp2_output_um.gfortran.kgo.$KGO_VERSION.nc.md5
# KGO: global model levels UM
cd ${DRIVER_DIR}
GDFILE='https://docs.google.com/uc?export=download&id=1kY1lRgzd0UhDiQef2u-VgTQql_iut3U2'
OUTPATH=data/outputs/UKMO/cosp2_output.um_global_model_levels.gfortran.kgo.$KGO_VERSION.nc.gz
wget --no-check-certificate $GDFILE -O $OUTPATH
gunzip ${OUTPATH}
cd data/outputs/UKMO
md5sum -c cosp2_output.um_global_model_levels.gfortran.kgo.$KGO_VERSION.nc.md5
