#/bin/bash

KGO_VERSION=v004
echo "KGO_VERSION=$KGO_VERSION"

kgo_links=( "https://docs.google.com/uc?export=download&id=1kC9RViPBdAsGcOivpYXxs3hHJ11Kv8IN" \
            "https://docs.google.com/uc?export=download&id=1X_oOzvY2lf-kyAR-D1E6JfuefkGg1idn" \
            "https://docs.google.com/uc?export=download&id=1c14qBf9VwYJWYVGCu-Cw35F-qqHx0mSg" )

out_type=( "cosp2_output_um.gfortran.kgo" \
           "cosp2_output.um_global.gfortran.kgo" \
           "cosp2_output.um_global_model_levels.gfortran.kgo" )

# Input: global UM
DRIVER_DIR=${PWD}
GDFILE='https://docs.google.com/uc?export=download&id=17eK4_DVEvFOE9Uf6siXJDpWZJKT1aqkU'
OUTPATH=data/inputs/UKMO/cosp_input.um_global.nc.gz
wget --no-check-certificate $GDFILE -O $OUTPATH
gunzip -f ${OUTPATH}
cd data/inputs/UKMO
md5sum -c cosp_input.um_global.nc.md5

for i in ${!kgo_links[@]}; do
  cd ${DRIVER_DIR}
  GDFILE=${kgo_links[$i]}
  OUTPATH=data/outputs/UKMO/${out_type[$i]}.$KGO_VERSION.nc.gz
  wget --no-check-certificate $GDFILE -O $OUTPATH
  gunzip -f ${OUTPATH}
  cd data/outputs/UKMO
  md5sum -c ${out_type[$i]}.$KGO_VERSION.nc.md5
done
