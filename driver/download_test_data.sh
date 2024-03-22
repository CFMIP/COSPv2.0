#/bin/bash

if [[ -n ${KGO_VERSION} ]] ; then
    echo "KGO_VERSION=$KGO_VERSION"
else
    export KGO_VERSION=v003
    echo "set KGO_VERSION=$KGO_VERSION"
fi


kgo_links=( "https://docs.google.com/uc?export=download&id=11dKcIL3EQr7s6jbo4f9GsoW0SufesGbq" \
            "https://docs.google.com/uc?export=download&id=1s5Ha6Hqnv_hWbRUs8KQpJ4Lxy8uvJDar" \
            "https://docs.google.com/uc?export=download&id=1kY1lRgzd0UhDiQef2u-VgTQql_iut3U2" )

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

if [[ ${KGO_VERSION} == "v002" ]] ; then
  for i in ${!kgo_links[@]}; do
    cd ${DRIVER_DIR}
    GDFILE=${kgo_links[$i]}
    OUTPATH=data/outputs/UKMO/${out_type[$i]}.$KGO_VERSION.nc.gz
    wget --no-check-certificate $GDFILE -O $OUTPATH
    gunzip -f ${OUTPATH}
    cd data/outputs/UKMO
    md5sum -c ${out_type[$i]}.$KGO_VERSION.nc.md5
  done
else
  echo "wrong KGO_VERSION supplied, try 'export KGO_VERSION=v002'"
fi