#/bin/bash

COMPILER=$1

case $COMPILER in

  "gfortran")
    KGO_VERSION=v005
    kgo_links=( "https://docs.google.com/uc?export=download&id=1olgrWzCSpFyPdimFnRfeFxkq5SoIfPDv" \
                "https://docs.google.com/uc?export=download&id=1B9Too1KSfT7fcWKJj2OeR7cytk56eZO2" \
                "https://docs.google.com/uc?export=download&id=1ddBYJsdRXy4z7kLrhRq54KvtuO6AjnTq" )
    out_type=( "cosp2_output_um.${COMPILER}.kgo" \
               "cosp2_output.um_global.${COMPILER}.kgo" \
               "cosp2_output.um_global_model_levels.${COMPILER}.kgo" )
    ;;

  "ifort")
    KGO_VERSION=v006
    kgo_links=( "https://docs.google.com/uc?export=download&id=121bSDuGNdbkb9WKhJLu9Pwe8KnhUMtDx" \
                "https://docs.google.com/uc?export=download&id=1JwWJKw8rO1MpDLXFvdWvrqKF5_HbKZAK" )
    out_type=( "cosp2_output_um.${COMPILER}.kgo" \
               "cosp2_output.um_global.${COMPILER}.kgo" )
    ;;

  "ifx")
    KGO_VERSION=v006
    kgo_links=( "https://docs.google.com/uc?export=download&id=1rby-uiuB8G9vVa4JvRUCYNYPPOG5TmLB" \
                "https://docs.google.com/uc?export=download&id=1SDxQodiDvrwcwfi3ZCMtQITkmHMXIwD-" )
    out_type=( "cosp2_output_um.${COMPILER}.kgo" \
               "cosp2_output.um_global.${COMPILER}.kgo" )
    ;;

  *)
    echo "ERROR: Invalid compiler"
    exit 1
    ;;
esac

echo "KGO_VERSION=$KGO_VERSION"

# Input: global UM
DRIVER_DIR=${PWD}
GDFILE='https://docs.google.com/uc?export=download&id=17eK4_DVEvFOE9Uf6siXJDpWZJKT1aqkU'
OUTPATH=data/inputs/UKMO/cosp_input.um_global.nc.gz
curl -sSfL -o $OUTPATH $GDFILE
gunzip -f ${OUTPATH}
cd data/inputs/UKMO
md5sum -c cosp_input.um_global.nc.md5

for i in ${!kgo_links[@]}; do
  cd ${DRIVER_DIR}
  GDFILE=${kgo_links[$i]}
  OUTPATH=data/outputs/UKMO/${out_type[$i]}.$KGO_VERSION.nc.gz
  curl -sSfL -o $OUTPATH $GDFILE
  gunzip -f ${OUTPATH}
  cd data/outputs/UKMO
  md5sum -c ${out_type[$i]}.$KGO_VERSION.nc.md5
done
