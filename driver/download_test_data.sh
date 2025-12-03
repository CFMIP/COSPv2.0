#/bin/bash

COMPILER=$1

case $COMPILER in

  "gfortran")
    KGO_VERSION=v006
    kgo_links=( "https://docs.google.com/uc?export=download&id=1HgtrBVI5-7ypQWxzzMPvAq_5IIE5q9-B" \
                "https://docs.google.com/uc?export=download&id=17Ah7z2oGiwwZ-yW1JX3Vg81_l7Q1zFf7" \
                "https://docs.google.com/uc?export=download&id=1IzBbnhnfYeaSPLRSvWsWCEzsL_D9_ZKq" \
                "https://docs.google.com/uc?export=download&id=1qqwHbMSuD38GcIVIbqsmJehHjSQEW2Mk" )
    out_type=( "cosp2_output_um.${COMPILER}.kgo" \
               "cosp2_output.um_global.${COMPILER}.kgo" \
               "cosp2_output.um_global_model_levels.${COMPILER}.kgo" \
               "cosp2_swath_output.um_global.${COMPILER}.kgo" )
    ;;

  "ifort")
    KGO_VERSION=v006
    kgo_links=( "https://docs.google.com/uc?export=download&id=121bSDuGNdbkb9WKhJLu9Pwe8KnhUMtDx" \
                "https://docs.google.com/uc?export=download&id=1JwWJKw8rO1MpDLXFvdWvrqKF5_HbKZAK" \
                "https://docs.google.com/uc?export=download&id=1bkLKYHQqZdskVZsCH_XiWIayf0PFzz5U" )
    out_type=( "cosp2_output_um.${COMPILER}.kgo" \
               "cosp2_output.um_global.${COMPILER}.kgo" \
               "cosp2_swath_output.um_global.${COMPILER}.kgo" )
    ;;

  "ifx")
    KGO_VERSION=v007
    kgo_links=( "https://docs.google.com/uc?export=download&id=1epHiss6ESy2rDuBsQwXn-8KAyjmAip6j" \
                "https://docs.google.com/uc?export=download&id=16KKAP_GFmY30d-_NR90Y_Lj0zlisYDRH" \
                "https://docs.google.com/uc?export=download&id=1llJMOnQRIws8DLkEs8svTO8GTM_29DF6" )
    out_type=( "cosp2_output_um.${COMPILER}.kgo" \
               "cosp2_output.um_global.${COMPILER}.kgo" \
               "cosp2_swath_output.um_global.${COMPILER}.kgo" )
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
