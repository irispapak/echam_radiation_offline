#! /bin/ksh
#
#         
#
set -ex
#
#
typeset -Z4 YY 
#
#
#
#
YY=1940      #Jahr
EXP=2xCO2_cld
#
#
#
SD=1
ED=365
#
DD=${SD}
#
while [ ${DD} -le ${ED} ]
do
 cdo setyear,${YY} -timmean -sellevel,1 ../output/flsh_${DD}.nc ../output/flsh_${DD}_lev1.nc
 cdo setyear,${YY} -timmean -sellevel,1 ../output/flth_${DD}.nc ../output/flth_${DD}_lev1.nc
 cdo cat ../output/flsh_${DD}_lev1.nc ../output/${EXP}_flsh_${YY}_lev1.nc
 cdo cat ../output/flth_${DD}_lev1.nc ../output/${EXP}_flth_${YY}_lev1.nc
 rm ../output/flsh_${DD}_lev1.nc
 rm ../output/flth_${DD}_lev1.nc
 rm ../output/flsh_${DD}.nc
 rm ../output/flth_${DD}.nc
    ((DD=DD+1))
done
#
mv ../output/${EXP}_flsh_${YY}_lev1.nc ../../Results_thorsten/
mv ../output/${EXP}_flth_${YY}_lev1.nc ../../Results_thorsten/
#
#
#
exit












