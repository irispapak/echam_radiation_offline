#! /bin/ksh
#
#         
#
set -ex
#
EXP=2xCO2_alb
#
cdo cat ../../Results_thorsten/${EXP}_flsh_1940_lev1.nc ../../Results_thorsten/${EXP}_flsh_lev1.nc
cdo cat ../../Results_thorsten/${EXP}_flsh_1941_lev1.nc ../../Results_thorsten/${EXP}_flsh_lev1.nc
cdo cat ../../Results_thorsten/${EXP}_flsh_1942_lev1.nc ../../Results_thorsten/${EXP}_flsh_lev1.nc
cdo cat ../../Results_thorsten/${EXP}_flsh_1943_lev1.nc ../../Results_thorsten/${EXP}_flsh_lev1.nc
cdo cat ../../Results_thorsten/${EXP}_flsh_1944_lev1.nc ../../Results_thorsten/${EXP}_flsh_lev1.nc
cdo cat ../../Results_thorsten/${EXP}_flsh_1945_lev1.nc ../../Results_thorsten/${EXP}_flsh_lev1.nc 
#
cdo cat ../../Results_thorsten/${EXP}_flth_1940_lev1.nc ../../Results_thorsten/${EXP}_flth_lev1.nc
cdo cat ../../Results_thorsten/${EXP}_flth_1941_lev1.nc ../../Results_thorsten/${EXP}_flth_lev1.nc
cdo cat ../../Results_thorsten/${EXP}_flth_1942_lev1.nc ../../Results_thorsten/${EXP}_flth_lev1.nc
cdo cat ../../Results_thorsten/${EXP}_flth_1943_lev1.nc ../../Results_thorsten/${EXP}_flth_lev1.nc
cdo cat ../../Results_thorsten/${EXP}_flth_1944_lev1.nc ../../Results_thorsten/${EXP}_flth_lev1.nc
cdo cat ../../Results_thorsten/${EXP}_flth_1945_lev1.nc ../../Results_thorsten/${EXP}_flth_lev1.nc
#
#
cdo settaxis,1940-01-01,18:00:00,1days -setcalendar,365days ../../Results_thorsten/${EXP}_flth_lev1.nc ../../Results_thorsten/${EXP}_hilf_flth.nc
cdo settaxis,1940-01-01,18:00:00,1days -setcalendar,365days ../../Results_thorsten/${EXP}_flsh_lev1.nc ../../Results_thorsten/${EXP}_hilf_flsh.nc 
#
cdo timmean ../../Results_thorsten/${EXP}_hilf_flsh.nc ../../Results_thorsten/${EXP}_flsh_lev1_timm.nc
cdo timmean ../../Results_thorsten/${EXP}_hilf_flth.nc ../../Results_thorsten/${EXP}_flth_lev1_timm.nc
#
cdo fldmean ../../Results_thorsten/${EXP}_flsh_lev1_timm.nc ../../Results_thorsten/${EXP}_flsh_lev1_timm_fldmean.nc
cdo fldmean ../../Results_thorsten/${EXP}_flth_lev1_timm.nc ../../Results_thorsten/${EXP}_flth_lev1_timm_fldmean.nc
#
rm ../../Results_thorsten/${EXP}_hilf_flsh.nc
rm ../../Results_thorsten/${EXP}_hilf_flth.nc
#
exit
