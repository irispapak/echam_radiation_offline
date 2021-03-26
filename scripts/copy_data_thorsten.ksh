#! /bin/ksh
#
#         
#
set -ex
#
#
typeset -Z4 YY EY SY
#
#
#
#
SY=1950      #Startjahr
EY=1950      #Endjahr
#
#
#
YY=${SY}
#
while [ ${YY} -le ${EY} ]
do
    scp blizzard.dkrz.de:/work/mh0730/m300057/experiments/feedback/CTRL/CTRL_${YY}_echam.grb /scratch/mpi/bm0546/m222079/Radiation/Radiation_standalone_2011_07_18/data/
    cd /scratch/mpi/bm0546/m222079/Radiation/Radiation_standalone_2011_07_18/data/
    cdo selcode,134 CTRL_${YY}_echam.grb code134_${YY}.grb
    cdo -f nc copy code134_${YY}.grb code134_${YY}.nc
    rm code134_${YY}.grb
#
    cdo selcode,153 CTRL_${YY}_echam.grb code153_${YY}.grb
    cdo -f nc copy code153_${YY}.grb code153_${YY}.nc
    rm code153_${YY}.grb
#
    cdo selcode,154 CTRL_${YY}_echam.grb code154_${YY}.grb
    cdo -f nc copy code154_${YY}.grb code154_${YY}.nc
    rm code154_${YY}.grb
#
    cdo selcode,129 CTRL_${YY}_echam.grb code129_${YY}.grb
    cdo -f nc copy code129_${YY}.grb code129_${YY}.nc
    rm code129_${YY}.grb
#
    cdo selcode,157 CTRL_${YY}_echam.grb code157_${YY}.grb
    cdo -f nc copy code157_${YY}.grb code157_${YY}.nc
    rm code157_${YY}.grb
#
    cdo selcode,175 CTRL_${YY}_echam.grb code175_${YY}.grb
    cdo -f nc copy code175_${YY}.grb code175_${YY}.nc
    rm code175_${YY}.grb
#
#cloud droplet number concentration fehlt leider
#    cdo selcode,136 CTRL_${YY}_echam.grb code136_${YY}.grb
#    cdo -f nc copy code137_${YY}.grb code136_${YY}.nc
#    rm code136_${YY}.grb
#
    cdo selcode,162 CTRL_${YY}_echam.grb code162_${YY}.grb
    cdo -f nc copy code162_${YY}.grb code162_${YY}.nc
    rm code162_${YY}.grb
#
    cdo selcode,169 CTRL_${YY}_echam.grb code169_${YY}.grb
    cdo -f nc copy code169_${YY}.grb code169_${YY}.nc
    rm code169_${YY}.grb
#
    cdo selcode,236 CTRL_${YY}_echam.grb code236_${YY}.grb
    cdo -f nc copy code236_${YY}.grb code236_${YY}.nc
    rm code236_${YY}.grb
#
    cdo selcode,237 CTRL_${YY}_echam.grb code237_${YY}.grb
    cdo -f nc copy code237_${YY}.grb code237_${YY}.nc
    rm code237_${YY}.grb
#
    cdo selcode,130 CTRL_${YY}_echam.grb code130_${YY}.grb
    cdo -f nc copy code130_${YY}.grb code130_${YY}.nc
    rm code130_${YY}.grb
#
    cdo selcode,133 CTRL_${YY}_echam.grb code133_${YY}.grb
    cdo -f nc copy code133_${YY}.grb code133_${YY}.nc
    rm code133_${YY}.grb  
#
    cdo merge code134_${YY}.nc code153_${YY}.nc code154_${YY}.nc code129_${YY}.nc code157_${YY}.nc code175_${YY}.nc code136_${YY}.nc code162_${YY}.nc code169_${YY}.nc code236_${YY}.nc code237_${YY}.nc code130_${YY}.nc code133_${YY}.nc CTRL_${YY}.nc
    rm CTRL_${YY}_echam.grb
    rm code134_${YY}.nc
    rm code153_${YY}.nc
    rm code154_${YY}.nc
    rm code129_${YY}.nc
    rm code157_${YY}.nc
    rm code175_${YY}.nc
    rm code162_${YY}.nc
    rm code169_${YY}.nc
    rm code236_${YY}.nc 
    rm code237_${YY}.nc
    rm code130_${YY}.nc
    rm code133_${YY}.nc    
    cdo sp2gp CTRL_${YY}.nc CTRL_${YY}_gg.nc
    rm CTRL_${YY}.nc 
    cdo -f nc -t echam5 -r copy CTRL_${YY}_gg.nc CTRL_${YY}_gg_named.nc
    cdo splitmon CTRL_${YY}_gg_named.nc CTRL_${YY} 
    rm CTRL_${YY}_gg.nc
    rm CTRL_${YY}_gg_named.nc  
    cdo splitday CTRL_${YY}01.nc CTRL_${YY}01
    rm CTRL_${YY}01.nc
    cdo splitday CTRL_${YY}02.nc CTRL_${YY}02
    rm CTRL_${YY}02.nc
    cdo splitday CTRL_${YY}03.nc CTRL_${YY}03
    rm CTRL_${YY}03.nc
    cdo splitday CTRL_${YY}04.nc CTRL_${YY}04
    rm CTRL_${YY}04.nc
    cdo splitday CTRL_${YY}05.nc CTRL_${YY}05
    rm CTRL_${YY}05.nc
    cdo splitday CTRL_${YY}06.nc CTRL_${YY}06
    rm CTRL_${YY}06.nc
    cdo splitday CTRL_${YY}07.nc CTRL_${YY}07
    rm CTRL_${YY}07.nc
    cdo splitday CTRL_${YY}08.nc CTRL_${YY}08
    rm CTRL_${YY}08.nc
    cdo splitday CTRL_${YY}09.nc CTRL_${YY}09
    rm CTRL_${YY}09.nc
    cdo splitday CTRL_${YY}10.nc CTRL_${YY}10
    rm CTRL_${YY}10.nc  
    cdo splitday CTRL_${YY}11.nc CTRL_${YY}11
    rm CTRL_${YY}11.nc
    cdo splitday CTRL_${YY}12.nc CTRL_${YY}12
    rm CTRL_${YY}12.nc
    ((YY=YY+1))
done
#
cd /scratch/mpi/bm0546/m222079/Radiation/Radiation_standalone_2011_07_18/scripts
#
#
exit












