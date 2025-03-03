#!/bin/bash
#------------------------------------------------------------------------------------------------------------------
mv SHiELD_SRC SHiELD_SRC_temp

# Set up SHiELD_build
echo "------------------------------------------------------------------------------------------"
git clone https://github.com/NOAA-GFDL/SHiELD_build.git
echo "------------------------------------------------------------------------------------------"

# initialize the mkmf submodule
cd SHiELD_build
git submodule update --init mkmf

./CHECKOUT_code

cd ..

rm -rf SHiELD_SRC/GFDL_atmos_cubed_sphere
mv SHiELD_SRC_temp/GFDL_atmos_cubed_sphere SHiELD_SRC/GFDL_atmos_cubed_sphere

rm -rf SHiELD_SRC/atmos_drivers
mv SHiELD_SRC_temp/atmos_drivers SHiELD_SRC/atmos_drivers
rm -rf SHiELD_SRC_temp


