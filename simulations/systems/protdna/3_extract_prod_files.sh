#!/bin/bash

if [ -z "$SIMULATIONS_DIR" ]; then
    echo "SIMULATIONS_DIR variable is not defined. Exiting."
    exit 1
fi

if [ -z "$NAME" ]; then
    echo "NAME variable is not defined. Exiting."
    exit 1
fi

folder="$SIMULATIONS_DIR/systems/$NAME"
folder_ions="$folder/3-adding-ions"
folder_equil="$folder/6-equilibration"

folder_prod="$folder/files-for-prod"

mkdir -p "$folder_prod"

cp -r "$FF" "$folder_prod"

cp "$folder_ions"/main.top "$folder_prod/top.top"

cp "$folder_ions"/atomtypes.itp "$folder_prod/atomtypes.itp"
cp "$folder_ions"/system1.itp "$folder_prod/system1.itp"
cp "$folder_ions"/system1_posre.itp "$folder_prod/system1_posre.itp"
cp "$folder_ions"/system2.itp "$folder_prod/system2.itp"
cp "$folder_ions"/system2_posre.itp "$folder_prod/system2_posre.itp"
cp "$folder_ions"/system3.itp "$folder_prod/system3.itp"
cp "$folder_ions"/system3_posre.itp "$folder_prod/system3_posre.itp"
cp "$folder_ions"/system4.itp "$folder_prod/system4.itp"
cp "$folder_ions"/system4_posre.itp "$folder_prod/system4_posre.itp"
cp "$folder_ions"/system5.itp "$folder_prod/system5.itp"
cp "$folder_ions"/system5_posre.itp "$folder_prod/system5_posre.itp"
cp "$folder_ions"/system6.itp "$folder_prod/system6.itp"
cp "$folder_ions"/system6_posre.itp "$folder_prod/system6_posre.itp"
cp "$folder_ions"/system7.itp "$folder_prod/system7.itp"
cp "$folder_ions"/system7_posre.itp "$folder_prod/system7_posre.itp"
cp "$folder_ions"/system8.itp "$folder_prod/system8.itp"
cp "$folder_ions"/system8_posre.itp "$folder_prod/system8_posre.itp"
cp "$folder_ions"/system11.itp "$folder_prod/system11.itp"
cp "$folder_ions"/system11_posre.itp "$folder_prod/system11_posre.itp"
cp "$folder_ions"/system12.itp "$folder_prod/system12.itp"
cp "$folder_ions"/system12_posre.itp "$folder_prod/system12_posre.itp"
cp "$folder_ions"/system13.itp "$folder_prod/system13.itp"
cp "$folder_ions"/system13_posre.itp "$folder_prod/system13_posre.itp"
cp "$folder_ions"/system14.itp "$folder_prod/system14.itp"
cp "$folder_ions"/system14_posre.itp "$folder_prod/system14_posre.itp"
cp "$folder_ions"/system5a.itp "$folder_prod/system5a.itp"
cp "$folder_ions"/system5a_posre.itp "$folder_prod/system5a_posre.itp"
cp "$folder_ions"/system15.itp "$folder_prod/system15.itp"
cp "$folder_ions"/system15_posre.itp "$folder_prod/system15_posre.itp"
cp "$folder_ions"/system7a.itp "$folder_prod/system7a.itp"
cp "$folder_ions"/system7a_posre.itp "$folder_prod/system7a_posre.itp"
cp "$folder_ions"/system16.itp "$folder_prod/system16.itp"
cp "$folder_ions"/system16_posre.itp "$folder_prod/system16_posre.itp"
cp "$folder_ions"/dna1.itp "$folder_prod/dna1.itp"
cp "$folder_ions"/dna1_posre.itp "$folder_prod/dna1_posre.itp"
cp "$folder_ions"/dna2.itp "$folder_prod/dna2.itp"
cp "$folder_ions"/dna2_posre.itp "$folder_prod/dna2_posre.itp"
cp "$folder_ions"/dna1a.itp "$folder_prod/dna1a.itp"
cp "$folder_ions"/dna1a_posre.itp "$folder_prod/dna1a_posre.itp"
cp "$folder_ions"/dna2a.itp "$folder_prod/dna2a.itp"
cp "$folder_ions"/dna2a_posre.itp "$folder_prod/dna2a_posre.itp"
cp "$folder_ions"/opc4.itp "$folder_prod/opc4.itp"

cp "$folder_equil"/state.cpt "$folder_prod/cpt.cpt"
cp "$folder_equil"/equil3.gro "$folder_prod/gro.gro"
