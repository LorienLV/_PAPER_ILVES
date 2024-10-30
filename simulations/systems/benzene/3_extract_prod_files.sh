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
folder_prep="$folder/1-preparation"
folder_equil="$folder/3-equilibration"

folder_prod="$folder/files-for-prod"
folder_tpr="$folder/files-for-tpr"

mkdir -p "$folder_prod"

cp -r "$FF" "$folder_prod"
cp "$folder_prep"/${NAME}_ua.itp "$folder_prod/${NAME}_ua.itp"
cp "$folder_prep"/$NAME.top "$folder_prod/top.top"
cp "$folder_equil"/state.cpt "$folder_prod/cpt.cpt"
cp "$folder_equil"/equil3.gro "$folder_prod/gro.gro"
