#!/bin/bash

if [ -z "$ROOT" ]; then
    echo "ROOT variable is not defined. Exiting."
    exit 1
fi

if [ -z "$NAME" ]; then
    echo "NAME variable is not defined. Exiting."
    exit 1
fi

folder="$ROOT/systems/$NAME"
folder_prep="$folder/1-preparation"
folder_equil="$folder/2-equilibration"

folder_prod="$folder/files-for-prod"

mkdir -p "$folder_prod"

cp -r "$FF" "$folder_prod"
cp "$folder_prep"/topol.top "$folder_prod/top.top"
cp "$folder_prep"/index.ndx "$folder_prod/ndx.ndx"
cp "$folder_equil"/step6.6_equilibration.cpt "$folder_prod/cpt.cpt"
cp "$folder_equil"/step6.6_equilibration.gro "$folder_prod/gro.gro"

cp -r "$folder_prep"/toppar "$folder_prod"