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
folder_ions="$folder/3-adding-ions"
folder_equil="$folder/6-equilibration"

folder_prod="$folder/files-for-prod"

mkdir -p "$folder_prod"

cp -r "$FF" "$folder_prod"
cp "$folder_ions"/topol.top "$folder_prod/top.top"
cp "$folder_equil"/state.cpt "$folder_prod/cpt.cpt"
cp "$folder_equil"/equil3.gro "$folder_prod/gro.gro"
