#!/bin/bash

motif_file="Nvec_predicted_motifs_pwm.txt"
cisbp_base="CisBP"
output_dir="PWMs"
missing_file="missing_motifs.txt"

mkdir -p "$output_dir"
> "$missing_file"  # Empty the file at the start

while read -r motif; do
    if [[ -z "$motif" ]]; then
        continue
    fi
    
    found=0
    for species_dir in "$cisbp_base"/*/pwms_all_motifs; do
        pwm_file="${species_dir}/${motif}.txt"
        if [[ -f "$pwm_file" ]]; then
            cp "$pwm_file" "${output_dir}/${motif}.txt"
            echo "✅ Found and copied $motif from $species_dir"
            found=1
            break
        fi
    done
    
    if [[ $found -eq 0 ]]; then
        echo "❌ Motif $motif not found"
        echo "$motif" >> "$missing_file"
    fi

done < "$motif_file"

echo "🎉 Done! Missing motifs are logged in $missing_file"

