#!/bin/bash
# Molecular docking via AutoDock Vina (config-file mode, Windows-safe)
# 4LXZ=HDAC2(SHH/SAHA)  5D7P=SIRT2(OCZ/EX-243)  6EN4=SF3b(BGZ/macrolide)
set -e
cd /e/KIMI/GBM_WEBS_Project/docking
OBABEL=/e/Miniconda3/envs/docking/Library/bin/obabel.exe
mkdir -p prep results logs

box_center() {  # $1=pdb file, $2=resname -> "x y z"
  awk -v res="$2" '$1=="HETATM" && substr($0,18,3) ~ res {x+=substr($0,31,8); y+=substr($0,39,8); z+=substr($0,47,8); n++}
       END {if(n>0) printf "%.2f %.2f %.2f", x/n, y/n, z/n; else print "0 0 0"}' "$1"
}

prep_receptor() {  # $1=pdb id
  awk '$1=="ATOM" || ($1=="HETATM" && substr($0,18,3)=="ZN ") || $1=="TER" || $1=="END"' pdb/$1.pdb > prep/$1_rec.pdb
  "$OBABEL" -ipdb prep/$1_rec.pdb -O prep/$1_rec.pdbqt -xr -h --partialcharge gasteiger >/dev/null 2>&1
}

prep_ligand() {  # $1=name
  "$OBABEL" ligands/$1_3d.sdf -O prep/$1.pdbqt -h --partialcharge gasteiger >/dev/null 2>&1
}

run_vina() {  # $1=target, $2=ligand name, $3=x $4=y $5=z
  cat > prep/cfg_$1_$2.txt <<EOF
receptor = prep\\$1_rec.pdbqt
ligand = prep\\$2.pdbqt
center_x = $3
center_y = $4
center_z = $5
size_x = 25
size_y = 25
size_z = 25
exhaustiveness = 16
num_modes = 9
out = results\\$1_$2.pdbqt
EOF
  ./vina.exe --config prep/cfg_$1_$2.txt > logs/$1_$2.stdout 2>&1
  aff=$(grep -m1 "REMARK VINA RESULT" results/$1_$2.pdbqt | awk '{print $4}')
  echo "$1,$2,$aff" >> results/binding_energies.csv
  echo "$1 <- $2 : $aff kcal/mol"
}

echo "target,ligand,best_affinity_kcal_mol" > results/binding_energies.csv

run_pair() {  # $1=pdb id, $2=cocrystal resname, $3=drug list
  local id=$1 lig=$2; shift 2
  local center=$(box_center pdb/$id.pdb $lig)
  echo "=== $id (box center: $center) ==="
  prep_receptor $id
  for drug in "$@"; do
    prep_ligand $drug
    run_vina $id $drug $center
  done
}

run_pair 4LXZ SHH vorinostat
run_pair 5D7P OCZ AGK2
run_pair 6EN4 BGZ H3B-8800 pladienolide_B

echo "=== ALL DONE ==="
cat results/binding_energies.csv
