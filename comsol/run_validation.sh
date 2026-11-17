#!/bin/bash
# ptpipe - run the COMSOL templates on the validation cases (macOS / Linux terminal).
#
#   bash run_validation.sh                  (quick set: sphere_R20 mim_s0 mim_device)
#   bash run_validation.sh all              (every case in configs/)
#   bash run_validation.sh rod_AR4 dimer_g05
#   bash run_validation.sh --cores 6 rod_AR4  (leave cores free to keep using the computer)
#
# Then, from the repository root:   python validation/validate.py
#
# COMSOL executable: set COMSOL=/path/to/comsol if it is not in the default location.

cd "$(dirname "$0")"
C="${COMSOL:-/Applications/COMSOL64/Multiphysics/bin/comsol}"
ROOT="$(cd .. && pwd)/"
mkdir -p ../results

template_of() {   # case -> template, from the keys of its case file (any case name works)
  f="configs/$1.cfg"
  [ -f "$f" ] || { echo ""; return; }
  if grep -qE '^[[:space:]]*template[[:space:]]*=' "$f"; then
    grep -E '^[[:space:]]*template[[:space:]]*=' "$f" | head -1 | sed -E 's/^[^=]*=[[:space:]]*([A-Za-z_]+).*/\1/'
  elif grep -qE '^[[:space:]]*(shape|cluster_file|R_np)[[:space:]]*=' "$f"; then echo PT_Particle
  elif grep -qE '^[[:space:]]*layers[[:space:]]*=' "$f"; then echo PT_Multilayer
  elif grep -qE '^[[:space:]]*A_abs[[:space:]]*=' "$f"; then echo PT_Device
  else echo ""; fi
}
cfg_of() {
  case "$1" in PT_Particle) echo particle.cfg ;; PT_Multilayer) echo multilayer.cfg ;; PT_Device) echo device.cfg ;; esac
}

NP=""
if [ "$1" = "--cores" ]; then NP="-np $2"; shift 2; fi

if [ "$1" = "figures" ]; then
  # export meshes and field maps of solved cases (no solving): bash run_validation.sh figures [case ...]
  shift
  if [ $# -gt 0 ]; then printf "%s\n" "$@" > figures.cfg
  else ls ../results/*_Model.mph 2>/dev/null | xargs -n1 basename | sed 's/_Model\.mph$//' > figures.cfg; fi
  mkdir -p ../results/build
  perl -pe "s|static String BASE = \".*\";|static String BASE = \"$ROOT\";|" PT_Figures.java > ../results/build/PT_Figures.java
  "$C" compile ../results/build/PT_Figures.java
  sleep 3
  "$C" batch $NP -inputfile ../results/build/PT_Figures.class -batchlog ../results/figures_batch.log
  echo "Done. Next: python validation/plot_fields.py \$(cat figures.cfg)"
  exit 0
fi

if [ $# -eq 0 ]; then CASES="sphere_R20 mim_s0 mim_device"
elif [ "$1" = "all" ]; then CASES=$(ls configs/*.cfg | xargs -n1 basename | sed 's/\.cfg$//')
else
  CASES=""
  for a in "$@"; do
    case "$a" in \#*) break ;; esac            # anything after a '#' is a comment
    CASES="$CASES $a"
  done
fi

KEEP=""; command -v caffeinate >/dev/null && KEEP="caffeinate -i"
COMPILED=" "
for c in $CASES; do
  T=$(template_of "$c")
  if [ -z "$T" ] || [ ! -f "configs/$c.cfg" ]; then
    echo "!!! skipping '$c': no configs/$c.cfg or template not recognised (available: $(ls configs | sed 's/\.cfg//' | tr '\n' ' '))"
    continue
  fi
  if [[ "$COMPILED" != *" $T "* ]]; then
    # build copy with the repository path written in (sources stay untouched), then compile once
    mkdir -p ../results/build
    perl -pe "s|static String BASE = \".*\";|static String BASE = \"$ROOT\";|" "$T.java" > "../results/build/$T.java"
    "$C" compile "../results/build/$T.java"
    sleep 3            # avoids EOF errors when the class is read back from slow (SD) disks
    COMPILED="$COMPILED$T "
  fi
  cp "configs/$c.cfg" "$(cfg_of $T)"
  echo "=== $c  ($T)  $(date)  ${NP:-all cores}"
  $KEEP "$C" batch $NP -inputfile "../results/build/$T.class" -outputfile "../results/$c.mph" -batchlog "../results/${c}_batch.log" \
    || echo "!!! $c failed: see results/${c}_batch.log and results/${c}_log.txt"
  sleep 3
done
echo "Done. Next: python validation/validate.py"
