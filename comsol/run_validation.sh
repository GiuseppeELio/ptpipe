#!/bin/bash
# ptpipe - run the COMSOL templates on the validation cases (macOS / Linux terminal).
#
#   ./run_validation.sh                 # quick set: sphere_R20 mim_s0 mim_device
#   ./run_validation.sh all             # every case in configs/
#   ./run_validation.sh rod_AR4 dimer_g05
#
# Then, from the repository root:   python validation/validate.py
#
# COMSOL executable: set COMSOL=/path/to/comsol if it is not in the default location.

cd "$(dirname "$0")"
C="${COMSOL:-/Applications/COMSOL64/Multiphysics/bin/comsol}"
ROOT="$(cd .. && pwd)/"
mkdir -p ../results

template_of() {   # case name -> template
  case "$1" in
    sphere_*|rod_*|dimer_*|vesicle_*|cluster_*) echo PT_Particle ;;
    mim_device*|device_*)                       echo PT_Device ;;
    mim_*|stack_*)                              echo PT_Multilayer ;;
    *) echo "" ;;
  esac
}
cfg_of() {
  case "$1" in PT_Particle) echo particle.cfg ;; PT_Multilayer) echo multilayer.cfg ;; PT_Device) echo device.cfg ;; esac
}

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
    echo "!!! skipping '$c': no configs/$c.cfg (available: $(ls configs | sed 's/\.cfg//' | tr '\n' ' '))"
    continue
  fi
  if [[ "$COMPILED" != *" $T "* ]]; then
    # write the repository path into the template, then compile (once per template)
    perl -pi -e "s|static String BASE = \".*\";|static String BASE = \"$ROOT\";|" "$T.java"
    "$C" compile "$T.java"
    sleep 3            # avoids EOF errors when the class is read back from slow (SD) disks
    COMPILED="$COMPILED$T "
  fi
  cp "configs/$c.cfg" "$(cfg_of $T)"
  echo "=== $c  ($T)  $(date)"
  $KEEP "$C" batch -inputfile "$T.class" -outputfile "../results/$c.mph" -batchlog "../results/${c}_batch.log" \
    || echo "!!! $c failed: see results/${c}_batch.log and results/${c}_log.txt"
  sleep 3
done
echo "Done. Next: python validation/validate.py"
