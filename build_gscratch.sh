#!/bin/bash
# build_gscratch.sh
#
# Builds all heliolinx binaries in gscratch to avoid home directory quota issues.
# Object files and final binaries land in BUILDDIR; nothing is written to src/.
#
# Usage:
#   bash build_gscratch.sh              # build all programs
#   bash build_gscratch.sh -j8          # parallel compile (pass -jN to g++)
#   BUILDDIR=/other/path bash build_gscratch.sh
#
set -euo pipefail

SRCDIR="$(cd "$(dirname "$0")/src" && pwd)"
BUILDDIR="${BUILDDIR:-/gscratch/astro/rstrau/heliolinx_build}"
JOBS="${1:-}"   # e.g. -j8

CXX="${CXX:-g++}"
CXXFLAGS="-O3 -std=c++11 -fopenmp -I${SRCDIR}"

PROGRAMS="make_tracklets make_outim heliolinc heliolinc_omp heliolinc_lowmem \
          heliolinc_lowmem_omp heliovane link_purify link_planarity \
          link_planarity_omp parse_clust2det_MPC80 parse_clust2det \
          modsplit_hlfile merge_tracklet_files make_trailed_tracklets \
          parse_trk2det calc_heliohypmat label_hldet helio_highgrade \
          helio_highgrade_omp analyze_linkage01a"

mkdir -p "$BUILDDIR"
cd "$BUILDDIR"

echo "=== Building libheliolinx.a ==="
$CXX $CXXFLAGS -c "${SRCDIR}/solarsyst_dyn_geo01.cpp" -o solarsyst_dyn_geo01.o
ar -rc libheliolinx.a solarsyst_dyn_geo01.o
echo "  libheliolinx.a done"

echo "=== Building programs ==="
compile_one() {
    local prog=$1
    $CXX $CXXFLAGS -c "${SRCDIR}/${prog}.cpp" -o "${prog}.o"
    $CXX "${prog}.o" -L"${BUILDDIR}" -lheliolinx -fopenmp -o "${prog}"
    echo "  ${prog} done"
}
export -f compile_one
export CXX CXXFLAGS BUILDDIR SRCDIR

if command -v parallel &>/dev/null && [[ "$JOBS" == -j* ]]; then
    N="${JOBS#-j}"
    echo "$PROGRAMS" | tr ' ' '\n' | parallel -j"${N}" compile_one {}
else
    for prog in $PROGRAMS; do
        compile_one "$prog"
    done
fi

echo ""
echo "=== Done. Binaries in ${BUILDDIR} ==="
ls -1 "${BUILDDIR}"
