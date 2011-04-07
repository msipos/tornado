#!/bin/bash

# Simple c-linkage script.  Calls c-linkage on a FASTA file. Rest of
# parameters are passed to c-linkage.
FASTA=$1
shift
REST=$@
ROOTDIR=`dirname $0`
SPLICER="python2 $ROOTDIR/splicer/splicer"
PSIDIST="$ROOTDIR/psi-distance/psi-distance"
CLINKAGE="$ROOTDIR/c-linkage/c-linkage"

if [ ! -x $PSIDIST ]; then
  echo "no $PSIDIST"
  cd $ROOTDIR/psi-distance
  ./build.sh
  cd -
  if [ ! -x $PSIDIST ]; then
    echo "still no $PSIDIST"
    exit 1
  fi
fi

if [ ! -x $CLINKAGE ]; then
  echo "no $CLINKAGE"
  cd $ROOTDIR/c-linkage
  ./build.sh
  cd -
  if [ ! -x $CLINKAGE ]; then
    echo "still no $CLINKAGE"
    exit 1
  fi
fi

$SPLICER dereplicate $FASTA
$PSIDIST $FASTA.derepl > $FASTA.derepl.dist
$CLINKAGE $FASTA.derepl.dist $FASTA.names $REST
rm $FASTA.derepl
rm $FASTA.derepl.dist
rm $FASTA.names
