#!/bin/sh
# Given a file with a list of links,
# downloads each one sequentially and saves it with a unique name:
# Example: SGDP_POP_OriginalFileHandle.vcf.gz.tbi
#
# AUTHOR: Rocky Aikens rockyaikens@gmail.com
# USAGE: wget_from_file.sh FILENAME POP

if [ $# -eq 0 ]
  then
    echo "No arguments supplied."
    echo "# Given a file with a list of links,
# downloads each one sequentially and saves it with a unique name:
# Example: SGDP_POP_OriginalFileHandle.vcf.gz.tbi
#
# AUTHOR: Rocky Aikens rockyaikens@gmail.com
# USAGE: wget_from_file.sh FILENAME POP"
    exit 1
fi

FILENAME=${1}
POP=${2}

while read LINK; do
    echo "Getting from link:"
    echo ${LINK}
    echo "Saving to file:"
    b=${LINK%%\?*}
    INFILE=${b##*/}
    OUTFILE="SGDP_${POP}_${INFILE}"
    echo ${OUTFILE}
    wget -O ${OUTFILE} "${LINK}"
    if [[ "$?" != 0 ]]; then
        echo "Error downloading ${INFILE}"
    else
        echo "Successfully downloaded ${INFILE} to ${OUTFILE}"
    fi
    
done <${FILENAME}
