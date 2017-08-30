#!/bin/sh
# sanity checks that all the gunzipped files in the working directory are in working order

echoerr() { echo "$@" 1>&2; }

for file in *.gz; do
	echoerr "${file}"
	bgzip -cd "${file}" | wc
done
