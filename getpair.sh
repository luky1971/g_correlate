#!/bin/bash

# Exit if there's an error
set -e

# Check arguments
if [ "$#" -lt 2 ]; then
    echo "ERROR: Need at least 2 arguments: filename, first pair, [last pair]"
    exit 1
fi

FILE="$1"
FIRST="$2"
LAST="$3"

if [ "$#" -lt 3 ]; then
    LAST="$FIRST"
fi

# Get filename without extra fluff
FNAME=$(echo "$FILE" | rev | cut -d'/' -f 1 | cut -d'.' -f 2 | rev)

awk -F '[, ]' -v file="$FNAME" -v first="$FIRST" -v last="$LAST" '/PAIR/{pair=$3}{if (pair >= first && pair <= last) print > (file "_p" pair ".dat")}' "$FILE"

echo "Separated files saved as "$FNAME"_p#.dat"
