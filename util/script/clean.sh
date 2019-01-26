#!/bin/sh

(cd init; make clean)
(cd src; make clean)
rm -fv DAT/* rxmd *.txt

<<<<<<< HEAD:util/script/clean.sh
find regtests/*/run/DAT -type f ! -name .gitignore | xargs rm -v
=======
find regtests/*/run/DAT -type f ! -name .gitignore | xargs rm -fv
>>>>>>> master:util/script/clean.sh
