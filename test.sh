@!/bin/bash

make clean

make -f revmakefile lib

read -p 'Press Enter '

python3 test.py

read -p 'Press Enter '

