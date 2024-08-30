gcc -g -Wall -D__USE_FIXED_PROTOTYPES__ -O2 -I klib   -c -o oligotm.o oligotm.c
ar rv liboligotm.a oligotm.o
gcc -g -Wall -D__USE_FIXED_PROTOTYPES__ -O2 -I klib -o oligotm oligotm_main.c liboligotm.a  -lm
./oligotm -fo 0.8 TAGCTAGCTAGCTAGCTATGCTATCG 
echo -------python out-------------------
python ./p3_f.py