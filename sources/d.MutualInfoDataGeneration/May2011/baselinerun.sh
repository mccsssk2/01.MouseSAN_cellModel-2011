# make directorries and run the san binary
rm san
# cc -o san mousesan_E4031.c -lm
cc -o san mousesan.c -lm
#
./san 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
