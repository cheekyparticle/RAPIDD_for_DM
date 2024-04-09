mkdir lib/build
cd lib/build
cmake ../
make -j
cd ../../
rm -f pyRAPIDD.c
LDFLAGS="-L$PWD/lib/build" python setup.py build_ext -i
cp lib/build/libRAPIDD.so .
