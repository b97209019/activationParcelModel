for COM in AS; do
#for COM in AS NACL HCL S N ; do 
    make clean; make FOO=__${COM}__; mv main main${COM}
done
