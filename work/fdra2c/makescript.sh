module purge

# option 1, using openmpi 1.8.8 and the default gcc (version 4.8.5). It compile
module load engaging/openmpi/1.8.8

# option 2, using openmpi/3.1.3 which requires gcc 6.3.0. It doesn't compile
#module load gcc/6.3.0
#module load openmpi/3.1.3

cd external 
git clone git@github.com:camcat/hmmvp.git
cd hmmvp
mkdir lib bin
make
cd ../..
make clean; 
make mode=p;
make mode=p add_bc

ls bin -lh
