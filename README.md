### Maggie Kerr

### Mount Allison University

### August 2021

### This repository includes an analysis module for kLambda events after being processed through the ECCE fun4all simulation. It is derived from previous modules created by Wenliang (Bill) Li and Stephen Kay.

### Set-up instructions are adapted from Stephen kay's ECCE_DEMP_Analysis repository.

### Within the VM (Virtual Machine):

cd ~/

source setup.sh

mkdir work

cd work

git clone https://github.com/ECCE-EIC/macros

git clone https://github.com/ECCE-EIC/tutorials

mkdir ECCE_kLambda

cd ECCE_kLambda

CreateSubsysRecoModule.pl ECCE_DEMP

cp ~/work/tutorials/CaloAna/src/*.[sa]* .

cp ~/work/ECCE_kLambda_Analysis/ECCE_kLambda_Ana/configure.ac ./

cp ~/work/ECCE_kLambda_Analysis/ECCE_kLambda_Ana/Makefile.am ./

chmod +x autogen.sh

mkdir build

cd build

../autogen.sh --prefix=$MYINSTALL

make install

### Now - copy in the latest versions of the actual script and rebuild the module

cp /home/fun4all/work/ECCE_kLambda_Analysis/ECCE_kLambda_Ana/ECCE_kLambda.* ./ && make install

### Do this whenever you want to actually re-build the plugin!

cd /home/fun4all/work/macros/detectors/EICDetector

cp ~/work/ECCE_kLambda_Analysis/main_macro/G4_User.C .

cp ~/work/ECCE_kLambda_Analysis/main_macro/Fun4all_reana.C .
