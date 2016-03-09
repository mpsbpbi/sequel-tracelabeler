# TROL installation guide

Here is a walkthrough for installing this package on pa-dev01.lab.nanofluidics.com

```bash
module load virtualenv VE
module load python/2.7.9
module unload composer_xe/2015.2.164
module load pacbio-thirdparty

sudo yum install python-matplotlib
sudo yum install libfreetype6-devel
sudo yum install freetype-devel

virtualenv VE
source VE/bin/activate


pip install cython
pip install numpy
pip install h5py
pip install pysam
pip install pbcore
pip install matplotlib
```
