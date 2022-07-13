To run MAESTROeX with a Pytorch ML model, you will first need to download the appropriate Pytorch library (libtorch) in the current directory. You can either choose to use the CPU host or GPU CUDA version of libtorch.

```shell
# navigate to the external folder in MAESTROeX main directory
cd <MAESTROeX_dir>/external

# CPU only
wget https://download.pytorch.org/libtorch/cpu/libtorch-cxx11-abi-shared-with-deps-1.9.0%2Bcpu.zip
unzip libtorch-cxx11-abi-shared-with-deps-1.9.0+cpu.zip

# CUDA 11.1
wget https://download.pytorch.org/libtorch/cu111/libtorch-cxx11-abi-shared-with-deps-1.9.0%2Bcu111.zip
unzip libtorch-cxx11-abi-shared-with-deps-1.9.0+cu111.zip

# CUDA 10.2
wget https://download.pytorch.org/libtorch/cu102/libtorch-cxx11-abi-shared-with-deps-1.9.0%2Bcu102.zip
unzip libtorch-cxx11-abi-shared-with-deps-1.9.0+cu102.zip
```
