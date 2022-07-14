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

To compile:
```shell
# with ML model
make -j

# without ML model (original Strang)
make -j USE_ML=FALSE
```


## Training a Pytorch model

1. Compile the code with `USE_ML=FALSE`.
2. Run executable with input file `input_2d_smallscale_traindata`.
3. Create a new directory named `/data`.
4. Move all output directories starting with prefix `flame_` and `react_` into the previously created `/data/` directory.
5. Navigate to `/python_code` directory.
6. Run `python3 test_igsimple.py` to start training a model for 10 epochs.
