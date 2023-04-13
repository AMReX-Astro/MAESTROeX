To run this flame test problem, you will first need to download the appropriate Pytorch library (libtorch) in `MAESTROeX/external` directory. We can choose to use either CPU or CUDA.

```shell
cd $MAESTROEX_HOME/external

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

To compile the problem with ML enabled, enter
`make -j`.

Then you can run the problem.

`mpiexec -n 1 ./Maestro2d.gnu.MPI.ML.ex inputs_2d_smallscale_ml`


<!---
Torchscript code based on https://pytorch.org/tutorials/advanced/cpp_export.html 
-->
