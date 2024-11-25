# JacobiSVD

```bash
.
├── BlockParallel # Block Parallel Jacobi
│   ├── BlockParallel.cpp # Include Jacobi Rotation and One-sided Jacobi SVD
│   └── BlockParallel.h 
├── CMakeLists.txt # CMake
├── common.h # Set some parameters
├── demo.cpp # demp function
├── main.cpp # main function
├── Matrix # Matrix with simple operation
│   ├── Matrix.cpp 
│   └── Matrix.h 
├── Parallel # Direct Parallel Jacobi
│   ├── Parallel.cpp # Include Jacobi Rotation and One-sided Jacobi SVD
│   └── Parallel.h 
├── Serial # Serial Jacobi
│   ├── Serial.cpp # Include Jacobi Rotation and One-sided Jacobi SVD
│   └── Serial.h 
└── build # where you should using cmake
```