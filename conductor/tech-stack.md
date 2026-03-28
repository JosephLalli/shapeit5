# Technology Stack

## Programming Language
- **C++17:** Core language for the project, chosen for performance and modern features.

## Libraries & Frameworks
- **HTSlib (v1.22):** The standard library for handling high-throughput sequencing data formats (SAM/BAM/CRAM/VCF/BCF).
- **Boost Libraries:**
    - `boost::iostreams`: For efficient input/output operations.
    - `boost::program_options`: For parsing command-line arguments.
    - `boost::serialization`: For serializing data structures.

## Performance Optimization
- **SIMD Intrinsics:** Heavy usage of Intel AVX2 and FMA (Fused Multiply-Add) intrinsics for accelerating HMM computations.

## Build System
- **Makefiles:** Hierarchical makefile system managing the build process across multiple sub-projects.

## Target Platform
- **Linux/Unix:** Optimized for high-performance computing environments.
