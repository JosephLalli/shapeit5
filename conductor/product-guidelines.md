# Product Guidelines

These guidelines define the development standards and philosophy for the SHAPEIT5 native multiallelic extension.

## Codebase Philosophy
- **Extension, Not Replacement:** The primary goal is to extend the original work of the SHAPEIT5 creators. New features should feel like a natural evolution of the existing codebase, preserving the core algorithmic "signature."
- **Strict Style Adherence:** All new code, including 8-bit genotype handling and updated readers, must strictly mimic the existing naming conventions, indentation patterns, and file structure of the original SHAPEIT5 project.
- **Transparency:** While processing should feel seamless, log messages should clearly indicate when the native multiallelic path is being utilized for a site to aid in debugging and validation.

## Documentation Standards
- **Technical Accuracy & Rationale:** Prioritize documenting the "why" behind algorithmic choices, especially regarding SIMD (AVX2/FMA) optimizations and HMM transition logic.
- **Implementation Clarity:** Provide detailed internal comments for all data structure modifications, particularly the transition from bit-packed genotypes to the 8-bit (1-byte) representation, so the original authors can easily reason about the changes.

## Quality & Performance
- **Zero Regression for Biallelic Sites:** This is the paramount requirement. Every change must be verified through rigorous parity testing, ensuring that output for biallelic-only datasets is identical to the legacy baseline.
- **Efficiency Targets:** The computational overhead for native multiallelic support should be kept within a 10-20% margin compared to the biallelic-only version.
