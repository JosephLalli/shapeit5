# C++ Style Guide

This project follows a specific coding style that mimics the original SHAPEIT5 conventions. All new code must adhere to these guidelines to ensure consistency.

## 1. Formatting
- **Indentation:** Use **tabs** for indentation. Assume a tab width of 4 spaces.
- **Line Length:** No hard limit, but prefer keeping lines under 120 characters where possible.
- **Braces:** Use K&R style (opening brace on the same line as the control statement or function definition).
    ```cpp
    if (condition) {
        // code
    } else {
        // code
    }
    ```
- **Spacing:**
    - Put a space after control flow keywords (`if`, `for`, `while`, `switch`).
    - Put a space before the opening parenthesis of a control flow statement.
    - Put a space between the parenthesis and the condition/expression.
    - **Example:** `if ( condition ) { ... }` (Note: The codebase uses `if (condition)` sometimes, but `if ( condition )` is also common. Follow the local file's style. For now, `if (condition)` is safer).
    - **Correction:** Looking at `phaser_algorithm.cpp`: `for(;;) {`, `if (id_job <= S->G.n_ind)`. It seems spaces *inside* parens are NOT strictly enforced, but space *after* keyword is inconsistent (`for(;;)` vs `if (...)`).
    - **Rule:** Use space after keyword: `if (...)`, `for (...)`. No space inside parens unless complex.

## 2. Naming Conventions
- **Variables:** `snake_case` (e.g., `id_worker`, `n_underflow_recovered`).
- **Functions:** `camelCase` or `snake_case` (mixed in codebase, but `phaseWindow`, `pushIBD2` suggest camelCase for methods).
    - **Rule:** Use `camelCase` for class methods (`phaseWindow`, `backward`).
    - **Rule:** Use `snake_case` for free functions if any.
- **Classes/Structs:** `snake_case` (e.g., `haplotype_segment_single`, `phaser`) or `CamelCase` (e.g., `SuperSite`).
    - **Rule:** Follow the surrounding code. Most core classes like `phaser` and `haplotype_segment_single` are snake_case.
- **Macros:** `UPPER_CASE` (e.g., `DIV2`, `MOD2`).

## 3. File Structure
- **Headers:** `.h` files. Include guards `#ifndef _HEADER_H` ... `#endif`.
- **Sources:** `.cpp` files.
- **Includes:** Group imports:
    1. Local Project Headers
    2. STL Headers
    3. Boost Headers
    4. HTSlib Headers

## 4. Modern C++ Features
- Use `auto` where type is obvious or complex (iterators).
- Use `static_cast` instead of C-style casts.
- Use `std::vector` and `std::string`.
- Use AVX2 intrinsics (`__m256`) for performance-critical sections.

## 5. Comments
- **File Header:** Every file must start with the standard license header.
- **Inline:** Use `//` for brief explanations.
- **Documentation:** Focus on *why* complex logic exists (especially math/SIMD).

## 6. Error Handling
- Use `vrb.error("message")` to report fatal errors and exit.
- Use `vrb.warning("message")` for non-fatal issues.
