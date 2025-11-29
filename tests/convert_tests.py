#!/usr/bin/env python3
"""
Automated test conversion script for SHAPEIT5 test suite.
Converts tests to use standardized test_reporting.h format.
"""

import re
import sys
from pathlib import Path

def extract_test_name(filepath):
    """Extract test name from filename."""
    return filepath.stem

def has_test_reporting(content):
    """Check if file already includes test_reporting.h."""
    return '#include "test_reporting.h"' in content or '#include "test_reporting.h' in content

def add_test_reporting_include(content):
    """Add test_reporting.h include after first #include block."""
    lines = content.split('\n')

    # Find first #include line
    first_include = -1
    for i, line in enumerate(lines):
        if line.strip().startswith('#include'):
            first_include = i
            break

    if first_include == -1:
        return content

    # Find end of include block
    last_include = first_include
    for i in range(first_include + 1, len(lines)):
        if lines[i].strip().startswith('#include'):
            last_include = i
        elif lines[i].strip() and not lines[i].strip().startswith('//'):
            break

    # Insert after includes but before #define
    insert_pos = last_include + 1
    # Skip blank lines
    while insert_pos < len(lines) and not lines[insert_pos].strip():
        insert_pos += 1

    # If next line is #define, insert before it
    if insert_pos < len(lines) and lines[insert_pos].strip().startswith('#define'):
        lines.insert(insert_pos, '')
        lines.insert(insert_pos, '#include "test_reporting.h"')
    else:
        lines.insert(insert_pos, '#include "test_reporting.h"')
        lines.insert(insert_pos, '')

    return '\n'.join(lines)

def add_test_init(content, test_name):
    """Add TEST_INIT at start of main()."""
    # Find int main() {
    pattern = r'(int\s+main\s*\([^)]*\)\s*\{)'

    def replacement(match):
        return match.group(1) + '\n    TEST_INIT("' + test_name + '");'

    return re.sub(pattern, replacement, content, count=1)

def add_test_summary(content):
    """Add TEST_SUMMARY() before final return 0."""
    # Find last return 0; in main
    lines = content.split('\n')

    # Find last return 0;
    last_return = -1
    for i in range(len(lines) - 1, -1, -1):
        if re.search(r'return\s+0\s*;', lines[i]):
            last_return = i
            break

    if last_return == -1:
        return content

    # Add TEST_SUMMARY before return
    indent = len(lines[last_return]) - len(lines[last_return].lstrip())
    lines.insert(last_return, ' ' * indent + 'TEST_SUMMARY();')

    return '\n'.join(lines)

def convert_simple_outputs(content):
    """Convert simple cout << "passed" patterns to TEST_PASS."""

    # Pattern: std::cout << "Test name: OK" or similar
    # This is conservative - only handles very simple cases
    content = re.sub(
        r'std::cout\s*<<\s*"([^"]+):\s*OK"\s*<<\s*std::endl\s*;',
        r'TEST_PASS("\1");  // was: OK',
        content
    )

    content = re.sub(
        r'std::cout\s*<<\s*"([^"]+):\s*PASS"\s*<<\s*std::endl\s*;',
        r'TEST_PASS("\1");  // was: PASS',
        content
    )

    return content

def convert_test_file(filepath, dry_run=False):
    """Convert a single test file to new format."""
    print(f"\nProcessing: {filepath.name}")

    content = filepath.read_text()

    # Skip if already converted
    if has_test_reporting(content):
        print(f"  ⏭  Already has test_reporting.h - skipping")
        return False

    # Skip test_toolbox (special case)
    if 'test_toolbox' in filepath.name:
        print(f"  ⏭  Toolbox file - skipping")
        return False

    test_name = extract_test_name(filepath)

    original = content

    # Step 1: Add include
    content = add_test_reporting_include(content)

    # Step 2: Add TEST_INIT
    content = add_test_init(content, test_name)

    # Step 3: Add TEST_SUMMARY
    content = add_test_summary(content)

    # Step 4: Convert simple patterns (conservative)
    content = convert_simple_outputs(content)

    if content == original:
        print(f"  ⚠  No changes made - may need manual conversion")
        return False

    if dry_run:
        print(f"  ✓  Would update (dry-run)")
        # Show diff summary
        orig_lines = original.split('\n')
        new_lines = content.split('\n')
        print(f"     Lines: {len(orig_lines)} → {len(new_lines)} (Δ{len(new_lines)-len(orig_lines):+d})")
    else:
        filepath.write_text(content)
        print(f"  ✓  Updated")

    return True

def main():
    tests_dir = Path('/mnt/d/shapeit5/tests/src')
    test_files = sorted(tests_dir.glob('test_*.cpp'))

    # Exclude already converted
    exclude = {'test_supersite_builder.cpp', 'test_supersite_expansion_epochs.cpp', 'test_toolbox.cpp'}
    test_files = [f for f in test_files if f.name not in exclude]

    print(f"Found {len(test_files)} test files to convert")
    print(f"(Excluding {len(exclude)} already converted/special files)")

    dry_run = '--dry-run' in sys.argv
    if dry_run:
        print("\n=== DRY RUN MODE ===\n")

    converted = 0
    skipped = 0

    for test_file in test_files:
        if convert_test_file(test_file, dry_run=dry_run):
            converted += 1
        else:
            skipped += 1

    print(f"\n{'='*60}")
    print(f"Summary:")
    print(f"  Converted: {converted}")
    print(f"  Skipped: {skipped}")
    print(f"  Total: {len(test_files)}")

    if dry_run:
        print(f"\nRun without --dry-run to apply changes")

if __name__ == '__main__':
    main()
