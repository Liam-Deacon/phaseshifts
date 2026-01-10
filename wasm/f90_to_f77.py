#!/usr/bin/env python3
"""
Convert Fortran 90/95 syntax to FORTRAN 77 for f2c compatibility.

This script converts F90 features that f2c doesn't understand:
- character(len=N) -> character*N
- intent(in/out/inout) declarations -> removed
- do while(...) -> labeled goto loops
- end do -> continue with labels
- trim() -> direct variable use (spaces are padded in F77)
- len_trim() -> custom function
- allocatable/allocate/deallocate -> static arrays (where possible)
- ! comments -> C comments
- end subroutine/function name -> end

Author: Generated for phaseshifts WASM build
"""

import re
import sys
from typing import List, Set, Tuple, Optional

# Constant for F90 comment prefix to avoid duplication
F90_COMMENT_PREFIX = "C F90: "


def scan_existing_labels(content: str) -> Set[int]:
    """Scan the source file for existing numeric labels to avoid collisions."""
    labels = set()
    for line in content.split("\n"):
        # Match labels at start of line (columns 1-5 in fixed-form Fortran)
        match = re.match(r"^\s*(\d+)\s+", line)
        if match:
            labels.add(int(match.group(1)))
        # Also check for goto/go to targets
        for goto_match in re.finditer(r"\b(?:goto|go\s+to)\s+(\d+)", line, re.IGNORECASE):
            labels.add(int(goto_match.group(1)))
        # Check for do loop labels: do 100 i=1,n
        do_match = re.match(r"\s*do\s+(\d+)\s+\w+\s*=", line, re.IGNORECASE)
        if do_match:
            labels.add(int(do_match.group(1)))
    return labels


def find_safe_label_start(existing_labels: Set[int], block_size: int = 1000) -> int:
    """Find a safe starting point for new labels that won't collide with existing ones."""
    if not existing_labels:
        return 9000

    max_label = max(existing_labels)
    # Start at next thousand boundary above max existing label
    start = ((max_label // block_size) + 1) * block_size
    return max(start, 9000)


def check_lentrim_exists(content: str) -> bool:
    """Check if a lentrim function definition already exists."""
    # Look for function definition
    pattern = r"^\s*(integer\s+)?function\s+lentrim\s*\("
    return bool(re.search(pattern, content, re.IGNORECASE | re.MULTILINE))


class LabelAllocator:
    """Manages label allocation to avoid collisions."""

    def __init__(self, existing_labels: Set[int]):
        self.existing_labels = existing_labels
        self.start = find_safe_label_start(existing_labels)
        self.counter = self.start
        self.allocated = set()

    def get_label(self) -> int:
        """Allocate a new unique label."""
        self.counter += 1
        while self.counter in self.existing_labels or self.counter in self.allocated:
            self.counter += 1
        self.allocated.add(self.counter)
        return self.counter


class LoopInfo:
    """Stores information about a loop for proper conversion."""

    def __init__(self, loop_type: str, start_label: int, exit_label: int, indent: str = "", condition: str = ""):
        self.loop_type = loop_type  # "while_true", "while_cond", "do_loop"
        self.start_label = start_label
        self.exit_label = exit_label
        self.indent = indent
        self.condition = condition


def convert_comments(line: str) -> Tuple[str, bool]:
    """Convert F90 comments to F77 style. Returns (converted_line, is_full_comment)."""
    # Full line comment
    if line.strip().startswith("!"):
        comment_text = line.strip()[1:]
        return "C" + comment_text, True

    # Inline comment - remove it
    if "!" in line:
        in_string = False
        quote_char = None
        for j, c in enumerate(line):
            if c in ('"', "'") and not in_string:
                in_string = True
                quote_char = c
            elif c == quote_char and in_string:
                in_string = False
                quote_char = None
            elif c == "!" and not in_string:
                return line[:j].rstrip(), False

    return line, False


def convert_character_decl(line: str) -> str:
    """Convert character(len=N) declarations to F77 style."""
    # Handle: character(len=255), intent(in) :: varname
    pattern = r"character\s*\(\s*len\s*=\s*(\d+|\*)\s*\)\s*,?\s*(intent\s*\([^)]*\)\s*)?::\s*"
    match = re.search(pattern, line, re.IGNORECASE)
    if match:
        length = match.group(1)
        if length == "*":
            length = "(*)"
        rest = line[match.end() :].strip()
        indent = line[: len(line) - len(line.lstrip())]
        return f"{indent}character*{length} {rest}"

    # Handle: character(len=N) varname (without ::)
    pattern = r"character\s*\(\s*len\s*=\s*(\d+|\*)\s*\)\s+(\w+)"
    match = re.search(pattern, line, re.IGNORECASE)
    if match:
        length = match.group(1)
        if length == "*":
            length = "(*)"
        varname = match.group(2)
        rest = line[match.end() :]
        prefix = line[: match.start()]
        return f"{prefix}character*{length} {varname}{rest}"

    return line


def convert_intent_and_colons(line: str) -> str:
    """Remove intent declarations and convert :: syntax."""
    # Remove intent declarations
    line = re.sub(r",\s*intent\s*\(\s*\w+\s*\)", "", line, flags=re.IGNORECASE)
    line = re.sub(r"intent\s*\(\s*\w+\s*\)\s*::", "::", line, flags=re.IGNORECASE)

    # Convert ":: varlist" to just "varlist" for declarations
    if "::" in line:
        decl_types = ["integer", "real", "double precision", "character", "logical"]
        is_decl = any(line.strip().lower().startswith(t) for t in decl_types)
        if is_decl:
            line = re.sub(r"\s*::\s*", " ", line)

    return line


def convert_string_functions(line: str) -> str:
    """Convert trim() and len_trim() to F77 equivalents."""
    # trim(var) -> var
    line = re.sub(r"\btrim\s*\(\s*(\w+)\s*\)", r"\1", line, flags=re.IGNORECASE)
    # len_trim(var) -> lentrim(var)
    line = re.sub(r"\blen_trim\s*\(\s*(\w+)\s*\)", r"lentrim(\1)", line, flags=re.IGNORECASE)
    return line


def convert_end_statements(line: str) -> str:
    """Convert F90 end statements to F77 style."""
    line = re.sub(r"\bend\s+subroutine\s+\w*", "end", line, flags=re.IGNORECASE)
    line = re.sub(r"\bend\s+function\s+\w*", "end", line, flags=re.IGNORECASE)
    line = re.sub(r"\bend\s+program\s+\w*", "end", line, flags=re.IGNORECASE)
    return line


def handle_allocatable(line: str, original_line: str) -> Optional[List[str]]:
    """Handle allocatable declarations. Returns list of output lines or None."""
    if "allocatable" not in line.lower():
        return None

    indent = line[: len(line) - len(line.lstrip())]

    # rpower(nrmax,0:15) - used in abinitio and pseudo
    if "rpower" in line.lower():
        return [F90_COMMENT_PREFIX + original_line.strip(), f"{indent}double precision rpower(4000,0:15)"]

    # Generic allocatable - just comment out
    return [F90_COMMENT_PREFIX + original_line.strip()]


def handle_allocate_deallocate(line: str) -> Optional[str]:
    """Handle allocate/deallocate statements. Returns comment or None."""
    if re.match(r"\s*allocate\s*\(", line, re.IGNORECASE):
        return F90_COMMENT_PREFIX + line.strip()
    if re.match(r"\s*deallocate\s*\(", line, re.IGNORECASE):
        return F90_COMMENT_PREFIX + line.strip()
    if "allocated(" in line.lower():
        return F90_COMMENT_PREFIX + line.strip()
    return None


def process_do_while_true(line: str, allocator: LabelAllocator, do_stack: List[LoopInfo]) -> Optional[str]:
    """Process 'do while(.true.)' statements."""
    match = re.match(r"(\s*)do\s+while\s*\(\s*\.true\.\s*\)", line, re.IGNORECASE)
    if not match:
        return None

    indent = match.group(1)
    start_label = allocator.get_label()
    exit_label = allocator.get_label()
    do_stack.append(LoopInfo("while_true", start_label, exit_label, indent))
    return f"{start_label:5d} continue"


def process_do_while_cond(line: str, allocator: LabelAllocator, do_stack: List[LoopInfo]) -> Optional[str]:
    """Process 'do while(condition)' statements."""
    match = re.match(r"(\s*)do\s+while\s*\(\s*(.+)\s*\)", line, re.IGNORECASE)
    if not match:
        return None

    indent = match.group(1)
    condition = match.group(2)
    start_label = allocator.get_label()
    exit_label = allocator.get_label()
    do_stack.append(LoopInfo("while_cond", start_label, exit_label, indent, condition))
    return f"{start_label:5d} if (.not.({condition})) goto {exit_label}"


# Pre-compiled regex patterns to reduce complexity
DO_LOOP_PATTERN = re.compile(r"(\s*)(\d+\s+)?do\s+(\w+)\s*=\s*([^,]+),([^,]+)(,.+)?$", re.IGNORECASE)


def process_do_loop(line: str, allocator: LabelAllocator, do_stack: List[LoopInfo]) -> Optional[str]:
    """Process simple 'do i=start,end' statements."""
    match = DO_LOOP_PATTERN.match(line)
    if not match or "while" in line.lower():
        return None

    indent = match.group(1)
    existing_label = match.group(2)
    var = match.group(3)
    start = match.group(4).strip()
    end = match.group(5).strip()
    step = match.group(6) if match.group(6) else ""

    loop_label = allocator.get_label()
    exit_label = allocator.get_label()
    do_stack.append(LoopInfo("do_loop", loop_label, exit_label, indent))

    if existing_label:
        return f"{existing_label.strip():>5s} do {loop_label} {var}={start},{end}{step}"
    return f"{indent}do {loop_label} {var}={start},{end}{step}"


def process_end_do(line: str, do_stack: List[LoopInfo]) -> Optional[List[str]]:
    """Process 'end do' statements."""
    match = re.match(r"(\s*)(\d+\s+)?end\s*do(\s*!.*)?$", line, re.IGNORECASE)
    if not match:
        return None

    existing_label = match.group(2)
    output = []

    if do_stack:
        loop_info = do_stack.pop()
        if loop_info.loop_type in ("while_true", "while_cond"):
            output.append(f"      goto {loop_info.start_label}")
            output.append(f"{loop_info.exit_label:5d} continue")
        elif loop_info.loop_type == "do_loop":
            output.append(f"{loop_info.start_label:5d} continue")
            output.append(f"{loop_info.exit_label:5d} continue")
    elif existing_label:
        output.append(f"{existing_label.strip():>5s} continue")
    else:
        output.append("      continue")

    return output


def process_exit(line: str, do_stack: List[LoopInfo]) -> Optional[str]:
    """Process standalone 'exit' statements."""
    match = re.match(r"(\s*)exit(\s*!.*)?$", line, re.IGNORECASE)
    if not match:
        return None

    if do_stack:
        loop_info = do_stack[-1]
        return f"      goto {loop_info.exit_label}"
    return F90_COMMENT_PREFIX + "exit (no matching loop)"


def process_inline_exit(line: str, do_stack: List[LoopInfo]) -> str:
    """Process inline 'if (cond) exit' statements."""
    match = re.search(r"\bexit\b(\s*!.*)?$", line, re.IGNORECASE)
    if match and do_stack:
        loop_info = do_stack[-1]
        line = re.sub(r"\bexit\b(\s*!.*)?$", f"goto {loop_info.exit_label}", line, flags=re.IGNORECASE)
    return line


def get_lentrim_helper() -> str:
    """Return the lentrim helper function code."""
    return """
C-----------------------------------------------------------------------
C     LENTRIM - F77 replacement for F90 LEN_TRIM intrinsic
C     Returns the length of string excluding trailing blanks
      integer function lentrim(str)
      character*(*) str
      integer i, n
      n = len(str)
      do 100 i = n, 1, -1
        if (str(i:i) .ne. ' ') then
          lentrim = i
          return
        endif
  100 continue
      lentrim = 0
      return
      end
"""


def find_last_end_position(lines: List[str]) -> int:
    """Find the position of the last 'end' statement (for inserting helper before it)."""
    for i in range(len(lines) - 1, -1, -1):
        if re.match(r"^\s*end\s*$", lines[i], re.IGNORECASE):
            return i
    return len(lines)


def process_line(
    line: str, original_line: str, allocator: LabelAllocator, do_stack: List[LoopInfo]
) -> Tuple[Optional[List[str]], bool]:
    """
    Process a single line of F90 code.

    Returns:
        Tuple of (output_lines, done) where:
        - output_lines: List of converted lines, or None to continue with default handling
        - done: If True, skip further processing of this line
    """
    # Convert comments
    line, is_full_comment = convert_comments(line)
    if is_full_comment:
        return [line], True

    # Convert character declarations
    line = convert_character_decl(line)

    # Remove intent and convert :: syntax
    line = convert_intent_and_colons(line)

    # Convert string functions
    line = convert_string_functions(line)

    # Handle allocatable
    alloc_result = handle_allocatable(line, original_line)
    if alloc_result is not None:
        return alloc_result, True

    # Handle allocate/deallocate
    alloc_stmt = handle_allocate_deallocate(line)
    if alloc_stmt is not None:
        return [alloc_stmt], True

    # Process do while(.true.)
    result = process_do_while_true(line, allocator, do_stack)
    if result is not None:
        return [result], True

    # Process do while(condition)
    result = process_do_while_cond(line, allocator, do_stack)
    if result is not None:
        return [result], True

    # Process simple do loop
    result = process_do_loop(line, allocator, do_stack)
    if result is not None:
        return [result], True

    # Process end do
    end_do_result = process_end_do(line, do_stack)
    if end_do_result is not None:
        return end_do_result, True

    # Process standalone exit
    exit_result = process_exit(line, do_stack)
    if exit_result is not None:
        return [exit_result], True

    # Process inline exit
    line = process_inline_exit(line, do_stack)

    # Convert end statements
    line = convert_end_statements(line)

    return [line], True


def convert_f90_to_f77(content: str) -> str:
    """Convert F90 code to F77 compatible code."""
    lines = content.split("\n")
    output_lines: List[str] = []

    # Scan for existing labels to avoid collisions
    existing_labels = scan_existing_labels(content)
    allocator = LabelAllocator(existing_labels)

    # Check if lentrim already exists
    lentrim_exists = check_lentrim_exists(content)

    # Stack to track do loops for end do conversion
    do_stack: List[LoopInfo] = []

    for line in lines:
        # Skip empty lines
        if not line.strip():
            output_lines.append(line)
            continue

        result, _ = process_line(line, line, allocator, do_stack)
        if result is not None:
            output_lines.extend(result)

    # Add lentrim helper if needed and not already present
    result = "\n".join(output_lines)
    if "lentrim(" in result.lower() and not lentrim_exists:
        # Find the last 'end' and insert before it
        last_end_pos = find_last_end_position(output_lines)
        helper = get_lentrim_helper()
        output_lines.insert(last_end_pos, helper)
        result = "\n".join(output_lines)

    return result


def main():
    if len(sys.argv) < 2:
        print("Usage: f90_to_f77.py <input.f> [output.f]")
        print("Converts Fortran 90 code to FORTRAN 77 for f2c compatibility")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else None

    # Read input file with error handling
    try:
        with open(input_file, "r", encoding="utf-8") as f:
            content = f.read()
    except FileNotFoundError as e:
        print(f"Error: Input file not found: {input_file}")
        print(f"  {e}")
        sys.exit(1)
    except PermissionError as e:
        print(f"Error: Permission denied reading file: {input_file}")
        print(f"  {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error: Failed to read input file: {input_file}")
        print(f"  {e}")
        sys.exit(1)

    converted = convert_f90_to_f77(content)

    if output_file:
        # Write output file with error handling
        try:
            with open(output_file, "w", encoding="utf-8") as f:
                f.write(converted)
            print(f"Converted {input_file} -> {output_file}")
        except FileNotFoundError as e:
            print(f"Error: Output directory not found for: {output_file}")
            print(f"  {e}")
            sys.exit(1)
        except PermissionError as e:
            print(f"Error: Permission denied writing file: {output_file}")
            print(f"  {e}")
            sys.exit(1)
        except Exception as e:
            print(f"Error: Failed to write output file: {output_file}")
            print(f"  {e}")
            sys.exit(1)
    else:
        print(converted)


if __name__ == "__main__":
    main()
