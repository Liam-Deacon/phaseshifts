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


def convert_f90_to_f77(content):
    """Convert F90 code to F77 compatible code."""
    lines = content.split("\n")
    output_lines = []

    # Counter for generating unique labels
    label_counter = [9000]
    do_stack = []  # Stack to track do loops for end do conversion

    def get_label():
        label_counter[0] += 1
        return label_counter[0]

    for i, line in enumerate(lines):
        original_line = line

        # Skip empty lines
        if not line.strip():
            output_lines.append(line)
            continue

        # Convert ! comments to C comments (for lines starting with !)
        # But preserve inline comments by converting them
        if line.strip().startswith("!"):
            # Full line comment - convert to C style
            indent = len(line) - len(line.lstrip())
            comment_text = line.strip()[1:]  # Remove the !
            output_lines.append("C" + comment_text)
            continue

        # Handle inline ! comments - convert to C style or remove
        if "!" in line and not line.strip().startswith("!"):
            # Check if ! is inside a string
            in_string = False
            quote_char = None
            comment_pos = -1
            for j, c in enumerate(line):
                if c in ('"', "'") and not in_string:
                    in_string = True
                    quote_char = c
                elif c == quote_char and in_string:
                    in_string = False
                    quote_char = None
                elif c == "!" and not in_string:
                    comment_pos = j
                    break
            if comment_pos > 0:
                # Remove inline comment for F77 compatibility
                line = line[:comment_pos].rstrip()

        # Convert character(len=N) to character*N
        # Handle: character(len=255), intent(in) :: varname
        pattern = r"character\s*\(\s*len\s*=\s*(\d+|\*)\s*\)\s*,?\s*(intent\s*\([^)]*\)\s*)?::\s*"
        match = re.search(pattern, line, re.IGNORECASE)
        if match:
            length = match.group(1)
            # For len=*, use (*) in F77
            if length == "*":
                length = "(*)"
            # Get variable names after ::
            rest = line[match.end() :]
            vars_part = rest.strip()
            # Reconstruct as F77 style
            indent = line[: len(line) - len(line.lstrip())]
            line = f"{indent}character*{length} {vars_part}"

        # Convert character(len=N) varname (without ::)
        pattern = r"character\s*\(\s*len\s*=\s*(\d+|\*)\s*\)\s+(\w+)"
        match = re.search(pattern, line, re.IGNORECASE)
        if match:
            length = match.group(1)
            # For len=*, use (*) in F77
            if length == "*":
                length = "(*)"
            varname = match.group(2)
            rest = line[match.end() :]
            prefix = line[: match.start()]
            line = f"{prefix}character*{length} {varname}{rest}"

        # Remove standalone intent declarations
        line = re.sub(r",\s*intent\s*\(\s*\w+\s*\)", "", line, flags=re.IGNORECASE)
        line = re.sub(r"intent\s*\(\s*\w+\s*\)\s*::", "::", line, flags=re.IGNORECASE)

        # Convert ":: varlist" to just "varlist" for declarations
        # But be careful - only do this after type declarations
        if "::" in line:
            # Check if this is a declaration line
            decl_types = ["integer", "real", "double precision", "character", "logical"]
            is_decl = any(line.strip().lower().startswith(t) for t in decl_types)
            if is_decl:
                line = re.sub(r"\s*::\s*", " ", line)

        # Convert trim(var) to var (F77 doesn't have trim, strings are space-padded)
        line = re.sub(r"\btrim\s*\(\s*(\w+)\s*\)", r"\1", line, flags=re.IGNORECASE)

        # Convert len_trim(var) to a custom function call or index function
        # For simplicity, we'll use INDEX approach or just use LEN
        # Actually, let's define a helper function and use it
        # For now, replace with a simpler check - len_trim returns length without trailing spaces
        # We'll replace len_trim(x) with LENTRIM(x) and define the function
        line = re.sub(r"\blen_trim\s*\(\s*(\w+)\s*\)", r"lentrim(\1)", line, flags=re.IGNORECASE)

        # Convert allocatable declarations - replace with static arrays
        # Example: double precision, allocatable :: rpower(:,:)
        # We need to make these static with appropriate sizes
        if "allocatable" in line.lower():
            # Common patterns we know about from libphsh.f
            indent = line[: len(line) - len(line.lstrip())]

            # rpower(nrmax,0:15) - used in abinitio and pseudo
            if "rpower" in line.lower():
                output_lines.append("C F90: " + original_line.strip())
                output_lines.append(f"{indent}double precision rpower(4000,0:15)")
                continue

            # Generic allocatable - just comment out
            output_lines.append("C F90: " + original_line.strip())
            continue

        # Convert allocate statements - comment them out
        if re.match(r"\s*allocate\s*\(", line, re.IGNORECASE):
            output_lines.append("C F90: " + line.strip())
            continue

        # Convert deallocate statements - comment them out
        if re.match(r"\s*deallocate\s*\(", line, re.IGNORECASE):
            output_lines.append("C F90: " + line.strip())
            continue

        # Convert "if (allocated(...))" - comment out
        if "allocated(" in line.lower():
            output_lines.append("C F90: " + line.strip())
            continue

        # Convert "do while(.true.)" to labeled infinite loop
        match = re.match(r"(\s*)do\s+while\s*\(\s*\.true\.\s*\)", line, re.IGNORECASE)
        if match:
            indent = match.group(1)
            label = get_label()
            do_stack.append(("while_true", label, indent))
            output_lines.append(f"{label:5d} continue")
            continue

        # Convert "do while(condition)" to labeled loop with if-goto
        match = re.match(r"(\s*)do\s+while\s*\(\s*(.+)\s*\)", line, re.IGNORECASE)
        if match:
            indent = match.group(1)
            condition = match.group(2)
            start_label = get_label()
            end_label = get_label()
            do_stack.append(("while_cond", start_label, end_label, indent, condition))
            output_lines.append(f"{start_label:5d} if (.not.({condition})) goto {end_label}")
            continue

        # Convert simple "do i=start,end" - track for end do
        # Also handles labeled do like "8 do IX=1,NGRID"
        match = re.match(r"(\s*)(\d+\s+)?do\s+(\w+)\s*=\s*(.+?)\s*,\s*(.+?)(\s*,\s*.+)?$", line, re.IGNORECASE)
        if match and "while" not in line.lower():
            indent = match.group(1)
            existing_label = match.group(2)
            var = match.group(3)
            start = match.group(4)
            end = match.group(5)
            step = match.group(6) if match.group(6) else ""
            label = get_label()
            do_stack.append(("do_loop", label, indent))
            if existing_label:
                # Preserve original label at start, use our label for continue
                output_lines.append(f"{existing_label.strip():>5s} do {label} {var}={start},{end}{step}")
            else:
                output_lines.append(f"{indent}do {label} {var}={start},{end}{step}")
            continue

        # Convert "end do" to "continue" with label
        # Also handle labeled end do like "2990   end do"
        match = re.match(r"(\s*)(\d+\s+)?end\s*do(\s*!.*)?$", line, re.IGNORECASE)
        if match:
            existing_label = match.group(2)
            if do_stack:
                loop_info = do_stack.pop()
                if loop_info[0] == "while_true":
                    # Jump back to start, add exit label after
                    output_lines.append(f"      goto {loop_info[1]}")
                    # Add exit label for any exit statements in this loop
                    exit_label = loop_info[1] + 5000  # offset for exit labels
                    output_lines.append(f"{exit_label:5d} continue")
                elif loop_info[0] == "while_cond":
                    # Jump back to condition check, then add end label
                    output_lines.append(f"      goto {loop_info[1]}")
                    output_lines.append(f"{loop_info[2]:5d} continue")
                elif loop_info[0] == "do_loop":
                    # Add continue with label, then exit label
                    output_lines.append(f"{loop_info[1]:5d} continue")
                    # Add exit label after the loop for exit statements
                    exit_label = loop_info[1] + 5000
                    output_lines.append(f"{exit_label:5d} continue")
            elif existing_label:
                # Has a label but no matching do in stack - keep the label
                output_lines.append(f"{existing_label.strip():>5s} continue")
            else:
                # No matching do - just output continue
                output_lines.append("      continue")
            continue

        # Convert "exit" to goto the end of the current loop
        # Handle both standalone "exit" and inline "if (cond) exit"
        # Standalone exit
        match = re.match(r"(\s*)exit(\s*!.*)?$", line, re.IGNORECASE)
        if match:
            if do_stack:
                # Get the exit label for current loop
                loop_info = do_stack[-1]  # peek, don't pop
                if loop_info[0] == "while_true":
                    exit_label = loop_info[1] + 5000
                    output_lines.append(f"      goto {exit_label}")
                elif loop_info[0] == "while_cond":
                    exit_label = loop_info[2]  # end label
                    output_lines.append(f"      goto {exit_label}")
                elif loop_info[0] == "do_loop":
                    # For do loops, we need to add an exit label after the continue
                    # Use the loop label + 5000 as convention
                    exit_label = loop_info[1] + 5000
                    output_lines.append(f"      goto {exit_label}")
            else:
                output_lines.append("C F90: exit (no matching loop)")
            continue

        # Inline exit: "if (condition) exit" -> "if (condition) goto label"
        match = re.search(r"\bexit\b(\s*!.*)?$", line, re.IGNORECASE)
        if match and do_stack:
            loop_info = do_stack[-1]  # peek
            if loop_info[0] == "while_true":
                exit_label = loop_info[1] + 5000
            elif loop_info[0] == "while_cond":
                exit_label = loop_info[2]
            elif loop_info[0] == "do_loop":
                exit_label = loop_info[1] + 5000
            else:
                exit_label = 99999  # fallback
            # Replace exit with goto
            line = re.sub(r"\bexit\b(\s*!.*)?$", f"goto {exit_label}", line, flags=re.IGNORECASE)

        # Convert "end subroutine name" to just "end"
        line = re.sub(r"\bend\s+subroutine\s+\w*", "end", line, flags=re.IGNORECASE)
        line = re.sub(r"\bend\s+function\s+\w*", "end", line, flags=re.IGNORECASE)
        line = re.sub(r"\bend\s+program\s+\w*", "end", line, flags=re.IGNORECASE)

        # Convert "exit" to goto (need context, skip for now - rare usage)
        # Convert "cycle" to goto (need context, skip for now)

        # Convert double precision function to older style if needed
        # (Usually OK as-is for f2c)

        output_lines.append(line)

    # Add helper function for lentrim at the end
    helper_functions = """
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

    result = "\n".join(output_lines)

    # Only add helper if lentrim is used
    if "lentrim(" in result.lower():
        result = result + helper_functions

    return result


def main():
    if len(sys.argv) < 2:
        print("Usage: f90_to_f77.py <input.f> [output.f]")
        print("Converts Fortran 90 code to FORTRAN 77 for f2c compatibility")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else None

    with open(input_file, "r") as f:
        content = f.read()

    converted = convert_f90_to_f77(content)

    if output_file:
        with open(output_file, "w") as f:
            f.write(converted)
        print(f"Converted {input_file} -> {output_file}")
    else:
        print(converted)


if __name__ == "__main__":
    main()
