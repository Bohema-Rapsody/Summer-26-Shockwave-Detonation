from pathlib import Path

# ============================================================
# USER SETTINGS
# ============================================================

input_file = "N2_Shock/dump/Ti_dump_300.xyz"
output_file = "N2_Shock/dump/Ti_dump_300_cleaned.xyz"

# Frame numbers to REMOVE.
# Frame numbering starts at 0.

frames_to_remove = set(range(48,74))


# ============================================================
# XYZ TRAJECTORY PROCESSING
# ============================================================

def remove_frames(input_file, output_file, frames_to_remove):

    input_path = Path(input_file)
    output_path = Path(output_file)

    frame_number = 0
    frames_read = 0
    frames_written = 0

    with input_path.open("r") as infile, output_path.open("w") as outfile:

        while True:

            # ------------------------------------------------
            # Read number of atoms
            # ------------------------------------------------
            line = infile.readline()

            if not line:
                break  # End of file

            line = line.strip()

            if not line:
                continue

            try:
                n_atoms = int(line)
            except ValueError:
                raise ValueError(
                    f"Could not read atom count at frame {frame_number}: "
                    f"{line!r}"
                )

            # ------------------------------------------------
            # Read comment/header line
            # ------------------------------------------------
            comment = infile.readline()

            if not comment:
                raise ValueError(
                    f"Unexpected end of file while reading "
                    f"header of frame {frame_number}"
                )

            # ------------------------------------------------
            # Read all atoms in this frame
            # ------------------------------------------------
            atoms = []

            for _ in range(n_atoms):

                atom_line = infile.readline()

                if not atom_line:
                    raise ValueError(
                        f"Unexpected end of file while reading "
                        f"frame {frame_number}"
                    )

                atoms.append(atom_line)

            frames_read += 1

            # ------------------------------------------------
            # Write frame unless it is in the removal list
            # ------------------------------------------------
            if frame_number not in frames_to_remove:

                outfile.write(f"{n_atoms}\n")
                outfile.write(comment)
                outfile.writelines(atoms)

                frames_written += 1

            else:
                print(f"Removing frame {frame_number}")

            frame_number += 1

    # ========================================================
    # SUMMARY
    # ========================================================

    print()
    print("Processing complete.")
    print("-------------------")
    print(f"Input file:       {input_path}")
    print(f"Output file:      {output_path}")
    print(f"Frames read:      {frames_read}")
    print(f"Frames removed:   {frames_read - frames_written}")
    print(f"Frames written:   {frames_written}")


# ============================================================
# RUN
# ============================================================

remove_frames(
    input_file,
    output_file,
    frames_to_remove
)