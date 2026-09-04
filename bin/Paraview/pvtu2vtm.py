#!/usr/bin/env python3
"""
Replace each "dataXXX.pvtu" manifest with a two-block "dataXXX.vtm"
multiblock file that points directly at the existing "dataXXX_b0.vtu"
(fluid) and "dataXXX_b1.vtu" (solid) pieces.

Usage:
    python3 pvtu_to_vtm.py [directory]

If no directory is given, the current directory is used.
"""
import re
import sys
from pathlib import Path

VTM_TEMPLATE = """<VTKFile type="vtkMultiBlockDataSet" version="1.0" byte_order="LittleEndian">
  <vtkMultiBlockDataSet>
    <DataSet index="0" name="fluid" file="{b0}"/>
    <DataSet index="1" name="solid" file="{b1}"/>
  </vtkMultiBlockDataSet>
</VTKFile>
"""

PVTU_PATTERN = re.compile(r"^data(\d+)\.pvtu$")


def convert(directory: Path, delete_pvtu: bool = False) -> None:
    pvtu_files = sorted(directory.glob("data*.pvtu"))

    if not pvtu_files:
        print(f"No dataXXX.pvtu files found in {directory}")
        return

    for pvtu_path in pvtu_files:
        match = PVTU_PATTERN.match(pvtu_path.name)
        if not match:
            continue

        xxx = match.group(1)
        b0_name = f"data{xxx}_b0.vtu"
        b1_name = f"data{xxx}_b1.vtu"
        b0_path = directory / b0_name
        b1_path = directory / b1_name

        if not b0_path.exists() or not b1_path.exists():
            missing = [p.name for p in (b0_path, b1_path) if not p.exists()]
            print(f"Skipping {pvtu_path.name}: missing block file(s) {missing}")
            continue

        vtm_path = directory / f"data{xxx}.vtm"
        vtm_path.write_text(VTM_TEMPLATE.format(b0=b0_name, b1=b1_name))
        print(f"Wrote {vtm_path.name}  (blocks: {b0_name}, {b1_name})")

        if delete_pvtu:
            pvtu_path.unlink()
            print(f"Removed {pvtu_path.name}")


if __name__ == "__main__":
    target_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")
    convert(target_dir, delete_pvtu=False)