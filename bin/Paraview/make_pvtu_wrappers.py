#!/usr/bin/env python3
"""
Write thin, single-piece .pvtu "wrapper" files for each dataXXX_bN.vtu file,
named so the timestep number is trailing (immediately before the extension)
-- the naming FieldView's transient-sequence dialog expects -- without
touching, renaming, or duplicating the actual (large) .vtu data files.

For a file like:
    data108_b0.vtu
this writes a wrapper:
    data_b0.108.pvtu
which points back at data108_b0.vtu via a single <Piece Source="..."/>.

The wrapper's PPointData/PCellData/PPoints schema is read directly from the
real .vtu file's own PointData DataArray tags, so it stays correct even if
the number/names of variables differ between blocks (e.g. liquid vs solid),
or changes in the future.

Usage:
    python3 make_pvtu_wrappers.py [directory]

If no directory is given, the current directory is used.
"""
import re
import sys
from pathlib import Path

VTU_PATTERN = re.compile(r"^(?P<label>.*?)(?P<num>\d+)_b(?P<block>\d+)\.vtu$")
DATAARRAY_PATTERN = re.compile(
    r'<DataArray\s+type="(?P<type>[^"]+)"\s+Name="(?P<name>[^"]+)"'
    r'(?:\s+NumberOfComponents="(?P<ncomp>\d+)")?'
)


def extract_point_data_schema(vtu_path: Path):
    """Read the real .vtu file's <PointData> block and return a list of
    (type, name, ncomp) tuples describing each DataArray it contains."""
    text = vtu_path.read_text()

    pd_match = re.search(r"<PointData[^>]*>(.*?)</PointData>", text, re.S)
    if not pd_match:
        return []

    schema = []
    for m in DATAARRAY_PATTERN.finditer(pd_match.group(1)):
        ncomp = m.group("ncomp") or "1"
        schema.append((m.group("type"), m.group("name"), ncomp))
    return schema


def write_pvtu_wrapper(wrapper_path: Path, source_name: str, schema) -> None:
    lines = []
    lines.append('<?xml version="1.0"?>')
    lines.append('<VTKFile type="PUnstructuredGrid" version="0.1" byte_order="LittleEndian">')
    lines.append('  <PUnstructuredGrid GhostLevel="0">')

    scalars_attr = f' Scalars="{schema[0][1]}"' if schema else ""
    lines.append(f'    <PPointData{scalars_attr}>')
    for dtype, name, ncomp in schema:
        ncomp_attr = "" if ncomp == "1" else f' NumberOfComponents="{ncomp}"'
        lines.append(f'      <PDataArray type="{dtype}" Name="{name}"{ncomp_attr}/>')
    lines.append('    </PPointData>')

    lines.append('    <PCellData>')
    lines.append('    </PCellData>')

    lines.append('    <PPoints>')
    lines.append('      <PDataArray type="Float32" NumberOfComponents="3"/>')
    lines.append('    </PPoints>')

    lines.append(f'    <Piece Source="{source_name}"/>')
    lines.append('  </PUnstructuredGrid>')
    lines.append('</VTKFile>')

    wrapper_path.write_text("\n".join(lines) + "\n")


def make_wrappers(directory: Path) -> None:
    vtu_files = sorted(directory.glob("*.vtu"))
    count = 0

    for vtu_path in vtu_files:
        match = VTU_PATTERN.match(vtu_path.name)
        if not match:
            continue

        label = match.group("label")
        num = match.group("num")
        block = match.group("block")

        schema = extract_point_data_schema(vtu_path)
        if not schema:
            print(f"Warning: no PointData arrays found in {vtu_path.name}, skipping")
            continue

        wrapper_name = f"b{block}_{label}{num}.pvtu"
        wrapper_path = directory / wrapper_name

        write_pvtu_wrapper(wrapper_path, vtu_path.name, schema)
        count += 1
        print(f"Wrote {wrapper_name}  -> {vtu_path.name}  "
              f"({len(schema)} variable(s): {', '.join(s[1] for s in schema)})")

    if count == 0:
        print(f"No files matching <label><NNN>_b<N>.vtu found in {directory}")


if __name__ == "__main__":
    target_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")
    make_wrappers(target_dir)