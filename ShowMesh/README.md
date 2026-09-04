# ShowMesh (SDL2 port)

This is a rewrite of `showmesh3.c` (the EasyMesh viewer originally by
Bojan Niceno) that replaces all the X11/Xlib calls with
[SDL2](https://www.libsdl.org/) + `SDL2_ttf`, so it builds and runs on
**both macOS and Linux** from the same source file (X11/Xlib itself
only exists on Linux/Unix, which is why the original wouldn't build on
a Mac at all).

All the mesh-file parsing and mesh math (`load_mesh`, the zoom/pan
transform, the Delaunay/Voronoi/boundary drawing logic) is unchanged
from the original. Only the windowing/graphics layer was replaced.

## What changed vs. the original

- **Windowing & drawing**: Xlib → SDL2. Buttons are drawn as SDL
  rectangles + text instead of X child windows.
- **Text**: X server-side fonts → `SDL2_ttf` rendering a TrueType font
  found on your system (or pointed to via `SHOWMESH_FONT`, see below).
  No font files are bundled.
- **Zoom/Move**: the original used a "click, then click again"
  interaction (driven by X's `XOR`-drawn rubber band). This port uses
  an ordinary click-drag-release, which is more natural with SDL and
  gives you live rubber-band feedback (drawn in red) while dragging.
- **Bug fix**: while porting, testing turned up a real bug in the
  original's line-visibility check (`DrawLine`'s culling test) — it
  could incorrectly hide a mesh edge that was still on-screen if its
  two endpoints stuck out past *opposite* edges of the view (e.g.
  right after a zoom). This port uses a correct bounding-box overlap
  test instead, so lines no longer vanish in that situation.
- Everything else — button layout, mutual exclusivity of
  Elements/Materials, Delaunay/Voronoi/Nodes/Sides/Boundary toggle
  rules, Fit, Quit — behaves the same as the original.

## Building

You need SDL2 and SDL2_ttf development packages, and `pkg-config`.

**macOS** (with [Homebrew](https://brew.sh)):
```
brew install sdl2 sdl2_ttf pkg-config
make
```

**Ubuntu/Debian:**
```
sudo apt install build-essential pkg-config libsdl2-dev libsdl2-ttf-dev
make
```

**Fedora:**
```
sudo dnf install gcc pkgconf-pkg-config SDL2-devel SDL2_ttf-devel
make
```

This produces a `showmesh` binary.

## Running

Same as the original:
```
./showmesh <mesh-file> [<mesh-file2> ...]
```

## Fonts

The program looks for a bold TrueType font in the usual system
locations (DejaVu/Liberation on Linux, Arial/Helvetica on macOS). If
it can't find one, it will print an error telling you to either
install a font package or set the `SHOWMESH_FONT` environment
variable to point at any `.ttf`/`.ttc` file you have, e.g.:

```
export SHOWMESH_FONT=/path/to/SomeBoldFont.ttf
./showmesh mymesh
```

On Linux, `sudo apt install fonts-dejavu-core` (or
`fonts-liberation`) is enough to satisfy this.

## Files

- `showmesh_sdl.c` — the full source.
- `Makefile` — cross-platform build via `pkg-config`.
