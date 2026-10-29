# Thompson_utilities

Tools for Thompson's group F.

## nf2tikz

Takes the normal form of an element of F and draws its reduced tree diagram in TikZ
(conventions of Belk, *Thompson's Group F*, Sections 1.2–1.3: domain tree on top,
range tree on the bottom).

It needs only Python 3; compiling the figures requires LaTeX with TikZ.

```bash
python3 nf2tikz.py "x0 x1 x4^-1 x2^-1" --standalone > g.tex
pdflatex g.tex
```

Write the word before the options. Forms such as `"x_0^2 x_1 x_3^{-1}"` are accepted too.

### Options

| Option | Effect |
|---|---|
| `--standalone` | outputs a compilable LaTeX document (otherwise just the `tikzpicture`) |
| `--scale S` | scale of the picture (default 0.4) |
| `--gap G` | half the vertical distance between the two rows of leaves |
| `--labels` | numbers the leaves 0, 1, …; `--labels "a,b,c"` gives custom labels |
| `--caret-labels` | numbers the carets (caret k is the one between leaves k and k+1); also accepts a list |
| `--name "g="` | name to the left of the diagram, in math mode |
| `--glue` | glues the leaves of the two trees together |

All labels are optional and can be combined, for example:

```bash
python3 nf2tikz.py "x0 x1 x4^-1 x2^-1" --name "g=" --glue --labels "a,b,c,d,e,f,g" --standalone
```

`esempio_con_etichette.pdf` shows a sample output.

## nf2knot

Takes the normal form of an element of F (generators `x`) or of the Brown–Thompson group
F_3 (generators `y`) and draws, in TikZ, the link produced by Jones' construction
(see V. Aiello, *An introduction to Thompson knot theory and to Jones subgroups*,
JKTR 34 (2025), Sections 2, 4, 5):

* an element of F is sent to F_3 by Jones' map ι (x_n ↦ y_{2n});
* the two ternary trees are glued at the leaves and the roots are joined by an edge
  around the left side;
* every 4-valent vertex becomes a crossing, with the strand joining the left and
  right children over the one joining the parent and the middle child.

With `--oriented` the element must belong to the oriented subgroup (vec F or vec F_3),
i.e. its Tait graph must be 2-colourable; the link is then oriented as the boundary of
the checkerboard surface (leftmost region coloured +, + regions oriented counterclockwise).

```bash
python3 nf2knot.py "x0 x2^2 x5 x6 x7^-1 x6^-1 x4^-1" --standalone > g.tex   # the figure-eight knot
python3 nf2knot.py "x0 x1" --oriented --standalone > h.tex
python3 nf2knot.py "y1^2" --oriented --jones
```

| Option | Effect |
|---|---|
| `--oriented` | oriented link (error if the element is not in vec F / vec F_3) |
| `--standalone` | outputs a compilable LaTeX document |
| `--scale S` | scale of the picture (default 0.6) |
| `--pd` | prints the PD code instead of the picture |
| `--jones` | prints the Jones polynomial instead of the picture |

Without `--oriented` the components are oriented arbitrarily for the PD code and the
Jones polynomial (for knots this makes no difference).
