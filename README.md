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
