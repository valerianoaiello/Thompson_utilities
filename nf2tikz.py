#!/usr/bin/env python3
"""
nf2tikz.py - From the normal form of an element of Thompson's group F
to its reduced tree diagram, drawn in TikZ.

Conventions (Belk, "Thompson's Group F", Sections 1.2-1.3):
  * domain tree on top, range tree on the bottom, leaves aligned;
  * [T] = positive element with top tree T and a right vine as bottom tree;
  * [T] x_n = [T ^ n]  (Prop. 1.3.3): x_n attaches a caret to leaf n of T;
  * f = p q^{-1}  =>  tree diagram  [T_p ; T_q].

Usage:
  python3 nf2tikz.py "x0 x1 x6 x4^-1 x2^-1"
  python3 nf2tikz.py "x_0^2 x_1 x_3^{-1} x_0^{-1}" --scale 0.5 --standalone
  python3 nf2tikz.py "x0 x1" --labels            (leaves numbered 0,1,...)
  python3 nf2tikz.py "x0" --labels "a,b,c"       (custom labels)
  python3 nf2tikz.py "x0 x1" --caret-labels      (numbered carets; also with a list "a,b,c")
  python3 nf2tikz.py "x0 x1" --name "g="         (name to the left of the diagram)
  python3 nf2tikz.py "x0 x1" --glue              (leaves of the two trees glued together)
  (all labels are optional and can be combined)
  (always write the word BEFORE the options)
"""
import re
import sys
import argparse

LEAF = None  # a leaf; an internal node is a pair (left, right)


# ---------- parsing the word ----------
def parse_word(s):
    """Accepts 'x0 x1^2 x4^-1', 'x_0x_1^{2}x_4^{-1}', etc. Returns [(index, exponent)]."""
    s = s.replace(" ", "").replace("{", "").replace("}", "").replace("_", "")
    tokens = re.findall(r"x(\d+)(?:\^(-?\d+))?", s)
    if not tokens or "".join(f"x{i}" + (f"^{e}" if e else "") for i, e in tokens) != s:
        raise ValueError(f"Parola non riconosciuta: {s!r}")
    return [(int(i), int(e) if e else 1) for i, e in tokens]


def split_pq(word):
    """Splits the word into a positive part p and a negative part q^{-1}; returns the generator lists of p and q."""
    p, q, seen_neg = [], [], False
    for n, e in word:
        if e > 0:
            if seen_neg:
                raise ValueError("La parola deve essere della forma (positivi)(negativi), come la forma normale.")
            p += [n] * e
        elif e < 0:
            seen_neg = True
            q += [n] * (-e)
    q.reverse()  # q^{-1} = x_n^{-b_n}...x_0^{-b_0}  =>  q = x_0^{b_0}...x_n^{b_n}
    return p, q


# ---------- binary trees ----------
def width(t):
    return 0 if t is LEAF else width(t[0]) + width(t[1]) + 1


def right_vine(w):
    t = LEAF
    for _ in range(w):
        t = (LEAF, t)
    return t


def attach(t, n):
    """Attaches a caret to leaf n (leaves numbered from 0 on the left)."""
    if t is LEAF:
        assert n == 0
        return (LEAF, LEAF)
    wl = width(t[0]) + 1  # number of leaves of the left subtree
    if n < wl:
        return (attach(t[0], n), t[1])
    return (t[0], attach(t[1], n - wl))


def positive_tree(gens):
    """Tree T such that [T] = x_{g1} x_{g2} ... (Prop. 1.3.3, extending the right stalk when needed)."""
    t = LEAF
    for n in gens:
        while width(t) <= n:  # need n < w: add a caret to the last leaf (does not change the element)
            t = attach(t, width(t))
        t = attach(t, n)
    return t


def leaf_carets(t, offset=0, out=None):
    """Indices i of the carets whose two children are the leaves i, i+1."""
    if out is None:
        out = set()
    if t is LEAF:
        return out
    if t[0] is LEAF and t[1] is LEAF:
        out.add(offset)
    leaf_carets(t[0], offset, out)
    leaf_carets(t[1], offset + width(t[0]) + 1, out)
    return out


def remove_caret(t, i, offset=0):
    """Replaces the caret with leaves i, i+1 by a single leaf."""
    if t[0] is LEAF and t[1] is LEAF and offset == i:
        return LEAF
    wl = width(t[0]) + 1
    if i < offset + wl:
        return (remove_caret(t[0], i, offset), t[1])
    return (t[0], remove_caret(t[1], i, offset + wl))


def tree_diagram(word):
    p, q = split_pq(word)
    T, U = positive_tree(p), positive_tree(q)
    while width(T) < width(U):
        T = attach(T, width(T))
    while width(U) < width(T):
        U = attach(U, width(U))
    # reduction: remove pairs of opposing carets (Def. 1.2.3)
    while True:
        common = leaf_carets(T) & leaf_carets(U)
        if not common:
            return T, U
        i = min(common)
        T, U = remove_caret(T, i), remove_caret(U, i)


# ---------- TikZ ----------
def edges(t, sign, gap, offset=0, out=None):
    """Returns (x_node, y_node, [segments]). All edges are at 45 degrees."""
    if out is None:
        out = []
    if t is LEAF:
        return offset, sign * gap, out
    wl = width(t[0]) + 1
    xl, yl, _ = edges(t[0], sign, gap, offset, out)
    xr, yr, _ = edges(t[1], sign, gap, offset + wl, out)
    lo, hi = offset, offset + width(t)
    x, y = (lo + hi) / 2, sign * (gap + (hi - lo) / 2)
    out.append(((x, y), (xl, yl)))
    out.append(((x, y), (xr, yr)))
    return x, y, out


def caret_positions(t, sign, gap, offset=0, out=None):
    """Caret positions, indexed like the leaves: caret k separates leaves k and k+1."""
    if out is None:
        out = {}
    if t is LEAF:
        return out
    lo, hi = offset, offset + width(t)
    out[offset + width(t[0])] = ((lo + hi) / 2, sign * (gap + (hi - lo) / 2))
    caret_positions(t[0], sign, gap, offset, out)
    caret_positions(t[1], sign, gap, offset + width(t[0]) + 1, out)
    return out


def fmt(v):
    return f"{v:g}"


def to_tikz(T, U, scale=0.4, gap=0.5, label=None, leaf_labels=None,
            caret_labels=None, name=None, glue=False):
    n = width(T)
    if leaf_labels is not None and len(leaf_labels) != n + 1:
        raise ValueError(f"Servono {n + 1} etichette (una per foglia), ne hai date {len(leaf_labels)}.")
    if caret_labels is not None and len(caret_labels) != n:
        raise ValueError(f"Servono {n} etichette (una per caret), ne hai date {len(caret_labels)}.")
    if glue:
        gap = 0
    lines = [f"\\begin{{tikzpicture}}[scale={scale}, line cap=round]"]
    if label:
        lines.append(f"  % {label}")
    for title, tree, sign in (("domain tree (top)", T, 1), ("range tree (bottom)", U, -1)):
        lines.append(f"  % {title}")
        _, _, segs = edges(tree, sign, gap)
        for (a, b) in segs:
            lines.append(f"  \\draw ({fmt(a[0])},{fmt(a[1])}) -- ({fmt(b[0])},{fmt(b[1])});")
    if not glue:
        lines.append("  % leaves")
        for i in range(n + 1):
            lines.append(f"  \\fill ({i},{fmt(gap)}) circle (2.5pt) ({i},{fmt(-gap)}) circle (2.5pt);")
    if leaf_labels is not None:
        if glue:
            lines.append("  % leaf labels (on the glued vertices)")
            style = "font=\\scriptsize, fill=white, inner sep=1pt"
        else:
            lines.append("  % leaf labels (between the two rows)")
            style = "font=\\scriptsize"
        for i, lab in enumerate(leaf_labels):
            lines.append(f"  \\node[{style}] at ({i},0) {{${lab}$}};")
    if caret_labels is not None:
        lines.append("  % caret labels")
        for tree, sign, anchor in ((T, 1, "south"), (U, -1, "north")):
            pos = caret_positions(tree, sign, gap)
            for k, lab in enumerate(caret_labels):
                x, y = pos[k]
                lines.append(f"  \\node[font=\\scriptsize, anchor={anchor}, inner sep=2pt] "
                             f"at ({fmt(x)},{fmt(y)}) {{${lab}$}};")
    if name:
        lines.append("  % name of the tree diagram")
        lines.append(f"  \\node[anchor=east] at (-0.7,0) {{${name}$}};")
    lines.append("\\end{tikzpicture}")
    return "\n".join(lines)


def parse_labels(value, count):
    """'' -> 0,1,...,count-1; otherwise a comma-separated list."""
    if value == "":
        return [str(i) for i in range(count)]
    return [s.strip() for s in value.split(",")]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("word", help='es. "x0 x1 x6 x4^-1 x2^-1"')
    ap.add_argument("--scale", type=float, default=0.4)
    ap.add_argument("--gap", type=float, default=0.5, help="meta' della distanza verticale tra le due file di foglie")
    ap.add_argument("--standalone", action="store_true", help="produce un documento LaTeX compilabile")
    ap.add_argument("--labels", nargs="?", const="", default=None, metavar="LISTA",
                    help="etichette delle foglie: senza valore numera 0,1,...,w; "
                         "oppure una lista separata da virgole, es. \"a,b,c,d\"")
    ap.add_argument("--caret-labels", nargs="?", const="", default=None, metavar="LISTA",
                    help="etichette dei caret (in entrambi gli alberi): senza valore numera 0,1,...; "
                         "il caret k e' quello tra le foglie k e k+1")
    ap.add_argument("--name", metavar="TESTO",
                    help="nome a sinistra del diagramma, in modo matematico, es. \"g=\"")
    ap.add_argument("--glue", action="store_true",
                    help="incolla le foglie dei due alberi (figura chiusa, senza spazio tra le file)")
    a = ap.parse_args()

    T, U = tree_diagram(parse_word(a.word))
    n = width(T)
    leaf_labels = None if a.labels is None else parse_labels(a.labels, n + 1)
    caret_labels = None if a.caret_labels is None else parse_labels(a.caret_labels, n)
    gap = a.gap
    if leaf_labels is not None and gap == 0.5:  # default gap: widen it to make room for the labels
        gap = 0.8
    tikz = to_tikz(T, U, a.scale, gap, label=a.word, leaf_labels=leaf_labels,
                   caret_labels=caret_labels, name=a.name, glue=a.glue)
    if a.standalone:
        tikz = ("\\documentclass[tikz,border=4pt]{standalone}\n\\begin{document}\n"
                + tikz + "\n\\end{document}")
    print(tikz)


if __name__ == "__main__":
    try:
        main()
    except ValueError as e:
        sys.exit(f"Errore: {e}")
