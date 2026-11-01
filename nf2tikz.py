#!/usr/bin/env python3
"""
nf2tikz.py - From the normal form of an element of Thompson's group F
(generators x) or of the Brown-Thompson group F_3 (generators y) to its
reduced tree diagram, drawn in TikZ.

Conventions (Belk, "Thompson's Group F", Sections 1.2-1.3):
  * domain tree on top, range tree on the bottom, leaves aligned;
  * [T] = positive element with top tree T and a right vine as bottom tree;
  * [T] x_n = [T ^ n]  (Prop. 1.3.3): x_n attaches a caret to leaf n of T;
  * f = p q^{-1}  =>  tree diagram  [T_p ; T_q].
The same conventions are used for F_3, with ternary carets: y_n attaches a
ternary caret to leaf n, as in nf2knot.py (so x_n -> y_{2n} under Jones' map iota).

Usage:
  python3 nf2tikz.py "x0 x1 x6 x4^-1 x2^-1"
  python3 nf2tikz.py "x_0^2 x_1 x_3^{-1} x_0^{-1}" --scale 0.5 --standalone
  python3 nf2tikz.py "x0 x1" --labels            (leaves numbered 0,1,...)
  python3 nf2tikz.py "x0" --labels "a,b,c"       (custom labels)
  python3 nf2tikz.py "x0 x1" --caret-labels      (numbered carets; also with a list "a,b,c")
  python3 nf2tikz.py "x0 x1" --name "g="         (name to the left of the diagram)
  python3 nf2tikz.py "x0 x1" --glue              (leaves of the two trees glued together)
  python3 nf2tikz.py "y0 y3 y1^-1" --labels      (element of F_3)
  (all labels are optional and can be combined)
  (always write the word BEFORE the options)
"""
import re
import sys
import argparse

LEAF = None  # a leaf; an internal node is a pair (left, right)


# ---------- parsing the word ----------
def parse_word(s):
    """Accepts 'x0 x1^2 x4^-1' or 'y_0y_1^{2}y_4^{-1}'. Returns (k, [(index, exponent)]), k = 2 for F, 3 for F_3."""
    s = s.replace(" ", "").replace("{", "").replace("}", "").replace("_", "")
    tokens = re.findall(r"([xy])(\d+)(?:\^(-?\d+))?", s)
    if not tokens or "".join(f"{c}{i}" + (f"^{e}" if e else "") for c, i, e in tokens) != s:
        raise ValueError(f"Word not recognized: {s!r}")
    letters = {c for c, _, _ in tokens}
    if len(letters) > 1:
        raise ValueError("Use only generators x (group F) or only generators y (group F_3).")
    return (2 if letters.pop() == "x" else 3), [(int(i), int(e) if e else 1) for _, i, e in tokens]


def split_pq(word):
    """Splits the word into a positive part p and a negative part q^{-1}; returns the generator lists of p and q."""
    p, q, seen_neg = [], [], False
    for n, e in word:
        if e > 0:
            if seen_neg:
                raise ValueError("The word must have the form (positive part)(negative part), like the normal form.")
            p += [n] * e
        elif e < 0:
            seen_neg = True
            q += [n] * (-e)
    q.reverse()  # q^{-1} = x_n^{-b_n}...x_0^{-b_0}  =>  q = x_0^{b_0}...x_n^{b_n}
    return p, q


# ---------- k-ary trees (k = 2: F, k = 3: F_3); an internal node is the tuple of its children ----------
def width(t):
    """Number of leaves minus 1."""
    return 0 if t is LEAF else sum(width(c) + 1 for c in t) - 1


def attach(t, n, k):
    """Attaches a k-ary caret to leaf n (leaves numbered from 0 on the left)."""
    if t is LEAF:
        assert n == 0
        return (LEAF,) * k
    out = []
    for c in t:
        w = width(c) + 1
        if 0 <= n < w:
            c = attach(c, n, k)
        n -= w
        out.append(c)
    return tuple(out)


def positive_tree(gens, k):
    """Tree T such that [T] = g1 g2 ... (Prop. 1.3.3, extending the right stalk when needed)."""
    t = LEAF
    for n in gens:
        while width(t) <= n:  # need n < w: add a caret to the last leaf (does not change the element)
            t = attach(t, width(t), k)
        t = attach(t, n, k)
    return t


def leaf_carets(t, offset=0, out=None):
    """First leaves i of the carets whose children are all leaves."""
    if out is None:
        out = set()
    if t is LEAF:
        return out
    if all(c is LEAF for c in t):
        out.add(offset)
    for c in t:
        leaf_carets(c, offset, out)
        offset += width(c) + 1
    return out


def remove_caret(t, i, offset=0):
    """Replaces the caret whose children are the leaves i, i+1, ... by a single leaf."""
    if all(c is LEAF for c in t) and offset == i:
        return LEAF
    out = []
    for c in t:
        w = width(c) + 1
        if c is not LEAF and offset <= i < offset + w:
            c = remove_caret(c, i, offset)
        offset += w
        out.append(c)
    return tuple(out)


def tree_diagram(word, k=2):
    p, q = split_pq(word)
    T, U = positive_tree(p, k), positive_tree(q, k)
    while width(T) < width(U):
        T = attach(T, width(T), k)
    while width(U) < width(T):
        U = attach(U, width(U), k)
    # reduction: remove pairs of opposing carets (Def. 1.2.3)
    while True:
        common = leaf_carets(T) & leaf_carets(U)
        if not common:
            return T, U
        i = min(common)
        T, U = remove_caret(T, i), remove_caret(U, i)


# ---------- TikZ ----------
def edges(t, sign, gap, offset=0, out=None):
    """Returns (x_node, y_node, [segments]). In binary trees all edges are at 45 degrees."""
    if out is None:
        out = []
    if t is LEAF:
        return offset, sign * gap, out
    lo, hi = offset, offset + width(t)
    x, y = (lo + hi) / 2, sign * (gap + (hi - lo) / 2)
    ends = []
    for c in t:
        xc, yc, _ = edges(c, sign, gap, offset, out)
        ends.append((xc, yc))
        offset += width(c) + 1
    out += [((x, y), e) for e in ends]
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
        raise ValueError(f"{n + 1} labels are needed (one per leaf), {len(leaf_labels)} were given.")
    if caret_labels is not None and len(caret_labels) != n:
        raise ValueError(f"{n} labels are needed (one per caret), {len(caret_labels)} were given.")
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
    ap.add_argument("word", help='normal form, e.g. "x0 x1 x6 x4^-1 x2^-1" (F) or "y0 y3 y1^-1" (F_3)')
    ap.add_argument("--scale", type=float, default=0.4)
    ap.add_argument("--gap", type=float, default=0.5, help="half the vertical distance between the two rows of leaves")
    ap.add_argument("--standalone", action="store_true", help="output a compilable LaTeX document")
    ap.add_argument("--labels", nargs="?", const="", default=None, metavar="LIST",
                    help="leaf labels: without a value the leaves are numbered 0,1,...,w; "
                         "otherwise a comma-separated list, e.g. \"a,b,c,d\"")
    ap.add_argument("--caret-labels", nargs="?", const="", default=None, metavar="LIST",
                    help="caret labels (in both trees): without a value the carets are numbered 0,1,...; "
                         "caret k is the one between leaves k and k+1")
    ap.add_argument("--name", metavar="TEXT",
                    help="name to the left of the diagram, in math mode, e.g. \"g=\"")
    ap.add_argument("--glue", action="store_true",
                    help="glue the leaves of the two trees (closed picture, no space between the rows)")
    a = ap.parse_args()

    k, word = parse_word(a.word)
    if k == 3 and a.caret_labels is not None:
        raise ValueError("--caret-labels is available only for elements of F.")
    T, U = tree_diagram(word, k)
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
        sys.exit(f"Error: {e}")
