#!/usr/bin/env python3
"""
nf2knot.py - From the normal form of an element of F or F_3 to the link
produced by Jones' construction, drawn in TikZ.

Construction (Aiello, "An introduction to Thompson knot theory and to Jones
subgroups", Sections 2, 4, 5; Jones, "On the construction of knots and links
from Thompson's groups"):
  * an element of F is first sent to F_3 by iota (Fig. 8): every binary caret
    gets a middle edge, so x_n -> y_{2n};
  * the two ternary trees are glued at the leaves and the two roots are joined
    by an edge going around the left side;
  * every 4-valent vertex becomes a crossing (Fig. 7): the strand joining the
    left and right children goes over, the one joining parent and middle child
    goes under (in both trees).
With --oriented the element must lie in the oriented subgroup (vec F or
vec F_3): the Tait graph Gamma (black regions of the checkerboard shading,
leftmost region black) must be 2-colourable. The leftmost vertex gets colour +,
and the boundary of a + region is oriented counterclockwise, that of a -
region clockwise.

Usage:
  python3 nf2knot.py "x0 x1"                     (unoriented link of x0 x1 in F)
  python3 nf2knot.py "x0 x1" --oriented           (oriented link, x0 x1 in vec F)
  python3 nf2knot.py "y0 y2 y1^-1" --oriented     (element of F_3)
  python3 nf2knot.py "x0 x2^2 x5 x6 x7^-1 x6^-1 x4^-1" --standalone
  python3 nf2knot.py "x0 x1" --pd --jones         (PD code and Jones polynomial)
  (always write the word BEFORE the options)
"""
import re
import sys
import argparse
from fractions import Fraction

LEAF = None  # a leaf; an internal node is a tuple of k children


# ---------- parsing the word ----------
def parse_word(s):
    """Accepts 'x0 x1^2 x4^-1' or 'y_0y_1^{2}y_4^{-1}'. Returns (letter, [(index, exponent)])."""
    s = s.replace(" ", "").replace("{", "").replace("}", "").replace("_", "")
    tokens = re.findall(r"([xy])(\d+)(?:\^(-?\d+))?", s)
    if not tokens or "".join(f"{c}{i}" + (f"^{e}" if e else "") for c, i, e in tokens) != s:
        raise ValueError(f"Word not recognized: {s!r}")
    letters = {c for c, _, _ in tokens}
    if len(letters) > 1:
        raise ValueError("Use only generators x (group F) or only generators y (group F_3).")
    return letters.pop(), [(int(i), int(e) if e else 1) for _, i, e in tokens]


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


# ---------- k-ary trees ----------
def leaves(t):
    return 1 if t is LEAF else sum(leaves(c) for c in t)


def attach(t, n, k):
    """Attaches a k-ary caret to leaf n (leaves numbered from 0 on the left)."""
    if t is LEAF:
        assert n == 0
        return (LEAF,) * k
    out = []
    for c in t:
        w = leaves(c)
        if 0 <= n < w:
            c = attach(c, n, k)
        n -= w
        out.append(c)
    return tuple(out)


def positive_tree(gens, k):
    """Tree T such that [T] = g1 g2 ..., the bottom tree being the right vine."""
    t = LEAF
    for n in gens:
        while leaves(t) - 1 <= n:  # leaf n must not be the last one: add a caret to the last leaf
            t = attach(t, leaves(t) - 1, k)
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
        offset += leaves(c)
    return out


def remove_caret(t, i, offset=0):
    """Replaces the caret whose children are the leaves i, i+1, ... by a single leaf."""
    if all(c is LEAF for c in t) and offset == i:
        return LEAF
    out = []
    for c in t:
        w = leaves(c)
        if c is not LEAF and offset <= i < offset + w:
            c = remove_caret(c, i, offset)
        offset += w
        out.append(c)
    return tuple(out)


def tree_diagram(word, k):
    p, q = split_pq(word)
    T, U = positive_tree(p, k), positive_tree(q, k)
    while leaves(T) < leaves(U):
        T = attach(T, leaves(T) - 1, k)
    while leaves(U) < leaves(T):
        U = attach(U, leaves(U) - 1, k)
    # reduction: remove pairs of opposing carets
    while True:
        common = leaf_carets(T) & leaf_carets(U)
        if not common:
            return T, U
        i = min(common)
        T, U = remove_caret(T, i), remove_caret(U, i)


def iota(t):
    """Jones' embedding F -> F_3: every binary caret gets a middle leaf (Fig. 8)."""
    return LEAF if t is LEAF else (iota(t[0]), LEAF, iota(t[1]))


# ---------- the 4-valent graph and the link ----------
# Slots of a crossing, in counterclockwise order. p = parent, l/m/r = left/middle/right child.
CCW = {"T": ["p", "l", "m", "r"], "B": ["p", "r", "m", "l"]}
PARTNER = {"p": "m", "m": "p", "l": "r", "r": "l"}  # l-r is the over strand, p-m the under one


class Diagram:
    def __init__(self, T, B, d=0.45):
        self.d = d
        self.m = leaves(T)
        self.nodes = {}  # id -> dict(tree, pos, lo, hi, kids=[(slot, child)], parent)
        self.leaf_parent = {}  # (tree, i) -> (node id, slot)
        self.roots = {"T": self._build(T, "T", 0, None), "B": self._build(B, "B", 0, None)}

    def _build(self, t, tree, lo, parent):
        if t is LEAF:
            return ("leaf", lo)
        nid = len(self.nodes)
        hi = lo + leaves(t) - 1
        sign = 1 if tree == "T" else -1
        node = dict(tree=tree, lo=lo, hi=hi, pos=((lo + hi) / 2, sign * (hi - lo) / 2), parent=parent, kids=[])
        self.nodes[nid] = node
        off = lo
        for slot, c in zip("lmr", t):
            child = self._build(c, tree, off, nid)
            node["kids"].append((slot, child))
            if child[0] == "leaf":
                self.leaf_parent[(tree, off)] = (nid, slot)
            off += leaves(c)
        return ("node", nid)

    # --- geometry ---
    def target(self, nid, slot):
        """Point at the other end of the edge leaving node nid through slot."""
        node = self.nodes[nid]
        if slot == "p":
            if node["parent"] is None:
                x, y = node["pos"]
                return (x, y + (1 if node["tree"] == "T" else -1))
            return self.nodes[node["parent"]]["pos"]
        child = dict(node["kids"])[slot]
        return (child[1], 0) if child[0] == "leaf" else self.nodes[child[1]]["pos"]

    def slot_point(self, nid, slot):
        (x, y), (tx, ty) = self.nodes[nid]["pos"], self.target(nid, slot)
        n = ((tx - x) ** 2 + (ty - y) ** 2) ** 0.5
        return (x + self.d * (tx - x) / n, y + self.d * (ty - y) / n)

    def toward(self, a, b):
        n = ((b[0] - a[0]) ** 2 + (b[1] - a[1]) ** 2) ** 0.5
        return (a[0] + self.d * (b[0] - a[0]) / n, a[1] + self.d * (b[1] - a[1]) / n)

    # --- arcs (edges of the 4-valent graph) ---
    def arcs(self):
        """Each arc: ends ((nid, slot), (nid, slot)), path (TikZ), arrow segment, regions (left, right)
        with respect to the direction from the first end to the second."""
        out = []
        for nid, node in self.nodes.items():
            for slot, child in node["kids"]:
                P = self.slot_point(nid, slot)
                if child[0] == "node":
                    c = self.nodes[child[1]]
                    Q = self.slot_point(child[1], "p")
                    if node["tree"] == "T":  # going down: the left side is the east side
                        regs = (c["hi"], c["lo"] - 1)
                    else:  # going up: the left side is the west side
                        regs = (c["lo"] - 1, c["hi"])
                    out.append(dict(ends=((nid, slot), (child[1], "p")), path=[("line", P, Q)], regs=regs))
                elif node["tree"] == "T":
                    i = child[1]
                    b, bslot = self.leaf_parent[("B", i)]
                    Q = self.slot_point(b, bslot)
                    L = (i, 0)
                    L1, L2 = self.toward(L, node["pos"]), self.toward(L, self.nodes[b]["pos"])
                    path = [("line", P, L1), ("curve", L1, L, L2), ("line", L2, Q)]
                    out.append(dict(ends=((nid, slot), (b, bslot)), path=path, regs=(i, i - 1)))
        tr, br = self.roots["T"][1], self.roots["B"][1]
        P, Q = self.slot_point(tr, "p"), self.slot_point(br, "p")
        Y = self.nodes[tr]["pos"][1] + 1
        path = [("loop", P, Q, Y)]
        out.append(dict(ends=((tr, "p"), (br, "p")), path=path, regs=(-1, self.m - 1)))
        return out

    def corner_pairs(self, nid):
        """The two pairs of opposite regions at a crossing."""
        node = self.nodes[nid]
        spans = []
        for _, child in node["kids"]:
            if child[0] == "leaf":
                spans.append((child[1], child[1]))
            else:
                c = self.nodes[child[1]]
                spans.append((c["lo"], c["hi"]))
        (lo, l_hi), (_, m_hi), (_, hi) = spans
        return (lo - 1, m_hi), (l_hi, hi)


def colouring(D):
    """2-colouring of the Tait graph (vertices = black regions = odd indices); None if impossible."""
    adj = {}
    for nid in D.nodes:
        for a, b in D.corner_pairs(nid):
            if a % 2 == 1:  # the black pair
                adj.setdefault(a, []).append(b)
                adj.setdefault(b, []).append(a)
    col = {-1: 1}
    stack = [-1]
    while stack:
        v = stack.pop()
        for w in adj.get(v, []):
            if w not in col:
                col[w] = -col[v]
                stack.append(w)
            elif col[w] == col[v]:
                return None
    return col


def orient(D, arcs, col):
    """Directs every arc so that + regions are on its left and - regions on its right."""
    for a in arcs:
        left, right = a["regs"]
        black, black_on_left = (left, True) if left % 2 == 1 else (right, False)
        a["forward"] = (col[black] == 1) == black_on_left


def traverse(D, arcs, oriented):
    """Labels the arcs 1..N along the components. Returns (labels, heads, number of components)."""
    by_end = {}
    for idx, a in enumerate(arcs):
        for e in a["ends"]:
            by_end[e] = idx
    label, head = {}, {}  # head[idx] = end of the arc where the traversal enters the next crossing
    n, comps = 0, 0
    for start in range(len(arcs)):
        if start in label:
            continue
        comps += 1
        idx = start
        a = arcs[idx]
        tail = a["ends"][0] if (a.get("forward", True) or not oriented) else a["ends"][1]
        while idx not in label:
            n += 1
            label[idx] = n
            e0, e1 = arcs[idx]["ends"]
            h = e1 if tail == e0 else e0
            head[idx] = h
            if oriented and arcs[idx]["forward"] != (h == e1):
                raise AssertionError("incoherent orientation")
            nxt = (h[0], PARTNER[h[1]])
            idx, tail = by_end[nxt], nxt
    return label, head, comps


def pd_code(D, arcs, label, head):
    """PD code X[i,j,k,l] (counterclockwise from the incoming under strand) and crossing signs."""
    slot_label, incoming = {}, set()
    for idx, a in enumerate(arcs):
        for e in a["ends"]:
            slot_label[e] = label[idx]
        incoming.add(head[idx])
    X, signs = [], []
    for nid, node in D.nodes.items():
        ccw = CCW[node["tree"]]
        start = ccw.index("p") if (nid, "p") in incoming else ccw.index("m")
        order = [ccw[(start + j) % 4] for j in range(4)]
        X.append([slot_label[(nid, s)] for s in order])
        signs.append(1 if (nid, order[3]) in incoming else -1)  # over strand going from l to j
    return X, signs


# ---------- Kauffman bracket and Jones polynomial ----------
def padd(p, q, c=1):
    r = dict(p)
    for e, v in q.items():
        r[e] = r.get(e, 0) + c * v
    return {e: v for e, v in r.items() if v}


def pmul(p, q):
    r = {}
    for e1, v1 in p.items():
        for e2, v2 in q.items():
            r[e1 + e2] = r.get(e1 + e2, 0) + v1 * v2
    return {e: v for e, v in r.items() if v}


DELTA = {2: -1, -2: -1}  # -A^2 - A^-2


def bracket(X):
    """Kauffman bracket <D> in A (dict exponent -> coefficient), with <O> = 1."""
    states = {frozenset(): {0: 1}}
    for (a, b, c, d) in X:
        new = {}
        for M, poly in states.items():
            for pairs, w in ((((a, b), (c, d)), 1), (((a, d), (b, c)), -1)):
                match, loops = dict(M), 0
                for p, q in pairs:
                    if p == q or match.get(p) == q:
                        match.pop(p, None)
                        match.pop(q, None)
                        loops += 1
                        continue
                    x = match.pop(p) if p in match else p
                    if x != p:
                        match.pop(x)
                    y = match.pop(q) if q in match else q
                    if y != q:
                        match.pop(y)
                    match[x], match[y] = y, x
                key = frozenset(match.items())
                term = {e + w: v for e, v in poly.items()}
                for _ in range(loops):
                    term = pmul(term, DELTA)
                new[key] = padd(new.get(key, {}), term)
        states = new
    total = states.get(frozenset(), {0: 1})
    # divide by delta once (the normalisation <O> = 1)
    q, r = {}, dict(total)
    while r:
        top = max(r)
        coef = -r[top]
        q[top - 2] = coef
        r = padd(r, pmul({top - 2: coef}, DELTA), -1)
    return q


def jones(X, signs):
    """Jones polynomial: dict Fraction(exponent of t) -> coefficient."""
    if not X:
        return {Fraction(0): 1}
    w = sum(signs)
    b = bracket(X)
    f = {3 * (-w): (-1) ** (w % 2)}  # (-A^3)^{-w}
    V = pmul(f, b)
    return {Fraction(-e, 4): v for e, v in V.items()}  # A = t^{-1/4}


def fmt_poly(V):
    terms = []
    for e in sorted(V):
        c = V[e]
        if e == 0:
            mono = ""
        elif e == 1:
            mono = "t"
        else:
            mono = "t^{" + (str(e.numerator) if e.denominator == 1 else f"{e.numerator}/{e.denominator}") + "}"
        coef = str(abs(c)) if (abs(c) != 1 or not mono) else ""
        terms.append(("-" if c < 0 else "+", coef + mono))
    if not terms:
        return "0"
    s = ("-" if terms[0][0] == "-" else "") + terms[0][1]
    for sg, t in terms[1:]:
        s += f" {sg} {t}"
    return s


# ---------- TikZ ----------
def fmt(v):
    return f"{v:.3f}".rstrip("0").rstrip(".")


def pt(p):
    return f"({fmt(p[0])},{fmt(p[1])})"


def path_tikz(path):
    s = ""
    for piece in path:
        if piece[0] == "line":
            s += (pt(piece[1]) if not s else "") + f" -- {pt(piece[2])}"
        elif piece[0] == "curve":
            s += (pt(piece[1]) if not s else "") + f" .. controls {pt(piece[2])} .. {pt(piece[3])}"
        else:  # the loop joining the two roots around the left side
            _, P, Q, Y = piece
            s += (f"{pt(P)} -- ({fmt(P[0])},{fmt(Y)}) -- (-1,{fmt(Y)}) -- (-1,{fmt(-Y)}) "
                  f"-- ({fmt(Q[0])},{fmt(-Y)}) -- {pt(Q)}")
    return s


def arrow(a):
    """Arrow head in the middle of the longest straight piece of the arc."""
    piece = a["path"][0]
    if piece[0] == "loop":
        _, P, Q, Y = piece
        A, B = (-1, Y), (-1, -Y)
    else:
        A, B = max((p for p in a["path"] if p[0] == "line"),
                   key=lambda p: (p[2][0] - p[1][0]) ** 2 + (p[2][1] - p[1][1]) ** 2)[1:]
    if not a["forward"]:
        A, B = B, A
    M = ((A[0] + B[0]) / 2, (A[1] + B[1]) / 2)
    n = ((B[0] - A[0]) ** 2 + (B[1] - A[1]) ** 2) ** 0.5
    E = (M[0] + 0.05 * (B[0] - A[0]) / n, M[1] + 0.05 * (B[1] - A[1]) / n)
    return f"  \\draw[->] {pt(M)} -- {pt(E)};"


def to_tikz(D, arcs, scale, oriented, comment):
    lines = [f"\\begin{{tikzpicture}}[scale={scale}, line cap=round, line join=round, >=stealth, thick]"]
    lines += [f"  % {c}" for c in comment]
    if not D.nodes:
        lines.append("  \\draw (0,0) circle (1);")
        if oriented:
            lines.append("  \\draw[->] (1,0) -- (1,0.05);")
        lines.append("\\end{tikzpicture}")
        return "\n".join(lines)
    lines.append("  % arcs")
    for a in arcs:
        lines.append(f"  \\draw[rounded corners=4pt] {path_tikz(a['path'])};")
    lines.append("  % under strands (parent - middle child)")
    for nid, node in D.nodes.items():
        P, M = D.slot_point(nid, "p"), D.slot_point(nid, "m")
        lines.append(f"  \\draw {pt(P)} .. controls {pt(node['pos'])} .. {pt(M)};")
    lines.append("  % over strands (left - right child)")
    for nid, node in D.nodes.items():
        L, R = D.slot_point(nid, "l"), D.slot_point(nid, "r")
        c = f"{pt(L)} .. controls {pt(node['pos'])} .. {pt(R)}"
        lines.append(f"  \\draw[white, line width=6pt, line cap=butt] {c};")
        lines.append(f"  \\draw {c};")
    if oriented:
        lines.append("  % orientation")
        lines += [arrow(a) for a in arcs]
    lines.append("\\end{tikzpicture}")
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("word", help='normal form, e.g. "x0 x1 x4^-1 x2^-1" (F) or "y0 y2 y1^-1" (F_3)')
    ap.add_argument("--oriented", action="store_true",
                    help="oriented link: the element must belong to vec F (or vec F_3)")
    ap.add_argument("--scale", type=float, default=0.6)
    ap.add_argument("--standalone", action="store_true", help="output a compilable LaTeX document")
    ap.add_argument("--pd", action="store_true", help="print the PD code instead of the picture")
    ap.add_argument("--jones", action="store_true", help="print the Jones polynomial instead of the picture")
    a = ap.parse_args()

    letter, word = parse_word(a.word)
    if letter == "x":
        T, U = tree_diagram(word, 2)
        T, U = iota(T), iota(U)
    else:
        T, U = tree_diagram(word, 3)
    D = Diagram(T, U)
    arcs = D.arcs() if D.nodes else []
    if a.oriented:
        col = colouring(D)
        if col is None:
            group = "vec F" if letter == "x" else "vec F_3"
            raise ValueError(f"the element is not in {group}: its Tait graph is not 2-colourable.")
        orient(D, arcs, col)
    label, head, comps = traverse(D, arcs, a.oriented) if arcs else ({}, {}, 1)
    X, signs = pd_code(D, arcs, label, head) if arcs else ([], [])

    if a.pd or a.jones:
        if a.pd:
            print("PD[" + ", ".join("X[" + ",".join(map(str, x)) + "]" for x in X) + "]")
        if a.jones:
            print(fmt_poly(jones(X, signs)))
        return
    comment = [a.word, f"{len(X)} crossings, {comps} component(s)"]
    if a.oriented:
        comment.append(f"writhe {sum(signs)}")
    tikz = to_tikz(D, arcs, a.scale, a.oriented, comment)
    if a.standalone:
        tikz = ("\\documentclass[tikz,border=4pt]{standalone}\n\\begin{document}\n"
                + tikz + "\n\\end{document}")
    print(tikz)


if __name__ == "__main__":
    try:
        main()
    except ValueError as e:
        sys.exit(f"Error: {e}")
