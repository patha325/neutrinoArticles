from pathlib import Path
import re
import shutil
import subprocess

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "src" / "first_draft.md"
WORK = ROOT / "work"
SUBMISSION = ROOT / "submission"

text = SOURCE.read_text(encoding="utf-8")
lines = text.splitlines()
title = lines[0][2:].strip()

# Keep the author line, but remove the working-draft label and status box.
lines = [line for line in lines if not line.startswith("**Working manuscript,")]
lines = [line for line in lines if line != f"**{'Sven-Patrik Hallsjö'}**  "]
lines = [line for line in lines if line != f"# {title}"]
text = "\n".join(lines) + "\n"
text = re.sub(r"\n> \*\*Status\.\*\*.*?(?=\n\n)", "", text, count=1, flags=re.S)
text = re.sub(r"commit `([0-9a-f]{40})`", r"commit \1", text)
text = re.sub(
    r"`P_ν,on = 0\.976 MW × \(L/1,000 km\)² × \(10 kt/M\) × \(θ/1 mrad\)² × \(R_slot/1 slot s⁻¹\)`\.\s*\(A1\)",
    r"$$\nP_{\\nu,\\mathrm{on}} = 0.976\\,\\mathrm{MW} \\times (L/1000\\,\\mathrm{km})^2 \\times (10\\,\\mathrm{kt}/M) \\times (\\theta/1\\,\\mathrm{mrad})^2 \\times (R_{\\mathrm{slot}}/1\\,\\mathrm{slot\\,s^{-1}}) \\tag{A1}\n$$",
    text,
)

# Point image references into the compact arXiv source tree and merge duplicate
# Markdown/standalone captions into one LaTeX figure caption.
text = re.sub(r"\]\(\.\./figures/0([1-7])_[^)]+\)", r"](figures/0\1.png)", text)
def merge_caption(match):
    alt, path, number, caption = match.groups()
    alt = re.sub(rf"^Figure {number}:\s*", "", alt)
    caption = re.sub(rf"^Figure {number}\.\s*", "", caption.strip())
    return f"![{alt} {caption}]({path})"
text = re.sub(
    r"!\[([^\]]*)\]\(([^)]+)\)\s*\n\n\*\*Figure (\d+)\.\*\* (.*?)(?=\n\n|\Z)",
    merge_caption,
    text,
    flags=re.S,
)

# Keep the oscillation illustration as two genuine inline equations; the source
# Markdown had wrapped an entire sentence (including English text) as code.
text = text.replace(
    "`Δm² = 2.5 × 10⁻³ eV²` and maximal mixing yields `P(νμ→νμ) ≈ 1 − sin²[1.267 Δm²(eV²)L(km)/E(GeV)] ≈ 0.29`",
    r"$\Delta m^2 = 2.5 \times 10^{-3}\,\mathrm{eV}^2$ and maximal mixing yields $P(\nu_\mu\to\nu_\mu) \approx 1 - \sin^2[1.267\,\Delta m^2(\mathrm{eV}^2)L(\mathrm{km})/E(\mathrm{GeV})] \approx 0.29$",
)
text = text.replace("(`0.1 kt`)", r"($0.1\\,\\mathrm{kt}$)")
text = text.replace("(`0.01 kt`)", r"($0.01\\,\\mathrm{kt}$)")

# Convert standalone numbered formula lines to display math. Inline expressions
# are converted below; filenames remain monospaced text.
def display_equation(match):
    formula, punctuation, number = match.groups()
    formula = formula.strip()
    if formula.endswith("."):
        formula = formula[:-1].rstrip()
    sup = {"⁰":"0", "¹":"1", "²":"2", "³":"3", "⁴":"4", "⁵":"5", "⁶":"6", "⁷":"7", "⁸":"8", "⁹":"9", "⁻":"-"}
    formula = re.sub(r"[⁻⁰¹²³⁴⁵⁶⁷⁸⁹]+", lambda m: "^{" + "".join(sup[c] for c in m.group()) + "}", formula)
    formula = formula.replace("½", r"\frac{1}{2}").replace("~", r"\sim")
    formula = re.sub(r"\bPoisson\b", r"\\operatorname{Poisson}", formula)
    formula = re.sub(r"\bexp\b", r"\\exp", formula)
    formula = re.sub(r"\bln\b", r"\\ln", formula)
    formula = re.sub(r" selected events per on symbol$", r" \\text{selected events per on symbol}", formula)
    return f"\n$$\n{formula} \\tag{{{number}}}\n$$\n"
text = re.sub(r"^`([^`]+)`([.,]?)\s*\((\d+|A\d+)\)\s*$", display_equation, text, flags=re.M)

math_symbols = {
    "λ̂": r"{\hat{\lambda}}", "λ": r"{\lambda}", "τ": r"{\tau}",
    "θ": r"{\theta}", "ν": r"{\nu}", "ε": r"{\varepsilon}",
    "μ": r"{\mu}", "σ": r"{\sigma}", "π": r"{\pi}",
    "Δ": r"{\Delta}", "≈": r"\approx", "≃": r"\simeq",
    "∝": r"\propto", "×": r"\times", "−": "-", "~": r"\sim",
    "₁": r"_{1}", "₂": r"_{2}", "₃": r"_{3}",
}
math_functions = {"exp", "ln", "sin", "cos", "log", "max", "min"}
math_pattern = re.compile(r"[=<>∝≈≃×−/|^_{}\[\]λτθν εμσπΔ₁₂₃⁻⁰¹²³⁴⁵⁶⁷⁸⁹]".replace(" ", ""))

def convert_inline(match):
    value = match.group(1)
    # Keep filenames and symbolic commit IDs as code.
    if re.search(r"\.(?:py|csv|md|txt)$", value) or re.fullmatch(r"[0-9a-f]{40}", value):
        return match.group(0)
    # A spaced code span is prose or a number-plus-unit unless it contains an
    # unambiguous relation; do not swallow adjacent English into math mode.
    if " " in value and not re.search(r"[=<>∝≈≃×/|]", value):
        return match.group(0)
    if not (math_pattern.search(value) or re.fullmatch(r"[A-Za-z]+", value) or re.search(r"\d", value)):
        return match.group(0)
    for old, new in math_symbols.items():
        value = value.replace(old, new)
    sup = {"⁰":"0", "¹":"1", "²":"2", "³":"3", "⁴":"4", "⁵":"5", "⁶":"6", "⁷":"7", "⁸":"8", "⁹":"9", "⁻":"-"}
    value = re.sub(r"[⁻⁰¹²³⁴⁵⁶⁷⁸⁹]+", lambda m: "^{" + "".join(sup[c] for c in m.group()) + "}", value)
    for name in math_functions:
        value = re.sub(rf"\b{name}\b", lambda _: "\\" + name, value)
    # Textual units in math mode should be upright.
    return f"${value}$"

Path("/tmp/pre-inline.md").write_text(text, encoding="utf-8")
text = re.sub(r"`([^`]+)`", convert_inline, text)
# Normalize escaped control sequences to one TeX command slash and preserve the
# repository script path as monospaced code.
text = text.replace("\\\\", "\\")
text = text.replace("$communication/analysis/solar_system_scaling.py`", "`communication/analysis/solar_system_scaling.py`")
# Some pre-existing Unicode superscripts can be mistaken for Markdown math
# boundaries; normalize any such nested markers generated in display math.
text = re.sub(r"\$\^\{([^}]+)\}\$", r"^{\1}", text)
# Plain English duration units must not be interpreted as the math operator
# \min when they occur in table cells.
text = re.sub(r"(?<=\d) \\min\b", " min", text)
text = text.replace("$communication/analysis/solar_system_scaling.py`", r"`communication/analysis/solar\_system\_scaling.py`")
# Make bare reference URLs breakable in the typeset bibliography.
def linkify_url(match):
    url = match.group(0)
    suffix = ""
    while url and url[-1] in ".,;":
        suffix = url[-1] + suffix
        url = url[:-1]
    return f"<{url}>{suffix}"
text = re.sub(r"https?://[^\s<>]+", linkify_url, text)
text = re.sub(r"^(#{2,6}) ", lambda m: "#" * (len(m.group(1)) - 1) + " ", text, flags=re.M)

WORK.mkdir(exist_ok=True)
SUBMISSION.mkdir(exist_ok=True)
(WORK / "paper.md").write_text(text, encoding="utf-8")
(WORK / "title.txt").write_text(title, encoding="utf-8")

shutil.rmtree(SUBMISSION / "figures", ignore_errors=True)
(SUBMISSION / "figures").mkdir(parents=True)
for i in range(1, 8):
    shutil.copy2(ROOT / "src" / "figures" / f"0{i}.png", SUBMISSION / "figures" / f"0{i}.png")

header = r'''\usepackage[margin=1in]{geometry}
\usepackage{caption}
\usepackage{float}
\usepackage{xurl}
\setlength{\emergencystretch}{3em}
\setcounter{secnumdepth}{-\maxdimen}
'''
(WORK / "header.tex").write_text(header, encoding="utf-8")

output = SUBMISSION / "main.tex"
subprocess.run(
    ["pandoc", str(WORK / "paper.md"), "--from=markdown+tex_math_dollars+pipe_tables+implicit_figures",
     "--to=latex", "--standalone", "--top-level-division=section",
     "--metadata", f"title={title}",
     "--metadata", "author=Sven-Patrik Hallsjö", "--metadata", "date=",
     "-H", str(WORK / "header.tex"), "-o", str(output)],
    check=True,
)

readme = '''Source for “%s”.\n\nCompile main.tex with XeLaTeX (run twice for references). The seven required PNG figures are in figures/. The author is formatted as in the author's earlier arXiv paper.\n''' % title
(SUBMISSION / "README.txt").write_text(readme, encoding="utf-8")
