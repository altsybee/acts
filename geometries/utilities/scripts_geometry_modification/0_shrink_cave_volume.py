from pathlib import Path

p = Path("o2sim_geometry.gdml")

# Shrink the cave volume to the tracker volume

old = '<box name="cave" x="4000" y="4000" z="6000" lunit="cm"/>'
new = '<box name="cave" x="100" y="100" z="800" lunit="cm"/>'

text = p.read_text(encoding="utf-8")
text = text.replace(old, new)

print(old, "-->", new, "replaced")

p.write_text(text, encoding="utf-8")
