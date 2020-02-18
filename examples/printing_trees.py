import matplotlib.pyplot as plt
import matplotlib

import networkx as nx
import pydot
from networkx.drawing.nx_pydot import graphviz_layout

T = nx.balanced_tree(2, 4)

plt.figure(figsize=(6,4))
pos = graphviz_layout(T, prog="dot")
pos[0] = (570,300)
nx.draw(T, pos, node_color="black", node_size=50)
plt.show()
matplotlib.use("pgf")
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
})
plt.savefig('tree___.pgf')
