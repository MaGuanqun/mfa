import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
import pandas as pd

name="persistent_diagram_1e_3.csv"
figure_name="persistent_diagram_1e_3.png"
persistent = pd.read_csv(name)

filtration_levels = np.linspace(0, 157, num=200)
max_counts = []

for level in filtration_levels:
    count = (persistent['Persistence']> level).sum()
    max_counts.append(count)

plt.plot(filtration_levels, max_counts, label='Max Count')
plt.ylabel('Maximum count')
plt.xlabel('Persistence')
plt.savefig(figure_name, dpi=600) # Saves the figure with high resolution
