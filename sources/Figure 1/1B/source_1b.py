import matplotlib.pyplot as plt
import numpy as np
import sys
import os

def find_project_root(target_folder='10ksgt6ss'):
    current = os.path.abspath(os.getcwd())
    while True:
        if os.path.basename(current) == target_folder:
            return current
        parent = os.path.dirname(current)
        if parent == current:
            raise FileNotFoundError(f"'{target_folder}' not found in path hierarchy.")
        current = parent

project_root = find_project_root('10ksgt6ss')
sources_path = os.path.join(project_root, 'sources')
sys.path.insert(0, sources_path)


from working_dfs import li1, li1i3, li1i3i4b, li2,li3, li3i4b, df10

xx = [li3,li1,li1i3,li2,li3i4b,li1i3i4b]
xxx = [len(x) for x in xx]
xxx.insert(0, df10.query('nei_c =="Orphan"').assembly.nunique())
values = [x/df10.assembly.nunique() for x in xxx]
colors = ['#8ec740', '#f26c50', '#fbaf5f', '#43b5e3', '#8a56a3', '#4c4b42', 'gold']
labels = ['Orphan', 'T6SS i3', 'T6SS i1', 'T6SSs i1+i3', 'T6SS i2','T6SSs i3+i4b', 'T6SSs i1+i3+i4b']
labels = [f'{x} {y}'for x,y in zip(labels,xxx)]

# Configuration
n = len(values)
thickness = 0.10
center_radius = 0.15
outermost_radius = center_radius + n * thickness
arrow_target_x = outermost_radius + 0.15
arrow_gap = 0.02

fig, ax = plt.subplots(figsize=(8, 8))
ax.axis('equal')
ax.axis('off')

#Central empty circle
theta = np.linspace(0, 2 * np.pi, 500)
x = center_radius * np.cos(theta)
y = center_radius * np.sin(theta)
ax.fill(x, y, color='white', edgecolor='black', linewidth=0.8, zorder=10)

#Rings and arrows
for i, (val, color, label) in enumerate(zip(values, colors, labels)):
    outer_radius = center_radius + (n - i) * thickness
    inner_radius = outer_radius - thickness
    angle_start = np.pi / 2  # 90°
    theta_full = np.linspace(angle_start, angle_start - 2 * np.pi, 500)

    #Ring background
    x_outer = outer_radius * np.cos(theta_full)
    y_outer = outer_radius * np.sin(theta_full)
    x_inner = inner_radius * np.cos(theta_full[::-1])
    y_inner = inner_radius * np.sin(theta_full[::-1])
    ax.fill(np.concatenate([x_outer, x_inner]),
            np.concatenate([y_outer, y_inner]),
            color='white', edgecolor='black', linewidth=0.8, zorder=1)

    #Arc colored
    theta_val = np.linspace(angle_start, angle_start - 2 * np.pi * val, 300)
    x_outer_val = outer_radius * np.cos(theta_val)
    y_outer_val = outer_radius * np.sin(theta_val)
    x_inner_val = inner_radius * np.cos(theta_val[::-1])
    y_inner_val = inner_radius * np.sin(theta_val[::-1])
    ax.fill(np.concatenate([x_outer_val, x_inner_val]),
            np.concatenate([y_outer_val, y_inner_val]),
            color=color, zorder=5)

    # ➤ Arrow (horizontal)
    arrow_radius = (outer_radius + inner_radius) / 2
    x_start = arrow_radius * np.cos(angle_start)
    y_start = arrow_radius * np.sin(angle_start)
    x_end = arrow_target_x
    y_end = y_start
    ax.annotate("",
                xy=(x_end, y_end), xytext=(x_start, y_start),
                arrowprops=dict(arrowstyle="->", color='black', lw=0.5),
                zorder=6)

    #Label
    ax.text(x_end + arrow_gap, y_end, label, va='center', ha='left', fontsize=10)

#Extend limits to show arrows and labels
ax.set_xlim(-outermost_radius, arrow_target_x)
ax.set_ylim(-outermost_radius, outermost_radius)

plt.tight_layout()
plt.savefig("Figure 1B.png", dpi=300, bbox_inches='tight')

