"""Render the conditional index summary from committed certificate counts.

Generated displays are not new spectral computations or independent proofs.
Run from any directory: python figures/plot_summary.py
"""
import json
from pathlib import Path
import runpy

import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

ROOT = Path(__file__).resolve().parents[1]
FIGURES = ROOT / 'figures'
plt.style.use(FIGURES / 'pub.mplstyle')
source = ROOT / 'scripts' / 'cases.py'
if not source.exists():
    source = ROOT / 'certificates' / 'cases.py'
case = runpy.run_path(str(source))['CASES'][4]
rows = []
for n in range(4, 9):
    if n == 4:
        counts = [case['K'][m] for m in range(2, 2 * n - 1)]
    else:
        cert = json.loads((ROOT / 'degrees5-8' / 'certificates' / f'lower_n{n}.json').read_text())
        counts = [cert['channels'][str(m)]['count'] for m in range(2, 2 * n - 1)]
    index = 2 * sum(counts)
    bound = n * (n - 1)
    rows.append(dict(degree=n, radial_counts=counts, conditional_real_index=index,
                     conjectured_upper_bound=bound, deficit=bound - index))
data = ROOT / 'data'
data.mkdir(exist_ok=True)
(data / 'morse-indices.json').write_text(json.dumps({
    'claim_label': 'CERTIFIED finite counts; full-plane indices conditional on analytic hypotheses',
    'source': 'scripts/cases.py (release: certificates/cases.py) and degrees5-8/certificates/lower_n*.json; upper matching verified by replay_all.py',
    'rows': rows,
}, indent=2) + '\n', encoding='utf-8', newline='\n')

fig, ax = plt.subplots(figsize=(6.3, 3.6))
n = [row['degree'] for row in rows]
index = [row['conditional_real_index'] for row in rows]
bound = [row['conjectured_upper_bound'] for row in rows]
ax.plot(n, bound, color='#D55E00', linestyle='--', marker='s',
        markerfacecolor='white', markersize=5, label=r'Conjectured bound $n(n-1)$')
ax.plot(n, index, color='#0072B2', linestyle='-', marker='o',
        markersize=4, label='Conditional real index')
for row in rows:
    if row['deficit']:
        x = row['degree']
        lo, hi = row['conditional_real_index'], row['conjectured_upper_bound']
        ax.plot([x, x], [lo, hi], color='0.35', linestyle=':', linewidth=0.9)
        ax.annotate(str(row['deficit']), (x, (lo + hi) / 2), xytext=(7, 0),
                    textcoords='offset points', va='center', fontsize=8)
ax.set(xlabel=r'Vortex degree $n$', ylabel='Real Morse index', xlim=(3.8, 8.35), ylim=(7, 61))
ax.set_xticks(n)
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.legend(loc='upper left')
fig.subplots_adjust(left=0.13, right=0.98, bottom=0.17, top=0.98)
fig.savefig(FIGURES / 'morse-indices.pdf', metadata={'Creator': 'plot_summary.py', 'CreationDate': None})
fig.savefig(FIGURES / 'morse-indices-preview.png', dpi=180)
plt.close(fig)
print(json.dumps(rows, indent=2))
