import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = json.load(open("/tmp/traj_data.json"))
ts, ed = d["ts"], d["ed"]
colors = {"D2": "tab:blue", "D3": "tab:green", "D4": "tab:red", "D5": "tab:purple"}

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
ax1.plot(ts, ed, "k-o", lw=2.5, ms=5, label="ED (exact, basis 743)", zorder=10)
for k in ["D2", "D3", "D4", "D5"]:
    if k not in d:
        continue
    ax1.plot(ts, d[k], "--o", color=colors[k], ms=3, alpha=0.85, label="iPEPS " + k.replace("D", "D="))
    ax2.plot(ts, [abs(a - b) for a, b in zip(d[k], ed)], "-o", color=colors[k], ms=3, label=k.replace("D", "D="))

for ax in (ax1, ax2):
    ax.axvline(2.6, color="gray", ls=":", alpha=0.6)
ax1.set_ylabel("n(t)  excitation density")
ax1.set_title("4x4 Neel PXP: exact-oracle iPEPS n(t) vs ED (rel_floor=1e-3)")
ax1.legend(ncol=2, fontsize=9)
ax2.set_yscale("log")
ax2.set_ylabel("|n_iPEPS(t) - n_ED(t)|")
ax2.set_xlabel("t")
ax2.legend(ncol=2, fontsize=9)
plt.tight_layout()
plt.savefig("/tmp/revival_trajectory.png", dpi=130)
print("SAVED /tmp/revival_trajectory.png")
