import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Milestone figure for the README: simple-update iPEPS reproduces the ED
# scar collapse-and-revival of the 2D PXP model. Data from
# scripts/dev_revival_trajectory.jl (exact 16-site oracle, no CTM).
d = json.load(open("/tmp/traj_data.json"))
ts, ed = d["ts"], d["ed"]

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(ts, ed, "k-o", lw=3, ms=6, label="ED (exact, basis 743)", zorder=10)
for k, c, lbl in [("D2", "tab:blue", "iPEPS D=2"),
                  ("D3", "tab:green", "iPEPS D=3"),
                  ("D4", "tab:red", "iPEPS D=4")]:
    if k in d:
        ax.plot(ts, d[k], "--o", color=c, ms=4, lw=1.8, alpha=0.9, label=lbl)

ax.annotate("collapse", xy=(1.2, 0.10), xytext=(1.15, 0.27), ha="center",
            fontsize=10, arrowprops=dict(arrowstyle="->", alpha=0.5))
ax.annotate("first revival", xy=(2.6, 0.483), xytext=(2.25, 0.39), ha="center",
            fontsize=10, arrowprops=dict(arrowstyle="->", alpha=0.5))
ax.set_xlabel("time  t")
ax.set_ylabel("excitation density  n(t)")
ax.set_title("2D PXP quench: simple-update iPEPS reproduces the ED scar revival\n"
             "(4x4 Neel, exact 16-site contraction, rel_floor=1e-3)")
ax.legend(loc="upper right", fontsize=9)
ax.grid(alpha=0.2)
plt.tight_layout()
plt.savefig("artifacts/revival_ipeps_vs_ed.png", dpi=140)
print("SAVED artifacts/revival_ipeps_vs_ed.png")
