import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D


# ============================================================
# Parameters you can edit
# ============================================================

alpha = 3.5
N = 20

# Which p_n constraints to show
p_start = 1
p_end = 13

# Plotting rectangle
dmin, dmax = 0.0, 5.0
gmin, gmax = -0.6, 0.0

# Grid resolution
nx, ny = 800, 800

# Transparency of the coloured violated regions
region_alpha = 0.1

# Whether to save PNG as well as display inline
save_png = False
save_path = "sausage-plot.png"


# ============================================================
# Grid
# ============================================================

D = np.linspace(dmin, dmax, nx)
Ga = np.linspace(gmin, gmax, ny)

DD, GG = np.meshgrid(D, Ga)


# ============================================================
# Generalised binomial coefficient
# ============================================================

def binom_real(a, j):
    """
    Compute generalized binomial coefficient (a choose j)
    for real a and nonnegative integer j.
    """
    out = 1.0
    for k in range(j):
        out *= (a - k) / (k + 1)
    return out


# ============================================================
# Recurrence for coefficients, up to a positive factor
# ============================================================
#
# G(t) = exp{-delta(1-t) - gamma(1-t)^alpha}
#      = exp(-delta-gamma) * exp(sum_{j>=1} b_j t^j).
#
# p_n = exp(-delta-gamma) q_n,
# so sign(p_n) = sign(q_n).
#
# q_0 = 1,
#
# n q_n = sum_{j=1}^n j b_j q_{n-j}.
#
# b_1 = delta + alpha gamma,
# b_j = -gamma (-1)^j binom(alpha,j),  j >= 2.
#
# ============================================================

b = [None] * (N + 1)

b[1] = DD + alpha * GG

for j in range(2, N + 1):
    b[j] = -GG * ((-1) ** j) * binom_real(alpha, j)


q = [np.ones_like(DD)]

# neg_masks[n] is True where p_n < 0, equivalently q_n < 0.
neg_masks = [np.zeros_like(DD, dtype=bool)]

for n in range(1, N + 1):
    s = np.zeros_like(DD)

    for j in range(1, n + 1):
        s += j * b[j] * q[n - j]

    qn = s / n
    q.append(qn)
    neg_masks.append(qn < 0)


# ============================================================
# Colours
# ============================================================

num_constraints = p_end - p_start + 1

# tab20 gives reasonably distinct colours for up to about 20 constraints.
colors = plt.cm.tab20(np.linspace(0, 1, num_constraints))


# ============================================================
# Plot
# ============================================================

fig, ax = plt.subplots(figsize=(9, 6), facecolor="white")
ax.set_facecolor("white")


# Draw semi-transparent solid regions and matching boundaries
for idx, n in enumerate(range(p_start, p_end + 1)):
    mask = neg_masks[n].astype(float)
    color = colors[idx]

    # Semi-transparent solid region where p_n < 0.
    ax.contourf(
        DD,
        GG,
        mask,
        levels=[0.5, 1.5],
        colors=[color],
        alpha=region_alpha
    )

    # Boundary p_n = 0, using same colour but opaque.
    vals = q[n]

    if np.nanmin(vals) < 0 < np.nanmax(vals):
        ax.contour(
            DD,
            GG,
            vals,
            levels=[0],
            colors=[color],
            linewidths=1.2,
            alpha=1
        )

ax.axhline(0, color="black", linewidth=0.8)
ax.axvline(0, color="black", linewidth=0.8)

ax.set_xlim(dmin, dmax)
ax.set_ylim(gmin, gmax)
ax.tick_params(labelsize=14)

ax.set_xlabel(r"$\delta$", fontsize=16)
ax.set_ylabel(r"$\gamma$", fontsize=16)

ax.set_title(
    rf"violated regions $p_n<0$, "
    rf"$n={p_start},\ldots,{p_end}$, "
    rf"$\alpha={alpha}$",
    fontsize=16
)

ax.grid(True, alpha=0.25)


# ============================================================
# Legend
# ============================================================

legend_handles = []

for idx, n in enumerate(range(p_start, p_end + 1)):
    color = colors[idx]

    legend_handles.append(
        Patch(
            facecolor=color,
            edgecolor=color,
            alpha=region_alpha,
            label=rf"$p_{ {n} } < 0$"
        )
    )

ax.legend(
    handles=legend_handles,
    loc="center left",
    bbox_to_anchor=(1.02, 0.5),
    frameon=True,
    facecolor="white",
    framealpha=0.95,
    fontsize=12
)

fig.tight_layout()


if save_png:
    fig.savefig(save_path, dpi=220, bbox_inches="tight", facecolor="white")

plt.show()