import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider


def grid_centers(n: int) -> np.ndarray:
    """Return n cell-center coordinates in [0, 1]."""
    return (np.arange(n) + 0.5) / n


def wall_azimuth(sx, sy, sz, wall: str, n: int) -> np.ndarray:
    """
    Compute local azimuth (0~360 deg) on a wall grid.
    Wall definitions:
      - 'x0': plane x=0, local axes (y, z)
      - 'y0': plane y=0, local axes (x, z)
      - 'z0': plane z=0, local axes (x, y)
    """
    c = grid_centers(n)

    if wall == "x0":
        yy, zz = np.meshgrid(c, c)
        az = np.degrees(np.arctan2(zz - sz, yy - sy))
        return (az + 360) % 360

    if wall == "y0":
        xx, zz = np.meshgrid(c, c)
        az = np.degrees(np.arctan2(zz - sz, xx - sx))
        return (az + 360) % 360

    if wall == "z0":
        xx, yy = np.meshgrid(c, c)
        az = np.degrees(np.arctan2(yy - sy, xx - sx))
        return (az + 360) % 360

    raise ValueError(f"Unknown wall: {wall}")


def rectangular_solid_angle(x1, x2, y1, y2, z) -> np.ndarray:
    """
    Exact solid angle of a rectangle seen from a point on its normal axis.
    Formula (corner-combination form):
      Omega = f(x2,y2) - f(x1,y2) - f(x2,y1) + f(x1,y1)
      f(x,y) = atan2(x*y, z*sqrt(x^2+y^2+z^2))
    Units: steradian (sr)
    """

    def f(x, y):
        return np.arctan2(x * y, z * np.sqrt(x * x + y * y + z * z))

    omega = f(x2, y2) - f(x1, y2) - f(x2, y1) + f(x1, y1)
    return np.clip(omega, 0.0, None)


def wall_solid_angle(sx, sy, sz, wall: str, n: int) -> np.ndarray:
    """Compute per-cell solid angle on one wall using exact rectangle formula."""
    e = np.linspace(0.0, 1.0, n + 1)
    u1 = e[:-1][None, :]
    u2 = e[1:][None, :]
    v1 = e[:-1][:, None]
    v2 = e[1:][:, None]

    if wall == "x0":
        # Plane x=0, local coordinates: u->y, v->z, distance zdist->sx
        return rectangular_solid_angle(u1 - sy, u2 - sy, v1 - sz, v2 - sz, sx)

    if wall == "y0":
        # Plane y=0, local coordinates: u->x, v->z, distance zdist->sy
        return rectangular_solid_angle(u1 - sx, u2 - sx, v1 - sz, v2 - sz, sy)

    if wall == "z0":
        # Plane z=0, local coordinates: u->x, v->y, distance zdist->sz
        return rectangular_solid_angle(u1 - sx, u2 - sx, v1 - sy, v2 - sy, sz)

    raise ValueError(f"Unknown wall: {wall}")


def draw_cube_walls(ax3d):
    c = np.linspace(0, 1, 2)
    yy, zz = np.meshgrid(c, c)
    xx, zz2 = np.meshgrid(c, c)
    xx2, yy2 = np.meshgrid(c, c)

    ax3d.plot_surface(np.zeros_like(yy), yy, zz, alpha=0.16, color="#4C78A8", edgecolor="none")
    ax3d.plot_surface(xx, np.zeros_like(xx), zz2, alpha=0.16, color="#F58518", edgecolor="none")
    ax3d.plot_surface(xx2, yy2, np.zeros_like(xx2), alpha=0.16, color="#54A24B", edgecolor="none")

    edges = [
        ((0, 0, 0), (1, 0, 0)),
        ((0, 0, 0), (0, 1, 0)),
        ((0, 0, 0), (0, 0, 1)),
        ((1, 0, 0), (1, 1, 0)),
        ((1, 0, 0), (1, 0, 1)),
        ((0, 1, 0), (1, 1, 0)),
        ((0, 1, 0), (0, 1, 1)),
        ((0, 0, 1), (1, 0, 1)),
        ((0, 0, 1), (0, 1, 1)),
        ((1, 1, 0), (1, 1, 1)),
        ((1, 0, 1), (1, 1, 1)),
        ((0, 1, 1), (1, 1, 1)),
    ]
    for p0, p1 in edges:
        ax3d.plot(
            [p0[0], p1[0]],
            [p0[1], p1[1]],
            [p0[2], p1[2]],
            color="gray",
            linewidth=0.9,
            alpha=0.75,
        )


def main():
    fig = plt.figure(figsize=(14, 8))
    gs = fig.add_gridspec(
        2,
        4,
        left=0.05,
        right=0.98,
        top=0.95,
        bottom=0.18,
        width_ratios=[1.25, 1.0, 1.0, 1.0],
        wspace=0.35,
        hspace=0.35,
    )

    ax3d = fig.add_subplot(gs[:, 0], projection="3d")
    ax_x0 = fig.add_subplot(gs[0, 1])
    ax_y0 = fig.add_subplot(gs[0, 2])
    ax_z0 = fig.add_subplot(gs[0, 3])
    ax_text = fig.add_subplot(gs[1, 1:])
    ax_text.axis("off")

    # Initial values
    init_pos = {"x": 0.50, "y": 0.50, "z": 0.50}
    init_grid = {"x0": 4, "y0": 4, "z0": 4}

    # Sliders
    s_ax_x = fig.add_axes([0.07, 0.12, 0.26, 0.03])
    s_ax_y = fig.add_axes([0.07, 0.08, 0.26, 0.03])
    s_ax_z = fig.add_axes([0.07, 0.04, 0.26, 0.03])
    s_ax_nx = fig.add_axes([0.40, 0.12, 0.17, 0.03])
    s_ax_ny = fig.add_axes([0.62, 0.12, 0.17, 0.03])
    s_ax_nz = fig.add_axes([0.84, 0.12, 0.12, 0.03])

    s_x = Slider(s_ax_x, "source x", 0.02, 0.98, valinit=init_pos["x"], valstep=0.01)
    s_y = Slider(s_ax_y, "source y", 0.02, 0.98, valinit=init_pos["y"], valstep=0.01)
    s_z = Slider(s_ax_z, "source z", 0.02, 0.98, valinit=init_pos["z"], valstep=0.01)
    s_nx = Slider(s_ax_nx, "x=0 grid", 2, 12, valinit=init_grid["x0"], valstep=1)
    s_ny = Slider(s_ax_ny, "y=0 grid", 2, 12, valinit=init_grid["y0"], valstep=1)
    s_nz = Slider(s_ax_nz, "z=0 grid", 2, 12, valinit=init_grid["z0"], valstep=1)

    cax1 = fig.add_axes([ax_x0.get_position().x1 + 0.005, ax_x0.get_position().y0, 0.008, ax_x0.get_position().height])
    cax2 = fig.add_axes([ax_y0.get_position().x1 + 0.005, ax_y0.get_position().y0, 0.008, ax_y0.get_position().height])
    cax3 = fig.add_axes([ax_z0.get_position().x1 + 0.005, ax_z0.get_position().y0, 0.008, ax_z0.get_position().height])

    # Static 3D scene (draw once)
    draw_cube_walls(ax3d)
    ax3d.set_xlim(0, 1)
    ax3d.set_ylim(0, 1)
    ax3d.set_zlim(0, 1)
    ax3d.set_box_aspect((1, 1, 1))
    ax3d.set_xlabel("X")
    ax3d.set_ylabel("Y")
    ax3d.set_zlabel("Z")
    ax3d.set_title("3D cube + point source")
    source_artist = ax3d.scatter(
        [init_pos["x"]], [init_pos["y"]], [init_pos["z"]], s=80, c="red", label="Point source"
    )
    ax3d.legend(loc="upper right")

    # Initial heatmaps (create once, update data later)
    om_x0_init = wall_solid_angle(init_pos["x"], init_pos["y"], init_pos["z"], "x0", init_grid["x0"])
    om_y0_init = wall_solid_angle(init_pos["x"], init_pos["y"], init_pos["z"], "y0", init_grid["y0"])
    om_z0_init = wall_solid_angle(init_pos["x"], init_pos["y"], init_pos["z"], "z0", init_grid["z0"])
    init_vmax = max(om_x0_init.max(), om_y0_init.max(), om_z0_init.max())

    im_x0 = ax_x0.imshow(om_x0_init, origin="lower", cmap="viridis", vmin=0, vmax=init_vmax, interpolation="nearest")
    im_y0 = ax_y0.imshow(om_y0_init, origin="lower", cmap="viridis", vmin=0, vmax=init_vmax, interpolation="nearest")
    im_z0 = ax_z0.imshow(om_z0_init, origin="lower", cmap="viridis", vmin=0, vmax=init_vmax, interpolation="nearest")

    ax_x0.set_title(f"x=0 wall solid angle ({init_grid['x0']}x{init_grid['x0']})")
    ax_y0.set_title(f"y=0 wall solid angle ({init_grid['y0']}x{init_grid['y0']})")
    ax_z0.set_title(f"z=0 wall solid angle ({init_grid['z0']}x{init_grid['z0']})")

    ax_x0.set_xlabel("u")
    ax_x0.set_ylabel("v")
    ax_y0.set_xlabel("u")
    ax_y0.set_ylabel("v")
    ax_z0.set_xlabel("u")
    ax_z0.set_ylabel("v")

    cb1 = fig.colorbar(im_x0, cax=cax1)
    cb2 = fig.colorbar(im_y0, cax=cax2)
    cb3 = fig.colorbar(im_z0, cax=cax3)
    cb1.set_label("sr")
    cb2.set_label("sr")
    cb3.set_label("sr")

    # Static + dynamic text (draw once, update dynamic part only)
    ax_text.text(
        0.01,
        0.86,
        "Solid angle per cell uses the exact rectangle formula:\n"
        "Omega = f(x2,y2) - f(x1,y2) - f(x2,y1) + f(x1,y1)\n"
        "f(x,y) = atan2(x*y, z*sqrt(x^2 + y^2 + z^2))\n"
        "Each cell uses its four boundaries, not the center-point approximation.",
        fontsize=11,
        va="top",
    )
    pos_text = ax_text.text(
        0.01,
        0.20,
        "",
        fontsize=12,
        weight="bold",
    )
    omega_text = ax_text.text(
        0.01,
        0.05,
        "",
        fontsize=10,
    )

    # Cache identical slider states to avoid duplicate redraw work
    last_state = {"value": None}

    def update(_=None):
        sx, sy, sz = s_x.val, s_y.val, s_z.val
        nx, ny, nz = int(s_nx.val), int(s_ny.val), int(s_nz.val)
        state = (sx, sy, sz, nx, ny, nz)
        if state == last_state["value"]:
            return
        last_state["value"] = state

        # Update point source position only
        source_artist._offsets3d = ([sx], [sy], [sz])

        # Update heatmap arrays and extent when grid size changes
        om_x0 = wall_solid_angle(sx, sy, sz, "x0", nx)
        om_y0 = wall_solid_angle(sx, sy, sz, "y0", ny)
        om_z0 = wall_solid_angle(sx, sy, sz, "z0", nz)
        vmax = max(om_x0.max(), om_y0.max(), om_z0.max())

        im_x0.set_data(om_x0)
        im_y0.set_data(om_y0)
        im_z0.set_data(om_z0)
        im_x0.set_clim(0, vmax)
        im_y0.set_clim(0, vmax)
        im_z0.set_clim(0, vmax)
        im_x0.set_extent((-0.5, nx - 0.5, -0.5, nx - 0.5))
        im_y0.set_extent((-0.5, ny - 0.5, -0.5, ny - 0.5))
        im_z0.set_extent((-0.5, nz - 0.5, -0.5, nz - 0.5))

        ax_x0.set_xlim(-0.5, nx - 0.5)
        ax_x0.set_ylim(-0.5, nx - 0.5)
        ax_y0.set_xlim(-0.5, ny - 0.5)
        ax_y0.set_ylim(-0.5, ny - 0.5)
        ax_z0.set_xlim(-0.5, nz - 0.5)
        ax_z0.set_ylim(-0.5, nz - 0.5)

        ax_x0.set_title(f"x=0 wall solid angle ({nx}x{nx})")
        ax_y0.set_title(f"y=0 wall solid angle ({ny}x{ny})")
        ax_z0.set_title(f"z=0 wall solid angle ({nz}x{nz})")

        pos_text.set_text(f"Current source position: ({sx:.2f}, {sy:.2f}, {sz:.2f})")
        omega_text.set_text(
            "Solid-angle totals (sr): "
            f"x=0 {om_x0.sum():.4f}, y=0 {om_y0.sum():.4f}, z=0 {om_z0.sum():.4f}"
        )
        fig.canvas.draw_idle()

    for s in (s_x, s_y, s_z, s_nx, s_ny, s_nz):
        s.on_changed(update)

    update()
    plt.show()


if __name__ == "__main__":
    main()
