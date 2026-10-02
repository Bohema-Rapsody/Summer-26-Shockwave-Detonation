import pandas as pd
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from math import floor

N_grid = 12
grid_size = 30


def chunk_id_to_coord(chunk_id):
    chunk_id -= 1

    Z_PLACE = chunk_id//(N_grid**2)
    Y_PLACE = (chunk_id%(N_grid**2))//N_grid
    X_PLACE = (chunk_id%(N_grid**2))%N_grid

    return (X_PLACE-N_grid/2+0.5)*grid_size,(Y_PLACE-N_grid/2+0.5)*grid_size,(Z_PLACE-N_grid/2+0.5)*grid_size


profiles = {}

with open("N2_shock/profiles/Ti_300_grid_properties.profile") as f:
    while True:

        line = f.readline()

        if not line:
            break

        if line.startswith("#"):
            continue

        words = line.split()

        if len(words) == 3:

            timestep = int(words[0])
            nchunks = int(words[1])

            # skip header
            #f.readline()

            data = []

            for i in range(nchunks):
                entry = list(map(float, f.readline().split()))

                if entry[1] > 100:
                    data.append(entry)

            for entry in data:

                coords = chunk_id_to_coord(entry[0])
                del(entry[0])
                entry.insert(0,coords[0])
                entry.insert(1,coords[1])
                entry.insert(2,coords[2])


            profiles[timestep] = pd.DataFrame(
                data,
                columns=[
                    "x",
                    "y",
                    "z",
                    "Ncount",
                    "temp",
                    "ke_atom",
                    "pe_atom"
                ],
            )


def _3D_plot(data,type):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    scatter = ax.scatter(
        data["x"],
        data["y"],
        data["z"],
        c=data[type],
        cmap="inferno"
    )

    fig.colorbar(scatter, label=type)

    ax.set_xlabel("x (Å)")
    ax.set_ylabel("y (Å)")
    ax.set_zlabel("z (Å)")

    plt.show()


def _2D_colour_map(data,type):
    data_grid = data.pivot(
        index="y",
        columns="x",
        values=type
    )


    plt.imshow(
        data_grid,
        origin="lower",
        aspect="equal",
        cmap="inferno"
    )

    plt.colorbar(label=type)
    plt.xlabel("x (Å)")
    plt.ylabel("y (Å)")
    plt.show()

def _2D_contours(data,type):

    data_grid = data.pivot(
        index="y",
        columns="x",
        values=type
    )

    plt.contourf(
        data_grid.columns,
        data_grid.index,
        data_grid.values,
        levels=20,
        cmap="inferno"
    )

    plt.contour(
        data_grid.columns,
        data_grid.index,
        data_grid.values,
        levels=20,
        colors="black",
        linewidths=0.5
    )

    plt.colorbar(label=type)
    plt.xlabel("x (Å)")
    plt.ylabel("y (Å)")
    plt.show()



def _3D_plot_animation(data,type):

    # ============================================================
    # Settings
    # ============================================================

    property_name = type

    # Dictionary keys = timesteps
    timesteps = sorted(data.keys())


    # ============================================================
    # Find global colour range
    # ============================================================

    global_min = min(
        df[property_name].min()
        for df in data.values()
    )

    global_max = max(
        df[property_name].max()
        for df in data.values()
    )


    # ============================================================
    # Initial timestep
    # ============================================================

    t0 = timesteps[0]
    df0 = data[t0]

    local_min = df0[property_name].min()
    local_max = df0[property_name].max()


    fig = go.Figure()


    # ------------------------------------------------------------
    # Trace 0: GLOBAL temperature scale
    # ------------------------------------------------------------

    fig.add_trace(
        go.Scatter3d(
            x=df0["x"],
            y=df0["y"],
            z=df0["z"],

            mode="markers",

            marker=dict(
                size=5,
                color=df0[property_name],
                colorscale="Inferno",

                cmin=global_min,
                cmax=global_max,

                colorbar=dict(
                    title=property_name,
                    x=1.05
                )
            ),

            name="Global"
        )
    )


    # ------------------------------------------------------------
    # Trace 1: LOCAL temperature scale
    # ------------------------------------------------------------

    fig.add_trace(
        go.Scatter3d(
            x=df0["x"],
            y=df0["y"],
            z=df0["z"],

            mode="markers",

            marker=dict(
                size=5,
                color=df0[property_name],
                colorscale="Inferno",

                cmin=local_min,
                cmax=local_max,

                colorbar=dict(
                    title=property_name,
                    x=1.05
                )
            ),

            name="Local",

            # Initially hidden
            visible=False
        )
    )


    # ============================================================
    # Create animation frames
    # ============================================================

    frames = []

    for timestep in timesteps:

        df = data[timestep]

        # Local range for THIS timestep
        local_min = df[property_name].min()
        local_max = df[property_name].max()


        # Two traces in every frame:
        #   trace 0 = global
        #   trace 1 = local

        frames.append(
            go.Frame(
                name=str(timestep),

                data=[

                    # --------------------------------------------
                    # GLOBAL
                    # --------------------------------------------

                    go.Scatter3d(
                        x=df["x"],
                        y=df["y"],
                        z=df["z"],

                        mode="markers",

                        marker=dict(
                            size=5,
                            color=df[property_name],
                            colorscale="Inferno",

                            cmin=global_min,
                            cmax=global_max
                        )
                    ),

                    # --------------------------------------------
                    # LOCAL
                    # --------------------------------------------

                    go.Scatter3d(
                        x=df["x"],
                        y=df["y"],
                        z=df["z"],

                        mode="markers",

                        marker=dict(
                            size=5,
                            color=df[property_name],
                            colorscale="Inferno",

                            cmin=local_min,
                            cmax=local_max
                        )
                    )
                ]
            )
        )


    fig.frames = frames


    # ============================================================
    # Timestep slider
    # ============================================================

    slider_steps = []

    for timestep in timesteps:

        slider_steps.append(
            dict(
                method="animate",

                args=[
                    [str(timestep)],

                    dict(
                        mode="immediate",

                        frame=dict(
                            duration=0,
                            redraw=True
                        ),

                        transition=dict(
                            duration=0
                        )
                    )
                ],

                label=str(timestep)
            )
        )


    # ============================================================
    # Layout
    # ============================================================

    fig.update_layout(

        scene=dict(
            xaxis_title="x (Å)",
            yaxis_title="y (Å)",
            zaxis_title="z (Å)"
        ),


        # --------------------------------------------------------
        # Timestep slider
        # --------------------------------------------------------

        sliders=[
            dict(
                active=0,

                currentvalue=dict(
                    prefix="Timestep: "
                ),

                steps=slider_steps
            )
        ],


        # --------------------------------------------------------
        # Dropdown + Play/Pause
        # --------------------------------------------------------

        updatemenus=[

            # ====================================================
            # Play / Pause
            # ====================================================

            dict(
                type="buttons",

                direction="left",

                x=0.0,
                y=1.15,

                buttons=[

                    dict(
                        label="▶ Play",

                        method="animate",

                        args=[
                            None,

                            dict(
                                frame=dict(
                                    duration=100,
                                    redraw=True
                                ),

                                transition=dict(
                                    duration=0
                                ),

                                fromcurrent=True
                            )
                        ]
                    ),

                    dict(
                        label="⏸ Pause",

                        method="animate",

                        args=[
                            [None],

                            dict(
                                frame=dict(
                                    duration=0,
                                    redraw=False
                                ),

                                mode="immediate"
                            )
                        ]
                    )
                ]
            ),


            # ====================================================
            # Temperature scale dropdown
            # ====================================================

            dict(
                type="dropdown",

                direction="down",

                x=0.25,
                y=1.15,

                buttons=[

                    # ------------------------------------------------
                    # GLOBAL
                    # ------------------------------------------------

                    dict(
                        label=f"Global {type} scale",

                        method="restyle",

                        args=[
                            {
                                "visible": [True, False]
                            }
                        ]
                    ),


                    # ------------------------------------------------
                    # LOCAL
                    # ------------------------------------------------

                    dict(
                        label=f"Local {type} scale",

                        method="restyle",

                        args=[
                            {
                                "visible": [False, True]
                            }
                        ]
                    )
                ]
            )
        ]
    )


    # ============================================================
    # Display
    # ============================================================
    #pio.renderers.default = "browser"
    #fig.show()
    fig.write_html(
        f"{type}_animation.html",
        auto_open=True
    )


coord_steps = [(i-N_grid/2+0.5)*grid_size for i in range(N_grid)]
print(f"Available coords: {coord_steps}")

profile = profiles[1000]
sliced_profile = profile[profile["z"] == coord_steps[floor(N_grid/2)]]

#_3D_plot(profile,"temp")
_3D_plot_animation(profiles,"temp")
_3D_plot_animation(profiles,"ke_atom")
_3D_plot_animation(profiles,"pe_atom")

