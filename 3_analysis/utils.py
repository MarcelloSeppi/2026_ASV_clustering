import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from typing import Literal


# makes a long 1d set of plots as a grid of plots/subplots for convenient visualization
def make_gridplot(
    N: int,
    max_ny: int = 5,
    type: Literal["ax", "fig"] = "ax",
    figscale = 1,
    **kwargs
):
    """
    Create a figure containing a grid of either Matplotlib axes or subfigures.

    The grid is automatically sized to accommodate ``N`` elements, using at
    most ``max_ny`` columns. The returned axes or subfigures are flattened into
    a one-dimensional NumPy array for convenient iteration.

    Parameters
    ----------
    N : int
        Total number of axes or subfigures to create.
    max_ny : int, default=5
        Maximum number of columns in the grid. The number of rows is computed
        automatically as ``ceil(N / max_ny)``.
    type : {"ax", "fig"}, default="ax"
        Type of objects to create:
        - ``"ax"``: return a grid of `matplotlib.axes.Axes`.
        - ``"fig"``: return a grid of `matplotlib.figure.SubFigure`.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The parent figure.
    items : numpy.ndarray
        A one-dimensional array containing either the created axes or
        subfigures, depending on the value of ``type``.

    Raises
    ------
    ValueError
        If ``type`` is not one of ``{"ax", "fig"}``.

    Examples
    --------
    Create six subplots:

    >>> fig, axs = make_gridplot(6)
    >>> for ax in axs[:6]:
    ...     ax.plot([0, 1], [0, 1])

    Create four subfigures:

    >>> fig, subfigs = make_gridplot(4, type="fig")
    >>> ax = subfigs[0].subplots()
    """

    
    ncols = min(N, max_ny)
    nrows = int(np.ceil(N / ncols))

    
    

    fig = plt.figure(figsize=(figscale* 4 * ncols, figscale * 2.8 * nrows),constrained_layout=True)

    if type == "ax":
        axs = fig.subplots(nrows, ncols,**kwargs)
        axs = np.atleast_1d(axs).ravel()
        return fig, axs

    elif type == "fig":
        subfigs = fig.subfigures(nrows, ncols)
        subfigs = np.atleast_1d(subfigs).ravel()
        return fig, subfigs

    else:
        raise ValueError(f"Unknown type: {type!r}")