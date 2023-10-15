import numpy as np
import matplotlib.pyplot as plt


def create_model(a=5, b=2):
    def model(x, t=0):
        V, W = x
        dV = V**2 - a*W/(b+V)
        dW = -W + V/(b+V)
        return np.array([dV, dW])
    model.a = a
    model.b = b
    return model


def phase_plane_plot(model, range_x=(-1, 1), range_y=None,
                     num_grid_points=50, show=False, points=None):
    '''
    Simple implementation of the phase plane plot in matplotlib.

    Input:
    -----
      *model* : function
        function that takes numpy.array as input with two elements
        representing two state variables
      *range_x* = (-1, 1) : tuple
        range of x axis
      *range_y* = None : tuple
        range of y axis; if None, the same range as *range_x*
      *num_grid_points* = 50 : int
        number of samples on grid
      *show* = False : bool
        if True it shows matplotlib plot
      *points* = None : np.array
        array of points to draw on plot
    '''
    if range_y is None:
        range_y = range_x
    x_ = np.linspace(range_x[0], range_x[1], num_grid_points)
    y_ = np.linspace(range_y[0], range_y[1], num_grid_points)

    x_step = x_[1]-x_[0]
    y_step = y_[1]-y_[0]
    diag_ = np.linalg.norm([x_step, y_step])

    grid = np.meshgrid(x_, y_)

    dfmat = np.zeros((num_grid_points, num_grid_points, 2))
    for nx in range(num_grid_points):
        for ny in range(num_grid_points):
            df = model([grid[0][nx, ny], grid[1][nx, ny]])
            df_len = np.linalg.norm(df)
            if 0.3*diag_ > df_len or df_len > 5*diag_:
                df *= 2*diag_/df_len
            dfmat[nx, ny, 0] = df[0]
            dfmat[nx, ny, 1] = df[1]

    plt.grid(True)
    plt.title(f'$\\alpha$={model.a}, $\\beta$={model.b}')
    plt.quiver(grid[0], grid[1], dfmat[:, :, 0], dfmat[:, :, 1])
    if points is not None:
        plt.scatter(points[:, 0], points[:, 1], s=[150])
    # plt.contour(grid[0], grid[1], dfmat[:, :, 0], [0], colors='r')
    # plt.contour(grid[0], grid[1], dfmat[:, :, 1], [0], colors='g')
    if show:
        plt.show()


if __name__ == "__main__":
    phase_plane_plot(create_model(5, 2),
                     range_x=(0.4, 1.),
                     range_y=(0., 0.6),
                     num_grid_points=20, show=True,
                     points=np.array([[.6906474480, .2566844826]]))
