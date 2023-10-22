import numpy as np

iteration = 0


def exit_condition(x, xt):
    return x == xt


def fill_D():
    D = ...  # TODO: initialize D with a random variable
    return D


# TODO: write input
A = ...
c = ...

# TODO: fill identity matrix
I = ...

while True:
    print("Iteration " + str(iteration))

    D = fill_D()

    At = A * D
    ct = D * c

    P = I - np.transpose(At) * np.invert(At * np.transpose(At)) * At
    cp = P * ct

    # TODO: initialize alpha and v
    alpha = ...
    v = ...

    row = ...
    xt = row + alpha / v * cp

    x = D * xt

    if exit_condition(x, xt):
        break
