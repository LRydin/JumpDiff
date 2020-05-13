## This is an implementation of a simple integrator for a generic
# jump-diffusion process. This includes for a Euler─Mayurama and a Milstein
# integration scheme. To use the Milstein scheme the derivative of the diffusion
# function needs to be given. Created by Leonardo Rydin Gorjão and Pedro G. Lind

import numpy as np

def jdprocess(time: float, delta_t: float, a: callable, b: callable,
        xi: float, lamb: float, init: float=None, solver: str='Euler',
        b_prime: callable=None) -> np.ndarray:
    """
    Integrates a jump-diffusion process with drift a(x), diffusion b(x), jump
    amplitude xi (ξ), and jump rate lamb (λ).

        dXₜ = a(x)dt + b(x)dW(t) + ξdJ(t), with J Poisson with jump rate λ.

    This integrator has both an Euler─Mayurama and a Milstein method of
    integration. For Milstein one has to introduce the derivative of the
    diffusion term 'b', denoted 'b_prime'.

    Parameters
    ----------
    time: float > 0
        Total integration time. Positive float or int.

    delta_t: float > 0
        Time sampling, the smaller the better.

    a: callable
        The drift function. Can be a function of a lambda. For an
        Ornstein─Uhlenbeck process with drift '-2x', a takes the form
            a =  lambda x: -2x

    b: callable
        The diffusion function. Can be a function of a lambda. For an
        Ornstein─Uhlenbeck process with diffusion 1, a takes the form
            b =  lambda x: 1

    xi: float > 0
        Variance of the jump amplitude, which will be turned into a normal
        distribution like N(0,√xi)

    lamb: float > 0
        Jump rate of the Poissonian jumps. This is implemented as the numpy
        function 'np.random.poisson(lam = lamb * delta_t)'

    init: float
        Initial conditions. If None given, generates a random value from a
        normal distribution ~ N(0,√delta_t)

    solver: 'Euler' or 'Milstein'
        The regular Euler─Maruyama solver (None) is the default, with an order
        of '√delta_t'. To employ a state-dependent diffusion, i.e., b(x) as a
        function of x, the Milstein scheme has an order of 'delta_t'. You must
        introduce as well the derivative of b(x), i.e., b'(x), as the argument
        'b_prime'.

    Returns
    -------
    X: np.array
        Timeseries of size int(time/delta_t).
    """

    # assert and conditions
    assert time > 0, "Total integration time must be positive"
    assert delta_t > 0, "Time sampling must be positive"
    if solver == 'Milstein':
        assert b_prime == None, "Introduce b'(x) to use the Milstein solver"
        assert callable(b_prime) == True, "b'(x) must be a function"

    assert callable(a) == True, "drift a(x) must be a function"
    assert callable(b) == True, "diffusion b(x) must be a function"
    assert type(lamb) == int or type(lamb) == float, ("'lamb' is not an int or "
            "float")
    assert type(xi) == int or type(xi) == float, "'xi' is not an int or float"


    # Define total length of timeseries
    length = int(time/delta_t)

    # Initialise the array X
    X = np.zeros(length)

    # randomise initial starting value or use given (after assert)
    if init is None:
        X[0] = np.random.normal(loc = 0, scale = np.sqrt(delta_t), size = 1)
    else:
        assert type(init) == int or type(init) == float, ("'init' is not an int"
            " or float")
        X[0] = float(init)

    # Generate the Gaussian noise
    dw = np.random.normal(loc = 0, scale = np.sqrt(delta_t), size = length)

    # Generate the Poissonian Jumps
    dJ = np.random.poisson(lam = lamb * delta_t, size = length)


    # Integration, either Euler
    if solver == 'Euler':
        for i in range(1, length):
            X[i] = X[i-1] + a(X[i-1]) * delta_t + b(X[i-1]) * dw[i]
            if dJ[i] > 0.:
                X[i] += dJ[i] * np.random.normal(loc = 0, scale = np.sqrt(xi))

    if solver == 'Milstein':
        # Generate corrective terms of the Milstein integration method
        dw_2 = (dw**2 - delta_t) * 0.5

        for i in range(1, length):
            X[i] = X[i-1] + a(X[i-1]) * delta_t + b(X[i-1]) * dw[i] \
                    + b(X[i-1]) * b_prime(X[i-1]) * dw_2[i]
            if dJ[i] > 0.:
                X[i] += dJ[i] * np.random.normal(loc = 0, scale = np.sqrt(xi))

    return X
