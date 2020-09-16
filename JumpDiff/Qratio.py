## This is based on  Anvari, M., Tabar, M. R. R., Peinke, J., Lehnertz, K.,
# 'Disentangling the stochastic behavior of complex time series.' Scientific
# Reports, 6, 35435, 2016. doi: 10.1038/srep35435, now implemented by Leonardo
# Rydin Gorjão and Pedro G. Lind for direct application.

import numpy as np
from .moments import moments

def Qratio(lag: np.ndarray, timeseries: np.ndarray, loc: int=None,
        correction: bool=False) -> np.ndarray:
    r"""
    Qratio method to distinguish pure diffusion from jump-diffusion timeseries,
    Given by the relation of the 4th and 6th Kramers─Moyal coefficient with
    increasing lag

    .. math::

        Q(x,\tau) = \frac{D_6(x,\tau)}{5 D_4(x,\tau)} = \left\{\begin{array}{ll}
            b(x)^2 \tau, & \text{diffusive}  \\ \sigma_\xi^2(x), & \text{jumpy}
        \end{array}\right.

    Parameters
    ----------
    lag: np.ndarray of ints
        An array with the time-lag to extract the Kramers–Moyal coefficient for
        different lags.

    timeseries: np.ndarray
        A 1-dimensional timeseries.

    loc: float (defaul ``None``)
        Use a particular point in space to calculate the ratio. If ``None``
        given, the maximum of the probability density function is taken.

    corrections: bool (defaul ``False``)
        Select whether to use corrective terms.

    Returns
    -------
    lag: np.ndarray of ints
        Same as input, but only lag > 0 and as ints.

    ratio: np.ndarray of len(lag)
        Ratio of the sixth-order over forth-order Kramers–Moyal coefficient.

    References
    ----------
    Anvari, M., Tabar, M. R. R., Peinke, J., Lehnertz, K., 'Disentangling the
    stochastic behavior of complex time series.' Scientific Reports, 6, 35435,
    2016. doi: 10.1038/srep35435.

    Lehnertz, K., Zabawa, L., and Tabar, M. R. R., 'Characterizing abrupt
    transitions in stochastic dynamics.' New Journal of Physics, 20(11):113043,
    2018. doi: 10.1088/1367-2630/aaf0d7.
    """

    # Force lag to be ints, ensure lag > order + 1, and removes duplicates
    lag = lag[lag > 0]
    lag = np.round(np.unique(lag)).astype(int)

    # Assert if timeseries is 1 dimensional
    if timeseries.ndim > 1:
        assert timeseries.shape[1] == 1, "Timeseries needs to be 1-dimensional"

    # Find maximum of distribution
    if loc == None:
        temp = moments(timeseries, power=0, bins=np.array([5000]))[1]
        loc = np.argmax(temp[0])

    temp = moments(timeseries, power=6, bins=np.array([5000]), lag = lag,
        correction = correction)[1]
    ratio = temp[6,loc,:]/(5 * temp[4,loc,:])

    return lag, ratio
