## Delevoped by Leonardo Rydin Gorjão and Pedro G. Lind.

from sympy import bell, symbols, factorial, simplify

def M_formula(power, tau = True):
    r"""
    Generate the formula for the conditional moments with second-order
    corrections based on the relation with the ordinary Bell polynomials

    .. math::

        M_n(x^{\prime},\tau) \sim (n!)\tau D_n(x^{\prime}) + \frac{(n!)\tau^2}{2}
        \sum_{m=1}^{n-1}   D_m(x^{\prime})  D_{n-m}(x^{\prime})

    Parameters
    ----------
    power: int
        Desired order of the formula.

    Returns
    -------
    term: sympy.symbols
        Expression up to given ``power``.
    """
    init_sym = symbols('D:'+str(int(power+1)))[1:]
    sym = ()
    for i in range(1,power + 1):
        sym += (factorial(i)*init_sym[i-1],)

    if tau == True:
        t = symbols('tau')

        term = t*bell(power, 1, sym) + t**2*bell(power, 2, sym)

    else:
        term = bell(power, 1, sym) + bell(power, 2, sym)

    return term

def F_formula(power):
    r"""
    Generate the formula for the conditional moments with second-order
    corrections based on the relation with the ordinary Bell polynomials

    .. math::

        D_n(x) &=  \frac{1}{\tau (n!)} \bigg[ \hat{B}_{n,1}
        \left(M_1(x,\tau),M_2(x,\tau),\ldots,M_{n}(x,\tau)\right) \\
        &\qquad  \left.-\frac{\tau}{2} \hat{B}_{n,2}
        \left(M_1(x,\tau),M_2(x,\tau),\ldots,M_{n-1}(x,\tau)\right)\right].

    Parameters
    ----------
    power: int
        Desired order of the formula.

    Returns
    -------
    term: sympy.symbols
        Expression up to given ``power``.
    """

    init_sym = symbols('D:'+str(int(power+1)))[1:]
    sym = ()
    for i in range(1,power + 1):
        sym += (factorial(i)*init_sym[i-1],)

    term = (symbols('M'+str(int(power))) - bell(power, 2, sym))/factorial(power)

    return term


def F_formula_solver(power):
    r"""
    Generate the reciprocal relation of the moments to the Kramers─Moyal
    coefficients by sequential iteration.

    .. math::

        D_n(x) &=  \frac{1}{\tau (n!)} \bigg[ \hat{B}_{n,1}
        \left(M_1(x,\tau),M_2(x,\tau),\ldots,M_{n}(x,\tau)\right) \\
        &\qquad  \left.-\frac{\tau}{2} \hat{B}_{n,2}
        \left(M_1(x,\tau),M_2(x,\tau),\ldots,M_{n-1}(x,\tau)\right)\right].

    Parameters
    ----------
    power: int
        Desired order of the formula.

    Returns
    -------
    term: sympy.symbols
        Expression up to given ``power``.
    """
    power = power + 1
    terms_to_sub = []
    for i in range(1, power):
        terms_to_sub += [F_formula(i)]
        for j in range(i):
            terms_to_sub[i-1] = terms_to_sub[i-1].subs('D'+str(j+1),terms_to_sub[j])

    # for i in range
    # term = (symbols('M'+str(int(power))) + bell(power, 2, sym))/factorial(power)

    return (terms_to_sub[power-2]*factorial(power-1)).expand(basic = True)
