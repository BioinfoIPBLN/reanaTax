#!/usr/bin/env python3
"""
sahmi_stats.py -- the small statistical kit the microbial filters and tests need.

Kept separate because both bin/minimizer_filter.py and bin/host_kmer_filter.py
import it, and because it is the one part of this pipeline that reimplements
something a library would normally provide: the abundance modules run in a plain
`python:3.12` container with no SciPy, and adding NumPy to pull in one
correlation test would be a strange trade.

Two deliberate choices, both matching what R's `cor.test(method = "spearman")`
does in the situation this pipeline is actually in:

  * Ranks are averaged over ties. Taxonomic count data is full of ties - whole
    blocks of taxa share a read count of 2 - and R switches to the asymptotic
    t approximation the moment ties are present, so that is what is used here
    rather than an exact permutation null that the ties would invalidate.

  * A perfect correlation is scored at its EXACT floor, 2/n!, not at the zero
    the t approximation tends to. That is not a cosmetic difference. With four
    samples the smallest attainable two-sided p-value is 2/4! = 0.083, so no
    taxon can clear p < 0.05 no matter how clean its evidence is; five samples
    give 2/5! = 0.017 and the test starts to have power. Encoding the floor
    makes that limit enforce itself instead of being a caveat in a docstring.
"""
import math

# Smallest sample count at which a two-sided Spearman test can reach p < 0.05.
MIN_N_FOR_SIGNIFICANCE = 5


def rank(values):
    """Ranks, 1-based, averaged within ties."""
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        average = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = average
        i = j + 1
    return ranks


def spearman(x, y):
    """(rho, two-sided p), or (None, None) where the statistic is undefined.

    Undefined means one of the vectors is constant: every rank is identical, the
    denominator is zero, and no monotonic relationship can be said to exist
    either way. That is reported as missing rather than as a failure, because a
    taxon with the same count in every sample has not been tested.
    """
    n = len(x)
    if n < 3:
        return None, None
    rx, ry = rank(x), rank(y)
    mx = sum(rx) / n
    my = sum(ry) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    sxx = sum((a - mx) ** 2 for a in rx)
    syy = sum((b - my) ** 2 for b in ry)
    if sxx <= 0 or syy <= 0:
        return None, None
    rho = max(-1.0, min(1.0, sxy / math.sqrt(sxx * syy)))
    if abs(rho) >= 1.0:
        # See the module docstring: the exact floor, not the approximation's 0.
        return rho, (2.0 / math.factorial(n) if n <= 12 else 0.0)
    df = n - 2
    t = rho * math.sqrt(df / (1.0 - rho * rho))
    return rho, betainc(df / 2.0, 0.5, df / (df + t * t))


def betainc(a, b, x):
    """Regularised incomplete beta I_x(a, b), by the standard continued fraction."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    front = math.exp(
        math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
        + a * math.log(x) + b * math.log1p(-x)
    )
    if x < (a + 1.0) / (a + b + 2.0):
        return front * _betacf(a, b, x) / a
    return 1.0 - front * _betacf(b, a, 1.0 - x) / b


def _betacf(a, b, x):
    max_iterations, epsilon, tiny = 300, 3.0e-16, 1.0e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < tiny:
        d = tiny
    d = 1.0 / d
    h = d
    for m in range(1, max_iterations + 1):
        m2 = 2 * m
        for numerator in (
            m * (b - m) * x / ((qam + m2) * (a + m2)),
            -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2)),
        ):
            d = 1.0 + numerator * d
            if abs(d) < tiny:
                d = tiny
            c = 1.0 + numerator / c
            if abs(c) < tiny:
                c = tiny
            d = 1.0 / d
            h *= d * c
        # The convergence test belongs on the second half-step of the pair.
        if abs(d * c - 1.0) < epsilon:
            break
    return h


def benjamini_hochberg(pvalues):
    """BH-adjusted q-values, in the input order. None passes through as None."""
    indexed = [(i, p) for i, p in enumerate(pvalues) if p is not None]
    if not indexed:
        return list(pvalues)
    count = len(indexed)
    indexed.sort(key=lambda item: item[1])
    adjusted = list(pvalues)
    running = 1.0
    for position in range(count, 0, -1):
        index, p = indexed[position - 1]
        running = min(running, p * count / position)
        adjusted[index] = running
    return adjusted


def holm(pvalues):
    """Holm-Bonferroni adjusted p-values, in the input order. None passes through.

    SAHMI's barcode-level step calls R's `p.adjust(p)` without naming a method,
    and R's default is "holm" - family-wise error rate, not FDR. Stricter than
    BH, which matters here because the test it guards is close to a tautology.
    """
    indexed = [(i, p) for i, p in enumerate(pvalues) if p is not None]
    if not indexed:
        return list(pvalues)
    count = len(indexed)
    indexed.sort(key=lambda item: item[1])
    adjusted = list(pvalues)
    running = 0.0
    for position, (index, p) in enumerate(indexed):
        running = max(running, min(1.0, (count - position) * p))
        adjusted[index] = running
    return adjusted


# ---------------------------------------------------------------------------
# Exact tests on 2x2 tables and their meta-analysis, for the single-cell branch.
#
# CSI-Microbes (Robinson et al., Sci Adv 2024) uses Fisher's exact test rather
# than chi-square, explicitly because "for the chi-square approximation to be
# valid, the expected frequency should be at least 5" and that fails for most
# sample x cell type x genus cells. Everything here is therefore exact or
# permutation-free, and computed in log space so the factorials do not overflow.
# ---------------------------------------------------------------------------


def _log_choose(n, k):
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hypergeom_sf(k, population, successes, draws):
    """P(X >= k) for a hypergeometric draw. Upper tail, inclusive of k."""
    low = max(0, draws - (population - successes))
    high = min(draws, successes)
    if k <= low:
        return 1.0
    if k > high:
        return 0.0
    denominator = _log_choose(population, draws)
    total = 0.0
    for value in range(int(k), high + 1):
        total += math.exp(
            _log_choose(successes, value)
            + _log_choose(population - successes, draws - value)
            - denominator
        )
    return min(1.0, total)


def fisher_exact_greater(a, b, c, d):
    """One-sided Fisher exact p for enrichment in row 1 of [[a, b], [c, d]].

    Same orientation as scipy's fisher_exact(..., alternative='greater'): the
    question is whether `a` is larger than the margins would predict.
    """
    population = a + b + c + d
    successes = a + c
    draws = a + b
    return hypergeom_sf(a, population, successes, draws)


def normal_cdf(x):
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def normal_ppf(p):
    """Inverse standard normal CDF (Acklam's rational approximation, refined).

    Accurate to well under 1e-9 after one Halley step, which is far more than
    a combined p-value needs.
    """
    if p <= 0.0:
        return float("-inf")
    if p >= 1.0:
        return float("inf")
    a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
         1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
    d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
         3.754408661907416e+00]
    low, high = 0.02425, 1 - 0.02425
    if p < low:
        q = math.sqrt(-2 * math.log(p))
        x = (((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) / ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1)
    elif p <= high:
        q = p - 0.5
        r = q * q
        x = (((((a[0]*r + a[1])*r + a[2])*r + a[3])*r + a[4])*r + a[5]) * q / (((((b[0]*r + b[1])*r + b[2])*r + b[3])*r + b[4])*r + 1)
    else:
        q = math.sqrt(-2 * math.log(1 - p))
        x = -(((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) / ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1)
    # One Halley refinement against the true CDF.
    error = normal_cdf(x) - p
    density = math.exp(-x * x / 2) / math.sqrt(2 * math.pi)
    if density > 0:
        step = error / density
        x -= step / (1 + x * step / 2)
    return x


def stouffer(pvalues, weights=None):
    """Weighted Stouffer's Z combination. Returns (z, combined p).

    The 0.9999 cap is not cosmetic: scipy's combine_pvalues returns (-inf, 1)
    the moment any input p is exactly 1 (scipy issue #8506), and a p of exactly
    1 is common here because a sample with no infected cells of a type gives
    Fisher p = 1. CSI-Microbes documents the same clamp.

    Weights are the EXPECTED count of infected cells per sample, so libraries
    with almost no signal cannot dominate the combination.
    """
    usable = [(p, w) for p, w in zip(pvalues, weights or [1.0] * len(pvalues)) if p is not None]
    if not usable:
        return None, None
    numerator = 0.0
    denominator = 0.0
    for p, weight in usable:
        z = normal_ppf(1.0 - min(p, 0.9999))
        numerator += weight * z
        denominator += weight * weight
    if denominator <= 0:
        return None, None
    z = numerator / math.sqrt(denominator)
    return z, 1.0 - normal_cdf(z)


def ranksums(x, y):
    """Wilcoxon rank-sum, normal approximation with tie correction. (z, two-sided p)."""
    n1, n2 = len(x), len(y)
    if not n1 or not n2:
        return None, None
    combined = list(x) + list(y)
    ranks = rank(combined)
    total = sum(ranks[:n1])
    mean = n1 * (n1 + n2 + 1) / 2.0
    counts = {}
    for value in combined:
        counts[value] = counts.get(value, 0) + 1
    ties = sum(c ** 3 - c for c in counts.values())
    n = n1 + n2
    variance = n1 * n2 / 12.0 * ((n + 1) - ties / float(n * (n - 1))) if n > 1 else 0.0
    if variance <= 0:
        return None, None
    z = (total - mean) / math.sqrt(variance)
    return z, 2.0 * (1.0 - normal_cdf(abs(z)))
