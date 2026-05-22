# BUGS_COMMANDS — CAS Command Bugs (cas.js Part 2)

Found via systematic code analysis + Node.js testing against the CAS engine.

---

## 🔴 Critical (Crashes / Wrong Results)

### BUG #1 — `_solve` crashes with TypeError on quintic/specific quartic polynomials; incomplete roots for degree ≥5
**Location:** `cas.js`, `_solve()` quartic handler (~line 4825) and `_durandKerner` (~line 8895)
**Description:** (1) `_durandKerner` always returned real numbers even for complex roots, causing `solve(x^n-1, x)` for n ≥ 5 to return only the real roots. (2) The quartic Ferrari solver's shift-back step could fail with `this.left.simplify is not a function` when `ySols` contained non-Expr objects.

**Status: FIXED**
- `_durandKerner` now returns `{re, im}` objects for complex roots, preserving imaginary parts
- `_solvePolynomialNumerically` updated to handle the new complex root format
- Filter added to quartic shift-back: `.filter(y => y instanceof Expr)`
- Ferrari failures now fall back to numerical companion-matrix method

**Affected cases (now fixed):**
```js
solve(x^5 - 1, x)   // was returning 5 real roots, now 5 correct (incl. complex)
solve(x^7 - 1, x)   // was returning only 5 roots, now returns all 7
solve(x^4 + x^2 + 1, x)  // was crashing, now returns 4 roots
```
**Correct behavior:** Should return the correct roots without crashing.

### BUG #1b — `solve(abs(x) = -3)` returns extraneous solutions
**Status: FIXED** — abs handler now checks RHS >= 0 and returns empty set otherwise.
const x = new Sym('x');
try {
  cas.evaluate(new Call('solve', [
    new Sub(new Sub(new Sub(new Sub(new Pow(x, new Num(5)), new Pow(x, new Num(4))),
      new Pow(x, new Num(3)), new Pow(x, new Num(2))), new Add(x, new Num(-1))
  ]), x]));
} catch(e) {
  console.log(e.message); // "this.left.simplify is not a function"
}
```

---

### BUG #2 — `solve(abs(x) = -3)` returns extraneous solutions (NOW FIXED)
**Location:** `cas.js`, `_solve()` abs handler (~line 4299)
**Description:** When solving `abs(x) = -3` (which has no solutions since `abs` is always non-negative), the CAS returns `{-3, 3}` — treating the equation as if `abs(x) = 3`. The code generates both `A = B` and `A = -B` but doesn't check that `B >= 0` for real solutions.
**Correct behavior:** `solve(abs(x) = -3, x)` → `{}` (empty set) or `false` or a single valid root if any exist.
**Test case:**
```js
cas.evaluate(new Call('solve', [
  new Eq(new Call('abs', [x]), new Num(-3)), x
]));
// Got: {-3, 3}
// Expected: {} (no solutions — abs(x) cannot equal a negative number)
```

---

## 🟡 Moderate (Incorrect Output)

### BUG #3 — `series(sqrt(1+x), x, -1, 3)` returns `0` instead of Laurent series
**Location:** `cas.js`, `_laurent()` / `series` handler
**Description:** `series(sqrt(1+x), x, -1, 3)` attempts to find the Laurent series of `sqrt(1+x)` around `x = -1`. The function has a branch-point singularity at `x = -1` (since `sqrt(0)`). The `_laurent` method detects the pole order `k=1` correctly but `_taylor` on `g = sqrt(1+x)*(x+1) = (x+1)^(3/2)` returns `0` because all derivatives at the branch point are `0` or `Infinity`. Dividing `0/(x+1) = 0`.
**Status: FIXED** — Added branch-point detection in `_laurent`: when `gSeries` is `0` at the base point, check if `f^2/(x-a)` simplifies to a nonzero constant. If so, `f` contains a `sqrt(x-a)` factor and is returned as the leading singular term.
**Test case:**
```js
cas.evaluate(new Call('series', [
  new Call('sqrt', [new Add(new Num(1), x)]), x, new Num(-1), new Num(3)
]));
// Got: sqrt(1+x) — leading singular term (FIXED)
// Previously: 0
```

---

## 🟢 Minor (Missing Features / Edge Cases)

### BUG #4 — `_sumSymbolic` doesn't handle `Div` expressions
**Location:** `cas.js`, `_sumSymbolic()` (~line 3209)
**Description:** `_sumSymbolic` has handlers for `Add`, `Sub`, and `Mul` but not `Div`. When the summand is a fraction like `a/k`, the expression isn't recognized as needing special handling and the geometric series ratio detection also fails.
**Status: FIXED** — Added Div handler to factor out constant numerators (e.g. `sum(a/k, k, 1, n)` → `a * sum(1/k, k, 1, n)`) and constant denominators.
**Test case:**
```js
// sum(a/k, k, 1, n) returns unevaluated
cas.evaluate(new Call('sum', [
  new Div(new Sym('a'), new Sym('k')), new Sym('k'), new Num(1), new Sym('n')
]));
// Was: sum((a / k), k, 1, n)  (unevaluated)
// Now: a * sum((1 / k), k, 1, n)  (constant factored out)
```

### BUG #5 — `_sumSymbolic` doesn't simplify telescoping series
**Location:** `cas.js`, `_sumSymbolic()` (~line 3209)
**Description:** `sum(1/(k*(k+1)), k, 1, n)` is a telescoping series that simplifies to `n/(n+1)`. The current code doesn't recognize this pattern.
**Status: FIXED** — Added telescoping series detection via quadratic factorization and partial fractions. Handles `sum(c/((k+a)*(k+b)), k, 1, n)` by decomposing to `(c/(b-a)) * (1/(1+a) - 1/(n+b))`.
**Test case:**
```js
// sum(1/(k*(k+1)), k, 1, n) = n/(n+1)
cas.evaluate(new Call('sum', [
  new Div(new Num(1), new Mul(new Sym('k'), new Add(new Sym('k'), new Num(1)))),
  new Sym('k'), new Num(1), new Sym('n')
]));
// Was: sum((1 / (k^2 + k)), k, 1, n)  (unevaluated)
// Now: (1 - (1 / (n + 1)))  which equals n/(n+1)
```

### BUG #6 — `_integrate` doesn't handle `piecewise` expressions symbolically
**Location:** `cas.js`, integrate handler
**Description:** `integrate(piecewise(x<0, x, x^2), x)` returns unevaluated, even though the indefinite integral can be computed piecewise: `-x^2/2` for `x<0`, `x^3/3` for `x≥0`.
**Correct behavior:** Should return `piecewise(x<0, -x^2/2, x^3/3)` or similar.
**Test case:**
```js
// integrate(piecewise(x<0, x, x^2), x) — returns unevaluated
// Expected: piecewise(x<0, -x^2/2, x^3/3)
```
Note: The definite integral `integrate(piecewise(x<0, 0, x^2), x, -1, 1)` correctly returns `1/3`.

---

## Summary

| # | Severity | Component | Summary | Status |
|---|----------|-----------|---------|--------|
| 1 | 🔴 Critical | `_solve` | Crashes with TypeError on quintic/specific quartic polynomials; incomplete roots for degree ≥5 | **FIXED** — filter added + Durand-Kerner returns complex roots |
| 2 | 🔴 Critical | `_solve` | `abs(x) = -3` gives extraneous solutions | **FIXED** — returns empty set when RHS < 0 |
| 3 | 🟡 Moderate | `_laurent` | `series(sqrt(1+x), x, -1, 3)` returns 0 instead of Laurent series | **FIXED** — branch-point detection returns leading singular term `sqrt(1+x)` |
| 4 | 🟢 Minor | `_sumSymbolic` | `Div` summand not recognized | **FIXED** — Div handler factors out constant numerators and denominators |
| 5 | 🟢 Minor | `_sumSymbolic` | Telescoping series not simplified | **FIXED** — Telescoping series via quadratic factorization and partial fractions |
| 6 | 🟢 Minor | integrate | `piecewise` indefinite integral not handled | UNFIXED |

---

## Additional Observations (Not Bugs — Expected Behavior)

The following produce correct or expected unevaluated results:

- `solve(tan(x)=1, x)` → `pi/4` (principal solution) ✓ **FIXED HEN-207** — was crashing before null-check fix
- `solve(x*exp(x)=1, x)` → `LambertW(1)` ✓ **FIXED HEN-207** — was crashing before null-check fix
- `limit((1+1/x)^x, x, Infinity)` → `e` ✓
- `limit(x^x, x, 0, right)` → `1` ✓
- `taylor(sin(x)/x, x, 0, 4)` → `1` ✓ (higher-order terms are 0)
- `product(k, k, 1, n)` → `gamma(n+1)` ✓
- `taylor(exp(x), x, 0, 10)` → correct series ✓
- `implicitDiff(y^2=x^3, y, x)` → `3x^2/(2y)` ✓
- `integrate(sinh(x), x)` → `cosh(x)` ✓
- `integrate(cosh(x), x)` → `sinh(x)` ✓
- `diff(piecewise(x<0, x, x^2), x)` → `piecewise(x<0, 1, 2x)` ✓
- `solve(abs(x)=2, x)` → `{-2, 2}` ✓
- `solve(x^4=16, x)` → all four roots (complex) ✓
- `limit(abs(x)/x, x, 0, right)` → `1` ✓
- `limit(abs(x)/x, x, 0, left)` → `-1` ✓
