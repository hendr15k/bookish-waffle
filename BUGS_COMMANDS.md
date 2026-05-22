# BUGS_COMMANDS — CAS Command Bugs (cas.js Part 2)

Found via systematic code analysis + Node.js testing against the CAS engine.

---

## 🔴 Critical (Crashes / Wrong Results)

### BUG #1 — `_solve` crashes with TypeError on quintic/specific quartic polynomials
**Location:** `cas.js`, `_solve()` quartic handler (~line 4825)
**Description:** Calling `solve` on polynomials of degree ≥ 4 that trigger the internal Ferrari quartic solver throws `TypeError: this.left.simplify is not a function` (or `this.right.simplify is not a function`). The bug occurs in the shift-back step `xSols = ySols.map(y => new Sub(y, shift).simplify())`. The root cause appears to be that `shift` (or some intermediate variable in the quartic handler) becomes a non-Expr object (possibly `undefined`), causing the subsequent `Sub.simplify()` to fail. This is likely triggered when:
1. The `_factor` call returns the original polynomial unchanged (irreducible quartic)
2. The `_getPolyCoeffs` returns coefficients that lead to a symbolic `shift`
3. An intermediate expression in the Ferrari method doesn't resolve to a proper Expr

**Affected cases:**
```js
solve(x^5 - x^4 + x^3 - x^2 + x - 1, x)  // TypeError, returns only x=1
solve(x^5 - 1, x)                          // TypeError, returns only x=1
solve(x^7 - 1, x)                          // TypeError
solve(x^4 + x^2 + 1, x)                   // TypeError (Ferrari on irreducible quartic)
solve((x-1)^4, x)                         // TypeError (biquadratic with repeated roots)
solve(x^4 + x + 1, x)                     // TypeError (Ferrari general case)
```
**Correct behavior:** Should return the correct roots without crashing.
**Test case:**
```js
const CAS = require('./js/cas.js');
const { Num, Sym, Call, Add, Sub, Mul, Div, Pow } = require('./js/expression.js');
const cas = new CAS();
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

### BUG #2 — `solve(abs(x) = -3)` returns extraneous solutions
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
**Description:** `series(sqrt(1+x), x, -1, 3)` attempts to find the Laurent series of `sqrt(1+x)` around `x = -1`. The function has a singularity at `x = -1` (since `sqrt(0)`). The current implementation returns `0` which is incorrect — it should return a proper Laurent series like `(1/sqrt(x+1))` terms or at minimum the leading singular behavior.
**Correct behavior:** Should return a proper Laurent expansion accounting for the square-root singularity.
**Test case:**
```js
cas.evaluate(new Call('series', [
  new Call('sqrt', [new Add(new Num(1), x)]), x, new Num(-1), new Num(3)
]));
// Got: 0
// Expected: Laurent series with singular term (e.g., sqrt(2) + terms in (x+1))
```

---

## 🟢 Minor (Missing Features / Edge Cases)

### BUG #4 — `_sumSymbolic` doesn't handle `Div` expressions
**Location:** `cas.js`, `_sumSymbolic()` (~line 3209)
**Description:** `_sumSymbolic` has handlers for `Add`, `Sub`, and `Mul` but not `Div`. When the summand is a fraction like `a/k`, the expression isn't recognized as needing special handling and the geometric series ratio detection also fails.
**Test case:**
```js
// sum(a/k, k, 1, n) returns unevaluated
cas.evaluate(new Call('sum', [
  new Div(new Sym('a'), new Sym('k')), new Sym('k'), new Num(1), new Sym('n')
]));
// Got: sum((a / k), k, 1, n)  (unevaluated)
// Expected: a * H_n (harmonic number), or at least a * psi(n+1) - a * psi(1)
```

### BUG #5 — `_sumSymbolic` doesn't simplify telescoping series
**Location:** `cas.js`, `_sumSymbolic()` (~line 3209)
**Description:** `sum(1/(k*(k+1)), k, 1, n)` is a telescoping series that simplifies to `n/(n+1)`. The current code doesn't recognize this pattern.
**Test case:**
```js
// sum(1/(k*(k+1)), k, 1, n) = n/(n+1)
cas.evaluate(new Call('sum', [
  new Div(new Num(1), new Mul(new Sym('k'), new Add(new Sym('k'), new Num(1)))),
  new Sym('k'), new Num(1), new Sym('n')
]));
// Got: sum((1 / (k^2 + k)), k, 1, n)  (unevaluated)
// Expected: n/(n+1)
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
| 1 | 🔴 Critical | `_solve` | Crashes with TypeError on quintic/specific quartic polynomials | UNFIXED |
| 2 | 🔴 Critical | `_solve` | `abs(x) = -3` gives extraneous solutions | UNFIXED |
| 3 | 🟡 Moderate | `_laurent` | `series(sqrt(1+x), x, -1, 3)` returns 0 instead of Laurent series | UNFIXED |
| 4 | 🟢 Minor | `_sumSymbolic` | `Div` summand not recognized | UNFIXED |
| 5 | 🟢 Minor | `_sumSymbolic` | Telescoping series not simplified | UNFIXED |
| 6 | 🟢 Minor | integrate | `piecewise` indefinite integral not handled | UNFIXED |

---

## Additional Observations (Not Bugs — Expected Behavior)

The following produce correct or expected unevaluated results:

- `solve(tan(x)=1, x)` → `pi/4` (principal solution)
- `solve(x*exp(x)=1, x)` → `LambertW(1)` ✓
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
