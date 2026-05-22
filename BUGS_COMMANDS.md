# BUGS_COMMANDS — cas.js Command Bugs

Found via Node.js testing against the CAS. All tests run from workspace root with correct module loading (`const { CAS } = require('./js/cas.js')`).

---

## 🔴 Critical

### BUG #1 — `solve` returns incomplete root set for degree ≥ 5 polynomials
**Location:** `cas.js`, `_solve()` polynomial branch (~line 4560-4840)
**Description:** `solve(x^5-1, x)` returns `1` instead of all 5 roots. `solve(x^5+1, x)` returns `-1` instead of all 5 roots. Only one root is returned for quintic equations. Similarly, `solve(x^3-3*x+2, x)` returns `{-2, 1}` instead of `{-2, 1, 1}` (double root `1` is missing).

**Root cause:** The polynomial solver uses numerical root-finding for degree ≥ 5 (Abel-Ruffini: no general formula). The numerical method likely fails to converge or misses roots for these specific polynomials. Also, duplicate root detection may be removing valid repeated roots.

**Correct behavior:**
- `solve(x^5-1, x)` → all 5 fifth roots of unity
- `solve(x^5+1, x)` → all 5 fifth roots of -1
- `solve(x^3-3*x+2, x)` → `{1, 1, -2}` (or `{-2, 1}` with deduplication but correct count)

**Test case:**
```js
const {Parser, Lexer} = require('./js/parser.js');
const { CAS } = require('./js/cas.js');
const cas = new CAS();
function parse(s) { return new Parser(new Lexer(s)).parse(); }
cas.evaluate(parse('solve(x^5-1,x)')).toString()
// Got: 1
// Expected: set of 5 complex roots

cas.evaluate(parse('solve(x^5+1,x)')).toString()
// Got: -1
// Expected: set of 5 complex roots

cas.evaluate(parse('solve(x^3-3*x+2,x)')).toString()
// Got: {-2, 1}
// Expected: {-2, 1, 1} (x^3-3*x+2 = (x+2)(x-1)^2)
```

**Severity:** 🔴 Critical — wrong mathematical results

---

### BUG #2 — `Sub.simplify` crashes with "this.left.simplify is not a function"
**Location:** `cas.js:4825` (inside `_solve()`) calling `Sub.simplify`
**Description:** When `_solve` processes certain polynomial equations (e.g., `solve(x^5+1,x)`, `solve(x^4+1,x)`), it creates `Sub` nodes where `this.left` is not an `Expr` object (likely `null` or `undefined`). Calling `.simplify()` on such a `Sub` crashes with `this.left.simplify is not a function`.

**Root cause:** `_solve` creates `Sub` nodes at line 4825 without ensuring both operands are valid `Expr` instances. The code `const newSub = new Sub(expr, ...)` is called with something that's not an `Expr` on the left side.

**Correct behavior:** `solve` should return valid root sets or handle edge cases gracefully.

**Test case:**
```js
cas.evaluate(parse('solve(x^5+1,x)')).toString()
// Got: TypeError: this.left.simplify is not a function
// Expected: set of 5 roots (or error message)
```

**Severity:** 🔴 Critical — crash, not graceful degradation

---

## 🟡 Moderate

### BUG #3 — `eval` doesn't evaluate numeric expressions
**Location:** `cas.js`, `_eval` command handler
**Description:** `eval` doesn't evaluate numeric expressions like `1/3`, `2^10`, or symbolic expressions with constants (`e`, `pi`). It just returns the expression unchanged.

**Correct behavior:**
- `eval(1/3)` → `0.33333...` or `0.3333`
- `eval(2^10)` → `1024`
- `eval(e)` → `2.71828...` (numeric approximation)
- `eval(pi)` → `3.14159...`

**Test case:**
```js
cas.evaluate(parse('eval(1/3,5)')).toString()
// Got: eval((1 / 3), 5)

cas.evaluate(parse('eval(2^10)')).toString()
// Got: eval(1024)

cas.evaluate(parse('eval(exp(1),5)')).toString()
// Got: eval(e, 5)
```

**Severity:** 🟡 Moderate — eval command barely functional

---

### BUG #4 — `series(log(1+x), x, 0, n)` produces incorrect output
**Location:** `cas.js`, `_series` command handler
**Description:** The Maclaurin series for `log(1+x)` is `x - x²/2 + x³/3 - x⁴/4 + ...`. The command produces deeply nested, incorrect expressions mixing `ln` and powers in wrong ways.

**Correct behavior:**
- `series(log(1+x), x, 0, 5)` → `x - x²/2 + x³/3 - x⁴/4 + x⁵/5`

**Test case:**
```js
cas.evaluate(parse('series(log(1+x),x,0,5)')).toString()
// Got: deeply nested incorrect expression with ln(100), ln(1000), etc.
// Expected: x - x²/2 + x³/3 - x⁴/4 + x⁵/5
```

**Note:** `series(exp(x)-1,x,0,5)` works correctly.

**Severity:** 🟡 Moderate — wrong mathematical result

---

### BUG #5 — `subst` command not supported (3-arg syntax)
**Location:** `cas.js`, CAS command handlers (~line 125)
**Description:** The CAS only handles `subs` and `substitute` commands, not `subst`. All three are standard aliases in CAS systems (Maple, Mathematica, etc.). Users calling `subst(expr, var, value)` get an unevaluated `subst(...)` node instead of the substituted result.

**Correct behavior:**
- `subst(x^2, x, 2)` → `4`
- `subs(x^2, x=2)` → `4` (current, works)

**Test case:**
```js
cas.evaluate(parse('subst(x^2,x,2)')).toString()
// Got: subst(x^2, x, 2)  ← unevaluated

cas.evaluate(parse('subs(x^2,x=2)')).toString()
// Got: 4  ← works with correct syntax
```

**Severity:** 🟡 Moderate — missing user-facing command

---

### BUG #6 — `simplify` doesn't simplify `abs(x^2)` → `x^2`
**Location:** `expression.js` (`Call.simplify()`)
**Description:** `abs(x)` of a non-negative expression should simplify. Since `x^2 ≥ 0` for all real `x`, `abs(x^2) = x^2`. Similarly `abs(x)^2 = x^2`.

**Correct behavior:**
- `simplify(abs(x^2))` → `x^2`
- `simplify(abs(x)^2)` → `x^2`
- `simplify(abs(x*y))` → `abs(x) * abs(y)` (optional, lower priority)

**Test case:**
```js
cas.evaluate(parse('simplify(abs(x^2))')).toString()
// Got: abs(x^2)

cas.evaluate(parse('simplify(abs(x)^2)')).toString()
// Got: abs(x)^2
```

**Severity:** 🟡 Moderate — incomplete simplification for known non-negative forms

---

### BUG #7 — `integrate(x^(-1/2), x)` not simplified
**Location:** `expression.js` (`Call.simplify()` numeric evaluation) + `cas.js` integrate handler
**Description:** `x^(-1/2)` is equivalent to `1/sqrt(x)`, and `integrate(1/sqrt(x), x) = 2*sqrt(x)`. But `integrate(x^(-1/2), x)` returns the unevaluated call instead of the same result.

**Correct behavior:** `integrate(x^(-1/2), x)` → `2 * sqrt(x)`

**Test case:**
```js
cas.evaluate(parse('integrate(x^(-1/2),x)')).toString()
// Got: integrate(x^((-1 / 2)), x)  ← unevaluated

cas.evaluate(parse('integrate(1/sqrt(x),x)')).toString()
// Got: (2 * sqrt(x))  ← correct
```

**Severity:** 🟡 Moderate — missing simplification before integration

---

### BUG #8 — `factor(x^4+1)` not factored
**Location:** `cas.js`, `_factor()` quadratic branch
**Description:** `x^4+1` is irreducible over ℝ but reducible over ℂ into `(x^2 + sqrt(2)*x + 1)*(x^2 - sqrt(2)*x + 1)`. Currently returns `x^4 + 1` unchanged.

**Correct behavior:**
- `factor(x^4+1)` → product of quadratics over ℝ (with sqrt(2))
- `factor(x^4+1)` over ℂ → four linear factors

**Test case:**
```js
cas.evaluate(parse('factor(x^4+1)')).toString()
// Got: (x^4 + 1)
```

**Severity:** 🟡 Moderate — incomplete factorization

---

## Summary

| # | Severity | Component | Summary | Status |
|---|----------|-----------|---------|--------|
| 1 | 🔴 Critical | `_solve` | Returns incomplete root set for degree ≥ 5 | Open |
| 2 | 🔴 Critical | `_solve` / `Sub.simplify` | Crash on certain polynomial solves | Open |
| 3 | 🟡 Moderate | `_eval` | eval doesn't evaluate numeric expressions | Open |
| 4 | 🟡 Moderate | `_series` | log(1+x) series gives wrong output | Open |
| 5 | 🟡 Moderate | CAS handlers | `subst` command not supported | Open |
| 6 | 🟡 Moderate | `Call.simplify()` | `abs(x^2)` not simplified to `x^2` | Open |
| 7 | 🟡 Moderate | `integrate` | `x^(-1/2)` not simplified before integration | Open |
| 8 | 🟡 Moderate | `_factor` | `x^4+1` not factored | Open |
