# BUGS_SIMPLIFY — Simplify, Algebraic Manipulation & Symbolic Computation

Found **22 bugs** in the CAS simplify/factor/expand/partfrac pipeline.

---

## 🔴 Critical (Incorrect Results or Unusable Output)

### BUG #1 — Add.simplify: like terms not collected
**Location:** `expression.js`, `Add.simplify()` (~line 488)
**Description:** When both operands of Add are `Mul(Num, Sym)`, the like terms are not combined. Only the pattern `x + x → 2x` works (via toString equality), but `2x + 3x` stays unreduced.
**Correct behavior:** `2*x + 3*x` → `(5 * x)`
**Test case:**
```js
simp("2*x + 3*x")  // Got: ((2 * x) + (3 * x))
simp("x + 2*x")    // Got: (x + (2 * x))
```
**Root cause:** Add.simplify checks `l.toString() === r.toString()` to combine into `2*l`, but doesn't extract numeric coefficients from Mul terms and combine them.
**Severity:** 🔴 — Core algebraic simplification failure

---

### BUG #2 — Add.simplify: like terms with symbolic coefficients not collected
**Location:** `expression.js`, `Add.simplify()` (~line 488)
**Description:** `a*x + b*x` should become `(a+b)*x` but Add.simplify has no logic to factor out common symbolic subexpressions from Mul terms.
**Correct behavior:** `a*x + b*x` → `((a + b) * x)`
**Test case:**
```js
simp("a*x + b*x")  // Got: ((a * x) + (b * x))
```
**Severity:** 🔴 — Fundamental algebraic simplification

---

### BUG #3 — Mul.simplify: Mul(c*x^n, d*x^m) not fully combined
**Location:** `expression.js`, `Mul.simplify()` (~line 799)
**Description:** The code handles `x * x^n → x^(n+1)` and `x^n * x^m → x^(n+m)` but doesn't handle when both operands are `Mul(Num, Pow(x,n))`. The numeric coefficients get combined (via scalar associativity) but the remaining `Pow(x,a) * x` or `Mul(c, x^a) * x` doesn't reduce.
**Correct behavior:** `2*x^2 * 3*x` → `(6 * x^3)`
**Test case:**
```js
simp("2*x^2 * 3*x")  // Got: ((6 * x^2) * x)
```
**Root cause:** After scalar associativity combines 2*3=6, the result is `Mul(Mul(6, x^2), x)`. The code then checks `x * x^n` but `Mul(6, x^2)` is not a Pow, so it doesn't match.
**Severity:** 🔴 — Breaks polynomial arithmetic

---

### BUG #4 — _partfrac: _quo returns unevaluated Call for polynomial division
**Location:** `cas.js`, `_quo()` (~line 7986)
**Description:** `_quo(a, b)` only handles the `Num / Num` case. For polynomial division (e.g., `(x^2+1) / (x-1)`), it returns `new Call('quo', [a, b])` — an unevaluated symbolic node. This breaks `_partfrac()` when degree(num) ≥ degree(den).
**Correct behavior:** `partfrac((x^2+1)/(x-1), x)` → `x + 1 + 2/(x-1)`
**Test case:**
```js
partfracCmd("(x^2+1)/(x-1)", "x")    // Got: (quo((x^2 + 1), (x - 1)) + (2 / (x - 1)))
partfracCmd("(x^3+1)/(x^2-1)", "x")  // Got: (quo(((1 - x) + x^2), (x - 1)) + (1 / (x - 1)))
```
**Severity:** 🔴 — `_quo` needs full polynomial long division implementation

---

### BUG #5 — _factor: factors polynomials over ℂ instead of ℝ
**Location:** `cas.js`, `_factor()` quadratic factoring step (~line 6890-6960)
**Description:** When factoring quadratics with complex roots (discriminant < 0), `_factor` calls `_solve` which returns complex roots, then creates linear factors `(x - root)` with complex numbers. Irreducible quadratics over ℝ like `x^2+1` should remain unfactored.
**Correct behavior:**
- `factor(x^2+1)` → `x^2 + 1` (irreducible over ℝ)
- `factor(x^4-1)` → `(x-1)(x+1)(x^2+1)`
- `factor(x^3-1)` → `(x-1)(x^2+x+1)`
- `factor(x^3-8)` → `(x-2)(x^2+2x+4)`
**Test case:**
```js
factorCmd("x^2+1")   // Got: ((x + i) * (x - i))
factorCmd("x^4-1")   // Got: (((x - 1) * (x + 1)) * ((x + i) * (x - i)))
factorCmd("x^3-1")   // Got: (x-1)(complex factors)
factorCmd("x^3-8")   // Got: (x-2)(complex factors)
```
**Root cause:** The quadratic factoring step solves for all roots (including complex) and creates linear factors. It should check if roots are complex and, if so, keep the quadratic as an irreducible factor.
**Severity:** 🔴 — Wrong factorization domain

---

### BUG #6 — _factor: multivariate expressions return unevaluated factor()
**Location:** `cas.js`, `_factor()` (~line 6706)
**Description:** When `_factor` is called on a multivariate expression, it falls through to `return new Call('factor', [n])` because the polynomial coefficient extraction only supports single-variable polynomials.
**Correct behavior:** `factor(x*y)` → `x * y`, `factor(a*x+b*x)` → `x*(a+b)` (after simplification)
**Test case:**
```js
factorCmd("x*y")     // Got: factor((x * y))
factorCmd("a*x+b*x") // Got: factor(((a * x) + (b * x)))
```
**Severity:** 🔴 — Factor returns unevaluated symbolic call

---

## 🟡 Moderate (Cosmetic or Partial Failures)

### BUG #7 — Sub.simplify: commutative equality not detected
**Location:** `expression.js`, `Sub.simplify()` (~line 638)
**Description:** Sub.simplify uses `toString()` comparison to detect cancellation. Commutative expressions like `(a+b)` and `(b+a)` produce different strings, so `(a+b)-(b+a)` is not simplified to 0.
**Correct behavior:** `(a+b)-(b+a)` → `0`, `(x+1)-(1+x)` → `0`
**Test case:**
```js
simp("(a+b)-(b+a)")  // Got: ((a + b) - (b + a))
simp("(x+1)-(1+x)")  // Got: ((x + 1) - (1 + x))
```
**Root cause:** `l.toString() === r.toString()` fails for commutative reorderings. Also affects `a*b - b*a`.
**Severity:** 🟡 — Common simplification miss

---

### BUG #8 — Sub.simplify: nested Add partial cancellation fails
**Location:** `expression.js`, `Sub.simplify()` (~line 638)
**Description:** When the left side is `Add(Add(a,b), c)` (representing `a+b+c`), subtracting `Add(a,b)` fails because the code only checks direct children of Add, not grandchildren.
**Correct behavior:** `(x+y+z)-(x+y)` → `z`
**Test case:**
```js
simp("(x+y+z)-(x+y)")  // Got: ((x + (y + z)) - (x + y))
```
**Root cause:** Sub checks `l.left.toString() === r.toString()` and `l.right.toString() === r.toString()` on direct children only. The nested structure `Add(x, Add(y, z))` doesn't match `Add(x, y)`.
**Severity:** 🟡 — Common pattern in algebra

---

### BUG #9 — Mul.simplify: commutative equality not detected
**Location:** `expression.js`, `Mul.simplify()` (~line 799)
**Description:** Like Sub, Mul uses `toString()` for equality. `a*b - b*a` doesn't simplify because `(a * b)` and `(b * a)` have different toString.
**Correct behavior:** `a*b - b*a` → `0`
**Test case:**
```js
simp("a*b - b*a")  // Got: ((a * b) - (b * a))
```
**Severity:** 🟡 — Affects all commutative cancellation

---

### BUG #10 — expand: result not simplified after distribution
**Location:** `expression.js`, `Mul.expand()` and `Pow.expand()` (~line 956, ~line 1880)
**Description:** The `expand()` method distributes multiplication over addition but doesn't call `simplify()` on the result. This leaves unsimplified expressions like `x^2 - x + x - 1` instead of `x^2 - 1`.
**Correct behavior:**
- `expand((x+1)*(x-1))` → `x^2 - 1`
- `expand((x+2)*(x-3))` → `x^2 - x - 6`
**Test case:**
```js
expandCmd("(x+1)*(x-1)")  // Got: ((x^2 - x) + (x - 1))
expandCmd("(x+2)*(x-3)")  // Got: ((x^2 - (3 * x)) + ((2 * x) - 6))
expandCmd("x*(x+1)*(x-1)") // Got: ((x^3 - x^2) + (x^2 - x)) — not simplified
```
**Root cause:** `Mul.expand()` creates `Add(Mul(a,c), Mul(b,c))` etc. without simplifying each sub-term or the overall result. Each expand step is correct structurally but unsimplified.
**Severity:** 🟡 — Common user-facing issue

---

### BUG #11 — _factor: repeated factors not combined into powers
**Location:** `cas.js`, `_factor()` final assembly (~line 6950-6970)
**Description:** After finding all factors, `_factor` multiplies them together as a flat `Mul` chain. Identical factors like `(x+1) * (x+1)` are not combined into `(x+1)^2`.
**Correct behavior:**
- `factor(2*x^2+4*x+2)` → `2 * (x+1)^2`
- `factor(3*x^2+6*x+3)` → `3 * (x+1)^2`
**Test case:**
```js
factorCmd("2*x^2+4*x+2")   // Got: ((2 * (x + 1)) * (x + 1))
factorCmd("3*x^2+6*x+3")   // Got: ((3 * (x + 1)) * (x + 1))
```
**Root cause:** The assembly loop at the end of `_factor` just does `result = new Mul(result, factors[i])` without grouping identical factors.
**Severity:** 🟡 — Factorization output quality

---

### BUG #12 — Pow.simplify: negative integer exponents not converted to fractions
**Location:** `expression.js`, `Pow.simplify()` (~line 1612)
**Description:** `x^(-n)` where n is a positive integer should simplify to `1/x^n` but stays as `x^-n`.
**Correct behavior:**
- `x^(-1)` → `(1 / x)`
- `x^(-2)` → `(1 / (x^2))`
**Test case:**
```js
simp("x^(-1)")  // Got: x^-1
simp("x^(-2)")  // Got: x^-2
```
**Severity:** 🟡 — User-facing display issue

---

### BUG #13 — Pow.simplify: (c*x)^(-1) gives fractional coefficient form
**Location:** `expression.js`, `Pow.simplify()` (~line 1612)
**Description:** `(2*x)^(-1)` simplifies to `(0.5 * x^-1)` via the `(c*x)^n → c^n * x^n` rule. This is mathematically correct but should ideally produce `(1 / (2 * x))` for consistency.
**Correct behavior:** `(2*x)^(-1)` → `(1 / (2 * x))`
**Test case:**
```js
simp("(2*x)^(-1)")  // Got: (0.5 * x^-1)
```
**Severity:** 🟡 — Inconsistent form

---

### BUG #14 — (-1)*(Add) only partially distributes
**Location:** `expression.js`, `Mul.simplify()` (~line 947)
**Description:** When multiplying -1 by a nested Add like `Add(a, Add(b, c))`, the rule `-1 * (A + B) → Sub(-A, B)` only distributes to the top level, leaving `b+c` undistributed.
**Correct behavior:** `(-1)*(a+b+c)` → `(-a - b - c)`
**Test case:**
```js
simp("(-1)*(a+b+c)")  // Got: ((-1 * a) - (b + c))
```
**Root cause:** The rule `if (l instanceof Num && l.value === -1 && r instanceof Add) return new Sub(new Mul(new Num(-1), r.left).simplify(), r.right).simplify()` only negates `r.left` and keeps `r.right` as-is.
**Severity:** 🟡 — Partial distribution

---

### BUG #15 — Div.simplify: trig simplification missing
**Location:** `expression.js`, `Div.simplify()` (~line 1063)
**Description:** Div.simplify doesn't recognize trig-based simplifications like `sin/cos → tan`, `1/sin → csc`, `cos/sin → cot`, `1/cos → sec`.
**Correct behavior:**
- `sin(x)/cos(x)` → `tan(x)`
- `1/sin(x)` → `csc(x)`
- `cos(x)/sin(x)` → `cot(x)`
- `1/cos(x)` → `sec(x)`
**Test case:**
```js
simp("sin(x)/cos(x)")  // Got: (sin(x) / cos(x))
simp("1/sin(x)")       // Got: (1 / sin(x))
```
**Note:** The `Mul.simplify` handles `tan*cos → sin`, `cot*sin → cos`, etc. but Div.simplify doesn't have the reciprocal rules.
**Severity:** 🟡 — Common trig simplification

---

### BUG #16 — _factor: multivariate difference of squares produces nested factor() calls
**Location:** `cas.js`, `_factor()` (~line 6758)
**Description:** When the difference-of-squares pattern is detected but involves a multivariate expression, the recursive `_factor` calls return unevaluated `Call('factor', ...)` nodes.
**Correct behavior:** `factor(x^2-(y+1)^2)` → `(x-y-1)*(x+y+1)`
**Test case:**
```js
factorCmd("x^2-(y+1)^2")  // Got: (factor((x - (y + 1))) * factor((x + (y + 1))))
```
**Severity:** 🟡 — Ugly output for valid algebraic identities

---

### BUG #17 — Div.simplify: multivariate polynomial GCD not supported
**Location:** `expression.js`, `Div.simplify()` polynomial GCD section (~line 1240)
**Description:** The polynomial GCD simplification in Div only supports single-variable polynomials. Multivariate expressions like `(x^2-y^2)/(x-y)` or `(x^3-y^3)/(x-y)` are not simplified.
**Correct behavior:**
- `(x^2-y^2)/(x-y)` → `(x + y)`
- `(x^3-y^3)/(x-y)` → `(x^2 + x*y + y^2)`
**Test case:**
```js
simp("(x^2-y^2)/(x-y)")  // Got: ((x^2 - y^2) / (x - y))
simp("(x^3-y^3)/(x-y)")  // Got: ((x^3 - y^3) / (x - y))
```
**Severity:** 🟡 — Missing feature for multivariate algebra

---

### BUG #18 — expand: Mul terms within expansion not simplified
**Location:** `expression.js`, `Pow.expand()` binomial expansion (~line 1880)
**Description:** The binomial expansion creates `Mul(Num(2), Mul(Sym(x), Sym(x)))` which doesn't get simplified to `Mul(Num(2), Pow(Sym(x), 2))` because expand doesn't call simplify on sub-terms.
**Correct behavior:** `expand((x+1)^2*(x-1))` should simplify internally
**Test case:**
```js
expandCmd("(x+1)^2*(x-1)")  // Got: ((x^3 - x^2) + ((((2 * x) * x) - (2 * x)) + (x - 1)))
// (2 * x) * x should be 2 * x^2
```
**Severity:** 🟡 — Ugly nested Mul in expanded output

---

## 🟢 Minor (Nice-to-Have Improvements)

### BUG #19 — Add.simplify: common factor not extracted
**Location:** `expression.js`, `Add.simplify()` (~line 488)
**Description:** `a*x + a*y` should factor to `a*(x+y)`. The Add.simplify doesn't have logic to extract common factors from Mul terms.
**Correct behavior:** `a*x + a*y` → `(a * (x + y))`
**Test case:**
```js
simp("a*x + a*y")  // Got: ((a * x) + (a * y))
```
**Severity:** 🟢 — Advanced factoring (CAS nice-to-have)

---

### BUG #20 — _factor: integer factoring display uses non-standard 'factored'/'power' wrappers
**Location:** `cas.js`, `_factor()` integer branch (~line 6738)
**Description:** Integer factorization returns `Call('factored', [...])` and `Call('power', [...])` which are non-standard display nodes. `factor(-12)` → `(-1 * factored(power(2, 2), 3))`.
**Correct behavior:** Should use standard math notation like `-1 * 2^2 * 3`
**Test case:**
```js
factorCmd("-12")  // Got: (-1 * factored(power(2, 2), 3))
```
**Severity:** 🟢 — Display convention

---

### BUG #21 — expand: expand of (a+b+c)^2 works but is deeply nested
**Location:** `expression.js`, `Pow.expand()` binomial expansion
**Description:** `(a+b+c)^2` expands correctly but produces deeply nested Add/Mul structure. Due to parser associativity, `a+b+c` is `Add(a, Add(b, c))`, and the binomial expansion handles it via recursive expand. The result is correct but hard to read.
**Correct behavior:** Clean output like `a^2 + b^2 + c^2 + 2ab + 2ac + 2bc`
**Test case:**
```js
expandCmd("(a+b+c)^2")  // Got: (a^2 + ((2 * (a * b)) + ((2 * (a * c)) + (b^2 + ((2 * (b * c)) + c^2)))))
```
**Severity:** 🟢 — Cosmetic

---

### BUG #22 — partfrac: non-monic denominator not handled cleanly
**Location:** `cas.js`, `_partfrac()` (~line 9961)
**Description:** `partfrac(1/(2*x+1), x)` gives `(1 / (2 * (x - (-1/2))))` instead of decomposing properly. The root-finding works but the coefficient extraction doesn't account for the leading coefficient of the denominator factor.
**Correct behavior:** `1/(2x+1)` → `1/(2x+1)` (single irreducible factor) or `(1/2) * 1/(x+1/2)`
**Test case:**
```js
partfracCmd("1/(2*x+1)", "x")  // Got: (1 / (2 * (x - (-1 / 2))))
```
**Severity:** 🟢 — Non-standard form

---

## Summary

| Severity | Count | Key Areas |
|----------|-------|-----------|
| 🔴 Critical | 6 | Like-term collection, Mul coefficient merging, polynomial division (_quo), complex factorization, multivariate factor |
| 🟡 Moderate | 12 | Commutative equality, expand simplification, negative powers, trig ratios, partial distribution, multivariate GCD |
| 🟢 Minor | 4 | Common factor extraction, integer display, nesting, non-monic partfrac |
| **Total** | **22** | |
