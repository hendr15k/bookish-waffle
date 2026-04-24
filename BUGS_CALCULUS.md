# Calculus Engine Bug Report

## Bug #1: `_partfrac()` treats transcendental denominators as polynomials
**Location:** `cas.js:_partfrac()` (~line 9961)
**Severity:** 🔴 Critical
**Description:** `_partfrac()` calls `_solve(den, varNode)` to find roots of the denominator, but never checks whether the denominator is actually a polynomial in `varNode`. For transcendental expressions like `e^x - 1`, it finds the root `x = 0` and constructs `1/x` as the partial fraction — which is completely wrong. Similarly for `e^x + 1`, it finds `x = ln(-1)` and produces `(-1)/(x - ln(-1))`.

This bug propagates to `integrate()`, which uses `_partfrac` as a fallback strategy, producing wrong antiderivatives.

**Correct behavior:** `_partfrac` should verify the denominator is a rational function (ratio of polynomials) before decomposing. If not, it should return the expression unchanged.

**Test cases:**
```js
// BUG: returns 1/x, should return 1/(e^x-1) unchanged
partfrac(1/(e^x-1), x)

// BUG: returns (-1/(x-ln(-1))), should return 1/(e^x+1) unchanged
partfrac(1/(e^x+1), x)

// BUG: integrate returns ln(abs(x)) — the antiderivative of 1/x, not 1/(e^x-1)
integrate(1/(e^x-1), x)
// Expected: ln(|1-e^(-x)|) or x-ln(e^x-1)

// BUG: integrate returns ln(1/(x-ln(-1))) — nonsensical
integrate(1/(e^x+1), x)
// Expected: x - ln(e^x+1) or -ln(e^(-x)+1)
```

---

## Bug #2: `_partfrac()` decomposes irreducible quadratics over ℂ instead of ℝ
**Location:** `cas.js:_partfrac()` (~line 9961)
**Severity:** 🔴 Critical
**Description:** `_partfrac` finds complex roots of irreducible quadratics (e.g., `x²+1` has roots `±i`) and decomposes over ℂ. This produces complex partial fractions like `(-1/(2i(x+i))) + (-i/(2(x-i)))` instead of keeping `1/(x²+1)` as-is. When integrated, this yields complex logarithms instead of `atan(x)`. The complex result also evaluates to `NaN` numerically.

Additionally, `partfrac(1/(x²+2), x)` returns `0` entirely — a complete failure where the complex residue computation yields zero.

**Correct behavior:** Over ℝ, irreducible quadratic factors should be kept intact. Partial fractions should produce terms of the form `(Ax+B)/(x²+bx+c)` for irreducible quadratics, not complex linear factors.

**Test cases:**
```js
// BUG: returns complex decomposition, should return 1/(x²+1)
partfrac(1/(x²+1), x)

// BUG: returns 0, should return 1/(x²+2)
partfrac(1/(x²+2), x)

// BUG: returns complex logs instead of atan(x)
integrate(x²/(x²+1), x)
// Expected: x - atan(x)
```

---

## Bug #3: `_integrateByParts()` infinite recursion on cyclic by-parts integrals
**Location:** `cas.js:_integrateByParts()` (~line 5675)
**Severity:** 🔴 Critical
**Description:** Integration by parts creates infinite recursion (stack overflow) for integrals where applying IBP twice returns to the original integral. The code calls `this.evaluate(new Call("integrate", [vdu, varNode]))` recursively without detecting cycles.

The canonical examples are `∫ eˣsin(x) dx` and `∫ sin(x)cos(x) dx`. For these, applying IBP produces an equation of the form `I = uv - I`, which should be solved algebraically as `I = uv/2`. Instead, the code recurses infinitely.

**Correct behavior:** After one round of IBP, check if the remaining integral is proportional to the original. If so, solve the algebraic equation `I = uv - kI` → `I = uv/(1+k)`.

**Test cases:**
```js
// BUG: Maximum call stack size exceeded
integrate(e^x*sin(x), x)
// Expected: (e^x*(sin(x)-cos(x)))/2

// BUG: Maximum call stack size exceeded
integrate(sin(x)*cos(x), x)
// Expected: sin(x)^2/2 or -cos(2x)/4
```

---

## Bug #4: `integrate(sin(x)², x)` fails — produces self-referential result
**Location:** `cas.js` integrate handler (~line 695-705) and `expression.js:Mul.integrate` / by-parts cycle
**Severity:** 🔴 Critical
**Description:** `integrate(sin(x)², x)` returns `((-1 * integrate(cos(x)², x)) + (integrate(sin(x)², x) + x)) / 2` — a self-referential expression that contains the original integral on both sides. The issue is that IBP on `sin²` creates a cycle (`sin² = sin·sin`), and the Pythagorean identity substitution `sin² = 1 - cos²` just shifts the problem to `cos²`. The `_linearizeTrig` function correctly produces `(-cos(2x)+1)/2`, but by the time it's tried, the earlier strategies have already produced the broken result.

**Correct behavior:** Should use power reduction: `sin(x)² = (1-cos(2x))/2`, integrate to `x/2 - sin(2x)/4`.

**Test cases:**
```js
// BUG: self-referential result
integrate(sin(x)^2, x)
// Expected: x/2 - sin(2*x)/4

// BUG: self-referential result
integrate(cos(x)^2, x)
// Expected: x/2 + sin(2*x)/4
```

---

## Bug #5: `_limit()` doesn't recognize `Sym('inf')` as Infinity
**Location:** `cas.js:_limit()` (~line 5783+)
**Severity:** 🔴 Critical
**Description:** The parser creates `Sym('inf')` for the token `inf`, but all `isInf()` helper functions inside `_limit()` check for `name === 'Infinity' || name === 'infinity'`. Since `'inf' !== 'Infinity'`, infinity is never recognized in limit calculations when the user types `inf`.

Using `Infinity` instead works correctly, but `inf` is the natural input and is what most users would type.

**Correct behavior:** `isInf()` checks should also match `name === 'inf'`. Alternatively, the parser should map `inf` to `Sym('Infinity')`.

**Test cases:**
```js
// BUG: returns inf² (should be inf)
limit(x², x, inf)

// BUG: returns 1/inf (should be 0)
limit(1/x, x, inf)

// BUG: returns e^inf (should be inf)
limit(e^x, x, inf)

// BUG: returns ((1+(1/inf)))^inf (should be e)
limit((1+1/x)^x, x, inf)

// OK with Infinity:
limit(x², x, Infinity)  // → Infinity ✓
limit(1/x, x, Infinity)  // → 0 ✓
```

---

## Bug #6: `_limit()` doesn't handle `∞/∞` L'Hôpital when using `inf` Sym
**Location:** `cas.js:_limit()` Div handler
**Severity:** 🟡 Major
**Description:** Related to Bug #5 — since `Sym('inf')` is not recognized as infinity, the `isInfinite()` check fails and L'Hôpital's rule is never triggered for `∞/∞` forms when the point is `inf`. With `Infinity` it works.

**Test cases:**
```js
// BUG: returns ln(inf)/inf (should be 0 via L'Hôpital)
limit(ln(x)/x, x, inf)

// BUG: returns inf/e^inf (should be 0)
limit(x/e^x, x, inf)

// These work with Infinity:
limit(ln(x)/x, x, Infinity)  // → 0 ✓
limit(x/e^x, x, Infinity)    // → 0 ✓
```

---

## Bug #7: `_limit()` doesn't detect oscillating limits at infinity
**Location:** `cas.js:_limit()` fallback path
**Severity:** 🟡 Major
**Description:** `limit(sin(x), x, Infinity)` returns `sin(Infinity)` instead of `undefined`. For bounded oscillating functions like `sin(x)` and `cos(x)`, the limit at infinity does not exist. The code has no detection for this case — it just substitutes and returns the unevaluated expression.

**Correct behavior:** If a function is bounded (like sin/cos with range [-1,1]) and the point is ±∞, and direct substitution gives an expression involving ∞ inside a periodic function, return `undefined` (or `NaN`/`"undefined"`).

**Test cases:**
```js
// BUG: returns sin(Infinity), should be undefined
limit(sin(x), x, Infinity)

// BUG: returns cos(Infinity), should be undefined
limit(cos(x), x, Infinity)
```

---

## Bug #8: `_limit()` doesn't handle `0·∞` form with `ln`
**Location:** `cas.js:_limit()` Mul handler
**Severity:** 🟡 Major
**Description:** `limit(x*ln(x), x, 0, 1)` returns `limit((x * ln(x)), x, 0)` unevaluated. This is a `0·(-∞)` indeterminate form. The Mul handler checks for `NaN` from `0 * Infinity`, but when the point is `0`, `ln(0)` returns `-Infinity` (a Mul node), and then `x * (-Infinity)` at `x=0` may not produce a recognizable `NaN`.

**Correct behavior:** Should rewrite as `ln(x)/(1/x)` and apply L'Hôpital → `d(ln(x))/d(1/x) = (1/x)/(-1/x²) = -x → 0`.

**Test case:**
```js
// BUG: returns unevaluated limit, should be 0
limit(x*ln(x), x, 0, 1)  // → 0 (by L'Hôpital)
```

---

## Bug #9: `diff()` doesn't recognize `arcsin`/`arccos`/`arctan` aliases
**Location:** `expression.js:Call.diff()` 
**Severity:** 🟡 Major
**Description:** The derivative rules recognize `asin`, `acos`, `atan` but NOT `arcsin`, `arccos`, `arctan`. Users commonly use the `arc`-prefix notation, and these return unevaluated `diff(arcsin(x), x)`.

**Correct behavior:** Both `asin`/`arcsin` etc. should produce the same derivative.

**Test cases:**
```js
// BUG: returns diff(arcsin(x), x) unevaluated
diff(arcsin(x), x)  // Expected: 1/sqrt(1-x²)

// BUG: returns diff(arccos(x), x) unevaluated
diff(arccos(x), x)  // Expected: -1/sqrt(1-x²)

// BUG: returns diff(arctan(x), x) unevaluated
diff(arctan(x), x)  // Expected: 1/(1+x²)

// OK:
diff(asin(x), x)  // → 1/sqrt(1-x²) ✓
diff(atan(x), x)  // → 1/(1+x²) ✓
```

---

## Bug #10: `solve()` returns duplicate roots for double roots
**Location:** `cas.js:_solve()`
**Severity:** 🟢 Minor
**Description:** `solve((x-1)², x)` returns `{1, 1}` (duplicate), while `solve((x-1)³, x)` returns just `1`. The deduplication is inconsistent. `_partfrac` actually deduplicates roots internally, but `_solve` itself does not, leading to inconsistent output.

**Correct behavior:** `solve()` should return each distinct root once. If multiplicity is needed, it should be represented differently (e.g., `{1: 2}` or separate multiplicity info).

**Test cases:**
```js
// BUG: returns {-1, -1}, should return {-1} or handle multiplicity
solve(x²+2*x+1, x)

// BUG: returns {0, 0, 0, 0}, should return {0}
solve(x⁴, x)

// Inconsistent: returns just 1 (no duplicates for triple root)
solve((x-1)³, x)
```

---

## Bug #11: `integrate(x²/(x²+1), x)` uses complex partial fractions instead of polynomial division
**Location:** `cas.js:_partfrac()` + integrate fallback
**Severity:** 🟡 Major
**Description:** When integrating `x²/(x²+1)`, the CAS should first do polynomial division: `x²/(x²+1) = 1 - 1/(x²+1)`, then integrate to `x - atan(x)`. Instead, `_partfrac` decomposes over ℂ (Bug #2) and the result contains `quo(x², (x²+1))` (unevaluated polynomial division) plus complex logarithms.

**Test case:**
```js
// BUG: returns integrate(quo(...)...) + complex logs
integrate(x²/(x²+1), x)
// Expected: x - atan(x)
```

---

## Bug #12: `integrate(1/sqrt(x²+1), x)` returns `asinh(x)` instead of `ln(x+sqrt(x²+1))`
**Location:** `cas.js:_integrateInverseHyperbolic()`
**Severity:** 🟢 Minor
**Description:** Returns `asinh(x)` which is mathematically correct but many users/CAS systems expect the logarithmic form `ln(x + sqrt(x²+1))`. Both are valid, this is more of a presentation preference.

**Test case:**
```js
// Returns asinh(x), could also return ln(x+sqrt(x²+1))
integrate(1/sqrt(x²+1), x)
```

---

## Bug #13: `integrate(1/sqrt(1-x²), x)` returns `asin(x)` instead of `arcsin(x)` 
**Location:** `expression.js:Div.integrate()` or CAS handler
**Severity:** 🟢 Minor  
**Description:** Returns `asin(x)` which is correct but inconsistent with user input if they use `arcsin`. Related to the alias issue in Bug #9.

---

## Summary

| # | Severity | Component | Summary |
|---|----------|-----------|---------|
| 1 | 🔴 | `_partfrac` | Treats transcendental denominators as polynomials |
| 2 | 🔴 | `_partfrac` | Decomposes irreducible quadratics over ℂ |
| 3 | 🔴 | `_integrateByParts` | Infinite recursion on cyclic IBP |
| 4 | 🔴 | `integrate` | sin²/cos² produce self-referential result |
| 5 | 🔴 | `_limit` | `inf` not recognized as Infinity |
| 6 | 🟡 | `_limit` | ∞/∞ L'Hôpital fails with `inf` |
| 7 | 🟡 | `_limit` | Oscillating limits at ∞ not detected |
| 8 | 🟡 | `_limit` | 0·∞ form with ln not handled |
| 9 | 🟡 | `diff` | arcsin/arccos/arctan aliases missing |
| 10 | 🟢 | `solve` | Duplicate roots for double roots |
| 11 | 🟡 | `integrate` | x²/(x²+1) uses complex PF instead of division |
| 12 | 🟢 | `integrate` | asinh vs ln form preference |
| 13 | 🟢 | `integrate` | asin vs arcsin inconsistency |
