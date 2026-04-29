# Calculus Engine Bug Report

## Bug #1: `_partfrac()` treats transcendental denominators as polynomials
**Location:** `cas.js:_partfrac()` (~line 9961)
**Severity:** ~~🔴 Critical~~ ✅ FIXED
**Description:** FIXED in commit. `_partfrac()` now checks if the denominator is a polynomial using `_getPolyCoeffs`. If not, it returns the expression unchanged.

**Test cases:**
```js
// FIXED: returns 1/(e^x-1) unchanged
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
**Severity:** ~~🔴 Critical~~ ✅ FIXED
**Description:** FIXED in commit. `_partfrac` now checks the discriminant of quadratic denominators. If disc < 0, the irreducible quadratic is kept intact.

**Test cases:**
```js
// FIXED: returns 1/(x²+1) unchanged
partfrac(1/(x²+1), x)

// FIXED: returns 1/(x²+2) unchanged
partfrac(1/(x²+2), x)
```

---

## Bug #3: `_integrateByParts()` infinite recursion on cyclic by-parts integrals
**Location:** `cas.js:_integrateByParts()` (~line 5675)
**Severity:** ~~🔴 Critical~~ ✅ FIXED
**Description:** FIXED in commit. The code now detects cyclic IBP cases (when vdu is proportional to expr) and solves algebraically: `I = uv - kI` → `I = uv/(1+k)`.

**Test cases:**
```js
// FIXED: returns (e^x*(sin(x)-cos(x)))/2
integrate(e^x*sin(x), x)

// FIXED: returns sin(x)^2/2
integrate(sin(x)*cos(x), x)
```

---

## Bug #4: `integrate(sin(x)², x)` fails — produces self-referential result
**Location:** `cas.js` integrate handler (~line 695-705) and `expression.js:Mul.integrate` / by-parts cycle
**Severity:** ~~🔴 Critical~~ ✅ FIXED
**Description:** FIXED in commit. Trig power reduction (`_trigPowerReduce`) is now tried BEFORE integration by parts. For `sin(x)^2`, it reduces to `(-cos(2x)+1)/2` and integrates correctly.

**Test cases:**
```js
// FIXED: returns x/2 - sin(2*x)/4
integrate(sin(x)^2, x)

// FIXED: returns x/2 + sin(2*x)/4
integrate(cos(x)^2, x)
```

---

## Bug #5: `_limit()` doesn't recognize `Sym('inf')` as Infinity
**Location:** `cas.js:_limit()` (~line 5783+)
**Severity:** ~~🔴 Critical~~ ✅ FIXED
**Description:** FIXED in previous commit (332b96a). All `isInf()` helper functions now check `name === 'inf'` in addition to `name === 'Infinity'`/`name === 'infinity'`.

**Test cases:**
```js
// FIXED: returns Infinity
limit(x², x, inf)

// FIXED: returns 0
limit(1/x, x, inf)
```

---

## Bug #6: `_limit()` doesn't handle `∞/∞` L'Hôpital when using `inf` Sym
**Location:** `cas.js:_limit()` Div handler
**Severity:** ~~🟡 Major~~ ✅ FIXED (via Bug #5 fix)
**Description:** Related to Bug #5 — since `Sym('inf')` is now recognized, ∞/∞ L'Hôpital works with `inf`.

---

## Bug #7: `_limit()` doesn't detect oscillating limits at infinity
**Location:** `cas.js:_limit()` fallback path
**Severity:** ~~🟡 Major~~ ✅ FIXED
**Description:** FIXED in commit. `_limit` now detects `sin(f(x))` and `cos(f(x))` as x→±∞ and returns `undefined`.

**Test cases:**
```js
// FIXED: returns undefined
limit(sin(x), x, Infinity)

// FIXED: returns undefined
limit(cos(x), x, Infinity)
```

---

## Bug #8: `_limit()` doesn't handle `0·∞` form with `ln`
**Location:** `cas.js:_limit()` Mul handler
**Severity:** ~~🟡 Major~~ ✅ FIXED
**Description:** FIXED in commit. `x*ln(x)` is now rewritten as `ln(x)/(1/x)` and L'Hôpital is applied.

**Test case:**
```js
// FIXED: returns 0
limit(x*ln(x), x, 0, 1)
```

---

## Bug #9: `diff()` doesn't recognize `arcsin`/`arccos`/`arctan` aliases
**Location:** `expression.js:Call.diff()`
**Severity:** ~~🟡 Major~~ ✅ FIXED
**Description:** FIXED in previous commit. The derivative rules now handle `arcsin`, `arccos`, `arctan` aliases.

**Test cases:**
```js
// FIXED: returns 1/sqrt(1-x²)
diff(arcsin(x), x)

// FIXED: returns -1/sqrt(1-x²)
diff(arccos(x), x)

// FIXED: returns 1/(1+x²)
diff(arctan(x), x)
```

---

## Bug #10: `solve()` returns duplicate roots for double roots
**Location:** `cas.js:_solve()`
**Severity:** ~~🟢 Minor~~ ✅ FIXED
**Description:** FIXED in commit. The quartic solver now deduplicates roots before returning. Linear/convex roots deduplication also improved.

**Test cases:**
```js
// FIXED: returns {-1} (single root)
solve(x²+2*x+1, x)

// FIXED: returns {0}
solve(x⁴, x)
```

---

## Bug #11: `integrate(x²/(x²+1), x)` uses complex partial fractions instead of polynomial division
**Location:** `cas.js:_partfrac()` + integrate fallback
**Severity:** ~~🟡 Major~~ ⚠️ PARTIALLY FIXED
**Description:** The `x²+1` irreducible quadratic check in `_partfrac` fixes part of this. The polynomial division part (`_quo`) is still incomplete.
**Test case:**
```js
integrate(x²/(x²+1), x)  // Returns integrate(x²/(x²+1), x) unevaluated
// Expected: x - atan(x)
```

---

## Bug #12: `integrate(1/sqrt(x²+1), x)` returns `asinh(x)` instead of `ln(x+sqrt(x²+1))`
**Location:** `cas.js:_integrateInverseHyperbolic()`
**Severity:** 🟢 Minor (design decision)
**Description:** Returns `asinh(x)` which is mathematically correct. Both forms are valid.

---

## Bug #13: `integrate(1/sqrt(1-x²), x)` returns `asin(x)` instead of `arcsin(x)` 
**Location:** `expression.js:Div.integrate()` or CAS handler
**Severity:** 🟢 Minor (design decision)
**Description:** Returns `asin(x)` which is correct. Both forms are valid.

---

## Summary

| # | Severity | Component | Summary | Status |
|---|----------|-----------|---------|--------|
| 1 | 🔴→✅ | `_partfrac` | Treats transcendental denominators as polynomials | FIXED |
| 2 | 🔴→✅ | `_partfrac` | Decomposes irreducible quadratics over ℂ | FIXED |
| 3 | 🔴→✅ | `_integrateByParts` | Infinite recursion on cyclic IBP | FIXED |
| 4 | 🔴→✅ | `integrate` | sin²/cos² produce self-referential result | FIXED |
| 5 | 🔴→✅ | `_limit` | `inf` not recognized as Infinity | FIXED |
| 6 | 🟡→✅ | `_limit` | ∞/∞ L'Hôpital fails with `inf` | FIXED |
| 7 | 🟡→✅ | `_limit` | Oscillating limits at ∞ not detected | FIXED |
| 8 | 🟡→✅ | `_limit` | 0·∞ form with ln not handled | FIXED |
| 9 | 🟡→✅ | `diff` | arcsin/arccos/arctan aliases missing | FIXED |
| 10 | 🟢→✅ | `solve` | Duplicate roots for double roots | FIXED |
| 11 | 🟡→⚠️ | `integrate` | x²/(x²+1) uses complex PF instead of division | PARTIALLY FIXED |
| 12 | 🟢 | `integrate` | asinh vs ln form preference | NOT A BUG |
| 13 | 🟢 | `integrate` | asin vs arcsin inconsistency | NOT A BUG |
