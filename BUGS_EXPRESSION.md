# BUGS_EXPRESSION — expression.js Bugs

Found via code analysis + Node.js testing. All tests run from workspace root.

---

## 🔴 Critical

### BUG #1 — `Mul.expand()` doesn't distribute when left side is Sym/Num
**Location:** `expression.js:1320-1327` (`Mul.expand()`)
**Description:** `Mul.expand()` distributes multiplication over addition only when the operand is already an `Add`/`Sub` node. It first calls `expand()` on both sides and stores results in `l` and `r`, then checks `l instanceof Add` and `r instanceof Add`. But if `l` and `r` are `Sym`/`Num` (not `Add`/`Sub`), both checks fail and it returns `new Mul(l, r)` — which is the unexpanded form. This means `Sym * Add` or `Num * Add` never distributes.

The second problem is that `Add.simplify()` at line 726-731 refactors `x*y + x*z` back into `x*(y+z)` via common factor extraction, completely defeating expand.

Together these two bugs make `expand()` fail for all expressions where the left operand is not an `Add`/`Sub`.

**Root cause:** Two-part bug:
1. `Mul.expand()` only checks `this.left instanceof Add` and `this.right instanceof Add` (original sides), not the expanded `l` and `r`.
2. `Add.simplify()` re-factors `x*y + x*z` → `x*(y+z)`.

**Correct behavior:**
- `expand(x*(y+z))` → `x*y + x*z`
- `expand(2*(x+y))` → `2*x + 2*y`
- `expand((x+1)*(y+2))` → `x*y + 2*x + y + 2`

**Test case:**
```js
const {Parser, Lexer} = require('./js/parser.js');
const { CAS } = require('./js/cas.js');
const cas = new CAS();
function parse(s) { return new Parser(new Lexer(s)).parse(); }

cas.evaluate(parse('expand(x*(y+z))')).toString();
// Got: (x * (y + z))
// Expected: (x * y) + (x * z)

cas.evaluate(parse('expand(2*(x+y))')).toString();
// Got: (2 * (x + y))
// Expected: (2 * x) + (2 * y)

cas.evaluate(parse('expand((x+1)*(y+2))')).toString();
// Got: ((y + 2) * (x + 1))
// Expected: (x * y) + (2 * x) + (y + 2)
```

**Severity:** 🔴 Critical — `expand()` is fundamentally broken for most expressions

---

## 🟡 Moderate

### BUG #2 — `Sub.simplify()` doesn't flatten `Sub(Num(0), x)` → `-x`
**Location:** `expression.js:638` and throughout (`Sub.simplify()`)
**Description:** When unary minus creates `Sub(Num(0), x)` (e.g., `(-1)*x` simplifies to `(0 - x)`), there's no simplification rule to convert this to `-x`. The pattern `if (r instanceof Num && r.value < 0) return new Add(l, new Num(-r.value)).simplify()` at line 953 only handles adding a positive number, not the `Sub(0, x)` pattern.

Also, `Mul.simplify()` creates `Sub(Num(0), expr)` for `(-1) * expr` (line 1195: `return new Sub(new Num(0), r)`), which is the correct structure but never gets cleaned up.

**Correct behavior:**
- `(-1)*x` → `-x` (or at least `-1*x` without nested Sub)
- `(0 - x)` → `-x`

**Test case:**
```js
simp("(-1)*x")  // Got: (0 - x)
simp("0-x")     // Got: (0 - x)
simp("-x")      // Got: (0 - x)
simp("--x")     // Got: x ← works correctly (Sub.simplify cancels)
```

**Severity:** 🟡 Moderate — ugly/verbose output, nested Sub structure

---

### BUG #3 — `Mul.simplify()` doesn't handle `Sub(0,x) * y` or `(-x) * y`
**Location:** `expression.js:1043-1328` (`Mul.simplify()`)
**Description:** When a `Sub(Num(0), x)` (representing `-x`) is on the left side of a multiplication, there's no simplification. `(-x)*y` stays as `((0 - x) * y)`. Similarly, `(-x)*(-y)` stays as `((0 - x) * (0 - y))` instead of simplifying to `x*y`.

**Correct behavior:**
- `(-x)*y` → `-(x*y)` or `-x*y`
- `(-x)*(-y)` → `x*y`
- `x/(-y)` → `-(x/y)` or `-x/y`
- `(-x)/y` → `-(x/y)` or `-x/y`

**Test case:**
```js
simp("(-x)*(-y)")  // Got: ((0 - x) * (0 - y)) — not simplified to x*y
simp("(-x)*(y)")   // Got: ((0 - x) * y)
simp("x/(-y)")     // Got: (x / (0 - y))
simp("(-x)/y")     // Got: ((0 - x) / y)
```

**Severity:** 🟡 Moderate — affects negative coefficient handling

---

### ~~BUG #4 — `Mul.simplify()` doesn't simplify power of powers `(a^x)^y` → `a^(x*y)`~~ ✅ FIXED
**Location:** `expression.js` (`Pow.simplify()`)
**Description:** Fixed in Pow.simplify() via the `(a^b)^c = a^(b*c)` rule at line ~2307. Now correctly simplifies `(a^x)^y` → `a^(x*y)`.

**Test case:**
```js
simp("(a^x)^y")   // Got: a^((x * y))  ✓
simp("(x^2)^3")   // Got: x^6          ✓
simp("((a^x)^y)^z") // Got: a^(((x * y) * z))  ✓
```

**Severity:** 🟡 Moderate — common symbolic exponent simplification
**Fixed:** HEN-90

---

### BUG #5 — `Mul.simplify()` doesn't combine `a^x * a^y` → `a^(x+y)`
**Location:** `expression.js` (`Mul.simplify()`)
**Description:** The code handles `x * x^n → x^(n+1)` and `x^n * x^m → x^(n+m)` at lines 1228-1241, but only for matching `toString()`. It doesn't handle `Pow(base, exp1) * Pow(base, exp2)` when `base` is a compound expression like `a`.

**Correct behavior:**
- `simplify(a^x * a^y)` → `a^(x+y)`
- `simplify(2^x * 2^y)` → `2^(x+y)`

**Test case:**
```js
simp("a^x*a^y")  // Got: (a^x * a^y)
simp("2^x*2^y")  // Got: (2^x * 2^y)
```

**Severity:** 🟡 Moderate — common exponent arithmetic

---

### BUG #6 — `Mul.simplify()` doesn't combine `a^x / a^y` → `a^(x-y)`
**Location:** `expression.js` (`Div.simplify()`)
**Description:** `Div.simplify()` doesn't have a rule for `Pow(base, exp1) / Pow(base, exp2)` → `Pow(base, Sub(exp1, exp2))`.

**Correct behavior:**
- `simplify(a^x/a^y)` → `a^(x-y)`
- `simplify(2^x/2^y)` → `2^(x-y)`

**Test case:**
```js
simp("a^x/a^y")  // Got: a^((x - y))
// Actually works! The Div.simplify handles Pow/Pow correctly.
```

**Note:** Actually this DOES work — Div.simplify correctly simplifies `a^x/a^y` → `a^(x-y)`. Confirmed working in tests.

**Severity:** 🟢 — Not actually a bug

---

## 🟢 Minor

### BUG #7 — `Call.simplify()` doesn't simplify `log(e^x)` → `x*ln(e)`
**Location:** `expression.js` (`Call.simplify()`)
**Description:** `log(exp(x))` doesn't simplify to `x` because there's no rule matching `log(Pow(e, x))`. The `ln(e)` simplifies to `1` but `log(exp(x))` stays as `log(exp(x))`.

**Correct behavior:** `log(exp(x))` → `x`

**Test case:**
```js
simp("log(exp(x))")  // Got: log(exp(x))
// Note: exp(ln(x)) correctly simplifies to x, but log(exp(x)) does not
```

**Severity:** 🟢 Minor — inconsistency between log∘exp and exp∘log

---

## Summary

| # | Severity | Component | Summary | Status |
|---|----------|-----------|---------|--------|
| 1 | 🔴 Critical | `Mul.expand()` + `Add.simplify()` | expand() doesn't distribute Mul*Add and simplify() re-factors result | Open |
| 2 | 🟡 Moderate | `Sub.simplify()` | `Sub(0, x)` not flattened to `-x` | Open |
| 3 | 🟡 Moderate | `Mul.simplify()` | `(-x)*y`, `(-x)*(-y)` etc. not simplified | Open |
| 4 | 🟡 Moderate | `Pow.simplify()` | `(a^x)^y` not simplified to `a^(x*y)` in simplify | ✅ Fixed (HEN-90) |
| 5 | 🟡 Moderate | `Mul.simplify()` | `a^x * a^y` not combined to `a^(x+y)` | Open |
| 6 | 🟢 | — | Not a bug — confirmed working | N/A |
| 7 | 🟢 Minor | `Call.simplify()` | `log(exp(x))` not simplified | Open |
