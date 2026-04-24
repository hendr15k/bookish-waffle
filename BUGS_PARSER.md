# Parser & LaTeX Engine Bugs

Found via code analysis + Node.js testing. All tests run from `/tmp/bookish-waffle/`.

---

## Bug #1: `toString()` roundtrip broken for all negative numbers
- **Location:** expression.js (Num.toString), parser.js (Lexer)
- **Description:** `Num(-1).toString()` outputs `-1`. When re-parsed, the Lexer sees `-` (unary minus) applied to `1`, creating `Mul(Num(-1), Num(1))` instead of `Num(-1)`. This means ANY expression containing a negative number cannot survive a toString → re-parse roundtrip. Each roundtrip multiplies the damage: `-x` → `(-1*x)` → `((-1*1)*x)` → `(((-1*1)*1)*x)` → ...
- **Correct behavior:** `(-1 * x)` should re-parse as `Mul(Num(-1), Sym('x'))`, not `Mul(Mul(Num(-1), Num(1)), Sym('x'))`
- **Test case:**
  ```js
  const {Num, Sym, Mul} = require('./js/expression.js');
  const {Parser, Lexer} = require('./js/parser.js');
  function parse(s) { return new Parser(new Lexer(s)).parse(); }

  // Direct: -1 * x roundtrip
  const original = '(-1 * x)';
  const s1 = parse(original).toString();  // (-1 * x)
  const s2 = parse(s1).toString();        // ((-1 * 1) * x) ← WRONG
  console.log(s1 === s2);  // false ← should be true

  // Cascading: (-1 * x^2 + 3 * x - 5)
  const expr = '(-1 * x^2 + 3 * x - 5)';
  const r1 = parse(expr).toString();
  const r2 = parse(r1).toString();
  console.log(r1 === r2);  // false ← gets worse every roundtrip
  ```
- **Severity:** 🔴 Critical — breaks the core CAS feedback loop (simplify → display → edit → re-parse)

---

## Bug #2: Mul.toLatex() missing parentheses for Mul×Mul
- **Location:** expression.js:1035 (Mul.toLatex)
- **Description:** When both left and right operands are `Mul`, no parentheses are rendered. `(a*b)*(c*d)` → `a b \cdot c d` which is visually ambiguous (looks like `a·b·c·d`). The code only wraps `Add` and `Sub` children in `\left(\right)`.
- **Correct behavior:** `(a*b)*(c*d)` should render as `\left(a b\right) \cdot \left(c d\right)` or at minimum `a b \cdot c d` with explicit grouping.
- **Test case:**
  ```js
  const {Sym, Mul} = require('./js/expression.js');
  const m = new Mul(new Mul(new Sym('a'), new Sym('b')), new Mul(new Sym('c'), new Sym('d')));
  console.log(m.toLatex());
  // Output:   a b \cdot c d
  // Expected: \left(a b\right) \cdot \left(c d\right)
  ```
- **Severity:** 🟡 Medium — ambiguous LaTeX rendering for nested multiplications

---

## Bug #3: `sec`, `csc`, `cot` duplicated in both standardFunctions and mapFunctions
- **Location:** expression.js:~3410 (Call.toLatex)
- **Description:** The `standardFunctions` array includes `'sec', 'csc', 'cot'` AND the `mapFunctions` object also has entries for them. Since `standardFunctions` is checked first with an early return, the `mapFunctions` entries for sec/csc/cot are dead code. This isn't causing wrong output currently, but it's a maintenance hazard — if someone edits only one location, the outputs will diverge silently.
- **Correct behavior:** Each function should be defined in exactly one place. sec/csc/cot should be removed from either `standardFunctions` or `mapFunctions`.
- **Test case:**
  ```js
  // The mapFunctions entries for sec/csc/cot are unreachable dead code:
  // standardFunctions check (line ~3410) catches them BEFORE mapFunctions check (line ~3425)
  // Verify: both produce the same output
  const {Call, Sym} = require('./js/expression.js');
  console.log(new Call('sec', [new Sym('x')]).toLatex());  // \sec\left(x\right) — from standardFunctions
  // mapFunctions['sec'] = '\\sec' is NEVER used
  ```
- **Severity:** 🟢 Low — no wrong output yet, but maintenance hazard / dead code

---

## Bug #4: `Num.toLatex()` produces invalid LaTeX for scientific notation
- **Location:** expression.js:410 (Num.toLatex)
- **Description:** Very small numbers like `1e-10` are rendered as literal `1e-10` in LaTeX, where `e` becomes a variable (italic *e*). Should render as `1 \times 10^{-10}`.
- **Correct behavior:** Scientific notation numbers should use proper LaTeX: `1 \times 10^{-10}`
- **Test case:**
  ```js
  const {Num} = require('./js/expression.js');
  console.log(new Num(1e-10).toLatex());  // "1e-10" — WRONG, 'e' renders as variable
  console.log(new Num(2.5e15).toLatex()); // "2500000000000000" — ugly but valid
  ```
- **Severity:** 🟡 Medium — wrong rendering for scientific notation numbers

---

## Bug #5: Mul.toLatex() missing parens for Mul as right operand
- **Location:** expression.js:1035 (Mul.toLatex)
- **Description:** `a*(b*c)` renders as `a \cdot b c` — the right Mul is not wrapped, making it ambiguous whether `c` multiplies `a·b` or just `b`. Similar issue for Sym×Mul: `a · b c` looks like three separate terms.
- **Correct behavior:** `a*(b*c)` should render as `a \cdot \left(b c\right)` to clarify grouping.
- **Test case:**
  ```js
  const {Sym, Mul} = require('./js/expression.js');
  const m = new Mul(new Sym('a'), new Mul(new Sym('b'), new Sym('c')));
  console.log(m.toLatex());
  // Output:   a \cdot b c
  // Expected: a \cdot \left(b c\right)
  ```
- **Severity:** 🟡 Medium — ambiguous LaTeX output

---

## Bug #6: `Mul(-1, ...).toLatex()` renders `-1` coefficient literally for non-simple right operands
- **Location:** expression.js:1035 (Mul.toLatex), handles `-1` coefficient at line ~1048
- **Description:** The `-1` coefficient handling only covers `Sym`, `Call`, `Pow(left=Sym)`, `Add`, and `Sub` on the right side. If the right side is a `Mul`, `Div`, or `Pow(left=Call)`, the output is `-1 \cdot ...` instead of `-(...)`. For example, `(-1)*(a/b)` → `-1 \cdot \frac{a}{b}` instead of `-\frac{a}{b}`.
- **Correct behavior:** `Mul(Num(-1), Div(Sym('a'), Sym('b')))` should render as `-\frac{a}{b}` or `-\left(\frac{a}{b}\right)`.
- **Test case:**
  ```js
  const {Num, Sym, Mul, Div} = require('./js/expression.js');
  const m = new Mul(new Num(-1), new Div(new Sym('a'), new Sym('b')));
  console.log(m.toLatex());
  // Output:   -1 \cdot \frac{a}{b}
  // Expected: -\frac{a}{b}
  ```
- **Severity:** 🟡 Medium — ugly/incorrect LaTeX for negated fractions and products

---

## Bug #7: `Pow.toLatex()` wraps Pow child in unnecessary `\left(\right)`
- **Location:** expression.js:1950 (Pow.toLatex)
- **Description:** Line 1950 includes `this.left instanceof Pow` in the check for wrapping with `\left(\right)`. This means `(x^2)^3` renders as `{\left({x}^{2}\right)}^{3}` — correct but visually heavy. More importantly, `x^2^3` (which should be a parse error or right-associative `x^{2^3}`) gets double-wrapped.
- **Correct behavior:** `Pow` children typically don't need extra parens since `{x^2}^3` is already unambiguous in LaTeX (superscripts apply left-to-right). The `\left(\right)` wrapping for Pow child is unnecessary.
- **Test case:**
  ```js
  const {Num, Sym, Pow} = require('./js/expression.js');
  const p = new Pow(new Pow(new Sym('x'), new Num(2)), new Num(3));
  console.log(p.toLatex());
  // Output:   {\left({x}^{2}\right)}^{3}
  // Expected: {{x}^{2}}^{3}  (no \left/\right needed)
  ```
- **Severity:** 🟢 Low — cosmetically heavy but technically correct

---

## Summary

| # | Bug | Severity | Impact |
|---|-----|----------|--------|
| 1 | toString roundtrip broken for negative numbers | 🔴 Critical | All simplified expressions with negatives can't be re-parsed |
| 2 | Mul.toLatex missing parens for Mul×Mul | 🟡 Medium | Ambiguous LaTeX for nested products |
| 3 | sec/csc/cot duplicated in standard+map | 🟢 Low | Dead code, maintenance hazard |
| 4 | Num.toLatex scientific notation invalid | 🟡 Medium | `e` renders as variable in `1e-10` |
| 5 | Mul.toLatex missing parens for right Mul | 🟡 Medium | Ambiguous grouping in products |
| 6 | Mul(-1,...) toLatex misses Div/Mul cases | 🟡 Medium | `-1 \cdot \frac{a}{b}` instead of `-\frac{a}{b}` |
| 7 | Pow.toLatex unnecessary \left(\right) on Pow | 🟢 Low | Excessive visual weight |
