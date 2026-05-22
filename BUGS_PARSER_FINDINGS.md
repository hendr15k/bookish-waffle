# Parser.js Bug Hunt Results

Searched `js/parser.js` for edge cases, incorrect token handling, associativity issues, and malformed AST creation. 11 bugs found (2 critical, 6 medium, 3 low). Test environment: Node.js, direct expression parsing.

---

## Bug #1: Chained comparisons silently truncated [🔴 Critical]

- **Location:** `parser.js` — `compExpr()` (line 823–847)
- **Description:** `compExpr()` parses exactly one comparison. After `0 < x`, it returns `Lt(0, x)` and the trailing `< 5` is silently discarded — no error.
- **Correct behavior:** `0<x<5` should parse as `And(Lt(0, x), Lt(x, 5))`.
- **Test:**
  ```js
  parse('0<x<5').toString()    // "0 < x"  ← "< 5" lost
  parse('a<b<c<d').toString() // "a < b"  ← "<c<d" lost
  parse('x>=1 && x<=10').toString() // "x >= 1"  ← broken
  ```
- **Severity:** 🔴 Critical — silent data loss.

---

## Bug #2: Subscript/array index on LHS of assignment ignored [🔴 Critical]

- **Location:** `parser.js` — `statement()` (line 896–1024)
- **Description:** `statement()` handles `name:=value` but not `expr:=value`. When it sees `a[0]:=5`, `factor()` returns `At(Sym('a'), Num(0))`, then the `:=` check fails because the LHS is an `At` node (not `Sym`). The assignment is discarded.
- **Correct behavior:** `a[0]:=5` should be `Assignment(At(Sym('a'), Num(0)), Num(5))`.
- **Test:**
  ```js
  parse('a[0]:=5').toString()   // "a[0]"  ← assignment lost
  parse('a[0] := 5').toString() // "a[0]"  ← same
  ```
- **Severity:** 🔴 Critical — array assignment silently does nothing.

---

## Bug #3: Double factorial `!!` parsed as single factorial [🟡 Medium]

- **Location:** `parser.js` — `power()` (line 727–739), `atom()` (line 729–731)
- **Description:** `power()` matches `!` in `factor()`. After `n!`, the second `!` is not consumed — `factor()` ends, `power()` ends, `unary()` ends. `n!!` becomes `fact(n)`.
- **Correct behavior:** `n!!` = `fact(fact(n))` (for integer n). At minimum, second `!` must not be silently lost.
- **Test:**
  ```js
  parse('n!!').toString()  // "fact(n)"  ← second ! lost
  parse('5!!').toString()  // "fact(5)"  // should be 120!
  ```
- **Severity:** 🟡 Medium — silently wrong result.

---

## Bug #4: Double pipe `||x||` fails with syntax error [🟡 Medium]

- **Location:** `parser.js` — `Lexer.getNextToken()` (line 350–356)
- **Description:** `||` is tokenized as `TOKEN_OR` (logical OR). After `||x|`, the lexer returns `TOKEN_PIPE`. The parser sees `|x|` where binary `|` is invalid, then another `|` with no matching open.
- **Correct behavior:** `||x||` should parse as `abs(abs(x))` (norm notation).
- **Test:**
  ```js
  parse('||x||')          // "Invalid syntax"
  parse('|x|').toString() // "abs(x)"  ← works
  // Tokens for "||x||": OR:||, IDENTIFIER:x, PIPE:|
  ```
- **Severity:** 🟡 Medium — norm notation broken.

---

## Bug #5: Subscript with braces loses operators [🟡 Medium]

- **Location:** `parser.js` — `atom()` subscript loop (line 425–440)
- **Description:** Inside `{...}`, only `IDENTIFIER` and `NUMBER` tokens are consumed. Operators like `+`, `-`, `*` are skipped. Result: `x_{a+b}` → `x_ab`.
- **Correct behavior:** Complex subscripts should be preserved or give an error.
- **Test:**
  ```js
  parse('x_{a+b}').toString()  // "x_ab"  ← "+" lost
  parse('x_{y+1}').toString()  // "x_y1"  ← "+" lost
  ```
- **Severity:** 🟡 Medium — data loss in subscript content.

---

## Bug #6: Sci notation followed by identifier merges unexpectedly [🟡 Medium]

- **Location:** `parser.js` — `Lexer.number()` (line 104–142), `Lexer.identifier()` (line 144–152)
- **Description:** After `1e10`, the lexer sees `e5` as a new token and parses it as identifier `e5`. With implicit multiplication, `1e10e5` becomes `(1e10) * (e5)`.
- **Correct behavior:** Either error or treat as `1e100000`.
- **Test:**
  ```js
  parse('1e10e5').toString()  // "(10000000000 * e5)"  ← wrong
  parse('1e1e2').toString()  // "(10 * e2)"            ← wrong
  parse('1e').toString()     // "(1 * e)"              ← also odd
  ```
- **Severity:** 🟡 Medium — confusing output.

---

## Bug #7: `++i` → `i`; `--i` → double negation [🟡 Medium]

- **Location:** `parser.js` — `Lexer` / `unary()` (line 710–725)
- **Description:** No `++`/`--` token types exist. `++i` parses as `unary()` → `+` → eat → `unary()` → `i`. Two unary plusses cancel. `--i` produces `(-1 * (-1 * i))`.
- **Correct behavior:** Should be self-assignments or errors, not identity/double-negation.
- **Test:**
  ```js
  parse('++i').toString()  // "i"                 ← confusing
  parse('--i').toString()  // "(-1 * (-1 * i))"  ← wrong for decrement
  ```
- **Severity:** 🟡 Medium — C-like semantics broken.

---

## Bug #8: Conditional expression `if(a,b,c)` fails [🟡 Medium]

- **Location:** `parser.js` — `atom()` (line 896–911)
- **Description:** `if` is only a control-flow statement (requires parens around condition, block body). Python-style `if(cond, trueVal, falseVal)` is not recognized as a `Call` expression.
- **Correct behavior:** Support `if(condition, trueVal, falseVal)` or give a clear "use block form" error.
- **Test:**
  ```js
  parse('if(x>0, x, -x)')                    // "Invalid syntax"
  parse('if(x>0) { x } else { -x }').toString() // works
  ```
- **Severity:** 🟡 Medium — common functional idiom not supported.

---

## Bug #9: For loop comma separator not supported [🟢 Low]

- **Location:** `parser.js` — `statement()` for-loop (line 921–945)
- **Description:** For loop requires semicolon-separated parts. Common comma syntax `for(i=0,i<5,i+1)` fails.
- **Test:**
  ```js
  parse('for(i=0,i<5,i+1)')    // "Invalid syntax"
  parse('for(i=0; i<5; i+1)').toString() // works
  ```
- **Severity:** 🟢 Low — may be intentional.

---

## Bug #10: `x^2!` parses as `x^fact(2)` [🟢 Low]

- **Location:** `parser.js` — `power()` (line 727–739), `atom()` (line 729–731)
- **Description:** `!` binds to the atom immediately to its left. `x^2!` = `Pow(x, fact(2))` = `x^2`. Ambiguous — could mean `(x^2)!`.
- **Test:**
  ```js
  parse('x^2!').toString()   // "x^fact(2)" = x^2
  parse('-x!').toString()     // "(-1 * fact(x))" = -(x!)
  ```
- **Severity:** 🟢 Low — precedence ambiguous but defensible.

---

## Bug #11: Compound assignment operators missing [🟢 Low]

- **Location:** `parser.js` — no `+=`, `-=`, `*=`, etc. tokens/handling
- **Description:** `i+=1`, `i-=1`, `i*=2` all fail.
- **Test:**
  ```js
  parse('i+=1')  // "Invalid syntax"
  ```
- **Severity:** 🟢 Low — uncommon in CAS expressions.

---

## Summary

| # | Bug | Severity | Impact |
|---|-----|----------|--------|
| 1 | Chained comparisons silently truncated | 🔴 Critical | `0<x<5` → `0<x`, data loss |
| 2 | Subscript LHS of assignment ignored | 🔴 Critical | `a[0]:=5` → `a[0]`, assignment lost |
| 3 | `!!` parsed as single `!` | 🟡 Medium | `n!!` → `fact(n)`, wrong |
| 4 | `||x||` fails | 🟡 Medium | Norm notation broken |
| 5 | Subscript braces lose operators | 🟡 Medium | `x_{a+b}` → `x_ab` |
| 6 | Sci notation + identifier merges | 🟡 Medium | `1e10e5` → `(1e10 * e5)` |
| 7 | `++i`/`--i` wrong | 🟡 Medium | Confusing C-like semantics broken |
| 8 | `if(a,b,c)` not supported | 🟡 Medium | Functional if-expr missing |
| 9 | For loop comma separator | 🟢 Low | Minor compatibility |
| 10 | `x^2!` → `x^fact(2)` | 🟢 Low | Precedence ambiguity |
| 11 | Compound assignment missing | 🟢 Low | Minor |
