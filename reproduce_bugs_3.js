
const { Lexer, Parser } = require('./js/parser.js');
const { CAS } = require('./js/cas.js');
const { Expr, Num, Sym, Add, Sub, Mul, Div, Pow, Call } = require('./js/expression.js');

const cas = new CAS();

function assert(condition, message) {
    if (!condition) {
        console.error(`FAIL: ${message}`);
        process.exit(1);
    } else {
        console.log(`PASS: ${message}`);
    }
}

// 1. Implicit Multiplication Precedence: 1/2x
// Should be 1/(2x) or (1/2)*x?
// Standard physics/engineering (and TI calculators) often treat implicit mult as higher precedence than division?
// Or equal?
// If equal: 1/2*x -> (1/2)*x.
// If higher: 1/2x -> 1/(2x).
// Most programming languages: (1/2)*x.
// Let's see what current parser does. It uses `term` loop for `*`, `/` and implicit.
// So `1 / 2 x` -> `(1/2) x` -> `((1/2)*x)`.
// This is standard for computer algebra usually (e.g. Python, JS), but users might expect `1/2x` to mean `1/(2x)`.
// Xcas treats `1/2x` as `1/(2*x)`.
// Let's verify current behavior.
try {
    const lex = new Lexer("1/2x");
    const par = new Parser(lex);
    const tree = par.parse();
    console.log("1/2x parsed as:", tree.toString());
} catch(e) {
    console.error(e);
}

// 2. Complex number parsing
// `1+2i`. `i` is identifier. `2i` is implicit mul.
// `1+2i` -> `1 + (2*i)`.
try {
    const lex = new Lexer("1+2i");
    const par = new Parser(lex);
    const tree = par.parse();
    console.log("1+2i parsed as:", tree.toString());
} catch(e) { console.error(e); }


// 3. sin^2(x) vs sin(x)^2
// Parser supports `identifier` followed by `(`.
// `sin^2(x)` would be parsed as `sin` (ID) -> wait, `power` calls `factor`.
// `factor` parses ID. If next is `^`, `power` handles it.
// So `sin^2` -> `Pow(Sym(sin), 2)`.
// Then `(x)` follows.
// `term` sees `Pow` then `(x)`. `(x)` is `LPAREN ...`. `factor` handles parens.
// So `term` sees `Pow` and `factor` (implicit mul).
// Result: `sin^2 * (x)`. Not `(sin(x))^2`.
// We need to support `sin^2(x)` syntax or at least check what happens.
try {
    const lex = new Lexer("sin^2(x)");
    const par = new Parser(lex);
    const tree = par.parse();
    console.log("sin^2(x) parsed as:", tree.toString());
} catch(e) { console.error(e); }

// 4. Roots
// `sqrt(x)` is supported. `x^(1/2)` should simplify to `sqrt(x)`?
const root = new Pow(new Sym('x'), new Div(new Num(1), new Num(2)));
console.log("x^(1/2) simplified:", root.simplify().toString());
// Currently `Pow` simplify doesn't do this.

// 5. Log properties
// log(a*b) -> log(a) + log(b) ?
// Usually expand() does this.
const logProd = new Call('log', [new Mul(new Sym('a'), new Sym('b'))]);
console.log("log(a*b) expanded:", logProd.expand().toString());
