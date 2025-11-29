
const { Expr, Num, Sym, Add, Sub, Mul, Div, Pow, Call } = require('./js/expression.js');
const { CAS } = require('./js/cas.js');

const cas = new CAS();

function assert(condition, message) {
    if (!condition) {
        console.error(`FAIL: ${message}`);
        process.exit(1);
    } else {
        console.log(`PASS: ${message}`);
    }
}

// 1. Test Power LaTeX precedence: (x^2)^3 vs x^(2^3)
const x = new Sym('x');
const p1 = new Pow(new Pow(x, new Num(2)), new Num(3)); // (x^2)^3
const p2 = new Pow(x, new Pow(new Num(2), new Num(3))); // x^(2^3)

console.log("p1 LaTeX:", p1.toLatex());
console.log("p2 LaTeX:", p2.toLatex());

// p1 should look like ({x}^{2})^{3} or similar, to distinguish from x^2^3
// Current impl: `{${lTex}}^{${this.right.toLatex()}}`
// If lTex is x^2, it becomes {x^2}^3 -> x^2^3 which is ambiguous/wrong precedence.
// Standard latex: (x^2)^3 -> \left(x^{2}\right)^{3}

// 2. Test 0^0
const zeroPowZero = new Pow(new Num(0), new Num(0));
const simp = zeroPowZero.simplify();
console.log("0^0 simplified:", simp.toString());
assert(simp.value === 1, "0^0 should be 1 by convention in combinatorics/series");


// 3. Test Differentiation of Log base 10
// d/dx log(x) = 1/(x ln(10))
const log10x = new Call('log', [x]);
const diffLog = log10x.diff(x).simplify();
console.log("diff(log(x)):", diffLog.toString());
// Expected: 1/(x * ln(10))
// Implementation: `new Div(u.diff(varName), new Mul(u, new Call('ln', [new Num(10)])))`
// This seems correct in code, let's verify output.

// 4. Test Limit of -1/x^2 at 0
// -1/x^2 -> -Infinity
// limit(-1/x^2, x, 0)
const exprLimit = new Div(new Num(-1), new Pow(x, new Num(2)));
const resLimit = cas._limit(exprLimit, x, new Num(0));
console.log("limit(-1/x^2, x, 0):", resLimit.toString());
assert(resLimit.toString() !== "Infinity", "Limit of negative function should be -Infinity (or handle sign)");
