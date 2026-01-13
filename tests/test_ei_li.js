
const { CAS } = require('../js/cas');
const { Parser, Lexer } = require('../js/parser');
require('../js/expression');

const cas = new CAS();

function assert(condition, message) {
    if (!condition) {
        console.error(`FAIL: ${message}`);
        process.exit(1);
    } else {
        console.log(`PASS: ${message}`);
    }
}

function eval(input) {
    const lexer = new Lexer(input);
    const parser = new Parser(lexer);
    const expr = parser.parse();
    return cas.evaluate(expr);
}

function approx(a, b, eps=1e-5) {
    return Math.abs(a - b) < eps;
}

// 1. Integration Tests
const intEi = eval("integrate(exp(x)/x, x)").toString();
assert(intEi === "Ei(x)", `integrate(exp(x)/x, x) should be Ei(x), got ${intEi}`);

const intLi = eval("integrate(1/ln(x), x)").toString();
assert(intLi === "Li(x)", `integrate(1/ln(x), x) should be Li(x), got ${intLi}`);

// 2. Numeric Tests
// Ei(1) approx 1.895117816
const ei1 = eval("Ei(1)").value;
assert(approx(ei1, 1.895117816), `Ei(1) approx 1.8951, got ${ei1}`);

// Ei(-1) approx -0.219383934
const ei_1 = eval("Ei(-1)").value;
assert(approx(ei_1, -0.219383934), `Ei(-1) approx -0.21938, got ${ei_1}`);

// Li(2) = Ei(ln(2)) approx 1.04516378
const li2 = eval("Li(2)").value;
assert(approx(li2, 1.04516378), `Li(2) approx 1.04516, got ${li2}`);

// 3. Differentiation Tests
// diff(Ei(x), x) = exp(x)/x
const diffEi = eval("diff(Ei(x), x)").toString();
// Expected: (exp(x) / x) or similar
assert(diffEi === "(exp(x) / x)" || diffEi === "exp(x) / x", `diff(Ei(x)) should be exp(x)/x, got ${diffEi}`);

// diff(Li(x), x) = 1/ln(x)
const diffLi = eval("diff(Li(x), x)").toString();
assert(diffLi === "(1 / ln(x))" || diffLi === "1 / ln(x)", `diff(Li(x)) should be 1/ln(x), got ${diffLi}`);

// 4. Integral of Ei and Li
// integrate(Ei(x), x) = x*Ei(x) - exp(x)
const intEi2 = eval("integrate(Ei(x), x)").toString();
// check if contains x*Ei(x) and -exp(x)
const targetEi2 = "((x * Ei(x)) - exp(x))";
assert(intEi2 === targetEi2 || intEi2 === "x * Ei(x) - exp(x)", `integrate(Ei(x)) failed, got ${intEi2}`);

// integrate(Li(x), x) = x*Li(x) - Ei(2*ln(x))
const intLi2 = eval("integrate(Li(x), x)").toString();
const targetLi2 = "((x * Li(x)) - Ei((2 * ln(x))))";
const targetLi3 = "((x * Li(x)) - Ei(ln(x^2)))";
assert(intLi2 === targetLi2 || intLi2 === targetLi3, `integrate(Li(x)) failed, got ${intLi2}`);

console.log("All Ei/Li tests passed!");
