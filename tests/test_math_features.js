
const { Expr, Num, Sym, BinaryOp, Add, Sub, Mul, Div, Pow, Call, Assignment, Eq, Vec, FunctionDef, Block, toExpr, And, Or, Xor, Implies, Iff, Not, Mod, Neq, Lt, Gt, Le, Ge, At, BooleanEq, If, While, For, Return, Break, Continue } = require('../js/expression');
const { Parser, Lexer } = require('../js/parser');
const { CAS } = require('../js/cas');

const cas = new CAS();

function test(exprStr, expectedStr) {
    try {
        const lexer = new Lexer(exprStr);
        const parser = new Parser(lexer);
        const expr = parser.parse();
        const res = cas.evaluate(expr);
        const resStr = res.toString();
        if (resStr === expectedStr) {
            console.log(`[PASS] ${exprStr} => ${resStr}`);
        } else {
            console.error(`[FAIL] ${exprStr} => ${resStr}, expected ${expectedStr}`);
        }
    } catch (e) {
        console.error(`[ERROR] ${exprStr}: ${e.message}`);
    }
}

console.log("--- Testing Math Features ---");

// Bernoulli (Existing feature check)
test("bernoulli(0)", "1");
test("bernoulli(1)", "(-1 / 2)");
test("bernoulli(2)", "(1 / 6)");
test("bernoulli(4)", "(-1 / 30)");

// Convolution
test("conv([1, 2], [3, 4])", "[3, 10, 8]");
test("conv([1, 1, 1], [1, 1])", "[1, 2, 2, 1]");

// Cross Correlation
test("xcorr([1, 2], [1, 2])", "[2, 5, 2]");

// Matrix Helpers
test("vandermonde([1, 2])", "[[1, 1], [1, 2]]");
test("hilbert(2)", "[[1, (1 / 2)], [(1 / 2), (1 / 3)]]");
test("toeplitz([1, 2], [1, 3])", "[[1, 3], [2, 1]]");

// Pade Approximation
try {
    const lexerPade = new Lexer("pade(exp(x), x, 1, 1)");
    const parserPade = new Parser(lexerPade);
    const padeRes = cas.evaluate(parserPade.parse());
    console.log("pade(exp(x), x, 1, 1) => " + padeRes.toString());
    // Also check taylor
    const lexerTaylor = new Lexer("taylor(exp(x), x, 0, 2)");
    const parserTaylor = new Parser(lexerTaylor);
    const taylorRes = cas.evaluate(parserTaylor.parse());
    console.log("taylor(exp(x), x, 0, 2) => " + taylorRes.toString());
} catch(e) { console.error(e); }

// Product Differentiation
// diff(product(x, k, 1, n), x) = product(x) * sum(1/x) = x^n * n/x = n*x^(n-1)
// Let's see if it produces symbolic form
try {
    const lexerDiff = new Lexer("diff(product(x, k, 1, n), x)");
    const parserDiff = new Parser(lexerDiff);
    const diffProd = cas.evaluate(parserDiff.parse());
    console.log("diff(product(x, k, 1, n), x) => " + diffProd.toString());
} catch(e) { console.error(e); }

console.log("--- Done ---");
