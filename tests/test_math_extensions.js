
const expression = require('../js/expression');
Object.assign(global, expression); // Make Expr, Num, etc. global for CAS

const { CAS } = require('../js/cas');
const { Parser } = require('../js/parser');
const { Lexer } = require('../js/parser');

// Helper to init CAS
function getCAS() {
    return new CAS();
}

function parse(str) {
    const lexer = new Lexer(str);
    const parser = new Parser(lexer);
    return parser.parse();
}

function testBernoulli() {
    console.log("Testing Bernoulli...");
    const cas = getCAS();

    // B0 = 1
    let res = cas.evaluate(parse("bernoulli(0)"));
    if (res.toString() !== "1") throw new Error("B0 failed: " + res);

    // B1 = -1/2
    res = cas.evaluate(parse("bernoulli(1)"));
    if (res.toString() !== "-0.5" && res.toString() !== "(-1 / 2)") throw new Error("B1 failed: " + res);

    // B2 = 1/6
    res = cas.evaluate(parse("bernoulli(2)"));
    // Expect 1/6 approx 0.1666... or rational
    console.log("B2: " + res.toString());

    // B3 = 0
    res = cas.evaluate(parse("bernoulli(3)"));
    if (res.value !== 0) throw new Error("B3 failed: " + res);

    // B4 = -1/30
    res = cas.evaluate(parse("bernoulli(4)"));
    console.log("B4: " + res.toString());
}

function testConv() {
    console.log("Testing Convolution...");
    const cas = getCAS();
    // conv([1, 2], [3, 4])
    // [1*3, 1*4+2*3, 2*4] = [3, 10, 8]
    const res = cas.evaluate(parse("conv([1, 2], [3, 4])"));
    console.log("conv([1,2], [3,4]): " + res.toString());
    if (res.toString() !== "[3, 10, 8]") throw new Error("Conv failed");
}

function testCond() {
    console.log("Testing Condition Number...");
    const cas = getCAS();
    // Identity matrix cond = 1
    let res = cas.evaluate(parse("cond(eye(2))"));
    if (res.toString() !== "1") throw new Error("Cond(I) failed");

    // Diagonal [1, 0; 0, 2] -> Singular -> Inf? Or max/min.
    // My implementation checks min singular value. If 0, Infinity.
    // [2, 0; 0, 0.5] -> cond = 2 / 0.5 = 4
    res = cas.evaluate(parse("cond([[2, 0], [0, 0.5]])"));
    console.log("Cond(diag): " + res.toString());
    if (Math.abs(res.evaluateNumeric() - 4) > 1e-9) throw new Error("Cond failed");
}

function testProductDiff() {
    console.log("Testing Product Differentiation...");
    const cas = getCAS();
    // diff(product(x, k, 1, 3), x)
    // product(x, k, 1, 3) = x*x*x = x^3. Diff = 3x^2.
    // Formula: product(...) * sum(diff(x,x)/x, ...)
    // = x^3 * sum(1/x, k, 1, 3) = x^3 * (1/x + 1/x + 1/x) = x^3 * 3/x = 3x^2.
    const res = cas.evaluate(parse("diff(product(x, k, 1, 3), x)"));
    const simplified = res.simplify();
    console.log("diff(product): " + simplified.toString());

    // Check numeric equivalence at x=2 -> 3*4 = 12
    const val = simplified.substitute(new Sym('x'), new Num(2)).evaluateNumeric();
    console.log("Val at x=2: " + val);
    if (Math.abs(val - 12) > 1e-9) throw new Error("Product diff failed");
}

function testPade() {
    console.log("Testing Pade Approximation...");
    const cas = getCAS();
    // pade(exp(x), x, 1, 1) -> (1 + x/2) / (1 - x/2) approx

    const res = cas.evaluate(parse("pade(exp(x), x, 1, 1)"));
    console.log("Pade(exp, 1, 1): " + res.toString());

    // Check value at x=1
    // (2+1)/(2-1) = 3
    const val = res.substitute(new Sym('x'), new Num(1)).evaluateNumeric();
    console.log("Value at x=1: " + val);
    if (Math.abs(val - 3) > 1e-9) throw new Error("Pade failed");
}

try {
    testBernoulli();
    testConv();
    testCond();
    testProductDiff();
    testPade();
    console.log("All tests passed!");
} catch (e) {
    console.error(e);
    process.exit(1);
}
