
const { CAS } = require('../js/cas');
const { Expr, Num, Sym, Call, Add, Sub, Mul, Div, Pow, Vec, Lt } = require('../js/expression');

// Mock GLOBAL classes for independent testing if needed, or rely on CAS
global.Expr = Expr;
global.Num = Num;
global.Sym = Sym;
global.Call = Call;
global.Add = Add;
global.Sub = Sub;
global.Mul = Mul;
global.Div = Div;
global.Pow = Pow;
global.Vec = Vec;

const cas = new CAS();

function test(name, exprStr, expectedStr) {
    try {
        // Since we don't have a parser here, we construct AST manually or assume 'evaluate' works if we pass AST
        // We will construct AST manually for precision.

        // Helper to evaluate
        const res = cas.evaluate(exprStr);
        const resStr = res.toString();

        if (resStr === expectedStr) {
            console.log(`[PASS] ${name}`);
        } else {
            console.error(`[FAIL] ${name}`);
            console.error(`  Expected: ${expectedStr}`);
            console.error(`  Actual:   ${resStr}`);
        }
    } catch (e) {
        console.error(`[ERROR] ${name}: ${e.message}`);
    }
}

// Helper to construct basic AST
const n = (v) => new Num(v);
const s = (v) => new Sym(v);
const call = (f, args) => new Call(f, args);

console.log("Testing New Features...");

// 1. Cyclotomic Polynomials
// Phi_1(x) = x - 1
test("Cyclotomic(1, x)", new Call('cyclotomic', [n(1), s('x')]), "(x - 1)");
// Phi_2(x) = x + 1
test("Cyclotomic(2, x)", new Call('cyclotomic', [n(2), s('x')]), "(x + 1)");
// Phi_3(x) = x^2 + x + 1
test("Cyclotomic(3, x)", new Call('cyclotomic', [n(3), s('x')]), "((x^2 + x) + 1)");
// Phi_4(x) = x^2 + 1
test("Cyclotomic(4, x)", new Call('cyclotomic', [n(4), s('x')]), "(x^2 + 1)");
// Phi_6(x) = x^2 - x + 1
test("Cyclotomic(6, x)", new Call('cyclotomic', [n(6), s('x')]), "((x^2 - x) + 1)");

// 2. Entropy
// H([0.5, 0.5]) = -0.5*log2(0.5) - 0.5*log2(0.5) = 0.5 + 0.5 = 1
test("Entropy([0.5, 0.5])", new Call('entropy', [new Vec([n(0.5), n(0.5)])]), "1");
// H([1, 0]) = 0
test("Entropy([1, 0])", new Call('entropy', [new Vec([n(1), n(0)])]), "0");

// 3. KL Divergence
// KL([0.5, 0.5], [0.5, 0.5]) = 0
test("KL([0.5, 0.5], [0.5, 0.5])", new Call('klDivergence', [new Vec([n(0.5), n(0.5)]), new Vec([n(0.5), n(0.5)])]), "0");

// 4. Sigmoid
// sigmoid(0) = 1/(1+1) = 0.5
test("sigmoid(0)", new Call('sigmoid', [n(0)]), "0.5");
// sigmoid(infinity) = 1
test("sigmoid(inf)", new Call('sigmoid', [s('infinity')]), "1");

// 5. ReLU
// relu(-5) = 0
test("relu(-5)", new Call('relu', [n(-5)]), "0");
// relu(5) = 5
test("relu(5)", new Call('relu', [n(5)]), "5");

// 6. Softplus
// softplus(0) = ln(1+1) = ln(2)
test("softplus(0)", new Call('softplus', [n(0)]), "ln(2)");

// 7. Inequality Solving
// solve(2x < 4, x) -> x < 2
// 2x - 4 < 0
test("solve(2x < 4)", new Call('solve', [new Lt(new Mul(n(2), s('x')), n(4)), s('x')]), "(x < 2)");

// 8. Polygamma
// polygamma(0, 1) = -gamma_const (approx -0.577)
// We check simplifier call
const resPoly = cas.evaluate(new Call('polygamma', [n(0), n(1)]));
if (Math.abs(resPoly.evaluateNumeric() + 0.57721) < 0.001) {
    console.log("[PASS] polygamma(0, 1) numeric");
} else {
    console.error("[FAIL] polygamma(0, 1) numeric: " + resPoly.evaluateNumeric());
}

// 9. Check diff(relu(x)) -> heaviside(x)
// relu(x) expands to x * H(x)
// diff(x*H(x)) = H(x) + x*delta(x) -> H(x) + 0 -> H(x)
// We need to check if diff works
const reluExpr = cas.evaluate(new Call('relu', [s('x')]));
const diffRelu = cas.evaluate(new Call('diff', [reluExpr, s('x')]));
// Expected: heaviside(x)
// Or equivalent
// diffRelu might come out as ((1 * heaviside(x)) + (x * dirac(x)))
// simplify should reduce it
test("diff(relu(x))", new Call('diff', [new Call('relu', [s('x')]), s('x')]), "heaviside(x)");

console.log("Tests Completed.");
