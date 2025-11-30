
const { Call, Add, Sub, Sym, Num, Mul } = require('../js/expression');

// Mock a simple test runner since the project uses node directly
function runTests() {
    let passed = 0;
    let failed = 0;

    const assertLatex = (expr, expected, message) => {
        const actual = expr.toLatex();
        if (actual === expected) {
            console.log(`[PASS] ${message}`);
            passed++;
        } else {
            console.log(`[FAIL] ${message}`);
            console.log(`  Expected: ${expected}`);
            console.log(`  Actual:   ${actual}`);
            failed++;
        }
    };

    console.log("--- Testing LaTeX Precedence for Binding Functions ---");

    // Test 1: Limit of a sum
    // limit(x + 1, x, 0)
    const limitSum = new Call('limit', [
        new Add(new Sym('x'), new Num(1)),
        new Sym('x'),
        new Num(0)
    ]);
    assertLatex(limitSum, '\\lim_{x \\to 0} \\left(x + 1\\right)', 'Limit of Sum');

    // Test 2: Integral of a sum
    // integrate(x + 1, x)
    const intSum = new Call('integrate', [
        new Add(new Sym('x'), new Num(1)),
        new Sym('x')
    ]);
    assertLatex(intSum, '\\int \\left(x + 1\\right) \\, dx', 'Indefinite Integral of Sum');

    // Test 3: Definite Integral of a sum
    // integrate(x - 1, x, 0, 1)
    const defIntSub = new Call('integrate', [
        new Sub(new Sym('x'), new Num(1)),
        new Sym('x'),
        new Num(0),
        new Num(1)
    ]);
    assertLatex(defIntSub, '\\int_{0}^{1} \\left(x - 1\\right) \\, dx', 'Definite Integral of Subtraction');

    // Test 4: Summation of a sum
    // sum(k + 1, k, 1, n)
    const sumSum = new Call('sum', [
        new Add(new Sym('k'), new Num(1)),
        new Sym('k'),
        new Num(1),
        new Sym('n')
    ]);
    assertLatex(sumSum, '\\sum_{k=1}^{n} \\left(k + 1\\right)', 'Summation of Sum');

    // Test 5: Product of a subtraction
    // product(k - 1, k, 1, n)
    const prodSub = new Call('product', [
        new Sub(new Sym('k'), new Num(1)),
        new Sym('k'),
        new Num(1),
        new Sym('n')
    ]);
    assertLatex(prodSub, '\\prod_{k=1}^{n} \\left(k - 1\\right)', 'Product of Subtraction');

    // Test 6: Check that single terms are NOT wrapped
    // limit(x, x, 0)
    const limitSingle = new Call('limit', [
        new Sym('x'),
        new Sym('x'),
        new Num(0)
    ]);
    assertLatex(limitSingle, '\\lim_{x \\to 0} x', 'Limit of Single Term');

    // Test 7: Check that multiplication is NOT wrapped (conventionally acceptable)
    // integrate(2*x, x)
    const intMul = new Call('integrate', [
        new Mul(new Num(2), new Sym('x')),
        new Sym('x')
    ]);
    // Mul.toLatex is '2x' (implicit) or '2 \cdot x'
    // Mul already handles its own parens if needed.
    // 2x should not be wrapped.
    assertLatex(intMul, '\\int 2x \\, dx', 'Integral of Multiplication');


    console.log(`\nTests Completed: ${passed} Passed, ${failed} Failed`);
    if (failed > 0) process.exit(1);
}

runTests();
