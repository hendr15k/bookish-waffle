
// Load expression classes globally
require('../js/expression.js');

// Load CAS
const { CAS } = require('../js/cas.js');

const cas = new CAS();

// Helpers for testing
const x = new Sym('x');
const y = new Sym('y');

function check(name, expr, expected) {
    try {
        const res = cas.evaluate(expr);
        const resStr = res.toString();
        console.log(`${name}: ${resStr}`);
        if (expected && resStr !== expected) {
            console.log(`FAIL: Expected ${expected}, got ${resStr}`);
        } else if (expected) {
            console.log(`PASS`);
        }
    } catch (e) {
        console.log(`${name}: Error - ${e.message}`);
    }
}

console.log("--- Feature 1: Radical Simplification ---");
// simplify(sqrt(8)) -> 2*sqrt(2)
check('simplify(sqrt(8))', new Call('sqrt', [new Num(8)]), '(2 * sqrt(2))');
check('simplify(sqrt(12))', new Call('sqrt', [new Num(12)]), '(2 * sqrt(3))');
check('simplify(sqrt(72))', new Call('sqrt', [new Num(72)]), '(6 * sqrt(2))');
// simplify(sqrt(2)) -> sqrt(2)
check('simplify(sqrt(2))', new Call('sqrt', [new Num(2)]), 'sqrt(2)');

console.log("\n--- Feature 2: Log/Exp Simplification ---");
// simplify(exp(2 * ln(x))) -> x^2
check('simplify(exp(2 * ln(x)))', new Call('exp', [new Mul(new Num(2), new Call('ln', [x]))]), 'x^2');
// simplify(exp(ln(x) * 2)) -> x^2
check('simplify(exp(ln(x) * 2))', new Call('exp', [new Mul(new Call('ln', [x]), new Num(2))]), 'x^2');

console.log("\n--- Feature 3: Negation Distribution ---");
// simplify(-1 * (x - 1)) -> 1 - x
check('simplify(-1 * (x - 1))', new Mul(new Num(-1), new Sub(x, new Num(1))), '(1 - x)');
// solve(sin(x) - 0.5, x) output check
// Previously: {(-1 * ((-1 * pi) + (pi / 6))), (pi / 6)} which is {-( -pi + pi/6 ), pi/6}
// Expected: {(pi - (pi / 6)), (pi / 6)} -> { (5pi/6), (pi/6) } approx.
// Note: pi - pi/6 simplifies? Sub(pi, pi/6).
// If simplify handles Sub(Num*pi, Num*pi), it should work.
check('solve(sin(x) - 0.5, x)', new Call('solve', [new Sub(new Call('sin', [x]), new Num(0.5)), x]));
