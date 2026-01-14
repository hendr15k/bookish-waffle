
const { Parser, Lexer } = require('../js/parser');
const { CAS } = require('../js/cas');
const { Num, Sym, Vec, Call, Add, Mul, Div, Pow } = require('../js/expression'); // Pre-load classes

// Polyfill for CAS environment
global.Num = Num;
global.Sym = Sym;
global.Vec = Vec;
global.Call = Call;
global.Add = Add;
global.Mul = Mul;
global.Div = Div;
global.Pow = Pow;

const cas = new CAS();

function test(input, expected) {
    const parser = new Parser(new Lexer(input));
    const expr = parser.parse();
    const result = cas.evaluate(expr);
    const resStr = result.toString();
    if (resStr !== expected) {
        console.error(`FAILED: ${input} => ${resStr} (Expected: ${expected})`);
    } else {
        console.log(`PASSED: ${input} => ${resStr}`);
    }
}

function testNumeric(input, expectedVal) {
    const parser = new Parser(new Lexer(input));
    const expr = parser.parse();
    const result = cas.evaluate(expr);
    const resVal = result.evaluateNumeric();
    if (Math.abs(resVal - expectedVal) > 1e-5) {
        console.error(`FAILED: ${input} => ${resVal} (Expected: ${expectedVal})`);
    } else {
        console.log(`PASSED: ${input} => ${resVal}`);
    }
}

console.log("--- Testing CRT ---");
// x = 2 mod 3, x = 3 mod 5 -> x = 8
// crt([2, 3], [3, 5])
test("crt([2, 3], [3, 5])", "8");
// x = 1 mod 2, x = 2 mod 3, x = 3 mod 5 -> 23
test("crt([1, 2, 3], [2, 3, 5])", "23");

console.log("--- Testing Lagrange ---");
// Points: (1, 1), (2, 4), (3, 9) -> x^2
// lagrange([[1,1], [2,4], [3,9]])
// Result might be expanded or not: 1 * ... + 4 * ...
// Simplify should handle it.
// Expected: x^2 (symbolic simplification might be tricky)
// Let's check value at 4 -> 16
let parser = new Parser(new Lexer("lagrange([[1,1], [2,4], [3,9]])"));
let poly = cas.evaluate(parser.parse());
console.log("Lagrange Poly:", poly.toString());
let valAt4 = poly.substitute(new Sym('x'), new Num(4)).simplify().evaluateNumeric();
if (Math.abs(valAt4 - 16) < 1e-9) console.log("PASSED: Lagrange value at 4 is 16");
else console.error("FAILED: Lagrange value at 4 is", valAt4);

console.log("--- Testing RK4 ---");
// y' = y, y(0) = 1. Solution y = e^t. At t=1, y=e.
// rk4(y, t, y, 0, 1, 0.1, 10)
// Last point should be approx (1, 2.718)
parser = new Parser(new Lexer("rk4(y, t, y, 0, 1, 0.1, 10)"));
let rkRes = cas.evaluate(parser.parse());
let lastPt = rkRes.elements[rkRes.elements.length - 1]; // [t, y]
let yEnd = lastPt.elements[1].evaluateNumeric();
if (Math.abs(yEnd - Math.E) < 1e-4) console.log("PASSED: RK4 y(1) approx e");
else console.error("FAILED: RK4 y(1) =", yEnd);

console.log("--- Testing LSQ ---");
// A = [[1, 1], [1, 2], [1, 3]]
// b = [6, 5, 7]
// x = [y-intercept, slope]
// lsq(A, b)
// x = [[1, 1], [1, 2], [1, 3]]
// y = [6, 5, 7]
// This is linear regression on (1,6), (2,5), (3,7).
// slope = 0.5, intercept = 5. (Wait, let's calculate)
// Mean x=2, Mean y=6.
// Sxy = (1-2)(6-6) + (2-2)(5-6) + (3-2)(7-6) = 0 + 0 + 1 = 1
// Sxx = 1 + 0 + 1 = 2
// Slope = 0.5.
// Intercept = 6 - 0.5*2 = 5.
// Result should be vector [5, 0.5] (col vector usually)
parser = new Parser(new Lexer("lsq([[1, 1], [1, 2], [1, 3]], [[6], [5], [7]])")); // b as col vector
let lsqRes = cas.evaluate(parser.parse());
console.log("LSQ Result:", lsqRes.toString());
// Expect [[5], [0.5]]
// Check numerical
let resVec = lsqRes.elements; // [[5], [0.5]]
if (Math.abs(resVec[0].elements[0].evaluateNumeric() - 5) < 1e-9 &&
    Math.abs(resVec[1].elements[0].evaluateNumeric() - 0.5) < 1e-9) {
    console.log("PASSED: LSQ");
} else {
    console.error("FAILED: LSQ result incorrect");
}

console.log("--- Testing Integrate Erf ---");
// integrate(erf(x), x)
test("integrate(erf(x), x)", "((x * erf(x)) + (exp((-1 * x^2)) / sqrt(pi)))");
