
const { CAS } = require('../js/cas');
const { Expr, Sym, Num, Call, Div, Pow, Add, Sub, Mul } = require('../js/expression');

const cas = new CAS();
const k = new Sym('k');
const n = new Sym('n');
const a = new Sym('a');

function assert(condition, message) {
    if (!condition) {
        throw new Error("Test failed: " + message);
    }
    console.log("PASS: " + message);
}

function assertStr(expr, expected, message) {
    const s = expr.toString();
    if (s !== expected) {
        console.error(`Expected: ${expected}`);
        console.error(`Actual:   ${s}`);
        throw new Error("Test failed: " + message);
    }
    console.log("PASS: " + message);
}

console.log("--- Testing HEN-91: Sum.simplify Div summands and telescoping series ---");

// Bug #4: sum(a/k, k, 1, n) should factor out constant numerator
console.log("\n--- Bug #4: Div summands ---");

const bug4 = cas.evaluate(new Call('sum', [
    new Div(a, k), k, new Num(1), n
]));
console.log("sum(a/k, k, 1, n) => " + bug4.toString());
assert(bug4.toString() === "(a * sum((1 / k), k, 1, n))",
    "sum(a/k, k, 1, n) factors out a");

const bug4b = cas.evaluate(new Call('sum', [
    new Div(new Num(5), k), k, new Num(1), n
]));
console.log("sum(5/k, k, 1, n) => " + bug4b.toString());
assert(bug4b.toString() === "(5 * sum((1 / k), k, 1, n))",
    "sum(5/k, k, 1, n) factors out 5");

const bug4c = cas.evaluate(new Call('sum', [
    new Div(new Num(3), new Add(k, new Num(1))),
    k, new Num(1), n
]));
console.log("sum(3/(k+1), k, 1, n) => " + bug4c.toString());
assert(bug4c.toString() === "(3 * sum((1 / (k + 1)), k, 1, n))",
    "sum(3/(k+1), k, 1, n) factors out 3");

// Bug #5: telescoping series
console.log("\n--- Bug #5: Telescoping series ---");

const bug5 = cas.evaluate(new Call('sum', [
    new Div(new Num(1), new Mul(k, new Add(k, new Num(1)))),
    k, new Num(1), n
]));
console.log("sum(1/(k*(k+1)), k, 1, n) => " + bug5.toString());
// Should be 1 - 1/(n+1) which equals n/(n+1)
const bug5str = bug5.toString();
assert(bug5str.includes("1 - (1 / (n + 1))") || bug5str.includes("n / (n + 1)") || bug5str.includes("(n / (n + 1))"),
    "sum(1/(k*(k+1)), k, 1, n) simplifies to n/(n+1)");

// Verify numerically
const numTest = cas.evaluate(new Call('sum', [
    new Div(new Num(1), new Mul(k, new Add(k, new Num(1)))),
    k, new Num(1), new Num(5)
]));
console.log("sum(1/(k*(k+1)), k, 1, 5) => " + numTest.toString());
assert(numTest.toString() === "(5 / 6)",
    "sum(1/(k*(k+1)), k, 1, 5) = 5/6");

// Additional telescoping tests
const test2 = cas.evaluate(new Call('sum', [
    new Div(new Num(1), new Mul(k, new Add(k, new Num(2)))),
    k, new Num(1), n
]));
console.log("sum(1/(k*(k+2)), k, 1, n) => " + test2.toString());
assert(test2.toString() === "((1 - (1 / (n + 2))) / 2)",
    "sum(1/(k*(k+2)), k, 1, n) = (1 - 1/(n+2))/2");

const test3 = cas.evaluate(new Call('sum', [
    new Div(new Num(2), new Mul(k, new Add(k, new Num(1)))),
    k, new Num(1), n
]));
console.log("sum(2/(k*(k+1)), k, 1, n) => " + test3.toString());
assert(test3.toString() === "(2 * (1 - (1 / (n + 1))))",
    "sum(2/(k*(k+1)), k, 1, n) = 2*(1 - 1/(n+1))");

// Regression tests
console.log("\n--- Regression tests ---");

const sumk = cas.evaluate(new Call('sum', [k, k, new Num(1), n]));
console.log("sum(k, k, 1, n) => " + sumk.toString());
assert(sumk.toString() === "((n * (n + 1)) / 2)",
    "sum(k, k, 1, n) = n(n+1)/2");

const sumk2 = cas.evaluate(new Call('sum', [new Pow(k, new Num(2)), k, new Num(1), n]));
console.log("sum(k^2, k, 1, n) => " + sumk2.toString());
assert(!sumk2.toString().includes("sum("),
    "sum(k^2, k, 1, n) evaluates (not unevaluated)");

const sum2k = cas.evaluate(new Call('sum', [new Mul(new Num(2), k), k, new Num(1), n]));
console.log("sum(2*k, k, 1, n) => " + sum2k.toString());
assert(sum2k.toString() === "(n * (n + 1))",
    "sum(2*k, k, 1, n) = n(n+1)");

const sumConst = cas.evaluate(new Call('sum', [new Num(3), k, new Num(1), n]));
console.log("sum(3, k, 1, n) => " + sumConst.toString());
assert(sumConst.toString() === "(3 * n)",
    "sum(3, k, 1, n) = 3n");

// Check sum(1/k, k, 1, n) stays unevaluated (harmonic number, no closed form)
const sumHarmonic = cas.evaluate(new Call('sum', [
    new Div(new Num(1), k), k, new Num(1), n
]));
console.log("sum(1/k, k, 1, n) => " + sumHarmonic.toString());
assert(sumHarmonic.toString().includes("sum("),
    "sum(1/k, k, 1, n) stays unevaluated (no closed form for harmonic number)");

console.log("\n--- All tests passed ---");
