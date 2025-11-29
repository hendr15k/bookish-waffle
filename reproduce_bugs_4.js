
const { CAS } = require('./js/cas.js');
const { Expr, Num, Sym, Add, Sub, Mul, Div, Pow, Call, Vec } = require('./js/expression.js');

const cas = new CAS();

function assert(condition, message) {
    if (!condition) {
        console.error(`FAIL: ${message}`);
        process.exit(1);
    } else {
        console.log(`PASS: ${message}`);
    }
}

// 1. Matrix Power
// A = [[1, 0], [0, 1]] (Identity)
// A^2 should be A
const I = new Vec([new Vec([new Num(1), new Num(0)]), new Vec([new Num(0), new Num(1)])]);
const I2 = new Pow(I, new Num(2)).simplify();
console.log("I^2:", I2.toString());
// Expected: [[1, 0], [0, 1]]
// Current implementation of Pow.simplify: `if (l instanceof Num && r instanceof Num) ... return new Pow(l, r)`
// It does NOT handle Matrix base.

// 2. Limit at Infinity
// limit(1/x, x, Infinity)
// Should be 0.
const expr = new Div(new Num(1), new Sym('x'));
const lim = cas._limit(expr, new Sym('x'), new Sym('Infinity'));
console.log("limit(1/x, x, Infinity):", lim.toString());
// Current `_limit` substitutes. `1/Infinity`.
// `Div.simplify` doesn't handle `Infinity` explicitly except maybe Num?
// `Sym` evaluateNumeric returns NaN for Infinity?
// Let's see what happens.

// 3. Limit x -> Infinity of x/(x+1)
const expr2 = new Div(new Sym('x'), new Add(new Sym('x'), new Num(1)));
const lim2 = cas._limit(expr2, new Sym('x'), new Sym('Infinity'));
console.log("limit(x/(x+1), x, Infinity):", lim2.toString());
// Expected: 1
// Current: Infinity / Infinity -> L'Hopital?
// `diff(Infinity)` ? No, `diff` is on expression.
// `substitute` returns `Infinity` symbols.
// `simplify` might not reduce `Infinity/Infinity`.

// 4. Matrix Determinant 3x3
// [[1,0,0],[0,1,0],[0,0,1]] -> 1
const I3 = new Vec([
    new Vec([new Num(1), new Num(0), new Num(0)]),
    new Vec([new Num(0), new Num(1), new Num(0)]),
    new Vec([new Num(0), new Num(0), new Num(1)])
]);
const det = cas._det(I3);
console.log("det(I3):", det.toString());
// Recursion should handle this if correct.
