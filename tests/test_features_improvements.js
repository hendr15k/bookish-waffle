
const fs = require('fs');
const path = require('path');
const vm = require('vm');

// Load files manually to simulate browser environment
const expressionCode = fs.readFileSync(path.join(__dirname, '../js/expression.js'), 'utf8');
const casCode = fs.readFileSync(path.join(__dirname, '../js/cas.js'), 'utf8');

const context = {
    console: console,
    Math: Math,
    Number: Number,
    parseFloat: parseFloat,
    parseInt: parseInt,
    isNaN: isNaN,
    isFinite: isFinite,
    module: { exports: {} },
    exports: {}
};
vm.createContext(context);

// Run expression.js
vm.runInContext(expressionCode, context);
vm.runInContext(casCode, context);

const { CAS, Call, Add, Sub, Mul, Div, Pow, Eq, Vec, Num, Sym } = context;

const cas = new CAS();

// 1. Test Cubic Solve (Irrational): x^3 - 3x + 1 = 0
const x = new Sym('x');
const cubicIrrational = new Add(
    new Sub(
        new Pow(x, new Num(3)),
        new Mul(new Num(3), x)
    ),
    new Num(1)
);

console.log("Cubic Poly (Irrational):", cubicIrrational.toString());
const rootsIrr = cas.evaluate(new Call('solve', [new Eq(cubicIrrational, new Num(0)), x]));
console.log("Roots found (Irrational):", rootsIrr.toString());

// 2. Test Limit: sin(x)/x at 0
const limExpr = new Div(new Call('sin', [x]), x);
const limitVal = cas.evaluate(new Call('limit', [limExpr, x, new Num(0)]));
console.log("Limit sin(x)/x at 0:", limitVal.toString());

// 3. Test Rank
const mat = new Vec([
    new Vec([new Num(1), new Num(2), new Num(3)]),
    new Vec([new Num(4), new Num(5), new Num(6)]),
    new Vec([new Num(7), new Num(8), new Num(9)])
]);
const rankVal = cas.evaluate(new Call('rank', [mat]));
console.log("Rank of [[1,2,3],[4,5,6],[7,8,9]]:", rankVal.toString());

// 4. Test Inverse Laplace: 1/(s^2 + 1) -> sin(t)
const s = new Sym('s');
const t = new Sym('t');
const lapExpr = new Div(new Num(1), new Add(new Pow(s, new Num(2)), new Num(1)));
const ilap = cas.evaluate(new Call('ilaplace', [lapExpr, s, t]));
console.log("Inverse Laplace of 1/(s^2+1):", ilap.toString());

// Check expectations
let pass = true;
// cubic roots should use cos for casus irreducibilis (3 real roots)
if (!rootsIrr.toString().includes('cos(')) {
    console.log("FAIL: Cubic roots should use trigonometric solution (cos) for 3 real roots");
    pass = false;
}
if (limitVal.toString() !== '1') {
    console.log("FAIL: Limit incorrect");
    pass = false;
}
if (rankVal.toString() !== '2') {
    console.log("FAIL: Rank incorrect");
    pass = false;
}
if (!ilap.toString().includes('sin(t)')) {
    console.log("FAIL: Inverse Laplace incorrect");
    pass = false;
}

if (pass) console.log("ALL TESTS PASSED");
else console.log("SOME TESTS FAILED");
