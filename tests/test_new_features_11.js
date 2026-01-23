const { CAS } = require('../js/cas');
const { Expr, Num, Sym, Add, Sub, Mul, Div, Pow, Call, Vec } = require('../js/expression');

const cas = new CAS();

console.log("--- Testing Geometric Sum ---");
// sum(2^k, k, 0, n)
// Ratio = 2^(k+1)/2^k = 2. Independent of k.
// Sum = 2^0 * (2^(n+1) - 1) / (2 - 1) = 2^(n+1) - 1.
const k = new Sym('k');
const n = new Sym('n');
const sumGeo = new Call('sum', [new Pow(new Num(2), k), k, new Num(0), n]);
const resGeo = cas.evaluate(sumGeo);
console.log(`sum(2^k, k, 0, n) = ${resGeo.toString()}`);

// sum(r^k, k, 0, n)
// Ratio = r.
// Sum = r^0 * (r^(n+1) - 1) / (r - 1) = (r^(n+1) - 1) / (r - 1).
const r = new Sym('r');
const sumGeoSym = new Call('sum', [new Pow(r, k), k, new Num(0), n]);
const resGeoSym = cas.evaluate(sumGeoSym);
console.log(`sum(r^k, k, 0, n) = ${resGeoSym.toString()}`);


console.log("\n--- Testing Elliptic Integrals ---");
// Numeric K(0.5)
const kVal = new Num(0.5);
const K = cas.evaluate(new Call('EllipticK', [kVal]));
console.log(`EllipticK(0.5) = ${K.toString()}`);
// Check value
const expectedK = 1.6857503548;
if (Math.abs(K.value - expectedK) < 1e-6) console.log("PASS: EllipticK(0.5) numeric");
else console.log("FAIL: EllipticK(0.5) numeric");

// Numeric E(0.5)
const E = cas.evaluate(new Call('EllipticE', [kVal]));
console.log(`EllipticE(0.5) = ${E.toString()}`);
// Expected: 1.4674622093
const expectedE = 1.4674622093;
if (Math.abs(E.value - expectedE) < 1e-6) console.log("PASS: EllipticE(0.5) numeric");
else console.log("FAIL: EllipticE(0.5) numeric");

// Derivative of K(k)
const diffK = cas.evaluate(new Call('diff', [new Call('EllipticK', [k]), k]));
console.log(`diff(K(k), k) = ${diffK.toString()}`);

console.log("\n--- Testing Interpolation ---");
// interpolate([[0,1], [1,3], [2,7]], x)
const points = new Vec([
    new Vec([new Num(0), new Num(1)]),
    new Vec([new Num(1), new Num(3)]),
    new Vec([new Num(2), new Num(7)])
]);
const x = new Sym('x');
const poly = cas.evaluate(new Call('interpolate', [points, x]));
console.log(`interpolate([[0,1],[1,3],[2,7]], x) = ${poly.toString()}`);

const polyExpanded = cas.evaluate(new Call('expand', [poly]));
console.log(`Expanded: ${polyExpanded.toString()}`);

// Verify values
const val0 = cas.evaluate(poly.substitute(x, new Num(0)));
const val1 = cas.evaluate(poly.substitute(x, new Num(1)));
const val2 = cas.evaluate(poly.substitute(x, new Num(2)));

console.log(`P(0)=${val0}, P(1)=${val1}, P(2)=${val2}`);

if (val0.value === 1 && val1.value === 3 && val2.value === 7) {
    console.log("PASS: Interpolation values correct");
} else {
    console.log("FAIL: Interpolation values incorrect");
}
