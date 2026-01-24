const { CAS } = require('../js/cas');
const { Expr, Num, Sym, Add, Sub, Mul, Div, Pow, Call, Vec } = require('../js/expression');

const cas = new CAS();

console.log("--- Testing Numeric Bessel Functions ---");
// J_0(0) = 1
const J0 = cas.evaluate(new Call('besselJ', [new Num(0), new Num(0)]));
console.log(`besselJ(0, 0) = ${J0.toString()}`);
if (Math.abs(J0.value - 1) < 1e-9) console.log("PASS: besselJ(0,0)");
else console.log("FAIL: besselJ(0,0)");

// J_1(1) ~ 0.44005
const J1 = cas.evaluate(new Call('besselJ', [new Num(1), new Num(1)]));
console.log(`besselJ(1, 1) = ${J1.toString()}`);
if (Math.abs(J1.value - 0.4400505857) < 1e-6) console.log("PASS: besselJ(1,1)");
else console.log("FAIL: besselJ(1,1)");

// Y_0(1) ~ 0.088256
const Y0 = cas.evaluate(new Call('besselY', [new Num(0), new Num(1)]));
console.log(`besselY(0, 1) = ${Y0.toString()}`);
if (Math.abs(Y0.value - 0.08825696) < 1e-5) console.log("PASS: besselY(0,1)");
else console.log("FAIL: besselY(0,1)");

// I_0(1) ~ 1.26606
const I0 = cas.evaluate(new Call('besselI', [new Num(0), new Num(1)]));
console.log(`besselI(0, 1) = ${I0.toString()}`);
if (Math.abs(I0.value - 1.2660658778) < 1e-5) console.log("PASS: besselI(0,1)");
else console.log("FAIL: besselI(0,1)");

console.log("\n--- Testing Numeric Hypergeometric hyp2f1 ---");
// hyp2f1(1, 1, 2, 0.5) = ln(1-0.5)/(-0.5) = ln(0.5)/-0.5 = -0.693/-0.5 = 1.38629
const H1 = cas.evaluate(new Call('hyp2f1', [new Num(1), new Num(1), new Num(2), new Num(0.5)]));
console.log(`hyp2f1(1,1,2,0.5) = ${H1.toString()}`);
const expectedH1 = Math.log(0.5) / -0.5;
if (Math.abs(H1.value - expectedH1) < 1e-5) console.log("PASS: hyp2f1(1,1,2,0.5)");
else console.log("FAIL: hyp2f1(1,1,2,0.5)");

console.log("\n--- Testing Symbolic hyp2f1 ---");
const z = new Sym('z');
const diffH = cas.evaluate(new Call('diff', [new Call('hyp2f1', [new Num(1), new Num(2), new Num(3), z]), z]));
console.log(`diff(hyp2f1(1,2,3,z), z) = ${diffH.toString()}`);
// Expected: 1*2/3 * hyp2f1(2,3,4,z)
// Result might be simplified float coeffs
if (diffH.toString().includes('hyp2f1(2, 3, 4, z)')) console.log("PASS: diff(hyp2f1)");
else console.log("FAIL: diff(hyp2f1)");

console.log("\n--- Testing Rational Simplification (polyGcd) ---");
const x = new Sym('x');
// (x^2 - 1) / (x + 1) -> x - 1
// x^2 - 1 = (x+1)(x-1)
const num1 = new Sub(new Pow(x, new Num(2)), new Num(1));
const den1 = new Add(x, new Num(1));
const rat1 = cas.evaluate(new Div(num1, den1)); // simplify() called inside evaluate
console.log(`(x^2 - 1) / (x + 1) = ${rat1.toString()}`);
// Should be x - 1
const valRat1 = rat1.substitute(x, new Num(5)).evaluateNumeric();
if (Math.abs(valRat1 - 4) < 1e-9 && rat1.toString().length < 15) console.log("PASS: (x^2-1)/(x+1) simplified");
else console.log("FAIL: (x^2-1)/(x+1)");

// (x^2 + 2x + 1) / (x + 1) -> x + 1
const num2 = new Add(new Add(new Pow(x, new Num(2)), new Mul(new Num(2), x)), new Num(1));
const den2 = new Add(x, new Num(1));
const rat2 = cas.evaluate(new Div(num2, den2));
console.log(`(x^2 + 2x + 1) / (x + 1) = ${rat2.toString()}`);
const valRat2 = rat2.substitute(x, new Num(5)).evaluateNumeric();
if (Math.abs(valRat2 - 6) < 1e-9 && rat2.toString().length < 15) console.log("PASS: (x^2+2x+1)/(x+1) simplified");
else console.log("FAIL: (x^2+2x+1)/(x+1)");

// (2x + 2) / 2 -> x + 1
const num3 = new Add(new Mul(new Num(2), x), new Num(2));
const den3 = new Num(2);
const rat3 = cas.evaluate(new Div(num3, den3));
console.log(`(2x + 2) / 2 = ${rat3.toString()}`);
if (rat3.toString().includes('x') && !rat3.toString().includes('/')) console.log("PASS: (2x+2)/2 simplified");
else console.log("FAIL: (2x+2)/2");

// (x^3 - 1) / (x - 1) -> x^2 + x + 1
const num4 = new Sub(new Pow(x, new Num(3)), new Num(1));
const den4 = new Sub(x, new Num(1));
const rat4 = cas.evaluate(new Div(num4, den4));
console.log(`(x^3 - 1) / (x - 1) = ${rat4.toString()}`);
const valRat4 = rat4.substitute(x, new Num(2)).evaluateNumeric();
// 7/1 = 7. 4+2+1 = 7.
if (Math.abs(valRat4 - 7) < 1e-9) console.log("PASS: (x^3-1)/(x-1) simplified");
else console.log("FAIL: (x^3-1)/(x-1)");
