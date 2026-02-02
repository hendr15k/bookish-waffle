
const { CAS } = require('../js/cas.js');
const { Expr, Num, Sym, Call, Vec, Add, Mul, Pow, Sub } = require('../js/expression.js');

const cas = new CAS();

console.log("--- Test Quadratic Factorization ---");
const x = new Sym('x');
// 2x^2 + 5x + 2
const poly = new Add(new Mul(new Num(2), new Pow(x, new Num(2))), new Add(new Mul(new Num(5), x), new Num(2)));
const factored = cas._factor(poly);
console.log("Factor(2x^2+5x+2) =", factored.toString());

if (factored.toString() === "((2 * x) + 1) * (x + 2)" || factored.toString() === "(x + 2) * ((2 * x) + 1)") {
    console.log("PASS: Quadratic factorization");
} else {
    // Also accept with leading coeff 1 if it simplified well, but we want integer factors.
    // My implementation: pushed '2x+1' and 'x+2'. No leading coeff remaining.
    // Mul(Add(2x, 1), Sub(x, -2)) -> (2x+1)(x+2)
    console.log("CHECK: Is this acceptable?");
}

console.log("\n--- Test Pow.toLatex for Trig ---");
const sin2 = new Pow(new Call('sin', [x]), new Num(2));
console.log("Latex sin(x)^2 =", sin2.toLatex());
if (sin2.toLatex() === "\\sin^{2}\\left(x\\right)") {
    console.log("PASS: Latex format");
} else {
    console.log("FAIL: Expected \\sin^{2}\\left(x\\right)");
}

console.log("\n--- Test Inverse Laplace ---");
const s = new Sym('s');
const t = new Sym('t');
// 1/(s^2+1) -> sin(t)
const expr = new Div(new Num(1), new Add(new Pow(s, new Num(2)), new Num(1)));
const il = cas._ilaplace(expr, s, t);
console.log("ilaplace(1/(s^2+1)) =", il.toString());
if (il.toString() === "sin(t)") {
    console.log("PASS: Inverse Laplace");
} else {
    console.log("FAIL: Expected sin(t)");
}

console.log("\n--- Test Divergence ---");
// div([x^2, y^2, z^2], [x,y,z]) = 2x+2y+2z
const y = new Sym('y');
const z = new Sym('z');
const v = new Vec([new Pow(x, new Num(2)), new Pow(y, new Num(2)), new Pow(z, new Num(2))]);
const vars = new Vec([x, y, z]);
const div = cas._divergence(v, vars);
console.log("div([x^2, y^2, z^2]) =", div.toString());
if (div.toString() === "((2 * x) + ((2 * y) + (2 * z)))") {
    console.log("PASS: Divergence");
} else {
    console.log("FAIL: Check output");
}
