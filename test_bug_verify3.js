const { CAS } = require('./js/cas.js');
const { Num, Sym, Call, Vec, Add, Sub, Mul, Div, Pow, Eq, Lt } = require('./js/expression.js');
const cas = new CAS();
globalThis.HELP_DATA = {};
const x = new Sym('x'), y = new Sym('y'), z = new Sym('z');

function cmd(name, ...args) { return cas.evaluate(new Call(name, args)); }

// Build x^5 - x^4 + x^3 - x^2 + x - 1 = 0 programmatically
// = ((x-1)*(x^4+x^2+1))
function buildPoly() {
  const x2 = new Pow(x, new Num(2));
  const x3 = new Pow(x, new Num(3));
  const x4 = new Pow(x, new Num(4));
  const x5 = new Pow(x, new Num(5));
  // x^5 - x^4 + x^3 - x^2 + x - 1
  return new Sub(new Sub(new Sub(new Sub(x5, x4), x3), x2), new Add(x, new Num(-1)));
}

console.log('\n== SOLVE: quintic that triggers quartic handler ==');
try {
  const result = cmd('solve', buildPoly(), x);
  console.log(`  solve(x^5-x^4+x^3-x^2+x-1, x): ${result.toString()}`);
} catch(e) {
  console.log(`  ERROR: ${e.message}`);
}

console.log('\n== SOLVE: simpler quartic (known good) ==');
try {
  // x^4 - 5*x^2 + 4 = 0 => (x^2-1)(x^2-4) => x = ±1, ±2
  const poly = new Sub(new Sub(new Pow(x, new Num(4)), new Mul(new Num(5), new Pow(x, new Num(2)))), new Num(4));
  const result = cmd('solve', poly, x);
  console.log(`  solve(x^4-5x^2+4, x): ${result.toString()}`);
} catch(e) {
  console.log(`  ERROR: ${e.message}`);
}

console.log('\n== SOLVE: another quintic ==');
try {
  // x^5 - 1 = 0 => 5th roots of unity
  const poly = new Sub(new Pow(x, new Num(5)), new Num(1));
  const result = cmd('solve', poly, x);
  console.log(`  solve(x^5-1, x): ${result.toString()}`);
} catch(e) {
  console.log(`  ERROR: ${e.message}`);
}

console.log('\n== SOLVE: sextic ==');
try {
  // (x^2-1)(x^4-1) = x^6 - x^4 - x^2 + 1 = 0 => x = ±1, ±1, ±1
  const poly = new Sub(new Sub(new Pow(x, new Num(6)), new Pow(x, new Num(4))), new Add(new Pow(x, new Num(2)), new Num(-1)));
  const result = cmd('solve', poly, x);
  console.log(`  solve(x^6-x^4-x^2+1, x): ${result.toString()}`);
} catch(e) {
  console.log(`  ERROR: ${e.message}`);
}

console.log('\n== SOLVE: another sextic (x^6+1) ==');
try {
  const poly = new Add(new Pow(x, new Num(6)), new Num(1));
  const result = cmd('solve', poly, x);
  console.log(`  solve(x^6+1, x): ${result.toString()}`);
} catch(e) {
  console.log(`  ERROR: ${e.message}`);
}

console.log('\n== SOLVE: sum(1/(k*(k+1)), k, 1, n) telescoping ==');
try {
  const k = new Sym('k');
  const n = new Sym('n');
  const telescoping = new Div(new Num(1), new Mul(k, new Add(k, new Num(1))));
  const result = cmd('sum', telescoping, k, new Num(1), n);
  console.log(`  sum(1/(k*(k+1)), k, 1, n): ${result.toString()}`);
} catch(e) {
  console.log(`  ERROR: ${e.message}`);
}

console.log('\n== INTEGRATE: piecewise definite integral ==');
try {
  const pw = new Call('piecewise', [
    new Lt(x, new Num(0)),
    new Num(0),
    new Pow(x, new Num(2))
  ]);
  const result = cmd('integrate', pw, x, new Num(-1), new Num(1));
  console.log(`  integrate(piecewise(x<0,0,x^2), x, -1, 1): ${result.toString()}`);
} catch(e) {
  console.log(`  ERROR: ${e.message}`);
}

console.log('\n== SOLVE: abs equation with negative RHS (BUG) ==');
try {
  const result = cmd('solve', new Eq(new Call('abs', [x]), new Num(-3)), x);
  console.log(`  solve(abs(x) = -3, x): ${result.toString()}`);
} catch(e) {
  console.log(`  ERROR: ${e.message}`);
}

console.log('\n== LIMIT: 0*Inf form x*(1/x^2) at 0+ ==');
try {
  const result = cmd('limit', new Mul(x, new Div(new Num(1), new Pow(x, new Num(2)))), x, new Num(0), new Num(1));
  console.log(`  limit(x/x^2, x, 0, right): ${result.toString()}`);
} catch(e) {
  console.log(`  ERROR: ${e.message}`);
}

console.log('\n== LIMIT: inf - inf with Div combination ==');
try {
  const result = cmd('limit', new Sub(new Div(new Num(1), x), new Div(new Num(1), new Sub(x, new Num(1)))), x, new Sym('Infinity'));
  console.log(`  limit(1/x - 1/(x-1), x, inf): ${result.toString()}`);
} catch(e) {
  console.log(`  ERROR: ${e.message}`);
}

console.log('\n== SERIES: sqrt(1+x) Laurent at x=-1 (pole) ==');
try {
  const result = cmd('series', new Call('sqrt', [new Add(new Num(1), x)]), x, new Num(-1), new Num(3));
  console.log(`  series(sqrt(1+x), x, -1, 3): ${result.toString()}`);
} catch(e) {
  console.log(`  ERROR: ${e.message}`);
}

console.log('\n== DONE ==');
