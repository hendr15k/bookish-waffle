const { CAS } = require('./js/cas.js');
const { Num, Sym, Call, Vec, Add, Sub, Mul, Div, Pow, Eq } = require('./js/expression.js');
const cas = new CAS();
globalThis.HELP_DATA = {};
const x = new Sym('x');

function cmd(name, ...args) { return cas.evaluate(new Call(name, args)); }

// Test various polynomials that might trigger the quartic handler crash
console.log('\n== Quartic with complex roots (triggers casus irreducibilis) ==');
try {
  // x^4 + x^2 + 1 = 0 -> factors as (x^2+x+1)(x^2-x+1), roots: e^(±iπ/3), e^(±2iπ/3)
  const poly = new Add(new Pow(x, new Num(4)), new Add(new Pow(x, new Num(2)), new Num(1)));
  console.log(`  solve(x^4+x^2+1, x): ${cmd('solve', poly, x).toString()}`);
} catch(e) { console.log(`  ERROR: ${e.message}`); }

console.log('\n== Quartic that triggers Ferrari method ==');
try {
  // x^4 - 4*x^3 + 6*x^2 - 4*x + 1 = (x-1)^4
  const poly = new Sub(new Add(new Add(new Pow(x, new Num(4)), new Mul(new Num(-4), new Pow(x, new Num(3)))), new Mul(new Num(6), new Pow(x, new Num(2)))), new Mul(new Num(-4), x)), new Num(1));
  console.log(`  solve((x-1)^4, x): ${cmd('solve', poly, x).toString()}`);
} catch(e) { console.log(`  ERROR: ${e.message}`); }

console.log('\n== Quartic with biquadratic resolvent (q=0 path) ==');
try {
  // x^4 - 5*x^2 + 4 = 0 -> factors as (x^2-1)(x^2-4)
  const poly = new Sub(new Sub(new Pow(x, new Num(4)), new Mul(new Num(5), new Pow(x, new Num(2)))), new Num(4));
  console.log(`  solve(x^4-5x^2+4, x): ${cmd('solve', poly, x).toString()}`);
} catch(e) { console.log(`  ERROR: ${e.message}`); }

console.log('\n== Quartic requiring Ferrari (not biquadratic) ==');
try {
  // x^4 + x^3 - x - 1 = (x^2+1)(x^2+x-1)
  const poly = new Sub(new Add(new Pow(x, new Num(4)), new Pow(x, new Num(3))), new Add(x, new Num(1)));
  console.log(`  solve(x^4+x^3-x-1, x): ${cmd('solve', poly, x).toString()}`);
} catch(e) { console.log(`  ERROR: ${e.message}`); }

console.log('\n== Pure quintic (cyclotomic x^5-1) ==');
try {
  const poly = new Sub(new Pow(x, new Num(5)), new Num(1));
  console.log(`  solve(x^5-1, x): ${cmd('solve', poly, x).toString()}`);
} catch(e) { console.log(`  ERROR: ${e.message}`); }

console.log('\n== Quintic (x-1)(x^2+x+1)(x^2-x+1) ==');
try {
  // (x-1)(x^4+x^2+1)
  const poly = new Mul(new Sub(x, new Num(1)), new Add(new Pow(x, new Num(4)), new Add(new Pow(x, new Num(2)), new Num(1))));
  console.log(`  solve((x-1)(x^4+x^2+1), x): ${cmd('solve', poly, x).toString()}`);
} catch(e) { console.log(`  ERROR: ${e.message}`); }

console.log('\n== High-degree that falls through to cubic handler ==');
try {
  // x^7 - 1
  const poly = new Sub(new Pow(x, new Num(7)), new Num(1));
  console.log(`  solve(x^7-1, x): ${cmd('solve', poly, x).toString()}`);
} catch(e) { console.log(`  ERROR: ${e.message}`); }

console.log('\n== Quartic irreducible: x^4 + 4*x^2 + 5 ==');
try {
  const poly = new Add(new Add(new Pow(x, new Num(4)), new Mul(new Num(4), new Pow(x, new Num(2)))), new Num(5));
  console.log(`  solve(x^4+4x^2+5, x): ${cmd('solve', poly, x).toString()}`);
} catch(e) { console.log(`  ERROR: ${e.message}`); }

console.log('\n== x^4 + x + 1 (no cubic resolvent path) ==');
try {
  const poly = new Add(new Add(new Pow(x, new Num(4)), x), new Num(1));
  console.log(`  solve(x^4+x+1, x): ${cmd('solve', poly, x).toString()}`);
} catch(e) { console.log(`  ERROR: ${e.message}`); }

console.log('\n== Factor on irreducible quartic ==');
try {
  const poly = new Add(new Add(new Pow(x, new Num(4)), new Pow(x, new Num(2))), new Num(1));
  console.log(`  factor(x^4+x^2+1): ${cmd('factor', poly).toString()}`);
} catch(e) { console.log(`  ERROR: ${e.message}`); }

console.log('\n== DONE ==');
