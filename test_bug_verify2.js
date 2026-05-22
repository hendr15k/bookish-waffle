const { CAS } = require('./js/cas.js');
const { Num, Sym, Call, Vec, Add, Sub, Mul, Div, Pow, Eq, Lt, Gt } = require('./js/expression.js');
const cas = new CAS();
globalThis.HELP_DATA = {};
const x = new Sym('x'), y = new Sym('y'), k = new Sym('k'), n = new Sym('n');
const num = v => new Num(v);
const cmd = (name, ...args) => cas.evaluate(new Call(name, args));

function test(name, fn) {
  try {
    const r = fn();
    console.log(`  ${name}: ${r && r.toString ? r.toString() : String(r)}`);
  } catch(e) {
    console.log(`  ${name}: ERROR ${e.message}`);
  }
}

console.log('\n== sum: div(a+b, k) ==');
test('sum((a+b)/k, k, 1, n)', () => cmd('sum', new Div(new Add(num(1), num(1)), k), k, num(1), n));
test('sum((x+k)/k, k, 1, n)', () => cmd('sum', new Div(new Add(x, k), k), k, num(1), n));

console.log('\n== sum: sum(1/(k*(k+1)), k, 1, n) ==');
const telescoping = new Div(num(1), new Mul(k, new Add(k, num(1))));
test('sum(1/(k*(k+1)), k, 1, n)', () => cmd('sum', telescoping, k, num(1), n));

console.log('\n== product: div symbolic ==');
test('product(1/k, k, 1, n)', () => cmd('product', new Div(num(1), k), k, num(1), n));
test('product((k+1)/k, k, 1, n)', () => cmd('product', new Div(new Add(k, num(1)), k), k, num(1), n));

console.log('\n== limit: 0 * Inf form ==');
test('limit(x*(1/x), x, 0, right)', () => cmd('limit', new Mul(x, new Div(num(1), x)), x, num(0), num(1)));
test('limit(x*(1/x^2), x, 0, right)', () => cmd('limit', new Mul(x, new Div(num(1), new Pow(x, num(2)))), x, num(0), num(1)));

console.log('\n== limit: 1^Inf form (0/0 form) ==');
test('limit((x)/x, x, 1)', () => cmd('limit', new Div(x, x), x, num(1)));

console.log('\n== limit: Inf - Inf ==');
test('limit(1/x - 1/(x-1), x, Infinity)', () => cmd('limit', new Sub(new Div(num(1), x), new Div(num(1), new Sub(x, num(1)))), x, new Sym('Infinity')));

console.log('\n== solve: very high degree ==');
test('solve(x^5 - x^4 + x^3 - x^2 + x - 1, x)', () => cmd('solve', new Sub(new Sub(new Sub(new Sub(new Pow(x, num(5)), new Pow(x, num(4))), new Pow(x, num(3))), new Pow(x, num(2))), x), x).toString());

console.log('\n== solve: with abs and zero RHS ==');
test('solve(abs(x) = 0, x)', () => cmd('solve', new Eq(new Call('abs', [x]), num(0)), x));

console.log('\n== solve: abs with positive RHS (valid) ==');
test('solve(abs(x) = 2, x)', () => cmd('solve', new Eq(new Call('abs', [x]), num(2)), x));

console.log('\n== diff: higher order of exp ==');
test('diff(exp(x), x, 10)', () => cmd('diff', new Call('exp', [x]), x, num(10)));

console.log('\n== taylor: higher order of exp ==');
test('taylor(exp(x), x, 0, 10)', () => cmd('taylor', new Call('exp', [x]), x, num(0), num(10)).toString());

console.log('\n== integrate: multiple integrals ==');
test('integrate(integrate(x*y, x), y)', () => cmd('integrate', cmd('integrate', new Mul(x, y), x), y).toString());

console.log('\n== series: sqrt(1+x) ==');
test('series(sqrt(1+x), x, 0, 5)', () => cmd('series', new Call('sqrt', [new Add(num(1), x)]), x, num(0), num(5)).toString());

console.log('\n== sum: nested sum ==');
test('sum(sum(k, k, 1, m), m, 1, n)', () => {
  const inner = cmd('sum', k, k, num(1), k);
  return cmd('sum', inner, k, num(1), n);
});

console.log('\n== limit: piecewise at discontinuity ==');
const pw = new Call('piecewise', [new Lt(x, num(0)), num(1), num(0)]);
test('limit(piecewise(x<0,1,0), x, 0)', () => cmd('limit', pw, x, num(0)));

console.log('\n== taylor: arcsin(x) ==');
test('taylor(asin(x), x, 0, 5)', () => cmd('taylor', new Call('asin', [x]), x, num(0), num(5)).toString());

console.log('\n== taylor: arctan(x) ==');
test('taylor(atan(x), x, 0, 5)', () => cmd('taylor', new Call('atan', [x]), x, num(0), num(5)).toString());

console.log('\n== solve: nonlinear system ==');
const sysEqs = new Vec([new Eq(new Add(x, y), num(3)), new Eq(new Sub(x, new Mul(num(2), y)), num(0))]);
const sysVars = new Vec([x, y]);
test('solve([x+y=3, x-2y=0], [x,y])', () => cmd('solve', sysEqs, sysVars).toString());

console.log('\n== done ==');
