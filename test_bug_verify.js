const { CAS } = require('./js/cas.js');
const { Num, Sym, Call, Vec, Add, Sub, Mul, Div, Pow, Eq, Lt, Gt, Le, Ge } = require('./js/expression.js');
const cas = new CAS();
globalThis.HELP_DATA = {};
const x = new Sym('x'), y = new Sym('y'), k = new Sym('k'), n = new Sym('n'), a = new Sym('a');
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

console.log('\n== solve: abs with negative RHS ==');
// abs(x) = -3 -> no solutions, but CAS returns {3, -3}
test('solve(abs(x) = -3)', () => cmd('solve', new Eq(new Call('abs', [x]), num(-3)), x));

console.log('\n== solve: quadratic inequalities ==');
// x^2 - 4 > 0: x < -2 or x > 2
const x2_gt_4 = new Gt(new Sub(new Pow(x, num(2)), num(4)), num(0));
test('solve(x^2 - 4 > 0)', () => cmd('solve', x2_gt_4, x));
// -x^2 + 4 > 0: -2 < x < 2
const neg_x2_gt_neg4 = new Gt(new Add(new Mul(num(-1), new Pow(x, num(2))), num(4)), num(0));
test('solve(-x^2 + 4 > 0)', () => cmd('solve', neg_x2_gt_neg4, x));
// -x^2 - 4 > 0: never
const neg_x2_gt_4 = new Gt(new Add(new Mul(num(-1), new Pow(x, num(2))), num(-4)), num(0));
test('solve(-x^2 - 4 > 0)', () => cmd('solve', neg_x2_gt_4, x));

console.log('\n== solve: even power negative RHS ==');
// x^2 + 4 = 0 -> no real solutions
test('solve(x^2 = -4)', () => cmd('solve', new Eq(new Add(new Pow(x, num(2)), num(4)), num(0)), x));
// x^4 = 16 -> x = 2 or x = -2
test('solve(x^4 = 16)', () => cmd('solve', new Eq(new Sub(new Pow(x, num(4)), num(16)), num(0)), x));

console.log('\n== limit: 0^0 form ==');
test('limit(0^x, x, 0)', () => cmd('limit', new Pow(num(0), x), x, num(0)));

console.log('\n== limit: 1^Infinity form ==');
test('limit((1+1/x)^x, x, Infinity)', () => cmd('limit', new Pow(new Add(num(1), new Div(num(1), x)), x), x, new Sym('Infinity')));
test('limit((1+2/x)^x, x, Infinity)', () => cmd('limit', new Pow(new Add(num(1), new Div(num(2), x)), x), x, new Sym('Infinity')));

console.log('\n== limit: sin(x)/x at 0 ==');
test('limit(sin(x)/x, x, 0)', () => cmd('limit', new Div(new Call('sin', [x]), x), x, num(0)));

console.log('\n== product: symbolic n (factorial-like) ==');
test('product(k, k, 1, n)', () => cmd('product', k, k, num(1), n));
test('product(k+1, k, 1, n)', () => cmd('product', new Add(k, num(1)), k, num(1), n));
test('product(2^k, k, 1, n)', () => cmd('product', new Pow(num(2), k), k, num(1), n));

console.log('\n== sum: a/k symbolic ==');
test('sum(a/k, k, 1, n)', () => cmd('sum', new Div(a, k), k, num(1), n));

console.log('\n== taylor: removable singularity sin(x)/x ==');
test('taylor(sin(x)/x, x, 0, 4)', () => cmd('taylor', new Div(new Call('sin', [x]), x), x, num(0), num(4)));

console.log('\n== integrate piecewise ==');
test('integrate(piecewise(x<0, x, x^2), x)', () => cmd('integrate', new Call('piecewise', [new Lt(x, num(0)), x, new Pow(x, num(2))]), x).toString());

console.log('\n== diff piecewise ==');
test('diff(piecewise(x<0, x, x^2), x)', () => cmd('diff', new Call('piecewise', [new Lt(x, num(0)), x, new Pow(x, num(2))]), x).toString());

console.log('\n== series essential singularity ==');
test('series(exp(1/x), x, 0, 3)', () => cmd('series', new Call('exp', [new Div(num(1), x)]), x, num(0), num(3)).toString());

console.log('\n== implicitDiff ==');
test('implicitDiff(y^2 = x^3, y, x)', () => cmd('implicitDiff', new Eq(new Pow(y, num(2)), new Pow(x, num(3))), y, x).toString());

console.log('\n== limit abs at discontinuity ==');
test('limit(abs(x)/x, x, 0, right)', () => cmd('limit', new Div(new Call('abs', [x]), x), x, num(0), num(1)));
test('limit(abs(x)/x, x, 0, left)', () => cmd('limit', new Div(new Call('abs', [x]), x), x, num(0), num(-1)));

console.log('\n== solve transcendental ==');
test('solve(tan(x)=1)', () => cmd('solve', new Eq(new Call('tan', [x]), num(1)), x));
test('solve(x*exp(x)=1)', () => cmd('solve', new Eq(new Mul(x, new Call('exp', [x])), num(1)), x));

console.log('\n== product: product(sin(k), k, 1, n) ==');
test('product(sin(k), k, 1, n)', () => cmd('product', new Call('sin', [k]), k, num(1), n));

console.log('\n== sum: sum(product(j,j,1,k), k, 1, n) ==');
const prodKJ = new Call('product', [k, k, num(1), k]);
test('sum(product(k,k,1,k), k, 1, n)', () => cmd('sum', prodKJ, k, num(1), n));

console.log('\n== limit x^x as x->0+ ==');
test('limit(x^x, x, 0, right)', () => cmd('limit', new Pow(x, x), x, num(0), num(1)));

console.log('\n== limit sin(sin(x))/x as x->0 ==');
test('limit(sin(sin(x))/x, x, 0)', () => cmd('limit', new Div(new Call('sin', [new Call('sin', [x])]), x), x, num(0)));

console.log('\n== solve cubic x^3 = 1 ==');
test('solve(x^3 = 1)', () => cmd('solve', new Eq(new Sub(new Pow(x, num(3)), num(1)), num(0)), x));

console.log('\n== integrate sinh, cosh ==');
test('integrate(sinh(x), x)', () => cmd('integrate', new Call('sinh', [x]), x));
test('integrate(cosh(x), x)', () => cmd('integrate', new Call('cosh', [x]), x));

console.log('\n== done ==');
