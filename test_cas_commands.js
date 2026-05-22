// BUG HUNT: CAS Commands - cas.js part 2
// Tests: solve, integrate, diff, limit, series, sum, product, plot
const { CAS } = require('./js/cas.js');
const { Num, Sym, Call, Vec, Add, Sub, Mul, Div, Pow, Eq, Lt, Gt, Le, Ge } = require('./js/expression.js');

const cas = new CAS();
globalThis.HELP_DATA = {};

const x = new Sym('x');
const y = new Sym('y');
const n = new Sym('n');
const k = new Sym('k');
const a = new Sym('a');
const b = new Sym('b');

function cmd(name, ...args) {
  return cas.evaluate(new Call(name, args));
}

function test(name, fn, expected) {
  try {
    const result = fn();
    const resultStr = result && result.toString ? result.toString() : String(result);
    let pass = expected === undefined ? true : (typeof expected === 'function' ? expected(result) : resultStr === expected.toString());
    console.log(`${pass ? 'PASS' : 'FAIL'}: ${name}`);
    if (!pass) {
      if (expected !== undefined) {
        console.log(`  Expected: ${expected}`);
      }
      console.log(`  Got: ${resultStr}`);
      console.log(`  Type: ${result && result.constructor ? result.constructor.name : typeof result}`);
    }
    return pass;
  } catch (e) {
    console.log(`ERROR: ${name}`);
    console.log(`  Exception: ${e.message}`);
    return false;
  }
}

function num(v) { return new Num(v); }
function sym(s) { return new Sym(s); }
function eq(l, r) { return new Eq(l, r); }
function lt(l, r) { return new Lt(l, r); }
function gt(l, r) { return new Gt(l, r); }

console.log('\n========== SOLVE bugs ==========');
{
  // Basic linear
  test('solve(2*x + 2*y = 0, x)', () => cmd('solve', eq(new Add(new Mul(num(2), x), new Mul(num(2), y)), num(0)), x));

  // Transcendental returns unevaluated (expected behavior)
  test('solve(sin(x) - 1, x)', () => cmd('solve', eq(new Sub(new Call('sin', [x]), num(1)), num(0)), x));

  // Constant equation
  test('solve(5 = 0, x)', () => cmd('solve', eq(num(5), num(0)), x).toString());

  // Solve with 0 = 0 (any x)
  test('solve(0*x = 0, x)', () => cmd('solve', eq(new Mul(num(0), x), num(0)), x).toString());

  // Solve with Eq object
  test('solve(x = 2, x)', () => cmd('solve', eq(x, num(2)), x).toString());

  // Cubic with multiple roots
  test('solve(x^3 - x, x)', () => cmd('solve', eq(new Sub(new Pow(x, num(3)), x), num(0)), x).toString());

  // Constant-only equation: 5 = 0 should give empty set or error
  test('solve(5, x)', () => cmd('solve', num(5), x).toString());

  // Quadratic with discriminant 0
  test('solve(x^2 - 2*x + 1 = 0, x)', () => cmd('solve', eq(new Sub(new Pow(x, num(2)), new Mul(num(2), x)), num(0)), x).toString());

  // Exponential equation
  test('solve(exp(x) - 1, x)', () => cmd('solve', eq(new Sub(new Call('exp', [x]), num(1)), num(0)), x).toString());

  // Inequality
  test('solve(x^2 < 4, x)', () => cmd('solve', lt(new Pow(x, num(2)), num(4)), x).toString());

  // Abs equation
  test('solve(abs(x) = 3, x)', () => cmd('solve', eq(new Call('abs', [x]), num(3)), x).toString());

  // Div equation P/Q = 0 -> P = 0
  test('solve((x-1)/(x-2) = 0, x)', () => cmd('solve', eq(new Div(new Sub(x, num(1)), new Sub(x, num(2))), num(0)), x).toString());

  // High degree via factoring
  test('solve(x^4 - 1, x)', () => cmd('solve', new Sub(new Pow(x, num(4)), num(1)), x).toString());

  // Mul: both sides depend on var
  test('solve(x * y, x)', () => cmd('solve', new Mul(x, y), x).toString());

  // Mul: only right depends
  test('solve(y * x, x)', () => cmd('solve', new Mul(y, x), x).toString());
}

console.log('\n========== INTEGRATE bugs ==========');
{
  // Already in BUGS_CALCULUS - we test new edge cases

  // sinh(x)
  test('integrate(sinh(x), x)', () => cmd('integrate', new Call('sinh', [x]), x).toString());

  // cosh(x)
  test('integrate(cosh(x), x)', () => cmd('integrate', new Call('cosh', [x]), x).toString());

  // sec(x)
  test('integrate(sec(x), x)', () => cmd('integrate', new Call('sec', [x]), x).toString());

  // csc(x)
  test('integrate(csc(x), x)', () => cmd('integrate', new Call('csc', [x]), x).toString());

  // x*ln(x) - integration by parts
  test('integrate(x*ln(x), x)', () => cmd('integrate', new Mul(x, new Call('ln', [x])), x).toString());

  // 1/ln(x) - non-elementary
  test('integrate(1/ln(x), x)', () => cmd('integrate', new Div(num(1), new Call('ln', [x])), x).toString());

  // tanh(x)
  test('integrate(tanh(x), x)', () => cmd('integrate', new Call('tanh', [x]), x).toString());

  // sech(x)
  test('integrate(sech(x), x)', () => cmd('integrate', new Call('sech', [x]), x).toString());

  // x*exp(x)
  test('integrate(x*exp(x), x)', () => cmd('integrate', new Mul(x, new Call('exp', [x])), x).toString());

  // x^2*exp(x)
  test('integrate(x^2*exp(x), x)', () => cmd('integrate', new Mul(new Pow(x, num(2)), new Call('exp', [x])), x).toString());

  // Definite integral
  test('integrate(x^2, x, 0, 1)', () => cmd('integrate', new Pow(x, num(2)), x, num(0), num(1)).toString());

  // Definite integral with infinite upper
  test('integrate(x*exp(-x), x, 0, Infinity)', () => cmd('integrate', new Mul(x, new Call('exp', new [new Mul(num(-1), x)])), x, num(0), new Sym('Infinity')).toString());

  // sin(x)*cos(x)
  test('integrate(sin(x)*cos(x), x)', () => cmd('integrate', new Mul(new Call('sin', [x]), new Call('cos', [x])), x).toString());

  // sec^2(x)
  test('integrate(sec(x)^2, x)', () => cmd('integrate', new Pow(new Call('sec', [x]), num(2)), x).toString());
}

console.log('\n========== DIFF bugs ==========');
{
  // sinh
  test('diff(sinh(x), x)', () => cmd('diff', new Call('sinh', [x]), x).toString());

  // cosh
  test('diff(cosh(x), x)', () => cmd('diff', new Call('cosh', [x]), x).toString());

  // abs
  test('diff(abs(x), x)', () => cmd('diff', new Call('abs', [x]), x).toString());

  // higher order
  test('diff(sin(x), x, 2)', () => cmd('diff', new Call('sin', [x]), x, num(2)).toString());

  // ln(ln(x))
  test('diff(ln(ln(x)), x)', () => cmd('diff', new Call('ln', [new Call('ln', [x])]), x).toString());

  // x^x
  test('diff(x^x, x)', () => cmd('diff', new Pow(x, x), x).toString());

  // log
  test('diff(log(x), x)', () => cmd('diff', new Call('log', [x]), x).toString());

  // asinh
  test('diff(asinh(x), x)', () => cmd('diff', new Call('asinh', [x]), x).toString());

  // acosh
  test('diff(acosh(x), x)', () => cmd('diff', new Call('acosh', [x]), x).toString());

  // atanh
  test('diff(atanh(x), x)', () => cmd('diff', new Call('atanh', [x]), x).toString());

  // multiple vars
  test('diff(x*y, x)', () => cmd('diff', new Mul(x, y), x).toString());

  // product of functions
  test('diff(sin(x)*cos(x), x)', () => cmd('diff', new Mul(new Call('sin', [x]), new Call('cos', [x])), x).toString());

  // piecewise
  test('diff(piecewise(x<0, -1, 1), x)', () => cmd('diff', new Call('piecewise', [lt(x, num(0)), num(-1), num(1)]), x).toString());
}

console.log('\n========== LIMIT bugs ==========');
{
  // sin(x)/x as x -> 0
  test('limit(sin(x)/x, x, 0)', () => cmd('limit', new Div(new Call('sin', [x]), x), x, num(0)).toString());

  // (1-cos(x))/x^2 as x -> 0
  test('limit((1-cos(x))/x^2, x, 0)', () => cmd('limit', new Div(new Sub(num(1), new Call('cos', [x])), new Pow(x, num(2))), x, num(0)).toString());

  // (exp(x)-1)/x as x -> 0
  test('limit((exp(x)-1)/x, x, 0)', () => cmd('limit', new Div(new Sub(new Call('exp', [x]), num(1)), x), x, num(0)).toString());

  // (1+1/x)^x as x -> Infinity
  test('limit((1+1/x)^x, x, Infinity)', () => cmd('limit', new Pow(new Add(num(1), new Div(num(1), x)), x), x, new Sym('Infinity')).toString());

  // ln(x)/x as x -> Infinity
  test('limit(ln(x)/x, x, Infinity)', () => cmd('limit', new Div(new Call('ln', [x]), x), x, new Sym('Infinity')).toString());

  // x*sin(1/x) as x -> 0
  test('limit(x*sin(1/x), x, 0)', () => cmd('limit', new Mul(x, new Call('sin', new [new Div(num(1), x)])), x, num(0)).toString());

  // directional limit right
  test('limit(1/x, x, 0, right)', () => cmd('limit', new Div(num(1), x), x, num(0), num(1)).toString());

  // directional limit left
  test('limit(1/x, x, 0, left)', () => cmd('limit', new Div(num(1), x), x, num(0), num(-1)).toString());

  // symbolic point
  test('limit(x^2, x, a)', () => cmd('limit', new Pow(x, num(2)), x, a).toString());

  // tan at pi/2 right
  test('limit(tan(x), x, new Mul(Div(num(1), num(2)), pi), right)', () => cmd('limit', new Call('tan', [x]), x, new Div(new Call('pi', []), num(2)), num(1)).toString());

  // sin(x)/x at Infinity
  test('limit(sin(x)/x, x, Infinity)', () => cmd('limit', new Div(new Call('sin', [x]), x), x, new Sym('Infinity')).toString());

  // sqrt with substitution
  test('limit((sqrt(x)-sqrt(a))/(x-a), x, a)', () => cmd('limit', new Div(new Sub(new Call('sqrt', [x]), new Call('sqrt', [a])), new Sub(x, a)), x, a).toString());

  // piecewise limit
  test('limit(piecewise(x<0, 0, x), x, 0)', () => cmd('limit', new Call('piecewise', [lt(x, num(0)), num(0), x]), x, num(0)).toString());

  // x*ln(x) -> 0
  test('limit(x*ln(x), x, 0, right)', () => cmd('limit', new Mul(x, new Call('ln', [x])), x, num(0), num(1)).toString());
}

console.log('\n========== SERIES/TAYLOR bugs ==========');
{
  // taylor sin
  test('taylor(sin(x), x, 0, 5)', () => cmd('taylor', new Call('sin', [x]), x, num(0), num(5)).toString());

  // taylor cos
  test('taylor(cos(x), x, 0, 4)', () => cmd('taylor', new Call('cos', [x]), x, num(0), num(4)).toString());

  // taylor exp
  test('taylor(exp(x), x, 0, 3)', () => cmd('taylor', new Call('exp', [x]), x, num(0), num(3)).toString());

  // taylor ln(1+x)
  test('taylor(ln(1+x), x, 0, 5)', () => cmd('taylor', new Call('ln', [new Add(num(1), x)]), x, num(0), num(5)).toString());

  // taylor 1/(1-x)
  test('taylor(1/(1-x), x, 0, 5)', () => cmd('taylor', new Div(num(1), new Sub(num(1), x)), x, num(0), num(5)).toString());

  // taylor at non-zero point
  test('taylor(sin(x), x, 1, 3)', () => cmd('taylor', new Call('sin', [x]), x, num(1), num(3)).toString());

  // taylor with singularity at point
  test('taylor(1/x, x, 1, 3)', () => cmd('taylor', new Div(num(1), x), x, num(1), num(3)).toString());

  // series of 1/x (Laurent)
  test('series(1/x, x, 0, 5)', () => cmd('series', new Div(num(1), x), x, num(0), num(5)).toString());

  // taylor ln(x) at 1
  test('taylor(ln(x), x, 1, 4)', () => cmd('taylor', new Call('ln', [x]), x, num(1), num(4)).toString());
}

console.log('\n========== SUM bugs ==========');
{
  // k^2 symbolic
  test('sum(k^2, k, 1, n)', () => cmd('sum', new Pow(k, num(2)), k, num(1), n).toString());

  // k^3 symbolic
  test('sum(k^3, k, 1, n)', () => cmd('sum', new Pow(k, num(3)), k, num(1), n).toString());

  // 1/k symbolic (harmonic)
  test('sum(1/k, k, 1, n)', () => cmd('sum', new Div(num(1), k), k, num(1), n).toString());

  // 2^k geometric
  test('sum(2^k, k, 1, n)', () => cmd('sum', new Pow(num(2), k), k, num(1), n).toString());

  // telescoping 1/(k*(k+1))
  test('sum(1/(k*(k+1)), k, 1, n)', () => cmd('sum', new Div(num(1), new Mul(k, new Add(k, num(1)))), k, num(1), n).toString());

  // k numerical
  test('sum(k^2, k, 1, 5)', () => cmd('sum', new Pow(k, num(2)), k, num(1), num(5)).toString());

  // sum over list
  test('sum([1,2,3,4,5])', () => cmd('sum', new Vec([num(1), num(2), num(3), num(4), num(5)])).toString());

  // symbolic lower bound
  test('sum(k^2, k, m, n)', () => cmd('sum', new Pow(k, num(2)), k, sym('m'), n).toString());

  // constant sum
  test('sum(5, k, 1, 10)', () => cmd('sum', num(5), k, num(1), num(10)).toString());

  // x^k
  test('sum(x^k, k, 1, n)', () => cmd('sum', new Pow(x, k), k, num(1), n).toString());
}

console.log('\n========== PRODUCT bugs ==========');
{
  // k factorial-like
  test('product(k, k, 1, n)', () => cmd('product', k, k, num(1), n).toString());

  // numerical
  test('product(k, k, 1, 5)', () => cmd('product', k, k, num(1), num(5)).toString());

  // 2^k
  test('product(2, k, 1, n)', () => cmd('product', num(2), k, num(1), n).toString());

  // over list
  test('product([1,2,3,4,5])', () => cmd('product', new Vec([num(1), num(2), num(3), num(4), num(5)])).toString());

  // k+1
  test('product(k+1, k, 0, n)', () => cmd('product', new Add(k, num(1)), k, num(0), n).toString());

  // 1 - k^2
  test('product(1 - k^2, k, 1, n)', () => cmd('product', new Sub(num(1), new Pow(k, num(2))), k, num(1), n).toString());

  // symbolic lower bound
  test('product(k, k, m, n)', () => cmd('product', k, k, sym('m'), n).toString());
}

console.log('\n========== PLOT bugs ==========');
{
  // Basic plot
  test('plot(sin(x), x)', () => { const r = cmd('plot', new Call('sin', [x]), x); return r.type === 'plot' ? 'plot_object' : r.toString(); }, 'plot_object');

  // Plot with range
  test('plot(sin(x), x, 0, 10)', () => { const r = cmd('plot', new Call('sin', [x]), x, num(0), num(10)); return r.type === 'plot' ? `min=${r.min}_max=${r.max}` : r.toString(); });

  // Multiple expressions
  test('plot([sin(x), cos(x)], x)', () => { const r = cmd('plot', new Vec([new Call('sin', [x]), new Call('cos', [x])]), x); return r.type === 'plot' ? `n=${r.expressions.length}` : r.toString(); });

  // nstep
  test('plot(sin(x), x, 0, 10, nstep=200)', () => { const r = cmd('plot', new Call('sin', [x]), x, num(0), num(10), new Eq(sym('nstep'), num(200))); return r.type === 'plot' ? `nstep=${r.nstep}` : r.toString(); });

  // color
  test('plot(sin(x), x, color=red)', () => { const r = cmd('plot', new Call('sin', [x]), x, new Eq(sym('color'), sym('red'))); return r.type === 'plot' ? `color=${r.color}` : r.toString(); });
}

console.log('\n========== ADDITIONAL COMMANDS ==========');
{
  // implicitDiff
  test('implicitDiff(x^2 + y^2 = 1, y, x)', () => cmd('implicitDiff', eq(new Add(new Pow(x, num(2)), new Pow(y, num(2))), num(1)), y, x).toString());

  // extrema
  test('extrema(x^3 - 3*x, x)', () => cmd('extrema', new Sub(new Pow(x, num(3)), new Mul(num(3), x)), x).toString());

  // stationary_points
  test('stationary_points(sin(x), x)', () => cmd('stationary_points', new Call('sin', [x]), x).toString());

  // asymptotes
  test('asymptotes(1/x, x)', () => cmd('asymptotes', new Div(num(1), x), x).toString());

  // completeSquare
  test('completeSquare(x^2 + 4*x + 5, x)', () => cmd('completeSquare', new Add(new Pow(x, num(2)), new Add(new Mul(num(4), x), num(5))), x).toString());

  // discriminant
  test('discriminant(x^2 + 3*x + 2, x)', () => cmd('discriminant', new Add(new Pow(x, num(2)), new Add(new Mul(num(3), x), num(2))), x).toString());

  // resultant
  test('resultant(x^2-1, x^2-4, x)', () => cmd('resultant', new Sub(new Pow(x, num(2)), num(1)), new Sub(new Pow(x, num(2)), num(4)), x).toString());

  // cumsum
  test('cumsum([1, 2, 3, 4])', () => cmd('cumsum', new Vec([num(1), num(2), num(3), num(4)])).toString());

  // flatten
  test('flatten([[1,2],[3,4]])', () => cmd('flatten', new Vec([new Vec([num(1), num(2)]), new Vec([num(3), num(4)])])).toString());

  // minimize
  test('minimize(x^2, x)', () => cmd('minimize', new Pow(x, num(2)), x).toString());

  // maximize
  test('maximize(-x^2, x)', () => cmd('maximize', new Mul(num(-1), new Pow(x, num(2))), x).toString());

  // nIntegrate (numeric)
  test('nIntegrate(x^2, x, 0, 1)', () => cmd('nIntegrate', new Pow(x, num(2)), x, num(0), num(1)).toString());
}

console.log('\n========== DONE ==========');
