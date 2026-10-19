// Regression tests for the October 2026 bug batch.
//
// Covered root causes:
//   1. Lexer rejected single-quoted strings (documented `'call'` / `'inch'` forms).
//   2. Lexer rejected a single '&' although cnf/dnf examples use it.
//   3. `mod_inverse` / `mod_pow` were lexed as the infix `mod` keyword.
//   4. Call.simplify() let `simplify()` pull a numeric coefficient back out of a
//      sum, which left expand() results only partially distributed.
//   5. Pow/Call.integrate() had no chain rule for linear arguments, so
//      integrate(x*exp(-x), x) recursed into a bogus series and the definite
//      form returned NaN; sqrt(x) stayed unevaluated; constant integrands
//      (sqrt(5), a Call with no variable) were not integrated at all.
//   6. Documented command forms threw: substitute(expr, var, val), quo/rem with
//      an explicit variable, div() as integer division, permutations(list),
//      del(f, vars), angle(z), toCylindrical/toSpherical(x, y, z),
//      lineEquation(p1, p2, vars), simplex/convolution/list forms.
//   7. ctrb/obsv/ss2tf/acker/isControllable/isObservable were dispatched but the
//      methods did not exist ("this._ctrb is not a function").
//   8. cfrac ignored its depth argument and crashed on a symbolic constant.
//   9. Unit conversion lacked the spelled-out unit names used by the docs.
const { ev, evalExpr } = require('./harness.js');

let failed = 0;
const num = (text) => {
    try { return evalExpr(text).evaluateNumeric(); } catch (e) { return NaN; }
};
const near = (a, b, tol = 1e-9) => Number.isFinite(a) && Number.isFinite(b) && Math.abs(a - b) <= tol * Math.max(1, Math.abs(a), Math.abs(b));

const checkStr = (desc, input, expected) => {
    let got;
    try { got = ev(input); } catch (e) { got = 'ERROR: ' + e.message; }
    if (got === expected) { console.log(`[PASS] ${desc}`); return; }
    failed++;
    console.error(`[FAIL] ${desc}\n  input:    ${input}\n  expected: ${expected}\n  got:      ${got}`);
};

const checkNum = (desc, input, expected, tol = 1e-9) => {
    const got = num(input);
    if (near(got, expected, tol)) { console.log(`[PASS] ${desc}`); return; }
    failed++;
    console.error(`[FAIL] ${desc}\n  input:    ${input}\n  expected: ${expected}\n  got:      ${got}`);
};

const checkNoThrow = (desc, input) => {
    try {
        const r = ev(input);
        if (/^ERROR:/.test(r)) throw new Error(r);
        console.log(`[PASS] ${desc}`);
    } catch (e) {
        failed++;
        console.error(`[FAIL] ${desc}\n  input:    ${input}\n  error:    ${e.message}`);
    }
};

// --- 1. single-quoted strings ------------------------------------------------
checkStr('single-quoted string literal', "'abc'", '"abc"');
checkNoThrow('blackScholes accepts the documented string argument', "blackScholes(100, 100, 1, 0.05, 0.2, 'call')");
checkNum('blackScholes call value', "blackScholes(100, 100, 1, 0.05, 0.2, 'call')", 10.450583572185565, 1e-6);
checkNum('unit conversion with spelled-out unit', "convert(1, 'inch', 'cm')", 2.54);
checkNum('unit conversion km -> m', "convert(1, 'km', 'm')", 1000);

// --- 2. single '&' as logical AND -------------------------------------------
checkStr('1 & 0', '1 & 0', '0');
checkStr('1 & 1', '1 & 1', '1');
checkStr('&& still works', '1 && 1', '1');

// --- 3. mod_inverse / mod_pow are commands, not the mod operator -------------
checkStr('mod_inverse', 'mod_inverse(3, 11)', '4');
checkStr('mod_pow', 'mod_pow(2, 10, 1000)', '24');
checkStr('mod still parses as a function', 'mod(10, 3)', '1');
checkStr('infix mod still works', '10 mod 3', '1');

// --- 4. expand() must be fully distributed ----------------------------------
checkStr('expand((x+1)^3) is fully expanded', 'expand((x+1)^3)', '(((x^3 + (3 * x^2)) + (3 * x)) + 1)');
checkStr('plain polynomial is not factored', '3*x^2+3*x', '((3 * x^2) + (3 * x))');
checkStr('expand((x+1)^5)', 'expand((x+1)^5)', '(((((x^5 + (5 * x^4)) + (10 * x^3)) + (10 * x^2)) + (5 * x)) + 1)');
checkStr('expand((x^2+x+1)^2)', 'expand((x^2+x+1)^2)', '((((x^4 + (2 * x^3)) + (3 * x^2)) + (2 * x)) + 1)');
checkStr('expand((a+b)*(a-b))', 'expand((a+b)*(a-b))', '(a^2 - b^2)');
checkStr('symbolic common factor still factorises', 'a*x+a*y', '(a * (x + y))');
checkStr('expanded (x+y+z)^2 matches the closed form', 'subs(expand((x+y+z)^2), [x=1, y=2, z=3])', '36');

// --- 5. integration ---------------------------------------------------------
checkNum('integrate(x*exp(-x), x) at a point', 'substitute(integrate(x*exp(-x), x), x, 1)', -2 / Math.E, 1e-9);
checkNum('definite integrate(x*exp(-x), x, 0, 2)', 'integrate(x*exp(-x), x, 0, 2)', 1 - 3 * Math.exp(-2), 1e-9);
checkNum('definite integrate(x*exp(-x), x, 1, 2)', 'integrate(x*exp(-x), x, 1, 2)', (2 * Math.exp(-1) - 3 * Math.exp(-2)), 1e-9);
checkNum('integrate(sqrt(x), x, 0, 4)', 'integrate(sqrt(x), x, 0, 4)', 16 / 3, 1e-9);
checkNum('integrate(sin(2*x), x)', 'integrate(sin(2*x), x, 0, 1)', (1 - Math.cos(2)) / 2, 1e-9);
checkNum('integrate(cos(3*x), x)', 'integrate(cos(3*x), x, 0, 1)', Math.sin(3) / 3, 1e-9);
checkNum('integrate((2*x+1)^3, x, 0, 1)', 'integrate((2*x+1)^3, x, 0, 1)', (81 - 1) / 8, 1e-9);
checkNum('constant Call integrand', 'integrate(sqrt(5), x, 0, 1)', Math.sqrt(5), 1e-9);
checkNum('integrate(1/sqrt(x), x, 0, 4)', 'integrate(1/sqrt(x), x, 0, 4)', 4, 1e-9);
checkNum('clamped-by-chain rule stays exact', 'integrate(exp(2*x+1), x, 0, 1)', (Math.exp(3) - Math.exp(1)) / 2, 1e-9);

// --- 6. documented command forms -------------------------------------------
checkStr('substitute(expr, var, value)', 'substitute(x+y, x, 3)', '(3 + y)');
checkStr('subs two-arg form still works', 'subs(x^2, x, 2)', '4');
checkStr('quo with explicit variable', 'quo(x^2-1, x-1, x)', '(x + 1)');
checkStr('rem with explicit variable', 'rem(x^2, x-1, x)', '1');
checkStr('div is integer division', 'div(10, 3)', '3');
checkStr('divergence is unchanged', 'divergence([x,y,z],[x,y,z])', '3');
checkStr('div still diverges vector fields', 'div([x,y,z],[x,y,z])', '3');
checkStr('permutations(list)', 'permutations([1,2])', '[[1, 2], [2, 1]]');
checkStr('del(f, vars) is the nabla operator', 'del(x^2*y, [x,y])', '[((2 * x) * y), x^2]');
checkNum('angle of a complex number', 'angle(1+i)', Math.PI / 4, 1e-9);
checkStr('toCylindrical with three scalars', 'toCylindrical(1, 1, 0)', '[sqrt(2), 0.785398163397, 0]');
checkStr('lineEquation honours the variable list', 'lineEquation([0,0], [1,1], [x,y])', 'y = x');
checkNoThrow('lagrangeMultipliers with a single constraint', 'lagrangeMultipliers(x^2+y^2, x+y-1, [x,y])');
checkStr('discrete convolution', 'convolution([1,2], [3,4])', '[3, 10, 8]');
checkStr('simplex with three arguments', 'simplex([1,1], [[1,2],[2,1]], [4,5])', '[[2, 1], 3]');
checkStr('map(list, expr)', 'map([1,2,3], x^2)', '[1, 4, 9]');

// --- 7. state-space commands ------------------------------------------------
checkStr('ctrb', 'ctrb([[0,1],[-2,-3]], [0,1])', '[[0, 1], [1, -3]]');
checkStr('obsv', 'obsv([[0,1],[-2,-3]], [[1,0]])', '[[1, 0], [0, 1]]');
checkStr('isControllable', 'isControllable([[0,1],[-2,-3]], [0,1])', '1');
checkStr('isObservable', 'isObservable([[0,1],[-2,-3]], [[1,0]])', '1');
checkStr('isControllable false case', 'isControllable([[1,0],[0,1]], [1,0])', '0');
checkStr('acker pole placement', 'acker([[0,1],[-2,-3]], [0,1], [-3,-4])', '[10, 4]');
// Closed loop A - B*K must have the requested characteristic polynomial s^2+7s+12
checkStr('acker closed loop', '[[0,1],[-2,-3]] - [[0],[1]]*[[10,4]]', '[[0, 1], [-12, -7]]');
checkNum('ss2tf at s=0', 'substitute(ss2tf([[0,1],[-2,-3]], [0,1], [[1,0]], 0), s, 0)', 0.5, 1e-9);
checkNum('ss2tf at s=1', 'substitute(ss2tf([[0,1],[-2,-3]], [0,1], [[1,0]], 0), s, 1)', 1 / 6, 1e-9);
checkStr('eigenvalues alias', 'eigenvalues([[1,2],[2,1]])', '[-1, 3]');

// --- 8. continued fractions -------------------------------------------------
checkStr('cfrac honours the depth argument', 'cfrac(pi, 5)', '[3, 7, 15, 1, 292]');
checkStr('continued_fraction on e', 'continued_fraction(e, 5)', '[2, 1, 2, 1, 1]');

// --- 9. polynomial cancellation must keep the sign ---------------------------
// polyGcd() normalised its result in place, which corrupted the caller's
// denominator coefficients and flipped the sign of the quotient.
checkStr('(x^2-4)/(2-x)', '(x^2-4)/(2-x)', '-(2 + x)');
checkStr('(x^2-9)/(3-x)', '(x^2-9)/(3-x)', '-(3 + x)');
checkNum('(x^2-4)/(2-x) at x=3', 'substitute((x^2-4)/(2-x), x, 3)', -5);
checkNum('(4-x^2)/(2-x) at x=3', 'substitute((4-x^2)/(2-x), x, 3)', 5);
checkStr('(x^2-4)/(x-2) unchanged', '(x^2-4)/(x-2)', '(2 + x)');
checkStr('gcd still cancels', '(2*x+4)/(x+2)', '2');

// --- 10. powers of powers ----------------------------------------------------
// (x^2)^(1/2) is |x|, so it must not collapse to x, and the printed form must
// re-parse as the same expression ((x^2)^(1/2), not x^2^(1/2)).
checkStr('fractional root of an even power stays symbolic', 'simplify((x^2)^(1/2))', '(x^2)^((1 / 2))');
checkNum('(x^2)^(1/2) at x=-3 is |x|', 'substitute(simplify((x^2)^(1/2)), x, -3)', 3);
checkStr('integer outer exponent still collapses', 'simplify((x^2)^3)', 'x^6');
checkStr('symbolic exponent still collapses', '(a^x)^y', 'a^((x * y))');
checkStr('sqrt(x^2) is still abs(x)', 'sqrt(x^2)', 'abs(x)');

// --- 11. numeric edge cases --------------------------------------------------
// Splitting (a/b)^n into a^n/b^n overflows for large n (360 months is normal):
// 121^360 and 120^360 are both Infinity, so the ratio came out NaN.
checkNum('(12.1/12)^(-360)', '(12.1/12)^(-360)', Math.pow(12.1 / 12, -360), 1e-9);
checkStr('small fractional powers stay exact', '(1/2)^2', '(1 / 4)');
// loan(P, annual_rate, years) treated the years as months
checkNum('loan() 30 years', 'loan(200000, 0.04, 30)', 200000 * (0.04 / 12) / (1 - Math.pow(1 + 0.04 / 12, -360)), 1e-9);
// `oo` was only normalised for limits, not for sums/products
checkStr('sum with oo bound', 'sum(1/n^2, n, 1, oo)', '(pi^2 / 6)');
checkStr('sum with inf bound unchanged', 'sum(1/n^2, n, 1, inf)', '(pi^2 / 6)');

// --- 9. every documented help example evaluates -----------------------------
{
    // help.js registers window.HELP_DATA; chemistry.js must live inside the
    // harness sandbox because the CAS looks up globalThis.Chemistry there.
    const fs = require('fs');
    const vm = require('vm');
    const path = require('path');
    const { sandbox } = require('./harness.js');
    vm.runInContext(fs.readFileSync(path.join(__dirname, '../js/chemistry.js'), 'utf8'), sandbox);
    global.window = global.window || {};
    require('../js/help.js');
    const HELP = global.window.HELP_DATA;
    let bad = [];
    for (const [name, info] of Object.entries(HELP)) {
        for (const ex of (info.examples || [])) {
            try {
                const out = evalExpr(ex).toString();
                if (/^(integrate|diff|limit|sum|product|solve)\(/.test(out)) bad.push(`${name}: ${ex} -> ${out}`);
            } catch (e) {
                bad.push(`${name}: ${ex} -> ${e.message}`);
            }
        }
    }
    if (bad.length === 0) {
        console.log('[PASS] all documented help examples evaluate');
    } else {
        failed++;
        console.error(`[FAIL] ${bad.length} documented help examples fail:\n  ` + bad.slice(0, 12).join('\n  '));
    }
}

if (failed > 0) {
    console.error(`\n${failed} regression test(s) failed`);
    process.exit(1);
}
console.log('\nAll bug-batch regression tests passed');
