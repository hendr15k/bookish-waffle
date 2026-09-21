// Regression tests for the September batch of CAS fixes:
//  - polynomial gcd/lcm (documented help example was unevaluated)
//  - e^(linear) integration (unlocking general u-substitution)
//  - limits at infinity: `oo` alias, L'Hôpital, (1+1/n)^n -> e
//  - variance/std and (a/b)^n numeric evaluation
//  - like-term collection across nested Sub nodes, display of 0 - a
const { ev, test } = require('./harness.js');

let failed = 0;
const check = (desc, input, expected) => {
    if (!test(desc, input, expected)) failed++;
};

// --- Algebra: polynomial gcd / lcm -----------------------------------------
check('gcd help example', 'gcd(x^2-1, x^2+2*x+1)', '(x + 1)');
check('gcd x^3-1, x^2-1', 'gcd(x^3-1, x^2-1)', '(x - 1)');
check('gcd x^2-1, x-1', 'gcd(x^2-1, x-1)', '(x - 1)');
check('gcd integers', 'gcd(12, 18)', '6');
check('lcm polynomial', 'lcm(x^2-1, x-1)', '(x^2 - 1)');
check('lcm integers', 'lcm(4, 6)', '12');

// --- Calculus: linear-exponent integration ----------------------------------
check('int e^(2x)', 'integrate(e^(2*x), x)', '(e^((2 * x)) / 2)');
check('int exp(3x)', 'integrate(exp(3*x), x)', '(exp((3 * x)) / 3)');
check('int e^(2x+1)', 'integrate(e^(2*x+1), x)', '(e^(((2 * x) + 1)) / 2)');
check('int 2^(3x)', 'integrate(2^(3*x), x)', '(2^((3 * x)) / ln(8))');

// --- Limits at infinity (incl. the `oo` alias) ------------------------------
check('lim sin(x)/x at oo', 'limit(sin(x)/x, x, oo)', '0');
check('lim ln(x)/x at oo', 'limit(ln(x)/x, x, oo)', '0');
check('lim 1/x at oo', 'limit(1/x, x, oo)', '0');
check('lim 1/x at -oo', 'limit(1/x, x, -oo)', '0');
check('lim sin(x)/x at inf', 'limit(sin(x)/x, x, inf)', '0');
check('lim (1+1/n)^n at oo', 'limit((1+1/n)^n, n, oo)', 'e');
check('lim e^x/x at oo', 'limit(e^x/x, x, oo)', 'Infinity');
check('lim sin(x)/x at 0 (unchanged)', 'limit(sin(x)/x, x, 0)', '1');
check('lim 1/x from the right (unchanged)', 'limit(1/x, x, 0, 1)', 'Infinity');

// --- Statistics / numeric simplification ------------------------------------
check('variance numeric', 'variance([1,2,3,4])', '(5 / 3)');
check('std numeric', 'std([1,2,3,4])', 'sqrt((5 / 3))');
check('power of a fraction', '(1/2)^2', '(1 / 4)');
check('power of a negative fraction', '(-3/2)^2', '(9 / 4)');
check('negative fraction, negative power', '(3/2)^(-2)', '(4 / 9)');

// --- Simplification / display -----------------------------------------------
check('collect across nested Sub', '(x^4 - x) - (x^4 - x^2)', '(x^2 - x)');
check('collect x - 2x', 'x - 2*x + 5', '(-x + 5)');
check('collect symbolic opposite', '(a+b)-(a-b)', '(2 * b)');
check('like terms still work', '2*x + 3*x', '(5 * x)');
check('trigReduce sin^2', 'trigReduce(sin(x)^2)', '((-cos((2 * x)) + 1) / 2)');

if (failed > 0) {
    console.error(`\n${failed} regression test(s) failed`);
    process.exit(1);
}
console.log('\nAll regression tests passed');
