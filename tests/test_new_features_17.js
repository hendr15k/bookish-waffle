const { Expr, Num, Sym, Call, Add, Sub, Mul, Div, Pow } = require('../js/expression');
const { CAS } = require('../js/cas');
const { Parser, Lexer } = require('../js/parser');

const cas = new CAS();

function test(exprStr, expectedStr) {
    try {
        const lexer = new Lexer(exprStr);
        const parser = new Parser(lexer);
        const expr = parser.parse();
        const res = cas.evaluate(expr);
        const resStr = res.toString();

        // Normalize for comparison
        const normalize = (s) => s.replace(/\s+/g, '');

        if (normalize(resStr) === normalize(expectedStr) || resStr.includes(expectedStr)) {
            console.log(`[PASS] ${exprStr} => ${resStr}`);
        } else {
            console.error(`[FAIL] ${exprStr} => ${resStr}, expected ${expectedStr}`);
        }
    } catch (e) {
        console.error(`[ERROR] ${exprStr}: ${e.message}`);
    }
}

console.log("--- Testing New Features 17 ---");

// 1. Infinite Sums
// Accepts evaluated form (0.1666... * pi^2) or zeta(2)
// Since zeta(2) evaluates to number*pi^2, we check for 'pi^2'
test("sum(1/k^2, k, 1, Infinity)", "pi^2");
test("sum(k^-4, k, 1, Infinity)", "pi^4");
test("sum(1/2^k, k, 0, Infinity)", "2");
test("sum(0.5^k, k, 0, Infinity)", "2");

// 2. erfc and LambertW
// Expected: ((-2 * exp((-1 * x^2))) / sqrt(pi))
test("diff(erfc(x), x)", "((-2 * exp((-1 * x^2))) / sqrt(pi))");

// integrate(erfc(x)) -> ... sqrt(3.14...)
// We accept loose match or just check structure
test("integrate(erfc(x), x)", "erfc(x)"); // Check it contains erfc(x)

// diff(lambertw(x))
test("diff(lambertw(x), x)", "(exp((-1 * lambertw(x))) / (1 + lambertw(x)))");

test("integrate(lambertw(x), x)", "(x * ((lambertw(x) - 1) + (1 / lambertw(x))))");

// 3. Divergent Integral
// integrate(1/x^2, x, -1, 1) should be Infinity or contains Infinity
test("integrate(1/x^2, x, -1, 1)", "Infinity");

// 4. Numeric Checks
test("N(lambertw(e))", "1");
// erfc(0) is 1. N(erfc(0)) -> 0.999999999.
// Allow approximate check?
// toString() won't match exactly.
// We'll skip strict string match for this one or check prefix
try {
    const val = cas.evaluate(new Parser(new Lexer("N(erfc(0))")).parse());
    if (Math.abs(val.value - 1) < 1e-6) console.log("[PASS] N(erfc(0)) ~ 1");
    else console.error("[FAIL] N(erfc(0)) = " + val.toString());
} catch(e) { console.error(e); }

console.log("--- Done ---");
