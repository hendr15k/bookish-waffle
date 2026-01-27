
const { CAS } = require('../js/cas.js');
const { Expr, Num, Sym, Call, Add, Div, Pow, Mul, Vec } = require('../js/expression.js');

// Mock loading - in real env, CAS uses global Expr classes.
// We need to setup global scope for CAS to work in Node environment test
global.Expr = Expr;
global.Num = Num;
global.Sym = Sym;
global.Call = Call;
global.Add = Add;
global.Div = Div;
global.Pow = Pow;
global.Mul = Mul;
global.Vec = Vec;

const cas = new CAS();

function runTest(name, exprStr, expectedStr) {
    try {
        console.log(`Testing ${name}...`);
        // We simulate parsing by constructing AST manually or using a simplified parser if available.
        // Here we construct AST manually for precision.

        let result = cas.evaluate(exprStr);
        let resStr = result.toString();

        if (resStr === expectedStr) {
            console.log(`PASS: ${resStr}`);
        } else {
            console.error(`FAIL: Expected ${expectedStr}, got ${resStr}`);
            process.exit(1);
        }
    } catch (e) {
        console.error(`ERROR in ${name}:`, e);
        process.exit(1);
    }
}

// 1. Partitions
// p(5) = 7
runTest('partitions(5)', new Call('partitions', [new Num(5)]), '7');
// p(10) = 42
runTest('partitions(10)', new Call('partitions', [new Num(10)]), '42');

// 2. Is Coprime
// 14, 15 -> 1
runTest('is_coprime(14, 15)', new Call('is_coprime', [new Num(14), new Num(15)]), '1');
// 14, 21 -> 0 (gcd 7)
runTest('is_coprime(14, 21)', new Call('is_coprime', [new Num(14), new Num(21)]), '0');

// 3. Egyptian Fraction
// 5/7 -> 1/2 + 1/5 + 1/70
// Output is a Vec [1/2, 1/5, 1/70]
// toString of Vec is [1/2, 1/5, 1/70] usually, assuming Num/Div toString is clean.
// Div toString is (1 / 2) usually.
// Let's verify string format.
// (1 / 2)
// Vec: [(1 / 2), (1 / 5), (1 / 70)]
const vecExpected = '[(1 / 2), (1 / 5), (1 / 70)]';
runTest('egyptian_fraction(5/7)', new Call('egyptian_fraction', [new Div(new Num(5), new Num(7))]), vecExpected);

// 4. Totient Sum
// Phi(10) = 32
runTest('totient_sum(10)', new Call('totient_sum', [new Num(10)]), '32');

console.log("All tests passed!");
