const { CAS } = require('../js/cas.js');
const { Sym, Num, Vec, Call, Div, Add, Sub, Mul, Pow } = require('../js/expression.js');

const test_new_features_5 = () => {
    const cas = new CAS();

    const assert = (condition, message) => {
        if (!condition) {
            console.error(`❌ ${message}`);
            throw new Error(message);
        } else {
            console.log(`✅ ${message}`);
        }
    };

    const assertEqual = (actual, expected, message) => {
        const sActual = actual.toString();
        const sExpected = expected.toString();
        if (sActual !== sExpected) {
            console.error(`❌ ${message}: Expected ${sExpected}, got ${sActual}`);
            throw new Error(`Expected ${sExpected}, got ${sActual}`);
        } else {
            console.log(`✅ ${message}`);
        }
    };

    console.log("--- Testing Vector Calculus (Line Integral) ---");

    // 1. Line Integral of F=[y, x] along r=[t, t] from 0 to 1
    // F(r(t)) = [t, t]. r'(t) = [1, 1]. Dot = t+t = 2t. Int(2t, 0, 1) = [t^2] = 1.
    const field = new Vec([new Sym('y'), new Sym('x')]);
    const vars = new Vec([new Sym('x'), new Sym('y')]);
    const curve = new Vec([new Sym('t'), new Sym('t')]);
    const res1 = cas._lineIntegral(field, vars, curve, new Sym('t'), new Num(0), new Num(1));
    assertEqual(res1, new Num(1), "Line Integral of [y,x] along t,t is 1");

    // 2. Line Integral of F=[1, 2] along r=[t, t^2] from 0 to 1
    // F=[1, 2]. r'=[1, 2t]. Dot = 1 + 4t. Int(1+4t) = t + 2t^2. From 0 to 1 = 3.
    const field2 = new Vec([new Num(1), new Num(2)]);
    const curve2 = new Vec([new Sym('t'), new Pow(new Sym('t'), new Num(2))]);
    const res2 = cas._lineIntegral(field2, vars, curve2, new Sym('t'), new Num(0), new Num(1));
    assertEqual(res2, new Num(3), "Line Integral of [1,2] along t,t^2 is 3");

    console.log("--- Testing Number Theory (Primitive Root) ---");

    // 3. primitiveRoot(7) should be 3.
    const pr7 = cas._primitiveRoot(new Num(7));
    assertEqual(pr7, new Num(3), "Primitive root of 7 is 3");

    // 4. primitiveRoot(1) -> 0
    assertEqual(cas._primitiveRoot(new Num(1)), new Num(0), "Primitive root of 1 is 0");

    console.log("--- Testing Statistics (Hypergeometric) ---");

    // 5. Hypergeometric PDF(k=1, M=5, n=2, N=10)
    // Draw 2 from 10, 5 are successes. Prob of 1 success.
    // comb(5,1)*comb(5,1) / comb(10,2) = 5*5 / 45 = 25/45 = 5/9.
    const hyp1 = cas._hypergeometricPDF(new Num(1), new Num(5), new Num(2), new Num(10));
    // 5/9 approx 0.55555...
    assert(Math.abs(hyp1.evaluateNumeric() - 5/9) < 1e-9, "Hypergeometric PDF numeric check");

    // 6. Hypergeometric CDF(k=2, M=5, n=2, N=10)
    // Sum P(X=0) + P(X=1) + P(X=2) = 1 (since max successes is 2)
    const hypCDF = cas._hypergeometricCDF(new Num(2), new Num(5), new Num(2), new Num(10));
    assert(Math.abs(hypCDF.evaluateNumeric() - 1) < 1e-9, "Hypergeometric CDF total prob is 1");

    console.log("--- All New Features Tests Passed ---");
};

test_new_features_5();
