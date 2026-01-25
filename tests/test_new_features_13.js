const exprLib = require('../js/expression.js');
Object.assign(global, exprLib); // Expose Expr classes globally for CAS

const { CAS } = require('../js/cas.js');

const cas = new CAS();

// Helpers
const n = (v) => new Num(v);
const s = (v) => new Sym(v);
const c = (f, args) => new Call(f, args);
const v = (arr) => new Vec(arr);

console.log("Running Tests for New Features 13 (Number Theory)");

// 1. Chinese Remainder Theorem
// x = 2 mod 3, x = 3 mod 5 -> x = 8 mod 15
try {
    const res = cas.evaluate(c('chineseRemainder', [v([n(2), n(3)]), v([n(3), n(5)])]));
    console.log(`chineseRemainder([2, 3], [3, 5]) = ${res.toString()}`);
    if (res.toString() !== '8') console.error("FAIL: expected 8");
} catch (e) {
    console.error("FAIL chineseRemainder 1:", e);
}

// x = 2 mod 3, x = 3 mod 4, x = 1 mod 5
// x = 11 mod 60
try {
    const res = cas.evaluate(c('crt', [v([n(2), n(3), n(1)]), v([n(3), n(4), n(5)])]));
    console.log(`crt([2, 3, 1], [3, 4, 5]) = ${res.toString()}`);
    if (res.toString() !== '11') console.error("FAIL: expected 11");
} catch (e) {
    console.error("FAIL chineseRemainder 2:", e);
}

// 2. Extended GCD
// xgcd(240, 46) -> [2, -9, 47]
try {
    const res = cas.evaluate(c('xgcd', [n(240), n(46)]));
    console.log(`xgcd(240, 46) = ${res.toString()}`);
    if (res.elements[0].value !== 2) console.error("FAIL: gcd != 2");
    if (res.elements[1].value !== -9) console.error("FAIL: x != -9");
    if (res.elements[2].value !== 47) console.error("FAIL: y != 47");
} catch (e) {
    console.error("FAIL xgcd:", e);
}

// 3. isPrimitiveRoot
try {
    // 3 is primitive root mod 7? Yes.
    const res1 = cas.evaluate(c('isPrimitiveRoot', [n(3), n(7)]));
    console.log(`isPrimitiveRoot(3, 7) = ${res1.toString()}`);
    if (res1.value !== 1) console.error("FAIL: expected 1 for primitive root");

    // 2 is primitive root mod 7? No (order 3).
    const res2 = cas.evaluate(c('isPrimitiveRoot', [n(2), n(7)]));
    console.log(`isPrimitiveRoot(2, 7) = ${res2.toString()}`);
    if (res2.value !== 0) console.error("FAIL: expected 0 for non-primitive root");
} catch (e) {
    console.error("FAIL isPrimitiveRoot:", e);
}

console.log("Tests Completed");
