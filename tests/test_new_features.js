const vm = require('vm');
const fs = require('fs');
const path = require('path');

// Load Core
const exprCode = fs.readFileSync(path.join(__dirname, '../js/expression.js'), 'utf8');
const casCode = fs.readFileSync(path.join(__dirname, '../js/cas.js'), 'utf8');

const sandbox = {
    console: console,
    Math: Math
};
vm.createContext(sandbox);
vm.runInContext(exprCode, sandbox);
vm.runInContext(casCode, sandbox);

// Helper
function test(name, fn) {
    try {
        fn();
        console.log(`[PASS] ${name}`);
    } catch (e) {
        console.error(`[FAIL] ${name}: ${e.message}`);
        process.exit(1);
    }
}

const cas = new sandbox.CAS();
const Num = sandbox.Num;

test('F-Distribution PDF', () => {
    // f(x, d1, d2)
    const x = new Num(1);
    const d1 = new Num(5);
    const d2 = new Num(10);
    const res = cas._fPDF(x, d1, d2);
    if (!(res instanceof Num)) throw new Error("Result is not a number");
    if (res.value <= 0) throw new Error("PDF should be positive");
});

test('F-Distribution CDF', () => {
    const x = new Num(1);
    const d1 = new Num(5);
    const d2 = new Num(10);
    const res = cas._fCDF(x, d1, d2);
    if (!(res instanceof Num)) throw new Error("Result is not a number");
    if (res.value < 0 || res.value > 1) throw new Error("CDF out of range");
});

test('Beta Function', () => {
    // B(1, 1) = 1
    const b = cas._beta(1, 1);
    if (Math.abs(b - 1) > 1e-9) throw new Error("B(1,1) != 1");
});
