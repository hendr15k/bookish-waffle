const fs = require('fs');
const vm = require('vm');
const path = require('path');

const exprCode = fs.readFileSync(path.join(__dirname, '../js/expression.js'), 'utf8');
const casCode = fs.readFileSync(path.join(__dirname, '../js/cas.js'), 'utf8');

const sandbox = { console, Math };
vm.createContext(sandbox);
vm.runInContext(exprCode, sandbox);
vm.runInContext(casCode, sandbox);

function runTest(name, func) {
    try {
        func();
        console.log(`PASS: ${name}`);
    } catch (e) {
        console.error(`FAIL: ${name}`);
        console.error(e.message);
        process.exit(1);
    }
}

runTest('Combinations (Numeric)', () => {
    const cas = new sandbox.CAS();
    // combinations(5, 2)
    const call = new sandbox.Call('combinations', [new sandbox.Num(5), new sandbox.Num(2)]);
    const res = cas.evaluate(call);
    if (res.value !== 10) throw new Error(`Expected 10, got ${res.toString()}`);
});

runTest('Combinations (List)', () => {
    const cas = new sandbox.CAS();
    // combinations([1, 2, 3], 2)
    const vec = new sandbox.Vec([new sandbox.Num(1), new sandbox.Num(2), new sandbox.Num(3)]);
    const call = new sandbox.Call('combinations', [vec, new sandbox.Num(2)]);
    const res = cas.evaluate(call);
    if (!(res instanceof sandbox.Vec)) throw new Error("Expected Vec");
    if (res.elements.length !== 3) throw new Error(`Expected length 3, got ${res.elements.length}`);
    if (res.elements[0].toString() !== "[1, 2]") throw new Error(`Expected [1, 2], got ${res.elements[0].toString()}`);
});

runTest('Permutations (List)', () => {
    const cas = new sandbox.CAS();
    // permutations([1, 2], 2)
    const vec = new sandbox.Vec([new sandbox.Num(1), new sandbox.Num(2)]);
    const call = new sandbox.Call('permutations', [vec, new sandbox.Num(2)]);
    const res = cas.evaluate(call);
    if (res.elements.length !== 2) throw new Error(`Expected length 2, got ${res.elements.length}`);
    if (res.toString() !== "[[1, 2], [2, 1]]") throw new Error(`Expected [[1, 2], [2, 1]], got ${res.toString()}`);
});

runTest('Unique', () => {
    const cas = new sandbox.CAS();
    // unique([1, 2, 2, 3])
    const vec = new sandbox.Vec([new sandbox.Num(1), new sandbox.Num(2), new sandbox.Num(2), new sandbox.Num(3)]);
    const call = new sandbox.Call('unique', [vec]);
    const res = cas.evaluate(call);
    if (res.toString() !== "[1, 2, 3]") throw new Error(`Expected [1, 2, 3], got ${res.toString()}`);
});

runTest('Integrate sinc(x)', () => {
    const cas = new sandbox.CAS();
    const x = new sandbox.Sym('x');
    // integrate(sinc(x), x)
    const call = new sandbox.Call('integrate', [new sandbox.Call('sinc', [x]), x]);
    const res = cas.evaluate(call);
    if (res.toString() !== "Si(x)") throw new Error(`Expected Si(x), got ${res.toString()}`);
});

runTest('exp2trig', () => {
    const cas = new sandbox.CAS();
    const x = new sandbox.Sym('x');
    const i = new sandbox.Sym('i');
    // exp2trig(exp(i*x))
    const arg = new sandbox.Call('exp', [new sandbox.Mul(i, x)]);
    const call = new sandbox.Call('exp2trig', [arg]);
    const res = cas.evaluate(call);
    // (cos(x) + i*sin(x))
    const s = res.toString();
    console.log("exp2trig result: " + s);
    if (!s.includes('cos(x)') || !s.includes('sin(x)')) throw new Error("Expected cos(x) + i*sin(x) form");
});

runTest('trig2exp', () => {
    const cas = new sandbox.CAS();
    const x = new sandbox.Sym('x');
    // trig2exp(sin(x))
    const call = new sandbox.Call('trig2exp', [new sandbox.Call('sin', [x])]);
    const res = cas.evaluate(call);
    const s = res.toString();
    console.log("trig2exp result: " + s);
    if (!s.includes('exp')) throw new Error("Expected exponential form");
});
