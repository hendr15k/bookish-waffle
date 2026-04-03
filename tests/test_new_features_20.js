
const { CAS } = require('../js/cas.js');
const { Expr, Sym, Num, Add, Sub, Mul, Div, Pow, Call, Vec } = require('../js/expression.js');

global.Expr = Expr;
global.Sym = Sym;
global.Num = Num;
global.Call = Call;
global.Add = Add;
global.Sub = Sub;
global.Mul = Mul;
global.Div = Div;
global.Pow = Pow;
global.Vec = Vec;

const cas = new CAS();

function test(name, expr, expectedStr) {
    console.log(`Testing ${name}...`);
    try {
        const res = cas.evaluate(expr);
        const resStr = res.toString();

        if (expectedStr && resStr !== expectedStr) {
            // Normalization for parens and spaces
            const normRes = resStr.replace(/\s/g, '').replace(/\(/g, '').replace(/\)/g, '');
            const normExp = expectedStr.replace(/\s/g, '').replace(/\(/g, '').replace(/\)/g, '');

            if (normRes === normExp) {
                 console.log(`  PASSED (ignoring format): ${resStr}`);
            } else {
                 console.error(`  FAILED: Expected ${expectedStr}, got ${resStr}`);
            }
        } else {
            console.log(`  PASSED: ${resStr}`);
        }
    } catch (e) {
        console.error(`  ERROR: ${e.message}`);
        console.error(e.stack);
    }
}

const x = new Sym('x');
const one = new Num(1);

// 1. Residue simple pole: 1/(x-1) at 1
const f1 = new Div(one, new Sub(x, one));
test("Residue simple pole", new Call('residue', [f1, x, one]), "1");

// 2. Residue order 2: 1/(x-1)^2 at 1
const f2 = new Div(one, new Pow(new Sub(x, one), new Num(2)));
test("Residue order 2", new Call('residue', [f2, x, one]), "0");

// 3. Singularities: 1/(x^2-1)
const f3 = new Div(one, new Sub(new Pow(x, new Num(2)), one));
test("Singularities", new Call('singularities', [f3, x]), "[-1, 1]");

// 4. Continuity sin(x) at 0
test("isContinuous sin(x)", new Call('isContinuous', [new Call('sin', [x]), x, new Num(0)]), "1");

// 5. Continuity 1/x at 0
test("isContinuous 1/x", new Call('isContinuous', [new Div(one, x), x, new Num(0)]), "0");

// 6. Differentiability x^2 at 0
test("isDifferentiable x^2", new Call('isDifferentiable', [new Pow(x, new Num(2)), x, new Num(0)]), "1");

// 7. Differentiability abs(x) at 0
test("isDifferentiable abs(x)", new Call('isDifferentiable', [new Call('abs', [x]), x, new Num(0)]), "0");

// 8. Routh Hurwitz s^3 + 2s^2 + 3s + 1
const poly = new Add(new Add(new Pow(x, new Num(3)), new Mul(new Num(2), new Pow(x, new Num(2)))), new Add(new Mul(new Num(3), x), one));
test("Routh Hurwitz", new Call('routhHurwitz', [poly, x]), "[[1, 3], [2, 1], [5/2, 0], [1, 0]]");

// 9. Complex parts
const c = new Add(one, new Mul(new Num(2), new Sym('i')));
test("re(1+2i)", new Call('re', [c]), "1");
test("im(1+2i)", new Call('im', [c]), "2");
test("conj(1+2i)", new Call('conj', [c]), "1 - 2 * i");

// 10. Plot Contour
const y = new Sym('y');
const plotExpr = new Add(new Pow(x, new Num(2)), new Pow(y, new Num(2)));
const callPlot = new Call('plotcontour', [plotExpr, x, y]);
console.log("Testing plotcontour...");
try {
    const resPlot = cas.evaluate(callPlot);
    if (resPlot.type === 'plot' && resPlot.subtype === 'contour') {
        console.log("  PASSED: plotcontour returned contour plot object");
    } else {
        console.error("  FAILED: plotcontour returned " + JSON.stringify(resPlot));
    }
} catch(e) {
    console.error("  ERROR: " + e.message);
}
