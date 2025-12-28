const { Parser, Lexer } = require('../js/parser.js');
const { CAS } = require('../js/cas.js');
const { Call, Sym, Num, Eq, Vec } = require('../js/expression.js'); // Ensure required classes are imported

// Mock globalThis for tests
globalThis.Parser = Parser;
globalThis.Lexer = Lexer;
globalThis.CAS = CAS;
globalThis.Call = Call;
globalThis.Sym = Sym;
globalThis.Num = Num;
globalThis.Eq = Eq;
globalThis.Vec = Vec;

const cas = new CAS();

function testPlotNstep() {
    console.log("Testing plot with nstep...");

    // Test case 1: plot(sin(x), x, -5, 5, nstep=100)
    const input1 = "plot(sin(x), x, -5, 5, nstep=100)";
    const lexer1 = new Lexer(input1);
    const parser1 = new Parser(lexer1);
    const tree1 = parser1.parse();
    const result1 = cas.evaluate(tree1);

    if (result1.type === 'plot' && result1.nstep === 100) {
        console.log("PASS: plot nstep=100 detected correctly.");
    } else {
        console.error("FAIL: plot nstep=100 failed. Got:", result1.nstep);
    }

    // Test case 2: plotparam(cos(t), sin(t), t, 0, 2*pi, nstep=50)
    const input2 = "plotparam(cos(t), sin(t), t, 0, 6.28, nstep=50)";
    const lexer2 = new Lexer(input2);
    const parser2 = new Parser(lexer2);
    const tree2 = parser2.parse();
    const result2 = cas.evaluate(tree2);

    if (result2.type === 'plot' && result2.subtype === 'parametric' && result2.nstep === 50) {
        console.log("PASS: plotparam nstep=50 detected correctly.");
    } else {
        console.error("FAIL: plotparam nstep=50 failed. Got:", result2.nstep);
    }

    // Test case 3: plotpolar(1, t, 0, 2*pi, nstep=20)
    const input3 = "plotpolar(1, t, 0, 6.28, nstep=20)";
    const lexer3 = new Lexer(input3);
    const parser3 = new Parser(lexer3);
    const tree3 = parser3.parse();
    const result3 = cas.evaluate(tree3);

    if (result3.type === 'plot' && result3.subtype === 'polar' && result3.nstep === 20) {
        console.log("PASS: plotpolar nstep=20 detected correctly.");
    } else {
        console.error("FAIL: plotpolar nstep=20 failed. Got:", result3.nstep);
    }

    // Test case 4: Default behavior (no nstep)
    const input4 = "plot(x^2, x)";
    const lexer4 = new Lexer(input4);
    const parser4 = new Parser(lexer4);
    const tree4 = parser4.parse();
    const result4 = cas.evaluate(tree4);

    if (result4.type === 'plot' && result4.nstep === null) {
        console.log("PASS: plot default nstep is null.");
    } else {
        console.error("FAIL: plot default nstep failed. Got:", result4.nstep);
    }
}

try {
    testPlotNstep();
} catch (e) {
    console.error("ERROR running tests:", e);
}
