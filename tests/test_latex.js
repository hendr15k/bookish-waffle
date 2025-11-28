
const fs = require('fs');
const vm = require('vm');

function loadFile(filePath) {
    return fs.readFileSync(filePath, 'utf8');
}

const expressionCode = loadFile('js/expression.js');
const parserCode = loadFile('js/parser.js');
const casCode = loadFile('js/cas.js');

const sandbox = {
    console: console,
    Math: Math,
    Number: Number,
    parseFloat: parseFloat,
    parseInt: parseInt
};

vm.createContext(sandbox);

// Helper to expose classes
const exposeClasses = `
    globalThis.Expr = Expr;
    globalThis.Num = Num;
    globalThis.Sym = Sym;
    globalThis.Add = Add;
    globalThis.Sub = Sub;
    globalThis.Mul = Mul;
    globalThis.Div = Div;
    globalThis.Pow = Pow;
    globalThis.Call = Call;
    globalThis.Assignment = Assignment;
    globalThis.Eq = Eq;
    globalThis.Vec = Vec;
    globalThis.FunctionDef = FunctionDef;

    globalThis.Lexer = Lexer;
    globalThis.Parser = Parser;

    globalThis.CAS = CAS;
`;

vm.runInContext(expressionCode + "\n" + parserCode + "\n" + casCode + "\n" + exposeClasses, sandbox);

const { Call, Sym, Num, Mul } = sandbox;

function checkDirectLatex(description, exprObj, expectedSubstring) {
    try {
        const latex = exprObj.toLatex();
        if (latex.includes(expectedSubstring) || latex === expectedSubstring) {
            console.log(`[PASS] ${description}: ${latex}`);
        } else {
            console.error(`[FAIL] ${description}`);
            console.error(`  Expected substring: ${expectedSubstring}`);
            console.error(`  Got:      ${latex}`);
        }
    } catch (e) {
        console.error(`[FAIL] ${description}`);
        console.error(`  Error:    ${e.message}`);
    }
}

console.log("--- Testing LaTeX Output ---");

// Greek letters
checkDirectLatex("Greek alpha", new Sym("alpha"), "\\alpha");
checkDirectLatex("Greek Gamma", new Sym("Gamma"), "\\Gamma");

// Calculus formatting - Test Call objects directly to avoid CAS evaluation/simplification
// limit(f(x), x, 0)
const limitCall = new Call("limit", [
    new Call("f", [new Sym("x")]),
    new Sym("x"),
    new Num(0)
]);
checkDirectLatex("Limit", limitCall, "\\lim_{x \\to 0}");

// integrate(f(x), x)
const integrateCall = new Call("integrate", [
    new Call("f", [new Sym("x")]),
    new Sym("x")
]);
checkDirectLatex("Integral", integrateCall, "\\int");

// integrate(f(x), x, 0, 1)
const defIntegrateCall = new Call("integrate", [
    new Call("f", [new Sym("x")]),
    new Sym("x"),
    new Num(0),
    new Num(1)
]);
checkDirectLatex("Definite Integral", defIntegrateCall, "\\int_{0}^{1}");

// sum(k, k, 1, n)
const sumCall = new Call("sum", [
    new Sym("k"),
    new Sym("k"),
    new Num(1),
    new Sym("n")
]);
checkDirectLatex("Sum", sumCall, "\\sum_{k=1}^{n}");

// diff(f(x), x)
const diffCall = new Call("diff", [
    new Call("f", [new Sym("x")]),
    new Sym("x")
]);
checkDirectLatex("Diff", diffCall, "\\frac{d}{d x}");

// Explicit Mul NumNum (only if we construct it directly and it doesn't simplify)
// NOTE: Mul constructor does NOT simplify automatically in the class definition in js/expression.js.
// `new Mul(...)` creates the object. `simplify()` is a method.
const mulCall = new Mul(new Num(2), new Num(3));
checkDirectLatex("Explicit Mul NumNum", mulCall, "2 \\cdot 3");

// Implicit Mul NumVar
const impMul1 = new Mul(new Num(2), new Sym("x"));
checkDirectLatex("Implicit Mul NumVar", impMul1, "2x");

// Implicit Mul Vars
const impMul2 = new Mul(new Sym("x"), new Sym("y"));
checkDirectLatex("Implicit Mul Vars", impMul2, "x y");

