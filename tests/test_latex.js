
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
    globalThis.Symbol = Symbol;
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

const cas = new sandbox.CAS();
const parser = (text) => new sandbox.Parser(new sandbox.Lexer(text));

function evalExpr(text) {
    const expr = parser(text).parse();
    // We want to test toLatex() of the expression tree as parsed, or evaluated?
    // User sees evaluated result.
    const result = cas.evaluate(expr);
    return result;
}

function checkLatex(description, input, expectedSubstring) {
    try {
        const result = evalExpr(input);
        const latex = result.toLatex();
        if (latex.includes(expectedSubstring) || latex === expectedSubstring) {
            console.log(`[PASS] ${description}: ${input} -> ${latex}`);
        } else {
            console.error(`[FAIL] ${description}`);
            console.error(`  Input: ${input}`);
            console.error(`  Expected substring: ${expectedSubstring}`);
            console.error(`  Got:      ${latex}`);
        }
    } catch (e) {
        console.error(`[FAIL] ${description}`);
        console.error(`  Input: ${input}`);
        console.error(`  Error:    ${e.message}`);
    }
}

console.log("--- Testing LaTeX Output ---");

// Greek letters
checkLatex("Greek alpha", "alpha", "\\alpha");
checkLatex("Greek Gamma", "Gamma", "\\Gamma");
checkLatex("Infinity", "limit(1/x, x, 0)", "\\infty"); // Assuming 1/x limit returns Infinity

// Calculus formatting
checkLatex("Integral", "integrate(f(x), x)", "\\int");
checkLatex("Definite Integral", "integrate(f(x), x, 0, 1)", "\\int_{0}^{1}");
checkLatex("Limit", "limit(f(x), x, 0)", "\\lim_{x \\to 0}");
checkLatex("Sum", "sum(k, k, 1, n)", "\\sum_{k=1}^{n}");
checkLatex("Product", "product(k, k, 1, n)", "\\prod_{k=1}^{n}");
checkLatex("Diff", "diff(f(x), x)", "\\frac{d}{d x}");

// Multiplication
checkLatex("Implicit Mul NumVar", "2 * x", "2x"); // Prefer 2x over 2 \cdot x
checkLatex("Implicit Mul Vars", "x * y", "x y"); // Prefer x y or x \cdot y. Let's aim for x y or just check it's not awful.
checkLatex("Explicit Mul NumNum", "2 * 3", "2 \\cdot 3");
