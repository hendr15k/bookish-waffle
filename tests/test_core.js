
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
    return cas.evaluate(expr);
}

function test(description, input, expectedOutput) {
    try {
        const result = evalExpr(input);
        const resultStr = result.toString();
        if (resultStr === expectedOutput) {
            console.log(`[PASS] ${description}`);
        } else {
            console.error(`[FAIL] ${description}`);
            console.error(`  Input: ${input}`);
            console.error(`  Expected: ${expectedOutput}`);
            console.error(`  Got:      ${resultStr}`);
        }
    } catch (e) {
        console.error(`[FAIL] ${description}`);
        console.error(`  Input: ${input}`);
        console.error(`  Error:    ${e.message}`);
    }
}

// Basic arithmetic
test("Addition", "1 + 1", "2");
test("Subtraction", "5 - 3", "2");
test("Multiplication", "4 * 2", "8");
test("Division", "10 / 2", "5");
test("Power", "2 ^ 3", "8");

// Algebra
test("Variable", "x", "x");
test("Symbolic Addition", "x + x", "(2 * x)");
test("Symbolic Mul", "x * x", "(x * x)");

// Vectors
test("Vector Addition", "[1, 2] + [3, 4]", "[4, 6]");
test("Dot Product", "[1, 2] * [3, 4]", "11");
test("Matrix * Vector", "[[1, 2], [3, 4]] * [5, 6]", "[17, 39]");
