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
    parseInt: parseInt,
    isNaN: isNaN
};

vm.createContext(sandbox);
vm.runInContext(expressionCode + "\n" + parserCode + "\n" + casCode + "\n" + `
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
`, sandbox);

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
            return true;
        } else {
            console.error(`[FAIL] ${description}`);
            console.error(`  Input: ${input}`);
            console.error(`  Expected: ${expectedOutput}`);
            console.error(`  Got:      ${resultStr}`);
            return false;
        }
    } catch (e) {
        console.error(`[FAIL] ${description}`);
        console.error(`  Input: ${input}`);
        console.error(`  Error:    ${e.message}`);
        return false;
    }
}

console.log('=== Testing BUG #6: multivariate factor ===');
test('factor(x*y)', 'factor(x*y)', '(x * y)');
test('factor(a*x+b*x)', 'factor(a*x+b*x)', '(x * (a + b))');
test('factor(x^2-y^2)', 'factor(x^2-y^2)', '((x - y) * (x + y))');

console.log('\n=== Testing that existing factor still works ===');
test('factor(x^2-1)', 'factor(x^2-1)', '((x - 1) * (x + 1))');
test('factor(x^2+1)', 'factor(x^2+1)', '(x^2 + 1)');
test('factor(x^4-1)', 'factor(x^4-1)', '(((x - 1) * (x + 1)) * (x^2 + 1))');
