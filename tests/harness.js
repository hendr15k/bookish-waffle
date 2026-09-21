// Reusable CAS test harness: loads engine into a sandbox and exposes helpers.
const fs = require('fs');
const path = require('path');
const vm = require('vm');

const root = path.join(__dirname, '..');

const sandbox = {
    console: console,
    Math: Math,
    Number: Number,
    parseFloat: parseFloat,
    parseInt: parseInt
};

vm.createContext(sandbox);

const exposeClasses = `
    globalThis.Expr = Expr;
    globalThis.Num = Num;
    globalThis.Sym = Sym;
    globalThis.Symbol = Sym;
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

vm.runInContext(
    fs.readFileSync(path.join(root, 'js/expression.js'), 'utf8') + '\n' +
    fs.readFileSync(path.join(root, 'js/parser.js'), 'utf8') + '\n' +
    fs.readFileSync(path.join(root, 'js/cas.js'), 'utf8') + '\n' +
    exposeClasses, sandbox);

const cas = new sandbox.CAS();

function evalExpr(text) {
    const expr = new sandbox.Parser(new sandbox.Lexer(text)).parse();
    return cas.evaluate(expr);
}

function ev(text) {
    return evalExpr(text).toString();
}

// evalExpr that keeps the session (variables persist) but resets between groups
function resetCas() { cas.variables['x'] = undefined; delete cas.variables['x']; }

function test(description, input, expectedOutput) {
    try {
        const result = evalExpr(input);
        const resultStr = result.toString();
        if (resultStr === expectedOutput) {
            console.log(`[PASS] ${description}`);
            return true;
        }
        console.error(`[FAIL] ${description}`);
        console.error(`  Input: ${input}`);
        console.error(`  Expected: ${expectedOutput}`);
        console.error(`  Got:      ${resultStr}`);
        return false;
    } catch (e) {
        console.error(`[FAIL] ${description}`);
        console.error(`  Input: ${input}`);
        console.error(`  Error:    ${e.message}`);
        return false;
    }
}

function testShouldThrow(description, input, expectedErrorFragment) {
    try {
        evalExpr(input);
        console.error(`[FAIL] ${description} - Expected error but got success.`);
        return false;
    } catch (e) {
        if (e.message.includes(expectedErrorFragment)) {
            console.log(`[PASS] ${description}`);
            return true;
        }
        console.error(`[FAIL] ${description} - Expected error containing "${expectedErrorFragment}" but got "${e.message}"`);
        return false;
    }
}

module.exports = { sandbox, cas, evalExpr, ev, test, testShouldThrow, resetCas };
