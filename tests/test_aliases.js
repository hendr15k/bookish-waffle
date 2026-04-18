// Alias tests - verifies that all command aliases resolve to the same result
// Uses the same vm-based approach as test_core.js (avoids JSDOM async script loading issues)
const fs = require('fs');
const path = require('path');
const vm = require('vm');

function loadFile(filePath) {
    return fs.readFileSync(filePath, 'utf8');
}

const expressionCode = loadFile(path.join(__dirname, '..', 'js', 'expression.js'));
const parserCode = loadFile(path.join(__dirname, '..', 'js', 'parser.js'));
const casCode = loadFile(path.join(__dirname, '..', 'js', 'cas.js'));

const sandbox = {
    console: console,
    Math: Math,
    Number: Number,
    parseFloat: parseFloat,
    parseInt: parseInt,
    path: path,
    __dirname: __dirname
};

// Expose classes to globalThis so evaluate() can find them
const exposeCode = `
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

vm.createContext(sandbox);
vm.runInContext(expressionCode + "\n" + parserCode + "\n" + casCode + "\n" + exposeCode, sandbox);

const cas = new sandbox.CAS();
const parser = (text) => new sandbox.Parser(new sandbox.Lexer(text));

function assertEqual(actual, expected, label) {
    if (actual !== expected) {
        throw new Error(`${label} failed: expected ${expected}, got ${actual}`);
    }
    console.log(`[PASS] ${label}`);
}

function evalExpr(text) {
    return cas.evaluate(parser(text).parse());
}

assertEqual(evalExpr('egcd(240,46)').toString(), evalExpr('xgcd(240,46)').toString(), 'egcd alias');
assertEqual(evalExpr('invmod(3,11)').toString(), evalExpr('modInverse(3,11)').toString(), 'invmod alias');
assertEqual(evalExpr('inverse(3,11)').toString(), evalExpr('modInverse(3,11)').toString(), 'inverse alias');
assertEqual(evalExpr('powmod(2,10,1000)').toString(), evalExpr('modPow(2,10,1000)').toString(), 'powmod alias');
assertEqual(evalExpr('gcd(12,18)').toString(), '6', 'gcd basic');
assertEqual(evalExpr('isPrime(17)').toString(), '1', 'isPrime basic');
assertEqual(evalExpr('nCr(5,2)').toString(), evalExpr('binomial(5,2)').toString(), 'nCr/binomial alias');
assertEqual(evalExpr('nPr(5,2)').toString(), evalExpr('perm(5,2)').toString(), 'nPr/perm alias');
assertEqual(evalExpr('det([[1,2],[3,4]])').toString(), '-2', 'det');
assertEqual(evalExpr('trans([[1,2],[3,4]])').toString(), evalExpr('transpose([[1,2],[3,4]])').toString(), 'trans/transpose alias');
assertEqual(evalExpr('eye(3)').toString(), evalExpr('identity(3)').toString(), 'eye/identity alias');

console.log('All alias tests passed');
