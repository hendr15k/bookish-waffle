
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
    Date: Date
};

vm.createContext(sandbox);

const exposeClasses = `
    globalThis.Expr = Expr;
    globalThis.Num = Num;
    globalThis.Sym = Sym;
    globalThis.Vec = Vec;
    globalThis.CAS = CAS;
    globalThis.Parser = Parser;
    globalThis.Lexer = Lexer;
    globalThis.Call = Call;
    globalThis.Add = Add;
    globalThis.Sub = Sub;
    globalThis.Mul = Mul;
    globalThis.Div = Div;
    globalThis.Pow = Pow;
`;

vm.runInContext(expressionCode + "\n" + parserCode + "\n" + casCode + "\n" + exposeClasses, sandbox);

const cas = new sandbox.CAS();

function checkSingularity() {
    console.log("Checking integrate(1/x^2, x, -1, 1)...");
    const parser = new sandbox.Parser(new sandbox.Lexer("integrate(1/x^2, x, -1, 1)"));
    const expr = parser.parse();
    const res = cas.evaluate(expr);
    console.log("Result:", res.toString());
}

checkSingularity();
