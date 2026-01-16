
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

function runPerfTest() {
    console.log("Running Inverse Performance Test...");
    const size = 6;
    console.log(`Matrix Size: ${size}x${size}`);

    // Create matrix
    const rows = [];
    for(let i=0; i<size; i++) {
        const row = [];
        for(let j=0; j<size; j++) {
            row.push(new sandbox.Num(Math.floor(Math.random() * 10)));
        }
        rows.push(new sandbox.Vec(row));
    }
    const matrix = new sandbox.Vec(rows);

    const start = Date.now();
    try {
        cas._inv(matrix);
        const end = Date.now();
        console.log(`Time: ${end - start}ms`);
    } catch(e) {
        console.log("Error:", e.message);
    }
}

function runLimitTest() {
    console.log("Running Limit Test...");
    // limit(1/x, x, Infinity) -> 0
    try {
        const parser = new sandbox.Parser(new sandbox.Lexer("limit(1/x, x, infinity)"));
        const expr = parser.parse();
        const res = cas.evaluate(expr);
        console.log("limit(1/x, x, infinity) =", res.toString());
    } catch(e) {
        console.log("Limit Error:", e.message);
    }
}

runPerfTest();
runLimitTest();
