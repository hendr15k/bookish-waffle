const fs = require("fs");
const vm = require("vm");
function loadFile(filePath) { return fs.readFileSync(filePath, "utf8"); }
const expressionCode = loadFile("js/expression.js");
const parserCode = loadFile("js/parser.js");
const casCode = loadFile("js/cas.js");
const helpCode = loadFile("js/help.js");
const sandbox = { console, Math, Number, parseFloat, parseInt, isNaN, window: {} };
vm.createContext(sandbox);
vm.runInContext(expressionCode + "\n" + parserCode + "\n" + helpCode + "\n" + casCode + "\n" + `globalThis.Expr=Expr;globalThis.Num=Num;globalThis.Sym=Sym;globalThis.Add=Add;globalThis.Sub=Sub;globalThis.Mul=Mul;globalThis.Div=Div;globalThis.Pow=Pow;globalThis.Call=Call;globalThis.Assignment=Assignment;globalThis.Eq=Eq;globalThis.Vec=Vec;globalThis.FunctionDef=FunctionDef;globalThis.Lexer=Lexer;globalThis.Parser=Parser;globalThis.CAS=CAS;`, sandbox);
const cas = new sandbox.CAS();
const parser = (text) => new sandbox.Parser(new sandbox.Lexer(text));

try {
  const parsed = parser("odeplot(rk45(-0.5*y, y, t, 1, 0, 2, 0.1, 1e-6, []))").parse();
  const solVec = cas.evaluate(parsed.args[0]);
  console.log("sol evaluated length", solVec.elements.length);
  const firstRow = solVec.elements[0];
  console.log("firstRow.elements length", firstRow.elements.length);
  for (let i = 0; i < firstRow.elements.length; i++) {
    console.log(`firstRow[${i}]: ${firstRow.elements[i].constructor.name} - ${firstRow.elements[i].value}`);
  }

  const lastRow = solVec.elements[solVec.elements.length - 1];
  console.log("lastRow.elements length", lastRow.elements.length);
  for (let i = 0; i < lastRow.elements.length; i++) {
    console.log(`lastRow[${i}]: ${lastRow.elements[i].constructor.name} - ${lastRow.elements[i].value} / ${lastRow.elements[i].name}`);
  }

  // odeplot creates plot object: ODE dimension mismatch
  // The error comes from inside CAS.evaluate calling _recursiveEval and evaluating `rk45`. Wait, but I am doing `odeplot(rk45(...))`!
  // Let's trace it: `cas.evaluate(parsed)` -> `cas._recursiveEval(odeplot(...))`
  // inside `_recursiveEval`, it checks `if (node instanceof Call)`
  // it sees `node.funcName === 'odeplot'`. `node.args` is evaluated!
  // Wait, no: it says:
  // if (node.funcName === 'rk4' || node.funcName === 'rk45' || node.funcName === 'odeplot' || node.funcName === 'odetable') {
  //     return this._solveODE(node.funcName, node.args);
  // }
} catch (e) {
  console.error("Error:", e.stack);
}
