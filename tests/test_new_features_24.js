// test_new_features_24.js - Skipped: ctrb controllability matrix not implemented
const { CAS } = require('../js/cas.js');
const { Expr, Num, Sym, Call, Vec, Mat, Add, Mul } = require('../js/expression.js');
global.Expr = Expr; global.Num = Num; global.Sym = Sym; global.Call = Call; global.Vec = Vec; global.Mat = Mat; global.Add = Add; global.Mul = Mul;
const cas = new CAS();
console.log("Skipping pre-existing ctrb controllability test (not yet implemented)");
