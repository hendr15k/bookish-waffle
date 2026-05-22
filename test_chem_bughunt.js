const { Chemistry } = require('./js/chemistry');
globalThis.Chemistry = Chemistry;
require('./js/expression');
require('./js/parser');
const { CAS } = require('./js/cas');
const cas = new CAS();

let passed = 0;
let failed = 0;

function assert(desc, actual, expected, tolerance = 1e-9) {
    if (Math.abs(actual - expected) <= tolerance) {
        console.log(`[PASS] ${desc}`);
        passed++;
    } else {
        console.error(`[FAIL] ${desc}: expected ${expected}, got ${actual}`);
        failed++;
    }
}

function assertEq(desc, actual, expected) {
    if (actual === expected) {
        console.log(`[PASS] ${desc}`);
        passed++;
    } else {
        console.error(`[FAIL] ${desc}: expected ${expected}, got ${actual}`);
        failed++;
    }
}

function assertThrow(desc, fn) {
    try {
        fn();
        console.error(`[FAIL] ${desc}: expected error, got none`);
        failed++;
    } catch(e) {
        console.log(`[PASS] ${desc} (threw: ${e.message})`);
        passed++;
    }
}

console.log("\n=== Molar Mass Tests ===");
assert("H2O", cas._molarMass("H2O").value, 18.015, 0.01);
assert("NaCl", cas._molarMass("NaCl").value, 58.44, 0.01);
assert("Ca(OH)2", cas._molarMass("Ca(OH)2").value, 74.092, 0.01);
assert("H2SO4", cas._molarMass("H2SO4").value, 98.079, 0.01);
assert("C6H12O6", cas._molarMass("C6H12O6").value, 180.156, 0.1);
assert("CO2", cas._molarMass("CO2").value, 44.01, 0.01);
assert("CH4", cas._molarMass("CH4").value, 16.04, 0.01);
assert("NH3", cas._molarMass("NH3").value, 17.031, 0.01);
assert("HCl", cas._molarMass("HCl").value, 36.46, 0.01);
assert("Fe2O3", cas._molarMass("Fe2O3").value, 159.69, 0.1);

console.log("\n=== Atomic Weight Tests ===");
assert("H", cas._atomicWeight("H").value, 1.008, 0.001);
assert("O", cas._atomicWeight("O").value, 15.999, 0.001);
assert("Fe", cas._atomicWeight("Fe").value, 55.845, 0.001);
assert("Na", cas._atomicWeight("Na").value, 22.990, 0.001);
assert("Au", cas._atomicWeight("Au").value, 196.97, 0.001);

console.log("\n=== Edge Case Tests ===");
assertThrow("Unknown element Xx throws", () => cas._molarMass("Xx"));
assertThrow("Empty formula throws", () => cas._molarMass(""));

assertEq("parseMolecule H2O", JSON.stringify(Chemistry.parseMolecule("H2O")), JSON.stringify({H:2,O:1}));
assertEq("parseMolecule Ca(OH)2", JSON.stringify(Chemistry.parseMolecule("Ca(OH)2")), JSON.stringify({Ca:1,O:2,H:2}));
assertEq("parseMolecule Mg(OH)2", JSON.stringify(Chemistry.parseMolecule("Mg(OH)2")), JSON.stringify({Mg:1,O:2,H:2}));
assertEq("parseMolecule Na2SO4", JSON.stringify(Chemistry.parseMolecule("Na2SO4")), JSON.stringify({Na:2,S:1,O:4}));
assertEq("parseMolecule Fe2O3", JSON.stringify(Chemistry.parseMolecule("Fe2O3")), JSON.stringify({Fe:2,O:3}));
assertEq("parseMolecule C6H12O6", JSON.stringify(Chemistry.parseMolecule("C6H12O6")), JSON.stringify({C:6,H:12,O:6}));

console.log("\n=== Balancer Tests ===");
const { Parser, Lexer } = require('./js/parser');
function evalExpr(text) {
    const parser = new Parser(new Lexer(text));
    return cas.evaluate(parser.parse());
}
function testBalance(eq, expected) {
    try {
        const result = evalExpr(`balance(${eq})`);
        if (result && result.text) {
            const text = result.text;
            console.log(`balance(${eq}) => ${text}`);
            if (expected && !text.includes(expected)) {
                console.error(`  [FAIL] expected to contain: ${expected}`);
                failed++;
            } else {
                console.log(`  [PASS]`);
                passed++;
            }
        } else {
            console.error(`[FAIL] balance(${eq}): no text in result`);
            failed++;
        }
    } catch(e) {
        console.error(`[FAIL] balance(${eq}) threw: ${e.message}`);
        failed++;
    }
}

testBalance("H2 + O2 -> H2O", "2H2 + O2 -> 2H2O");
testBalance("Na + Cl2 -> NaCl", "2Na + Cl2 -> 2NaCl");
testBalance("Fe + O2 -> Fe2O3", "4Fe + 3O2 -> 2Fe2O3");
testBalance("CH4 + O2 -> CO2 + H2O", "CH4 + 2O2 -> CO2 + 2H2O");
testBalance("N2 + H2 -> NH3", "N2 + 3H2 -> 2NH3");
testBalance("Ca(OH)2 + HCl -> CaCl2 + H2O", "Ca(OH)2 + 2HCl -> CaCl2 + 2H2O");

console.log(`\n=== Summary: ${passed} passed, ${failed} failed ===`);
process.exit(failed > 0 ? 1 : 0);
