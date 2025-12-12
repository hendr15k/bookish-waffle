const { CAS } = require('../js/cas');
const { Chemistry } = require('../js/chemistry');
// Mock globalThis.Chemistry for CAS
globalThis.Chemistry = Chemistry;

// We need to import Expr and Num classes manually if they are not exposed by cas.js
// Wait, cas.js doesn't export Expr/Num?
// It does: if (typeof globalThis !== 'undefined') Object.assign(globalThis, exports);
// But we need to load expression.js first.
require('../js/expression');
require('../js/parser');

const cas = new CAS();

function test_molarMass() {
    const formulas = [
        ['H2O', 18.015],
        ['C6H12O6', 180.156],
        ['NaCl', 58.44], // Na=22.99, Cl=35.45
        ['Ca(OH)2', 74.093], // Ca=40.078, O=15.999, H=1.008 -> 40.078 + 2*(15.999+1.008) = 74.092
    ];

    for (const [f, expected] of formulas) {
        // We test _molarMass direct call first
        try {
            const mass = cas._molarMass(f);
            if (Math.abs(mass.value - expected) > 0.1) {
                console.error(`FAIL: Molar Mass ${f}. Expected ${expected}, got ${mass.value}`);
            } else {
                console.log(`[PASS] Molar Mass ${f}: ${mass.value}`);
            }
        } catch(e) {
            console.error(`FAIL: Molar Mass ${f} threw error: ${e.message}`);
        }
    }
}

function test_atomicWeight() {
    const els = [
        ['H', 1.008],
        ['O', 15.999],
        ['Au', 196.97]
    ];

    for (const [el, expected] of els) {
        try {
            const w = cas._atomicWeight(el);
            if (Math.abs(w.value - expected) > 0.001) {
                console.error(`FAIL: Atomic Weight ${el}. Expected ${expected}, got ${w.value}`);
            } else {
                console.log(`[PASS] Atomic Weight ${el}: ${w.value}`);
            }
        } catch(e) {
            console.error(`FAIL: Atomic Weight ${el} threw error: ${e.message}`);
        }
    }
}

test_molarMass();
test_atomicWeight();
