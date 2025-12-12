
const PeriodicTable = {
    'H': { name: 'Hydrogen', mass: 1.008, num: 1 },
    'He': { name: 'Helium', mass: 4.0026, num: 2 },
    'Li': { name: 'Lithium', mass: 6.94, num: 3 },
    'Be': { name: 'Beryllium', mass: 9.0122, num: 4 },
    'B': { name: 'Boron', mass: 10.81, num: 5 },
    'C': { name: 'Carbon', mass: 12.011, num: 6 },
    'N': { name: 'Nitrogen', mass: 14.007, num: 7 },
    'O': { name: 'Oxygen', mass: 15.999, num: 8 },
    'F': { name: 'Fluorine', mass: 18.998, num: 9 },
    'Ne': { name: 'Neon', mass: 20.180, num: 10 },
    'Na': { name: 'Sodium', mass: 22.990, num: 11 },
    'Mg': { name: 'Magnesium', mass: 24.305, num: 12 },
    'Al': { name: 'Aluminium', mass: 26.982, num: 13 },
    'Si': { name: 'Silicon', mass: 28.085, num: 14 },
    'P': { name: 'Phosphorus', mass: 30.974, num: 15 },
    'S': { name: 'Sulfur', mass: 32.06, num: 16 },
    'Cl': { name: 'Chlorine', mass: 35.45, num: 17 },
    'K': { name: 'Potassium', mass: 39.098, num: 19 },
    'Ar': { name: 'Argon', mass: 39.948, num: 18 },
    'Ca': { name: 'Calcium', mass: 40.078, num: 20 },
    'Sc': { name: 'Scandium', mass: 44.956, num: 21 },
    'Ti': { name: 'Titanium', mass: 47.867, num: 22 },
    'V': { name: 'Vanadium', mass: 50.942, num: 23 },
    'Cr': { name: 'Chromium', mass: 51.996, num: 24 },
    'Mn': { name: 'Manganese', mass: 54.938, num: 25 },
    'Fe': { name: 'Iron', mass: 55.845, num: 26 },
    'Ni': { name: 'Nickel', mass: 58.693, num: 28 },
    'Co': { name: 'Cobalt', mass: 58.933, num: 27 },
    'Cu': { name: 'Copper', mass: 63.546, num: 29 },
    'Zn': { name: 'Zinc', mass: 65.38, num: 30 },
    'Ga': { name: 'Gallium', mass: 69.723, num: 31 },
    'Ge': { name: 'Germanium', mass: 72.63, num: 32 },
    'As': { name: 'Arsenic', mass: 74.922, num: 33 },
    'Se': { name: 'Selenium', mass: 78.96, num: 34 },
    'Br': { name: 'Bromine', mass: 79.904, num: 35 },
    'Kr': { name: 'Krypton', mass: 83.798, num: 36 },
    'Rb': { name: 'Rubidium', mass: 85.468, num: 37 },
    'Sr': { name: 'Strontium', mass: 87.62, num: 38 },
    'Y': { name: 'Yttrium', mass: 88.906, num: 39 },
    'Zr': { name: 'Zirconium', mass: 91.224, num: 40 },
    'Nb': { name: 'Niobium', mass: 92.906, num: 41 },
    'Mo': { name: 'Molybdenum', mass: 95.95, num: 42 },
    'Tc': { name: 'Technetium', mass: 98, num: 43 },
    'Ru': { name: 'Ruthenium', mass: 101.07, num: 44 },
    'Rh': { name: 'Rhodium', mass: 102.91, num: 45 },
    'Pd': { name: 'Palladium', mass: 106.42, num: 46 },
    'Ag': { name: 'Silver', mass: 107.87, num: 47 },
    'Cd': { name: 'Cadmium', mass: 112.41, num: 48 },
    'In': { name: 'Indium', mass: 114.82, num: 49 },
    'Sn': { name: 'Tin', mass: 118.71, num: 50 },
    'Sb': { name: 'Antimony', mass: 121.76, num: 51 },
    'I': { name: 'Iodine', mass: 126.90, num: 53 },
    'Te': { name: 'Tellurium', mass: 127.60, num: 52 },
    'Xe': { name: 'Xenon', mass: 131.29, num: 54 },
    'Cs': { name: 'Cesium', mass: 132.91, num: 55 },
    'Ba': { name: 'Barium', mass: 137.33, num: 56 },
    'La': { name: 'Lanthanum', mass: 138.91, num: 57 },
    'Ce': { name: 'Cerium', mass: 140.12, num: 58 },
    'Pr': { name: 'Praseodymium', mass: 140.91, num: 59 },
    'Nd': { name: 'Neodymium', mass: 144.24, num: 60 },
    'Pm': { name: 'Promethium', mass: 145, num: 61 },
    'Sm': { name: 'Samarium', mass: 150.36, num: 62 },
    'Eu': { name: 'Europium', mass: 151.96, num: 63 },
    'Gd': { name: 'Gadolinium', mass: 157.25, num: 64 },
    'Tb': { name: 'Terbium', mass: 158.93, num: 65 },
    'Dy': { name: 'Dysprosium', mass: 162.50, num: 66 },
    'Ho': { name: 'Holmium', mass: 164.93, num: 67 },
    'Er': { name: 'Erbium', mass: 167.26, num: 68 },
    'Tm': { name: 'Thulium', mass: 168.93, num: 69 },
    'Yb': { name: 'Ytterbium', mass: 173.05, num: 70 },
    'Lu': { name: 'Lutetium', mass: 174.97, num: 71 },
    'Hf': { name: 'Hafnium', mass: 178.49, num: 72 },
    'Ta': { name: 'Tantalum', mass: 180.95, num: 73 },
    'W': { name: 'Tungsten', mass: 183.84, num: 74 },
    'Re': { name: 'Rhenium', mass: 186.21, num: 75 },
    'Os': { name: 'Osmium', mass: 190.23, num: 76 },
    'Ir': { name: 'Iridium', mass: 192.22, num: 77 },
    'Pt': { name: 'Platinum', mass: 195.08, num: 78 },
    'Au': { name: 'Gold', mass: 196.97, num: 79 },
    'Hg': { name: 'Mercury', mass: 200.59, num: 80 },
    'Tl': { name: 'Thallium', mass: 204.38, num: 81 },
    'Pb': { name: 'Lead', mass: 207.2, num: 82 },
    'Bi': { name: 'Bismuth', mass: 208.98, num: 83 },
    'Po': { name: 'Polonium', mass: 209, num: 84 },
    'At': { name: 'Astatine', mass: 210, num: 85 },
    'Rn': { name: 'Radon', mass: 222, num: 86 },
    'Fr': { name: 'Francium', mass: 223, num: 87 },
    'Ra': { name: 'Radium', mass: 226, num: 88 },
    'Ac': { name: 'Actinium', mass: 227, num: 89 },
    'Th': { name: 'Thorium', mass: 232.04, num: 90 },
    'Pa': { name: 'Protactinium', mass: 231.04, num: 91 },
    'U': { name: 'Uranium', mass: 238.03, num: 92 },
    'Np': { name: 'Neptunium', mass: 237, num: 93 },
    'Pu': { name: 'Plutonium', mass: 244, num: 94 },
    'Am': { name: 'Americium', mass: 243, num: 95 },
    'Cm': { name: 'Curium', mass: 247, num: 96 },
    'Bk': { name: 'Berkelium', mass: 247, num: 97 },
    'Cf': { name: 'Californium', mass: 251, num: 98 },
    'Es': { name: 'Einsteinium', mass: 252, num: 99 },
    'Fm': { name: 'Fermium', mass: 257, num: 100 },
    'Md': { name: 'Mendelevium', mass: 258, num: 101 },
    'No': { name: 'Nobelium', mass: 259, num: 102 },
    'Lr': { name: 'Lawrencium', mass: 262, num: 103 },
    'Rf': { name: 'Rutherfordium', mass: 267, num: 104 },
    'Db': { name: 'Dubnium', mass: 268, num: 105 },
    'Sg': { name: 'Seaborgium', mass: 269, num: 106 },
    'Bh': { name: 'Bohrium', mass: 270, num: 107 },
    'Hs': { name: 'Hassium', mass: 269, num: 108 },
    'Mt': { name: 'Meitnerium', mass: 278, num: 109 },
    'Ds': { name: 'Darmstadtium', mass: 281, num: 110 },
    'Rg': { name: 'Roentgenium', mass: 281, num: 111 },
    'Cn': { name: 'Copernicium', mass: 285, num: 112 },
    'Nh': { name: 'Nihonium', mass: 284, num: 113 },
    'Fl': { name: 'Flerovium', mass: 289, num: 114 },
    'Mc': { name: 'Moscovium', mass: 288, num: 115 },
    'Lv': { name: 'Livermorium', mass: 293, num: 116 },
    'Ts': { name: 'Tennessine', mass: 294, num: 117 },
    'Og': { name: 'Oganesson', mass: 294, num: 118 }
};

class Chemistry {
    static getElement(symbol) {
        // Handle case sensitivity loosely
        let sym = symbol;
        if (PeriodicTable[sym]) return PeriodicTable[sym];

        // Try strict case: First Upper, rest lower
        sym = sym.charAt(0).toUpperCase() + sym.slice(1).toLowerCase();
        if (PeriodicTable[sym]) return PeriodicTable[sym];

        return null;
    }

    static calculateMolarMass(formula) {
        // Simple parser for chemical formulas e.g. H2O, C6H12O6, Ca(OH)2
        // Tokenize
        // Groups: Element, Number, (, )

        let tokens = [];
        let i = 0;
        const len = formula.length;

        while(i < len) {
            const char = formula[i];

            if (char === '(' || char === ')') {
                tokens.push({ type: 'bracket', val: char });
                i++;
                continue;
            }

            if (/[A-Z]/.test(char)) {
                // Element start
                let elem = char;
                i++;
                while(i < len && /[a-z]/.test(formula[i])) {
                    elem += formula[i];
                    i++;
                }
                tokens.push({ type: 'element', val: elem });
                continue;
            }

            if (/[0-9]/.test(char)) {
                let num = char;
                i++;
                while(i < len && /[0-9]/.test(formula[i])) {
                    num += formula[i];
                    i++;
                }
                tokens.push({ type: 'number', val: parseInt(num) });
                continue;
            }

            // Unknown char (ignore or throw?)
            i++;
        }

        // Parse tokens
        // Stack for handling brackets
        // We accumulate mass in current scope

        let stack = [{ mass: 0 }];

        for(let j=0; j<tokens.length; j++) {
            const token = tokens[j];

            if (token.type === 'element') {
                const el = this.getElement(token.val);
                if (!el) throw new Error("Unknown element: " + token.val);

                let count = 1;
                // Check if next is number
                if (j+1 < tokens.length && tokens[j+1].type === 'number') {
                    count = tokens[j+1].val;
                    j++; // Skip number
                }

                stack[stack.length-1].mass += el.mass * count;
            } else if (token.type === 'bracket') {
                if (token.val === '(') {
                    stack.push({ mass: 0 });
                } else {
                    if (stack.length === 1) throw new Error("Unmatched bracket");
                    const group = stack.pop();

                    let count = 1;
                    if (j+1 < tokens.length && tokens[j+1].type === 'number') {
                        count = tokens[j+1].val;
                        j++;
                    }

                    stack[stack.length-1].mass += group.mass * count;
                }
            }
        }

        if (stack.length !== 1) throw new Error("Unmatched bracket");
        return stack[0].mass;
    }
}

// Export
(function() {
    const exports = { PeriodicTable, Chemistry };
    if (typeof globalThis !== 'undefined') {
        Object.assign(globalThis, exports);
    }
    if (typeof module !== 'undefined' && module.exports) {
        module.exports = exports;
    }
})();
