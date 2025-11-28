
class Expr {
    add(other) { return new Add(this, toExpr(other)); }
    sub(other) { return new Sub(this, toExpr(other)); }
    mul(other) { return new Mul(this, toExpr(other)); }
    div(other) { return new Div(this, toExpr(other)); }
    pow(other) { return new Pow(this, toExpr(other)); }

    simplify() { return this; }
    evaluateNumeric() { return NaN; }
    diff(varName) { throw new Error("diff not implemented"); }
    integrate(varName) { throw new Error("integrate not implemented for " + this.constructor.name); }
    expand() { return this; }
    substitute(varName, value) { return this; }
    toLatex() { return this.toString(); }
}

function toExpr(other) {
    if (other instanceof Expr) return other;
    if (typeof other === 'number') return new Num(other);
    throw new Error("Cannot convert to Expr: " + other);
}

class Num extends Expr {
    constructor(value) {
        super();
        this.value = value;
    }
    toString() { return this.value.toString(); }
    simplify() { return this; }
    evaluateNumeric() { return this.value; }
    diff(varName) { return new Num(0); }
    integrate(varName) {
        if (this.value === 0) return new Num(0);
        return new Mul(this, varName);
    }
    substitute(varName, value) { return this; }
    toLatex() { return this.value.toString(); }
    equals(other) { return other instanceof Num && this.value === other.value; }
}

class Symbol extends Expr {
    constructor(name) {
        super();
        this.name = name;
    }
    toString() { return this.name; }
    simplify() { return this; }
    evaluateNumeric() {
        if (this.name === 'pi') return Math.PI;
        if (this.name === 'e') return Math.E;
        if (this.name === 'i') return NaN; // Or complex object?
        return NaN; // Cannot evaluate generic symbol
    }
    diff(varName) {
        return (this.name === varName.name) ? new Num(1) : new Num(0);
    }
    integrate(varName) {
        if (this.name === varName.name) {
            return new Div(new Pow(this, new Num(2)), new Num(2));
        } else {
            return new Mul(this, varName);
        }
    }
    substitute(varName, value) {
        return (this.name === varName.name) ? value : this;
    }
    toLatex() {
        const greek = {
            'alpha': '\\alpha', 'beta': '\\beta', 'gamma': '\\gamma', 'delta': '\\delta',
            'epsilon': '\\epsilon', 'zeta': '\\zeta', 'eta': '\\eta', 'theta': '\\theta',
            'iota': '\\iota', 'kappa': '\\kappa', 'lambda': '\\lambda', 'mu': '\\mu',
            'nu': '\\nu', 'xi': '\\xi', 'omicron': 'o', 'pi': '\\pi', 'rho': '\\rho',
            'sigma': '\\sigma', 'tau': '\\tau', 'upsilon': '\\upsilon', 'phi': '\\phi',
            'chi': '\\chi', 'psi': '\\psi', 'omega': '\\omega',
            'Alpha': 'A', 'Beta': 'B', 'Gamma': '\\Gamma', 'Delta': '\\Delta',
            'Epsilon': 'E', 'Zeta': 'Z', 'Eta': 'H', 'Theta': '\\Theta',
            'Iota': 'I', 'Kappa': 'K', 'Lambda': '\\Lambda', 'Mu': 'M',
            'Nu': 'N', 'Xi': '\\Xi', 'Omicron': 'O', 'Pi': '\\Pi', 'Rho': 'P',
            'Sigma': '\\Sigma', 'Tau': 'T', 'Upsilon': '\\Upsilon', 'Phi': '\\Phi',
            'Chi': 'X', 'Psi': '\\Psi', 'Omega': '\\Omega',
            'Infinity': '\\infty'
        };
        if (greek.hasOwnProperty(this.name)) return greek[this.name];
        return this.name;
    }
    equals(other) { return other instanceof Symbol && this.name === other.name; }
}

class BinaryOp extends Expr {
    constructor(left, right) {
        super();
        this.left = left;
        this.right = right;
    }
    substitute(varName, value) {
        return new this.constructor(this.left.substitute(varName, value), this.right.substitute(varName, value));
    }
}

class Add extends BinaryOp {
    toString() { return `(${this.left} + ${this.right})`; }
    simplify() {
        const l = this.left.simplify();
        const r = this.right.simplify();

        // Vector addition
        if (l instanceof Vec && r instanceof Vec) {
            if (l.elements.length !== r.elements.length) throw new Error("Vector length mismatch in addition");
            const newElements = [];
            for (let i = 0; i < l.elements.length; i++) {
                newElements.push(new Add(l.elements[i], r.elements[i]).simplify());
            }
            return new Vec(newElements);
        }

        if (l instanceof Num && r instanceof Num) return new Num(l.value + r.value);
        if (l instanceof Num && l.value === 0) return r;
        if (r instanceof Num && r.value === 0) return l;
        if (r instanceof Num && r.value < 0) return new Sub(l, new Num(-r.value)).simplify();
        if (l instanceof Num && l.value < 0) return new Sub(r, new Num(-l.value)).simplify(); // (-a) + b -> b - a
        if (l.toString() === r.toString()) return new Mul(new Num(2), l);
        return new Add(l, r);
    }
    evaluateNumeric() { return this.left.evaluateNumeric() + this.right.evaluateNumeric(); }
    diff(varName) { return new Add(this.left.diff(varName), this.right.diff(varName)); }
    integrate(varName) { return new Add(this.left.integrate(varName), this.right.integrate(varName)); }
    expand() { return new Add(this.left.expand(), this.right.expand()); }
    toLatex() { return `${this.left.toLatex()} + ${this.right.toLatex()}`; }
}

class Sub extends BinaryOp {
    toString() { return `(${this.left} - ${this.right})`; }
    simplify() {
        const l = this.left.simplify();
        const r = this.right.simplify();

        // Vector subtraction
        if (l instanceof Vec && r instanceof Vec) {
            if (l.elements.length !== r.elements.length) throw new Error("Vector length mismatch in subtraction");
            const newElements = [];
            for (let i = 0; i < l.elements.length; i++) {
                newElements.push(new Sub(l.elements[i], r.elements[i]).simplify());
            }
            return new Vec(newElements);
        }

        if (l instanceof Num && r instanceof Num) return new Num(l.value - r.value);
        if (r instanceof Num && r.value === 0) return l;
        if (l.toString() === r.toString()) return new Num(0);
        return new Sub(l, r);
    }
    evaluateNumeric() { return this.left.evaluateNumeric() - this.right.evaluateNumeric(); }
    diff(varName) { return new Sub(this.left.diff(varName), this.right.diff(varName)); }
    integrate(varName) { return new Sub(this.left.integrate(varName), this.right.integrate(varName)); }
    expand() { return new Sub(this.left.expand(), this.right.expand()); }
    toLatex() { return `${this.left.toLatex()} - ${this.right.toLatex()}`; }
}

class Mul extends BinaryOp {
    toString() { return `(${this.left} * ${this.right})`; }
    simplify() {
        const l = this.left.simplify();
        const r = this.right.simplify();

        // Complex number rule: i * i = -1
        if (l instanceof Symbol && l.name === 'i' && r instanceof Symbol && r.name === 'i') {
            return new Num(-1);
        }

        // Matrix multiplication logic
        if (l instanceof Vec && r instanceof Vec) {
            // Check if Matrix * Matrix or Matrix * Vector or Vector * Vector (Dot)
            const lIsMatrix = l.elements.length > 0 && l.elements[0] instanceof Vec;
            const rIsMatrix = r.elements.length > 0 && r.elements[0] instanceof Vec;

            if (!lIsMatrix && !rIsMatrix) {
                if (l.elements.length !== r.elements.length) throw new Error("Vector length mismatch in dot product");
                let sum = new Num(0);
                for (let i = 0; i < l.elements.length; i++) {
                    sum = new Add(sum, new Mul(l.elements[i], r.elements[i])).simplify();
                }
                return sum;
            } else if (!lIsMatrix && rIsMatrix) {
                // Vector * Matrix (Row Vector * Matrix)
                const colsVector = l.elements.length;
                const rowsMatrix = r.elements.length;
                const colsMatrix = r.elements[0].elements.length;

                if (colsVector !== rowsMatrix) throw new Error(`Vector-Matrix dimension mismatch: ${colsVector} * ${rowsMatrix}x${colsMatrix}`);

                const newElements = [];
                for (let j = 0; j < colsMatrix; j++) {
                    let sum = new Num(0);
                    for (let k = 0; k < colsVector; k++) {
                        sum = new Add(sum, new Mul(l.elements[k], r.elements[k].elements[j])).simplify();
                    }
                    newElements.push(sum);
                }
                return new Vec(newElements);
            } else if (lIsMatrix && rIsMatrix) {
                const rowsA = l.elements.length;
                const colsA = l.elements[0].elements.length;
                const rowsB = r.elements.length;
                const colsB = r.elements[0].elements.length;

                if (colsA !== rowsB) throw new Error(`Matrix dimension mismatch: ${rowsA}x${colsA} * ${rowsB}x${colsB}`);

                const resultRows = [];
                for (let i = 0; i < rowsA; i++) {
                    const row = [];
                    for (let j = 0; j < colsB; j++) {
                        let sum = new Num(0);
                        for (let k = 0; k < colsA; k++) {
                            sum = new Add(sum, new Mul(l.elements[i].elements[k], r.elements[k].elements[j])).simplify();
                        }
                        row.push(sum);
                    }
                    resultRows.push(new Vec(row));
                }
                return new Vec(resultRows);
            } else if (lIsMatrix && !rIsMatrix) {
                // Matrix * Vector
                const rowsA = l.elements.length;
                const colsA = l.elements[0].elements.length;
                const rowsB = r.elements.length;

                if (colsA !== rowsB) throw new Error(`Matrix-Vector dimension mismatch: ${rowsA}x${colsA} * ${rowsB}`);

                const newElements = [];
                for (let i = 0; i < rowsA; i++) {
                    let sum = new Num(0);
                    for (let k = 0; k < colsA; k++) {
                        sum = new Add(sum, new Mul(l.elements[i].elements[k], r.elements[k])).simplify();
                    }
                    newElements.push(sum);
                }
                return new Vec(newElements);
            }
        }

        // Scalar * Vector
        if (l instanceof Num && r instanceof Vec) {
            return new Vec(r.elements.map(e => new Mul(l, e).simplify()));
        }
        if (l instanceof Vec && r instanceof Num) {
            return new Vec(l.elements.map(e => new Mul(e, r).simplify()));
        }

        if (l instanceof Num && r instanceof Num) return new Num(l.value * r.value);
        if (l instanceof Num && l.value === 0) return new Num(0);
        if (r instanceof Num && r.value === 0) return new Num(0);
        if (l instanceof Num && l.value === 1) return r;
        if (r instanceof Num && r.value === 1) return l;

        // Simplify Mul with Div: a * (b / c) -> (a * b) / c
        if (r instanceof Div) {
            return new Div(new Mul(l, r.left), r.right).simplify();
        }
        if (l instanceof Div) {
            return new Div(new Mul(l.left, r), l.right).simplify();
        }

        return new Mul(l, r);
    }
    evaluateNumeric() { return this.left.evaluateNumeric() * this.right.evaluateNumeric(); }
    diff(varName) {
        return new Add(new Mul(this.left.diff(varName), this.right), new Mul(this.left, this.right.diff(varName)));
    }
    integrate(varName) {
        if (this.left instanceof Num) return new Mul(this.left, this.right.integrate(varName));
        if (this.right instanceof Num) return new Mul(this.right, this.left.integrate(varName));
        return new Call("integrate", [this, varName]);
    }
    expand() {
        const l = this.left.expand();
        const r = this.right.expand();
        if (l instanceof Add) return new Add(new Mul(l.left, r).expand(), new Mul(l.right, r).expand());
        if (r instanceof Add) return new Add(new Mul(l, r.left).expand(), new Mul(l, r.right).expand());
        return new Mul(l, r);
    }
    toLatex() {
        let lTex = this.left.toLatex();
        let rTex = this.right.toLatex();
        let op = " \\cdot ";

        // Implicit multiplication formatting rules
        if (this.left instanceof Num && this.right instanceof Symbol) op = ""; // 2x
        if (this.left instanceof Num && this.right instanceof Call) op = ""; // 2sin(x)
        if (this.left instanceof Symbol && this.right instanceof Symbol) op = " "; // x y
        if (this.left instanceof Symbol && this.right instanceof Call) op = ""; // x sin(y)
        if (this.left instanceof Num && this.right instanceof Pow && this.right.left instanceof Symbol) op = ""; // 2x^2

        if (this.left instanceof Add || this.left instanceof Sub) lTex = `\\left(${lTex}\\right)`;
        if (this.right instanceof Add || this.right instanceof Sub) rTex = `\\left(${rTex}\\right)`;
        return `${lTex}${op}${rTex}`;
    }
}

class Div extends BinaryOp {
    toString() { return `(${this.left} / ${this.right})`; }
    simplify() {
        const l = this.left.simplify();
        const r = this.right.simplify();
        if (l instanceof Num && r instanceof Num) {
            if (r.value === 0) throw new Error("Division by zero");
            if (l.value % r.value === 0) return new Num(l.value / r.value);
            // Simplify signs: a/-b -> -a/b, -a/-b -> a/b
            if (r.value < 0) {
                return new Div(new Num(-l.value), new Num(-r.value)).simplify();
            }
            return new Div(l, r);
        }
        // Vector division by scalar: Vec / Num
        if (l instanceof Vec && r instanceof Num) {
            return new Vec(l.elements.map(e => new Div(e, r).simplify()));
        }

        if (l instanceof Num && l.value === 0) return new Num(0);
        if (r instanceof Num && r.value === 1) return l;
        if (l.toString() === r.toString()) return new Num(1);

        // Cancellation: (a * b) / a -> b
        if (l instanceof Mul) {
            if (l.left.toString() === r.toString()) return l.right;
            if (l.right.toString() === r.toString()) return l.left;
        }
        // Cancellation: a / (a * b) -> 1/b (Not full implementation but basic)

        return new Div(l, r);
    }
    evaluateNumeric() { return this.left.evaluateNumeric() / this.right.evaluateNumeric(); }
    diff(varName) {
        const num = new Sub(new Mul(this.left.diff(varName), this.right), new Mul(this.left, this.right.diff(varName)));
        const den = new Pow(this.right, new Num(2));
        return new Div(num, den);
    }
    integrate(varName) {
        if (this.left instanceof Num && this.left.value === 1 && this.right instanceof Symbol && this.right.name === varName.name) {
            return new Call("ln", [varName]);
        }
        if (this.right instanceof Num) {
            return new Mul(new Div(new Num(1), this.right), this.left.integrate(varName));
        }
        return new Call("integrate", [this, varName]);
    }
    expand() {
        const l = this.left.expand();
        const r = this.right.expand();
        if (l instanceof Add) return new Add(new Div(l.left, r).expand(), new Div(l.right, r).expand());
        return new Div(l, r);
    }
    toLatex() { return `\\frac{${this.left.toLatex()}}{${this.right.toLatex()}}`; }
}

class Pow extends BinaryOp {
    toString() { return `(${this.left}^${this.right})`; }
    simplify() {
        const l = this.left.simplify();
        const r = this.right.simplify();

        // i^2 = -1
        if (l instanceof Symbol && l.name === 'i' && r instanceof Num && r.value === 2) {
            return new Num(-1);
        }

        if (r instanceof Num && r.value === 0) return new Num(1);
        if (r instanceof Num && r.value === 1) return l;
        if (l instanceof Num && r instanceof Num) return new Num(Math.pow(l.value, r.value));
        return new Pow(l, r);
    }
    evaluateNumeric() { return Math.pow(this.left.evaluateNumeric(), this.right.evaluateNumeric()); }
    diff(varName) {
        if (this.right instanceof Num) {
            const n = this.right;
            const u = this.left;
            return new Mul(new Mul(n, new Pow(u, new Num(n.value - 1))), u.diff(varName));
        }
        return new Mul(this, new Mul(this.right, new Call("ln", [this.left])).diff(varName));
    }
    integrate(varName) {
        if (this.left instanceof Symbol && this.left.name === varName.name) {
            if (this.right instanceof Num) {
                if (this.right.value === -1) return new Call("ln", [varName]);
                const n = this.right.value;
                return new Div(new Pow(this.left, new Num(n + 1)), new Num(n + 1));
            }
        }
        return new Call("integrate", [this, varName]);
    }
    expand() {
        const l = this.left.expand();
        const r = this.right.expand();
        if (r instanceof Num && l instanceof Add && r.value === 2) {
            const a = l.left;
            const b = l.right;
            return new Add(new Add(new Pow(a, new Num(2)), new Mul(new Num(2), new Mul(a, b))), new Pow(b, new Num(2))).expand();
        }
        return new Pow(l, r);
    }
    toLatex() {
        let lTex = this.left.toLatex();
        if (this.left instanceof Add || this.left instanceof Sub || this.left instanceof Mul || this.left instanceof Div) lTex = `\\left(${lTex}\\right)`;
        return `{${lTex}}^{${this.right.toLatex()}}`;
    }
}

class Call extends Expr {
    constructor(funcName, args) {
        super();
        this.funcName = funcName;
        this.args = args;
    }
    toString() { return `${this.funcName}(${this.args.join(", ")})`; }
    simplify() {
        const simpleArgs = this.args.map(a => a.simplify());

        if (this.funcName === 'sqrt') {
            const arg = simpleArgs[0];
            if (arg instanceof Num) {
                if (arg.value >= 0) {
                    const sqrtVal = Math.sqrt(arg.value);
                    if (Number.isInteger(sqrtVal)) return new Num(sqrtVal);
                    // Return symbolic sqrt if not integer?
                    // But we must avoid infinite loop if we return Call('sqrt') and simplify calls this.
                    // We can just return the Call object (this) if we haven't changed args,
                    // or a new Call with simplified args.
                    // Since simpleArgs are simplified, we just return new Call.
                    return new Call('sqrt', simpleArgs);
                }
                // sqrt(-x)
                const absVal = Math.abs(arg.value);
                const sqrtVal = Math.sqrt(absVal);
                if (Number.isInteger(sqrtVal)) {
                    return new Mul(new Symbol('i'), new Num(sqrtVal));
                }
                return new Mul(new Symbol('i'), new Call('sqrt', [new Num(absVal)]));
            }
            // sqrt(x) -> x^0.5
            // Maybe keep as sqrt(x) for display?
            // Xcas uses sqrt(x). Pow(0.5) is for calculus often.
            // Let's keep sqrt(x) unless we want to normalize.
            return new Call('sqrt', simpleArgs);
            // return new Pow(arg, new Num(0.5));
        }

        if (this.funcName === 'sin') {
            const arg = simpleArgs[0];
            if (arg instanceof Num && arg.value === 0) return new Num(0);
            if (arg instanceof Num && arg.value === Math.PI) return new Num(0);
        }
        if (this.funcName === 'cos') {
            const arg = simpleArgs[0];
            if (arg instanceof Num && arg.value === 0) return new Num(1);
        }
        if (this.funcName === 'tan') {
            const arg = simpleArgs[0];
            if (arg instanceof Num && arg.value === 0) return new Num(0);
        }
        if (this.funcName === 'asin') {
            const arg = simpleArgs[0];
            if (arg instanceof Num && arg.value === 0) return new Num(0);
            if (arg instanceof Num && arg.value === 1) return new Num(Math.PI / 2); // Approximate? Num stores float.
            if (arg instanceof Num && arg.value === -1) return new Num(-Math.PI / 2);
        }
        if (this.funcName === 'acos') {
            const arg = simpleArgs[0];
            if (arg instanceof Num && arg.value === 1) return new Num(0);
            if (arg instanceof Num && arg.value === 0) return new Num(Math.PI / 2);
            if (arg instanceof Num && arg.value === -1) return new Num(Math.PI);
        }
        if (this.funcName === 'atan') {
            const arg = simpleArgs[0];
            if (arg instanceof Num && arg.value === 0) return new Num(0);
        }
        if (this.funcName === 'sinh') {
            const arg = simpleArgs[0];
            if (arg instanceof Num && arg.value === 0) return new Num(0);
        }
        if (this.funcName === 'cosh') {
            const arg = simpleArgs[0];
            if (arg instanceof Num && arg.value === 0) return new Num(1);
        }
        if (this.funcName === 'tanh') {
            const arg = simpleArgs[0];
            if (arg instanceof Num && arg.value === 0) return new Num(0);
        }
        if (this.funcName === 'ln') {
            const arg = simpleArgs[0];
            if (arg instanceof Num && arg.value === 1) return new Num(0);
            if (arg instanceof Symbol && arg.name === 'e') return new Num(1);
        }
        if (this.funcName === 'log') {
            const arg = simpleArgs[0];
            if (arg instanceof Num && arg.value === 1) return new Num(0);
            if (arg instanceof Num && arg.value === 10) return new Num(1);
            // log(10^x) -> x
            if (arg instanceof Pow && arg.left instanceof Num && arg.left.value === 10) return arg.right;
            // log(100) -> 2
            if (arg instanceof Num) {
                const val = Math.log10(arg.value);
                if (Number.isInteger(val)) return new Num(val);
            }
        }
        if (this.funcName === 'exp') {
            const arg = simpleArgs[0];
            if (arg instanceof Num && arg.value === 0) return new Num(1);
        }

        return new Call(this.funcName, simpleArgs);
    }
    evaluateNumeric() {
        const argsVal = this.args.map(a => a.evaluateNumeric());
        if (this.funcName === 'sin') return Math.sin(argsVal[0]);
        if (this.funcName === 'cos') return Math.cos(argsVal[0]);
        if (this.funcName === 'tan') return Math.tan(argsVal[0]);
        if (this.funcName === 'asin') return Math.asin(argsVal[0]);
        if (this.funcName === 'acos') return Math.acos(argsVal[0]);
        if (this.funcName === 'atan') return Math.atan(argsVal[0]);
        if (this.funcName === 'sinh') return Math.sinh(argsVal[0]);
        if (this.funcName === 'cosh') return Math.cosh(argsVal[0]);
        if (this.funcName === 'tanh') return Math.tanh(argsVal[0]);
        if (this.funcName === 'exp') return Math.exp(argsVal[0]);
        if (this.funcName === 'ln') return Math.log(argsVal[0]);
        if (this.funcName === 'log') return Math.log10(argsVal[0]);
        if (this.funcName === 'sqrt') return Math.sqrt(argsVal[0]);
        if (this.funcName === 'abs') return Math.abs(argsVal[0]);
        return NaN; // Unknown
    }
    diff(varName) {
        const u = this.args[0];
        if (this.funcName === 'sin') return new Mul(new Call('cos', [u]), u.diff(varName));
        if (this.funcName === 'cos') return new Mul(new Mul(new Num(-1), new Call('sin', [u])), u.diff(varName));
        if (this.funcName === 'tan') return new Mul(new Div(new Num(1), new Pow(new Call('cos', [u]), new Num(2))), u.diff(varName));

        if (this.funcName === 'asin') {
            // 1/sqrt(1-u^2) * u'
            return new Mul(new Div(new Num(1), new Call('sqrt', [new Sub(new Num(1), new Pow(u, new Num(2)))])), u.diff(varName));
        }
        if (this.funcName === 'acos') {
            // -1/sqrt(1-u^2) * u'
            return new Mul(new Div(new Num(-1), new Call('sqrt', [new Sub(new Num(1), new Pow(u, new Num(2)))])), u.diff(varName));
        }
        if (this.funcName === 'atan') {
            // 1/(1+u^2) * u'
            return new Mul(new Div(new Num(1), new Add(new Num(1), new Pow(u, new Num(2)))), u.diff(varName));
        }
        if (this.funcName === 'sinh') return new Mul(new Call('cosh', [u]), u.diff(varName));
        if (this.funcName === 'cosh') return new Mul(new Call('sinh', [u]), u.diff(varName));
        if (this.funcName === 'tanh') return new Mul(new Div(new Num(1), new Pow(new Call('cosh', [u]), new Num(2))), u.diff(varName));

        if (this.funcName === 'exp') return new Mul(this, u.diff(varName));
        if (this.funcName === 'ln') return new Div(u.diff(varName), u);
        if (this.funcName === 'log') {
            // d/dx log10(u) = u' / (u * ln(10))
            return new Div(u.diff(varName), new Mul(u, new Call('ln', [new Num(10)])));
        }
        if (this.funcName === 'sqrt') return new Div(u.diff(varName), new Mul(new Num(2), new Call('sqrt', [u])));
        // Default to symbolic diff
        return new Call('diff', [this, varName]);
    }
    integrate(varName) {
        if (this.funcName === 'sin' && this.args[0].toString() === varName.toString()) return new Mul(new Num(-1), new Call('cos', [varName]));
        if (this.funcName === 'cos' && this.args[0].toString() === varName.toString()) return new Call('sin', [varName]);
        if (this.funcName === 'exp' && this.args[0].toString() === varName.toString()) return this;
        return new Call("integrate", [this, varName]);
    }
    substitute(varName, value) {
        return new Call(this.funcName, this.args.map(a => a.substitute(varName, value)));
    }
    toLatex() {
        const argsTex = this.args.map(a => a.toLatex());

        if (this.funcName === 'sqrt') return `\\sqrt{${argsTex[0]}}`;

        const standardFunctions = ['sin', 'cos', 'tan', 'sinh', 'cosh', 'tanh', 'exp', 'ln', 'log', 'det', 'gcd', 'sec', 'csc', 'cot'];
        if (standardFunctions.includes(this.funcName)) {
            return `\\${this.funcName}\\left(${argsTex.join(", ")}\\right)`;
        }

        const mapFunctions = {
            'asin': '\\arcsin',
            'acos': '\\arccos',
            'atan': '\\arctan',
            'lcm': '\\operatorname{lcm}'
        };
        if (mapFunctions.hasOwnProperty(this.funcName)) {
            return `${mapFunctions[this.funcName]}\\left(${argsTex.join(", ")}\\right)`;
        }

        if (this.funcName === 'inv' && argsTex.length === 1) {
             return `\\left(${argsTex[0]}\\right)^{-1}`;
        }
        if (this.funcName === 'trans' && argsTex.length === 1) {
             return `\\left(${argsTex[0]}\\right)^{T}`;
        }
        if (this.funcName === 'cross' && argsTex.length === 2) {
             return `\\left(${argsTex[0]} \\times ${argsTex[1]}\\right)`;
        }
        if (this.funcName === 'abs' && argsTex.length === 1) {
             return `\\left|${argsTex[0]}\\right|`;
        }

        if (this.funcName === 'limit') {
            // limit(expr, var, val)
            if (argsTex.length === 3) {
                return `\\lim_{${argsTex[1]} \\to ${argsTex[2]}} ${argsTex[0]}`;
            }
        }

        if (this.funcName === 'integrate') {
            // integrate(expr, var) or integrate(expr, var, a, b)
            if (argsTex.length === 2) {
                return `\\int ${argsTex[0]} \\, d${argsTex[1]}`;
            }
            if (argsTex.length === 4) {
                return `\\int_{${argsTex[2]}}^{${argsTex[3]}} ${argsTex[0]} \\, d${argsTex[1]}`;
            }
        }

        if (this.funcName === 'sum') {
            // sum(expr, var, start, end)
            if (argsTex.length === 4) {
                return `\\sum_{${argsTex[1]}=${argsTex[2]}}^{${argsTex[3]}} ${argsTex[0]}`;
            }
        }

        if (this.funcName === 'product') {
            // product(expr, var, start, end)
            if (argsTex.length === 4) {
                return `\\prod_{${argsTex[1]}=${argsTex[2]}}^{${argsTex[3]}} ${argsTex[0]}`;
            }
        }

        if (this.funcName === 'diff') {
            // diff(expr, var)
            if (argsTex.length === 2) {
                // If expr is complex, maybe \frac{d}{dx} (expr)
                return `\\frac{d}{d ${argsTex[1]}} \\left( ${argsTex[0]} \\right)`;
            }
        }

        if (this.funcName === 'factored') {
            return argsTex.join(" \\cdot ");
        }

        return `\\text{${this.funcName}}\\left(${argsTex.join(", ")}\\right)`;
    }
}

class Assignment extends Expr {
    constructor(target, value) {
        super();
        this.target = target;
        this.value = value;
    }
    toString() { return `${this.target} := ${this.value}`; }
    simplify() { return new Assignment(this.target, this.value.simplify()); }
    evaluateNumeric() { return this.value.evaluateNumeric(); }
    diff(varName) { throw new Error("Cannot differentiate assignment"); }
    integrate(varName) { throw new Error("Cannot integrate assignment"); }
    expand() { return new Assignment(this.target, this.value.expand()); }
    substitute(varName, value) { return new Assignment(this.target, this.value.substitute(varName, value)); }
    toLatex() { return `${this.target.toLatex()} \\coloneqq ${this.value.toLatex()}`; }
}

class Eq extends Expr {
    constructor(left, right) {
        super();
        this.left = left;
        this.right = right;
    }
    toString() { return `${this.left} = ${this.right}`; }
    simplify() { return new Eq(this.left.simplify(), this.right.simplify()); }
    evaluateNumeric() { return NaN; }
    diff(varName) { return new Eq(this.left.diff(varName), this.right.diff(varName)); }
    integrate(varName) { return new Eq(this.left.integrate(varName), this.right.integrate(varName)); }
    expand() { return new Eq(this.left.expand(), this.right.expand()); }
    substitute(varName, value) { return new Eq(this.left.substitute(varName, value), this.right.substitute(varName, value)); }
    toLatex() { return `${this.left.toLatex()} = ${this.right.toLatex()}`; }
}

class Vec extends Expr {
    constructor(elements) {
        super();
        this.elements = elements;
    }
    toString() { return `[${this.elements.join(", ")}]`; }
    simplify() { return new Vec(this.elements.map(e => e.simplify())); }
    evaluateNumeric() { return NaN; } // Or array of numbers?
    diff(varName) { return new Vec(this.elements.map(e => e.diff(varName))); }
    integrate(varName) { return new Vec(this.elements.map(e => e.integrate(varName))); }
    substitute(varName, value) { return new Vec(this.elements.map(e => e.substitute(varName, value))); }
    toLatex() {
        // Heuristic: check if elements are vectors (Matrix)
        if (this.elements.length > 0 && this.elements[0] instanceof Vec) {
            // Matrix
            const rows = this.elements.map(row => row.elements.map(e => e.toLatex()).join(" & ")).join(" \\\\ ");
            return `\\begin{bmatrix} ${rows} \\end{bmatrix}`;
        }
        return `\\begin{bmatrix} ${this.elements.map(e => e.toLatex()).join(" \\\\ ")} \\end{bmatrix}`; // Column vector default? Or row?
    }
}

class FunctionDef extends Expr {
    constructor(name, params, body) {
        super();
        this.name = name;
        this.params = params;
        this.body = body;
    }
    toString() { return `${this.name}(${this.params.join(", ")}) := ${this.body}`; }
    simplify() { return new FunctionDef(this.name, this.params, this.body.simplify()); }
    evaluateNumeric() { return NaN; }
    toLatex() {
        const paramsTex = this.params.join(", ");
        return `\\text{${this.name}}(${paramsTex}) \\coloneqq ${this.body.toLatex()}`;
    }
}
