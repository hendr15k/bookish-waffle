
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
    dependsOn(varName) {
        if (this instanceof Sym) return this.name === varName.name;
        if (this instanceof Num) return false;
        if (this instanceof BinaryOp) return this.left.dependsOn(varName) || this.right.dependsOn(varName);
        if (this instanceof Call) return this.args.some(a => a.dependsOn(varName));
        if (this instanceof Vec) return this.elements.some(e => e.dependsOn(varName));
        if (this instanceof Pow) return this.left.dependsOn(varName) || this.right.dependsOn(varName);
        if (this instanceof Assignment) return this.value.dependsOn(varName);
        if (this instanceof Eq) return this.left.dependsOn(varName) || this.right.dependsOn(varName);
        if (this instanceof Not) return this.arg.dependsOn(varName);
        if (this instanceof At) return this.obj.dependsOn(varName) || this.index.dependsOn(varName);
        return false;
    }
    expand() {
        // Expand properties like log(a*b) -> log(a) + log(b)
        // Only if user explicitly calls expand()
        if (this.args) {
            const simpleArgs = this.args.map(a => a.expand());

            if (this.funcName === 'log' || this.funcName === 'ln') {
                const arg = simpleArgs[0];
                if (arg instanceof Mul) {
                    // log(a*b) = log(a) + log(b)
                    return new Add(new Call(this.funcName, [arg.left]), new Call(this.funcName, [arg.right])).expand();
                }
                if (arg instanceof Div) {
                    // log(a/b) = log(a) - log(b)
                    return new Sub(new Call(this.funcName, [arg.left]), new Call(this.funcName, [arg.right])).expand();
                }
                if (arg instanceof Pow) {
                    // log(a^b) = b * log(a)
                    return new Mul(arg.right, new Call(this.funcName, [arg.left])).expand();
                }
            }

            // Trigonometric Expansion
            if (this.funcName === 'sin') {
                const arg = simpleArgs[0];
                // sin(2*x) -> 2*sin(x)*cos(x)
                if (arg instanceof Mul) {
                    if (arg.left instanceof Num && arg.left.value === 2) {
                        return new Mul(new Num(2), new Mul(new Call('sin', [arg.right]), new Call('cos', [arg.right]))).expand();
                    }
                    if (arg.right instanceof Num && arg.right.value === 2) {
                        return new Mul(new Num(2), new Mul(new Call('sin', [arg.left]), new Call('cos', [arg.left]))).expand();
                    }
                }
                // sin(a+b) -> sin(a)cos(b) + cos(a)sin(b)
                if (arg instanceof Add) {
                    const a = arg.left;
                    const b = arg.right;
                    return new Add(
                        new Mul(new Call('sin', [a]), new Call('cos', [b])),
                        new Mul(new Call('cos', [a]), new Call('sin', [b]))
                    ).expand();
                }
                // sin(a-b) -> sin(a)cos(b) - cos(a)sin(b)
                if (arg instanceof Sub) {
                    const a = arg.left;
                    const b = arg.right;
                    return new Sub(
                        new Mul(new Call('sin', [a]), new Call('cos', [b])),
                        new Mul(new Call('cos', [a]), new Call('sin', [b]))
                    ).expand();
                }
            }

            if (this.funcName === 'cos') {
                const arg = simpleArgs[0];
                // cos(2*x) -> cos(x)^2 - sin(x)^2
                if (arg instanceof Mul) {
                    if (arg.left instanceof Num && arg.left.value === 2) {
                        return new Sub(new Pow(new Call('cos', [arg.right]), new Num(2)), new Pow(new Call('sin', [arg.right]), new Num(2))).expand();
                    }
                    if (arg.right instanceof Num && arg.right.value === 2) {
                        return new Sub(new Pow(new Call('cos', [arg.left]), new Num(2)), new Pow(new Call('sin', [arg.left]), new Num(2))).expand();
                    }
                }
                // cos(a+b) -> cos(a)cos(b) - sin(a)sin(b)
                if (arg instanceof Add) {
                    const a = arg.left;
                    const b = arg.right;
                    return new Sub(
                        new Mul(new Call('cos', [a]), new Call('cos', [b])),
                        new Mul(new Call('sin', [a]), new Call('sin', [b]))
                    ).expand();
                }
                // cos(a-b) -> cos(a)cos(b) + sin(a)sin(b)
                if (arg instanceof Sub) {
                    const a = arg.left;
                    const b = arg.right;
                    return new Add(
                        new Mul(new Call('cos', [a]), new Call('cos', [b])),
                        new Mul(new Call('sin', [a]), new Call('sin', [b]))
                    ).expand();
                }
            }

            if (this.funcName === 'tan') {
                const arg = simpleArgs[0];
                // tan(2*x) -> 2tan(x)/(1-tan(x)^2)
                if (arg instanceof Mul) {
                    if ((arg.left instanceof Num && arg.left.value === 2) || (arg.right instanceof Num && arg.right.value === 2)) {
                        const x = (arg.left instanceof Num) ? arg.right : arg.left;
                        return new Div(
                            new Mul(new Num(2), new Call('tan', [x])),
                            new Sub(new Num(1), new Pow(new Call('tan', [x]), new Num(2)))
                        ).expand();
                    }
                }
                // tan(a+b) -> (tan(a)+tan(b))/(1-tan(a)tan(b))
                if (arg instanceof Add) {
                    const a = arg.left;
                    const b = arg.right;
                    return new Div(
                        new Add(new Call('tan', [a]), new Call('tan', [b])),
                        new Sub(new Num(1), new Mul(new Call('tan', [a]), new Call('tan', [b])))
                    ).expand();
                }
                // tan(a-b) -> (tan(a)-tan(b))/(1+tan(a)tan(b))
                if (arg instanceof Sub) {
                    const a = arg.left;
                    const b = arg.right;
                    return new Div(
                        new Sub(new Call('tan', [a]), new Call('tan', [b])),
                        new Add(new Num(1), new Mul(new Call('tan', [a]), new Call('tan', [b])))
                    ).expand();
                }
            }

            return new Call(this.funcName, simpleArgs);
        }
        return this;
    }
    substitute(varName, value) { return this; }
    toLatex() { return this.toString(); }
}

function toExpr(other) {
    if (other instanceof Expr) return other;
    if (typeof other === 'number') return new Num(other);
    throw new Error("Cannot convert to Expr: " + other);
}

function getPolyCoeffs(expr, varNode) {
    if (!varNode) return null;
    const x = varNode;
    const getTermInfo = (term) => {
        if (term instanceof Num) return { deg: 0, coeff: term.value };
        if (term instanceof Sym && term.name === x.name) return { deg: 1, coeff: 1 };
        if (term instanceof Mul && term.left instanceof Num && term.right instanceof Sym && term.right.name === x.name) return { deg: 1, coeff: term.left.value };
        if (term instanceof Mul && term.right instanceof Num && term.left instanceof Sym && term.left.name === x.name) return { deg: 1, coeff: term.right.value };
        if (term instanceof Pow && term.left instanceof Sym && term.left.name === x.name && term.right instanceof Num && Number.isInteger(term.right.value) && term.right.value >= 0) return { deg: term.right.value, coeff: 1 };
        if (term instanceof Mul && term.left instanceof Num && term.right instanceof Pow && term.right.left instanceof Sym && term.right.left.name === x.name && term.right.right instanceof Num && Number.isInteger(term.right.right.value)) {
            return { deg: term.right.right.value, coeff: term.left.value };
        }
        if (term instanceof Mul && term.left instanceof Num && term.left.value === -1 && term.right instanceof Sym && term.right.name === x.name) return { deg: 1, coeff: -1 };
        return null;
    };

    const terms = [];
    const collectTerms = (e, sign = 1) => {
        if (e instanceof Add) { collectTerms(e.left, sign); collectTerms(e.right, sign); }
        else if (e instanceof Sub) { collectTerms(e.left, sign); collectTerms(e.right, -sign); }
        else terms.push({ expr: e, sign });
    };
    collectTerms(expr);

    const coeffs = {};
    let maxDeg = 0;
    for(const item of terms) {
        const info = getTermInfo(item.expr);
        if (info) {
            const deg = info.deg;
            const val = info.coeff * item.sign;
            coeffs[deg] = (coeffs[deg] || 0) + val;
            if (deg > maxDeg) maxDeg = deg;
        } else {
            return null; // Not a polynomial in x
        }
    }

    const P = [];
    for(let i=0; i<=maxDeg; i++) P[i] = coeffs[i] || 0;
    while(P.length > 0 && Math.abs(P[P.length-1]) < 1e-15) P.pop();
    return P;
}

function polyRem(P, Q) {
    if (Q.length === 0) throw new Error("Divide by zero polynomial");
    let R = [...P];
    const m = Q.length - 1;
    const lc = Q[m];
    while (R.length - 1 >= m) {
        const degR = R.length - 1;
        const factor = R[degR] / lc;
        const degDiff = degR - m;
        for(let i=0; i<=m; i++) {
            R[i + degDiff] -= factor * Q[i];
        }
        while(R.length > 0 && Math.abs(R[R.length-1]) < 1e-15) R.pop();
    }
    return R;
}

function polyQuotient(P, Q) {
    if (Q.length === 0) throw new Error("Divide by zero");
    let R = [...P];
    let Quot = new Array(Math.max(0, R.length - Q.length + 1)).fill(0);
    const m = Q.length - 1;
    const lc = Q[m];
    while(R.length - 1 >= m) {
        const degR = R.length - 1;
        const factor = R[degR] / lc;
        const degDiff = degR - m;
        Quot[degDiff] = factor;
        for(let i=0; i<=m; i++) {
            R[i + degDiff] -= factor * Q[i];
        }
        while(R.length > 0 && Math.abs(R[R.length-1]) < 1e-15) R.pop();
    }
    while(Quot.length > 0 && Math.abs(Quot[Quot.length-1]) < 1e-15) Quot.pop();
    return Quot;
}

function polyGcd(P, Q) {
    let a = P;
    let b = Q;
    while (b.length > 0 && !(b.length === 1 && Math.abs(b[0]) < 1e-15)) {
        const r = polyRem(a, b);
        a = b;
        b = r;
    }
    if (a.length > 0) {
        const lc = a[a.length-1];
        if (Math.abs(lc) > 1e-15) {
            for(let i=0; i<a.length; i++) a[i] /= lc;
        }
    }
    return a;
}

function polyFromCoeffs(coeffs, varNode) {
    if (coeffs.length === 0) return new Num(0);
    let res = new Num(0);
    for(let i=0; i<coeffs.length; i++) {
        const c = coeffs[i];
        if (Math.abs(c) < 1e-15) continue;
        const intVal = Math.round(c);
        const numNode = (Math.abs(c - intVal) < 1e-10) ? new Num(intVal) : new Num(c);

        let term;
        if (i === 0) term = numNode;
        else if (i === 1) term = (Math.abs(c-1)<1e-10) ? varNode : ((Math.abs(c+1)<1e-10) ? new Mul(new Num(-1), varNode) : new Mul(numNode, varNode));
        else {
            const pow = new Pow(varNode, new Num(i));
            if (Math.abs(c - 1) < 1e-10) term = pow;
            else if (Math.abs(c + 1) < 1e-10) term = new Mul(new Num(-1), pow);
            else term = new Mul(numNode, pow);
        }

        if (res instanceof Num && res.value === 0) res = term;
        else res = new Add(res, term);
    }
    return res;
}

function polyDiv(numerator, denominator) {
    // Basic synthetic division for P(x) / (x - c)
    // First, extract coefficients of numerator.
    // We assume numerator is a polynomial in one variable.
    // Denominator is (x - c).

    let x, c;
    if (denominator instanceof Sub) {
        x = denominator.left;
        c = denominator.right.value;
    } else if (denominator instanceof Add) {
        x = denominator.left;
        c = -denominator.right.value;
    } else {
        return new Div(numerator, denominator);
    }

    // Helper to get degree and coefficient of a term
    const getTermInfo = (term) => {
        if (term instanceof Num) return { deg: 0, coeff: term.value };
        if (term instanceof Sym && term.name === x.name) return { deg: 1, coeff: 1 };
        if (term instanceof Mul && term.left instanceof Num && term.right instanceof Sym && term.right.name === x.name) return { deg: 1, coeff: term.left.value };
        if (term instanceof Mul && term.right instanceof Num && term.left instanceof Sym && term.left.name === x.name) return { deg: 1, coeff: term.right.value }; // x*2
        if (term instanceof Pow && term.left instanceof Sym && term.left.name === x.name && term.right instanceof Num) return { deg: term.right.value, coeff: 1 };
        if (term instanceof Mul && term.left instanceof Num && term.right instanceof Pow && term.right.left instanceof Sym && term.right.left.name === x.name && term.right.right instanceof Num) {
            return { deg: term.right.right.value, coeff: term.left.value };
        }
        // Handle negative terms -x -> -1*x
        if (term instanceof Mul && term.left instanceof Num && term.left.value === -1 && term.right instanceof Sym && term.right.name === x.name) return { deg: 1, coeff: -1 };
        return null; // Not a simple monomial
    };

    // Flatten Add/Sub tree to list of terms
    const terms = [];
    const collectTerms = (expr, sign = 1) => {
        if (expr instanceof Add) {
            collectTerms(expr.left, sign);
            collectTerms(expr.right, sign);
        } else if (expr instanceof Sub) {
            collectTerms(expr.left, sign);
            collectTerms(expr.right, -sign);
        } else {
            terms.push({ expr, sign });
        }
    };
    collectTerms(numerator);

    const coeffs = {};
    let maxDeg = 0;

    for (const item of terms) {
        const info = getTermInfo(item.expr);
        if (info) {
             const deg = info.deg;
             const val = info.coeff * item.sign;
             coeffs[deg] = (coeffs[deg] || 0) + val;
             if (deg > maxDeg) maxDeg = deg;
        } else {
            // Found non-polynomial term, abort
            return new Div(numerator, denominator);
        }
    }

    // Convert to array [deg0, deg1, ...]
    const P = [];
    for(let i=0; i<=maxDeg; i++) P[i] = coeffs[i] || 0;

    // Synthetic division by c
    // P(x) / (x - c)
    // Coefficients are P[maxDeg], P[maxDeg-1], ... P[0]
    // We want Q(x) coeffs.
    // Q coeffs: q[n-1] = p[n]
    // q[k] = p[k+1] + c * q[k+1]
    // Remainder r = p[0] + c * q[0] (should be 0)

    const Q = new Array(maxDeg).fill(0);
    let carry = 0;
    // Correct loop
    for (let i = maxDeg - 1; i >= 0; i--) {
        const p_next = P[i+1];
        if (i === maxDeg - 1) {
            Q[i] = p_next;
        } else {
            Q[i] = P[i+1] + c * Q[i+1];
        }
    }

    // Construct Result Polynomial from Q
    let result = new Num(0);
    for (let i = 0; i < maxDeg; i++) {
        const coeff = Q[i];
        if (coeff === 0) continue;

        let term;
        if (i === 0) term = new Num(coeff);
        else if (i === 1) {
             if (coeff === 1) term = x;
             else if (coeff === -1) term = new Mul(new Num(-1), x);
             else term = new Mul(new Num(coeff), x);
        } else {
             const pow = new Pow(x, new Num(i));
             if (coeff === 1) term = pow;
             else if (coeff === -1) term = new Mul(new Num(-1), pow);
             else term = new Mul(new Num(coeff), pow);
        }

        result = new Add(result, term);
    }

    return result.simplify();
}

class Num extends Expr {
    constructor(value) {
        super();
        this.value = value;
    }
    toString() {
        if (Number.isInteger(this.value)) return this.value.toString();
        const precision = 12;
        return parseFloat(this.value.toPrecision(precision)).toString();
    }
    simplify() { return this; }
    evaluateNumeric() { return this.value; }
    diff(varName) { return new Num(0); }
    integrate(varName) {
        if (this.value === 0) return new Num(0);
        return new Mul(this, varName);
    }
    substitute(varName, value) { return this; }
    toLatex() {
        if (Number.isInteger(this.value)) return this.value.toString();
        const precision = 12;
        return parseFloat(this.value.toPrecision(precision)).toString();
    }
    equals(other) { return other instanceof Num && this.value === other.value; }
}

class Sym extends Expr {
    constructor(name) {
        super();
        this.name = name;
    }
    toString() { return this.name; }
    simplify() { return this; }
    evaluateNumeric() {
        if (this.name === 'pi') return Math.PI;
        if (this.name === 'e') return Math.E;
        if (this.name === 'Infinity' || this.name === 'infinity') return Infinity;
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
        if (this.name === varName.name) return value;
        // Check for deep equality if varName is complex (e.g. replacing sin(x) with u)
        if (varName instanceof Sym && this.name === varName.name) return value;
        // If varName is an expression, we can't substitute a symbol with it unless the symbol IS the expression (already checked)
        // But what if we want to substitute 'x' in 'x' with 'u'? (handled)

        // If we try to substitute a sub-expression, it is handled by the caller traversing the tree.
        // But if `varName` is `Pow(x, 2)` and `this` is `Sym('x')`, we don't match.
        return this;
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
    equals(other) { return other instanceof Sym && this.name === other.name; }
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

        // Associativity: (A + B) + C -> A + (B + C)
        // This helps combining terms like (1 - 2x) + 2x -> 1 + (-2x + 2x) -> 1
        if (l instanceof Add) {
            return new Add(l.left, new Add(l.right, r).simplify()).simplify();
        }

        // Vector addition
        if (l instanceof Vec && r instanceof Vec) {
            if (l.elements.length !== r.elements.length) throw new Error("Vector length mismatch in addition");
            const newElements = [];
            for (let i = 0; i < l.elements.length; i++) {
                newElements.push(new Add(l.elements[i], r.elements[i]).simplify());
            }
            return new Vec(newElements);
        }

        // Vector + Scalar
        if (l instanceof Vec && r instanceof Num) {
            const newElements = l.elements.map(e => new Add(e, r).simplify());
            return new Vec(newElements);
        }
        if (l instanceof Num && r instanceof Vec) {
            const newElements = r.elements.map(e => new Add(l, e).simplify());
            return new Vec(newElements);
        }

        if (l instanceof Num && r instanceof Num) {
            const sum = l.value + r.value;
            if (!isFinite(sum) && isFinite(l.value) && isFinite(r.value)) return new Sym("Infinity"); // Or handle signs
            // Handle Infinity
            if (l.value === Infinity && r.value === Infinity) return new Sym("Infinity");
            if (l.value === -Infinity && r.value === -Infinity) return new Mul(new Num(-1), new Sym("Infinity"));
            if (!isFinite(sum)) return new Sym("NaN"); // inf - inf
            return new Num(sum);
        }

        const isInf = (n) => n instanceof Sym && (n.name === 'Infinity' || n.name === 'infinity');
        const isNegInf = (n) => n instanceof Mul && n.left instanceof Num && n.left.value === -1 && isInf(n.right);
        const checkInf = (n) => isInf(n) ? 1 : (isNegInf(n) ? -1 : 0);

        const infL = checkInf(l);
        const infR = checkInf(r);

        if (infL !== 0 || infR !== 0) {
            if (infL === 1 && infR === 1) return new Sym("Infinity");
            if (infL === -1 && infR === -1) return new Mul(new Num(-1), new Sym("Infinity"));
            if (infL !== 0 && infR === 0) return l;
            if (infL === 0 && infR !== 0) return r;
            // inf - inf
            return new Sym("NaN");
        }

        // Symbolic Fraction addition
        // Combine if same denominator: A/C + B/C -> (A+B)/C
        const getDenom = (e) => (e instanceof Div) ? e.right : new Num(1);
        const getNum = (e) => (e instanceof Div) ? e.left : e;

        if (l instanceof Div || r instanceof Div) {
            const d1 = getDenom(l);
            const d2 = getDenom(r);

            // Check for identical denominators
            if (d1.toString() === d2.toString()) {
                const n1 = getNum(l);
                const n2 = getNum(r);
                return new Div(new Add(n1, n2).simplify(), d1).simplify();
            }

            // Combine if both are purely numeric fractions (existing logic)
            const isNumericDiv = (expr) => (expr instanceof Div && expr.left instanceof Num && expr.right instanceof Num);
            const isNumericNum = (expr) => (expr instanceof Num);

            if ((isNumericDiv(l) || isNumericNum(l)) && (isNumericDiv(r) || isNumericNum(r))) {
                 const n1Val = (l instanceof Div) ? l.left.value : l.value;
                 const d1Val = (l instanceof Div) ? l.right.value : 1;
                 const n2Val = (r instanceof Div) ? r.left.value : r.value;
                 const d2Val = (r instanceof Div) ? r.right.value : 1;

                 if (d1Val !== 0 && d2Val !== 0) {
                      const newNum = n1Val * d2Val + n2Val * d1Val;
                      const newDen = d1Val * d2Val;
                      return new Div(new Num(newNum), new Num(newDen)).simplify();
                 }
            }
        }

        if (l instanceof Num && l.value === 0) return r;
        if (r instanceof Num && r.value === 0) return l;
        if (r instanceof Num && r.value < 0) return new Sub(l, new Num(-r.value)).simplify();
        if (l instanceof Num && l.value < 0) return new Sub(r, new Num(-l.value)).simplify(); // (-a) + b -> b - a

        // x + (-y) -> x - y
        if (r instanceof Mul && r.left instanceof Num && r.left.value === -1) {
            return new Sub(l, r.right).simplify();
        }

        // (-x) + x -> 0
        if (l instanceof Mul && l.left instanceof Num && l.left.value === -1 && l.right.toString() === r.toString()) {
            return new Num(0);
        }

        if (l.toString() === r.toString()) return new Mul(new Num(2), l);

        // Logarithmic Combination: ln(a) + ln(b) -> ln(a*b)
        if (l instanceof Call && r instanceof Call && l.funcName === r.funcName && (l.funcName === 'ln' || l.funcName === 'log')) {
             return new Call(l.funcName, [new Mul(l.args[0], r.args[0])]).simplify();
        }

        // Pythagorean Identity: sin(x)^2 + cos(x)^2 -> 1
        if ((l instanceof Pow && r instanceof Pow) || (l instanceof Pow && r instanceof Mul) || (l instanceof Mul && r instanceof Pow)) {
             // Helper to check for sin(x)^2
             const isTrigSq = (expr, func) => {
                 return expr instanceof Pow &&
                        expr.left instanceof Call &&
                        expr.left.funcName === func &&
                        expr.right instanceof Num &&
                        expr.right.value === 2;
             };

             if ((isTrigSq(l, 'sin') && isTrigSq(r, 'cos')) ||
                 (isTrigSq(l, 'cos') && isTrigSq(r, 'sin'))) {
                 // Check if arguments match
                 if (l.left.args[0].toString() === r.left.args[0].toString()) {
                     return new Num(1);
                 }
             }
        }

        // (A - B) + B -> A
        if (l instanceof Sub) {
            if (l.right.toString() === r.toString()) return l.left;
        }
        // A + (B - A) -> B
        if (r instanceof Sub) {
            if (r.right.toString() === l.toString()) return r.left;
        }

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

        // Cancellation: (A+B)-A -> B, (A+B)-B -> A
        if (l instanceof Add) {
            if (l.left.toString() === r.toString()) return l.right;
            if (l.right.toString() === r.toString()) return l.left;
        }
        // (A-B)-A -> -B
        if (l instanceof Sub) {
            if (l.left.toString() === r.toString()) return new Mul(new Num(-1), l.right).simplify();
        }

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

        const isInf = (n) => n instanceof Sym && (n.name === 'Infinity' || n.name === 'infinity');
        const isNegInf = (n) => n instanceof Mul && n.left instanceof Num && n.left.value === -1 && isInf(n.right);
        const checkInf = (n) => isInf(n) ? 1 : (isNegInf(n) ? -1 : 0);

        const infL = checkInf(l);
        const infR = checkInf(r);

        if (infL !== 0 || infR !== 0) {
            if (infL === 1 && infR === -1) return new Sym("Infinity"); // Inf - (-Inf) = Inf
            if (infL === -1 && infR === 1) return new Mul(new Num(-1), new Sym("Infinity")); // -Inf - Inf = -Inf
            if (infL !== 0 && infR === 0) return l;
            if (infL === 0 && infR !== 0) return new Mul(new Num(-1), r).simplify();
            // Inf - Inf
            return new Sym("NaN");
        }

        if (l instanceof Num && l.value === 0) return new Mul(new Num(-1), r).simplify();

        // Symbolic Fraction subtraction
        // Combine if same denominator: A/C - B/C -> (A-B)/C
        const getDenom = (e) => (e instanceof Div) ? e.right : new Num(1);
        const getNum = (e) => (e instanceof Div) ? e.left : e;

        if (l instanceof Div || r instanceof Div) {
            const d1 = getDenom(l);
            const d2 = getDenom(r);

            // Check for identical denominators
            if (d1.toString() === d2.toString()) {
                const n1 = getNum(l);
                const n2 = getNum(r);
                return new Div(new Sub(n1, n2).simplify(), d1).simplify();
            }

            // Combine if both are purely numeric fractions
            const isNumericDiv = (expr) => (expr instanceof Div && expr.left instanceof Num && expr.right instanceof Num);
            const isNumericNum = (expr) => (expr instanceof Num);

            if ((isNumericDiv(l) || isNumericNum(l)) && (isNumericDiv(r) || isNumericNum(r))) {
                 const n1Val = (l instanceof Div) ? l.left.value : l.value;
                 const d1Val = (l instanceof Div) ? l.right.value : 1;
                 const n2Val = (r instanceof Div) ? r.left.value : r.value;
                 const d2Val = (r instanceof Div) ? r.right.value : 1;

                 if (d1Val !== 0 && d2Val !== 0) {
                      const newNum = n1Val * d2Val - n2Val * d1Val;
                      const newDen = d1Val * d2Val;
                      return new Div(new Num(newNum), new Num(newDen)).simplify();
                 }
            }
        }

        if (r instanceof Num && r.value === 0) return l;
        if (r instanceof Num && r.value < 0) return new Add(l, new Num(-r.value)).simplify();

        // x - (-y) -> x + y
        if (r instanceof Mul && r.left instanceof Num && r.left.value < 0) {
            const newR = new Mul(new Num(-r.left.value), r.right).simplify();
            return new Add(l, newR).simplify();
        }

        if (l.toString() === r.toString()) return new Num(0);

        // Logarithmic Combination: ln(a) - ln(b) -> ln(a/b)
        if (l instanceof Call && r instanceof Call && l.funcName === r.funcName && (l.funcName === 'ln' || l.funcName === 'log')) {
             return new Call(l.funcName, [new Div(l.args[0], r.args[0])]).simplify();
        }

        // Pythagorean Identities
        // 1 - sin(x)^2 -> cos(x)^2
        // 1 - cos(x)^2 -> sin(x)^2
        if (l instanceof Num && l.value === 1 && r instanceof Pow && r.right instanceof Num && r.right.value === 2) {
             if (r.left instanceof Call) {
                 if (r.left.funcName === 'sin') return new Pow(new Call('cos', r.left.args), new Num(2));
                 if (r.left.funcName === 'cos') return new Pow(new Call('sin', r.left.args), new Num(2));
             }
        }
        // sec(x)^2 - 1 -> tan(x)^2
        // csc(x)^2 - 1 -> cot(x)^2
        if (l instanceof Pow && l.right instanceof Num && l.right.value === 2 && r instanceof Num && r.value === 1) {
             if (l.left instanceof Call) {
                 if (l.left.funcName === 'sec') return new Pow(new Call('tan', l.left.args), new Num(2));
                 if (l.left.funcName === 'csc') return new Pow(new Call('cot', l.left.args), new Num(2));
             }
        }

        // sec^2(x) - tan^2(x) -> 1
        // csc^2(x) - cot^2(x) -> 1
        // tan^2(x) - sec^2(x) -> -1
        // cot^2(x) - csc^2(x) -> -1
        // cosh^2(x) - sinh^2(x) -> 1
        // sinh^2(x) - cosh^2(x) -> -1
        if (l instanceof Pow && r instanceof Pow && l.right instanceof Num && l.right.value === 2 && r.right instanceof Num && r.right.value === 2) {
             if (l.left instanceof Call && r.left instanceof Call && l.left.args[0].toString() === r.left.args[0].toString()) {
                 const f1 = l.left.funcName;
                 const f2 = r.left.funcName;
                 if (f1 === 'sec' && f2 === 'tan') return new Num(1);
                 if (f1 === 'csc' && f2 === 'cot') return new Num(1);
                 if (f1 === 'tan' && f2 === 'sec') return new Num(-1);
                 if (f1 === 'cot' && f2 === 'csc') return new Num(-1);
                 if (f1 === 'cosh' && f2 === 'sinh') return new Num(1);
                 if (f1 === 'sinh' && f2 === 'cosh') return new Num(-1);
             }
        }

        // (A + B) - A -> B
        if (l instanceof Add) {
            if (l.left.toString() === r.toString()) return l.right;
            if (l.right.toString() === r.toString()) return l.left;
        }
        // A - (A + B) -> -B
        if (r instanceof Add) {
             if (r.left.toString() === l.toString()) return new Mul(new Num(-1), r.right).simplify();
             if (r.right.toString() === l.toString()) return new Mul(new Num(-1), r.left).simplify();
        }
        // (A - B) - A -> -B
        if (l instanceof Sub) {
             if (l.left.toString() === r.toString()) return new Mul(new Num(-1), l.right).simplify();
        }
        // A - (A - B) -> B
        if (r instanceof Sub) {
             if (r.left.toString() === l.toString()) return r.right;
        }

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
        if (l instanceof Sym && l.name === 'i' && r instanceof Sym && r.name === 'i') {
            return new Num(-1);
        }

        // Complex identity: i * 1 -> i
        if (l instanceof Sym && l.name === 'i' && r instanceof Num && r.value === 1) return l;
        if (r instanceof Sym && r.name === 'i' && l instanceof Num && l.value === 1) return r;

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

        const isInf = (n) => n instanceof Sym && (n.name === 'Infinity' || n.name === 'infinity');
        const isNegInf = (n) => n instanceof Mul && n.left instanceof Num && n.left.value === -1 && isInf(n.right);
        const checkInf = (n) => isInf(n) ? 1 : (isNegInf(n) ? -1 : 0);

        const infL = checkInf(l);
        const infR = checkInf(r);

        if (infL !== 0 || infR !== 0) {
            // Check for 0 * Inf
            if ((l instanceof Num && l.value === 0) || (r instanceof Num && r.value === 0)) return new Sym("NaN");

            // Determine sign
            let sign = (infL !== 0 ? infL : 1) * (infR !== 0 ? infR : 1);
            if (l instanceof Num && l.value < 0) sign *= -1;
            if (r instanceof Num && r.value < 0) sign *= -1;

            if (sign > 0) return new Sym("Infinity");
            return new Mul(new Num(-1), new Sym("Infinity"));
        }

        if (l instanceof Num && l.value === 0) return new Num(0);
        if (r instanceof Num && r.value === 0) return new Num(0);
        if (l instanceof Num && l.value === 1) return r;
        if (r instanceof Num && r.value === 1) return l;

        // Scalar Associativity: c1 * (c2 * x) -> (c1*c2) * x
        if (l instanceof Num && r instanceof Mul && r.left instanceof Num) {
            return new Mul(new Num(l.value * r.left.value), r.right).simplify();
        }

        // Simplify Mul with Div: a * (b / c) -> (a * b) / c
        if (r instanceof Div) {
            return new Div(new Mul(l, r.left), r.right).simplify();
        }
        if (l instanceof Div) {
            return new Div(new Mul(l.left, r), l.right).simplify();
        }

        // Commutativity with Number: x * 2 -> 2 * x
        if (r instanceof Num && !(l instanceof Num)) {
            return new Mul(r, l).simplify();
        }

        // Distribute (-1) over Subtraction: -1 * (A - B) -> B - A
        if (l instanceof Num && l.value === -1 && r instanceof Sub) {
            return new Sub(r.right, r.left).simplify();
        }
        // Distribute (-1) over Addition: -1 * (A + B) -> (-A - B)
        if (l instanceof Num && l.value === -1 && r instanceof Add) {
            return new Sub(new Mul(new Num(-1), r.left).simplify(), r.right).simplify();
        }

        // x * x -> x^2
        if (l.toString() === r.toString()) {
            return new Pow(l, new Num(2));
        }

        // x * x^n -> x^(n+1)
        if (r instanceof Pow && r.left.toString() === l.toString() && r.right instanceof Num) {
            return new Pow(l, new Num(r.right.value + 1)).simplify();
        }
        if (l instanceof Pow && l.left.toString() === r.toString() && l.right instanceof Num) {
            return new Pow(r, new Num(l.right.value + 1)).simplify();
        }

        // Logarithm Condense: n * ln(x) -> ln(x^n)
        if (l instanceof Num && r instanceof Call && (r.funcName === 'ln' || r.funcName === 'log')) {
             return new Call(r.funcName, [new Pow(r.args[0], l).simplify()]);
        }
        // ln(x) * n -> ln(x^n)
        if (r instanceof Num && l instanceof Call && (l.funcName === 'ln' || l.funcName === 'log')) {
             return new Call(l.funcName, [new Pow(l.args[0], r).simplify()]);
        }

        // Trigonometric Simplifications
        if (l instanceof Call && r instanceof Call) {
             const is = (n, name) => n.funcName === name;
             // Check argument equality
             if (l.args.length === 1 && r.args.length === 1 && l.args[0].toString() === r.args[0].toString()) {
                  const arg = l.args[0];
                  if ((is(l, 'tan') && is(r, 'cos')) || (is(l, 'cos') && is(r, 'tan'))) return new Call('sin', [arg]);
                  if ((is(l, 'cot') && is(r, 'sin')) || (is(l, 'sin') && is(r, 'cot'))) return new Call('cos', [arg]);
                  if ((is(l, 'sec') && is(r, 'cos')) || (is(l, 'cos') && is(r, 'sec'))) return new Num(1);
                  if ((is(l, 'csc') && is(r, 'sin')) || (is(l, 'sin') && is(r, 'csc'))) return new Num(1);
                  if ((is(l, 'tan') && is(r, 'cot')) || (is(l, 'cot') && is(r, 'tan'))) return new Num(1);
             }
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

        // Check for constant factors
        if (!this.left.dependsOn(varName)) {
            return new Mul(this.left, this.right.integrate(varName));
        }
        if (!this.right.dependsOn(varName)) {
            return new Mul(this.right, this.left.integrate(varName));
        }

        // If varName is not provided or invalid, default to generic Call
        if (!varName) return new Call("integrate", [this]);
        return new Call("integrate", [this, varName]);
    }
    expand() {
        const l = this.left.expand();
        const r = this.right.expand();
        if (l instanceof Add) return new Add(new Mul(l.left, r).expand(), new Mul(l.right, r).expand());
        if (l instanceof Sub) return new Sub(new Mul(l.left, r).expand(), new Mul(l.right, r).expand());
        if (r instanceof Add) return new Add(new Mul(l, r.left).expand(), new Mul(l, r.right).expand());
        if (r instanceof Sub) return new Sub(new Mul(l, r.left).expand(), new Mul(l, r.right).expand());
        return new Mul(l, r);
    }
    toLatex() {
        let lTex = this.left.toLatex();
        let rTex = this.right.toLatex();
        let op = " \\cdot ";

        // Implicit multiplication formatting rules
        if (this.left instanceof Num && this.right instanceof Sym) op = ""; // 2x
        if (this.left instanceof Num && this.right instanceof Call) op = ""; // 2sin(x)
        if (this.left instanceof Sym && this.right instanceof Sym) op = " "; // x y
        if (this.left instanceof Sym && this.right instanceof Call) op = ""; // x sin(y)
        if (this.left instanceof Num && this.right instanceof Pow && this.right.left instanceof Sym) op = ""; // 2x^2

        // Handle -1 coefficient
        if (this.left instanceof Num && this.left.value === -1) {
            if (this.right instanceof Sym || this.right instanceof Call || (this.right instanceof Pow && this.right.left instanceof Sym)) {
                return `-${this.right.toLatex()}`;
            }
            if (this.right instanceof Add || this.right instanceof Sub) {
                return `-\\left(${this.right.toLatex()}\\right)`;
            }
        }

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
            if (r.value === 0) {
                 if (l.value === 0) return new Sym("NaN"); // 0/0
                 if (l.value < 0) return new Mul(new Num(-1), new Sym("Infinity")); // -1/0 -> -Infinity
                 return new Sym("Infinity"); // 1/0
            }
            if (l.value % r.value === 0) return new Num(l.value / r.value);
            // Simplify signs: a/-b -> -a/b, -a/-b -> a/b
            if (r.value < 0) {
                return new Div(new Num(-l.value), new Num(-r.value)).simplify();
            }
            // GCD simplification for fractions
            const gcd = (a, b) => !b ? a : gcd(b, a % b);
            const common = gcd(Math.abs(l.value), Math.abs(r.value));
            if (common > 1) {
                return new Div(new Num(l.value / common), new Num(r.value / common));
            }

            return new Div(l, r);
        }

        const isInf = (n) => n instanceof Sym && (n.name === 'Infinity' || n.name === 'infinity');
        const isNegInf = (n) => n instanceof Mul && n.left instanceof Num && n.left.value === -1 && isInf(n.right);
        const checkInf = (n) => isInf(n) ? 1 : (isNegInf(n) ? -1 : 0);

        const infL = checkInf(l);
        const infR = checkInf(r);

        if (infL !== 0 && infR === 0) {
            // Inf / Finite
            if (r instanceof Num) {
                const sign = (r.value > 0 ? 1 : -1) * infL;
                return sign > 0 ? new Sym("Infinity") : new Mul(new Num(-1), new Sym("Infinity"));
            }
        }
        if (infL === 0 && infR !== 0) {
            // Finite / Inf -> 0
            return new Num(0);
        }
        if (infL !== 0 && infR !== 0) {
            // Inf / Inf
            return new Sym("NaN");
        }

        // Vector division by scalar: Vec / Scalar
        if (l instanceof Vec && !(r instanceof Vec)) {
            return new Vec(l.elements.map(e => new Div(e, r).simplify()));
        }

        // Generic sign simplification for negative denominator
        if (r instanceof Num && r.value < 0) {
            return new Div(new Mul(new Num(-1), l).simplify(), new Num(-r.value)).simplify();
        }

        if (l instanceof Num && l.value === 0) return new Num(0);
        if (r instanceof Num && r.value === 1) return l;
        if (l.toString() === r.toString()) return new Num(1);

        // Complex 1/i -> -i
        if (l instanceof Num && l.value === 1 && r instanceof Sym && r.name === 'i') {
            return new Mul(new Num(-1), new Sym('i'));
        }
        // Complex 1/-i -> i
        if (l instanceof Num && l.value === 1 && r instanceof Mul && r.left instanceof Num && r.left.value === -1 && r.right instanceof Sym && r.right.name === 'i') {
            return new Sym('i');
        }

        // Cancellation: (a * b) / a -> b
        if (l instanceof Mul) {
            if (l.left.toString() === r.toString()) return l.right;
            if (l.right.toString() === r.toString()) return l.left;

            // Handle scalar denominator vs term in numerator
            // e.g. (2 * x) / 2 -> x (where r is Num(2) and l.left is Num(2))
            if (r instanceof Num && l.left instanceof Num && l.left.value === r.value) return l.right;
            if (r instanceof Num && l.right instanceof Num && l.right.value === r.value) return l.left;
        }
        // Cancellation: a / (a * b) -> 1/b
        if (r instanceof Mul) {
            if (r.left.toString() === l.toString()) return new Div(new Num(1), r.right).simplify();
            if (r.right.toString() === l.toString()) return new Div(new Num(1), r.left).simplify();
        }

        // Cancellation: (a * b) / (a * c) -> b / c
        if (l instanceof Mul && r instanceof Mul) {
            if (l.left.toString() === r.left.toString()) return new Div(l.right, r.right).simplify();
            if (l.left.toString() === r.right.toString()) return new Div(l.right, r.left).simplify();
            if (l.right.toString() === r.left.toString()) return new Div(l.left, r.right).simplify();
            if (l.right.toString() === r.right.toString()) return new Div(l.left, r.left).simplify();
        }

        // Div of Divs: (a/b) / (c/d) -> (a*d) / (b*c)
        if (l instanceof Div && r instanceof Div) {
            return new Div(new Mul(l.left, r.right).simplify(), new Mul(l.right, r.left).simplify()).simplify();
        }
        if (l instanceof Div) { // (a/b) / c -> a / (b*c)
            return new Div(l.left, new Mul(l.right, r).simplify()).simplify();
        }
        if (r instanceof Div) { // a / (b/c) -> (a*c) / b
            return new Div(new Mul(l, r.right).simplify(), r.left).simplify();
        }

        // Recursive Reduction for Mul terms
        const isSameDiv = (res, n, d) => (res instanceof Div && res.left.toString() === n.toString() && res.right.toString() === d.toString());

        // Distribute division over Add/Sub if it simplifies
        if (l instanceof Add || l instanceof Sub) {
            const d1 = new Div(l.left, r).simplify();
            const d2 = new Div(l.right, r).simplify();

            // Check for clean cancellation
            if (r instanceof Num) {
                const hasFraction = (e) => {
                    if (e instanceof Div && e.right instanceof Num) return true;
                    if (e instanceof Mul) return hasFraction(e.left) || hasFraction(e.right);
                    if (e instanceof Add || e instanceof Sub) return hasFraction(e.left) || hasFraction(e.right);
                    return false;
                };

                if (!hasFraction(d1) && !hasFraction(d2)) {
                    if (l instanceof Add) return new Add(d1, d2).simplify();
                    if (l instanceof Sub) return new Sub(d1, d2).simplify();
                }
            } else {
                const effective = (res) => {
                    if (isSameDiv(res, l.left, r)) return false;
                    if (res instanceof Div && res.right.toString() === r.toString()) return false;
                    return true;
                };

                if (effective(d1) && effective(d2)) {
                    if (l instanceof Add) return new Add(d1, d2).simplify();
                    if (l instanceof Sub) return new Sub(d1, d2).simplify();
                }
            }
        }

        // (A * B) / C -> (A / C) * B or (B / C) * A
        if (l instanceof Mul) {
            const divLeft = new Div(l.left, r).simplify();
            if (!isSameDiv(divLeft, l.left, r)) {
                return new Mul(divLeft, l.right).simplify();
            }
            const divRight = new Div(l.right, r).simplify();
            if (!isSameDiv(divRight, l.right, r)) {
                return new Mul(divRight, l.left).simplify();
            }
        }

        // C / (A * B) -> (C / A) / B or (C / B) / A
        if (r instanceof Mul) {
             const divLeft = new Div(l, r.left).simplify();
             if (!isSameDiv(divLeft, l, r.left)) {
                 return new Div(divLeft, r.right).simplify();
             }
             const divRight = new Div(l, r.right).simplify();
             if (!isSameDiv(divRight, l, r.right)) {
                 return new Div(divRight, r.left).simplify();
             }
        }

        // Simplify Powers in Division: x^a / x^b -> x^(a-b)
        let baseL = l;
        let expL = new Num(1);
        if (l instanceof Pow) {
            baseL = l.left;
            expL = l.right;
        }

        let baseR = r;
        let expR = new Num(1);
        if (r instanceof Pow) {
            baseR = r.left;
            expR = r.right;
        }

        if (baseL.toString() === baseR.toString()) {
            const newExp = new Sub(expL, expR).simplify();
            return new Pow(baseL, newExp).simplify();
        }

        // Polynomial GCD Simplification
        const getVar = (node) => {
            if (node instanceof Sym && node.name !== 'pi' && node.name !== 'e' && node.name !== 'i') return node;
            if (node instanceof Pow && node.left instanceof Sym) return node.left;
            if (node instanceof BinaryOp) return getVar(node.left) || getVar(node.right);
            return null;
        };
        const varNode = getVar(l) || getVar(r);

        if (varNode) {
            const P = getPolyCoeffs(l, varNode);
            const Q = getPolyCoeffs(r, varNode);
            if (P && Q) {
                // Check if denominator constant (handled by coeff extraction but verify)
                if (Q.length <= 1) {
                    // Constant denominator, already handled by Num/Div rules usually, but coeff division?
                    // e.g. (2x+2)/2 -> x+1
                    // Only apply if division is exact to avoid floats (e.g. x/6 -> 0.166x)
                    if (Q.length === 1 && Q[0] !== 0) {
                        const allExact = P.every(c => Math.abs(c % Q[0]) < 1e-10);
                        if (allExact) {
                            const resCoeffs = P.map(c => c / Q[0]);
                            return polyFromCoeffs(resCoeffs, varNode);
                        }
                    }
                } else {
                    const G = polyGcd(P, Q);
                    // If gcd degree > 0 or (degree=0 and value != 1)
                    if (G.length > 1 || (G.length === 1 && Math.abs(G[0] - 1) > 1e-9)) {
                        // Divide num and den by G
                        const newNumCoeffs = polyQuotient(P, G);
                        const newDenCoeffs = polyQuotient(Q, G);

                        const newNum = polyFromCoeffs(newNumCoeffs, varNode);
                        const newDen = polyFromCoeffs(newDenCoeffs, varNode);

                        if (newDenCoeffs.length === 1 && Math.abs(newDenCoeffs[0] - 1) < 1e-9) return newNum;
                        return new Div(newNum, newDen);
                    }
                }
            }
        }

        // Factorial/Gamma Simplification
        if (l instanceof Call && r instanceof Call && l.funcName === r.funcName) {
            if (l.funcName === 'factorial') {
                const a = l.args[0];
                const b = r.args[0];
                const diff = new Sub(a, b).simplify();
                if (diff instanceof Num && Number.isInteger(diff.value)) {
                    const k = diff.value;
                    if (k > 0) {
                        let res = new Num(1);
                        for(let i=1; i<=k; i++) {
                            res = new Mul(res, new Add(b, new Num(i)));
                        }
                        return res.simplify();
                    } else if (k < 0) {
                        let res = new Num(1);
                        for(let i=1; i<=-k; i++) {
                            res = new Mul(res, new Add(a, new Num(i)));
                        }
                        return new Div(new Num(1), res).simplify();
                    }
                }
            }
            if (l.funcName === 'gamma') {
                const a = l.args[0];
                const b = r.args[0];
                const diff = new Sub(a, b).simplify();
                if (diff instanceof Num && Number.isInteger(diff.value)) {
                    const k = diff.value;
                    if (k > 0) {
                        let res = new Num(1);
                        for(let i=0; i<k; i++) {
                            res = new Mul(res, new Add(b, new Num(i)));
                        }
                        return res.simplify();
                    } else if (k < 0) {
                        const m = -k;
                        let res = new Num(1);
                        for(let i=1; i<=m; i++) {
                            res = new Mul(res, new Sub(b, new Num(i)));
                        }
                        return new Div(new Num(1), res).simplify();
                    }
                }
            }
        }

        return new Div(l, r);
    }
    evaluateNumeric() { return this.left.evaluateNumeric() / this.right.evaluateNumeric(); }
    diff(varName) {
        const num = new Sub(new Mul(this.left.diff(varName), this.right), new Mul(this.left, this.right.diff(varName)));
        const den = new Pow(this.right, new Num(2));
        return new Div(num, den);
    }
    integrate(varName) {
        // sin(x)/x -> Si(x)
        if (this.right.toString() === varName.toString()) {
            if (this.left instanceof Call && this.left.funcName === 'sin' && this.left.args[0].toString() === varName.toString()) {
                return new Call('Si', [varName]);
            }
            // cos(x)/x -> Ci(x)
            if (this.left instanceof Call && this.left.funcName === 'cos' && this.left.args[0].toString() === varName.toString()) {
                return new Call('Ci', [varName]);
            }
            // exp(x)/x -> Ei(x)
            if (this.left instanceof Call && this.left.funcName === 'exp' && this.left.args[0].toString() === varName.toString()) {
                return new Call('Ei', [varName]);
            }
        }

        // 1/ln(x) -> Li(x)
        if (this.left instanceof Num && this.left.value === 1 &&
            this.right instanceof Call && (this.right.funcName === 'ln' || this.right.funcName === 'log') &&
            this.right.args[0].toString() === varName.toString()) {
            return new Call('Li', [varName]);
        }

        const isX2 = (node) => node instanceof Pow && node.left.toString() === varName.toString() && node.right instanceof Num && node.right.value === 2;
        const isConst = (node) => (node instanceof Num) || (node instanceof Sym && node.name !== varName.name);

        if (this.left instanceof Num) {
             const c = this.left;
             const den = this.right;

             // 1/(x^2+a)
             if (den instanceof Add) {
                 let xPart, aPart;
                 if (isX2(den.left) && isConst(den.right)) { xPart = den.left; aPart = den.right; }
                 else if (isConst(den.left) && isX2(den.right)) { xPart = den.right; aPart = den.left; }

                 if (xPart) {
                     // int(c/(x^2+a)) = c * (1/sqrt(a)) * atan(x/sqrt(a))
                     const aVal = aPart.evaluateNumeric();
                     if (aVal > 0) {
                          const sqrtA = new Call('sqrt', [aPart]).simplify();
                          return new Mul(new Div(c, sqrtA), new Call('atan', [new Div(varName, sqrtA)]));
                     }
                 }

                 // Completing the square: x^2 + bx + c
                 // Look for Add(Add(Pow(x,2), Mul(b, x)), c)
                 // Or flat Add structure? This parser/simplifier usually nests Adds. ((x^2 + bx) + c)

                 // Helper to match x^2 + bx
                 const isQuadPart = (n) => {
                     if (n instanceof Add) {
                         if (isX2(n.left) && n.right instanceof Mul && n.right.right.toString() === varName.toString() && isConst(n.right.left)) {
                             return { b: n.right.left };
                         }
                         if (isX2(n.right) && n.left instanceof Mul && n.left.right.toString() === varName.toString() && isConst(n.left.left)) {
                             return { b: n.left.left };
                         }
                         // Handle x^2 + x (b=1)
                         if (isX2(n.left) && n.right.toString() === varName.toString()) return { b: new Num(1) };
                         if (isX2(n.right) && n.left.toString() === varName.toString()) return { b: new Num(1) };
                     }
                     return null;
                 };

                 // Check den = (x^2+bx) + c
                 const qp = isQuadPart(den.left);
                 if (qp && isConst(den.right)) {
                     const b = qp.b;
                     const cVal = den.right;
                     // Form: x^2 + bx + c = (x + b/2)^2 + (c - b^2/4)
                     // Let u = x + b/2. du = dx.
                     // Den = u^2 + K. K = c - b^2/4.
                     const bOver2 = new Div(b, new Num(2)).simplify();
                     const bSqOver4 = new Div(new Pow(b, new Num(2)), new Num(4)).simplify();
                     const K = new Sub(cVal, bSqOver4).simplify();

                     const KVal = K.evaluateNumeric();
                     if (KVal > 0) {
                         const sqrtK = new Call('sqrt', [K]).simplify();
                         const u = new Add(varName, bOver2).simplify();
                         // Result: (1/sqrtK) * atan(u/sqrtK)
                         return new Mul(c, new Mul(new Div(new Num(1), sqrtK), new Call('atan', [new Div(u, sqrtK)])));
                     }
                 }
             }

             // 1/sqrt(a - x^2) -> asin(x/sqrt(a))
             if (den instanceof Call && den.funcName === 'sqrt') {
                 const arg = den.args[0];
                 if (arg instanceof Sub) {
                     // a - x^2
                     if (isConst(arg.left) && isX2(arg.right)) {
                         const aPart = arg.left;
                         const aVal = aPart.evaluateNumeric();
                         if (aVal > 0) {
                              const sqrtA = new Call('sqrt', [aPart]).simplify();
                              return new Mul(c, new Call('asin', [new Div(varName, sqrtA)]));
                         }
                     }
                 }
             }
        }

        if (this.left instanceof Num && this.left.value === 1 && this.right instanceof Sym && this.right.name === varName.name) {
            return new Call("ln", [new Call("abs", [varName])]);
        }
        if (this.right instanceof Num) {
            return new Mul(new Div(new Num(1), this.right), this.left.integrate(varName));
        }

        // Helper to check if expression depends on variable
        const dependsOn = (expr, v) => {
            if (expr instanceof Sym) return expr.name === v.name;
            if (expr instanceof Num) return false;
            if (expr instanceof BinaryOp) return dependsOn(expr.left, v) || dependsOn(expr.right, v);
            if (expr instanceof Call) return expr.args.some(a => dependsOn(a, v));
            return false;
        };

        // integrate(c / (a*x+b)) -> c/a * ln(a*x+b)
        // Only if denominator is linear in x

        if (!dependsOn(this.left, varName)) {
             const c = this.left;
             const den = this.right;

             // Check if den is linear: a*x + b
             // diff(den, x) should be independent of x
             // We can't use full simplification here easily without circular dependencies or complex logic.
             // But we can check structural forms.

             // 1. Simple: x +/- b (b is constant)
             if ((den instanceof Add || den instanceof Sub) && den.left instanceof Sym && den.left.name === varName.name && !dependsOn(den.right, varName)) {
                 // x + b -> ln(x+b)
                 // x - b -> ln(x-b)
                 return new Mul(c, new Call("ln", [den]));
             }

             // 2. Linear: a*x +/- b
             // den.diff(varName) -> a
             // If den.diff(varName) is constant (does not depend on varName) and non-zero
             try {
                 const diff = den.diff(varName).simplify();
                 // Ensure we are not dividing by zero symbolic or numeric
                 let isZero = false;
                 if (diff instanceof Num && diff.value === 0) isZero = true;

                 if (!isZero && !dependsOn(diff, varName)) {
                      // Result: c/a * ln(den)
                      const coeff = new Div(c, diff);
                      return new Mul(coeff, new Call("ln", [den]));
                 }
             } catch(e) {}
        }

        // integrate(1/x^n, x)
        if (this.left instanceof Num && this.left.value === 1 &&
            this.right instanceof Pow &&
            this.right.left instanceof Sym && this.right.left.name === varName.name &&
            this.right.right instanceof Num) {

            const n = this.right.right.value;
            if (n !== 1) {
                // x^(-n) -> x^(-n+1) / (-n+1)
                const newExp = -n + 1;
                return new Div(new Pow(this.right.left, new Num(newExp)), new Num(newExp));
            } else {
                 return new Call("ln", [new Call("abs", [varName])]);
            }
        }

        // integrate(1/sqrt(x), x) -> 2sqrt(x)
        if (this.left instanceof Num && this.left.value === 1 && this.right instanceof Call && this.right.funcName === 'sqrt' && this.right.args[0].name === varName.name) {
             return new Mul(new Num(2), new Call('sqrt', [varName]));
        }

        // integrate(x / (x^2 + a)) -> 0.5 * ln(x^2 + a)
        // More generally: integrate(f'(x)/f(x)) = ln(f(x))
        // Check if numerator is derivative of denominator (up to constant)
        const denomDiff = this.right.diff(varName).simplify();
        // Check if num = k * denomDiff
        // k = num / denomDiff
        const ratio = new Div(this.left, denomDiff).simplify();

        const dependsOn2 = (expr, v) => {
            if (expr instanceof Sym) return expr.name === v.name;
            if (expr instanceof Num) return false;
            if (expr instanceof BinaryOp) return dependsOn2(expr.left, v) || dependsOn2(expr.right, v);
            if (expr instanceof Call) return expr.args.some(a => dependsOn2(a, v));
            return false;
        };

        if (!dependsOn2(ratio, varName)) {
             // result = ratio * ln(denom)
             return new Mul(ratio, new Call("ln", [this.right]));
        }

        return new Call("integrate", [this, varName]);
    }
    expand() {
        const l = this.left.expand();
        const r = this.right.expand();
        if (l instanceof Add) return new Add(new Div(l.left, r).expand(), new Div(l.right, r).expand());
        return new Div(l, r);
    }
    toLatex() {
        // Pull out negative sign
        if (this.left instanceof Num && this.left.value < 0 && this.right instanceof Num && this.right.value > 0) {
             return `-\\frac{${Math.abs(this.left.value)}}{${this.right.value}}`;
        }
        return `\\frac{${this.left.toLatex()}}{${this.right.toLatex()}}`;
    }
}

class Pow extends BinaryOp {
    toString() {
        const leftStr = this.left.toString();
        // If left is a function call like cos(x), don't wrap in parens if unnecessary
        // but current logic wraps everything in parens.
        // Users want cos(x)^2 not (cos(x)^2).
        // Check if left is Call?
        // Actually, just removing outer parens might be unsafe if context requires them (like 2*(x^2)).
        // But the user complained about simplifier output "Got (cos(x)^2)".
        // If Pow.toString returns "(...)", then Simplify returns Pow, so it prints "(...)".
        // Let's remove outer parens. Precedence handling usually belongs to the parent (Mul, Add)
        // or we need a precedence argument in toString(precedence).
        // Since we don't have precedence passing, removing outer parens is risky for "2*x^2" -> "2*x^2" (fine)
        // "x^2 + 1" -> "x^2 + 1" (fine).
        // "(x+1)^2" -> "(x+1)^2" (left handles its own parens if needed).
        // So `return "${this.left}^${this.right}"` might be better, but need to check left/right precedence.

        let l = this.left.toString();
        let r = this.right.toString();

        // Wrap left if it's an operator with lower precedence than Pow
        if (this.left instanceof Add || this.left instanceof Sub || this.left instanceof Mul || this.left instanceof Div) {
            l = `(${l})`;
        }
        // Also wrap if left is negative number?
        if (this.left instanceof Num && this.left.value < 0) {
            l = `(${l})`;
        }

        // Wrap right if complex? usually right is evaluated first?
        // x^(a+b) needs parens
        if (this.right instanceof Add || this.right instanceof Sub || this.right instanceof Mul || this.right instanceof Div) {
            r = `(${r})`;
        }

        return `${l}^${r}`;
    }
    substitute(varName, value) {
        // If we are substituting the whole power expression
        if (this.toString() === varName.toString()) return value; // Heuristic check
        // Better equality check
        if (varName instanceof Pow && this.left.toString() === varName.left.toString() && this.right.toString() === varName.right.toString()) {
             // This is weak if toString() is not unique or order matters, but for standard forms it helps.
             // Ideal would be .equals() method on all Expr.
             return value;
        }

        return new Pow(this.left.substitute(varName, value), this.right.substitute(varName, value));
    }
    simplify() {
        const l = this.left.simplify();
        const r = this.right.simplify();

        // Complex number rules for i^n
        if (l instanceof Sym && l.name === 'i' && r instanceof Num && Number.isInteger(r.value)) {
            let n = r.value;
            // Normalize n to [0, 3] if positive, or handle negative
            // i^-1 = -i, i^-2 = -1, i^-3 = i, i^-4 = 1
            // Generalized: i^n = i^(n mod 4)
            // Need correct mod for negative numbers: ((n % 4) + 4) % 4
            const rem = ((n % 4) + 4) % 4;
            if (rem === 0) return new Num(1);
            if (rem === 1) return l; // i
            if (rem === 2) return new Num(-1);
            if (rem === 3) return new Mul(new Num(-1), l); // -i
        }

        // Matrix Power M^n
        if (l instanceof Vec && r instanceof Num && Number.isInteger(r.value) && r.value >= 0) {
             const rows = l.elements.length;
             // Check square matrix
             if (rows > 0 && l.elements[0] instanceof Vec && l.elements[0].elements.length === rows) {
                 const n = r.value;
                 // Identity
                 const idRows = [];
                 for(let i=0; i<rows; i++) {
                     const row = [];
                     for(let j=0; j<rows; j++) row.push(new Num(i===j?1:0));
                     idRows.push(new Vec(row));
                 }
                 let res = new Vec(idRows);
                 if (n === 0) return res;

                 let base = l;
                 // Binary Exponentiation or just simple loop
                 for(let i=0; i<n; i++) {
                     res = new Mul(res, base).simplify();
                 }
                 return res;
             }
        }

        if (r instanceof Num && r.value === 0) return new Num(1);
        if (r instanceof Num && r.value === 1) return l;
        if (l instanceof Num && r instanceof Num) return new Num(Math.pow(l.value, r.value));

        // Rational Exponents: n ^ (a/b)
        if (l instanceof Num && r instanceof Div && r.left instanceof Num && r.right instanceof Num) {
             const base = l.value;
             const num = r.left.value;
             const den = r.right.value;
             // Try to evaluate
             // Handle negative base for odd roots
             if (base < 0 && den % 2 !== 0) {
                 const val = Math.pow(-base, num/den);
                 if (Math.abs(val - Math.round(val)) < 1e-9) {
                     const res = Math.round(val);
                     // If num is even, result is positive. If num is odd, result is negative.
                     return new Num((num % 2 === 0) ? res : -res);
                 }
             } else if (base >= 0) {
                 const val = Math.pow(base, num/den);
                 if (Math.abs(val - Math.round(val)) < 1e-9) {
                     return new Num(Math.round(val));
                 }
             }
        }

        // Infinity Powers
        const isInf = (n) => n instanceof Sym && (n.name === 'Infinity' || n.name === 'infinity');
        const isNegInf = (n) => n instanceof Mul && n.left instanceof Num && n.left.value === -1 && isInf(n.right);
        const checkInf = (n) => isInf(n) ? 1 : (isNegInf(n) ? -1 : 0);

        const infL = checkInf(l);
        if (infL !== 0) {
            if (r instanceof Num) {
                if (r.value === 0) return new Sym("NaN"); // Inf^0 indeterminate
                if (r.value > 0) {
                    if (infL === 1) return new Sym("Infinity");
                    // (-Inf)^n
                    if (Number.isInteger(r.value)) {
                        return (r.value % 2 === 0) ? new Sym("Infinity") : new Mul(new Num(-1), new Sym("Infinity"));
                    }
                    return new Sym("NaN"); // Complex?
                }
                // r < 0: 1/Inf -> 0
                return new Num(0);
            }
        }

        // Simplify (-x)^even -> x^even
        if (l instanceof Mul && l.left instanceof Num && l.left.value === -1 && r instanceof Num && Number.isInteger(r.value) && r.value % 2 === 0) {
             return new Pow(l.right, r).simplify();
        }

        // (c * x)^n -> c^n * x^n (if c, n are numbers)
        if (l instanceof Mul && l.left instanceof Num && r instanceof Num) {
             const c = l.left;
             const x = l.right;
             const cn = new Pow(c, r).simplify();
             const xn = new Pow(x, r).simplify();
             return new Mul(cn, xn).simplify();
        }

        // Simplify roots: (x^a)^(1/b) -> x^(a/b) ?
        // Simplification rule: (a^b)^c = a^(b*c)
        if (l instanceof Pow) {
             const base = l.left;
             const exp1 = l.right;
             const exp2 = r;
             return new Pow(base, new Mul(exp1, exp2).simplify()).simplify();
        }

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
        if (this.left instanceof Sym && this.left.name === varName.name) {
            if (this.right instanceof Num) {
                if (this.right.value === -1) return new Call("ln", [varName]);
                const n = this.right.value;
                return new Div(new Pow(this.left, new Num(n + 1)), new Num(n + 1));
            }
        }
        // a^x -> a^x / ln(a)
        if (this.right instanceof Sym && this.right.name === varName.name) {
             const isSimpleConst = (node) => {
                 if (node instanceof Num) return true;
                 if (node instanceof Sym && node.name !== varName.name) return true;
                 return false;
             };
             if (isSimpleConst(this.left)) {
                 if (this.left instanceof Sym && this.left.name === 'e') return this;
                 return new Div(this, new Call('ln', [this.left]));
             }
        }
        return new Call("integrate", [this, varName]);
    }
    expand() {
        const l = this.left.expand();
        const r = this.right.expand();

        // (a * b)^n -> a^n * b^n
        if (l instanceof Mul && r instanceof Num && Number.isInteger(r.value)) {
            return new Mul(new Pow(l.left, r).expand(), new Pow(l.right, r).expand()).simplify();
        }

        // (a / b)^n -> a^n / b^n
        if (l instanceof Div && r instanceof Num && Number.isInteger(r.value)) {
            return new Div(new Pow(l.left, r).expand(), new Pow(l.right, r).expand()).simplify();
        }

        // (a + b)^n
        if (r instanceof Num && Number.isInteger(r.value) && r.value > 1) {
            const n = r.value;

            if (l instanceof Add || l instanceof Sub) {
                const a = l.left;
                const b = l.right;
                const isSub = (l instanceof Sub);

                // Binomial expansion: sum( nCk * a^(n-k) * b^k )
                // For sub: sum( nCk * a^(n-k) * (-b)^k ) => alternating signs if k is odd

                // Helper to compute factorial
                const factorial = (num) => {
                    if (num <= 1) return 1;
                    let res = 1;
                    for(let i=2; i<=num; i++) res *= i;
                    return res;
                };

                // Helper to compute nCk
                const nCr = (n, k) => {
                     return factorial(n) / (factorial(k) * factorial(n - k));
                };

                let result = null;

                for(let k=0; k<=n; k++) {
                    const coeffVal = nCr(n, k);
                    const coeff = new Num(coeffVal);

                    // Term: coeff * a^(n-k) * b^k
                    // Simplify powers: x^0 = 1, x^1 = x

                    let termA = null;
                    if (n - k === 0) termA = new Num(1);
                    else if (n - k === 1) termA = a;
                    else termA = new Pow(a, new Num(n - k));

                    let termB = null;
                    if (k === 0) termB = new Num(1);
                    else if (k === 1) termB = b;
                    else termB = new Pow(b, new Num(k));

                    // If subtraction and k is odd, term is negative
                    let isNeg = false;
                    if (isSub && k % 2 !== 0) isNeg = true;

                    let term = new Mul(coeff, new Mul(termA, termB));

                    if (result === null) {
                         if (isNeg) result = new Mul(new Num(-1), term);
                         else result = term;
                    } else {
                        if (isNeg) result = new Sub(result, term);
                        else result = new Add(result, term);
                    }
                }

                return result.expand().simplify();
            }
        }

        return new Pow(l, r);
    }
    toLatex() {
        let lTex = this.left.toLatex();
        if (this.left instanceof Add || this.left instanceof Sub || this.left instanceof Mul || this.left instanceof Div || this.left instanceof Pow) lTex = `\\left(${lTex}\\right)`;
        return `{${lTex}}^{${this.right.toLatex()}}`;
    }
}

class Call extends Expr {
    constructor(funcName, args) {
        super();
        this.funcName = funcName;
        this.args = args;
    }
    toString() {
        if (this.funcName === 'set') {
            return `{${this.args.join(", ")}}`;
        }
        return `${this.funcName}(${this.args.join(", ")})`;
    }
    simplify() {
        const simpleArgs = this.args.map(a => a.simplify());

        if (this.funcName === 'sqrt') {
            const arg = simpleArgs[0];
            // sqrt(x^2) -> abs(x)
            if (arg instanceof Pow && arg.right instanceof Num && arg.right.value === 2) {
                return new Call('abs', [arg.left]);
            }
            if (arg instanceof Num) {
                if (arg.value >= 0) {
                    const sqrtVal = Math.sqrt(arg.value);
                    if (Number.isInteger(sqrtVal)) return new Num(sqrtVal);

                    // Simplify radical for integer: sqrt(8) -> 2*sqrt(2)
                    if (Number.isInteger(arg.value)) {
                        const simp = math_intSqrtSimplify(arg.value);
                        if (simp && simp.coeff !== 1) {
                            return new Mul(new Num(simp.coeff), new Call('sqrt', [new Num(simp.radical)])).simplify();
                        }
                    }

                    return new Call('sqrt', simpleArgs);
                }
                // sqrt(-x)
                const absVal = Math.abs(arg.value);
                const sqrtVal = Math.sqrt(absVal);
                if (Number.isInteger(sqrtVal)) {
                    return new Mul(new Sym('i'), new Num(sqrtVal)).simplify();
                }
                return new Mul(new Sym('i'), new Call('sqrt', [new Num(absVal)])).simplify();
            }
            return new Call('sqrt', simpleArgs);
        }

        if (this.funcName === 'sin') {
            const arg = simpleArgs[0];
            // sin(asin(x)) -> x
            if (arg instanceof Call && arg.funcName === 'asin') return arg.args[0];
            const val = arg.evaluateNumeric();
            if (!isNaN(val)) {
                if (Math.abs(val) < 1e-10) return new Num(0);
                // Check if multiple of pi/2
                const multiple = val / (Math.PI / 2);
                const eps = 1e-10;
                if (Math.abs(multiple - Math.round(multiple)) < eps) {
                    const k = Math.round(multiple);
                    // k * pi/2
                    // k % 4: 0 -> 0, 1 -> 1, 2 -> 0, 3 -> -1
                    const rem = ((k % 4) + 4) % 4;
                    if (rem === 0) return new Num(0);
                    if (rem === 1) return new Num(1);
                    if (rem === 2) return new Num(0);
                    if (rem === 3) return new Num(-1);
                }
            }
        }
        if (this.funcName === 'cos') {
            const arg = simpleArgs[0];
            // cos(acos(x)) -> x
            if (arg instanceof Call && arg.funcName === 'acos') return arg.args[0];
            const val = arg.evaluateNumeric();
            if (!isNaN(val)) {
                if (Math.abs(val) < 1e-10) return new Num(1);
                const multiple = val / (Math.PI / 2);
                const eps = 1e-10;
                if (Math.abs(multiple - Math.round(multiple)) < eps) {
                    const k = Math.round(multiple);
                    const rem = ((k % 4) + 4) % 4;
                    if (rem === 0) return new Num(1);
                    if (rem === 1) return new Num(0);
                    if (rem === 2) return new Num(-1);
                    if (rem === 3) return new Num(0);
                }
            }
        }
        if (this.funcName === 'tan') {
            const arg = simpleArgs[0];
            // tan(atan(x)) -> x
            if (arg instanceof Call && arg.funcName === 'atan') return arg.args[0];
            const val = arg.evaluateNumeric();
            if (!isNaN(val)) {
                if (Math.abs(val) < 1e-10) return new Num(0);
                const multiple = val / (Math.PI / 4);
                const eps = 1e-10;
                if (Math.abs(multiple - Math.round(multiple)) < eps) {
                    const k = Math.round(multiple);
                    const rem = ((k % 8) + 8) % 8;
                    if (rem === 0) return new Num(0);
                    if (rem === 1) return new Num(1);
                    if (rem === 2) return new Sym("Infinity");
                    if (rem === 3) return new Num(-1);
                    if (rem === 4) return new Num(0);
                    if (rem === 5) return new Num(1);
                    if (rem === 6) return new Sym("Infinity");
                    if (rem === 7) return new Num(-1);
                }
            }
        }
        if (this.funcName === 'asin') {
            const arg = simpleArgs[0];
            // asin(sin(x)) -> x (simplified, ignoring domain issues for now)
            if (arg instanceof Call && arg.funcName === 'sin') return arg.args[0];

            const val = arg.evaluateNumeric();
            if (!isNaN(val)) {
                if (Math.abs(val) < 1e-9) return new Num(0);
                if (Math.abs(val - 1) < 1e-9) return new Div(new Sym('pi'), new Num(2)).simplify();
                if (Math.abs(val + 1) < 1e-9) return new Div(new Mul(new Num(-1), new Sym('pi')), new Num(2)).simplify();
                if (Math.abs(val - 0.5) < 1e-9) return new Div(new Sym('pi'), new Num(6)).simplify();
                if (Math.abs(val + 0.5) < 1e-9) return new Div(new Mul(new Num(-1), new Sym('pi')), new Num(6)).simplify();
                if (Math.abs(val - Math.sqrt(2)/2) < 1e-9) return new Div(new Sym('pi'), new Num(4)).simplify();
                if (Math.abs(val + Math.sqrt(2)/2) < 1e-9) return new Div(new Mul(new Num(-1), new Sym('pi')), new Num(4)).simplify();
                if (Math.abs(val - Math.sqrt(3)/2) < 1e-9) return new Div(new Sym('pi'), new Num(3)).simplify();
                if (Math.abs(val + Math.sqrt(3)/2) < 1e-9) return new Div(new Mul(new Num(-1), new Sym('pi')), new Num(3)).simplify();
            }
        }
        if (this.funcName === 'acos') {
            const arg = simpleArgs[0];
            // acos(cos(x)) -> x
            if (arg instanceof Call && arg.funcName === 'cos') return arg.args[0];

            const val = arg.evaluateNumeric();
            if (!isNaN(val)) {
                if (Math.abs(val - 1) < 1e-9) return new Num(0);
                if (Math.abs(val) < 1e-9) return new Div(new Sym('pi'), new Num(2)).simplify();
                if (Math.abs(val + 1) < 1e-9) return new Sym('pi');
                if (Math.abs(val - 0.5) < 1e-9) return new Div(new Sym('pi'), new Num(3)).simplify();
                if (Math.abs(val + 0.5) < 1e-9) return new Div(new Mul(new Num(2), new Sym('pi')), new Num(3)).simplify();
                if (Math.abs(val - Math.sqrt(2)/2) < 1e-9) return new Div(new Sym('pi'), new Num(4)).simplify();
                if (Math.abs(val + Math.sqrt(2)/2) < 1e-9) return new Div(new Mul(new Num(3), new Sym('pi')), new Num(4)).simplify();
                if (Math.abs(val - Math.sqrt(3)/2) < 1e-9) return new Div(new Sym('pi'), new Num(6)).simplify();
                if (Math.abs(val + Math.sqrt(3)/2) < 1e-9) return new Div(new Mul(new Num(5), new Sym('pi')), new Num(6)).simplify();
            }
        }
        if (this.funcName === 'atan') {
            const arg = simpleArgs[0];
            // atan(tan(x)) -> x
            if (arg instanceof Call && arg.funcName === 'tan') return arg.args[0];
            if (arg instanceof Num && arg.value === 0) return new Num(0);
            if (arg instanceof Num && arg.value === 1) return new Div(new Sym('pi'), new Num(4)).simplify();
            if (arg instanceof Num && arg.value === -1) return new Div(new Mul(new Num(-1), new Sym('pi')), new Num(4)).simplify();
            // sqrt(3) -> pi/3, 1/sqrt(3) -> pi/6
            if (arg instanceof Num && Math.abs(arg.value - Math.sqrt(3)) < 1e-9) return new Div(new Sym('pi'), new Num(3)).simplify();
            if (arg instanceof Num && Math.abs(arg.value - 1/Math.sqrt(3)) < 1e-9) return new Div(new Sym('pi'), new Num(6)).simplify();
            if (arg instanceof Num && Math.abs(arg.value - -Math.sqrt(3)) < 1e-9) return new Div(new Mul(new Num(-1), new Sym('pi')), new Num(3)).simplify();
            if (arg instanceof Num && Math.abs(arg.value - -1/Math.sqrt(3)) < 1e-9) return new Div(new Mul(new Num(-1), new Sym('pi')), new Num(6)).simplify();
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
        if (this.funcName === 'sec') {
            const arg = simpleArgs[0];
            const val = arg.evaluateNumeric();
            if (!isNaN(val)) {
                // sec(x) = 1/cos(x)
                if (Math.abs(val) < 1e-10) return new Num(1);
                const cosVal = Math.cos(val);
                if (Math.abs(cosVal) > 1e-10) return new Num(1 / cosVal);
                // undefined
            }
        }
        if (this.funcName === 'csc') {
             const arg = simpleArgs[0];
             const val = arg.evaluateNumeric();
             if (!isNaN(val)) {
                 // csc(x) = 1/sin(x)
                 const sinVal = Math.sin(val);
                 if (Math.abs(sinVal) > 1e-10) return new Num(1 / sinVal);
                 // undefined
             }
        }
        if (this.funcName === 'cot') {
             const arg = simpleArgs[0];
             const val = arg.evaluateNumeric();
             if (!isNaN(val)) {
                 // cot(x) = 1/tan(x) = cos(x)/sin(x)
                 if (Math.abs(val - Math.PI/2) < 1e-10) return new Num(0);
                 const tanVal = Math.tan(val);
                 if (Math.abs(tanVal) > 1e-10) return new Num(1 / tanVal);
                 // undefined or infinity?
             }
        }

        if (this.funcName === 'ln') {
            const arg = simpleArgs[0];
            // ln(exp(x)) -> x
            if (arg instanceof Call && arg.funcName === 'exp') return arg.args[0];
            // ln(e^x) -> x
            if (arg instanceof Pow && arg.left instanceof Sym && arg.left.name === 'e') return arg.right;
            // ln(e) -> 1
            if (arg instanceof Sym && arg.name === 'e') return new Num(1);
            if (arg instanceof Num && arg.value === 1) return new Num(0);
            if (arg instanceof Sym && (arg.name === 'Infinity' || arg.name === 'infinity')) return new Sym('Infinity');
            if (arg instanceof Num && arg.value === 0) return new Mul(new Num(-1), new Sym('Infinity'));
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
        if (this.funcName === 'log2') {
             const arg = simpleArgs[0];
             if (arg instanceof Num && arg.value === 1) return new Num(0);
             if (arg instanceof Num && arg.value === 2) return new Num(1);
             if (arg instanceof Num) {
                 const val = Math.log2(arg.value);
                 if (Number.isInteger(val)) return new Num(val);
             }
        }
        if (this.funcName === 'exp') {
            const arg = simpleArgs[0];
            // exp(ln(x)) -> x
            if (arg instanceof Call && arg.funcName === 'ln') return arg.args[0];
            // exp(n * ln(x)) -> x^n
            if (arg instanceof Mul) {
                if (arg.right instanceof Call && arg.right.funcName === 'ln') {
                    // n * ln(x) -> x^n
                    return new Pow(arg.right.args[0], arg.left).simplify();
                }
                if (arg.left instanceof Call && arg.left.funcName === 'ln') {
                    // ln(x) * n -> x^n
                    return new Pow(arg.left.args[0], arg.right).simplify();
                }
            }
            if (arg instanceof Num && arg.value === 0) return new Num(1);
            if (arg instanceof Num && arg.value === 1) return new Sym('e');
            if (arg instanceof Sym && (arg.name === 'Infinity' || arg.name === 'infinity')) return new Sym('Infinity');
            // exp(-Infinity) -> 0
            if (arg instanceof Mul && arg.left instanceof Num && arg.left.value === -1 && (arg.right.name === 'Infinity' || arg.right.name === 'infinity')) return new Num(0);
            if (arg instanceof Num && arg.value === -Infinity) return new Num(0);
        }
        if (this.funcName === 'sign') {
            const arg = simpleArgs[0];
            if (arg instanceof Num) return new Num(Math.sign(arg.value));
            // sign(c * x) -> sign(x) if c > 0, -sign(x) if c < 0
            if (arg instanceof Mul && arg.left instanceof Num) {
                const c = arg.left;
                if (c.value > 0) return new Call('sign', [arg.right]).simplify();
                if (c.value < 0) return new Mul(new Num(-1), new Call('sign', [arg.right])).simplify();
            }
            // sign(x^even) -> 1 (assuming x != 0)
            if (arg instanceof Pow && arg.right instanceof Num && arg.right.value % 2 === 0) {
                return new Num(1);
            }
            // sign(abs(x)) -> 1
            if (arg instanceof Call && arg.funcName === 'abs') {
                return new Num(1);
            }
        }

        if (this.funcName === 'Ei') {
            const arg = simpleArgs[0];
            if (arg instanceof Num) {
                return new Num(math_Ei(arg.value));
            }
        }
        if (this.funcName === 'Li') {
            const arg = simpleArgs[0];
            if (arg instanceof Num) {
                const val = arg.value;
                if (val > 0 && val !== 1) {
                    return new Num(math_Ei(Math.log(val)));
                }
            }
        }

        if (this.funcName === 'EllipticK') {
            const arg = simpleArgs[0];
            if (arg instanceof Num) {
                return new Num(math_EllipticK(arg.value));
            }
        }
        if (this.funcName === 'EllipticE') {
            const arg = simpleArgs[0];
            if (arg instanceof Num) {
                return new Num(math_EllipticE(arg.value));
            }
        }

        if (this.funcName === 'erf') {
            const arg = simpleArgs[0];
            if (arg instanceof Num) {
                return new Num(math_erf(arg.value));
            }
            if (arg instanceof Mul && arg.left instanceof Num && arg.left.value < 0) {
                // erf(-x) = -erf(x)
                return new Mul(new Num(-1), new Call('erf', [new Mul(new Num(-arg.left.value), arg.right).simplify()])).simplify();
            }
            if (arg instanceof Num && arg.value === 0) return new Num(0);
        }

        if (this.funcName === 'erfc') {
            const arg = simpleArgs[0];
            if (arg instanceof Num) {
                // erfc(0) = 1
                if (arg.value === 0) return new Num(1);
                return new Num(1 - math_erf(arg.value));
            }
            // erfc(-x) = 2 - erfc(x)
            if (arg instanceof Mul && arg.left instanceof Num && arg.left.value < 0) {
                 // erfc(-x) = 2 - erfc(x)
                 return new Sub(new Num(2), new Call('erfc', [new Mul(new Num(-arg.left.value), arg.right).simplify()])).simplify();
            }
            if (arg instanceof Num && arg.value === 0) return new Num(1);
        }

        if (this.funcName === 'erfi') {
            const arg = simpleArgs[0];
            if (arg instanceof Num && arg.value === 0) return new Num(0);
            if (arg instanceof Num && arg.value < 0) {
                 return new Mul(new Num(-1), new Call('erfi', [new Num(-arg.value)])).simplify();
            }
            if (arg instanceof Mul && arg.left instanceof Num && arg.left.value < 0) {
                 return new Mul(new Num(-1), new Call('erfi', [new Mul(new Num(-1), arg).simplify()])).simplify();
            }
            // erfi(i*x) = i*erf(x)
            if (arg instanceof Mul && arg.left instanceof Sym && arg.left.name === 'i') {
                 return new Mul(new Sym('i'), new Call('erf', [arg.right])).simplify();
            }
            if (arg instanceof Sym && arg.name === 'i') {
                 return new Mul(new Sym('i'), new Call('erf', [new Num(1)])).simplify();
            }
        }

        if (this.funcName === 'erfinv') {
            const arg = simpleArgs[0];
            if (arg instanceof Num && arg.value === 0) return new Num(0);
            if (arg instanceof Num) {
                return new Num(math_erfinv(arg.value));
            }
            // erfinv(-x) = -erfinv(x)
            if (arg instanceof Mul && arg.left instanceof Num && arg.left.value < 0) {
                return new Mul(new Num(-1), new Call('erfinv', [new Mul(new Num(-arg.left.value), arg.right).simplify()])).simplify();
            }
        }

        if (this.funcName === 'power') {
             // Do not simplify to Pow or Num, keep as 'power' to preserve factorization structure
             return new Call('power', simpleArgs);
        }

        if (this.funcName === 'floor') {
            const arg = simpleArgs[0];
            if (arg instanceof Num) return new Num(Math.floor(arg.value));
        }
        if (this.funcName === 'ceil') {
            const arg = simpleArgs[0];
            if (arg instanceof Num) return new Num(Math.ceil(arg.value));
        }
        if (this.funcName === 'round') {
            const arg = simpleArgs[0];
            const decimals = simpleArgs.length > 1 ? simpleArgs[1] : null;
            if (arg instanceof Num) {
                if (decimals instanceof Num) {
                     const factor = Math.pow(10, decimals.value);
                     return new Num(Math.round(arg.value * factor) / factor);
                }
                if (!decimals) return new Num(Math.round(arg.value));
            }
        }
        if (this.funcName === 'min' || this.funcName === 'max') {
            const allNumeric = simpleArgs.every(a => a instanceof Num);
            if (allNumeric) {
                const vals = simpleArgs.map(a => a.value);
                if (this.funcName === 'min') return new Num(Math.min(...vals));
                if (this.funcName === 'max') return new Num(Math.max(...vals));
            }
            // Flatten nested min/max if same type? min(min(a,b), c) -> min(a,b,c)
            const flatArgs = [];
            for(const a of simpleArgs) {
                if (a instanceof Call && a.funcName === this.funcName) {
                    flatArgs.push(...a.args);
                } else {
                    flatArgs.push(a);
                }
            }
            if (flatArgs.length !== simpleArgs.length) {
                return new Call(this.funcName, flatArgs).simplify();
            }
        }
        if (this.funcName === 'piecewise') {
            // piecewise(cond1, val1, cond2, val2, ..., [default])
            for(let i=0; i<simpleArgs.length; i+=2) {
                if (i + 1 >= simpleArgs.length) {
                    // Default case (last argument)
                    return simpleArgs[i];
                }
                const cond = simpleArgs[i];
                const val = simpleArgs[i+1];
                // Check if cond is true/false
                if (cond instanceof Num) {
                    if (cond.value !== 0) return val; // True -> return val
                    // False -> continue to next case
                } else if (cond instanceof Sym && (cond.name === 'true' || cond.name === 'false')) {
                    if (cond.name === 'true') return val;
                } else {
                    // Symbolic condition, cannot simplify further branches easily without assuming false
                    // But we can reconstruct the piecewise with remaining args
                    return new Call('piecewise', simpleArgs.slice(i));
                }
            }
            return new Num(NaN); // No cases matched
        }
        if (this.funcName === 'abs') {
            const arg = simpleArgs[0];
            if (arg instanceof Num) return new Num(Math.abs(arg.value));
        }
        if (this.funcName === 'real') {
            const arg = simpleArgs[0];
            if (arg instanceof Num) return arg;
            if (arg instanceof Sym && arg.name === 'i') return new Num(0);
        }
        if (this.funcName === 'imag') {
             const arg = simpleArgs[0];
             if (arg instanceof Num) return new Num(0);
             if (arg instanceof Sym && arg.name === 'i') return new Num(1);
        }
        if (this.funcName === 'conj') {
             const arg = simpleArgs[0];
             if (arg instanceof Num) return arg;
             if (arg instanceof Sym && arg.name === 'i') return new Mul(new Num(-1), arg);
        }

        if (this.funcName === 'lambertw') {
            const arg = simpleArgs[0];
            if (arg instanceof Num) {
                if (arg.value === 0) return new Num(0);
                if (Math.abs(arg.value - Math.E) < 1e-9) return new Num(1);
                if (Math.abs(arg.value + 1/Math.E) < 1e-9) return new Num(-1);
            }
            if (arg instanceof Sym && arg.name === 'e') return new Num(1);
        }

        if (this.funcName === 'psi' || this.funcName === 'digamma') {
            const val = simpleArgs[0].evaluateNumeric();
            if (!isNaN(val)) {
                return new Num(math_psi(val));
            }
        }

        if (this.funcName === 'zeta') {
            const arg = simpleArgs[0];
            if (arg instanceof Num) {
                if (arg.value === 0) return new Num(-0.5);
                if (arg.value === 1) return new Sym('Infinity');
                if (Number.isInteger(arg.value)) {
                    const n = arg.value;
                    if (n > 0 && n % 2 === 0) {
                        // zeta(2n) = (-1)^(n+1) * B_2n * (2pi)^(2n) / (2 * (2n)!)
                        const k = n / 2;
                        const B = getBernoulliExpr(n);
                        const piPow = new Pow(new Sym('pi'), new Num(n));
                        const twoPow = new Pow(new Num(2), new Num(n - 1));

                        // Factorial (2n)!
                        let fact = 1;
                        for(let i=2; i<=n; i++) fact *= i;

                        // (-1)^(k+1) * B * 2^(2n-1) * pi^(2n) / (2n)!
                        // Formula simplified: zeta(2k) = |B_2k| * (2pi)^(2k) / (2 * (2k)!)
                        // = |B_2k| * 2^(2k) * pi^(2k) / 2 / (2k)!
                        // = |B_2k| * 2^(2k-1) * pi^(2k) / (2k)!

                        // zeta(2k) is always positive.
                        // Formula: zeta(2k) = (-1)^(k+1) * B_2k * (2pi)^2k / 2(2k)!
                        // We use simpler: |B_2k| * (2pi)^2k / 2(2k)!
                        // If B_2k is Num, we can just take abs. If symbolic (Div), we handle sign.

                        let absB = B;
                        if (B instanceof Num) absB = new Num(Math.abs(B.value));
                        else if (B instanceof Div && B.left instanceof Num) absB = new Div(new Num(Math.abs(B.left.value)), B.right).simplify();
                        else absB = new Call('abs', [B]);

                        const num = new Mul(absB, new Mul(new Pow(new Num(2), new Num(n-1)), piPow));
                        const den = new Num(fact);
                        return new Div(num, den).simplify();
                    }
                    if (n < 0 && Number.isInteger(n)) {
                        // zeta(-n) = (-1)^n * B_(n+1) / (n+1)
                        // B_n is non-zero only for even n (except B1)
                        // If n is even negative integer (e.g. -2, -4), n+1 is odd => B_odd=0.
                        // So zeta(-2k) = 0.
                        if (n % 2 === 0) return new Num(0);

                        // If n is odd negative integer (e.g. -1, -3), n+1 is even => B_even != 0.
                        // zeta(-1) = -1 * B_2 / 2 = -1/12.
                        // zeta(-3) = -1 * B_4 / 4 = -(-1/30)/4 = 1/120.
                        const k = -n; // Positive
                        const B = getBernoulliExpr(k + 1);
                        const sign = (k % 2 === 0) ? 1 : -1; // (-1)^(-k) = (-1)^k
                        return new Div(new Mul(new Num(sign), B), new Num(k + 1)).simplify();
                    }
                }
            }
        }

        if (this.funcName === 'heaviside') {
            const arg = simpleArgs[0];
            if (arg instanceof Num) {
                if (arg.value < 0) return new Num(0);
                if (arg.value > 0) return new Num(1);
                // Heaviside(0) depends on convention. Standard in many CAS is 1/2 or undefined.
                // We return 0.5
                return new Num(0.5);
            }
        }

        if (this.funcName === 'dirac') {
            const arg = simpleArgs[0];
            if (arg instanceof Num) {
                if (arg.value !== 0) return new Num(0);
                return new Sym("Infinity");
            }
        }

        if (this.funcName === 'sinc') {
            // sinc(x) = sin(x)/x
            const arg = simpleArgs[0];
            if (arg instanceof Num) {
                if (arg.value === 0) return new Num(1); // Limit
                return new Num(Math.sin(arg.value) / arg.value);
            }
        }

        if (this.funcName === 'rect') {
            // rect(x) = 1 if |x|<0.5, 0.5 if |x|=0.5, 0 else
            const arg = simpleArgs[0];
            if (arg instanceof Num) {
                const abs = Math.abs(arg.value);
                if (abs < 0.5) return new Num(1);
                if (Math.abs(abs - 0.5) < 1e-15) return new Num(0.5);
                return new Num(0);
            }
        }

        if (this.funcName === 'sigmoid') {
            const arg = simpleArgs[0];
            if (arg instanceof Num) return new Num(1 / (1 + Math.exp(-arg.value)));
        }
        if (this.funcName === 'relu') {
            const arg = simpleArgs[0];
            if (arg instanceof Num) return new Num(Math.max(0, arg.value));
        }
        if (this.funcName === 'softplus') {
            const arg = simpleArgs[0];
            if (arg instanceof Num) return new Num(Math.log(1 + Math.exp(arg.value)));
        }

        if (this.funcName === 'gamma') {
            const arg = simpleArgs[0];
            if (arg instanceof Div && arg.left instanceof Num && arg.right instanceof Num) {
                // gamma(1/2) = sqrt(pi)
                if (arg.left.value === 1 && arg.right.value === 2) {
                    return new Call('sqrt', [new Sym('pi')]);
                }
                // gamma(3/2) = 0.5 * sqrt(pi)
                if (arg.left.value === 3 && arg.right.value === 2) {
                    return new Mul(new Num(0.5), new Call('sqrt', [new Sym('pi')])).simplify();
                }
                // gamma(5/2) = 0.75 * sqrt(pi)
                if (arg.left.value === 5 && arg.right.value === 2) {
                    return new Mul(new Num(0.75), new Call('sqrt', [new Sym('pi')])).simplify();
                }
            }
            // gamma(n) = (n-1)! for integer n
            if (arg instanceof Num && Number.isInteger(arg.value) && arg.value > 0) {
                // Factorial logic is usually in _factorial, but simplify can do small ones
                if (arg.value <= 10) {
                    let res = 1;
                    for(let i=1; i<arg.value; i++) res *= i;
                    return new Num(res);
                }
            }
        }

        if (['besselJ', 'besselY', 'besselI', 'besselK', 'hyp2f1'].includes(this.funcName)) {
             const allNum = simpleArgs.every(a => a instanceof Num);
             if (allNum) {
                 const val = new Call(this.funcName, simpleArgs).evaluateNumeric();
                 if (!isNaN(val)) return new Num(val);
             }
        }

        return new Call(this.funcName, simpleArgs);
    }
    evaluateNumeric() {
        const argsVal = this.args.map(a => a.evaluateNumeric());
        if (this.funcName === 'heaviside') {
            const v = argsVal[0];
            if (v < 0) return 0;
            if (v > 0) return 1;
            return 0.5;
        }
        if (this.funcName === 'dirac') {
            const v = argsVal[0];
            if (v === 0) return Infinity;
            return 0;
        }
        if (this.funcName === 'sinc') {
            const v = argsVal[0];
            if (Math.abs(v) < 1e-15) return 1;
            return Math.sin(v) / v;
        }
        if (this.funcName === 'rect') {
            const v = Math.abs(argsVal[0]);
            if (v < 0.5) return 1;
            if (Math.abs(v - 0.5) < 1e-15) return 0.5;
            return 0;
        }
        if (this.funcName === 'sin') return Math.sin(argsVal[0]);
        if (this.funcName === 'cos') return Math.cos(argsVal[0]);
        if (this.funcName === 'tan') return Math.tan(argsVal[0]);
        if (this.funcName === 'asin') return Math.asin(argsVal[0]);
        if (this.funcName === 'acos') return Math.acos(argsVal[0]);
        if (this.funcName === 'atan') return Math.atan(argsVal[0]);
        if (this.funcName === 'asinh') return Math.asinh(argsVal[0]);
        if (this.funcName === 'acosh') return Math.acosh(argsVal[0]);
        if (this.funcName === 'atanh') return Math.atanh(argsVal[0]);
        if (this.funcName === 'sinh') return Math.sinh(argsVal[0]);
        if (this.funcName === 'cosh') return Math.cosh(argsVal[0]);
        if (this.funcName === 'tanh') return Math.tanh(argsVal[0]);
        if (this.funcName === 'exp') return Math.exp(argsVal[0]);
        if (this.funcName === 'sec') return 1 / Math.cos(argsVal[0]);
        if (this.funcName === 'csc') return 1 / Math.sin(argsVal[0]);
        if (this.funcName === 'cot') return 1 / Math.tan(argsVal[0]);
        if (this.funcName === 'asec') return Math.acos(1 / argsVal[0]);
        if (this.funcName === 'acsc') return Math.asin(1 / argsVal[0]);
        if (this.funcName === 'acot') return Math.atan(1 / argsVal[0]);
        if (this.funcName === 'sech') return 1 / Math.cosh(argsVal[0]);
        if (this.funcName === 'csch') return 1 / Math.sinh(argsVal[0]);
        if (this.funcName === 'coth') return 1 / Math.tanh(argsVal[0]);
        if (this.funcName === 'asinh') return Math.asinh(argsVal[0]); // Already there but confirming
        if (this.funcName === 'acosh') return Math.acosh(argsVal[0]);
        if (this.funcName === 'atanh') return Math.atanh(argsVal[0]);
        if (this.funcName === 'asech') return Math.acosh(1 / argsVal[0]);
        if (this.funcName === 'acsch') return Math.asinh(1 / argsVal[0]);
        if (this.funcName === 'acoth') return Math.atanh(1 / argsVal[0]);

        if (this.funcName === 'ln') return Math.log(argsVal[0]);
        if (this.funcName === 'log') return Math.log10(argsVal[0]);
        if (this.funcName === 'log2') return Math.log2(argsVal[0]);
        if (this.funcName === 'sqrt') return Math.sqrt(argsVal[0]);
        if (this.funcName === 'abs') return Math.abs(argsVal[0]);
        if (this.funcName === 'sign') return Math.sign(argsVal[0]);
        if (this.funcName === 'floor') return Math.floor(argsVal[0]);
        if (this.funcName === 'ceil') return Math.ceil(argsVal[0]);
        if (this.funcName === 'round') {
            if (argsVal.length > 1 && !isNaN(argsVal[1])) {
                const factor = Math.pow(10, argsVal[1]);
                return Math.round(argsVal[0] * factor) / factor;
            }
            return Math.round(argsVal[0]);
        }
        if (this.funcName === 'min') return Math.min(...argsVal);
        if (this.funcName === 'max') return Math.max(...argsVal);
        if (this.funcName === 'piecewise') {
             for(let i=0; i<argsVal.length; i+=2) {
                 if (i + 1 >= argsVal.length) return argsVal[i]; // Default
                 if (argsVal[i] !== 0) return argsVal[i+1];
             }
             return NaN;
        }
        if (this.funcName === 'erf') {
             return math_erf(argsVal[0]);
        }
        if (this.funcName === 'erfc') {
             return 1 - math_erf(argsVal[0]);
        }
        if (this.funcName === 'Q') {
             return 0.5 * (1 - math_erf(argsVal[0] / Math.SQRT2));
        }
        if (this.funcName === 'erfinv') {
             return math_erfinv(argsVal[0]);
        }
        if (this.funcName === 'Si') {
            return math_Si(argsVal[0]);
        }
        if (this.funcName === 'Ci') {
            return math_Ci(argsVal[0]);
        }
        if (this.funcName === 'FresnelS') {
            return math_fresnelS(argsVal[0]);
        }
        if (this.funcName === 'FresnelC') {
            return math_fresnelC(argsVal[0]);
        }
        if (this.funcName === 'EllipticK') {
            return math_EllipticK(argsVal[0]);
        }
        if (this.funcName === 'EllipticE') {
            return math_EllipticE(argsVal[0]);
        }
        if (this.funcName === 'Ei') {
            return math_Ei(argsVal[0]);
        }
        if (this.funcName === 'Li') {
            const val = argsVal[0];
            if (val <= 0 || val === 1) return NaN;
            return math_Ei(Math.log(val));
        }
        if (this.funcName === 'gamma') {
            return math_gamma(argsVal[0]);
        }
        if (this.funcName === 'factorial') {
            return math_gamma(argsVal[0] + 1);
        }
        if (this.funcName === 'beta') {
            return math_beta(argsVal[0], argsVal[1]);
        }
        if (this.funcName === 'zeta') {
            return math_zeta(argsVal[0]);
        }
        if (this.funcName === 'lambertw') {
            return math_lambertw(argsVal[0]);
        }
        if (this.funcName === 'psi' || this.funcName === 'digamma') {
            return math_psi(argsVal[0]);
        }
        if (this.funcName === 'polygamma') {
            // polygamma(n, z)
            if (argsVal.length === 2 && !isNaN(argsVal[0]) && !isNaN(argsVal[1])) {
                return math_polygamma(argsVal[0], argsVal[1]);
            }
        }
        if (this.funcName === 'besselJ') return math_besselJ(argsVal[0], argsVal[1]);
        if (this.funcName === 'besselY') return math_besselY(argsVal[0], argsVal[1]);
        if (this.funcName === 'besselI') return math_besselI(argsVal[0], argsVal[1]);
        if (this.funcName === 'besselK') return math_besselK(argsVal[0], argsVal[1]);
        if (this.funcName === 'hyp2f1') return math_hyp2f1(argsVal[0], argsVal[1], argsVal[2], argsVal[3]);

        if (this.funcName === 'sigmoid') return 1 / (1 + Math.exp(-argsVal[0]));
        if (this.funcName === 'relu') return Math.max(0, argsVal[0]);
        if (this.funcName === 'softplus') return Math.log(1 + Math.exp(argsVal[0]));

        return NaN; // Unknown
    }
    diff(varName) {
        const u = this.args[0];
        if (this.funcName === 'gamma') {
            // diff(gamma(u)) = gamma(u) * psi(u) * u'
            return new Mul(new Mul(this, new Call('psi', [u])), u.diff(varName));
        }
        if (this.funcName === 'psi' || this.funcName === 'digamma') {
            // diff(psi(u)) = polygamma(1, u) * u'
            return new Mul(new Call('polygamma', [new Num(1), u]), u.diff(varName));
        }
        if (this.funcName === 'polygamma') {
            // diff(polygamma(n, x)) = polygamma(n+1, x)
            if (this.args.length === 2) {
                const n = this.args[0];
                const x = this.args[1];
                // Assuming n is constant w.r.t x for standard polygamma derivative
                if (n instanceof Num) {
                    return new Mul(new Call('polygamma', [new Num(n.value + 1), x]), x.diff(varName));
                }
                return new Mul(new Call('polygamma', [new Add(n, new Num(1)), x]), x.diff(varName));
            }
        }
        if (this.funcName === 'heaviside') {
            // diff(heaviside(u)) = dirac(u) * u'
            return new Mul(new Call('dirac', [u]), u.diff(varName));
        }
        if (this.funcName === 'dirac') {
            // diff(dirac(u)) = dirac(u, 1) * u' ? Or symbolic
            // We don't support distributional derivatives fully yet.
            // Standard: dirac'(x)
            return new Call('diff', [this, varName]);
        }
        if (this.funcName === 'sign') {
            // diff(sign(u)) = 2*dirac(u) * u'
            return new Mul(new Mul(new Num(2), new Call('dirac', [u])), u.diff(varName));
        }
        if (['floor', 'ceil', 'round'].includes(this.funcName)) {
            return new Num(0);
        }
        if (this.funcName === 'Ei') {
            // Ei'(u) = exp(u)/u * u'
            return new Mul(new Div(new Call('exp', [u]), u), u.diff(varName));
        }
        if (this.funcName === 'Li') {
            // Li'(u) = 1/ln(u) * u'
            return new Mul(new Div(new Num(1), new Call('ln', [u])), u.diff(varName));
        }
        if (this.funcName === 'Si') {
            // Si'(u) = sin(u)/u * u'
            return new Mul(new Div(new Call('sin', [u]), u), u.diff(varName));
        }
        if (this.funcName === 'Ci') {
            // Ci'(u) = cos(u)/u * u'
            return new Mul(new Div(new Call('cos', [u]), u), u.diff(varName));
        }
        if (this.funcName === 'sinc') {
            // diff(sin(x)/x) = (x*cos(x) - sin(x)) / x^2
            const num = new Sub(new Mul(u, new Call('cos', [u])), new Call('sin', [u]));
            const den = new Pow(u, new Num(2));
            return new Mul(new Div(num, den), u.diff(varName));
        }
        if (this.funcName === 'rect') {
            // diff(rect(t)) = dirac(t+0.5) - dirac(t-0.5)
            // rect(t) = H(t+0.5) - H(t-0.5)
            const t = u;
            const term1 = new Call('dirac', [new Add(t, new Num(0.5))]);
            const term2 = new Call('dirac', [new Sub(t, new Num(0.5))]);
            return new Mul(new Sub(term1, term2), u.diff(varName));
        }

        if (this.funcName === 'lambertw') {
            // W'(x) = W(x) / (x * (1 + W(x)))
            // if x=0, limit is 1.
            // But if u is 0? W(0)=0. denom = 0 * 1 = 0.
            // Diff of W(x) at x=0 is 1. (x=0 -> W=0 -> limit W/x = 1)
            // But formula gives 0/0.
            // Using exp(W)*W' = 1 -> W' = exp(-W)/(1+W)
            // Or W' = W / (x(1+W))
            // The formula W' = exp(-W)/(1+W) is better for x=0.
            // W(0) = 0. exp(0)/(1+0) = 1. Correct.
            // W = lambertw(u)
            const W = new Call('lambertw', [u]);
            const num = new Call('exp', [new Mul(new Num(-1), W)]);
            const den = new Add(new Num(1), W);
            return new Mul(new Div(num, den), u.diff(varName));
        }
        if (this.funcName === 'sigmoid') {
            const s = new Call('sigmoid', [u]);
            return new Mul(new Mul(s, new Sub(new Num(1), s)), u.diff(varName));
        }
        if (this.funcName === 'relu') {
            return new Mul(new Call('heaviside', [u]), u.diff(varName));
        }
        if (this.funcName === 'softplus') {
            return new Mul(new Call('sigmoid', [u]), u.diff(varName));
        }
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
        if (this.funcName === 'asinh') {
             // 1/sqrt(u^2 + 1) * u'
             return new Mul(new Div(new Num(1), new Call('sqrt', [new Add(new Pow(u, new Num(2)), new Num(1))])), u.diff(varName));
        }
        if (this.funcName === 'acosh') {
             // 1/sqrt(u^2 - 1) * u'
             return new Mul(new Div(new Num(1), new Call('sqrt', [new Sub(new Pow(u, new Num(2)), new Num(1))])), u.diff(varName));
        }
        if (this.funcName === 'atanh') {
             // 1/(1-u^2) * u'
             return new Mul(new Div(new Num(1), new Sub(new Num(1), new Pow(u, new Num(2)))), u.diff(varName));
        }
        if (this.funcName === 'sinh') return new Mul(new Call('cosh', [u]), u.diff(varName));
        if (this.funcName === 'cosh') return new Mul(new Call('sinh', [u]), u.diff(varName));
        if (this.funcName === 'tanh') return new Mul(new Div(new Num(1), new Pow(new Call('cosh', [u]), new Num(2))), u.diff(varName));
        if (this.funcName === 'sech') return new Mul(new Mul(new Num(-1), new Mul(new Call('sech', [u]), new Call('tanh', [u]))), u.diff(varName));
        if (this.funcName === 'csch') return new Mul(new Mul(new Num(-1), new Mul(new Call('csch', [u]), new Call('coth', [u]))), u.diff(varName));
        if (this.funcName === 'coth') return new Mul(new Mul(new Num(-1), new Pow(new Call('csch', [u]), new Num(2))), u.diff(varName));

        if (this.funcName === 'sec') return new Mul(new Mul(new Call('sec', [u]), new Call('tan', [u])), u.diff(varName));
        if (this.funcName === 'csc') return new Mul(new Mul(new Num(-1), new Mul(new Call('csc', [u]), new Call('cot', [u]))), u.diff(varName));
        if (this.funcName === 'cot') return new Mul(new Mul(new Num(-1), new Pow(new Call('csc', [u]), new Num(2))), u.diff(varName));

        if (this.funcName === 'asech') {
             // -1 / (x * sqrt(1-x^2))
             return new Mul(new Div(new Num(-1), new Mul(u, new Call('sqrt', [new Sub(new Num(1), new Pow(u, new Num(2))) ]))), u.diff(varName));
        }
        if (this.funcName === 'acsch') {
             // -1 / (|x| * sqrt(1+x^2))
             return new Mul(new Div(new Num(-1), new Mul(new Call('abs', [u]), new Call('sqrt', [new Add(new Num(1), new Pow(u, new Num(2))) ]))), u.diff(varName));
        }
        if (this.funcName === 'acoth') {
             // 1 / (1-x^2)
             return new Mul(new Div(new Num(1), new Sub(new Num(1), new Pow(u, new Num(2)))), u.diff(varName));
        }

        if (this.funcName === 'exp') return new Mul(this, u.diff(varName));
        if (this.funcName === 'ln') return new Div(u.diff(varName), u);
        if (this.funcName === 'log') {
            // d/dx log10(u) = u' / (u * ln(10))
            return new Div(u.diff(varName), new Mul(u, new Call('ln', [new Num(10)])));
        }
        if (this.funcName === 'log2') {
             return new Div(u.diff(varName), new Mul(u, new Call('ln', [new Num(2)])));
        }
        if (this.funcName === 'sqrt') return new Div(u.diff(varName), new Mul(new Num(2), new Call('sqrt', [u])));
        if (this.funcName === 'abs') return new Mul(new Call('sign', [u]), u.diff(varName));
        if (this.funcName === 'min' || this.funcName === 'max') {
             // Derivative of min(f, g) is f' * step(g-f) + g' * step(f-g) roughly (at intersections undefined)
             // We can return derivative of the evaluated piecewise if possible, but structure is dynamic.
             // Simplest: piecewise(diff(f), f<g, diff(g), g<=f)
             // For now, symbolic.
             return new Call('diff', [this, varName]);
        }
        if (this.funcName === 'piecewise') {
            const newArgs = [];
            for(let i=0; i<this.args.length; i++) {
                // Odd indices are values, Even are conditions.
                // deriv(piecewise) = piecewise(cond1, diff(val1), cond2, diff(val2))
                // Conditions remain same (ignoring Dirac deltas at boundaries)
                if (i % 2 === 1 || (i === this.args.length - 1 && this.args.length % 2 === 1)) {
                    newArgs.push(this.args[i].diff(varName));
                } else {
                    newArgs.push(this.args[i]);
                }
            }
            return new Call('piecewise', newArgs).simplify();
        }
        if (this.funcName === 'erf') {
            // d/dx erf(u) = 2/sqrt(pi) * e^(-u^2) * u'
            const coeff = new Div(new Num(2), new Call('sqrt', [new Sym('pi')]));
            const exp = new Call('exp', [new Mul(new Num(-1), new Pow(u, new Num(2)))]);
            return new Mul(new Mul(coeff, exp), u.diff(varName));
        }
        if (this.funcName === 'erfc') {
            // d/dx erfc(u) = -2/sqrt(pi) * e^(-u^2) * u'
            const coeff = new Div(new Num(-2), new Call('sqrt', [new Sym('pi')]));
            const exp = new Call('exp', [new Mul(new Num(-1), new Pow(u, new Num(2)))]);
            return new Mul(new Mul(coeff, exp), u.diff(varName));
        }
        if (this.funcName === 'erfi') {
            // d/dx erfi(u) = 2/sqrt(pi) * e^(u^2) * u'
            const coeff = new Div(new Num(2), new Call('sqrt', [new Sym('pi')]));
            const exp = new Call('exp', [new Pow(u, new Num(2))]);
            return new Mul(new Mul(coeff, exp), u.diff(varName));
        }

        if (this.funcName === 'erfinv') {
            // d/dx erfinv(u) = sqrt(pi)/2 * exp(erfinv(u)^2) * u'
            const coeff = new Div(new Call('sqrt', [new Sym('pi')]), new Num(2));
            const exp = new Call('exp', [new Pow(new Call('erfinv', [u]), new Num(2))]);
            return new Mul(new Mul(coeff, exp), u.diff(varName));
        }

        if (this.funcName === 'FresnelS') {
            // S'(x) = sin(pi/2 x^2)
            const arg = new Mul(new Div(new Sym('pi'), new Num(2)), new Pow(u, new Num(2)));
            return new Mul(new Call('sin', [arg]), u.diff(varName));
        }
        if (this.funcName === 'FresnelC') {
            // C'(x) = cos(pi/2 x^2)
            const arg = new Mul(new Div(new Sym('pi'), new Num(2)), new Pow(u, new Num(2)));
            return new Mul(new Call('cos', [arg]), u.diff(varName));
        }

        if (this.funcName === 'EllipticK') {
            // dK/dk = E(k)/(k(1-k^2)) - K(k)/k
            // = (E(k) - (1-k^2)K(k)) / (k(1-k^2))
            const k = u;
            const E = new Call('EllipticE', [k]);
            const K = new Call('EllipticK', [k]);
            const k2 = new Pow(k, new Num(2));
            const oneMinusK2 = new Sub(new Num(1), k2);

            const num = new Sub(E, new Mul(oneMinusK2, K));
            const den = new Mul(k, oneMinusK2);
            return new Mul(new Div(num, den), k.diff(varName));
        }

        if (this.funcName === 'EllipticE') {
            // dE/dk = (E(k) - K(k)) / k
            const k = u;
            const E = new Call('EllipticE', [k]);
            const K = new Call('EllipticK', [k]);
            const num = new Sub(E, K);
            return new Mul(new Div(num, k), k.diff(varName));
        }

        // Bessel Functions Derivatives
        // J_v'(x) = 0.5 * (J_{v-1}(x) - J_{v+1}(x))
        if (this.funcName === 'besselJ') {
            if (this.args.length === 2) {
                const v = this.args[0];
                const x = this.args[1];
                // diff wrt x (argument)
                // If v depends on varName, it's much more complex. Assuming order v is constant.
                const term1 = new Call('besselJ', [new Sub(v, new Num(1)), x]);
                const term2 = new Call('besselJ', [new Add(v, new Num(1)), x]);
                const deriv = new Mul(new Num(0.5), new Sub(term1, term2));
                return new Mul(deriv, x.diff(varName));
            }
        }
        // Y_v'(x) = 0.5 * (Y_{v-1}(x) - Y_{v+1}(x))
        if (this.funcName === 'besselY') {
            if (this.args.length === 2) {
                const v = this.args[0];
                const x = this.args[1];
                const term1 = new Call('besselY', [new Sub(v, new Num(1)), x]);
                const term2 = new Call('besselY', [new Add(v, new Num(1)), x]);
                const deriv = new Mul(new Num(0.5), new Sub(term1, term2));
                return new Mul(deriv, x.diff(varName));
            }
        }
        // I_v'(x) = 0.5 * (I_{v-1}(x) + I_{v+1}(x))
        if (this.funcName === 'besselI') {
            if (this.args.length === 2) {
                const v = this.args[0];
                const x = this.args[1];
                const term1 = new Call('besselI', [new Sub(v, new Num(1)), x]);
                const term2 = new Call('besselI', [new Add(v, new Num(1)), x]);
                const deriv = new Mul(new Num(0.5), new Add(term1, term2));
                return new Mul(deriv, x.diff(varName));
            }
        }
        // K_v'(x) = -0.5 * (K_{v-1}(x) + K_{v+1}(x))
        if (this.funcName === 'besselK') {
            if (this.args.length === 2) {
                const v = this.args[0];
                const x = this.args[1];
                const term1 = new Call('besselK', [new Sub(v, new Num(1)), x]);
                const term2 = new Call('besselK', [new Add(v, new Num(1)), x]);
                const deriv = new Mul(new Num(-0.5), new Add(term1, term2));
                return new Mul(deriv, x.diff(varName));
            }
        }
        if (this.funcName === 'hyp2f1') {
            if (this.args.length === 4) {
                const a = this.args[0];
                const b = this.args[1];
                const c = this.args[2];
                const z = this.args[3];
                // d/dz hyp2f1(a,b,c,z) = a*b/c * hyp2f1(a+1, b+1, c+1, z) * z'
                const coeff = new Div(new Mul(a, b), c);
                const nextF = new Call('hyp2f1', [
                    new Add(a, new Num(1)),
                    new Add(b, new Num(1)),
                    new Add(c, new Num(1)),
                    z
                ]);
                return new Mul(new Mul(coeff, nextF), z.diff(varName));
            }
        }
        // Airy Functions
        // Ai'(x) -> derivative of Ai is usually denoted Ai'(x).
        // Ai''(x) = x Ai(x).
        // We can't express Ai'(x) in terms of Ai(x) simply without defining AiPrime.
        // Let's leave it symbolic Call('diff', ...) for now unless we introduce airyAiPrime.

        if (this.funcName === 'sum') {
            // diff(sum(f, i, a, b), x) = sum(diff(f, x), i, a, b)
            // Assuming limits a, b do not depend on x
            if (this.args.length === 4) {
                const expr = this.args[0];
                const i = this.args[1];
                const a = this.args[2];
                const b = this.args[3];
                // Check if bounds depend on varName?
                // Standard CAS assumption: differentiation under sum sign
                // If bounds depend on x, Leibniz rule applies.
                // For now, assume independent bounds or handle inner term
                return new Call('sum', [expr.diff(varName), i, a, b]);
            }
        }

        if (this.funcName === 'product') {
            // diff(product(f, i, a, b), x) = product(f) * sum(diff(f)/f, i, a, b)
            if (this.args.length === 4) {
                const expr = this.args[0];
                const i = this.args[1];
                const a = this.args[2];
                const b = this.args[3];
                const term = new Div(expr.diff(varName), expr);
                return new Mul(this, new Call('sum', [term, i, a, b]));
            }
        }

        // Default to symbolic diff
        return new Call('diff', [this, varName]);
    }
    integrate(varName) {
        if (this.args[0].toString() === varName.toString()) {
            if (this.funcName === 'sin') return new Mul(new Num(-1), new Call('cos', [varName]));
            if (this.funcName === 'cos') return new Call('sin', [varName]);
            if (this.funcName === 'tan') return new Mul(new Num(-1), new Call('ln', [new Call('abs', [new Call('cos', [varName])])]));
            if (this.funcName === 'cot') return new Call('ln', [new Call('abs', [new Call('sin', [varName])])]);
            if (this.funcName === 'sec') return new Call('ln', [new Call('abs', [new Add(new Call('sec', [varName]), new Call('tan', [varName]))])]);
            if (this.funcName === 'csc') return new Mul(new Num(-1), new Call('ln', [new Call('abs', [new Add(new Call('csc', [varName]), new Call('cot', [varName]))])]));
            if (this.funcName === 'sinh') return new Call('cosh', [varName]);
            if (this.funcName === 'cosh') return new Call('sinh', [varName]);
            if (this.funcName === 'tanh') return new Call('ln', [new Call('cosh', [varName])]);
            if (this.funcName === 'sech') return new Call('atan', [new Call('sinh', [varName])]);
            if (this.funcName === 'csch') return new Call('ln', [new Call('tanh', [new Div(varName, new Num(2))])]);
            if (this.funcName === 'coth') return new Call('ln', [new Call('sinh', [varName])]);
            if (this.funcName === 'acoth') {
                // x*acoth(x) + 0.5*ln(x^2-1)
                const term1 = new Mul(varName, new Call('acoth', [varName]));
                const term2 = new Mul(new Num(0.5), new Call('ln', [new Sub(new Pow(varName, new Num(2)), new Num(1))]));
                return new Add(term1, term2);
            }
            if (this.funcName === 'exp') return this;
            if (this.funcName === 'ln') return new Sub(new Mul(varName, new Call('ln', [varName])), varName);
            if (this.funcName === 'log') return new Div(new Sub(new Mul(varName, new Call('ln', [varName])), varName), new Call('ln', [new Num(10)]));

            // Inverse Trig Integrals
            // integrate(asin(x)) = x*asin(x) + sqrt(1-x^2)
            if (this.funcName === 'asin') {
                return new Add(
                    new Mul(varName, new Call('asin', [varName])),
                    new Call('sqrt', [new Sub(new Num(1), new Pow(varName, new Num(2)))])
                );
            }
            // integrate(acos(x)) = x*acos(x) - sqrt(1-x^2)
            if (this.funcName === 'acos') {
                return new Sub(
                    new Mul(varName, new Call('acos', [varName])),
                    new Call('sqrt', [new Sub(new Num(1), new Pow(varName, new Num(2)))])
                );
            }
            // integrate(atan(x)) = x*atan(x) - 0.5*ln(1+x^2)
            if (this.funcName === 'atan') {
                return new Sub(
                    new Mul(varName, new Call('atan', [varName])),
                    new Mul(new Num(0.5), new Call('ln', [new Add(new Num(1), new Pow(varName, new Num(2)))]))
                );
            }
            if (this.funcName === 'dirac') {
                return new Call('heaviside', [varName]);
            }
            if (this.funcName === 'heaviside') {
                // x * H(x) (Ramp function)
                return new Mul(varName, new Call('heaviside', [varName]));
            }
            if (this.funcName === 'abs') {
                // x*abs(x)/2
                return new Div(new Mul(varName, new Call('abs', [varName])), new Num(2));
            }
            if (this.funcName === 'sign') {
                // abs(x)
                return new Call('abs', [varName]);
            }
            if (this.funcName === 'Si') {
                // x*Si(x) + cos(x)
                return new Add(new Mul(varName, new Call('Si', [varName])), new Call('cos', [varName]));
            }
            if (this.funcName === 'Ci') {
                // x*Ci(x) - sin(x)
                return new Sub(new Mul(varName, new Call('Ci', [varName])), new Call('sin', [varName]));
            }
            if (this.funcName === 'Ei') {
                // x*Ei(x) - exp(x)
                return new Sub(new Mul(varName, new Call('Ei', [varName])), new Call('exp', [varName]));
            }
            if (this.funcName === 'Li') {
                // x*Li(x) - Ei(2*ln(x))
                const term2 = new Call('Ei', [new Mul(new Num(2), new Call('ln', [varName]))]);
                return new Sub(new Mul(varName, new Call('Li', [varName])), term2);
            }
            if (this.funcName === 'erf') {
                // x*erf(x) + exp(-x^2)/sqrt(pi)
                const term1 = new Mul(varName, new Call('erf', [varName]));
                const term2 = new Div(new Call('exp', [new Mul(new Num(-1), new Pow(varName, new Num(2)))]), new Call('sqrt', [new Sym('pi')]));
                return new Add(term1, term2);
            }
            if (this.funcName === 'erfc') {
                // x*erfc(x) - exp(-x^2)/sqrt(pi)
                const term1 = new Mul(varName, new Call('erfc', [varName]));
                const term2 = new Div(new Call('exp', [new Mul(new Num(-1), new Pow(varName, new Num(2)))]), new Call('sqrt', [new Sym('pi')]));
                return new Sub(term1, term2);
            }
            if (this.funcName === 'erfi') {
                // x*erfi(x) - exp(x^2)/sqrt(pi)
                const term1 = new Mul(varName, new Call('erfi', [varName]));
                const term2 = new Div(new Call('exp', [new Pow(varName, new Num(2))]), new Call('sqrt', [new Sym('pi')]));
                return new Sub(term1, term2);
            }
            if (this.funcName === 'FresnelS') {
                // x S(x) + 1/pi cos(pi/2 x^2)
                const term1 = new Mul(varName, new Call('FresnelS', [varName]));
                const arg = new Mul(new Div(new Sym('pi'), new Num(2)), new Pow(varName, new Num(2)));
                const term2 = new Div(new Call('cos', [arg]), new Sym('pi'));
                return new Add(term1, term2);
            }
            if (this.funcName === 'FresnelC') {
                // x C(x) - 1/pi sin(pi/2 x^2)
                const term1 = new Mul(varName, new Call('FresnelC', [varName]));
                const arg = new Mul(new Div(new Sym('pi'), new Num(2)), new Pow(varName, new Num(2)));
                const term2 = new Div(new Call('sin', [arg]), new Sym('pi'));
                return new Sub(term1, term2);
            }
            if (this.funcName === 'erfinv') {
                // integral(erfinv(x)) = x*erfinv(x) + exp(-erfinv(x)^2)/sqrt(pi) is WRONG.
                // Correct: -exp(-erfinv(x)^2)/sqrt(pi) (plus constant) because d/dx(-1/sqrt(pi) * e^(-y^2)) = ... = y.
                // Let y = erfinv(x). x = erf(y). dx = 2/sqrt(pi) e^-y^2 dy.
                // int y dx = int y * 2/sqrt(pi) e^-y^2 dy = -1/sqrt(pi) e^-y^2.
                const arg = new Call('erfinv', [varName]);
                const expTerm = new Call('exp', [new Mul(new Num(-1), new Pow(arg, new Num(2)))]);
                const res = new Div(expTerm, new Call('sqrt', [new Sym('pi')]));
                return new Mul(new Num(-1), res);
            }
            if (this.funcName === 'sigmoid') {
                return new Call('softplus', [varName]);
            }
            if (this.funcName === 'relu') {
                return new Mul(new Num(0.5), new Mul(varName, new Call('relu', [varName])));
            }
            if (this.funcName === 'softplus') {
                // -Li2(-e^x)
                return new Mul(new Num(-1), new Call('polylog', [new Num(2), new Mul(new Num(-1), new Call('exp', [varName]))]));
            }
        }
        return new Call("integrate", [this, varName]);
    }
    substitute(varName, value) {
        // Check if we are substituting the whole function call
        if (this.toString() === varName.toString()) return value; // Heuristic
        // Better equality check
        if (varName instanceof Call && this.funcName === varName.funcName && this.args.length === varName.args.length) {
             let allMatch = true;
             for(let i=0; i<this.args.length; i++) {
                 if (this.args[i].toString() !== varName.args[i].toString()) {
                     allMatch = false; break;
                 }
             }
             if (allMatch) return value;
        }

        // Handle variable binding functions (integrate, sum, product, limit, diff, solve, etc.)
        // These functions typically take the bound variable as the second argument (index 1).
        const bindingFunctions = ['integrate', 'sum', 'product', 'limit', 'diff', 'solve', 'plot', 'taylor'];
        if (bindingFunctions.includes(this.funcName) && this.args.length >= 2) {
            const boundVar = this.args[1];
            if (boundVar instanceof Sym && boundVar.name === varName.name) {
                // The variable being substituted is the bound variable.
                // Do NOT substitute in the body (arg 0) or the bound variable definition (arg 1).
                // However, limits/ranges (args 2+) might still need substitution if they are distinct from the bound variable scope.
                const newArgs = [...this.args];
                for (let i = 2; i < this.args.length; i++) {
                    newArgs[i] = newArgs[i].substitute(varName, value);
                }
                return new Call(this.funcName, newArgs);
            }
        }
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
            'sec': '\\sec',
            'csc': '\\csc',
            'cot': '\\cot',
            'asec': '\\operatorname{arcsec}',
            'acsc': '\\operatorname{arccsc}',
            'acot': '\\operatorname{arccot}',
            'asinh': '\\operatorname{asinh}',
            'acosh': '\\operatorname{acosh}',
            'atanh': '\\operatorname{atanh}',
            'sech': '\\operatorname{sech}',
            'csch': '\\operatorname{csch}',
            'coth': '\\operatorname{coth}',
            'lcm': '\\operatorname{lcm}',
            'floor': '\\lfloor ' + argsTex[0] + ' \\rfloor',
            'ceil': '\\lceil ' + argsTex[0] + ' \\rceil',
            'trace': '\\operatorname{tr}',
            'real': '\\Re',
            'imag': '\\Im',
            'conj': '\\overline{' + argsTex[0] + '}',
            'sign': '\\operatorname{sgn}',
            'erf': '\\operatorname{erf}',
            'psi': '\\psi',
            'digamma': '\\psi',
            'Si': '\\operatorname{Si}',
            'Ci': '\\operatorname{Ci}',
            'Ei': '\\operatorname{Ei}',
            'Li': '\\operatorname{Li}',
            'erfi': '\\operatorname{erfi}',
            'erfinv': '\\operatorname{erf}^{-1}',
            'FresnelS': 'S',
            'FresnelC': 'C',
            'EllipticK': 'K',
            'EllipticE': 'E'
        };

        if (this.funcName === 'besselJ' && argsTex.length === 2) return `J_{${argsTex[0]}}\\left(${argsTex[1]}\\right)`;
        if (this.funcName === 'besselY' && argsTex.length === 2) return `Y_{${argsTex[0]}}\\left(${argsTex[1]}\\right)`;
        if (this.funcName === 'besselI' && argsTex.length === 2) return `I_{${argsTex[0]}}\\left(${argsTex[1]}\\right)`;
        if (this.funcName === 'besselK' && argsTex.length === 2) return `K_{${argsTex[0]}}\\left(${argsTex[1]}\\right)`;
        if (this.funcName === 'airyAi') return `\\operatorname{Ai}\\left(${argsTex[0]}\\right)`;
        if (this.funcName === 'airyBi') return `\\operatorname{Bi}\\left(${argsTex[0]}\\right)`;

        if (this.funcName === 'hyp2f1' && argsTex.length === 4) {
             return `{}_2F_1\\left(${argsTex[0]}, ${argsTex[1]}; ${argsTex[2]}; ${argsTex[3]}\\right)`;
        }

        if (this.funcName === 'polygamma' && argsTex.length === 2) {
            return `\\psi^{(${argsTex[0]})}\\left(${argsTex[1]}\\right)`;
        }

        if (this.funcName === 'floor') return `\\lfloor ${argsTex[0]} \\rfloor`;
        if (this.funcName === 'ceil') return `\\lceil ${argsTex[0]} \\rceil`;
        if (this.funcName === 'conj') return `\\overline{${argsTex[0]}}`;

        if (mapFunctions.hasOwnProperty(this.funcName)) {
            return `${mapFunctions[this.funcName]}\\left(${argsTex.join(", ")}\\right)`;
        }

        if (this.funcName === 'inv' && argsTex.length === 1) {
             return `\\left(${argsTex[0]}\\right)^{-1}`;
        }
        if (this.funcName === 'log2' && argsTex.length === 1) {
             return `\\log_2\\left(${argsTex[0]}\\right)`;
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
                let exprTex = argsTex[0];
                if (this.args[0] instanceof Add || this.args[0] instanceof Sub) {
                    exprTex = `\\left(${exprTex}\\right)`;
                }
                return `\\lim_{${argsTex[1]} \\to ${argsTex[2]}} ${exprTex}`;
            }
        }

        if (this.funcName === 'integrate') {
            let exprTex = argsTex[0];
            if (this.args[0] instanceof Add || this.args[0] instanceof Sub) {
                exprTex = `\\left(${exprTex}\\right)`;
            }
            // integrate(expr, var) or integrate(expr, var, a, b)
            if (argsTex.length === 2) {
                return `\\int ${exprTex} \\, d${argsTex[1]}`;
            }
            if (argsTex.length === 4) {
                return `\\int_{${argsTex[2]}}^{${argsTex[3]}} ${exprTex} \\, d${argsTex[1]}`;
            }
        }

        if (this.funcName === 'sum') {
            // sum(expr, var, start, end)
            if (argsTex.length === 4) {
                let exprTex = argsTex[0];
                if (this.args[0] instanceof Add || this.args[0] instanceof Sub) {
                    exprTex = `\\left(${exprTex}\\right)`;
                }
                return `\\sum_{${argsTex[1]}=${argsTex[2]}}^{${argsTex[3]}} ${exprTex}`;
            }
        }

        if (this.funcName === 'product') {
            // product(expr, var, start, end)
            if (argsTex.length === 4) {
                let exprTex = argsTex[0];
                if (this.args[0] instanceof Add || this.args[0] instanceof Sub) {
                    exprTex = `\\left(${exprTex}\\right)`;
                }
                return `\\prod_{${argsTex[1]}=${argsTex[2]}}^{${argsTex[3]}} ${exprTex}`;
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

        if (this.funcName === 'set') {
            return `\\left\\{ ${argsTex.join(", ")} \\right\\}`;
        }

        if (this.funcName === 'nCr' && argsTex.length === 2) {
            return `\\binom{${argsTex[0]}}{${argsTex[1]}}`;
        }

        if (this.funcName === 'nPr' && argsTex.length === 2) {
             return `{}_{${argsTex[0]}}P_{${argsTex[1]}}`;
        }

        if (this.funcName === 'power') {
             return `{${argsTex[0]}}^{${argsTex[1]}}`;
        }

        if (this.funcName === 'grad' && argsTex.length === 2) {
             return `\\nabla ${argsTex[0]}`;
        }
        if (this.funcName === 'curl' && argsTex.length === 2) {
             return `\\nabla \\times ${argsTex[0]}`;
        }
        if ((this.funcName === 'divergence' || this.funcName === 'div') && argsTex.length === 2) {
             return `\\nabla \\cdot ${argsTex[0]}`;
        }
        if (this.funcName === 'norm' && argsTex.length === 1) {
             return `\\left\\|${argsTex[0]}\\right\\|`;
        }
        if (this.funcName === 'size' && argsTex.length === 1) {
             return `\\left|${argsTex[0]}\\right|`;
        }
        if (this.funcName === 'arg' && argsTex.length === 1) {
             return `\\arg\\left(${argsTex[0]}\\right)`;
        }
        if (this.funcName === 'approx') {
             return argsTex[0];
        }

        if (this.funcName === 'piecewise') {
            let s = "\\begin{cases} ";
            for(let i=0; i<argsTex.length; i+=2) {
                if (i + 1 >= argsTex.length) {
                    s += argsTex[i] + " & \\text{otherwise}";
                } else {
                    s += argsTex[i+1] + " & \\text{if } " + argsTex[i] + " \\\\ ";
                }
            }
            s += " \\end{cases}";
            return s;
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
    evaluateNumeric() {
        const l = this.left.evaluateNumeric();
        const r = this.right.evaluateNumeric();
        if (!isNaN(l) && !isNaN(r)) return (Math.abs(l - r) < 1e-10) ? 1 : 0;
        return NaN;
    }
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
        const isMatrix = this.elements.length > 0 && this.elements.every(e => e instanceof Vec);

        if (isMatrix) {
            // Matrix
            const rows = this.elements.map(row => {
                if (row instanceof Vec) {
                    return row.elements.map(e => e.toLatex()).join(" & ");
                }
                return row.toLatex();
            }).join(" \\\\ ");
            return `\\begin{bmatrix} ${rows} \\end{bmatrix}`;
        }
        // Vector (Row or Column? Xcas [1,2] is row-like but usually displayed as column in LaTeX if we want vertical?)
        // Standard convention for [a, b, c] is often row vector or just list.
        // But for linear algebra tools, column vectors are standard for v.
        // If we use \\\\ it becomes a column vector.
        return `\\begin{bmatrix} ${this.elements.map(e => e.toLatex()).join(" \\\\ ")} \\end{bmatrix}`;
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

class Block extends Expr {
    constructor(statements) {
        super();
        this.statements = statements;
    }
    toString() { return this.statements.map(s => s.toString()).join("; "); }
    simplify() { return new Block(this.statements.map(s => s.simplify())); }
    evaluateNumeric() { return NaN; }
    toLatex() {
        return this.statements.map(s => s.toLatex()).join("; \\; ");
    }
}

class And extends BinaryOp {
    toString() { return `(${this.left} and ${this.right})`; }
    simplify() {
        const l = this.left.simplify();
        const r = this.right.simplify();
        if (l instanceof Num && r instanceof Num) {
            return new Num((l.value && r.value) ? 1 : 0);
        }
        // Short-circuit logic: false && x -> false
        if (l instanceof Num && l.value === 0) return new Num(0);
        if (r instanceof Num && r.value === 0) return new Num(0);
        // true && x -> x
        if (l instanceof Num && l.value !== 0) return r;
        if (r instanceof Num && r.value !== 0) return l;

        // Idempotence: A and A -> A
        if (l.toString() === r.toString()) return l;

        // Absorption: A and (A or B) -> A
        if (r instanceof Or) {
            if (r.left.toString() === l.toString() || r.right.toString() === l.toString()) return l;
        }
        // Absorption: (A or B) and A -> A
        if (l instanceof Or) {
            if (l.left.toString() === r.toString() || l.right.toString() === r.toString()) return r;
        }

        // Complementarity: A and not(A) -> false
        if (l instanceof Not && l.arg.toString() === r.toString()) return new Num(0);
        if (r instanceof Not && r.arg.toString() === l.toString()) return new Num(0);

        // Distributivity / Factoring for AND of ORs
        // (A or B) and (A or C) -> A or (B and C)
        if (l instanceof Or && r instanceof Or) {
            let common = null;
            let term1 = null;
            let term2 = null;

            if (l.left.toString() === r.left.toString()) { common = l.left; term1 = l.right; term2 = r.right; }
            else if (l.left.toString() === r.right.toString()) { common = l.left; term1 = l.right; term2 = r.left; }
            else if (l.right.toString() === r.left.toString()) { common = l.right; term1 = l.left; term2 = r.right; }
            else if (l.right.toString() === r.right.toString()) { common = l.right; term1 = l.left; term2 = r.left; }

            if (common) {
                const combined = new And(term1, term2).simplify();
                return new Or(common, combined).simplify();
            }
        }

        return new And(l, r);
    }
    evaluateNumeric() { return (this.left.evaluateNumeric() && this.right.evaluateNumeric()) ? 1 : 0; }
    toLatex() { return `${this.left.toLatex()} \\land ${this.right.toLatex()}`; }
}

class Or extends BinaryOp {
    toString() { return `(${this.left} or ${this.right})`; }
    simplify() {
        const l = this.left.simplify();
        const r = this.right.simplify();
        if (l instanceof Num && r instanceof Num) {
            return new Num((l.value || r.value) ? 1 : 0);
        }
        // Short-circuit logic: true || x -> true
        if (l instanceof Num && l.value !== 0) return new Num(1);
        if (r instanceof Num && r.value !== 0) return new Num(1);
        // false || x -> x
        if (l instanceof Num && l.value === 0) return r;
        if (r instanceof Num && r.value === 0) return l;

        // Idempotence: A or A -> A
        if (l.toString() === r.toString()) return l;

        // Absorption: A or (A and B) -> A
        if (r instanceof And) {
            if (r.left.toString() === l.toString() || r.right.toString() === l.toString()) return l;
        }
        // Absorption: (A and B) or A -> A
        if (l instanceof And) {
            if (l.left.toString() === r.toString() || l.right.toString() === r.toString()) return r;
        }

        // Complementarity: A or not(A) -> true
        if (l instanceof Not && l.arg.toString() === r.toString()) return new Num(1);
        if (r instanceof Not && r.arg.toString() === l.toString()) return new Num(1);

        // Distributivity / Factoring
        // (A and B) or (A and C) -> A and (B or C)
        // If (B or C) simplifies to true (e.g. B is not C), then result is A.
        if (l instanceof And && r instanceof And) {
            let common = null;
            let term1 = null;
            let term2 = null;

            if (l.left.toString() === r.left.toString()) { common = l.left; term1 = l.right; term2 = r.right; }
            else if (l.left.toString() === r.right.toString()) { common = l.left; term1 = l.right; term2 = r.left; }
            else if (l.right.toString() === r.left.toString()) { common = l.right; term1 = l.left; term2 = r.right; }
            else if (l.right.toString() === r.right.toString()) { common = l.right; term1 = l.left; term2 = r.left; }

            if (common) {
                const combined = new Or(term1, term2).simplify();
                return new And(common, combined).simplify();
            }
        }

        return new Or(l, r);
    }
    evaluateNumeric() { return (this.left.evaluateNumeric() || this.right.evaluateNumeric()) ? 1 : 0; }
    toLatex() { return `${this.left.toLatex()} \\lor ${this.right.toLatex()}`; }
}

class Xor extends BinaryOp {
    toString() { return `(${this.left} xor ${this.right})`; }
    simplify() {
        const l = this.left.simplify();
        const r = this.right.simplify();
        if (l instanceof Num && r instanceof Num) {
            return new Num((!!l.value !== !!r.value) ? 1 : 0);
        }
        // xor(0, x) -> x, xor(1, x) -> not x
        if (l instanceof Num) {
            if (l.value === 0) return r;
            if (l.value !== 0) return new Not(r).simplify();
        }
        if (r instanceof Num) {
            if (r.value === 0) return l;
            if (r.value !== 0) return new Not(l).simplify();
        }

        if (l.toString() === r.toString()) return new Num(0);
        if (l instanceof Not && l.arg.toString() === r.toString()) return new Num(1);
        if (r instanceof Not && r.arg.toString() === l.toString()) return new Num(1);

        return new Xor(l, r);
    }
    evaluateNumeric() { return (!!this.left.evaluateNumeric() !== !!this.right.evaluateNumeric()) ? 1 : 0; }
    toLatex() { return `${this.left.toLatex()} \\oplus ${this.right.toLatex()}`; }
}

class Implies extends BinaryOp {
    toString() { return `(${this.left} -> ${this.right})`; }
    simplify() {
        const l = this.left.simplify();
        const r = this.right.simplify();
        if (l instanceof Num && r instanceof Num) {
            // !l || r
            const lv = l.value !== 0;
            const rv = r.value !== 0;
            return new Num((!lv || rv) ? 1 : 0);
        }
        // false -> x  => true
        if (l instanceof Num && l.value === 0) return new Num(1);
        // true -> x   => x
        if (l instanceof Num && l.value !== 0) return r;
        // x -> true   => true
        if (r instanceof Num && r.value !== 0) return new Num(1);
        // x -> false  => not x
        if (r instanceof Num && r.value === 0) return new Not(l).simplify();
        // x -> x      => true
        if (l.toString() === r.toString()) return new Num(1);

        return new Implies(l, r);
    }
    evaluateNumeric() {
        const lv = this.left.evaluateNumeric();
        const rv = this.right.evaluateNumeric();
        return (!lv || rv) ? 1 : 0;
    }
    toLatex() { return `${this.left.toLatex()} \\implies ${this.right.toLatex()}`; }
}

class Iff extends BinaryOp {
    toString() { return `(${this.left} <-> ${this.right})`; }
    simplify() {
        const l = this.left.simplify();
        const r = this.right.simplify();
        if (l instanceof Num && r instanceof Num) {
            const lv = !!l.value;
            const rv = !!r.value;
            return new Num((lv === rv) ? 1 : 0);
        }
        // true <-> x  => x
        if (l instanceof Num && l.value !== 0) return r;
        if (r instanceof Num && r.value !== 0) return l;
        // false <-> x => not x
        if (l instanceof Num && l.value === 0) return new Not(r).simplify();
        if (r instanceof Num && r.value === 0) return new Not(l).simplify();
        // x <-> x     => true
        if (l.toString() === r.toString()) return new Num(1);

        return new Iff(l, r);
    }
    evaluateNumeric() {
        const lv = !!this.left.evaluateNumeric();
        const rv = !!this.right.evaluateNumeric();
        return (lv === rv) ? 1 : 0;
    }
    toLatex() { return `${this.left.toLatex()} \\iff ${this.right.toLatex()}`; }
}

class Not extends Expr {
    constructor(arg) {
        super();
        this.arg = arg;
    }
    toString() { return `not(${this.arg})`; }
    simplify() {
        const a = this.arg.simplify();
        if (a instanceof Num) {
            // Check for explicit 0 or non-zero
            return new Num(a.value === 0 ? 1 : 0);
        }
        // Double negation: not(not(x)) -> x
        if (a instanceof Not) {
            return a.arg;
        }

        // De Morgan's Laws
        // not(A or B) -> not A and not B
        if (a instanceof Or) {
            return new And(new Not(a.left), new Not(a.right)).simplify();
        }
        // not(A and B) -> not A or not B
        if (a instanceof And) {
            return new Or(new Not(a.left), new Not(a.right)).simplify();
        }

        return new Not(a);
    }
    evaluateNumeric() { return (!this.arg.evaluateNumeric()) ? 1 : 0; }
    toLatex() { return `\\neg ${this.arg.toLatex()}`; }
    substitute(varName, value) { return new Not(this.arg.substitute(varName, value)); }
}

class Mod extends BinaryOp {
    toString() { return `(${this.left} mod ${this.right})`; }
    simplify() {
        const l = this.left.simplify();
        const r = this.right.simplify();
        if (l instanceof Num && r instanceof Num && r.value !== 0) {
            return new Num(l.value % r.value);
        }
        return new Mod(l, r);
    }
    evaluateNumeric() { return this.left.evaluateNumeric() % this.right.evaluateNumeric(); }
    toLatex() { return `${this.left.toLatex()} \\pmod{${this.right.toLatex()}}`; }
}

class Neq extends BinaryOp {
    toString() { return `${this.left} != ${this.right}`; }
    simplify() {
        const l = this.left.simplify();
        const r = this.right.simplify();
        if (l instanceof Num && r instanceof Num) {
            return new Num((l.value !== r.value) ? 1 : 0);
        }
        return new Neq(l, r);
    }
    evaluateNumeric() { return (this.left.evaluateNumeric() !== this.right.evaluateNumeric()) ? 1 : 0; }
    toLatex() { return `${this.left.toLatex()} \\neq ${this.right.toLatex()}`; }
}

class Lt extends BinaryOp {
    toString() { return `${this.left} < ${this.right}`; }
    simplify() {
        const l = this.left.simplify();
        const r = this.right.simplify();
        if (l instanceof Num && r instanceof Num) {
            return new Num((l.value < r.value) ? 1 : 0);
        }
        return new Lt(l, r);
    }
    evaluateNumeric() { return (this.left.evaluateNumeric() < this.right.evaluateNumeric()) ? 1 : 0; }
    diff(varName) { return new Num(0); }
    toLatex() { return `${this.left.toLatex()} < ${this.right.toLatex()}`; }
    diff(varName) { return new Num(0); } // Treat as constant for differentiation safety
}

class Gt extends BinaryOp {
    toString() { return `${this.left} > ${this.right}`; }
    simplify() {
        const l = this.left.simplify();
        const r = this.right.simplify();
        if (l instanceof Num && r instanceof Num) {
            return new Num((l.value > r.value) ? 1 : 0);
        }
        return new Gt(l, r);
    }
    evaluateNumeric() { return (this.left.evaluateNumeric() > this.right.evaluateNumeric()) ? 1 : 0; }
    diff(varName) { return new Num(0); }
    toLatex() { return `${this.left.toLatex()} > ${this.right.toLatex()}`; }
    diff(varName) { return new Num(0); }
}

class Le extends BinaryOp {
    toString() { return `${this.left} <= ${this.right}`; }
    simplify() {
        const l = this.left.simplify();
        const r = this.right.simplify();
        if (l instanceof Num && r instanceof Num) {
            return new Num((l.value <= r.value) ? 1 : 0);
        }
        return new Le(l, r);
    }
    evaluateNumeric() { return (this.left.evaluateNumeric() <= this.right.evaluateNumeric()) ? 1 : 0; }
    diff(varName) { return new Num(0); }
    toLatex() { return `${this.left.toLatex()} \\leq ${this.right.toLatex()}`; }
    diff(varName) { return new Num(0); }
}

class Ge extends BinaryOp {
    toString() { return `${this.left} >= ${this.right}`; }
    simplify() {
        const l = this.left.simplify();
        const r = this.right.simplify();
        if (l instanceof Num && r instanceof Num) {
            return new Num((l.value >= r.value) ? 1 : 0);
        }
        return new Ge(l, r);
    }
    evaluateNumeric() { return (this.left.evaluateNumeric() >= this.right.evaluateNumeric()) ? 1 : 0; }
    diff(varName) { return new Num(0); }
    toLatex() { return `${this.left.toLatex()} \\geq ${this.right.toLatex()}`; }
    diff(varName) { return new Num(0); }
}

class At extends Expr {
    constructor(obj, index) {
        super();
        this.obj = obj;
        this.index = index;
    }
    toString() { return `${this.obj}[${this.index}]`; }
    simplify() {
        const o = this.obj.simplify();
        const i = this.index.simplify();
        if (o instanceof Vec && i instanceof Num && Number.isInteger(i.value)) {
            if (i.value >= 0 && i.value < o.elements.length) {
                return o.elements[i.value];
            }
            // Index out of bounds? Return symbolic or throw?
            // Xcas might return undef. Let's return At(o, i).
        }
        return new At(o, i);
    }
    evaluateNumeric() {
        const o = this.obj.simplify();
        const i = this.index.evaluateNumeric();
        if (o instanceof Vec && !isNaN(i) && Number.isInteger(i)) {
            if (i >= 0 && i < o.elements.length) {
                return o.elements[i].evaluateNumeric();
            }
        }
        return NaN;
    }
    diff(varName) {
        // Derivative of indexed element?
        // If index is constant, it's diff(obj[i]).
        return new At(this.obj.diff(varName), this.index);
    }
    integrate(varName) {
        return new At(this.obj.integrate(varName), this.index);
    }
    substitute(varName, value) {
        return new At(this.obj.substitute(varName, value), this.index.substitute(varName, value));
    }
    toLatex() {
        return `${this.obj.toLatex()}_{${this.index.toLatex()}}`;
    }
}

class BooleanEq extends BinaryOp {
    toString() { return `${this.left} == ${this.right}`; }
    simplify() {
        const l = this.left.simplify();
        const r = this.right.simplify();
        if (l instanceof Num && r instanceof Num) {
            return new Num((Math.abs(l.value - r.value) < 1e-10) ? 1 : 0);
        }
        if (l.toString() === r.toString()) return new Num(1);
        return new BooleanEq(l, r);
    }
    evaluateNumeric() {
        const l = this.left.evaluateNumeric();
        const r = this.right.evaluateNumeric();
        if (!isNaN(l) && !isNaN(r)) return (Math.abs(l - r) < 1e-10) ? 1 : 0;
        return NaN;
    }
    toLatex() { return `${this.left.toLatex()} == ${this.right.toLatex()}`; }
    substitute(varName, value) { return new BooleanEq(this.left.substitute(varName, value), this.right.substitute(varName, value)); }
}

class If extends Expr {
    constructor(condition, trueBlock, falseBlock) {
        super();
        this.condition = condition;
        this.trueBlock = trueBlock;
        this.falseBlock = falseBlock;
    }
    toString() {
        let s = `if (${this.condition}) { ${this.trueBlock} }`;
        if (this.falseBlock) s += ` else { ${this.falseBlock} }`;
        return s;
    }
    simplify() {
        return new If(this.condition.simplify(), this.trueBlock.simplify(), this.falseBlock ? this.falseBlock.simplify() : null);
    }
    toLatex() {
        let s = `\\text{if } ${this.condition.toLatex()} \\text{ then } \\left\\{ ${this.trueBlock.toLatex()} \\right\\}`;
        if (this.falseBlock) s += ` \\text{ else } \\left\\{ ${this.falseBlock.toLatex()} \\right\\}`;
        return s;
    }
}

class While extends Expr {
    constructor(condition, body) {
        super();
        this.condition = condition;
        this.body = body;
    }
    toString() { return `while (${this.condition}) { ${this.body} }`; }
    simplify() { return new While(this.condition.simplify(), this.body.simplify()); }
    toLatex() { return `\\text{while } ${this.condition.toLatex()} \\text{ do } \\left\\{ ${this.body.toLatex()} \\right\\}`; }
}

class For extends Expr {
    constructor(init, condition, step, body) {
        super();
        this.init = init;
        this.condition = condition;
        this.step = step;
        this.body = body;
    }
    toString() { return `for (${this.init}; ${this.condition}; ${this.step}) { ${this.body} }`; }
    simplify() {
        return new For(
            this.init ? this.init.simplify() : null,
            this.condition ? this.condition.simplify() : null,
            this.step ? this.step.simplify() : null,
            this.body.simplify()
        );
    }
    toLatex() { return `\\text{for } ...`; }
}

class Return extends Expr {
    constructor(value) {
        super();
        this.value = value;
    }
    toString() { return `return ${this.value}`; }
    simplify() { return new Return(this.value.simplify()); }
    toLatex() { return `\\text{return } ${this.value.toLatex()}`; }
}

class Break extends Expr {
    toString() { return `break`; }
    simplify() { return this; }
    toLatex() { return `\\text{break}`; }
}

class Continue extends Expr {
    toString() { return `continue`; }
    simplify() { return this; }
    toLatex() { return `\\text{continue}`; }
}

// --- Numeric Helpers ---

function math_gamma(z) {
    if (z === 0 || (z < 0 && Number.isInteger(z))) return NaN; // Poles
    if (z < 0.5) {
        return Math.PI / (Math.sin(Math.PI * z) * math_gamma(1 - z));
    }
    return Math.exp(math_logGamma(z));
}

function math_logGamma(z) {
    if (z < 0.5) {
        return Math.log(Math.PI / (Math.sin(Math.PI * z) * math_gamma(1 - z)));
    }
    const g = 7;
    const C = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
    let x = C[0];
    for (let i = 1; i < g + 2; i++) x += C[i] / (z - 1 + i);
    const t = z - 1 + g + 0.5;
    return 0.5 * Math.log(2 * Math.PI) + (z - 0.5) * Math.log(t) - t + Math.log(x);
}

function math_beta(a, b) {
    // exp(lgamma(a) + lgamma(b) - lgamma(a+b)) to avoid overflow
    return Math.exp(math_logGamma(a) + math_logGamma(b) - math_logGamma(a + b));
}

function math_zeta(s) {
    if (s === 1) return Infinity;
    if (s < 0) {
        // Reflection: zeta(s) = 2^s * pi^(s-1) * sin(pi*s/2) * gamma(1-s) * zeta(1-s)
        return Math.pow(2, s) * Math.pow(Math.PI, s - 1) * Math.sin(Math.PI * s / 2) * math_gamma(1 - s) * math_zeta(1 - s);
    }

    // Simple series for s > 1
    if (s > 1) {
        let sum = 0;
        let limit = 100000; // Increased precision
        for (let i = 1; i <= limit; i++) {
            sum += Math.pow(i, -s);
        }
        return sum;
    }
    return NaN; // Not implemented for critical strip yet
}

function math_lambertw(x) {
    // Halley's method
    if (x < -1/Math.E) return NaN;
    let w = Math.log(x + 1); // Initial guess
    if (x > 1) w = Math.log(x) - Math.log(Math.log(x));

    for (let i = 0; i < 10; i++) {
        const ew = Math.exp(w);
        const wew = w * ew;
        const wewx = wew - x;
        w = w - wewx / (ew * (w + 1) - (w + 2) * wewx / (2 * w + 2));
    }
    return w;
}

function math_psi(x) {
    // Digamma function approximation
    if (x <= 0) {
        // Reflection: psi(1-x) - psi(x) = pi * cot(pi*x)
        // psi(x) = psi(1-x) - pi * cot(pi*x)
        // Ensure not integer (pole)
        if (Number.isInteger(x)) return NaN; // Pole
        return math_psi(1 - x) - Math.PI / Math.tan(Math.PI * x);
    }

    // Use asymptotic expansion for large x
    // psi(x) ~ ln(x) - 1/2x - 1/12x^2 + 1/120x^4 - 1/252x^6 ...
    if (x < 10) {
        // Shift up
        return math_psi(x + 1) - 1/x;
    }

    let res = Math.log(x) - 1/(2*x);
    const x2 = x * x;
    const x4 = x2 * x2;
    res -= 1/(12*x2);
    res += 1/(120*x4);
    res -= 1/(252*x2*x4);
    return res;
}

function math_polygamma(n, z) {
    // n: integer order (0 = psi)
    // z: argument
    if (n < 0 || !Number.isInteger(n)) return NaN;
    if (n === 0) return math_psi(z);

    // Reflection Formula: (-1)^n * psi_n(1-z) - psi_n(z) = pi * d^n/dz^n (cot(pi*z))
    // This is complex. We stick to z > 0 recursion.
    if (z <= 0) {
        if (Number.isInteger(z)) return NaN; // Pole
        // Use recurrence relation: psi_n(z+1) = psi_n(z) + (-1)^n * n! / z^(n+1)
        // -> psi_n(z) = psi_n(z+1) - (-1)^n * n! / z^(n+1)
        // We shift UP until positive.
        // Wait, for negative non-integers, we can shift up to positive region.
        // Formula: psi_n(z) = psi_n(z+k) - (-1)^n * n! * sum_{j=0}^{k-1} (z+j)^-(n+1)

        // Simple recursive shift
        // (-1)^n * n! is constant
        const factN = (function f(k){ return k<=1?1:k*f(k-1); })(n);
        const term = (n % 2 === 0 ? 1 : -1) * factN * Math.pow(z, -(n+1));
        return math_polygamma(n, z+1) - term;
    }

    // Recurrence for small positive z
    if (z < 15) {
        const factN = (function f(k){ return k<=1?1:k*f(k-1); })(n);
        const term = (n % 2 === 0 ? 1 : -1) * factN * Math.pow(z, -(n+1));
        return math_polygamma(n, z + 1) - term;
    }

    // Asymptotic Expansion for large z
    // psi_n(z) ~ (-1)^(n-1) [ (n-1)! / z^n  +  n! / (2 z^(n+1)) + ... ]
    // Check Trigamma (n=1): 1/z + 1/2z^2 + 1/6z^3 ...
    // Formula: psi_n(z) = (-1)^(n+1) * n! * [ 1/(n*z^n) + 1/(2*z^(n+1)) + sum B_2k / (2k)! * (n+2k-1)!/(n!) * z^-(n+2k) ] ?
    // Let's use the explicit summation form
    // (-1)^(n+1) * [ (n-1)!/z^n + n!/(2z^(n+1)) + sum B_2k * (n+2k-1)!/(2k)! * z^-(n+2k) ]
    // Wait, let's verify Trigamma n=1. (-1)^2 [ 0!/z + 1!/2z^2 + B2 * 2!/2! * z^-3 ] = 1/z + 1/2z^2 + 1/6z^3. Correct.

    // Factorial helper
    const fact = (k) => { let r=1; for(let i=2; i<=k; i++) r*=i; return r; };
    const factN = fact(n);
    const factN_1 = fact(n-1);

    let res = factN_1 / Math.pow(z, n);
    res += factN / (2 * Math.pow(z, n+1));

    // Bernoulli numbers B2, B4, ...
    const B = [1/6, -1/30, 1/42, -1/30, 5/66, -691/2730, 7/6];
    // terms z^-(n+2), z^-(n+4)...

    for(let k=1; k<=B.length; k++) {
        const b2k = B[k-1];
        const exponent = n + 2*k;
        const coeff = b2k * fact(n + 2*k - 1) / fact(2*k); // (n+2k-1)! / (2k)!
        res += coeff * Math.pow(z, -exponent);
    }

    return (n % 2 === 0 ? -1 : 1) * res; // (-1)^(n+1)
}

function math_Si(x) {
    if (x === 0) return 0;
    // Series approximation
    // Si(x) = sum_{k=0} (-1)^k x^(2k+1) / ((2k+1)(2k+1)!)
    let sum = 0;
    let term = x; // k=0: x / (1*1!) = x
    let k = 0;
    let fact = 1; // (2k+1)!
    // x term already added
    sum = x;

    // Iterate
    for(k=1; k<50; k++) {
        // term k: (-1)^k * x^(2k+1) / ((2k+1)*(2k+1)!)
        // prev term (k-1): (-1)^(k-1) * x^(2k-1) / ((2k-1)*(2k-1)!)

        if (Math.abs(x) > 10) {
            // Asymptotic: pi/2 - cos(x)/x - sin(x)/x^2 ... (for x>0)
            const sign = x > 0 ? 1 : -1;
            return sign * Math.PI/2 - Math.cos(x)/x - Math.sin(x)/(x*x);
        }

        // Update factorial from (2k-1)! to (2k+1)!
        // (2k+1)! = (2k+1)(2k)(2k-1)!
        const twoK = 2*k;
        fact *= twoK * (twoK + 1);

        const num = Math.pow(x, 2*k + 1);
        const den = (2*k + 1) * fact;
        const add = (k % 2 === 0 ? 1 : -1) * num / den;

        sum += add;
        if (Math.abs(add) < 1e-15) break;
    }
    return sum;
}

function math_Ci(x) {
    if (x <= 0) return NaN; // undefined for x<=0
    // Ci(x) = gamma + ln(x) + sum_{k=1} (-1)^k x^2k / (2k * (2k)!)
    const EULER = 0.5772156649;
    let sum = 0;
    let fact = 1; // (2k)!

    if (x > 10) {
        // Asymptotic: sin(x)/x - cos(x)/x^2 ...
        return Math.sin(x)/x - Math.cos(x)/(x*x);
    }

    for(let k=1; k<50; k++) {
        const twoK = 2*k;
        fact *= (twoK-1) * twoK;

        const num = Math.pow(x, twoK);
        const den = twoK * fact;
        const add = (k % 2 === 0 ? 1 : -1) * num / den;

        sum += add;
        if (Math.abs(add) < 1e-15) break;
    }
    return EULER + Math.log(x) + sum;
}

function math_Ei(x) {
    if (x === 0) return -Infinity;
    // Ei(x) = gamma + ln|x| + sum_{k=1} x^k / (k * k!)
    // For large positive x, Asymptotic: exp(x)/x * (1 + 1/x + 2!/x^2 + ...)

    if (x > 20) {
        // Asymptotic
        let sum = 1;
        let term = 1;
        for (let k = 1; k < 10; k++) {
            term *= k / x;
            sum += term;
        }
        return Math.exp(x) / x * sum;
    }

    const EULER = 0.5772156649;
    let sum = 0;
    let fact = 1; // k!

    for (let k = 1; k < 100; k++) {
        fact *= k;
        const num = Math.pow(x, k);
        const den = k * fact;
        const add = num / den;
        sum += add;
        if (Math.abs(add) < 1e-15) break;
    }

    return EULER + Math.log(Math.abs(x)) + sum;
}

function math_erf(val) {
    if (isNaN(val)) return NaN;
    const sign = (val >= 0) ? 1 : -1;
    const x = Math.abs(val);
    const a1 =  0.254829592;
    const a2 = -0.284496736;
    const a3 =  1.421413741;
    const a4 = -1.453152027;
    const a5 =  1.061405429;
    const p  =  0.3275911;
    const t = 1.0 / (1.0 + p * x);
    const y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-x * x);
    return sign * y;
}

function math_erfinv(x) {
    // erfinv(x) = invNorm((x+1)/2) / sqrt(2)
    // Domain (-1, 1)
    if (x <= -1 || x >= 1) {
        if (x === 1) return Infinity;
        if (x === -1) return -Infinity;
        return NaN;
    }
    const p = (x + 1.0) / 2.0;
    const invNorm = math_invNormStandard(p);
    return invNorm / Math.SQRT2;
}

function math_invNormStandard(p) {
    // Rational approximation for standard normal quantile
    // Abramowitz and Stegun 26.2.23
    if (p <= 0) return -Infinity;
    if (p >= 1) return Infinity;
    if (Math.abs(p - 0.5) < 1e-15) return 0;

    const c0 = 2.515517;
    const c1 = 0.802853;
    const c2 = 0.010328;
    const d1 = 1.432788;
    const d2 = 0.189269;
    const d3 = 0.001308;

    let t, num, den, xp;

    let q = (p < 0.5) ? p : 1.0 - p;
    t = Math.sqrt(-2.0 * Math.log(q));
    num = c0 + c1 * t + c2 * t * t;
    den = 1.0 + d1 * t + d2 * t * t + d3 * t * t * t;

    xp = t - num / den;

    if (p < 0.5) xp = -xp;
    return xp;
}

function math_fresnelS(x) {
    // Power series for small x
    // S(x) = sum_{n=0} (-1)^n (pi/2)^(2n+1) x^(4n+3) / ((2n+1)!(4n+3))
    if (x === 0) return 0;
    if (x < 0) return -math_fresnelS(-x);

    // For large x, asymptotic expansion or limit 0.5
    if (x > 10) return 0.5;

    let sum = 0;
    let fact = 1;
    const pi2 = Math.PI / 2;
    for (let n = 0; n < 20; n++) {
        // (2n+1)!
        if (n > 0) fact *= (2*n) * (2*n + 1);
        const num = Math.pow(-1, n) * Math.pow(pi2, 2*n + 1) * Math.pow(x, 4*n + 3);
        const den = fact * (4*n + 3);
        sum += num / den;
    }
    return sum;
}

function math_fresnelC(x) {
    // C(x) = sum_{n=0} (-1)^n (pi/2)^(2n) x^(4n+1) / ((2n)!(4n+1))
    if (x === 0) return 0;
    if (x < 0) return -math_fresnelC(-x);
    if (x > 10) return 0.5;

    let sum = 0;
    let fact = 1;
    const pi2 = Math.PI / 2;
    for (let n = 0; n < 20; n++) {
        if (n > 0) fact *= (2*n - 1) * (2*n);
        const num = Math.pow(-1, n) * Math.pow(pi2, 2*n) * Math.pow(x, 4*n + 1);
        const den = fact * (4*n + 1);
        sum += num / den;
    }
    return sum;
}

function math_EllipticK(k) {
    if (Math.abs(k) >= 1) return Infinity;
    if (k === 0) return Math.PI / 2;
    let a = 1;
    let b = Math.sqrt(1 - k*k);
    for(let i=0; i<100; i++) {
        let an = (a + b) / 2;
        let bn = Math.sqrt(a * b);
        if (Math.abs(a - b) < 1e-15) { a = an; break; }
        a = an;
        b = bn;
    }
    return Math.PI / (2 * a);
}

function math_EllipticE(k) {
    if (Math.abs(k) > 1) return NaN;
    if (Math.abs(k) === 1) return 1;
    if (k === 0) return Math.PI / 2;

    let a = 1;
    let b = Math.sqrt(1 - k*k);
    let sum = 0.5 * k * k; // n=0 term
    let p = 1; // 2^0

    // First step of EllipticK to get K value at end?
    // We can run AGM loop and sum simultaneously

    for(let i=0; i<100; i++) {
        let an = (a + b) / 2;
        let bn = Math.sqrt(a * b);

        let c_sq = Math.pow((a - b)/2, 2);
        sum += p * c_sq;

        p *= 2;

        if (Math.abs(a - b) < 1e-15) { a = an; break; }
        a = an;
        b = bn;
    }

    const K = Math.PI / (2 * a);
    return K * (1 - sum);
}

function math_besselJ(v, x) {
    if (x === 0) return (v === 0) ? 1 : 0;
    // Simple series for small x
    if (Math.abs(x) < 10) {
        let sum = 0;
        for(let k=0; k<50; k++) {
            const num = Math.pow(-1, k) * Math.pow(x/2, 2*k + v);
            const den = math_gamma(k + 1) * math_gamma(k + v + 1); // factorial(k)*gamma
            sum += num / den;
            if (Math.abs(num/den) < 1e-15 * Math.abs(sum)) break;
        }
        return sum;
    }
    // Asymptotic for large x
    return Math.sqrt(2 / (Math.PI * x)) * Math.cos(x - v * Math.PI / 2 - Math.PI / 4);
}

function math_besselY(v, x) {
    if (Number.isInteger(v)) {
        const eps = 1e-10;
        return (math_besselJ(v + eps, x) * Math.cos((v + eps) * Math.PI) - math_besselJ(-(v + eps), x)) / Math.sin((v + eps) * Math.PI);
    }
    return (math_besselJ(v, x) * Math.cos(v * Math.PI) - math_besselJ(-v, x)) / Math.sin(v * Math.PI);
}

function math_besselI(v, x) {
    if (Math.abs(x) < 20) {
        let sum = 0;
        for(let k=0; k<50; k++) {
            const num = Math.pow(x/2, 2*k + v);
            const den = math_gamma(k + 1) * math_gamma(k + v + 1);
            sum += num / den;
            if (Math.abs(num/den) < 1e-15 * Math.abs(sum)) break;
        }
        return sum;
    }
    return Math.exp(x) / Math.sqrt(2 * Math.PI * x);
}

function math_besselK(v, x) {
    if (Number.isInteger(v)) {
        const eps = 1e-10;
        return Math.PI / 2 * (math_besselI(-(v+eps), x) - math_besselI(v+eps, x)) / Math.sin((v+eps) * Math.PI);
    }
    return Math.PI / 2 * (math_besselI(-v, x) - math_besselI(v, x)) / Math.sin(v * Math.PI);
}

function math_hyp2f1(a, b, c, z) {
    if (Math.abs(z) >= 1) return NaN;
    let sum = 1;
    let term = 1;
    for(let n=0; n<100; n++) {
        term *= (a + n) * (b + n) / ((c + n) * (n + 1)) * z;
        sum += term;
        if (Math.abs(term) < 1e-15 * Math.abs(sum)) break;
    }
    return sum;
}

function getBernoulliExpr(n) {
    const B = {
        0: [1, 1], 1: [-1, 2], 2: [1, 6], 4: [-1, 30], 6: [1, 42],
        8: [-1, 30], 10: [5, 66], 12: [-691, 2730], 14: [7, 6], 16: [-3617, 510],
        18: [43867, 798], 20: [-174611, 330]
    };
    if (B[n]) {
        const val = B[n];
        if (val[1] === 1) return new Num(val[0]);
        return new Div(new Num(val[0]), new Num(val[1]));
    }
    // Fallback for n > 20: return symbolic
    return new Call('bernoulli', [new Num(n)]);
}

function math_intSqrtSimplify(n) {
    if (n < 0) return null;
    if (n === 0) return { coeff: 1, radical: 0 };
    if (n === 1) return { coeff: 1, radical: 1 };

    let coeff = 1;
    let radical = n;

    // Check 2
    while (radical % 4 === 0) {
        coeff *= 2;
        radical /= 4;
    }

    // Check odd factors
    let limit = Math.floor(Math.sqrt(radical));
    for (let k = 3; k <= limit; k += 2) {
        let k2 = k*k;
        if (k2 > radical) break;
        while (radical % k2 === 0) {
            coeff *= k;
            radical /= k2;
        }
    }
    return { coeff, radical };
}

// Export classes for Global/CommonJS environments
(function() {
    const exports = {
        Expr, Num, Sym, BinaryOp, Add, Sub, Mul, Div, Pow, Call, Assignment, Eq, Vec, FunctionDef, Block, toExpr,
        And, Or, Xor, Implies, Iff, Not, Mod, Neq, Lt, Gt, Le, Ge, At, BooleanEq,
        If, While, For, Return, Break, Continue
    };
    if (typeof globalThis !== 'undefined') {
        Object.assign(globalThis, exports);
    }
    if (typeof module !== 'undefined' && module.exports) {
        module.exports = exports;
    }
})();
