
const TOKEN_NUMBER = 'NUMBER';
const TOKEN_IDENTIFIER = 'IDENTIFIER';
const TOKEN_PLUS = 'PLUS';
const TOKEN_MINUS = 'MINUS';
const TOKEN_STAR = 'STAR';
const TOKEN_SLASH = 'SLASH';
const TOKEN_CARET = 'CARET';
const TOKEN_LPAREN = 'LPAREN';
const TOKEN_RPAREN = 'RPAREN';
const TOKEN_LBRACKET = 'LBRACKET';
const TOKEN_RBRACKET = 'RBRACKET';
const TOKEN_COMMA = 'COMMA';
const TOKEN_ASSIGN = 'ASSIGN';
const TOKEN_EQ = 'EQ';
const TOKEN_NEQ = 'NEQ';
const TOKEN_LT = 'LT';
const TOKEN_GT = 'GT';
const TOKEN_LE = 'LE';
const TOKEN_GE = 'GE';
const TOKEN_AND = 'AND';
const TOKEN_OR = 'OR';
const TOKEN_NOT = 'NOT';
const TOKEN_XOR = 'XOR';
const TOKEN_MOD = 'MOD';
const TOKEN_SEMI = 'SEMI';
const TOKEN_EOF = 'EOF';

class Token {
    constructor(type, value) {
        this.type = type;
        this.value = value;
    }
}

class Lexer {
    constructor(text) {
        this.text = text;
        this.pos = 0;
        this.currentChar = this.text.length > 0 ? this.text[0] : null;
    }

    advance() {
        this.pos++;
        if (this.pos < this.text.length) {
            this.currentChar = this.text[this.pos];
        } else {
            this.currentChar = null;
        }
    }

    peek() {
        const peekPos = this.pos + 1;
        return peekPos < this.text.length ? this.text[peekPos] : null;
    }

    skipWhitespace() {
        while (this.currentChar !== null && /\s/.test(this.currentChar)) {
            this.advance();
        }
    }

    skipComment() {
        while (this.currentChar !== null && this.currentChar !== '\n') {
            this.advance();
        }
    }

    number() {
        let result = '';
        let dotCount = 0;
        while (this.currentChar !== null && (/\d/.test(this.currentChar) || this.currentChar === '.')) {
            if (this.currentChar === '.') {
                dotCount++;
                if (dotCount > 1) {
                    throw new Error("Invalid number format: multiple decimal points");
                }
            }
            result += this.currentChar;
            this.advance();
        }

        // Scientific notation: 1e5, 1.2e-3
        // Lookahead to ensure it is actually scientific notation and not an identifier starting with e
        if (this.currentChar === 'e' || this.currentChar === 'E') {
            const peek1 = this.peek();
            const peek2 = (peek1 === '+' || peek1 === '-') ? this.text[this.pos + 2] : null;

            // Condition: 'e' followed by digit OR 'e' followed by +/- and then digit
            const isSci = (peek1 !== null && /\d/.test(peek1)) ||
                          ((peek1 === '+' || peek1 === '-') && peek2 !== undefined && /\d/.test(peek2));

            if (isSci) {
                result += this.currentChar;
                this.advance();
                if (this.currentChar === '+' || this.currentChar === '-') {
                    result += this.currentChar;
                    this.advance();
                }
                while (this.currentChar !== null && /\d/.test(this.currentChar)) {
                    result += this.currentChar;
                    this.advance();
                }
            }
        }

        return parseFloat(result);
    }

    identifier() {
        let result = '';
        while (this.currentChar !== null && (/[a-zA-Z0-9_]/.test(this.currentChar))) {
            result += this.currentChar;
            this.advance();
        }
        return result;
    }

    getNextToken() {
        while (this.currentChar !== null) {
            if (/\s/.test(this.currentChar)) {
                this.skipWhitespace();
                continue;
            }
            if (/\d/.test(this.currentChar) || this.currentChar === '.') {
                return new Token(TOKEN_NUMBER, this.number());
            }
            if (/[a-zA-Z]/.test(this.currentChar)) {
                const id = this.identifier();
                const lowerId = id.toLowerCase();
                if (lowerId === 'and') return new Token(TOKEN_AND, id);
                if (lowerId === 'or') return new Token(TOKEN_OR, id);
                if (lowerId === 'not') return new Token(TOKEN_NOT, id);
                if (lowerId === 'xor') return new Token(TOKEN_XOR, id);
                if (lowerId === 'mod') return new Token(TOKEN_MOD, id);
                return new Token(TOKEN_IDENTIFIER, id);
            }
            if (this.currentChar === ':') {
                if (this.peek() === '=') {
                    this.advance();
                    this.advance();
                    return new Token(TOKEN_ASSIGN, ':=');
                }
                throw new Error("Invalid character: :");
            }
            if (this.currentChar === '!') {
                if (this.peek() === '=') {
                    this.advance();
                    this.advance();
                    return new Token(TOKEN_NEQ, '!=');
                }
                // Factorial is handled in Parser, but here we might want to return a token for it?
                // The current Parser doesn't seem to explicitly handle '!' factorial token in lexer except maybe as char?
                // Let's check Parser.factor/atom.
                // Ah, the lexer threw Error for '!'.
                // If it's just '!', return something or error?
                // Xcas uses '!' for factorial.
                // Let's assume for now '!' is factorial if not followed by '='.
                // But wait, "Invalid character: !" was thrown before.
                // So we should return a token for '!'. Let's call it TOKEN_BANG.
                // But I didn't define TOKEN_BANG above. Let's stick to Xcas syntax: `!=` is NEQ. `!` is factorial.
                this.advance();
                return new Token('BANG', '!');
            }
            if (this.currentChar === '<') {
                if (this.peek() === '=') {
                    this.advance();
                    this.advance();
                    return new Token(TOKEN_LE, '<=');
                }
                this.advance();
                return new Token(TOKEN_LT, '<');
            }
            if (this.currentChar === '>') {
                if (this.peek() === '=') {
                    this.advance();
                    this.advance();
                    return new Token(TOKEN_GE, '>=');
                }
                this.advance();
                return new Token(TOKEN_GT, '>');
            }
            if (this.currentChar === '=') {
                if (this.peek() === '=') {
                    // Support == as equality too (Xcas style boolean eq, or standard)
                    this.advance();
                    this.advance();
                    return new Token(TOKEN_EQ, '==');
                }
                this.advance();
                return new Token(TOKEN_EQ, '=');
            }
            if (this.currentChar === '+') { this.advance(); return new Token(TOKEN_PLUS, '+'); }
            if (this.currentChar === '-') { this.advance(); return new Token(TOKEN_MINUS, '-'); }
            if (this.currentChar === '*') { this.advance(); return new Token(TOKEN_STAR, '*'); }
            if (this.currentChar === '/') {
                if (this.peek() === '/') {
                    this.skipComment();
                    this.skipWhitespace();
                    continue;
                }
                this.advance();
                return new Token(TOKEN_SLASH, '/');
            }
            if (this.currentChar === '%') { this.advance(); return new Token(TOKEN_MOD, '%'); }
            if (this.currentChar === '^') { this.advance(); return new Token(TOKEN_CARET, '^'); }
            if (this.currentChar === '&') {
                if (this.peek() === '&') {
                    this.advance(); this.advance();
                    return new Token(TOKEN_AND, '&&');
                }
                throw new Error("Invalid character: &");
            }
            if (this.currentChar === '|') {
                if (this.peek() === '|') {
                    this.advance(); this.advance();
                    return new Token(TOKEN_OR, '||');
                }
                throw new Error("Invalid character: |");
            }
            if (this.currentChar === '(') { this.advance(); return new Token(TOKEN_LPAREN, '('); }
            if (this.currentChar === ')') { this.advance(); return new Token(TOKEN_RPAREN, ')'); }
            if (this.currentChar === '[') { this.advance(); return new Token(TOKEN_LBRACKET, '['); }
            if (this.currentChar === ']') { this.advance(); return new Token(TOKEN_RBRACKET, ']'); }
            if (this.currentChar === ',') { this.advance(); return new Token(TOKEN_COMMA, ','); }
            if (this.currentChar === ';') { this.advance(); return new Token(TOKEN_SEMI, ';'); }

            throw new Error(`Invalid character: ${this.currentChar}`);
        }
        return new Token(TOKEN_EOF, null);
    }
}

class Parser {
    constructor(lexer) {
        this.lexer = lexer;
        this.currentToken = this.lexer.getNextToken();
    }

    error() {
        throw new Error("Invalid syntax");
    }

    eat(tokenType) {
        if (this.currentToken.type === tokenType) {
            this.currentToken = this.lexer.getNextToken();
        } else {
            this.error();
        }
    }

    factor() {
        let node = this.atom();
        // Check for postfix indexing: A[0]
        while (this.currentToken.type === TOKEN_LBRACKET) {
             // We need to differentiate A[0] (index) from A [0] (implicit mul).
             // However, implicit mul usually requires spaces or specific context.
             // Standard Xcas syntax A[0] means index.
             // If I have `x [1]`, is it x * [1] or x[1]?
             // In Xcas, `x[1]` is index. `x*[1]` is mul.
             // If we parse `[` as index here, `x [1]` will be index.
             // If we want `x [1]` to be mul, user must use `*`.
             // But my implicitMul logic allows `x [1]`.
             // If I consume `[` here, implicitMul won't see it.
             // Is there any case where `Expr [...]` should be multiplication?
             // Only if Expr is a number? `2[1]` -> `2 * [1]`.
             // If Expr is a vector? `[1][0]` -> `1`.
             // If Expr is identifier? `A[0]` -> Index.
             // If Expr is `(a+b)[0]` -> Index.

             // So, basically `[` immediately following an atom (without space?)
             // Lexer doesn't tell us about space.
             // But typical parser logic: Postfix binds tighter than implicit mul.
             // So `A[0]` is index. `2[0]` is index (undefined for number, but syntax wise).
             // If user wants `2 * [0]`, they should type `2 * [0]` or `2 [0]`?
             // If `2[0]` parses as `At(2, 0)`, then evaluation fails or returns something.
             // If we want `2 [0]` to be mul, we have a problem without whitespace info.

             // However, `At` node can decide. If `obj` is Num, maybe `At` acts as mul?
             // No, `At` is index.
             // Let's assume `[` is ALWAYS index if it follows a factor directly.
             // Users usually write `2*[1]` or `2 [1]`.
             // If I consume it here, `implicitMul` loop won't see it.

             // Wait, `implicitMul` logic:
             // `node = this.unary()`
             // `while (isImplicitMulStart)`...

             // `unary` calls `power`, `power` calls `factor`.
             // If `factor` consumes `[`, `unary` returns `At(...)`.
             // Then `implicitMul` sees next token.
             // If input is `A[0]`, `factor` eats `A` then `[0]`. Returns `At`. Next is EOF.
             // If input is `A [0]`. `factor` eats `A`. `[0]` is eaten. Returns `At`.

             // This breaks `x [1]` as implicit mul.
             // But `x(1)` is Call, not mul.
             // So `x[1]` being Index is consistent.

             // What about `2 [1]`? `2(1)` is implicit mul in some systems, but my parser handles `2(1)`?
             // `factor` handles `LPAREN` after `NUMBER`?
             // `factor` code:
             // if NUMBER -> return Num.
             // It does NOT check for `(` after number.
             // So `2(1)` -> `Num(2)`. `implicitMul` sees `(`. `isImplicitMulStart` is true.
             // So `2(1)` becomes `Mul(2, 1)`.

             // But `factor` handles `IDENTIFIER` then `(`.
             // So `x(1)` is `Call`.

             // If I put `[` handling in `atom` (or `factor` bottom), I need to be careful.
             // If `token` is `NUMBER`, `factor` returns `Num`. It does NOT look ahead.
             // So `2[1]` -> `Num(2)`. `unary` returns `Num(2)`.
             // `implicitMul` sees `[`. `isImplicitMulStart` is true.
             // So `2[1]` becomes `Mul(2, Vec([1]))`. Correct for `2 [1]`.

             // But `IDENTIFIER` block consumes `(` for Call.
             // It explicitly checks `LPAREN`.
             // I should add `LBRACKET` check there too.

             // BUT, `(a+b)[0]`?
             // `factor` handles `(` ... `)`. Returns `Add`.
             // If I simply add logic to `IDENTIFIER` block, `(a+b)[0]` won't work.

             // I should restructure `factor` to:
             // `node = primary()`.
             // `while (LBRACKET)` -> `node = At(node, ...)`
             // Then return node.

             // Let's rename `factor` logic (the long if/else) to `atom` or `primary`.
             // And make `factor` wrap it.

             // Wait, `factor` currently handles `sin^2(x)` logic which is weirdly tied to Identifier.
             // Let's refactor carefully.

             this.eat(TOKEN_LBRACKET);
             const index = this.statement(); // Allow expression index
             this.eat(TOKEN_RBRACKET);
             node = new At(node, index);
        }
        return node;
    }

    atom() {
        const token = this.currentToken;
        if (token.type === TOKEN_NUMBER) {
            this.eat(TOKEN_NUMBER);
            return new Num(token.value);
        } else if (token.type === TOKEN_IDENTIFIER) {
            const name = token.value;
            this.eat(TOKEN_IDENTIFIER);

            // Handle sin^2(x)
            if (this.currentToken.type === TOKEN_CARET) {
                // Check if it's a known function that permits exponent syntax
                const trigFunctions = ['sin', 'cos', 'tan', 'sec', 'csc', 'cot', 'sinh', 'cosh', 'tanh'];
                if (trigFunctions.includes(name)) {
                    this.eat(TOKEN_CARET);
                    const exponent = this.unary(); // Use unary to handle sin^-1(x)

                    if (this.currentToken.type === TOKEN_LPAREN) {
                        this.eat(TOKEN_LPAREN);
                        const arg = this.statement();
                        this.eat(TOKEN_RPAREN);
                        return new Pow(new Call(name, [arg]), exponent);
                    } else if (this.isImplicitMulStart(this.currentToken)) {
                         const arg = this.term();
                         return new Pow(new Call(name, [arg]), exponent);
                    }
                } else {
                     // For non-trig functions, `f^2(x)` is ambiguous. `(f^2)(x)`? or `f(x)^2`?
                     // Usually treated as symbol power if f is var.
                     // But we already ate identifier. We need to return Pow(Sym(name), exp)
                     // If followed by LPAREN, it becomes implicit mul? `x^2(y)` -> `x^2 * y`?
                     // Let's fallback to standard power logic if not special trig syntax.
                     // But we are inside factor(). We processed ID.
                     // We need to return Sym(name) and let `power()` handle the caret.
                     // BUT we already ate the caret inside this if block if we matched trig.
                     // If we are here, we matched caret.
                     // If not trig, we construct Pow(Sym, exp).

                     // Wait, `power()` calls `factor()`. `factor()` consumes ID. `power()` consumes caret.
                     // So we should NOT consume caret here unless we are sure it is `sin^2(x)`.
                     // If we are here, we ate ID.
                     // If we see caret, we return Sym(name) and let `power()` handle it?
                     // NO, `power` calls `factor` then checks caret.
                     // If `factor` consumes caret, `power` won't see it.
                     // So `sin^2(x)` logic must be inside `factor` OR `power`.

                     // Actually `sin^n(x)` is tricky because `power` binds tighter than `call`?
                     // If `sin` is handled as `Call` in `factor`, we never return to `power` with just `Sym`.
                     // The loop in `factor` for `(` handles explicit calls.
                }
            }

            if (this.currentToken.type === TOKEN_LPAREN) {
                this.eat(TOKEN_LPAREN);
                const args = [];
                if (this.currentToken.type !== TOKEN_RPAREN) {
                    args.push(this.statement());
                    while (this.currentToken.type === TOKEN_COMMA) {
                        this.eat(TOKEN_COMMA);
                        args.push(this.statement());
                    }
                }
                this.eat(TOKEN_RPAREN);
                return new Call(name, args);
            } else {
                // Check for implicit function call: sin x, log 10
                const knownFunctions = [
                    'sin', 'cos', 'tan', 'asin', 'acos', 'atan',
                    'sinh', 'cosh', 'tanh', 'asinh', 'acosh', 'atanh',
                    'sec', 'csc', 'cot', 'asec', 'acsc', 'acot',
                    'sqrt', 'log', 'ln', 'exp', 'abs', 'sign',
                    'floor', 'ceil', 'round', 'fact', 'factorial', 'gamma'
                ];
                if (knownFunctions.includes(name)) {
                     // If next is a factor start (implicit arg)
                     if (this.isImplicitMulStart(this.currentToken)) {
                         const arg = this.term();
                         return new Call(name, [arg]);
                     }
                }

                return new Sym(name);
            }
        } else if (token.type === TOKEN_LPAREN) {
            this.eat(TOKEN_LPAREN);
            const node = this.statement();
            this.eat(TOKEN_RPAREN);
            return node;
        } else if (token.type === TOKEN_LBRACKET) {
            this.eat(TOKEN_LBRACKET);
            const elements = [];
            if (this.currentToken.type !== TOKEN_RBRACKET) {
                elements.push(this.statement());
                while (this.currentToken.type === TOKEN_COMMA) {
                    this.eat(TOKEN_COMMA);
                    elements.push(this.statement());
                }
            }
            this.eat(TOKEN_RBRACKET);
            return new Vec(elements);
        }
        this.error();
    }

    unary() {
        if (this.currentToken.type === TOKEN_PLUS) {
            this.eat(TOKEN_PLUS);
            return this.unary();
        } else if (this.currentToken.type === TOKEN_MINUS) {
            this.eat(TOKEN_MINUS);
            return new Mul(new Num(-1), this.unary());
        } else if (this.currentToken.type === TOKEN_NOT) {
            this.eat(TOKEN_NOT);
            return new Not(this.unary());
        }
        return this.power();
    }

    power() {
        let node = this.factor();
        // Handle factorial '!' here if we had it as postfix
        // But wait, the lexer returns 'BANG' for '!'.
        if (this.currentToken.type === 'BANG') {
            this.eat('BANG');
            node = new Call('fact', [node]);
        }

        if (this.currentToken.type === TOKEN_CARET) {
            this.eat(TOKEN_CARET);
            node = new Pow(node, this.unary());
        }
        return node;
    }

    term() {
        let node = this.implicitMul();
        while (this.currentToken.type === TOKEN_STAR ||
               this.currentToken.type === TOKEN_SLASH ||
               this.currentToken.type === TOKEN_MOD) {

            if (this.currentToken.type === TOKEN_STAR) {
                this.eat(TOKEN_STAR);
                node = new Mul(node, this.implicitMul());
            } else if (this.currentToken.type === TOKEN_SLASH) {
                this.eat(TOKEN_SLASH);
                node = new Div(node, this.implicitMul());
            } else if (this.currentToken.type === TOKEN_MOD) {
                this.eat(TOKEN_MOD);
                node = new Mod(node, this.implicitMul());
            }
        }
        return node;
    }

    implicitMul() {
        let node = this.unary();
        // Check for implicit multiplication start
        while (this.isImplicitMulStart(this.currentToken)) {
             // Implicit multiplication has higher precedence than * /
             // so it binds tight here.
             const right = this.unary();
             node = new Mul(node, right);
        }
        return node;
    }

    isImplicitMulStart(token) {
        // Tokens that can start a power/factor:
        // NUMBER, IDENTIFIER, LPAREN, LBRACKET
        // (MINUS/PLUS are handled in expr, so implicit mul like `x -y` isn't standard, usually `x * -y` needs parens or explicit `*`)
        // However, `2x` works. `x y` works. `x(y)` is parsed in factor, so if we are here, we saw a factor and next is something else.
        // e.g. `x` (Sym) then `y` (Sym).
        return token.type === TOKEN_NUMBER ||
               token.type === TOKEN_IDENTIFIER ||
               token.type === TOKEN_LPAREN ||
               token.type === TOKEN_LBRACKET ||
               token.type === TOKEN_NOT; // not x -> mul? No.
               // If we have `x not y`, that's invalid syntax.
               // So NOT is NOT a start of implicit mul if it's an operator.
               // But `x not(y)`?
               // Let's exclude NOT from implicit mul start for now.
    }

    arithExpr() {
        let node = this.term();
        while (this.currentToken.type === TOKEN_PLUS || this.currentToken.type === TOKEN_MINUS) {
            const token = this.currentToken;
            if (token.type === TOKEN_PLUS) {
                this.eat(TOKEN_PLUS);
                node = new Add(node, this.term());
            } else if (token.type === TOKEN_MINUS) {
                this.eat(TOKEN_MINUS);
                node = new Sub(node, this.term());
            }
        }
        return node;
    }

    compExpr() {
        let node = this.arithExpr();
        if (this.currentToken.type === TOKEN_EQ) {
            this.eat(TOKEN_EQ);
            node = new Eq(node, this.arithExpr());
        } else if (this.currentToken.type === TOKEN_NEQ) {
            this.eat(TOKEN_NEQ);
            node = new Neq(node, this.arithExpr());
        } else if (this.currentToken.type === TOKEN_LT) {
            this.eat(TOKEN_LT);
            node = new Lt(node, this.arithExpr());
        } else if (this.currentToken.type === TOKEN_GT) {
            this.eat(TOKEN_GT);
            node = new Gt(node, this.arithExpr());
        } else if (this.currentToken.type === TOKEN_LE) {
            this.eat(TOKEN_LE);
            node = new Le(node, this.arithExpr());
        } else if (this.currentToken.type === TOKEN_GE) {
            this.eat(TOKEN_GE);
            node = new Ge(node, this.arithExpr());
        }
        return node;
    }

    andExpr() {
        let node = this.compExpr();
        while (this.currentToken.type === TOKEN_AND) {
            this.eat(TOKEN_AND);
            node = new And(node, this.compExpr());
        }
        return node;
    }

    logicExpr() {
        let node = this.andExpr();
        while (this.currentToken.type === TOKEN_OR || this.currentToken.type === TOKEN_XOR) {
            const token = this.currentToken;
            if (token.type === TOKEN_OR) {
                this.eat(TOKEN_OR);
                node = new Or(node, this.andExpr());
            } else if (token.type === TOKEN_XOR) {
                this.eat(TOKEN_XOR);
                node = new Xor(node, this.andExpr());
            }
        }
        return node;
    }

    equation() {
        // In Xcas, assignment is :=, equation is =.
        // We already have compExpr handling = (EQ).
        // However, assignment logic is in statement().
        // So equation() here really means the top-level expression.
        return this.logicExpr();
    }

    statement() {
        // Look ahead hack for Assignment: ID := ... or ID(args) := ...
        if (this.currentToken.type === TOKEN_IDENTIFIER) {
            // Save state
            const savedPos = this.lexer.pos;
            const savedChar = this.lexer.currentChar;
            const savedToken = this.currentToken;

            const name = this.currentToken.value;
            this.eat(TOKEN_IDENTIFIER);

            // Check for simple assignment: ID :=
            if (this.currentToken.type === TOKEN_ASSIGN) {
                this.eat(TOKEN_ASSIGN);
                const value = this.statement();
                return new Assignment(new Sym(name), value);
            }

            // Check for function definition: ID(args) :=
            if (this.currentToken.type === TOKEN_LPAREN) {
                this.eat(TOKEN_LPAREN);
                const args = [];
                let isFuncDef = true;

                // Parse args list, but must be identifiers for definition
                if (this.currentToken.type === TOKEN_IDENTIFIER) {
                    args.push(this.currentToken.value);
                    this.eat(TOKEN_IDENTIFIER);
                    while (this.currentToken.type === TOKEN_COMMA) {
                        this.eat(TOKEN_COMMA);
                        if (this.currentToken.type === TOKEN_IDENTIFIER) {
                            args.push(this.currentToken.value);
                            this.eat(TOKEN_IDENTIFIER);
                        } else {
                            isFuncDef = false;
                            break;
                        }
                    }
                } else if (this.currentToken.type !== TOKEN_RPAREN) {
                    isFuncDef = false;
                }

                if (isFuncDef && this.currentToken.type === TOKEN_RPAREN) {
                    this.eat(TOKEN_RPAREN);
                    if (this.currentToken.type === TOKEN_ASSIGN) {
                        this.eat(TOKEN_ASSIGN);
                        const body = this.statement();
                        return new FunctionDef(name, args, body);
                    }
                }
            }

            // Restore state if not assignment or function def
            this.lexer.pos = savedPos;
            this.lexer.currentChar = savedChar;
            this.currentToken = savedToken;
        }
        return this.equation();
    }

    parse() {
        const statements = [];
        statements.push(this.statement());
        while (this.currentToken.type === TOKEN_SEMI) {
            this.eat(TOKEN_SEMI);
            if (this.currentToken.type !== TOKEN_EOF) {
                statements.push(this.statement());
            }
        }
        if (statements.length === 1) {
            return statements[0];
        }
        return new Block(statements);
    }
}

// Export classes for Global/CommonJS environments
(function() {
    const exports = {
        Token, Lexer, Parser
    };
    if (typeof globalThis !== 'undefined') {
        Object.assign(globalThis, exports);
    }
    if (typeof module !== 'undefined' && module.exports) {
        module.exports = exports;
    }
})();
