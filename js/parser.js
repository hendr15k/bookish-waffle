
const TOKEN_NUMBER = 'NUMBER';
const TOKEN_IDENTIFIER = 'IDENTIFIER';
const TOKEN_PLUS = 'PLUS';
const TOKEN_MINUS = 'MINUS';
const TOKEN_STAR = 'STAR';
const TOKEN_SLASH = 'SLASH';
const TOKEN_CARET = 'CARET';
const TOKEN_UNDERSCORE = 'UNDERSCORE';
const TOKEN_LPAREN = 'LPAREN';
const TOKEN_RPAREN = 'RPAREN';
const TOKEN_LBRACKET = 'LBRACKET';
const TOKEN_RBRACKET = 'RBRACKET';
const TOKEN_COMMA = 'COMMA';
const TOKEN_ASSIGN = 'ASSIGN';
const TOKEN_EQ = 'EQ';
const TOKEN_BOOL_EQ = 'BOOL_EQ';
const TOKEN_NEQ = 'NEQ';
const TOKEN_LT = 'LT';
const TOKEN_GT = 'GT';
const TOKEN_LE = 'LE';
const TOKEN_GE = 'GE';
const TOKEN_AND = 'AND';
const TOKEN_OR = 'OR';
const TOKEN_NOT = 'NOT';
const TOKEN_XOR = 'XOR';
const TOKEN_IMPLIES = 'IMPLIES';
const TOKEN_IFF = 'IFF';
const TOKEN_MOD = 'MOD';
const TOKEN_RANGE = 'RANGE';
const TOKEN_SEMI = 'SEMI';
const TOKEN_EOF = 'EOF';
const TOKEN_LBRACE = 'LBRACE';
const TOKEN_RBRACE = 'RBRACE';
const TOKEN_PIPE = 'PIPE';
const TOKEN_LATEX_FRAC = 'LATEX_FRAC';
const TOKEN_LATEX_SQRT = 'LATEX_SQRT';
const TOKEN_LATEX_SUM = 'LATEX_SUM';
const TOKEN_LATEX_PROD = 'LATEX_PROD';
const TOKEN_LATEX_INT = 'LATEX_INT';
const TOKEN_LATEX_LIMIT = 'LATEX_LIMIT';

// Control Flow Tokens
const TOKEN_IF = 'IF';
const TOKEN_ELSE = 'ELSE';
const TOKEN_WHILE = 'WHILE';
const TOKEN_FOR = 'FOR';
const TOKEN_RETURN = 'RETURN';
const TOKEN_BREAK = 'BREAK';
const TOKEN_CONTINUE = 'CONTINUE';

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
                if (this.peek() === '.') {
                    break;
                }
                dotCount++;
                if (dotCount > 1) {
                    throw new Error("Invalid number format: multiple decimal points");
                }
            }
            result += this.currentChar;
            this.advance();
        }

        if (this.currentChar === 'e' || this.currentChar === 'E') {
            const peek1 = this.peek();
            const peek2 = (peek1 === '+' || peek1 === '-') ? this.text[this.pos + 2] : null;
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
        // REMOVED '_' from identifier to support x_1 as x subscript 1
        while (this.currentChar !== null && (/[a-zA-Z0-9]/.test(this.currentChar))) {
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

            if (this.currentChar === '\\') {
                this.advance();
                if (this.currentChar === '{') { this.advance(); return new Token(TOKEN_LBRACKET, '{'); }
                if (this.currentChar === '}') { this.advance(); return new Token(TOKEN_RBRACKET, '}'); }
                if (this.currentChar === ',') { this.advance(); continue; }
                if (this.currentChar === ';') { this.advance(); continue; }
                if (this.currentChar === ':') { this.advance(); continue; }
                if (this.currentChar === ' ') { this.advance(); continue; }
                if (this.currentChar === '|') { this.advance(); return new Token(TOKEN_PIPE, '|'); }

                let cmd = '';
                while (this.currentChar !== null && /[a-zA-Z]/.test(this.currentChar)) {
                    cmd += this.currentChar;
                    this.advance();
                }

                if (cmd === 'frac') return new Token(TOKEN_LATEX_FRAC, cmd);
                if (cmd === 'sqrt') return new Token(TOKEN_LATEX_SQRT, cmd);
                if (cmd === 'sum') return new Token(TOKEN_LATEX_SUM, cmd);
                if (cmd === 'prod') return new Token(TOKEN_LATEX_PROD, cmd);
                if (cmd === 'int') return new Token(TOKEN_LATEX_INT, cmd);
                if (cmd === 'lim') return new Token(TOKEN_LATEX_LIMIT, cmd);

                if (cmd === 'left' || cmd === 'right') {
                    this.skipWhitespace();
                    if (this.currentChar === '.') {
                        this.advance();
                    }
                    continue;
                }

                if (cmd === 'cdot') return new Token(TOKEN_STAR, '*');
                if (cmd === 'times') return new Token(TOKEN_STAR, '*');
                if (cmd === 'div') return new Token(TOKEN_SLASH, '/');
                if (cmd === 'le') return new Token(TOKEN_LE, '<=');
                if (cmd === 'ge') return new Token(TOKEN_GE, '>=');
                if (cmd === 'ne' || cmd === 'neq') return new Token(TOKEN_NEQ, '!=');
                if (cmd === 'to' || cmd === 'rightarrow') return new Token(TOKEN_IDENTIFIER, 'to');

                if (cmd === 'pi') return new Token(TOKEN_IDENTIFIER, 'pi');
                if (cmd === 'infty') return new Token(TOKEN_IDENTIFIER, 'infinity');

                return new Token(TOKEN_IDENTIFIER, cmd);
            }

            if (/\d/.test(this.currentChar)) {
                return new Token(TOKEN_NUMBER, this.number());
            }
            if (this.currentChar === '.') {
                if (this.peek() === '.') {
                    this.advance(); this.advance();
                    return new Token(TOKEN_RANGE, '..');
                }
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

                // Control Flow
                if (lowerId === 'if') return new Token(TOKEN_IF, id);
                if (lowerId === 'else') return new Token(TOKEN_ELSE, id);
                if (lowerId === 'while') return new Token(TOKEN_WHILE, id);
                if (lowerId === 'for') return new Token(TOKEN_FOR, id);
                if (lowerId === 'return') return new Token(TOKEN_RETURN, id);
                if (lowerId === 'break') return new Token(TOKEN_BREAK, id);
                if (lowerId === 'continue') return new Token(TOKEN_CONTINUE, id);

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
                this.advance();
                return new Token('BANG', '!');
            }
            if (this.currentChar === '<') {
                if (this.peek() === '=') {
                    this.advance();
                    this.advance();
                    if (this.currentChar === '>') { // <=>
                        this.advance();
                        return new Token(TOKEN_IFF, '<=>');
                    }
                    return new Token(TOKEN_LE, '<=');
                }
                if (this.peek() === '-') {
                    const next = this.text[this.pos + 2];
                    if (next === '>') { // <->
                        this.advance();
                        this.advance();
                        this.advance();
                        return new Token(TOKEN_IFF, '<->');
                    }
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
                    this.advance();
                    this.advance();
                    return new Token(TOKEN_BOOL_EQ, '==');
                }
                if (this.peek() === '>') {
                    this.advance();
                    this.advance();
                    return new Token(TOKEN_IMPLIES, '=>');
                }
                this.advance();
                return new Token(TOKEN_EQ, '=');
            }
            if (this.currentChar === '+') { this.advance(); return new Token(TOKEN_PLUS, '+'); }
            if (this.currentChar === '-') {
                if (this.peek() === '>') {
                    this.advance();
                    this.advance();
                    return new Token(TOKEN_IMPLIES, '->');
                }
                this.advance();
                return new Token(TOKEN_MINUS, '-');
            }
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
            if (this.currentChar === '_') { this.advance(); return new Token(TOKEN_UNDERSCORE, '_'); }
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
                this.advance();
                return new Token(TOKEN_PIPE, '|');
            }
            if (this.currentChar === '(') { this.advance(); return new Token(TOKEN_LPAREN, '('); }
            if (this.currentChar === ')') { this.advance(); return new Token(TOKEN_RPAREN, ')'); }
            if (this.currentChar === '[') { this.advance(); return new Token(TOKEN_LBRACKET, '['); }
            if (this.currentChar === ']') { this.advance(); return new Token(TOKEN_RBRACKET, ']'); }
            if (this.currentChar === '{') { this.advance(); return new Token(TOKEN_LBRACE, '{'); }
            if (this.currentChar === '}') { this.advance(); return new Token(TOKEN_RBRACE, '}'); }
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
        this.absDepth = 0;
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
        while (this.currentToken.type === TOKEN_LBRACKET) {
             this.eat(TOKEN_LBRACKET);
             const index = this.statement();
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
            let name = token.value;
            this.eat(TOKEN_IDENTIFIER);

            // Handle identifier subscript (name merging)
            while (this.currentToken.type === TOKEN_UNDERSCORE) {
                this.eat(TOKEN_UNDERSCORE);
                let suffix = '';
                if (this.currentToken.type === TOKEN_NUMBER) {
                    suffix = this.currentToken.value.toString();
                    this.eat(TOKEN_NUMBER);
                } else if (this.currentToken.type === TOKEN_IDENTIFIER) {
                    suffix = this.currentToken.value;
                    this.eat(TOKEN_IDENTIFIER);
                } else if (this.currentToken.type === TOKEN_LBRACE) {
                    this.eat(TOKEN_LBRACE);
                    // Consume everything inside braces as part of the name
                    // This is a simplification to allow x_{foo} -> x_foo
                    while (this.currentToken.type !== TOKEN_RBRACE && this.currentToken.type !== TOKEN_EOF) {
                        if (this.currentToken.type === TOKEN_IDENTIFIER ||
                            this.currentToken.type === TOKEN_NUMBER) {
                            suffix += this.currentToken.value.toString();
                        } else {
                             // Skip other chars? or error?
                        }
                        this.eat(this.currentToken.type);
                    }
                    this.eat(TOKEN_RBRACE);
                }
                name += '_' + suffix;
            }

            // After merging subscripts, check for function calls
            // Handle sin^2(x)
            if (this.currentToken.type === TOKEN_CARET) {
                const trigFunctions = ['sin', 'cos', 'tan', 'sec', 'csc', 'cot', 'sinh', 'cosh', 'tanh'];
                if (trigFunctions.includes(name)) {
                    this.eat(TOKEN_CARET);
                    const exponent = this.unary();

                    if (this.currentToken.type === TOKEN_LPAREN) {
                        this.eat(TOKEN_LPAREN);
                        const arg = this.statement();
                        this.eat(TOKEN_RPAREN);
                        return new Pow(new Call(name, [arg]), exponent);
                    } else if (this.isImplicitMulStart(this.currentToken)) {
                         const arg = this.term();
                         return new Pow(new Call(name, [arg]), exponent);
                    }
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
                const knownFunctions = [
                    'sin', 'cos', 'tan', 'asin', 'acos', 'atan',
                    'sinh', 'cosh', 'tanh', 'asinh', 'acosh', 'atanh',
                    'sec', 'csc', 'cot', 'asec', 'acsc', 'acot',
                    'sqrt', 'log', 'ln', 'exp', 'abs', 'sign',
                    'floor', 'ceil', 'round', 'fact', 'factorial', 'gamma'
                ];
                if (knownFunctions.includes(name)) {
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
        } else if (token.type === TOKEN_LBRACE) {
            this.eat(TOKEN_LBRACE);
            if (this.currentToken.type === TOKEN_RBRACE) {
                // Empty block?
                this.eat(TOKEN_RBRACE);
                return new Block([]);
            }
            // Block parsing handled by statement list if block is not expression?
            // But { a; b; c } is a block.
            // Our parse() loop handles the top level.
            // Here we want to parse statements inside braces.
            const stmts = [];
            while (this.currentToken.type !== TOKEN_RBRACE && this.currentToken.type !== TOKEN_EOF) {
                stmts.push(this.statement());
                while (this.currentToken.type === TOKEN_SEMI) {
                    this.eat(TOKEN_SEMI);
                }
            }
            this.eat(TOKEN_RBRACE);
            return new Block(stmts);
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
        } else if (token.type === TOKEN_PIPE) {
            this.eat(TOKEN_PIPE);
            this.absDepth++;
            const node = this.statement();
            this.absDepth--;
            this.eat(TOKEN_PIPE);
            return new Call('abs', [node]);
        } else if (token.type === TOKEN_LATEX_FRAC) {
            this.eat(TOKEN_LATEX_FRAC);
            const num = this.atom();
            const den = this.atom();
            return new Div(num, den);
        } else if (token.type === TOKEN_LATEX_SQRT) {
            this.eat(TOKEN_LATEX_SQRT);
            if (this.currentToken.type === TOKEN_LBRACKET) {
                this.eat(TOKEN_LBRACKET);
                const index = this.statement();
                this.eat(TOKEN_RBRACKET);
                const arg = this.atom();
                return new Pow(arg, new Div(new Num(1), index));
            } else {
                const arg = this.atom();
                return new Call('sqrt', [arg]);
            }
        } else if (token.type === TOKEN_LATEX_SUM || token.type === TOKEN_LATEX_PROD) {
            const isSum = token.type === TOKEN_LATEX_SUM;
            const funcName = isSum ? 'sum' : 'product';
            this.eat(token.type);

            let variable = new Sym('i');
            let start = new Num(0);
            let end = new Sym('n');
            let hasLimits = false;

            if (this.currentToken.type === TOKEN_UNDERSCORE) {
                this.eat(TOKEN_UNDERSCORE);
                hasLimits = true;
                if (this.currentToken.type === TOKEN_LBRACE) {
                    this.eat(TOKEN_LBRACE);
                    const parsed = this.statement();
                    this.eat(TOKEN_RBRACE);

                    if (parsed.constructor.name === 'Eq' || parsed.constructor.name === 'Assignment') {
                        variable = parsed.left;
                        start = parsed.right;
                    }
                } else {
                     variable = this.atom();
                }
            }

            if (this.currentToken.type === TOKEN_CARET) {
                this.eat(TOKEN_CARET);
                if (this.currentToken.type === TOKEN_LBRACE) {
                    this.eat(TOKEN_LBRACE);
                    end = this.statement();
                    this.eat(TOKEN_RBRACE);
                } else {
                    end = this.atom();
                }
            }

            const term = this.term();
            if (hasLimits) {
                return new Call(funcName, [term, variable, start, end]);
            } else {
                return new Call(funcName, [term]);
            }

        } else if (token.type === TOKEN_LATEX_INT) {
            this.eat(TOKEN_LATEX_INT);
            let lower = null;
            let upper = null;

            if (this.currentToken.type === TOKEN_UNDERSCORE) {
                this.eat(TOKEN_UNDERSCORE);
                if (this.currentToken.type === TOKEN_LBRACE) {
                     this.eat(TOKEN_LBRACE);
                     lower = this.statement();
                     this.eat(TOKEN_RBRACE);
                } else {
                     lower = this.atom();
                }
            }

            if (this.currentToken.type === TOKEN_CARET) {
                this.eat(TOKEN_CARET);
                if (this.currentToken.type === TOKEN_LBRACE) {
                     this.eat(TOKEN_LBRACE);
                     upper = this.statement();
                     this.eat(TOKEN_RBRACE);
                } else {
                     upper = this.atom();
                }
            }

            let term = this.term();

            // Heuristic to strip 'dx'
            let variable = new Sym('x');
            let integrand = term;

            // Check if term is Mul(..., d<var>)
            if (term.constructor.name === 'Mul') {
                const operands = [];
                function collect(node) {
                    if (node.constructor.name === 'Mul') {
                        collect(node.left);
                        collect(node.right);
                    } else {
                        operands.push(node);
                    }
                }
                collect(term);

                if (operands.length > 0) {
                    const last = operands[operands.length - 1];
                    if (last.constructor.name === 'Sym' && last.name.startsWith('d') && last.name.length > 1) {
                        const varName = last.name.substring(1);
                        variable = new Sym(varName);
                        // Rebuild term without last operand
                        operands.pop();
                        // Reduce operands back to Mul
                        if (operands.length === 0) {
                            integrand = new Num(1);
                        } else {
                            integrand = operands[0];
                            for (let i = 1; i < operands.length; i++) {
                                integrand = new Mul(integrand, operands[i]);
                            }
                        }
                    }
                }
            } else if (term.constructor.name === 'Sym' && term.name.startsWith('d') && term.name.length > 1) {
                 // \int dx case -> 1 dx
                 const varName = term.name.substring(1);
                 variable = new Sym(varName);
                 integrand = new Num(1);
            }

            if (lower && upper) {
                return new Call('integrate', [integrand, variable, lower, upper]);
            }
            return new Call('integrate', [integrand, variable]);

        } else if (token.type === TOKEN_LATEX_LIMIT) {
             this.eat(TOKEN_LATEX_LIMIT);
             let variable = new Sym('x');
             let target = new Num(0);

             if (this.currentToken.type === TOKEN_UNDERSCORE) {
                 this.eat(TOKEN_UNDERSCORE);
                 if (this.currentToken.type === TOKEN_LBRACE) {
                     this.eat(TOKEN_LBRACE);

                     if (this.currentToken.type === TOKEN_IDENTIFIER) {
                         variable = new Sym(this.currentToken.value);
                         this.eat(TOKEN_IDENTIFIER);
                     }

                     // Consume 'to'
                     if (this.currentToken.type === TOKEN_IDENTIFIER && this.currentToken.value === 'to') {
                         this.eat(TOKEN_IDENTIFIER);
                     } else if (this.currentToken.type === TOKEN_MINUS) { // ->
                         this.eat(TOKEN_MINUS);
                         if (this.currentToken.type === TOKEN_GT) {
                             this.eat(TOKEN_GT);
                         }
                     }

                     target = this.statement();
                     this.eat(TOKEN_RBRACE);
                 }
             }
             const term = this.atom();
             return new Call('limit', [term, variable, target]);
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
        while (this.isImplicitMulStart(this.currentToken)) {
             if (this.currentToken.type === TOKEN_PIPE && this.absDepth > 0) {
                 break;
             }
             const right = this.unary();
             node = new Mul(node, right);
        }
        return node;
    }

    isImplicitMulStart(token) {
        return token.type === TOKEN_NUMBER ||
               token.type === TOKEN_IDENTIFIER ||
               token.type === TOKEN_LPAREN ||
               token.type === TOKEN_LBRACKET ||
               token.type === TOKEN_LBRACE ||
               token.type === TOKEN_PIPE ||
               token.type === TOKEN_LATEX_FRAC ||
               token.type === TOKEN_LATEX_SQRT ||
               token.type === TOKEN_LATEX_SUM ||
               token.type === TOKEN_LATEX_PROD ||
               token.type === TOKEN_LATEX_INT ||
               token.type === TOKEN_LATEX_LIMIT ||
               // New control flow keywords should NOT trigger implicit mul (e.g. "x if" -> Mul(x, if) is wrong)
               // But they are identifiers so Lexer might return them as keywords.
               // We need to ensure IF, ELSE etc are not treated as implicit mul.
               // Token types are specific.
               token.type === TOKEN_IF || token.type === TOKEN_WHILE || token.type === TOKEN_FOR || token.type === TOKEN_RETURN;
               // Wait, if I type "x if", implicitMul calls unary() which calls factor() which calls atom().
               // atom() expects NUMBER, IDENTIFIER, etc.
               // IF is a token type TOKEN_IF. atom() does not handle TOKEN_IF.
               // So unary() will fail.
               // This is correct behavior. "x if" is syntax error.
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

    rangeExpr() {
        let node = this.arithExpr();
        if (this.currentToken.type === TOKEN_RANGE) {
            this.eat(TOKEN_RANGE);
            const right = this.arithExpr();
            return new Call('range', [node, new Add(right, new Num(1)), new Num(1)]);
        }
        return node;
    }

    compExpr() {
        let node = this.rangeExpr();
        if (this.currentToken.type === TOKEN_EQ) {
            this.eat(TOKEN_EQ);
            node = new Eq(node, this.rangeExpr());
        } else if (this.currentToken.type === TOKEN_BOOL_EQ) {
            this.eat(TOKEN_BOOL_EQ);
            node = new BooleanEq(node, this.rangeExpr());
        } else if (this.currentToken.type === TOKEN_NEQ) {
            this.eat(TOKEN_NEQ);
            node = new Neq(node, this.rangeExpr());
        } else if (this.currentToken.type === TOKEN_LT) {
            this.eat(TOKEN_LT);
            node = new Lt(node, this.rangeExpr());
        } else if (this.currentToken.type === TOKEN_GT) {
            this.eat(TOKEN_GT);
            node = new Gt(node, this.rangeExpr());
        } else if (this.currentToken.type === TOKEN_LE) {
            this.eat(TOKEN_LE);
            node = new Le(node, this.rangeExpr());
        } else if (this.currentToken.type === TOKEN_GE) {
            this.eat(TOKEN_GE);
            node = new Ge(node, this.rangeExpr());
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

    impliesExpr() {
        let node = this.logicExpr();
        if (this.currentToken.type === TOKEN_IMPLIES) {
            this.eat(TOKEN_IMPLIES);
            node = new Implies(node, this.impliesExpr());
        }
        return node;
    }

    iffExpr() {
        let node = this.impliesExpr();
        while (this.currentToken.type === TOKEN_IFF) {
            this.eat(TOKEN_IFF);
            node = new Iff(node, this.impliesExpr());
        }
        return node;
    }

    equation() {
        return this.iffExpr();
    }

    statement() {
        if (this.currentToken.type === TOKEN_IF) {
            this.eat(TOKEN_IF);
            this.eat(TOKEN_LPAREN);
            const condition = this.equation();
            this.eat(TOKEN_RPAREN);

            // Allow single statement or block
            const trueBlock = this.statement();

            let falseBlock = null;
            if (this.currentToken.type === TOKEN_ELSE) {
                this.eat(TOKEN_ELSE);
                falseBlock = this.statement();
            }
            return new If(condition, trueBlock, falseBlock);

        } else if (this.currentToken.type === TOKEN_WHILE) {
            this.eat(TOKEN_WHILE);
            this.eat(TOKEN_LPAREN);
            const condition = this.equation();
            this.eat(TOKEN_RPAREN);
            const body = this.statement();
            return new While(condition, body);

        } else if (this.currentToken.type === TOKEN_FOR) {
            this.eat(TOKEN_FOR);
            this.eat(TOKEN_LPAREN);
            // for (init; cond; step)

            let init = null;
            if (this.currentToken.type !== TOKEN_SEMI) {
                init = this.statement(); // This consumes assignment
            }
            this.eat(TOKEN_SEMI);

            let cond = null;
            if (this.currentToken.type !== TOKEN_SEMI) {
                cond = this.equation();
            }
            this.eat(TOKEN_SEMI);

            let step = null;
            if (this.currentToken.type !== TOKEN_RPAREN) {
                step = this.statement(); // Assignment usually
            }
            this.eat(TOKEN_RPAREN);

            const body = this.statement();
            return new For(init, cond, step, body);

        } else if (this.currentToken.type === TOKEN_RETURN) {
            this.eat(TOKEN_RETURN);
            let value = null;
            if (this.currentToken.type !== TOKEN_SEMI && this.currentToken.type !== TOKEN_EOF && this.currentToken.type !== TOKEN_RBRACE) {
                value = this.equation();
            }
            // Semi is consumed by block parser or here?
            // Usually statements don't consume their own semi if inside a block loop, but here it's cleaner to leave it?
            // Our block parser `parse()` and `LBRACE` handler consume SEMI.
            return new Return(value);

        } else if (this.currentToken.type === TOKEN_BREAK) {
            this.eat(TOKEN_BREAK);
            return new Break();

        } else if (this.currentToken.type === TOKEN_CONTINUE) {
            this.eat(TOKEN_CONTINUE);
            return new Continue();

        } else if (this.currentToken.type === TOKEN_LBRACE) {
            // Forward to atom handling which calls block?
            // atom() handles LBRACE for simple expressions or blocks.
            // But if we are in statement, and see LBRACE, it's a block statement.
            return this.atom();
        }

        if (this.currentToken.type === TOKEN_IDENTIFIER) {
            const savedPos = this.lexer.pos;
            const savedChar = this.lexer.currentChar;
            const savedToken = this.currentToken;

            const name = this.currentToken.value;
            this.eat(TOKEN_IDENTIFIER);

            if (this.currentToken.type === TOKEN_ASSIGN) {
                this.eat(TOKEN_ASSIGN);
                const value = this.statement();
                return new Assignment(new Sym(name), value);
            }

            if (this.currentToken.type === TOKEN_LPAREN) {
                this.eat(TOKEN_LPAREN);
                const args = [];
                let isFuncDef = true;

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

(function() {
    const exports = {
        Token, Lexer, Parser
    };
    if (typeof globalThis !== 'undefined') {
        Object.assign(globalThis, exports);
    }
    if (typeof module !== 'undefined' && module.exports) {
        if (typeof globalThis.Call === 'undefined') {
             try {
                 const expr = require('./expression.js');
                 Object.assign(globalThis, expr);
             } catch (e) {}
        }
        module.exports = exports;
    }
})();
