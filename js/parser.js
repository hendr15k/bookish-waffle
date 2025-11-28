
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
            if (/\d/.test(this.currentChar)) {
                return new Token(TOKEN_NUMBER, this.number());
            }
            if (/[a-zA-Z]/.test(this.currentChar)) {
                return new Token(TOKEN_IDENTIFIER, this.identifier());
            }
            if (this.currentChar === ':') {
                if (this.peek() === '=') {
                    this.advance();
                    this.advance();
                    return new Token(TOKEN_ASSIGN, ':=');
                }
                throw new Error("Invalid character: :");
            }
            if (this.currentChar === '=') {
                this.advance();
                return new Token(TOKEN_EQ, '=');
            }
            if (this.currentChar === '+') { this.advance(); return new Token(TOKEN_PLUS, '+'); }
            if (this.currentChar === '-') { this.advance(); return new Token(TOKEN_MINUS, '-'); }
            if (this.currentChar === '*') { this.advance(); return new Token(TOKEN_STAR, '*'); }
            if (this.currentChar === '/') { this.advance(); return new Token(TOKEN_SLASH, '/'); }
            if (this.currentChar === '^') { this.advance(); return new Token(TOKEN_CARET, '^'); }
            if (this.currentChar === '(') { this.advance(); return new Token(TOKEN_LPAREN, '('); }
            if (this.currentChar === ')') { this.advance(); return new Token(TOKEN_RPAREN, ')'); }
            if (this.currentChar === '[') { this.advance(); return new Token(TOKEN_LBRACKET, '['); }
            if (this.currentChar === ']') { this.advance(); return new Token(TOKEN_RBRACKET, ']'); }
            if (this.currentChar === ',') { this.advance(); return new Token(TOKEN_COMMA, ','); }

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
        const token = this.currentToken;
        if (token.type === TOKEN_NUMBER) {
            this.eat(TOKEN_NUMBER);
            return new Num(token.value);
        } else if (token.type === TOKEN_IDENTIFIER) {
            const name = token.value;
            this.eat(TOKEN_IDENTIFIER);
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
                return new Symbol(name);
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
        } else if (token.type === TOKEN_PLUS) {
            this.eat(TOKEN_PLUS);
            return this.factor();
        } else if (token.type === TOKEN_MINUS) {
            this.eat(TOKEN_MINUS);
            return new Mul(new Num(-1), this.factor());
        }
        this.error();
    }

    power() {
        let node = this.factor();
        if (this.currentToken.type === TOKEN_CARET) {
            this.eat(TOKEN_CARET);
            node = new Pow(node, this.power());
        }
        return node;
    }

    term() {
        let node = this.power();
        while (this.currentToken.type === TOKEN_STAR ||
               this.currentToken.type === TOKEN_SLASH ||
               this.isImplicitMulStart(this.currentToken)) {

            if (this.currentToken.type === TOKEN_STAR) {
                this.eat(TOKEN_STAR);
                node = new Mul(node, this.power());
            } else if (this.currentToken.type === TOKEN_SLASH) {
                this.eat(TOKEN_SLASH);
                node = new Div(node, this.power());
            } else {
                // Implicit multiplication: no operator consumed
                node = new Mul(node, this.power());
            }
        }
        return node;
    }

    isImplicitMulStart(token) {
        // Tokens that can start a power/factor:
        // NUMBER, IDENTIFIER, LPAREN, LBRACKET
        // (MINUS/PLUS are handled in expr, so implicit mul like `x -y` isn't standard, usually `x * -y` needs parens or explicit `*`)
        // However, `2x` works. `x y` works. `x(y)` is parsed in factor, so if we are here, we saw a factor and next is something else.
        // e.g. `x` (Symbol) then `y` (Symbol).
        return token.type === TOKEN_NUMBER ||
               token.type === TOKEN_IDENTIFIER ||
               token.type === TOKEN_LPAREN ||
               token.type === TOKEN_LBRACKET;
    }

    expr() {
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

    equation() {
        let node = this.expr();
        if (this.currentToken.type === TOKEN_EQ) {
            this.eat(TOKEN_EQ);
            node = new Eq(node, this.expr());
        }
        return node;
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
                return new Assignment(new Symbol(name), value);
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
        return this.statement();
    }
}
