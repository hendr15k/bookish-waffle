
class GeminiAdapter {
    constructor(casInstance) {
        this.cas = casInstance;
        // If casInstance is not provided, try to find global one
        if (!this.cas && typeof cas !== 'undefined') {
            this.cas = cas;
        }
    }

    /**
     * Returns the tool schema compatible with Gemini/OpenAI function calling.
     */
    getTools() {
        return [
            {
                name: "calculate",
                description: "Evaluates a mathematical expression or performs a calculation.",
                parameters: {
                    type: "object",
                    properties: {
                        expression: {
                            type: "string",
                            description: "The mathematical expression to evaluate (e.g., '1+1', 'sin(pi/4)', 'diff(x^2, x)')"
                        }
                    },
                    required: ["expression"]
                }
            },
            {
                name: "solve_equation",
                description: "Solves an equation for a specific variable.",
                parameters: {
                    type: "object",
                    properties: {
                        equation: {
                            type: "string",
                            description: "The equation to solve (e.g., 'x^2 - 4 = 0')"
                        },
                        variable: {
                            type: "string",
                            description: "The variable to solve for (e.g., 'x')"
                        }
                    },
                    required: ["equation", "variable"]
                }
            },
            {
                name: "plot_function",
                description: "Generates a 2D plot of a function.",
                parameters: {
                    type: "object",
                    properties: {
                        expression: {
                            type: "string",
                            description: "The function expression to plot (e.g., 'sin(x)')"
                        },
                        variable: {
                            type: "string",
                            description: "The independent variable (e.g., 'x')"
                        },
                        min: {
                            type: "number",
                            description: "Minimum x value"
                        },
                        max: {
                            type: "number",
                            description: "Maximum x value"
                        }
                    },
                    required: ["expression", "variable"]
                }
            },
            {
                name: "run_script",
                description: "Runs a multi-line CAS script.",
                parameters: {
                    type: "object",
                    properties: {
                        script: {
                            type: "string",
                            description: "The script content with statements separated by semicolons."
                        }
                    },
                    required: ["script"]
                }
            }
        ];
    }

    /**
     * Executes a tool by name with the given parameters.
     * @param {string} name - The name of the tool.
     * @param {object} params - The parameters for the tool.
     * @returns {object} - The result of the execution.
     */
    executeTool(name, params) {
        if (!this.cas) throw new Error("CAS instance not initialized in GeminiAdapter");

        // Helper to evaluate string command
        const evalCommand = (cmd) => {
            try {
                // We need Lexer and Parser. In browser they are global.
                // In Node they might need to be imported, but we assume global for now based on project structure.
                let LexerClass, ParserClass;

                if (typeof Lexer !== 'undefined') LexerClass = Lexer;
                else if (typeof globalThis.Lexer !== 'undefined') LexerClass = globalThis.Lexer;

                if (typeof Parser !== 'undefined') ParserClass = Parser;
                else if (typeof globalThis.Parser !== 'undefined') ParserClass = globalThis.Parser;

                if (!LexerClass || !ParserClass) {
                    throw new Error("Lexer/Parser not found");
                }

                const lexer = new LexerClass(cmd);
                const parser = new ParserClass(lexer);
                const tree = parser.parse();
                const result = this.cas.evaluate(tree);

                // Format result
                return {
                    status: "success",
                    result: result.toString(),
                    latex: (result.toLatex) ? result.toLatex() : null,
                    type: result.type || 'expression',
                    raw: (result.type === 'plot') ? result : undefined // Include raw object for plots
                };
            } catch (e) {
                return {
                    status: "error",
                    error: e.message
                };
            }
        };

        if (name === "calculate") {
            return evalCommand(params.expression);
        }

        if (name === "solve_equation") {
            return evalCommand(`solve(${params.equation}, ${params.variable})`);
        }

        if (name === "plot_function") {
            const min = params.min !== undefined ? params.min : -10;
            const max = params.max !== undefined ? params.max : 10;
            return evalCommand(`plot(${params.expression}, ${params.variable}, ${min}, ${max})`);
        }

        if (name === "run_script") {
            // Normalize script to avoid double semicolons or empty statements which parser dislikes
            // 1. Remove comments if needed (Lexer handles them, but simple cleanup helps)
            // 2. Normalize newlines to semicolons, but avoid ;;
            let script = params.script;
            // Replace newlines with semicolons if not preceded by semicolon
            script = script.replace(/([^;])\n/g, '$1;');
            // Remove newlines
            script = script.replace(/\n/g, '');
            // Remove double semicolons
            script = script.replace(/;+/g, ';');
            // Remove leading/trailing semicolons
            script = script.replace(/^;+/, '').replace(/;+$/, '');

            return evalCommand(script);
        }

        return { status: "error", error: `Unknown tool: ${name}` };
    }
}

// Export
(function() {
    const exports = {
        GeminiAdapter
    };
    if (typeof globalThis !== 'undefined') {
        Object.assign(globalThis, exports);
    }
    if (typeof module !== 'undefined' && module.exports) {
        module.exports = exports;
    }
})();
