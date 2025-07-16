"""
    get_qsymbols(expr)

Extract all quantum operators from an expression.

Recursively traverses a symbolic expression and returns all quantum operators
(objects that are neither numbers nor symbolic parameters).

# Arguments
- `expr`: A symbolic expression containing quantum operators and/or parameters

# Returns
- `Set`: Set of quantum operators found in the expression

"""
function get_qsymbols(expr)
    if Symbolics.iscall(expr)
        return union(get_qsymbols.(Symbolics.arguments(expr))...)
    else
        return (expr isa Number) || (expr isa SymbolicUtils.Symbolic) ? Set() : Set([expr])
    end
end

"""
    get_numsymbols(expr)

Extract all symbolic parameters from an expression.

Recursively traverses a symbolic expression and returns all symbolic parameters
(SymbolicUtils.Symbolic objects like those created with @cnumbers).

# Arguments
- `expr`: A symbolic expression containing quantum operators and/or parameters

# Returns
- `Set`: Set of symbolic parameters found in the expression

"""
function get_numsymbols(expr)
    if Symbolics.iscall(expr)
        return union(get_numsymbols.(Symbolics.arguments(expr))...)
    else
        return (expr isa SymbolicUtils.Symbolic) ? Set([expr]) : Set() 
    end
end

"""
    get_additive_terms(expr)

Extract additive terms from a quantum operator expression.

Takes an expression containing quantum operators and returns a list of terms
that contain no addition, only multiplication. Summing all returned terms
results in the original expression.

# Arguments
- `expr`: A symbolic expression containing quantum operators

# Returns
- `Vector`: List of terms without addition operators

"""
function get_additive_terms(expr)
    if Symbolics.iscall(expr)
        op = Symbolics.operation(expr)
        args = Symbolics.arguments(expr)
        
        if op === (+)
            # If this is an addition, recursively get terms from each argument
            terms = []
            for arg in args
                append!(terms, get_additive_terms(arg))
            end
            return terms
        else
            # If this is not an addition (multiplication, function call, etc.), 
            # return the whole expression as a single term
            return [expr]
        end
    else
        # If not a tree (atom), return as single term
        return [expr]
    end
end
