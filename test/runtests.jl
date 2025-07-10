using SystemLevelHamiltonian
using Test
using SecondQuantizedAlgebra
using Symbolics

@testset "SystemLevelHamiltonian.jl" begin
    
    @testset "get_additive_terms" begin
        # Create test Hilbert space and operators
        hf = FockSpace(:cavity)
        @qnumbers a::Destroy(hf)
        
        # Create symbolic parameters
        @cnumbers g ω κ
        
        # Test 1: Single term (no addition)
        single_term = g * a
        terms = get_additive_terms(single_term)
        @test length(terms) == 1
        @test isequal(terms[1], g * a)
        @test isequal(sum(terms), single_term)
        
        # Test 2: Simple addition of two terms
        two_terms = g * a + ω * a'
        terms = get_additive_terms(two_terms)
        @test length(terms) == 2
        @test isequal(sum(terms), two_terms)
        
        # Test 3: Multiple terms with different operators
        multi_terms = g * a + ω * a' + κ * a' * a
        terms = get_additive_terms(multi_terms)
        @test length(terms) == 3
        @test isequal(sum(terms), multi_terms)
        
        # Test 4: Nested additions
        nested = (g * a + ω * a') + κ * a' * a
        terms = get_additive_terms(nested)
        @test length(terms) == 3
        @test isequal(sum(terms), nested)
        
        # Test 5: Single number
        single_num = 5
        terms = get_additive_terms(single_num)
        @test length(terms) == 1
        @test terms[1] == 5
        @test isequal(sum(terms), single_num)
        
        # Test 6: Addition with numbers (numeric terms grouped together)
        num_terms = 3 + 2*g + ω * a
        terms = get_additive_terms(num_terms)
        @test length(terms) == 3
        @test iszero(simplify(sum(terms) - num_terms))
        
        # Test 7: Complex multiplication term (commutes to two terms due to [a,a']=1)
        complex_mult = g * ω * a * a'
        terms = get_additive_terms(complex_mult)
        @test length(terms) == 2  # g*ω*(a′*a) + g*ω due to commutation
        @test isequal(sum(terms), complex_mult)
        
        # Test 8: Zero term
        zero_term = 0 * a
        terms = get_additive_terms(zero_term)
        @test length(terms) == 1
        @test isequal(sum(terms), zero_term)
    end
    
    @testset "get_qsymbols" begin
        # Create test Hilbert space and operators
        hf = FockSpace(:cavity)
        @qnumbers a::Destroy(hf)
        
        # Create symbolic parameters
        @cnumbers g ω κ
        
        # Test 1: Single operator
        qsyms = get_qsymbols(a)
        @test length(qsyms) == 1
        @test a in qsyms
        
        # Test 2: Multiple operators
        expr = g * a + ω * a'
        qsyms = get_qsymbols(expr)
        @test length(qsyms) == 2
        @test a in qsyms
        @test a' in qsyms
        
        # Test 3: Just parameters (no quantum operators)
        param_expr = g + ω * κ
        qsyms = get_qsymbols(param_expr)
        @test length(qsyms) == 0
        
        # Test 4: Complex expression
        complex_expr = g * a' * a + κ * a
        qsyms = get_qsymbols(complex_expr)
        @test length(qsyms) == 2
        @test a in qsyms
        @test a' in qsyms
    end
    
    @testset "get_numsymbols" begin
        # Create test Hilbert space and operators
        hf = FockSpace(:cavity)
        @qnumbers a::Destroy(hf)
        
        # Create symbolic parameters
        @cnumbers g ω κ
        
        # Test 1: Single parameter
        numsyms = get_numsymbols(g)
        @test length(numsyms) == 1
        @test g in numsyms
        
        # Test 2: Multiple parameters
        expr = g * ω + κ
        numsyms = get_numsymbols(expr)
        @test length(numsyms) == 3
        @test g in numsyms
        @test ω in numsyms
        @test κ in numsyms
        
        # Test 3: Parameters with operators
        mixed_expr = g * a + ω * a'
        numsyms = get_numsymbols(mixed_expr)
        @test length(numsyms) == 2
        @test g in numsyms
        @test ω in numsyms
        
        # Test 4: Just operators (no parameters)
        op_expr = a + a'
        numsyms = get_numsymbols(op_expr)
        @test length(numsyms) == 0
    end
    
    
    @testset "SLH struct and basic operations" begin
        # Create test system
        hf = FockSpace(:cavity)
        @qnumbers a::Destroy(hf)
        @cnumbers ω κ g
        
        H = ω * a' * a
        L = [√κ * a]
        S = [1]
        sys = SLH(:test, [:in], [:out], S, L, H)
        
        # Test 1: SLH construction
        @test sys.name == :test
        @test sys.inputs == [:in]
        @test sys.outputs == [:out]
        @test sys.S == [1]
        @test length(sys.L) == 1
        
        
        # Test 3: operators function
        ops = operators(sys)
        @test length(ops) == 2
        @test a in ops
        @test a' in ops
        
        # Test 4: parameters function
        params = parameters(sys)
        @test ω in params
        @test κ in params
    end
    
    @testset "cascade" begin
        # Create two test systems
        hf1 = FockSpace(:cavity1)
        hf2 = FockSpace(:cavity2)
        @qnumbers a::Destroy(hf1) b::Destroy(hf2)
        @cnumbers ω1 ω2 κ1 κ2
        
        H1 = ω1 * a' * a
        L1 = [√κ1 * a]
        S1 = [1]
        sys1 = SLH(:sys1, [:in1], [:out1], S1, L1, H1)
        
        H2 = ω2 * b' * b
        L2 = [√κ2 * b]
        S2 = [1]
        sys2 = SLH(:sys2, [:in2], [:out2], S2, L2, H2)
        
        # Test concatenation
        combined = concatenate(:combined, [sys1, sys2])
        
        @test combined.name == :combined
        @test length(combined.inputs) == 2
        @test length(combined.outputs) == 2
        @test size(combined.S) == (2, 2)
        @test length(combined.L) == 2
        
        # Check that inputs and outputs are properly renamed
        @test :in1_sys1 in combined.inputs
        @test :in2_sys2 in combined.inputs
        @test :out1_sys1 in combined.outputs
        @test :out2_sys2 in combined.outputs
    
        cascaded = feedbackreduce(sys,:out1,:in2)
    end

    
end
