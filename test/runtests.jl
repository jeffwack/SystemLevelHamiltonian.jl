using SystemLevelHamiltonian
using Test
using QuantumCumulants
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
    
end
