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
    
    @testset "check_linear" begin
        # Create test Hilbert space and operators
        hf = FockSpace(:cavity)
        @qnumbers a::Destroy(hf)
        
        # Create symbolic parameters
        @cnumbers g ω κ
        
        # Test 1: Single operator (linear)
        single_op = a
        @test check_linear(single_op) == true
        
        # Test 2: Single multiplication term with two operators (linear)
        # Note: We need to test the actual multiplicative terms, not the expanded expressions
        two_ops_term = get_additive_terms(a * a')[1]  # Get the (a′*a) term
        @test check_linear(two_ops_term) == true
        
        # Test 3: Single multiplication term with three operators (nonlinear)
        three_ops_term = get_additive_terms(a * a' * a)[1]  # Get the (a′*a*a) term
        @test check_linear(three_ops_term) == false
        
        # Test 4: Parameter times operator (linear)
        param_op = g * a
        @test check_linear(param_op) == true
        
        # Test 5: Parameter times two operators (linear)
        param_two_ops_term = get_additive_terms(g * a * a')[1]  # Get the multiplicative term
        @test check_linear(param_two_ops_term) == true
        
        # Test 6: Parameter times three operators (nonlinear)
        param_three_ops_term = get_additive_terms(g * a * a' * a)[1]  # Get the multiplicative term
        @test check_linear(param_three_ops_term) == false
        
        # Test 7: Just a number (linear)
        just_number = 5
        @test check_linear(just_number) == true
        
        # Test 8: Just a parameter (linear)
        just_param = g
        @test check_linear(just_param) == true
    end
    
    @testset "count_quantum_operators" begin
        # Create test Hilbert space and operators
        hf = FockSpace(:cavity)
        @qnumbers a::Destroy(hf)
        
        # Create symbolic parameters
        @cnumbers g ω κ
        
        # Test 1: Single operator
        @test count_quantum_operators(a) == 1
        @test count_quantum_operators(a') == 1
        
        # Test 2: Two operators in multiplicative term
        two_ops_term = get_additive_terms(a * a')[1]  # Get the (a′*a) term
        @test count_quantum_operators(two_ops_term) == 2
        
        # Test 3: Three operators in multiplicative term
        three_ops_term = get_additive_terms(a * a' * a)[1]  # Get the (a′*a*a) term
        @test count_quantum_operators(three_ops_term) == 3
        
        # Test 4: Parameter times operator
        @test count_quantum_operators(g * a) == 1
        
        # Test 5: Parameter times two operators in multiplicative term
        param_two_ops_term = get_additive_terms(g * a * a')[1]  # Get the multiplicative term
        @test count_quantum_operators(param_two_ops_term) == 2
        
        # Test 6: Just a number
        @test count_quantum_operators(5) == 0
        
        # Test 7: Just a parameter
        @test count_quantum_operators(g) == 0
        
        # Test 8: Complex parameter expression
        @test count_quantum_operators(g * ω + κ) == 0
    end
    
    @testset "get_additive_terms and check_linear integration" begin
        # Create test Hilbert space and operators
        hf = FockSpace(:cavity)
        @qnumbers a::Destroy(hf)
        
        # Create symbolic parameters
        @cnumbers g ω κ
        
        # Test 1: Linear Hamiltonian with additive terms
        H_linear = g * a + ω * a' + κ * a' * a
        terms = get_additive_terms(H_linear)
        @test length(terms) == 3
        @test all(check_linear.(terms))
        
        # Test 2: Nonlinear Hamiltonian with one nonlinear term
        H_nonlinear = g * a + ω * a' + κ * a' * a * a
        terms = get_additive_terms(H_nonlinear)
        @test length(terms) == 3
        linearity_results = check_linear.(terms)
        @test sum(linearity_results) == 2  # Two linear terms, one nonlinear
        @test !linearity_results[3]  # The third term (a' * a * a) should be nonlinear
        
        # Test 3: Mixed linear and nonlinear terms
        H_mixed = g * a' * a + ω * a' * a * a + κ * a
        terms = get_additive_terms(H_mixed)
        @test length(terms) == 3
        linearity_results = check_linear.(terms)
        @test sum(linearity_results) == 2  # Two linear terms, one nonlinear
        
        # Test 4: All linear terms
        H_all_linear = g * a + ω * a' + κ * a' * a + (g + ω) * a
        terms = get_additive_terms(H_all_linear)
        @test all(check_linear.(terms))
        
        # Test 5: All nonlinear terms (use terms that won't simplify)
        H_all_nonlinear = g * a * a * a + ω * a' * a' * a' + κ * a' * a * a
        terms = get_additive_terms(H_all_nonlinear)
        @test all(.!check_linear.(terms))  # All terms should be nonlinear
    end
    
    @testset "ordered_ops" begin
        # Test 1: Single mode system
        hf = FockSpace(:cavity)
        @qnumbers a::Destroy(hf)
        @cnumbers ω κ
        
        H_single = ω * a' * a
        L_single = [√κ * a]
        S_single = [1]
        sys_single = SLH(:single, [:in], [:out], S_single, L_single, H_single)
        
        ops = ordered_ops(sys_single)
        @test length(ops) == 2
        @test ops[1] isa Destroy  # annihilation operator first
        @test ops[2] isa Create   # creation operator second
        @test ops[1].name == :a
        @test ops[2].name == :a
        
        # Test 2: Two mode system (optomechanical from slh_to_abcd.md)
        hilb = FockSpace(:optical) ⊗ FockSpace(:mechanical)
        a_opt = Destroy(hilb, :a, 1)
        b_mech = Destroy(hilb, :b, 2)
        @cnumbers ω Δ g κ
        
        H_optomech = ω * b_mech' * b_mech + Δ * a_opt' * a_opt + g * (a_opt' + a_opt) * (b_mech' + b_mech)
        L_optomech = [√κ * a_opt]
        S_optomech = [1]
        sys_optomech = SLH(:optomech, [:in], [:out], S_optomech, L_optomech, H_optomech)
        
        ops = ordered_ops(sys_optomech)
        @test length(ops) == 4  # 2 modes × 2 operators per mode
        
        # Check ordering: [a, a†, b, b†] (sorted by name)
        @test ops[1] isa Destroy && ops[1].name == :a
        @test ops[2] isa Create && ops[2].name == :a
        @test ops[3] isa Destroy && ops[3].name == :b
        @test ops[4] isa Create && ops[4].name == :b
        
        # Test 3: System with only coupling operators (no H)
        H_empty = 0
        L_test = [√κ * a_opt, √κ * b_mech]
        S_test = [1]
        sys_coupling = SLH(:coupling, [:in], [:out], S_test, L_test, H_empty)
        
        ops = ordered_ops(sys_coupling)
        @test length(ops) == 2  # Only annihilation operators from coupling
        @test ops[1] isa Destroy && ops[1].name == :a
        @test ops[2] isa Destroy && ops[2].name == :b
        
        # Test 4: System with mixed operators in H and L
        hf1 = FockSpace(:mode1)
        hf2 = FockSpace(:mode2)
        hilb_mixed = hf1 ⊗ hf2
        c = Destroy(hilb_mixed, :c, 1)
        d = Destroy(hilb_mixed, :d, 2)
        
        H_mixed = ω * c' * c + g * c' * d  # c appears in both H and L
        L_mixed = [√κ * c, √κ * d]
        S_mixed = [1]
        sys_mixed = SLH(:mixed, [:in], [:out], S_mixed, L_mixed, H_mixed)
        
        ops = ordered_ops(sys_mixed)
        @test length(ops) == 3  # c, c', d (no d' operator in system)
        # Check that operators are properly deduplicated and sorted
        @test ops[1] isa Destroy && ops[1].name == :c
        @test ops[2] isa Create && ops[2].name == :c
        @test ops[3] isa Destroy && ops[3].name == :d
        
        # Test 5: Empty system (no operators)
        H_zero = 0
        L_zero = [0]
        S_zero = [1]
        sys_zero = SLH(:zero, [:in], [:out], S_zero, L_zero, H_zero)
        
        ops = ordered_ops(sys_zero)
        @test length(ops) == 0
    end
    
end
