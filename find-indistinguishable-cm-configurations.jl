using Combinatorics: combinations

# -------------------------------------------

nextpermutation(swap::Function) = ([1], 1)

@inline function nextpermutation(swap::Function, val, state)
    c, i = state

    while true
        if  c[i] < i
            if isodd(i)
                val = swap(1, i, val)
            else
                val = swap(c[i], i, val)
            end

            c[i] += 1
            i = 1

            return val, (c, i)
        else
            c[i] = 1
            i += 1
            while i > length(c)
                push!(c, 1)
            end
        end
    end
end

nextpermutation!(a::AbstractVector, args...) = nextpermutation((i,j) -> ((a[i], a[j]) = (a[j], a[i])), args...)

function for_permutations_with_form(f, a, form)
    k = size(form, 2)
    for b in combinations(a, k)
        @inline function swap(i, j, val)
            val += (form[1, j] - form[1, i]) * (b[i] - b[j])
            b[i], b[j] = b[j], b[i]
            return val
        end

        val = form * b
        state = nextpermutation(swap)
        f(b, val)
        for _ in 2:factorial(k)
            val, state = nextpermutation(swap, val, state)
            f(b, val)
        end
    end
end

# -------------------------------------------

using RootSystems

function RootSystems.coefficients_on_simple_roots(Φ::RootSystem, ::Val{:TOP})
    return sum(coefficients_on_simple_roots(Φ, α) for α in positive_roots(Φ))
end

function extend_simple_root_values_to_positive_roots(Φ::RootSystem, simple_root_values)
    @assert length(simple_root_values) == rank(Φ)

    map(enumerate(positive_roots(Φ))) do I
        i, α = I
        sum(c * v for (c, v) in zip(coefficients_on_simple_roots(Φ, α), simple_root_values))
    end
end

@doc """
Try to find an additive assignment of positive values to the positive
roots in `root_system`, which have more than the expected number (2 in
most cases) of ways of assigning these numbers.
"""
function find_indistinguishable_configurations(Φ::RootSystem, searchrange)
    expected_number = length(dynkin_diagram_automorphisms(Φ))
    top_coeffs = coefficients_on_simple_roots(Φ, Val(:TOP))

    form = transpose(top_coeffs)

    #for simple_values in ((1, 1, 4, 2, 3),)
    #for simple_values in ((5, 4, 3, 11, 2),)
    #for simple_values in ((7, 10, 9, 3, 2),)
    #for simple_values in ((1, 1, 1, 4, 3, 3),)
    #for simple_values in ((3, 2, 1, 3, 4, 1),)
    while true
        simple_values = ntuple(_ -> rand(searchrange), Val(rank(Φ)))


        all_values = extend_simple_root_values_to_positive_roots(Φ, simple_values)

        max_possible_simple_value = begin
            top_coeffs_sorted = sort(top_coeffs)
            v = sort(all_values)

            min_coeff, coeffs = top_coeffs_sorted[1], top_coeffs_sorted[end-1:-1:1]

            min_total_rest = sum(a*b for (a,b) in zip(coeffs, @view v[1:length(coeffs)]))

            (sum(all_values) - min_total_rest)//min_coeff
        end

        subset_that_can_be_assigned_to_simple_root = filter(v -> v <= max_possible_simple_value, all_values)
        number_of_loops = factorial(big(length(subset_that_can_be_assigned_to_simple_root))) ÷ factorial(big(length(subset_that_can_be_assigned_to_simple_root) - rank(Φ) + 1))
        if number_of_loops * 100 > 30E9
            continue
        end
        #@info "Trying configuration" simple_values number_of_loops

        possible_assignments = Set()
        sorted_all_values = sort(all_values)
        total_of_values = sum(all_values)

        partialform, f = form[:,1:end-1], form[1,end]
        for_permutations_with_form(subset_that_can_be_assigned_to_simple_root, partialform) do other_assignment, other_total_of_values
            delta = total_of_values - other_total_of_values

            last_assignment, remainder = divrem(delta, f)
            iszero(remainder) || return
            last_assignment in subset_that_can_be_assigned_to_simple_root || return

            other_assignment = [other_assignment; last_assignment]

            other_all_values = extend_simple_root_values_to_positive_roots(Φ, other_assignment)

            if sorted_all_values == sort(other_all_values)
                push!(possible_assignments, copy(other_assignment))
            end
        end
        if length(possible_assignments) > expected_number
            @show possible_assignments
            @show sort(all_values)
        end
    end
end

find_indistinguishable_configurations(D(5), 1:10)
