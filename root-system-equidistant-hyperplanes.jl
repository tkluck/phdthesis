using RootSystems

normalized(h) = h .// h[findfirst(!iszero, h)]

function orbit(Φ, h)
    G = map(coordinates, dynkin_diagram_automorphisms(Φ))
    res = Set{Any}([h])
    while true
        newres = Set(normalized(reflectalong(coordinates(δ), r)) for δ in simple_roots(Φ) for r in res)
        union!(newres, Set(normalized(g * r) for r in res for g in G))
        prevlength = length(res)
        union!(res, newres)
        length(res) == prevlength && break
    end
    return res
end

function hyperplaneorbits(Φ)
    R = map(coordinates, roots(Φ))
    P = [(R[i], R[j]) for i in 1:length(R) for j in i+1:length(R)]
    H = Set(normalized(a-b) for (a,b) in P)
    orbits = []
    while !isempty(H)
        o = orbit(Φ, pop!(H))
        setdiff!(H, o)
        push!(orbits, o)
    end
    @info "Equidistant hyperplane orbits for $Φ" hyperplanes=sum(length, orbits) orbits=length(orbits)
    for o in orbits
        @info "Orbit" repres=tuple(first(o)...) length=length(o) multiplicity=count(((a,b),) -> (normalized(a-b) == first(o)), P) preimage=filter(((a,b),) -> (normalized(a-b) == first(o)), P)
    end
end

#hyperplaneorbits(A(2))
#hyperplaneorbits(A(3))
#hyperplaneorbits(A(4))
#hyperplaneorbits(A(5))

#hyperplaneorbits(D(5))
#hyperplaneorbits(D(6))
#hyperplaneorbits(D(7))
#
hyperplaneorbits(E(6))
