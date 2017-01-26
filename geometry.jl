typealias Float Float64
    
immutable Point
    x::Float
    y::Float
end


Base.hash(p::Point, h::UInt) = hash(hash(p.x, hash(p.y)), h)

Base.isequal(p::Point, q::Point) = isequal(p.x, q.x) && isequal(p.y, q.y)

Base.isapprox(p::Point, q::Point; rtol::Real = sqrt(eps(Float)), atol::Real = 0) = 
    isapprox(p.x, q.x; rtol = rtol, atol = atol) && isapprox(p.y, q.y; rtol = rtol, atol = atol)

Base.show(io::IO, p::Point) = Base.print(io, "($(p.x), $(p.y))")


Base.:+(p::Point, q::Point)::Point = Point(p.x + q.x, p.y + q.y)
Base.:-(p::Point, q::Point)::Point = Point(p.x - q.x, p.y - q.y)
Base.:*(alpha::Real, q::Point)::Point = Point(alpha * q.x, alpha * q.y)
Base.:*(p::Point, α::Real)::Point = α * p
Base.:/(p::Point, α::Real)::Point = inv(α) * p

Base.norm(p::Point)::Float = sqrt(p.x ^ 2 + p.y ^ 2)

distance(p::Point, q::Point)::Float = norm(p - q)



iscollinear(::Point) = true

iscollinear(::Point, ::Point) = true

function iscollinear(x::Point, y::Point, z::Point)::Bool
    # http://math.stackexchange.com/a/59248/31127
    M = ones(Float, 3, 3)
    for (i, p) in enumerate([x, y, z])
        M[i, 2:3] = [p.x, p.y]
    end
    
    area = det(M) / 2
    area ≈ 0
end

"""
    iscollinear(ps::Point...)

Determine if all points in `ps` are collinear.  No point counts as collinear.
"""
function iscollinear(ps::Point...)::Bool
    triplets = zip(ps, ps[2:end], ps[3:end])
    all(iscollinear(x, y, z) for (x, y, z) in triplets)
end



immutable Rectangle
    lower_left::Point
    upper_right::Point
end
