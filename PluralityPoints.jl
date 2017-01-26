module PluralityPoints

# we need these for linear programming the tangent intersections
using JuMP
using GLPKMathProgInterface

# ENV["PYTHON"] = "/home/philipp/anaconda3/bin/python"
using PyPlot
using Distributions

import Base
include("geometry.jl")

export Point, Rectangle
export verify_candidates, internal_tangent_intersection, tukeydepth
export plurality_points
    

"""
    select(xs::Array, i::Integer, parts::Int = 5)::T

Select the `i`-th largest element from `xs` in linear time, using the algorithm
from "Cormen, T.H., Leiserson, C.E., Rivest, R.L., Stein, C.: Introduction to Algorithms".
`parts` determines the recursive splitting factor.
"""
function select{T}(xs::Vector{T}, k::Integer, parts::Int = 5)::T
    # see also: http://stackoverflow.com/a/28089259/1346276

    const n = length(xs)
    @assert 1 <= k <= n
    
    # if !(1 <= k <= n)
    #     @show k
    #     @show xs
    #     error("Invalid k")
    # end

    const median_median = cld(cld(n, parts), 2)

    if n <= parts
        return sort(xs)[k]
    end

    const groups = [xs[p:min(p + parts, end)] for p in 1:parts:length(xs)]
    const medians = [sort(g)[cld(length(g), 2)] for g in groups]
    const median_of_medians = select(medians, median_median, parts)
    
    const l = filter(s -> s < median_of_medians, xs)
    if k <= length(l)
        return select(l, k)
    end

    const m = filter(s -> s ≈ median_of_medians, xs)
    const r = filter(s -> s > median_of_medians, xs)
    if k > length(l) + length(m)
        return select(r, k - length(l) - length(m), parts)
    else
        return median_of_medians
    end
end


"""
    tukeydepth(θ::Point, points::Vector{Point})::Int

Calculate the Tukey depth of a point with respect to a point set.
From: Rousseuw, P.J., and Ruts, I. (1996), AS 307 : Bivariate location
depth, Applied Statistics (JRRSS-C), vol.45, 516-526
"""
function tukeydepth(θ::Point, points::Vector{Point})::Int
    # mercylessly translated from Fortran original; don't complain.
    
    const n = length(points)
    const ε = 0.00000001

    nt = 0
    α = zeros(Float, n)

    # construct α
    for i in 1:n
        d = distance(θ, points[i])

        if d <= ε
            nt += 1
        else
            u = (points[i] - θ) / d

            if abs(u.x) > abs(u.y)
                if points[i].x >= θ.x
                    α[i - nt] = asin(u.y)
                    if α[i - nt] < 0
                        α[i - nt] += 2π
                    end
                else
                    α[i - nt] = π - asin(u.y)
                end
            else
                if points[i].y >= θ.y
                    α[i - nt] = acos(u.x)
                else
                    α[i - nt] = 2π - acos(u.x)
                end
            end

            if α[i - nt] >= (2π - ε)
                α[i - nt] = 0
            end
        end
    end

    nn = n - nt

    if nn <= 1
        return nt
    end

    sort!(α)


    # check wheter θ lies outside points
    angle = α[1] - α[nn] + 2π

    for i in 2:nn
        angle = max(angle, α[i] - α[i - 1])
    end

    if angle > (π + ε)
        return nt
    end


    # compute nu, which is #(αs < pi)
    angle = α[1]
    nu = 0

    for i in 1:nn
        α[i] -= angle
        if α[i] < (pi - ε)
            nu += 1
        end
    end

    if nu >= nn
        return nt
    end
    

    # mergesort the α with antipodal angles; update i, f
    ja = 1
    jb = 1
    αₖ = α[1]
    βₖ = α[nu + 1] - pi
    i = nu
    nf = nn
    f = zeros(Float, n)

    for j in 1:2nn
        if (αₖ + ε) < βₖ
            nf += 1

            if ja < nn
                ja += 1
                αₖ = α[ja]
            else
                αₖ = 2π + 1
            end
        else
            i += 1

            if i == nn + 1
                i = 1
                nf -= nn
            end

            f[i] = nf

            if jb < nn
                jb += 1
                if (jb + nu) <= nn
                    βₖ = α[jb + nu] - π
                else
                    βₖ = α[jb + nu - nn] + π
                end
            else
                βₖ = 2π + 1
            end
        end
    end


    # compute numh
    gi = 0
    ja = 1
    angle = α[1]
    numh = min(f[1], nn - f[1])

    for i in 2:nn
        if α[i] <= (angle + ε)
            ja += 1
        else
            gi += ja
            ja = 1
            angle = α[i]
        end

        ki = f[i] - gi
        numh = min(numh, ki, nn - ki)
    end

    return numh + nt
end


function internal_tangent(A, B, direction::Symbol)::Nullable{Tuple{Float, Float}}
    # ATTENTION: A needs to be "left" of B!
    
    if direction == :Min
        A, B = B, A
    elseif direction != :Max
        error("Invalid direction")
    end
    
    m = Model(solver = GLPKSolverLP())
    
    @variable(m, k)
    @variable(m, d)
    
    for v in A
        @constraint(m, v.y <= k * v.x + d)
    end
    for v in B
        @constraint(m, v.y >= k * v.x + d)
    end

    @objective(m, direction, k)
    
    if solve(m) == :Optimal
        return Nullable((getvalue(k), getvalue(d)))
    else
        return Nullable()
    end
end


"""
    internal_tangent_intersection(A::Vector{Point}, B::Vector{Point})::Nullable{Point}

Compute intersection of the internal tangents of (disjoint and separable) point sets `A` and `B`.
"""
function intersect_lines(l1::Tuple{Float, Float}, l2::Tuple{Float, Float})::Nullable{Point}
    k₁, d₁ = l1
    k₂, d₂ = l2
    f = (k₁ - k₂)

    if f != 0
        x = (d₂ - d₁) / f
        y = k₁ * x + d₁
        return Nullable(Point(x, y))
    else
        return Nullable{Point}()
    end
end


"""
    internal_tangent_intersection(A::Vector{Point}, B::Vector{Point})::Nullable{Point}

Compute intersection of the internal tangents of (disjoint and separable) point sets `A` and `B`.
"""
function internal_tangent_intersection(A, B)::Nullable{Point}
    l1 = internal_tangent(A, B, :Max)
    l2 = internal_tangent(A, B, :Min)
    
    if !isnull(l1) && !isnull(l2)
        return intersect_lines(get(l1), get(l2))
    else
        return Nullable{Point}()
    end
end


verify_candidates(C, V) = filter(p -> tukeydepth(p, V) >= length(V) / 2, C)


function plurality_points(V::Vector{Point})
    n = length(V)

    V_x = [v.x for v in V]
    V_y = [v.y for v in V]
    x_h = select(V_x, cld(n + 1, 2))
    x_l = select(V_x, cld(n, 2))
    y_h = select(V_y, cld(n + 1, 2))
    y_l = select(V_y, cld(n, 2))
    
    if iscollinear(V...)
        return Rectangle(Point(x_l, y_l), Point(x_h, y_h))
    else
        C = Set([Point(x_h, y_h), Point(x_h, y_l), Point(x_l, y_h), Point(x_l, y_l)])
        
        if length(C) == 1
            return verify_candidates(C, V)
        else
            P = verify_candidates(intersect(V, C), V)
            
            if isempty(P)
                # not a corner ⇒ Δ ∉ V
                if x_h ≈ x_l
                    V_a = Set(v for v in V if v.y <= y_l)
                    V_b = Set(v for v in V if v.y >= y_h)
                else
                    V_a = Set(v for v in V if v.x <= x_l)
                    V_b = Set(v for v in V if v.x >= x_h)
                end
                
                p = internal_tangent_intersection(V_a, V_b)
                if isnull(p)
                    error("Something went wrong with intersection...")
                else
                    return verify_candidates(Set([get(p)]), V)
                end
            else
                # Δ ∈ V
                return P
            end
        end
    end
end


####################################################################################################
# TEST FUNCTIONS
####################################################################################################


function test_select(;trials = 100, size = 20)
    @assert all(1:trials) do nix
        x = rand(1:10, size)
        sort(x) == [select(x, k) for k in 1:size]
    end
end

function test_tukey(;samples = 10, plot = false)
    xs, ys = randn(Float, samples), randn(Float, samples)
    testpoint = Point(rand(Float), rand(Float))

    if plot
        scatter(xs, ys, color = "r")
        scatter(testpoint.x, testpoint.y, color = "g")
    end
    
    tukeydepth(testpoint, map(Point, xs, ys))
end

function test_tukey2(;samples = 10, plot = false)
    # generate points on a circle, and test one point within
    const dist = MvNormal(eye(2))
    coords = mapslices(normalize, rand(dist, samples), 1)
    xs, ys = coords[1, :], coords[2, :]
    testpoint = Point(0, 0)

    if plot
        scatter(xs, ys, color = "r")
        scatter(testpoint.x, testpoint.y, color = "g")
    end        
    
    tukeydepth(testpoint, map(Point, xs, ys))
end

# fun fact: if the direction of points is uniform, the tukey depth grows linearly with the
# number of points:
# ms = [mean(PluralityPoints.test_tukey2(samples = s) for _ in 1:100) for s in 1:100]
# ms2 = [mean(PluralityPoints.test_tukey(samples = s) for _ in 1:100) for s in 1:100]
# plot(1:100, ms2)


function test_tangent_intersection(;samples = 10, plot = false)
    xs, ys = sort(rand(Float, samples)), rand(Float, samples)
    i = get(internal_tangent_intersection(map(Point, xs[1:nₗ], ys[1:nₗ]),
                                          map(Point, xs[nₗ+1:end], ys[nₗ+1:end])))

    if plot
        scatter(xs[1:nₗ], ys[1:nₗ], color = "r")
        scatter(xs[nₗ+1:end], ys[nₗ+1:end], color = "y")
        scatter(i.x, i.y, color = "b")
    end

    return i
end


function plurality_points_animated(V::Vector{Point})
    n = length(V)

    V_x = [v.x for v in V]
    V_y = [v.y for v in V]
    x_h = select(V_x, cld(n + 1, 2))
    x_l = select(V_x, cld(n, 2))
    y_h = select(V_y, cld(n + 1, 2))
    y_l = select(V_y, cld(n, 2))
    
    scatter(V_x, V_y, marker = "o", color = "b")
    plot([x_l, x_h, x_h, x_l, x_l], [y_l, y_l, y_h, y_h, y_l], color = "k")
    
    if iscollinear(V...)
        info("Collinear")
    else
        C = Set([Point(x_h, y_h), Point(x_h, y_l), Point(x_l, y_h), Point(x_l, y_l)])

        if length(C) == 1
            unique_candidate = collect(C)[1]
            info("Unique median")
            info("Tukey depth $(tukeydepth(unique_candidate, V)), should be ≥ $(length(V) / 2)")
            scatter(unique_candidate.x, unique_candidate.y, marker = "x", color = "y")
            for p in verify_candidates(C, V)
                scatter(p.x, p.y, marker = "*", color = "r")
            end
        else
            P = verify_candidates(intersect(V, C), V)
            
            if isempty(P)
                info("Inside C")
                # not a corner ⇒ Δ ∉ V
                if x_h ≈ x_l
                    V_a = Set(v for v in V if v.y <= y_l)
                    V_b = Set(v for v in V if v.y >= y_h)
                    info("x_h ≈ x_l")
                else
                    V_a = Set(v for v in V if v.x <= x_l)
                    V_b = Set(v for v in V if v.x >= x_h)
                end
                
                l1 = get(internal_tangent(V_a, V_b, :Max))
                l2 = get(internal_tangent(V_a, V_b, :Min))
                candidate = get(intersect_lines(l1, l2))
                scatter(candidate.x, candidate.y, marker = "x", color = "y")
                
                info("Tukey depth $(tukeydepth(candidate, V)), should be ≥ $(length(V) / 2)")
                
                xmin, xmax, ymin, ymax = axis()
                xaxis = linspace(xmin, xmax, 100)
                plot(xaxis, l1[1] * xaxis + l1[2], color = "k")
                plot(xaxis, l2[1] * xaxis + l2[2], color = "k")
                
                for p in verify_candidates([candidate], V)
                    scatter(p.x, p.y, marker = "*", color = "r")
                end
            else
                # Δ ∈ V
                info("Median corner and in V")
                for p in P
                    scatter(p.x, p.y, marker = "*", color = "r")
                end
            end
        end
    end
end

function test_plurality_points1(;samples = 10)
    xs, ys = rand(Float64, samples), rand(Float64, samples)
    points = [Point(x, y) for (x, y) in zip(xs, ys)]
    xlim(0, 1)
    ylim(0, 1)
    plurality_points_animated(points)
    plurality_points(points)
end

function test_plurality_points2(;samples = 10, plot = false)
    xs, ys = rand(Float64, samples), rand(Float64, samples)
    points = [Point(x, y) for (x, y) in zip(xs, ys)]
    pp = plurality_points(points)

    if plot
        scatter(xs, ys, color = "r")
        scatter([p.x for p in pp], [p.y for p in pp], color = "g")
    end
    
    return pp
end

function test_plurality_points3()
    plurality_points_animated([Point(0.0, 0.0), Point(-1.0, 0.0), Point(-0.5, 0.0),
                               Point(1.0, 0.0), Point(0.0, -1.0), Point(0.0, 1.0)])
end





function findexample(samples = 5, maxiterations = 100)
    xs, ys = rand(Float64, samples), rand(Float64, samples)
    testpoints = [Point(x, y) for (x, y) in zip(xs, ys)]            
    pp = plurality_points(testpoints)
    i = 1
    
    while isempty(pp)
        i += 1
        xs, ys = rand(Float64, samples), rand(Float64, samples)
        testpoints = [Point(x, y) for (x, y) in zip(xs, ys)]            
        pp = plurality_points(testpoints)

        if i > maxiterations
            warn("No example found")
            return
        end
    end

    info("Example found after ", i, " iterations")
    pp = plurality_points(testpoints; debug = true)
    
    scatter(xs, ys, marker = "o", color = "b")
    scatter([p.x for p in pp], [p.y for p in pp], marker = "o", color = "r")
end


end # module
