
"""Assert that `select` does the right thing by using it to sort `trials` sequences."""
function test_select(;trials = 100, size = 20)
    @assert all(1:trials) do nix
        x = rand(1:10, size)
        sort(x) == [select(x, k) for k in 1:size]
    end
end

"""Generate `samples` uniform random points, and calculate the Tukey depth of another one."""
function test_tukey(;samples = 10, showplot = false)
    xs, ys = randn(Float, samples), randn(Float, samples)
    testpoint = Point(rand(Float), rand(Float))

    if showplot
        scatter(xs, ys, color = "r")
        scatter(testpoint.x, testpoint.y, color = "g")
    end
    
    tukeydepth(testpoint, map(Point, xs, ys))
end

"""Generate `samples` random points on a circle, and calculate the Tukey depth of the center."""
function test_tukey_circle(;samples = 10, showplot = false)
    # generate points on a circle, and test one point within
    const dist = MvNormal(eye(2))
    coords = mapslices(normalize, rand(dist, samples), 1)
    xs, ys = coords[1, :], coords[2, :]
    testpoint = Point(0, 0)

    if showplot
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


"""Calculate the internal tangent intersection of a random set of points, randomly divided"""
function test_tangent_intersection(;samples = 10, showplot = false)
    @assert samples > 4
    xs, ys = sort(rand(Float, samples)), rand(Float, samples)
    nₗ = rand(2:(samples-2))
    left = map(Point, xs[1:nₗ], ys[1:nₗ])
    right = map(Point, xs[nₗ+1:end], ys[nₗ+1:end])
    
    l1 = get(internal_tangent(left, right, :Max))
    l2 = get(internal_tangent(left, right, :Min))
    p = get(intersect_lines(l1, l2))
    
    if showplot
        xlim(-0.5, 1.5)
        ylim(-0.5, 1.5)

        scatter(xs[1:nₗ], ys[1:nₗ], color = "r")
        scatter(xs[nₗ+1:end], ys[nₗ+1:end], color = "y")
        scatter(p.x, p.y, color = "b")

        xmin, xmax, ymin, ymax = axis()
        xaxis = linspace(xmin, xmax, 100)
        plot(xaxis, l1[1] * xaxis + l1[2], color = "k")
        plot(xaxis, l2[1] * xaxis + l2[2], color = "k")
    end

    return p
end


"""Calculate and animate the plurality point of a random set."""
function test_plurality_points1(;samples = 4)
    points = map(Point, rand(Float64, samples), rand(Float64, samples))
    plurality_points_animated(points)
    plurality_points(points)
end


function test_plurality_points2()
    plurality_points_animated([Point(0.0, 0.0), Point(-1.0, 0.0), Point(-0.5, 0.0),
                               Point(1.0, 0.0), Point(0.0, -1.0), Point(0.0, 1.0)])
end


"""Generate random point sets until one possessing a plurality point is found, and animate
that. """
function findexample(;samples = 4, maxiterations = 100)
    for i = 1:maxiterations
        testpoints = map(Point, rand(Float64, samples), rand(Float64, samples))
        pp = plurality_points(testpoints)

        if !isempty(pp)
            info("Example found after ", i, " iterations")
            plurality_points_animated(testpoints)
            return pp
        end
    end

    info("No example found after ", maxiterations, " iterations")
end
