function ols_regression(X::Matrix, y::Vector)
    # Add a column of ones to X to account for the intercept
    X = hcat(ones(size(X, 1)), X)

    # Compute the OLS solution: (X'X)^(-1) X'y
    β = inv(X' * X) * X' * y

    # Predicted values
    ŷ = X * β

    # Residuals
    residuals = y - ŷ

    # Return regression coefficients, predicted values, and residuals
    return β, ŷ, residuals
end


function get_slope(nretics::Vector{Int}, neg_plls::Vector{Float64})
    X = zeros(length(nretics), 1)
    X .= nretics
    ols_coefs, _, _ = ols_regression(X, neg_plls)
    return ols_coefs
end


function select_model(nretics::Vector{Int}, nplls::Vector{Float64})
    all(nplls .<= 0) || error("All negative log-liks should be <= 0")
    length(nretics) == length(nplls) || error("Supplied $(length(nretics)) retic numbers but $(length(nplls)) -logliks.")

    ols_int, ols_slope = get_slope(nretics, nplls)
    penalized_scores = [.-nplls[i] + 2*ols_slope*nretics[i] for i = 1:length(nretics)]
    _, _h_star = findmin(penalized_scores)
    return nretics[_h_star], penalized_scores, ols_int, ols_slope
end


function plot_model_selection(h::AbstractVector{Int}, nlls::Vector{Float64})
    for j = 2:length(nlls)
        nlls[j] = max(nlls[j-1], nlls[j])
    end

    best_h, h_scores, ols_int, ols_slope = select_model(h, nlls)
    p = plot(h, nlls, label="Neg. Loglik")
    plot!(p, h, .-h_scores, label="Neg. Pen. Scores")
    vline!(p, [best_h], label="Best h", lc=:red, alpha=0.25, linestyle=:dash)
    return p
end

