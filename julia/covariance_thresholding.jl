using LinearAlgebra

#Soft-thresholding operators
function soft_thresholding(a::Real; τ::Real=0.1)
    if a > τ
        return a - τ
    elseif a < -τ
        return a + τ
    else
        return 0
    end
end
function soft_thresholding(A::Union{Array,Matrix}; τ::Real=0.1)
    AS = zeros(size(A))
    for i in eachindex(A)
        AS[i] = soft_thresholding(A[i], τ=τ)
    end
    return AS
end

function covarianceThresholding_adaptive(S1, r::Int, klist::Array{Int}; ν::Real = 0.6745, νp::Real = 2)
    p = size(S1,1)
    #Scaling for the thresholding parameter based on problem characteristics
    ν0 = if p > sum(klist)^2 
        sqrt(log(p/(sum(klist)^2)))
    elseif p > mean(klist)^2 
        sqrt(log(p/(mean(klist)^2)))
    else
        sqrt(log(p/(mean(klist))))
    end

    α_grid = collect(1.0:0.2:3.0)
    n_grid = p .* (2.0.^collect(-2:1:5))
    best_a = 0; best_n = 0; best_obj = -Inf
    for a in α_grid
        for n in n_grid
            τ = a*ν0/sqrt(n)
            η1 = soft_thresholding(S1; τ=τ)

            F =  svd(η1, full=true) #SVD on soft-thresholded version
            V = F.U[:,1:r]

            for t in 1:r

                support_PC = sortperm(abs.(V[:,t]), rev=true)[1:klist[t]] #Tie-robust implementation
                V[setdiff(1:p,support_PC),t] .= 0 
                if norm(V[:,t]) > 0
                    V[:,t] ./= norm(V[:,t])
                end

                @assert sum(abs.(V[:,t]) .> 0) <= klist[t]
            end

            obj = dot(S1*V,V) / tr(S1)
            if obj > best_obj 
                best_obj = obj
                best_a = a
                best_n = n
            end
        end
    end

    τ = best_a*ν0/sqrt(best_n)
    η1 = soft_thresholding(S1; τ=τ)

    F =  svd(η1, full=true) #SVD on soft-thresholded version
    V = F.U[:,1:r]

    for t in 1:r

        support_PC = sortperm(abs.(V[:,t]), rev=true)[1:klist[t]] #Tie-robust implementation
        V[setdiff(1:p,support_PC),t] .= 0 

        if norm(V[:,t]) > 0
            V[:,t] ./= norm(V[:,t])
        end

        @assert sum(abs.(V[:,t]) .> 0) <= klist[t]

    end

    return V

end

function covarianceThresholding(S1, r::Int, n::Int, klist::Array{Int}; ν::Real = 0.6745, νp::Real = 2)
    p = size(S1,1)
    τ = if p > sum(klist)^2 
        νp*sqrt(log(p/(sum(klist)^2))/n)
    elseif p > mean(klist)^2 
        νp*sqrt(log(p/(mean(klist)^2))/n)
    else
        νp*sqrt(log(p/(mean(klist)))/n)
    end
    η1 = soft_thresholding(S1; τ=τ)

    F =  svd(η1, full=true) #SVD on soft-thresholded version
    V = F.U[:,1:r]

    for t in 1:r
        τ = sort(abs.(V[:,t]), rev=true)[klist[t]]
        
        V[:,t] .*= (abs.(V[:,t]) .>= τ)
        if norm(V[:,t]) > 0
            V[:,t] ./= norm(V[:,t])
        end
    end

    return V
end

