function rouwenhorst(μ,ρ,σ,n)
    n==1 && (return [μ],[1.0])

	mu_eps=0
	q = (ρ+1)/2
	nu = ((n-1)/(1-ρ^2))^(1/2) * σ
	P = [q 1-q; 1-q q]
	for i = 2:n-1
		P = q*[P zeros(i,1);zeros(1,i+1)] + (1-q)*[zeros(i,1) P; zeros(1,i+1)]+ (1-q)*[zeros(1,i+1); P zeros(i,1)] + q*[zeros(1,i+1); zeros(i,1) P]
		P[2:i,:] = P[2:i,:]/2
	end
	x = [linspace(mu_eps/(1-ρ).-nu,mu_eps/(1-ρ).+nu,n);].+μ
	return x,P
end

function cdf_normal(x) :inline
    c = 0.5 * erfc(-x/sqrt(2))
end

function tauchen(μ,ρ,σ,N,m=3)
	if N==1
		return [μ],[1.0]
	end
	Z     = zeros(N)
	Zprob = zeros(N,N)
	a     = (1-ρ)*μ

	Z[N]  = m * sqrt(σ^2 / (1 - ρ^2))
	Z[1]  = -Z[N]
	zstep = (Z[N] - Z[1]) / (N - 1)

	for i=2:(N-1)
	    Z[i] = Z[1] + zstep * (i - 1)
	end

	Z = Z .+ a / (1-ρ)

	for j = 1:N
	    for k = 1:N
	        if k == 1
	            Zprob[j,k] = cdf_normal((Z[1] - a - ρ * Z[j] + zstep / 2) / σ)
	        elseif k == N
	            Zprob[j,k] = 1 - cdf_normal((Z[N] - a - ρ * Z[j] - zstep / 2) / σ)
	        else
	            Zprob[j,k] = cdf_normal((Z[k] - a - ρ * Z[j] + zstep / 2) / σ) - cdf_normal((Z[k] - a - ρ * Z[j] - zstep / 2) / σ)
	        end
	    end
	end
	return Z,Zprob
end
