function reglsq( ;
                Γ = nothing,
                R = nothing,
                z = nothing,
                η0 = nothing,
                τ = nothing,
                tol = sqrt(τ^2 - η0^2),
                abstol = 1e-5,
                reltol = 1e-2,
                λinit = 1e-3,
                λmax = 1e2,
                verbose = true)
   @assert size(Γ, 2) == size(R, 2)
   @assert length(z) == size(R, 1)
   Qc, Rc = qr(Γ)
   Y = [y; zeros(size(Γ, 1))]
   verbose && @info("`reglsq` : solve regularised least squares")
   function _solve(λ)
      x = qr([ R; λ2 * Rc ]) \ Y
      return x, norm(R * x - z)
   end
   # Bracketing
   λ1 = 0.0
   λ2 = 1e-3
   while true
      x, η = _solve(λ2)
      if η > tol
         break
      end
      if λ2 > λmax
         @error("""reglsq reached λmax = $λmax in the bracketing stage.
                  The returned coeffcients are in the interior of the constraint
                  set ||Rx-z|| ≦ tol. If stronger regularisation is required,
                  then please increase the λmax parameter.""")
         return x
      end
      λ1 = λ2
      λ2 *= 10
   end
   verbose && @info("found bracket, starting bisection")
   # Bisection
   while true
      λ = 0.5 * (λ1 + λ2)
      x, η = _solve(λ)
      if tol / (1+reltol) <= η <= tol * (1+reltol)
         # success
         break
      end
      if η > tol
         λ2 = λ
      else
         λ1 = λ
      end
      if abs(λ1 - λ2) < abstol
         @error("""bracket has become too small; returning current solutions,
                   if it is unsatisfactory, then please reduce the
                   `abstol` parameter.""")
         return x
      end
   end
   verbose && @info("found a solution")
   return x
end