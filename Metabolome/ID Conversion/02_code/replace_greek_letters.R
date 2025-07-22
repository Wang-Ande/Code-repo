replace_greek_letters <- function(x) {
  greek_map <- c(
    "α" = "alpha", "β" = "beta", "γ" = "gamma", "δ" = "delta",
    "ε" = "epsilon", "ζ" = "zeta", "η" = "eta", "θ" = "theta",
    "κ" = "kappa", "λ" = "lambda", "μ" = "mu", "ν" = "nu",
    "ξ" = "xi", "π" = "pi", "ρ" = "rho", "σ" = "sigma", "ς" = "sigma",
    "τ" = "tau", "φ" = "phi", "χ" = "chi", "ψ" = "psi", "ω" = "omega"
  )
  
  for (greek in names(greek_map)) {
    x <- gsub(greek, greek_map[[greek]], x)
  }
  return(x)
}

