# Collapsed annotation-alpha inference

Status: research prototype; not a maintained `sblr` model.

This workspace defines and qualifies the continuous-alpha collapsed target before any production implementation. It makes no novelty claim and does not alter standard SBayesRC or the separately specified SBayesRC-S selection model.

## Oracle provenance

The computational oracle originated as the inbox file `research/local_reference/sbayesrc_alpha_solution_v5.py`. It was copied to `python/sbayesrc_alpha_solution_v5.py`, then its shared probit numerical core was corrected after independent review so the evaluated target and gradient both use exact stable tail probabilities. The ignored inbox copy remains unchanged.

- Original inbox Python v5 SHA-256: `e3b4a12f29bbf3c4a314a164773c071b583db93e1c54132fb33a4252e7ad1332`
- Corrected research-copy SHA-256: `ec378f802b9b6aea7d856052bce312f9089890b55e370fae45c8f9e317ceaf64`
- Mathematical findings: `research/local_reference/sbayesrc_annotation_findings.md`
- Findings SHA-256: `e9e01a0c635971859d911286434ea962883d7219d87de258474aa90c47aaceb8`

The oracle uses NumPy, SciPy, pandas, and matplotlib. It is qualification evidence, not an `sblr` runtime dependency.

## Statistical model

For markers $j=1,\ldots,M$, the observed independent-summary model is

$$
\widehat\beta_j\mid\beta_j\sim
N(\beta_j,\mathrm{SE}_j^2).
$$

There are $K$ ordered mixture components. Conditional on allocation $z_j=k$,

$$
\beta_j\mid z_j=k\sim N(0,V_\beta\gamma_k),
\qquad k=0,\ldots,K-1.
$$

The null convention is $\gamma_0=0$; every active $\gamma_k$ is positive. This prototype conditions on $V_\beta$, $\gamma$, and every $\mathrm{SE}_j$. It does not learn variance parameters.

### Annotation design and continuation link

Let $x_j\in\mathbb R^P$ contain the non-intercept annotations and let $c_s\in\mathbb R^P$ be a fixed centre for stick $s=1,\ldots,K-1$. The continuation-centred coefficient is

$$
\theta_s=(\theta_{0s},\theta_{1s},\ldots,\theta_{Ps})^\mathsf T,
$$

with

$$
\eta_{js}=\theta_{0s}+(x_j-c_s)^\mathsf T\theta_{-0,s},
\qquad q_{js}=\Phi(\eta_{js}).
$$

The v5 oracle uses the probit link $\Phi$, not the logistic link. Its simulated continuous and binary columns are both centred and divided by their sample standard deviation before the intercept is added. The R research fixtures use the same convention. This differs from some historical production SBayesRC preprocessing profiles and is part of this prototype's declared model.

Centres are estimated once, then held fixed. The oracle fits an intercept-only collapsed model, forms posterior component responsibilities, and uses the implied probability of reaching each stick as weights for annotation means. This is a coordinate transformation, not an iteration-dependent target change.

The corresponding raw coefficient matrix $\alpha$ satisfies

$$
\alpha_{-0,s}=\theta_{-0,s},\qquad
\alpha_{0s}=\theta_{0s}-c_s^\mathsf T\theta_{-0,s}.
$$

Thus $[1,x_j^\mathsf T]\alpha_s=[1,(x_j-c_s)^\mathsf T]\theta_s$ exactly.

### Stick-breaking probabilities

The component order is null first and then increasing active mixture scale. With $K=4$,

$$
\begin{aligned}
\pi_{j0}&=1-q_{j1},\\
\pi_{j1}&=q_{j1}(1-q_{j2}),\\
\pi_{j2}&=q_{j1}q_{j2}(1-q_{j3}),\\
\pi_{j3}&=q_{j1}q_{j2}q_{j3}.
\end{aligned}
$$

For general $K$, component $k<K-1$ is the probability of continuing through sticks $1,\ldots,k$ and stopping at $k+1$; the last component continues through every stick. The scientific target is not epsilon-clipped. Likelihood calculations use the exact stable representations

$$
\log q_{js}=\log\Phi(\eta_{js}),\qquad
\log(1-q_{js})=\log\Phi(-\eta_{js}),
$$

evaluated with log-CDF and log-survival functions. Component log probabilities and log weights are constructed directly from these quantities. Probabilities produced only for reporting may round to exactly zero or one in finite-precision arithmetic while the log-space target remains usable.

### Intercept and global mixture

The intercept is not independent of a declared global mixture. For a marker population,

$$
\bar\pi_k=M^{-1}\sum_j\pi_{jk}(\alpha).
$$

Changing slopes, annotation prevalence, or annotation correlation generally changes $\bar\pi$ unless intercepts change too. The oracle simulation calibrates each raw intercept sequentially so the eligibility-weighted mean continuation equals the configured target $(0.12,0.50,0.45)$. In inference, the centred intercept prior mean is the probit transform of that target. This prior centre is not a constraint forcing the posterior global mixture to remain fixed.

## Exact collapsed likelihood

Integrating $\beta_j$ within each component gives

$$
\ell_{jk}=\log\phi\!\left(
\widehat\beta_j;0,\mathrm{SE}_j^2+V_\beta\gamma_k
\right).
$$

Integrating the allocation gives the exact marker likelihood

$$
p(\widehat\beta_j\mid\alpha,\Theta)=
\sum_{k=0}^{K-1}\pi_{jk}(\alpha)
\phi\!\left(
\widehat\beta_j;0,\mathrm{SE}_j^2+V_\beta\gamma_k
\right),
$$

where $\Theta=(V_\beta,\gamma,\mathrm{SE},c)$ is fixed here. Marker contributions are evaluated with log-sum-exp.

### Prior and log posterior

The centred coefficients have independent proper Gaussian priors,

$$
\theta_{sr}\sim N(m_{sr},\tau_{sr}^2).
$$

The complete log posterior is

$$
\log p(\theta\mid\widehat\beta,\Theta)=
\sum_j\log\sum_k\exp\{\log\pi_{jk}(\theta)+\ell_{jk}\}
-\frac12\sum_{s,r}\left[
\frac{(\theta_{sr}-m_{sr})^2}{\tau_{sr}^2}
+\log(2\pi\tau_{sr}^2)
\right].
$$

The copied Python class reports the same expression with the prior normalizing constant omitted. The R reference exposes both the normalized posterior and this oracle-compatible kernel.

## Analytic gradient

Define normalized component responsibilities

$$
r_{jk}=\frac{\pi_{jk}\exp(\ell_{jk})}
{\sum_h\pi_{jh}\exp(\ell_{jh})}.
$$

At stick $s$, let

$$
R^+_{js}=\sum_{k=s}^{K-1}r_{jk},\qquad
R^-_{js}=r_{j,s-1},
$$

where zero-based components make $R^-_{j1}=r_{j0}$. With centred design row $d_{js}=(1,(x_j-c_s)^\mathsf T)^\mathsf T$,

$$
\nabla_{\theta_s}\log p(\theta\mid-)=
\sum_j d_{js}\phi(\eta_{js})
\left\{
\frac{R^+_{js}}{q_{js}}-
\frac{R^-_{js}}{1-q_{js}}
\right\}
-\operatorname{diag}(\tau_s^{-2})(\theta_s-m_s).
$$

The gradient is a $(K-1)\times(P+1)$ matrix and is flattened stick first, coefficient second. Positive score means increasing continuation at that stick improves the posterior. The tail derivatives are evaluated from the same log representation as the target,

$$
\frac{d\log q}{d\eta}=\exp\{\log\phi(\eta)-\log\Phi(\eta)\},\qquad
\frac{d\log(1-q)}{d\eta}=-\exp\{\log\phi(\eta)-\log\Phi(-\eta)\},
$$

and are independently checked with five-point central finite differences.

## Laplace geometry and whitening

The R prototype maximizes the collapsed log posterior with analytic gradients. Failure or a nonfinite optimum is an error. At the mode $\widehat\theta$, it evaluates the observed negative-log-posterior Hessian $G$. If $G=U\Lambda U^\mathsf T$, the regularized precision is

$$
G_\mathrm{reg}=U\operatorname{diag}{\max(\lambda_d,\lambda_*)\}U^\mathsf T,
$$

where $\lambda_*=\max(10^{-6},\lambda_{\max}/10^8)$. The number of altered eigenvalues, floor, and maximum added precision are reported. No failed metric is silently replaced.

Let $\Sigma=G_\mathrm{reg}^{-1}=LL^\mathsf T$ with lower Cholesky factor $L$. Whitening is

$$
z=L^{-1}(\theta-\widehat\theta),\qquad
\theta=\widehat\theta+Lz,
$$

and the whitened gradient is $L^\mathsf T\nabla_\theta\log p(\theta\mid-)$. The Python oracle uses a BFGS inverse-Hessian approximation with eigenvalues clipped to $[10^{-6},10^3]$; optimizer-level metric differences are recorded separately from same-metric parity.

## Metropolis-corrected HMC

The HMC position is the whitened $z\in\mathbb R^{(K-1)(P+1)}$ and momentum is $p\sim N(0,I)$. For step size $\epsilon$ and fixed leapfrog length $L_\mathrm{frog}$,

$$
\begin{aligned}
p&\leftarrow p+\tfrac12\epsilon\nabla_z\log p(z),\\
z&\leftarrow z+\epsilon p,\\
p&\leftarrow p+\epsilon\nabla_z\log p(z)
\end{aligned}
$$

with the final momentum update halved. The Hamiltonian is

$$
H(z,p)=-\log p(z\mid-)+\tfrac12p^\mathsf Tp.
$$

The proposal is accepted with

$$
\min[1,\exp\{H(z,p)-H(z',p')\}],
$$

so both posterior and kinetic-energy changes enter the decision. Nonfinite proposals are rejected. Warmup adapts only $\epsilon$ toward a declared acceptance target; leapfrog length remains fixed and sampling uses the frozen final step size. Reported diagnostics are acceptance rate/probability, invalid proposals, mean and maximum absolute energy error, split R-hat, and a compact autocorrelation ESS. The accepted state, not the proposal, is retained.

## Future partially collapsed production transition

The following is a future transition requiring a separate invariance proof and qualification; it is not implemented here:

1. condition on current mixture scales, variance parameters, and summary-likelihood state;
2. update continuous $\alpha$ from the exact collapsed conditional appropriate to that likelihood;
3. reconstruct marker probabilities from the accepted $\alpha$;
4. draw allocations and marker effects conditional on those probabilities;
5. update the remaining BLR parameters without conditioning on any stale value that was integrated out in step 2.

An LD-block version would require a newly derived block likelihood. The markerwise expression cannot simply be inserted into an appreciable-LD sampler and called exact.

## Distinct targets and algorithms

- **Continuous-alpha collapsed HMC:** integrates allocations and marker effects, then samples every continuous coefficient from the resulting smooth posterior.
- **Hard-allocation Gibbs:** conditions alpha regression on sampled allocations or latent continuation outcomes; it targets a joint augmented model but can suffer allocation-alpha feedback.
- **Frozen PIP/component regression:** fits alpha to plug-in posterior summaries and is modular/two-stage, not an exact joint update.
- **Soft-feedback ECM/MCEM:** alternates expected allocation information with optimization or approximate updates; it is not this posterior sampler.
- **Collapsed SBayesRC-S selection:** adds discrete coefficient-inclusion states and a different spike-and-slab posterior. It is explicitly outside this task.

## Scope and identifiability

The factorized likelihood is exact for independent, or effectively independent, summary statistics under the declared mixture model. With appreciable LD it is a composite likelihood. This workspace implements no LD-block likelihood and makes no claim for one.

Good HMC mixing does not imply that raw coefficients are identified. Limitations include:

- shrinking later-stick risk sets and event counts;
- rare annotations with few informative markers;
- correlated active annotations and interchangeable proxies;
- saturated continuation surfaces with small probit derivatives;
- dependence of raw intercepts/slopes on centering and scale.

Induced continuation probabilities, component probabilities, global mixture averages, and annotation contrasts may be better identified than individual raw $\alpha$ values and are retained as primary research diagnostics.

## Qualification reductions and invariants

The prototype checks:

1. zero slopes reduce to intercept-only marker probabilities;
2. $\alpha\leftrightarrow\theta$ and whitening transformations round-trip;
3. every continuation probability lies in $(0,1)$;
4. every component-probability and responsibility row sums to one;
5. component marginal densities use variance $\mathrm{SE}_j^2+V_\beta\gamma_k$, including the null reduction to $\mathrm{SE}_j^2$;
6. analytic and central finite-difference gradients agree at several interior points;
7. R and Python agree on deterministic likelihood, prior-kernel, gradient, same-metric whitening, leapfrog, Hamiltonian, and acceptance calculations;
8. leapfrog integration is reversible within floating-point tolerance;
9. the Metropolis ratio uses the complete Hamiltonian difference;
10. identical R seeds reproduce identical draws;
11. null, sparse-independent, correlated-proxy, later-stick, and two-sided tail fixtures remain finite and normalized;
12. marker-by-stick and marker-by-component axes, including their IDs, remain explicit when either the marker or stick dimension is one.
