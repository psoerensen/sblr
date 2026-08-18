#!/usr/bin/env python3
"""SBayesRC annotation-alpha inference laboratory, version 5.

This script isolates the joint allocation--alpha problem identified by the
SBayesRC/SBayesRC-S studies.  It simulates mixed binary/continuous annotations
and four-component continuation-probit architectures, then compares:

  * hard_gibbs: the usual sampled-allocation -> probit-alpha update;
  * collapsed_ess: allocations and marker effects are integrated out of the
    alpha update and the resulting smooth posterior is sampled by elliptical
    slice sampling (ESS);
  * collapsed_hmc: the same collapsed posterior sampled after Laplace
    whitening by an exact Metropolis-corrected Hamiltonian update;
  * collapsed_lowrank: the collapsed likelihood with shared and later-stick
    annotation factors, which pools information across continuation sticks;
  * collapsed_selection: the same collapsed likelihood with exact
    annotation-by-stick spike/slab toggles and ESS coefficient updates.

For independent marker summaries, the collapsed likelihood is exact:

  b_j | alpha ~ sum_k pi_jk(alpha) N(0, se_j^2 + vb * gamma_k).

This makes the experiment a clean test of computation versus information.  A
robustness mode also generates major-polygenic and MAF-dependent effects; those
are deliberately misspecified relative to the fitted BayesR mixture.

Dependencies: numpy, scipy, pandas, matplotlib.  No R installation is needed.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import platform
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.optimize import brentq, minimize
from scipy.special import log_ndtr, logsumexp, ndtr, ndtri


EPS = 1e-12
GAMMA = np.array([0.0, 0.01, 0.10, 1.0], dtype=float)
STICK_NAMES = ("stick1", "stick2", "stick3")


@dataclass
class Config:
    preset: str = "smoke"
    seed: int = 20260821
    m: int = 240
    p_annotations: int = 12
    n_eff: int = 10000
    n_per_cell: int = 1
    annotation_scenarios: Tuple[str, ...] = (
        "null",
        "sparse_independent",
        "correlated_proxy",
        "later_stick_only",
    )
    genetic_scenarios: Tuple[str, ...] = ("bayesr_exact",)
    n_chains: int = 3
    n_iter: int = 1400
    burn: int = 600
    thin: int = 2
    slope_prior_sd: float = 0.75
    intercept_prior_sd: float = 1.50
    inclusion_prior: float = 0.15
    spike_sd: float = 0.05
    slab_sd: float = 0.75
    selection_toggle_every: int = 5
    correlation_swap_threshold: float = 0.85
    run_selection: bool = True
    run_lowrank: bool = False
    lowrank_factor_sd: float = 0.75
    run_hmc: bool = True
    hmc_leapfrog_steps: int = 12
    hmc_initial_step_size: float = 0.15
    hmc_target_accept: float = 0.75
    h2: float = 0.40
    target_continuation: Tuple[float, float, float] = (0.12, 0.50, 0.45)
    annotation_strength: float = 1.0
    binary_fraction: float = 0.40
    proxy_correlation: float = 0.95
    maf_s: float = -0.50
    major_share: float = 0.35
    max_slice_steps: int = 1000


@dataclass
class SimulatedData:
    A: np.ndarray
    annotation_names: List[str]
    annotation_types: List[str]
    maf: np.ndarray
    alpha_true: np.ndarray
    q_true: np.ndarray
    pi_true: np.ndarray
    component: np.ndarray
    beta: np.ndarray
    bhat: np.ndarray
    se2: np.ndarray
    vb_fit: float
    gamma: np.ndarray
    annotation_scenario: str
    genetic_scenario: str
    replicate: int
    seed: int


@dataclass
class ChainResult:
    theta: np.ndarray
    alpha: np.ndarray
    loglik: np.ndarray
    delta: Optional[np.ndarray] = None


def clip_probability(x: np.ndarray) -> np.ndarray:
    return np.clip(x, EPS, 1.0 - EPS)


def standardize(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    sd = x.std(axis=0, ddof=1)
    if np.any(~np.isfinite(sd)) or np.any(sd <= 1e-10):
        raise ValueError("Cannot standardize a constant annotation")
    return (x - x.mean(axis=0)) / sd


def make_annotation_panel(
    m: int,
    p_annotations: int,
    rng: np.random.Generator,
    binary_fraction: float,
    proxy_correlation: float,
) -> Tuple[np.ndarray, List[str], List[str], np.ndarray]:
    if p_annotations < 6:
        raise ValueError("p_annotations must be at least 6")
    p_binary = max(2, int(round(p_annotations * binary_fraction)))
    p_cont = p_annotations - p_binary
    maf = rng.uniform(0.03, 0.50, size=m)
    rare_score = standardize(
        (-np.log(2.0 * maf * (1.0 - maf)))[:, None]
    )[:, 0]

    continuous = rng.normal(size=(m, p_cont))
    continuous[:, 0] = rare_score + rng.normal(scale=0.25, size=m)
    continuous[:, 1] = (
        proxy_correlation * continuous[:, 0]
        + math.sqrt(1.0 - proxy_correlation**2) * rng.normal(size=m)
    )
    continuous = standardize(continuous)

    binary = np.zeros((m, p_binary), dtype=float)
    prevalence = np.linspace(0.05, 0.45, p_binary)
    for j, prev in enumerate(prevalence):
        latent = rng.normal(size=m)
        if j == 0:
            latent = rare_score + rng.normal(scale=0.40, size=m)
        binary[:, j] = (latent > np.quantile(latent, 1.0 - prev)).astype(float)
    binary = standardize(binary)

    X = np.column_stack((continuous, binary))
    A = np.column_stack((np.ones(m), X))
    names = ["intercept"] + [f"continuous{j + 1}" for j in range(p_cont)] + [
        f"binary{j + 1}" for j in range(p_binary)
    ]
    types = ["intercept"] + ["continuous"] * p_cont + ["binary"] * p_binary
    return A, names, types, maf


def make_alpha_slopes(
    annotation_names: Sequence[str], scenario: str, strength: float
) -> np.ndarray:
    p = len(annotation_names) - 1
    slopes = np.zeros((p, 3), dtype=float)
    cont = [i - 1 for i, x in enumerate(annotation_names) if x.startswith("continuous")]
    binary = [i - 1 for i, x in enumerate(annotation_names) if x.startswith("binary")]

    if scenario == "null":
        pass
    elif scenario == "sparse_independent":
        slopes[cont[0], :] = (0.90, 0.55, 0.35)
        slopes[binary[min(1, len(binary) - 1)], :] = (0.70, 0.35, 0.20)
    elif scenario == "correlated_proxy":
        slopes[cont[0], :] = (0.90, 0.55, 0.35)
    elif scenario == "rare_binary":
        slopes[binary[0], :] = (1.30, 0.75, 0.45)
    elif scenario == "later_stick_only":
        slopes[cont[0], 1:] = (0.80, 0.55)
        slopes[binary[0], 1:] = (0.60, 0.40)
    elif scenario == "dense":
        active = cont[: min(4, len(cont))] + binary[: min(4, len(binary))]
        base = np.linspace(0.45, 0.15, len(active))
        slopes[active, 0] = base
        slopes[active, 1] = 0.60 * base
        slopes[active, 2] = 0.40 * base
    else:
        raise ValueError(f"Unknown annotation scenario: {scenario}")
    return slopes * strength


def calibrate_alpha(
    A: np.ndarray, slopes: np.ndarray, target: Sequence[float]
) -> Tuple[np.ndarray, np.ndarray]:
    alpha = np.vstack((np.zeros(3), slopes))
    q = np.zeros((A.shape[0], 3), dtype=float)
    eligibility = np.ones(A.shape[0], dtype=float)
    for stick in range(3):
        offset = A[:, 1:] @ slopes[:, stick]

        def objective(intercept: float) -> float:
            return float(
                np.sum(eligibility * ndtr(intercept + offset))
                / np.sum(eligibility)
                - target[stick]
            )

        alpha[0, stick] = brentq(objective, -20.0, 20.0, xtol=1e-12)
        q[:, stick] = ndtr(A @ alpha[:, stick])
        eligibility *= q[:, stick]
    return alpha, q


def stick_to_pi(q: np.ndarray) -> np.ndarray:
    q = np.asarray(q, dtype=float)
    if q.ndim != 2 or q.shape[0] < 1 or q.shape[1] < 1:
        raise ValueError("q must be a marker-by-stick matrix")
    if np.any(~np.isfinite(q)) or np.any((q < 0.0) | (q > 1.0)):
        raise ValueError("q must contain probabilities in [0, 1]")
    with np.errstate(divide="ignore"):
        log_q = np.log(q)
        log_1mq = np.log1p(-q)
    log_pi = log_stick_probabilities(log_q, log_1mq)
    pi = np.exp(log_pi)
    return pi / pi.sum(axis=1, keepdims=True)


def categorical_rows(probability: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    u = rng.random(probability.shape[0])
    return (u[:, None] > np.cumsum(probability, axis=1)).sum(axis=1)


def simulate_dataset(
    cfg: Config,
    annotation_scenario: str,
    genetic_scenario: str,
    replicate: int,
    seed: int,
) -> SimulatedData:
    rng = np.random.default_rng(seed)
    A, names, types, maf = make_annotation_panel(
        cfg.m,
        cfg.p_annotations,
        rng,
        cfg.binary_fraction,
        cfg.proxy_correlation,
    )
    slopes = make_alpha_slopes(names, annotation_scenario, cfg.annotation_strength)
    alpha_true, q_true = calibrate_alpha(A, slopes, cfg.target_continuation)
    pi_true = stick_to_pi(q_true)
    component = categorical_rows(pi_true, rng)
    if not np.any(component > 0):
        component[np.argmax(1.0 - pi_true[:, 0])] = 1

    expected_energy = np.sum(pi_true @ GAMMA)
    vb = cfg.h2 / max(expected_energy, EPS)
    beta = np.zeros(cfg.m, dtype=float)
    active = component > 0
    beta[active] = rng.normal(
        scale=np.sqrt(vb * GAMMA[component[active]]), size=active.sum()
    )

    if genetic_scenario == "bayesr_exact":
        pass
    elif genetic_scenario == "major_polygenic":
        idx = np.flatnonzero(active)
        if idx.size >= 2:
            major = idx[np.argmax(np.abs(beta[idx]))]
            remainder = idx[idx != major]
            signs = 1.0 if beta[major] >= 0 else -1.0
            beta[major] = signs * math.sqrt(cfg.h2 * cfg.major_share)
            norm = np.linalg.norm(beta[remainder])
            if norm <= EPS:
                beta[remainder] = rng.normal(size=remainder.size)
                norm = np.linalg.norm(beta[remainder])
            beta[remainder] *= math.sqrt(cfg.h2 * (1.0 - cfg.major_share)) / norm
    elif genetic_scenario == "maf_dependent":
        beta *= (2.0 * maf * (1.0 - maf)) ** (cfg.maf_s / 2.0)
        norm = np.linalg.norm(beta)
        if norm > EPS:
            beta *= math.sqrt(cfg.h2) / norm
    else:
        raise ValueError(f"Unknown genetic scenario: {genetic_scenario}")

    se2 = np.full(cfg.m, (1.0 - cfg.h2) / cfg.n_eff, dtype=float)
    bhat = beta + rng.normal(scale=np.sqrt(se2), size=cfg.m)
    return SimulatedData(
        A=A,
        annotation_names=names,
        annotation_types=types,
        maf=maf,
        alpha_true=alpha_true,
        q_true=q_true,
        pi_true=pi_true,
        component=component,
        beta=beta,
        bhat=bhat,
        se2=se2,
        vb_fit=vb,
        gamma=GAMMA.copy(),
        annotation_scenario=annotation_scenario,
        genetic_scenario=genetic_scenario,
        replicate=replicate,
        seed=seed,
    )


def marker_log_densities(data: SimulatedData) -> np.ndarray:
    variance = data.se2[:, None] + data.vb_fit * data.gamma[None, :]
    return -0.5 * (
        np.log(2.0 * np.pi * variance) + data.bhat[:, None] ** 2 / variance
    )


def build_designs(A: np.ndarray, centers: np.ndarray) -> List[np.ndarray]:
    X = A[:, 1:]
    return [
        np.column_stack((np.ones(A.shape[0]), X - centers[stick]))
        for stick in range(centers.shape[0])
    ]


def theta_to_alpha(theta: np.ndarray, centers: np.ndarray) -> np.ndarray:
    theta = np.asarray(theta, dtype=float)
    alpha = theta.copy()
    for stick in range(theta.shape[0]):
        alpha[stick, 0] = theta[stick, 0] - centers[stick] @ theta[stick, 1:]
    return alpha.T


def alpha_to_theta(alpha: np.ndarray, centers: np.ndarray) -> np.ndarray:
    alpha = np.asarray(alpha, dtype=float)
    theta = alpha.T.copy()
    for stick in range(theta.shape[0]):
        theta[stick, 0] = alpha[0, stick] + centers[stick] @ alpha[1:, stick]
    return theta


def log_stick_probabilities(
    log_q: np.ndarray, log_1mq: np.ndarray
) -> np.ndarray:
    log_q = np.asarray(log_q, dtype=float)
    log_1mq = np.asarray(log_1mq, dtype=float)
    if log_q.ndim != 2 or log_q.shape != log_1mq.shape:
        raise ValueError("log continuation arrays must be marker-by-stick matrices")
    n_marker, n_stick = log_q.shape
    log_pi = np.empty((n_marker, n_stick + 1), dtype=float)
    log_remaining = np.zeros(n_marker, dtype=float)
    for stick in range(n_stick):
        log_pi[:, stick] = log_remaining + log_1mq[:, stick]
        log_remaining = log_remaining + log_q[:, stick]
    log_pi[:, n_stick] = log_remaining
    return log_pi


def log_pi_from_theta(
    theta: np.ndarray, designs: Sequence[np.ndarray]
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    theta = np.asarray(theta, dtype=float)
    if theta.ndim != 2 or len(designs) != theta.shape[0]:
        raise ValueError("theta and centered designs must share the stick axis")
    eta = np.empty((designs[0].shape[0], theta.shape[0]), dtype=float)
    for stick in range(theta.shape[0]):
        eta[:, stick] = designs[stick] @ theta[stick]
    if np.any(~np.isfinite(eta)):
        raise ValueError("continuation linear predictors must be finite")
    q = ndtr(eta)
    logq = log_ndtr(eta)
    log1mq = log_ndtr(-eta)
    log_pi = log_stick_probabilities(logq, log1mq)
    return log_pi, q, eta


def responsibilities(
    theta: np.ndarray, designs: Sequence[np.ndarray], log_density: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    log_pi, q, eta = log_pi_from_theta(theta, designs)
    log_joint = log_pi + log_density
    normalizer = logsumexp(log_joint, axis=1)
    r = np.exp(log_joint - normalizer[:, None])
    return r, q, eta, float(normalizer.sum())


def estimate_fixed_centers(A: np.ndarray, log_density: np.ndarray) -> np.ndarray:
    """Estimate one fixed eligible-set center per continuation stick.

    Only a three-intercept null gating model is used.  The centers are computed
    once and held fixed in every chain, so this is preconditioning rather than
    an iteration-dependent reparameterization.
    """

    one_design = [np.ones((A.shape[0], 1)) for _ in range(3)]
    baseline = np.array([ndtri(0.12), ndtri(0.50), ndtri(0.45)])

    def objective(x: np.ndarray) -> Tuple[float, np.ndarray]:
        theta = x[:, None]
        r, _, eta, ll = responsibilities(theta, one_design, log_density)
        gradient = np.zeros(3)
        log_phi = -0.5 * eta**2 - 0.5 * math.log(2.0 * math.pi)
        log_q = log_ndtr(eta)
        log_1mq = log_ndtr(-eta)
        for s in range(3):
            success = r[:, (s + 1) :].sum(axis=1)
            failure = r[:, s]
            gradient[s] = np.sum(
                success * np.exp(log_phi[:, s] - log_q[:, s])
                - failure * np.exp(log_phi[:, s] - log_1mq[:, s])
            )
        prior_sd = 2.0
        ll_prior = -0.5 * np.sum(((x - baseline) / prior_sd) ** 2)
        gradient -= (x - baseline) / prior_sd**2
        return -(ll + ll_prior), -gradient

    fit = minimize(
        lambda x: objective(x)[0],
        baseline,
        jac=lambda x: objective(x)[1],
        method="L-BFGS-B",
    )
    theta = fit.x[:, None]
    r, _, _, _ = responsibilities(theta, one_design, log_density)
    trials = np.column_stack(
        (
            np.ones(A.shape[0]),
            r[:, 1:].sum(axis=1),
            r[:, 2:].sum(axis=1),
        )
    )
    X = A[:, 1:]
    centers = np.vstack(
        [np.average(X, axis=0, weights=np.maximum(trials[:, s], EPS)) for s in range(3)]
    )
    return centers


class CollapsedModel:
    def __init__(
        self,
        designs: Sequence[np.ndarray],
        log_density: np.ndarray,
        prior_mean: np.ndarray,
        prior_sd: np.ndarray,
    ) -> None:
        self.designs = list(designs)
        self.log_density = np.asarray(log_density)
        self.prior_mean = np.asarray(prior_mean)
        self.prior_sd = np.asarray(prior_sd)
        self.shape = self.prior_mean.shape

    def loglik(self, theta: np.ndarray) -> float:
        return responsibilities(theta, self.designs, self.log_density)[3]

    def logpost_and_grad(self, flat: np.ndarray) -> Tuple[float, np.ndarray]:
        theta = flat.reshape(self.shape)
        r, q, eta, ll = responsibilities(theta, self.designs, self.log_density)
        grad = np.zeros_like(theta)
        log_q = log_ndtr(eta)
        log_1mq = log_ndtr(-eta)
        log_phi = -0.5 * eta**2 - 0.5 * math.log(2.0 * math.pi)
        for s in range(theta.shape[0]):
            success = r[:, (s + 1) :].sum(axis=1)
            failure = r[:, s]
            score = (
                success * np.exp(log_phi[:, s] - log_q[:, s])
                - failure * np.exp(log_phi[:, s] - log_1mq[:, s])
            )
            grad[s] = self.designs[s].T @ score
        standardized = (theta - self.prior_mean) / self.prior_sd
        logprior = -0.5 * np.sum(standardized**2)
        grad -= (theta - self.prior_mean) / self.prior_sd**2
        return ll + logprior, grad.ravel()

    def map(self, start: Optional[np.ndarray] = None) -> np.ndarray:
        if start is None:
            start = self.prior_mean

        def f(x: np.ndarray) -> Tuple[float, np.ndarray]:
            value, gradient = self.logpost_and_grad(x)
            return -value, -gradient

        fit = minimize(
            lambda x: f(x)[0],
            np.asarray(start).ravel(),
            jac=lambda x: f(x)[1],
            method="L-BFGS-B",
            options={"maxiter": 2000, "ftol": 1e-12, "gtol": 1e-7},
        )
        if not fit.success and not np.isfinite(fit.fun):
            raise RuntimeError(f"Collapsed MAP failed: {fit.message}")
        return fit.x.reshape(self.shape)

    def laplace(self) -> Tuple[np.ndarray, np.ndarray]:
        """Posterior mode and a regularized BFGS covariance approximation."""

        def f(x: np.ndarray) -> Tuple[float, np.ndarray]:
            value, gradient = self.logpost_and_grad(x)
            return -value, -gradient

        fit = minimize(
            lambda x: f(x)[0],
            self.prior_mean.ravel(),
            jac=lambda x: f(x)[1],
            method="BFGS",
            options={"maxiter": 2000, "gtol": 1e-6},
        )
        if not np.isfinite(fit.fun):
            raise RuntimeError(f"Collapsed Laplace fit failed: {fit.message}")
        covariance = np.asarray(fit.hess_inv, dtype=float)
        covariance = 0.5 * (covariance + covariance.T)
        eigenvalue, eigenvector = np.linalg.eigh(covariance)
        eigenvalue = np.clip(eigenvalue, 1e-6, 1e3)
        covariance = (eigenvector * eigenvalue) @ eigenvector.T
        return fit.x.reshape(self.shape), covariance


class LowRankCollapsedModel:
    """Collapsed likelihood with a two-factor continuation-slope surface."""

    def __init__(
        self,
        designs: Sequence[np.ndarray],
        log_density: np.ndarray,
        loadings: np.ndarray,
        prior_mean: np.ndarray,
        prior_sd: np.ndarray,
    ) -> None:
        self.designs = list(designs)
        self.log_density = np.asarray(log_density)
        self.loadings = np.asarray(loadings)
        self.prior_mean = np.asarray(prior_mean)
        self.prior_sd = np.asarray(prior_sd)
        self.p = designs[0].shape[1] - 1
        self.n_factor = self.loadings.shape[1]

    def to_theta(self, phi: np.ndarray) -> np.ndarray:
        phi = np.asarray(phi)
        theta = np.zeros((3, self.p + 1), dtype=float)
        theta[:, 0] = phi[:3]
        factor = phi[3:].reshape(self.n_factor, self.p)
        theta[:, 1:] = self.loadings @ factor
        return theta

    def loglik(self, phi: np.ndarray) -> float:
        theta = self.to_theta(phi)
        return responsibilities(theta, self.designs, self.log_density)[3]


def lowrank_loadings() -> np.ndarray:
    """Orthonormal shared-effect and later-stick-effect directions."""

    shared = np.array([1.0, 0.60, 0.40])
    later = np.array([0.0, 1.0, 0.70])
    first = shared / np.linalg.norm(shared)
    later = later - first * (first @ later)
    second = later / np.linalg.norm(later)
    return np.column_stack((first, second))


def elliptical_slice_step(
    theta: np.ndarray,
    prior_mean: np.ndarray,
    prior_sd: np.ndarray,
    loglik,
    rng: np.random.Generator,
    max_steps: int,
) -> Tuple[np.ndarray, float, int]:
    centered = theta - prior_mean
    nu = rng.normal(scale=prior_sd, size=theta.shape)
    current_ll = float(loglik(theta))
    threshold = current_ll + math.log(rng.random())
    angle = rng.uniform(0.0, 2.0 * math.pi)
    angle_min = angle - 2.0 * math.pi
    angle_max = angle
    for step in range(1, max_steps + 1):
        proposal = prior_mean + centered * math.cos(angle) + nu * math.sin(angle)
        proposal_ll = float(loglik(proposal))
        if np.isfinite(proposal_ll) and proposal_ll >= threshold:
            return proposal, proposal_ll, step
        if angle < 0.0:
            angle_min = angle
        else:
            angle_max = angle
        angle = rng.uniform(angle_min, angle_max)
    raise RuntimeError("Elliptical slice sampler exceeded max_slice_steps")


def elliptical_slice_stick_step(
    theta: np.ndarray,
    stick: int,
    prior_mean: np.ndarray,
    prior_sd: np.ndarray,
    loglik,
    rng: np.random.Generator,
    max_steps: int,
) -> Tuple[np.ndarray, float, int]:
    """ESS update of one continuation-stick coefficient block."""

    current = theta[stick] - prior_mean[stick]
    nu = rng.normal(scale=prior_sd[stick], size=current.shape)
    current_ll = float(loglik(theta))
    threshold = current_ll + math.log(rng.random())
    angle = rng.uniform(0.0, 2.0 * math.pi)
    angle_min = angle - 2.0 * math.pi
    angle_max = angle
    for step in range(1, max_steps + 1):
        proposal = theta.copy()
        proposal[stick] = (
            prior_mean[stick]
            + current * math.cos(angle)
            + nu * math.sin(angle)
        )
        proposal_ll = float(loglik(proposal))
        if np.isfinite(proposal_ll) and proposal_ll >= threshold:
            return proposal, proposal_ll, step
        if angle < 0.0:
            angle_min = angle
        else:
            angle_max = angle
        angle = rng.uniform(angle_min, angle_max)
    raise RuntimeError("Blocked elliptical slice sampler exceeded max_slice_steps")


def elliptical_slice_vector_block(
    state: np.ndarray,
    index: np.ndarray,
    prior_mean: np.ndarray,
    prior_sd: np.ndarray,
    loglik,
    rng: np.random.Generator,
    max_steps: int,
) -> Tuple[np.ndarray, float, int]:
    """ESS update for arbitrary coordinates under an independent normal prior."""

    current = state[index] - prior_mean[index]
    nu = rng.normal(scale=prior_sd[index], size=len(index))
    current_ll = float(loglik(state))
    threshold = current_ll + math.log(rng.random())
    angle = rng.uniform(0.0, 2.0 * math.pi)
    angle_min = angle - 2.0 * math.pi
    angle_max = angle
    for step in range(1, max_steps + 1):
        proposal = state.copy()
        proposal[index] = (
            prior_mean[index]
            + current * math.cos(angle)
            + nu * math.sin(angle)
        )
        proposal_ll = float(loglik(proposal))
        if np.isfinite(proposal_ll) and proposal_ll >= threshold:
            return proposal, proposal_ll, step
        if angle < 0.0:
            angle_min = angle
        else:
            angle_max = angle
        angle = rng.uniform(angle_min, angle_max)
    raise RuntimeError("Vector-block elliptical slice sampler exceeded max_slice_steps")


def run_collapsed_chain(
    model: CollapsedModel,
    centers: np.ndarray,
    cfg: Config,
    seed: int,
) -> ChainResult:
    rng = np.random.default_rng(seed)
    # Dispersed prior initialization is intentional: convergence is not
    # manufactured by starting every chain at the common mode.
    theta = rng.normal(model.prior_mean, model.prior_sd)
    keep = []
    ll_keep = []
    for iteration in range(cfg.n_iter):
        # The stick blocks have strong within-block intercept/slope dependence,
        # while cross-stick dependence is weaker.  Updating all three blocks
        # every iteration and adding an occasional global ellipse is a simple
        # interweaving scheme that preserves the same collapsed posterior.
        ll = model.loglik(theta)
        for stick in rng.permutation(3):
            theta, ll, _ = elliptical_slice_stick_step(
                theta,
                int(stick),
                model.prior_mean,
                model.prior_sd,
                model.loglik,
                rng,
                cfg.max_slice_steps,
            )
        if iteration % 5 == 0:
            theta, ll, _ = elliptical_slice_step(
                theta,
                model.prior_mean,
                model.prior_sd,
                model.loglik,
                rng,
                cfg.max_slice_steps,
            )
        if iteration >= cfg.burn and (iteration - cfg.burn) % cfg.thin == 0:
            keep.append(theta.copy())
            ll_keep.append(ll)
    theta_draws = np.asarray(keep)
    alpha_draws = np.asarray([theta_to_alpha(x, centers) for x in theta_draws])
    return ChainResult(theta=theta_draws, alpha=alpha_draws, loglik=np.asarray(ll_keep))


def run_collapsed_hmc_chain(
    model: CollapsedModel,
    centers: np.ndarray,
    mode: np.ndarray,
    covariance: np.ndarray,
    cfg: Config,
    seed: int,
) -> ChainResult:
    """Exact HMC on a Laplace-whitened collapsed alpha posterior."""

    rng = np.random.default_rng(seed)
    dimension = mode.size
    jitter = 1e-10
    for _ in range(8):
        try:
            root = np.linalg.cholesky(covariance + np.eye(dimension) * jitter)
            break
        except np.linalg.LinAlgError:
            jitter *= 10.0
    else:
        raise RuntimeError("Could not factorize Laplace covariance")

    mode_flat = mode.ravel()

    def target(z: np.ndarray) -> Tuple[float, np.ndarray, np.ndarray]:
        theta_flat = mode_flat + root @ z
        value, gradient_theta = model.logpost_and_grad(theta_flat)
        gradient_z = root.T @ gradient_theta
        return value, gradient_z, theta_flat.reshape(model.shape)

    z = rng.normal(scale=1.25, size=dimension)
    value, gradient, theta = target(z)
    step_size = cfg.hmc_initial_step_size
    keep = []
    ll_keep = []
    for iteration in range(cfg.n_iter):
        momentum0 = rng.normal(size=dimension)
        z_new = z.copy()
        momentum = momentum0.copy()
        gradient_new = gradient.copy()
        momentum += 0.5 * step_size * gradient_new
        n_step = max(
            2,
            cfg.hmc_leapfrog_steps + int(rng.integers(-2, 3)),
        )
        proposed_value = -np.inf
        proposed_theta = theta
        for leapfrog in range(n_step):
            z_new += step_size * momentum
            proposed_value, gradient_new, proposed_theta = target(z_new)
            if not np.isfinite(proposed_value) or np.any(~np.isfinite(gradient_new)):
                proposed_value = -np.inf
                break
            if leapfrog != n_step - 1:
                momentum += step_size * gradient_new
        if np.isfinite(proposed_value):
            momentum += 0.5 * step_size * gradient_new
            momentum = -momentum
            log_accept = (
                proposed_value
                - 0.5 * np.dot(momentum, momentum)
                - value
                + 0.5 * np.dot(momentum0, momentum0)
            )
        else:
            log_accept = -np.inf
        accept_probability = min(1.0, math.exp(min(0.0, log_accept)))
        if math.log(rng.random()) < log_accept:
            z = z_new
            value = proposed_value
            gradient = gradient_new
            theta = proposed_theta

        if iteration < cfg.burn:
            learning_rate = 0.05 / math.sqrt(iteration + 1.0)
            log_step = math.log(step_size) + learning_rate * (
                accept_probability - cfg.hmc_target_accept
            )
            step_size = float(np.clip(math.exp(log_step), 0.005, 0.50))
        if iteration >= cfg.burn and (iteration - cfg.burn) % cfg.thin == 0:
            keep.append(theta.copy())
            ll_keep.append(model.loglik(theta))

    theta_draws = np.asarray(keep)
    alpha_draws = np.asarray([theta_to_alpha(x, centers) for x in theta_draws])
    return ChainResult(theta=theta_draws, alpha=alpha_draws, loglik=np.asarray(ll_keep))


def run_lowrank_collapsed_chain(
    model: LowRankCollapsedModel,
    centers: np.ndarray,
    cfg: Config,
    seed: int,
) -> ChainResult:
    rng = np.random.default_rng(seed)
    phi = rng.normal(model.prior_mean, model.prior_sd)
    keep_theta = []
    ll_keep = []
    blocks = [np.arange(3)] + [
        np.arange(3 + factor * model.p, 3 + (factor + 1) * model.p)
        for factor in range(model.n_factor)
    ]
    all_index = np.arange(phi.size)
    for iteration in range(cfg.n_iter):
        ll = model.loglik(phi)
        for block in rng.permutation(len(blocks)):
            phi, ll, _ = elliptical_slice_vector_block(
                phi,
                blocks[int(block)],
                model.prior_mean,
                model.prior_sd,
                model.loglik,
                rng,
                cfg.max_slice_steps,
            )
        if iteration % 5 == 0:
            phi, ll, _ = elliptical_slice_vector_block(
                phi,
                all_index,
                model.prior_mean,
                model.prior_sd,
                model.loglik,
                rng,
                cfg.max_slice_steps,
            )
        if iteration >= cfg.burn and (iteration - cfg.burn) % cfg.thin == 0:
            keep_theta.append(model.to_theta(phi))
            ll_keep.append(ll)
    theta_draws = np.asarray(keep_theta)
    alpha_draws = np.asarray([theta_to_alpha(x, centers) for x in theta_draws])
    return ChainResult(theta=theta_draws, alpha=alpha_draws, loglik=np.asarray(ll_keep))


def draw_truncated_normal(
    mean: np.ndarray, success: np.ndarray, rng: np.random.Generator
) -> np.ndarray:
    boundary = clip_probability(ndtr(-mean))
    u = rng.random(mean.size)
    probability = np.where(success, boundary + u * (1.0 - boundary), u * boundary)
    return mean + ndtri(clip_probability(probability))


def run_hard_gibbs_chain(
    designs: Sequence[np.ndarray],
    centers: np.ndarray,
    log_density: np.ndarray,
    prior_mean: np.ndarray,
    prior_sd: np.ndarray,
    cfg: Config,
    seed: int,
) -> ChainResult:
    rng = np.random.default_rng(seed)
    theta = rng.normal(prior_mean, prior_sd)
    keep = []
    ll_keep = []
    for iteration in range(cfg.n_iter):
        r, _, _, ll = responsibilities(theta, designs, log_density)
        component = categorical_rows(r, rng)
        for stick in range(3):
            eligible = component >= stick
            success = component[eligible] >= (stick + 1)
            D = designs[stick][eligible]
            eta = D @ theta[stick]
            latent = draw_truncated_normal(eta, success, rng)
            precision = D.T @ D + np.diag(1.0 / prior_sd[stick] ** 2)
            covariance = np.linalg.inv(precision)
            rhs = D.T @ latent + prior_mean[stick] / prior_sd[stick] ** 2
            mean = covariance @ rhs
            theta[stick] = rng.multivariate_normal(mean, covariance)
        if iteration >= cfg.burn and (iteration - cfg.burn) % cfg.thin == 0:
            keep.append(theta.copy())
            ll_keep.append(ll)
    theta_draws = np.asarray(keep)
    alpha_draws = np.asarray([theta_to_alpha(x, centers) for x in theta_draws])
    return ChainResult(theta=theta_draws, alpha=alpha_draws, loglik=np.asarray(ll_keep))


def run_collapsed_selection_chain(
    designs: Sequence[np.ndarray],
    centers: np.ndarray,
    log_density: np.ndarray,
    prior_mean: np.ndarray,
    cfg: Config,
    seed: int,
) -> ChainResult:
    """Collapsed SBayesRC-S-style spike/slab sampler.

    Each annotation-by-stick coefficient has its own inclusion indicator.  A
    toggle proposes the coefficient jointly with its destination state.  Prior
    and proposal densities cancel, leaving a simple exact MH likelihood ratio
    plus inclusion-odds ratio.  ESS then updates all coefficients conditional
    on the current inclusion states.  This is particularly important for the
    later-stick-only scenario, which a shared three-stick indicator can miss.
    """

    rng = np.random.default_rng(seed)
    p = prior_mean.shape[1] - 1
    delta = rng.binomial(1, cfg.inclusion_prior, size=(3, p)).astype(int)

    def make_sd(state: np.ndarray) -> np.ndarray:
        sd = np.full_like(prior_mean, cfg.intercept_prior_sd, dtype=float)
        slope_sd = np.where(state > 0, cfg.slab_sd, cfg.spike_sd)
        sd[:, 1:] = slope_sd
        return sd

    prior_sd = make_sd(delta)
    theta = rng.normal(prior_mean, prior_sd)
    current_ll = responsibilities(theta, designs, log_density)[3]
    annotation_correlation = np.corrcoef(designs[0][:, 1:], rowvar=False)
    correlated_pairs = [
        (left, right)
        for left in range(p)
        for right in range(left + 1, p)
        if abs(annotation_correlation[left, right]) >= cfg.correlation_swap_threshold
    ]
    keep = []
    delta_keep = []
    ll_keep = []
    log_prior_odds = math.log(cfg.inclusion_prior) - math.log1p(-cfg.inclusion_prior)

    for iteration in range(cfg.n_iter):
        if iteration % cfg.selection_toggle_every == 0:
            coordinate = [(stick, annotation) for stick in range(3) for annotation in range(p)]
            for position in rng.permutation(len(coordinate)):
                stick, annotation = coordinate[int(position)]
                proposed_delta = 1 - delta[stick, annotation]
                proposed_theta = theta.copy()
                proposed_sd = cfg.slab_sd if proposed_delta else cfg.spike_sd
                proposed_theta[stick, annotation + 1] = rng.normal(scale=proposed_sd)
                proposed_ll = responsibilities(
                    proposed_theta, designs, log_density
                )[3]
                log_ratio = proposed_ll - current_ll
                log_ratio += log_prior_odds if proposed_delta else -log_prior_odds
                if math.log(rng.random()) < log_ratio:
                    theta = proposed_theta
                    delta[stick, annotation] = proposed_delta
                    current_ll = proposed_ll
            # Exact symmetric mode-swap moves help chains traverse alternative
            # allocations of one signal to nearly exchangeable annotations.
            # They do not make the individual coefficients identifiable; they
            # make the posterior ambiguity visible as split PIPs.
            for stick in rng.permutation(3):
                for position in rng.permutation(len(correlated_pairs)):
                    left, right = correlated_pairs[int(position)]
                    proposal = theta.copy()
                    proposal[stick, [left + 1, right + 1]] = proposal[
                        stick, [right + 1, left + 1]
                    ]
                    proposed_ll = responsibilities(
                        proposal, designs, log_density
                    )[3]
                    if math.log(rng.random()) < proposed_ll - current_ll:
                        theta = proposal
                        delta[stick, [left, right]] = delta[stick, [right, left]]
                        current_ll = proposed_ll
            prior_sd = make_sd(delta)

        loglik = lambda x: responsibilities(x, designs, log_density)[3]
        for stick in rng.permutation(3):
            theta, current_ll, _ = elliptical_slice_stick_step(
                theta,
                int(stick),
                prior_mean,
                prior_sd,
                loglik,
                rng,
                cfg.max_slice_steps,
            )
        if iteration % 5 == 0:
            theta, current_ll, _ = elliptical_slice_step(
                theta,
                prior_mean,
                prior_sd,
                loglik,
                rng,
                cfg.max_slice_steps,
            )
        if iteration >= cfg.burn and (iteration - cfg.burn) % cfg.thin == 0:
            keep.append(theta.copy())
            delta_keep.append(delta.copy())
            ll_keep.append(current_ll)

    theta_draws = np.asarray(keep)
    alpha_draws = np.asarray([theta_to_alpha(x, centers) for x in theta_draws])
    return ChainResult(
        theta=theta_draws,
        alpha=alpha_draws,
        loglik=np.asarray(ll_keep),
        delta=np.asarray(delta_keep),
    )


def split_rhat(chains: np.ndarray) -> float:
    """Classical split R-hat for shape (chain, draw)."""

    x = np.asarray(chains, dtype=float)
    if x.ndim != 2 or x.shape[0] < 2 or x.shape[1] < 8:
        return np.nan
    half = x.shape[1] // 2
    x = np.concatenate((x[:, :half], x[:, -half:]), axis=0)
    n = x.shape[1]
    within = np.mean(np.var(x, axis=1, ddof=1))
    if not np.isfinite(within) or within <= 0:
        return np.nan
    between = n * np.var(np.mean(x, axis=1), ddof=1)
    return float(np.sqrt(((n - 1.0) / n * within + between / n) / within))


def one_chain_ess(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    n = x.size
    if n < 8 or np.std(x, ddof=1) <= 0:
        return np.nan
    centered = x - x.mean()
    fft_n = 1 << (2 * n - 1).bit_length()
    fft = np.fft.rfft(centered, n=fft_n)
    acov = np.fft.irfft(fft * np.conjugate(fft), n=fft_n)[:n]
    acov /= np.arange(n, 0, -1)
    rho = acov / acov[0]
    pair = rho[1::2][: (n - 1) // 2] + rho[2::2][: (n - 1) // 2]
    positive = pair[pair > 0]
    if positive.size < pair.size:
        positive = pair[: np.argmax(pair <= 0)]
    tau = 1.0 + 2.0 * positive.sum()
    return float(min(n / max(tau, 1.0), n))


def combined_ess(chains: np.ndarray) -> float:
    values = [one_chain_ess(x) for x in np.asarray(chains)]
    return float(np.nansum(values)) if np.any(np.isfinite(values)) else np.nan


def annotation_correlation_groups(
    A: np.ndarray, names: Sequence[str], threshold: float
) -> Tuple[np.ndarray, np.ndarray]:
    """Connected components of the absolute-correlation threshold graph."""

    X = np.asarray(A[:, 1:], dtype=float)
    p = X.shape[1]
    parent = np.arange(p)

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = int(parent[x])
        return x

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    correlation = np.corrcoef(X, rowvar=False)
    for left in range(p):
        for right in range(left + 1, p):
            if abs(correlation[left, right]) >= threshold:
                union(left, right)
    roots = [find(j) for j in range(p)]
    unique = {root: i + 1 for i, root in enumerate(dict.fromkeys(roots))}
    groups = np.array([unique[root] for root in roots], dtype=int)
    return groups, correlation


def posterior_surfaces(
    chains: Sequence[ChainResult], designs: Sequence[np.ndarray]
) -> Tuple[np.ndarray, np.ndarray]:
    q_sum = np.zeros((designs[0].shape[0], 3), dtype=float)
    n = 0
    for chain in chains:
        for theta in chain.theta:
            _, q, _ = log_pi_from_theta(theta, designs)
            q_sum += q
            n += 1
    q_mean = q_sum / n
    return q_mean, stick_to_pi(q_mean)


def summarize_functionals(
    method: str,
    chains: Sequence[ChainResult],
    data: SimulatedData,
    designs: Sequence[np.ndarray],
) -> pd.DataFrame:
    """Diagnose mixing for probability-scale quantities, not only raw alpha."""

    marker_order = np.argsort(data.A[:, 1])
    landmarks = marker_order[
        np.linspace(0, len(marker_order) - 1, min(12, len(marker_order))).astype(int)
    ]
    values_by_chain: List[Dict[str, np.ndarray]] = []
    for chain in chains:
        q_draws = []
        pi_draws = []
        for theta in chain.theta:
            _, q, _ = log_pi_from_theta(theta, designs)
            q_draws.append(q)
            pi_draws.append(stick_to_pi(q))
        q_array = np.asarray(q_draws)
        pi_array = np.asarray(pi_draws)
        values: Dict[str, np.ndarray] = {"collapsed_loglik": chain.loglik}
        for stick in range(3):
            values[f"q_mean_stick{stick + 1}"] = q_array[:, :, stick].mean(axis=1)
            values[f"q_sd_stick{stick + 1}"] = q_array[:, :, stick].std(axis=1)
            for landmark_number, marker in enumerate(landmarks, start=1):
                values[f"q_stick{stick + 1}_landmark{landmark_number:02d}"] = (
                    q_array[:, marker, stick]
                )
        for component in range(4):
            values[f"pi_mean_component{component}"] = pi_array[:, :, component].mean(axis=1)
        values_by_chain.append(values)

    rows = []
    for functional in values_by_chain[0]:
        samples = np.stack([x[functional] for x in values_by_chain], axis=0)
        rows.append(
            {
                "replicate": data.replicate,
                "annotation_scenario": data.annotation_scenario,
                "genetic_scenario": data.genetic_scenario,
                "method": method,
                "functional": functional,
                "posterior_mean": samples.mean(),
                "posterior_sd": samples.reshape(-1).std(ddof=1),
                "rhat": split_rhat(samples),
                "ess": combined_ess(samples),
            }
        )
    return pd.DataFrame(rows)


def summarize_method(
    method: str,
    chains: Sequence[ChainResult],
    data: SimulatedData,
    designs: Sequence[np.ndarray],
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    alpha_stack = np.stack([x.alpha for x in chains], axis=0)
    # chain x draw x coefficient x stick
    posterior_mean = alpha_stack.mean(axis=(0, 1))
    posterior_sd = alpha_stack.reshape(-1, *alpha_stack.shape[2:]).std(axis=0, ddof=1)
    selection_pip = None
    if chains[0].delta is not None:
        # draw x stick x annotation
        selection_pip = np.concatenate([x.delta for x in chains], axis=0).mean(axis=0)

    coefficient_rows = []
    for coefficient, name in enumerate(data.annotation_names):
        for stick in range(3):
            samples = alpha_stack[:, :, coefficient, stick]
            truth = data.alpha_true[coefficient, stick]
            pip = 1.0 if coefficient == 0 else (
                selection_pip[stick, coefficient - 1]
                if selection_pip is not None
                else np.nan
            )
            coefficient_rows.append(
                {
                    "replicate": data.replicate,
                    "annotation_scenario": data.annotation_scenario,
                    "genetic_scenario": data.genetic_scenario,
                    "method": method,
                    "coefficient": name,
                    "annotation_type": data.annotation_types[coefficient],
                    "stick": stick + 1,
                    "truth": truth,
                    "posterior_mean": posterior_mean[coefficient, stick],
                    "posterior_sd": posterior_sd[coefficient, stick],
                    "bias": posterior_mean[coefficient, stick] - truth,
                    "active_truth": coefficient > 0 and abs(truth) > 1e-12,
                    "rhat": split_rhat(samples),
                    "ess": combined_ess(samples),
                    "annotation_pip": pip,
                }
            )

    q_mean, pi_mean = posterior_surfaces(chains, designs)
    functionals = summarize_functionals(method, chains, data, designs)
    surface = pd.DataFrame(
        [
            {
                "replicate": data.replicate,
                "annotation_scenario": data.annotation_scenario,
                "genetic_scenario": data.genetic_scenario,
                "method": method,
                "q_rmse": np.sqrt(np.mean((q_mean - data.q_true) ** 2)),
                "pi_rmse": np.sqrt(np.mean((pi_mean - data.pi_true) ** 2)),
                "q_correlation": np.corrcoef(q_mean.ravel(), data.q_true.ravel())[0, 1],
                "pi_correlation": np.corrcoef(pi_mean.ravel(), data.pi_true.ravel())[0, 1],
                "max_rhat": np.nanmax([x["rhat"] for x in coefficient_rows]),
                "min_ess": np.nanmin([x["ess"] for x in coefficient_rows]),
                "max_functional_rhat": functionals["rhat"].max(),
                "min_functional_ess": functionals["ess"].min(),
            }
        ]
    )
    selection_rows = []
    if selection_pip is not None:
        delta_draws = np.concatenate([x.delta for x in chains], axis=0)
        groups, _ = annotation_correlation_groups(
            data.A, data.annotation_names[1:], threshold=0.85
        )
        truth_active = np.abs(data.alpha_true[1:, :]) > 1e-12
        for stick in range(3):
            for j, name in enumerate(data.annotation_names[1:]):
                members = np.flatnonzero(groups == groups[j])
                group_pip = np.mean(np.any(delta_draws[:, stick, members] > 0, axis=1))
                selection_rows.append(
                    {
                        "replicate": data.replicate,
                        "annotation_scenario": data.annotation_scenario,
                        "genetic_scenario": data.genetic_scenario,
                        "method": method,
                        "coefficient": name,
                        "stick": stick + 1,
                        "active_truth": truth_active[j, stick],
                        "annotation_pip": selection_pip[stick, j],
                        "correlation_group": f"group{groups[j]:03d}",
                        "group_active_truth": bool(np.any(truth_active[members, stick])),
                        "group_pip": group_pip,
                    }
                )
    return (
        pd.DataFrame(coefficient_rows),
        surface,
        pd.DataFrame(selection_rows),
        functionals,
    )


def fit_one_dataset(
    data: SimulatedData, cfg: Config, seed: int
) -> Tuple[
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
]:
    log_density = marker_log_densities(data)
    centers = estimate_fixed_centers(data.A, log_density)
    designs = build_designs(data.A, centers)
    d = data.A.shape[1]
    prior_mean = np.zeros((3, d), dtype=float)
    prior_mean[:, 0] = ndtri(np.asarray(cfg.target_continuation))
    prior_sd = np.full((3, d), cfg.slope_prior_sd, dtype=float)
    prior_sd[:, 0] = cfg.intercept_prior_sd
    model = CollapsedModel(designs, log_density, prior_mean, prior_sd)
    hmc_mode = hmc_covariance = None
    if cfg.run_hmc:
        hmc_mode, hmc_covariance = model.laplace()
    loadings = lowrank_loadings()
    lowrank_prior_mean = np.zeros(3 + loadings.shape[1] * (d - 1), dtype=float)
    lowrank_prior_mean[:3] = ndtri(np.asarray(cfg.target_continuation))
    lowrank_prior_sd = np.full_like(
        lowrank_prior_mean, cfg.lowrank_factor_sd, dtype=float
    )
    lowrank_prior_sd[:3] = cfg.intercept_prior_sd
    lowrank_model = LowRankCollapsedModel(
        designs,
        log_density,
        loadings,
        lowrank_prior_mean,
        lowrank_prior_sd,
    )

    methods: Dict[str, List[ChainResult]] = {"hard_gibbs": [], "collapsed_ess": []}
    if cfg.run_hmc:
        methods["collapsed_hmc"] = []
    if cfg.run_lowrank:
        methods["collapsed_lowrank"] = []
    if cfg.run_selection:
        methods["collapsed_selection"] = []

    for chain in range(cfg.n_chains):
        chain_seed = seed + 10007 * (chain + 1)
        methods["hard_gibbs"].append(
            run_hard_gibbs_chain(
                designs,
                centers,
                log_density,
                prior_mean,
                prior_sd,
                cfg,
                chain_seed + 101,
            )
        )
        methods["collapsed_ess"].append(
            run_collapsed_chain(model, centers, cfg, chain_seed + 211)
        )
        if cfg.run_hmc:
            methods["collapsed_hmc"].append(
                run_collapsed_hmc_chain(
                    model,
                    centers,
                    hmc_mode,
                    hmc_covariance,
                    cfg,
                    chain_seed + 239,
                )
            )
        if cfg.run_lowrank:
            methods["collapsed_lowrank"].append(
                run_lowrank_collapsed_chain(
                    lowrank_model, centers, cfg, chain_seed + 263
                )
            )
        if cfg.run_selection:
            methods["collapsed_selection"].append(
                run_collapsed_selection_chain(
                    designs,
                    centers,
                    log_density,
                    prior_mean,
                    cfg,
                    chain_seed + 307,
                )
            )

    coefficient = []
    surface = []
    selection = []
    functional = []
    for method, chains in methods.items():
        c, s, z, f = summarize_method(method, chains, data, designs)
        coefficient.append(c)
        surface.append(s)
        functional.append(f)
        if not z.empty:
            selection.append(z)

    center_rows = []
    for stick in range(3):
        for j, name in enumerate(data.annotation_names[1:]):
            center_rows.append(
                {
                    "replicate": data.replicate,
                    "annotation_scenario": data.annotation_scenario,
                    "genetic_scenario": data.genetic_scenario,
                    "stick": stick + 1,
                    "coefficient": name,
                    "fixed_center": centers[stick, j],
                }
            )
    groups, correlation = annotation_correlation_groups(
        data.A, data.annotation_names[1:], cfg.correlation_swap_threshold
    )
    group_rows = []
    for j, name in enumerate(data.annotation_names[1:]):
        members = np.flatnonzero(groups == groups[j])
        within = [
            abs(correlation[j, k]) for k in members if k != j
        ]
        group_rows.append(
            {
                "replicate": data.replicate,
                "annotation_scenario": data.annotation_scenario,
                "genetic_scenario": data.genetic_scenario,
                "coefficient": name,
                "correlation_group": f"group{groups[j]:03d}",
                "group_size": len(members),
                "maximum_within_group_abs_correlation": max(within) if within else np.nan,
            }
        )
    return (
        pd.concat(coefficient, ignore_index=True),
        pd.concat(surface, ignore_index=True),
        pd.concat(selection, ignore_index=True) if selection else pd.DataFrame(),
        pd.DataFrame(center_rows),
        pd.DataFrame(group_rows),
        pd.concat(functional, ignore_index=True),
    )


def truth_tables(data: SimulatedData) -> Tuple[pd.DataFrame, pd.DataFrame]:
    alpha_rows = []
    for coefficient, name in enumerate(data.annotation_names):
        for stick in range(3):
            alpha_rows.append(
                {
                    "replicate": data.replicate,
                    "annotation_scenario": data.annotation_scenario,
                    "genetic_scenario": data.genetic_scenario,
                    "coefficient": name,
                    "annotation_type": data.annotation_types[coefficient],
                    "stick": stick + 1,
                    "alpha_true": data.alpha_true[coefficient, stick],
                    "active_truth": coefficient > 0
                    and abs(data.alpha_true[coefficient, stick]) > 1e-12,
                }
            )
    global_row = {
        "replicate": data.replicate,
        "annotation_scenario": data.annotation_scenario,
        "genetic_scenario": data.genetic_scenario,
        "n_markers": data.A.shape[0],
        "n_annotations": data.A.shape[1] - 1,
        "realized_active_markers": int(np.sum(data.component > 0)),
        "realized_component0": float(np.mean(data.component == 0)),
        "realized_component1": float(np.mean(data.component == 1)),
        "realized_component2": float(np.mean(data.component == 2)),
        "realized_component3": float(np.mean(data.component == 3)),
        "realized_beta_energy": float(np.sum(data.beta**2)),
        "expected_beta_energy": float(data.vb_fit * np.sum(data.pi_true @ data.gamma)),
        "vb_fit": data.vb_fit,
    }
    return pd.DataFrame(alpha_rows), pd.DataFrame([global_row])


def aggregate_recovery(coefficient: pd.DataFrame, surface: pd.DataFrame) -> pd.DataFrame:
    slopes = coefficient[coefficient["coefficient"] != "intercept"].copy()
    rows = []
    for method, frame in slopes.groupby("method", sort=False):
        active = frame[frame["active_truth"]]
        inactive = frame[~frame["active_truth"]]
        null = frame[frame["annotation_scenario"] == "null"]
        s = surface[surface["method"] == method]
        rows.append(
            {
                "method": method,
                "active_alpha_rmse": np.sqrt(np.mean(active["bias"] ** 2))
                if len(active)
                else np.nan,
                "active_alpha_bias": np.mean(active["bias"]) if len(active) else np.nan,
                "inactive_alpha_rmse": np.sqrt(np.mean(inactive["bias"] ** 2))
                if len(inactive)
                else np.nan,
                "null_alpha_rmse": np.sqrt(np.mean(null["bias"] ** 2))
                if len(null)
                else np.nan,
                "mean_q_rmse": s["q_rmse"].mean(),
                "mean_pi_rmse": s["pi_rmse"].mean(),
                "maximum_rhat": frame["rhat"].max(),
                "minimum_ess": frame["ess"].min(),
                "maximum_functional_rhat": s["max_functional_rhat"].max(),
                "minimum_functional_ess": s["min_functional_ess"].min(),
            }
        )
    return pd.DataFrame(rows)


def selection_summary(selection: pd.DataFrame) -> pd.DataFrame:
    if selection.empty:
        return pd.DataFrame()
    return (
        selection.groupby(["method", "active_truth"], as_index=False)
        .agg(
            mean_pip=("annotation_pip", "mean"),
            median_pip=("annotation_pip", "median"),
            n=("annotation_pip", "size"),
        )
        .sort_values(["method", "active_truth"], ascending=[True, False])
    )


def selection_group_tables(selection: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if selection.empty:
        return pd.DataFrame(), pd.DataFrame()
    columns = [
        "replicate",
        "annotation_scenario",
        "genetic_scenario",
        "method",
        "stick",
        "correlation_group",
        "group_active_truth",
        "group_pip",
    ]
    recovery = selection[columns].drop_duplicates().reset_index(drop=True)
    summary = (
        recovery.groupby(["method", "group_active_truth"], as_index=False)
        .agg(
            mean_group_pip=("group_pip", "mean"),
            median_group_pip=("group_pip", "median"),
            n=("group_pip", "size"),
        )
        .sort_values(["method", "group_active_truth"], ascending=[True, False])
    )
    return recovery, summary


def exactness_checks(
    truth_alpha: pd.DataFrame,
    global_truth: pd.DataFrame,
    surface: pd.DataFrame,
    coefficient: pd.DataFrame,
) -> pd.DataFrame:
    checks = {
        "nonfinite_truth_alpha": int((~np.isfinite(truth_alpha["alpha_true"])).sum()),
        "beta_energy_expectation_error_exact_only": np.nan,
        "invalid_surface_metrics": int(
            (~np.isfinite(surface[["q_rmse", "pi_rmse"]].to_numpy())).sum()
        ),
        "invalid_probability_rmse": int(
            ((surface[["q_rmse", "pi_rmse"]].to_numpy() < 0)
             | (surface[["q_rmse", "pi_rmse"]].to_numpy() > 1)).sum()
        ),
        "invalid_posterior_means": int((~np.isfinite(coefficient["posterior_mean"])).sum()),
    }
    exact = global_truth[global_truth["genetic_scenario"] == "bayesr_exact"]
    if len(exact):
        checks["beta_energy_expectation_error_exact_only"] = float(
            np.mean(exact["realized_beta_energy"] - exact["expected_beta_energy"])
        )
    return pd.DataFrame(
        [{"check": key, "value": value} for key, value in checks.items()]
    )


def make_validation_pdf(
    recovery: pd.DataFrame,
    surface: pd.DataFrame,
    coefficient: pd.DataFrame,
    selection: pd.DataFrame,
    output_file: Path,
) -> None:
    os.environ.setdefault("MPLCONFIGDIR", str(output_file.parent / ".matplotlib"))
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    colors = {
        "hard_gibbs": "#c44e52",
        "collapsed_ess": "#4c72b0",
        "collapsed_hmc": "#64b5cd",
        "collapsed_lowrank": "#8172b2",
        "collapsed_selection": "#55a868",
    }
    with PdfPages(output_file) as pdf:
        fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))
        metrics = [
            ("active_alpha_rmse", "Active alpha RMSE"),
            ("inactive_alpha_rmse", "Inactive alpha RMSE"),
            ("mean_pi_rmse", "Mixture-surface RMSE"),
            ("maximum_rhat", "Maximum split R-hat"),
        ]
        for ax, (metric, title) in zip(axes.ravel(), metrics):
            x = np.arange(len(recovery))
            ax.bar(
                x,
                recovery[metric],
                color=[colors.get(m, "grey") for m in recovery["method"]],
            )
            ax.set_xticks(x, recovery["method"], rotation=20, ha="right")
            ax.set_title(title)
            ax.grid(axis="y", alpha=0.25)
            if metric == "maximum_rhat":
                ax.axhline(1.01, color="black", linestyle="--", linewidth=1)
        fig.suptitle("SBayesRC alpha solution v5: overall recovery")
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        fig, axes = plt.subplots(1, 2, figsize=(11, 5))
        for method, frame in surface.groupby("method", sort=False):
            axes[0].scatter(
                frame["q_rmse"],
                frame["pi_rmse"],
                label=method,
                alpha=0.75,
                color=colors.get(method),
            )
            axes[1].scatter(
                frame["max_rhat"],
                frame["min_ess"],
                label=method,
                alpha=0.75,
                color=colors.get(method),
            )
        axes[0].set(xlabel="q RMSE", ylabel="pi RMSE", title="Probability recovery")
        axes[1].set(xlabel="maximum split R-hat", ylabel="minimum ESS", title="Mixing")
        axes[0].legend()
        axes[1].axvline(1.01, color="black", linestyle="--", linewidth=1)
        for ax in axes:
            ax.grid(alpha=0.25)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        active = coefficient[
            (coefficient["coefficient"] != "intercept") & coefficient["active_truth"]
        ]
        fig, ax = plt.subplots(figsize=(8.5, 7))
        for method, frame in active.groupby("method", sort=False):
            ax.scatter(
                frame["truth"],
                frame["posterior_mean"],
                label=method,
                alpha=0.65,
                color=colors.get(method),
            )
        bounds = [
            min(active["truth"].min(), active["posterior_mean"].min()),
            max(active["truth"].max(), active["posterior_mean"].max()),
        ] if len(active) else [-1, 1]
        ax.plot(bounds, bounds, color="black", linestyle="--")
        ax.set(xlabel="True active alpha", ylabel="Posterior mean", title="Active coefficient recovery")
        ax.legend()
        ax.grid(alpha=0.25)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        if not selection.empty:
            fig, ax = plt.subplots(figsize=(8.5, 5.5))
            inactive = selection.loc[~selection["active_truth"], "annotation_pip"]
            active_pip = selection.loc[selection["active_truth"], "annotation_pip"]
            ax.hist(inactive, bins=np.linspace(0, 1, 21), alpha=0.65, label="inactive")
            ax.hist(active_pip, bins=np.linspace(0, 1, 21), alpha=0.65, label="active")
            ax.set(xlabel="Posterior inclusion probability", ylabel="Count", title="Collapsed structured selection")
            ax.legend()
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)


def write_readme(output_dir: Path, cfg: Config, elapsed: float) -> None:
    text = f"""# SBayesRC alpha solution v5 output

This run used the `{cfg.preset}` preset and completed in {elapsed:.1f} seconds.

The central test is `hard_gibbs` versus `collapsed_hmc`. Both target the same
independent-summary model and use the same continuation-probit surface and
Gaussian priors. The only substantive computational difference is whether the
latent marker allocations/effects are sampled before alpha (`hard_gibbs`) or
integrated out of the alpha likelihood (`collapsed_hmc`). `collapsed_ess`
targets the latter posterior with a less geometry-aware sampler.

`collapsed_selection` adds one spike/slab indicator per annotation-by-stick
coefficient. Its joint coefficient/indicator toggle is a valid
Metropolis-Hastings step; it is not the plug-in EM selection approximation used
in v4.2.

Start with:

* `method_comparison_v5.csv`
* `surface_recovery_v5.csv`
* `coefficient_recovery_v5.csv`
* `selection_summary_v5.csv`
* `selection_group_summary_v5.csv`
* `functional_diagnostics_v5.csv`
* `validation_v5.pdf`

Interpretation guardrails:

* `bayesr_exact` is the decisive computation test because the fitted collapsed
  likelihood is correctly specified.
* `major_polygenic` and `maf_dependent` are robustness tests of misspecification.
* Good recovery in a very small run is not sufficient evidence for a package
  change. Use the `full` preset before drawing final operating-characteristic
  conclusions.
* The independent-summary likelihood is exact for negligible LD. Production
  use with LD requires blockwise likelihoods or a clearly labelled composite
  likelihood; multiplying marginal marker likelihoods is not exact under LD.
"""
    (output_dir / "README_results_v5.md").write_text(text, encoding="utf-8")


def apply_preset(cfg: Config, preset: str) -> Config:
    cfg.preset = preset
    if preset == "smoke":
        cfg.m = 160
        cfg.p_annotations = 8
        cfg.annotation_scenarios = ("null", "sparse_independent")
        cfg.genetic_scenarios = ("bayesr_exact",)
        cfg.n_per_cell = 1
        cfg.n_chains = 2
        cfg.n_iter = 500
        cfg.burn = 200
        cfg.thin = 2
    elif preset == "focused":
        cfg.m = 300
        cfg.p_annotations = 12
        cfg.annotation_scenarios = (
            "null",
            "sparse_independent",
            "correlated_proxy",
            "later_stick_only",
        )
        cfg.genetic_scenarios = ("bayesr_exact",)
        cfg.n_per_cell = 2
        cfg.n_chains = 3
        cfg.n_iter = 1400
        cfg.burn = 600
        cfg.thin = 2
    elif preset == "robustness":
        cfg.m = 300
        cfg.p_annotations = 12
        cfg.annotation_scenarios = ("sparse_independent", "correlated_proxy")
        cfg.genetic_scenarios = ("bayesr_exact", "major_polygenic", "maf_dependent")
        cfg.n_per_cell = 3
        cfg.n_chains = 3
        cfg.n_iter = 1600
        cfg.burn = 700
        cfg.thin = 2
    elif preset == "full":
        cfg.m = 600
        cfg.p_annotations = 30
        cfg.annotation_scenarios = (
            "null",
            "sparse_independent",
            "correlated_proxy",
            "rare_binary",
            "later_stick_only",
            "dense",
        )
        cfg.genetic_scenarios = ("bayesr_exact", "major_polygenic", "maf_dependent")
        cfg.n_per_cell = 20
        cfg.n_chains = 4
        cfg.n_iter = 2400
        cfg.burn = 1000
        cfg.thin = 2
    else:
        raise ValueError("preset must be smoke, focused, robustness, or full")
    return cfg


def run(cfg: Config, output_dir: Path) -> Dict[str, pd.DataFrame]:
    output_dir.mkdir(parents=True, exist_ok=True)
    coefficient_all = []
    surface_all = []
    selection_all = []
    center_all = []
    group_all = []
    functional_all = []
    alpha_truth_all = []
    global_truth_all = []
    manifest = []
    total = (
        len(cfg.annotation_scenarios)
        * len(cfg.genetic_scenarios)
        * cfg.n_per_cell
    )
    index = 0
    start_time = time.time()
    for annotation_scenario in cfg.annotation_scenarios:
        for genetic_scenario in cfg.genetic_scenarios:
            for replicate_in_cell in range(1, cfg.n_per_cell + 1):
                index += 1
                seed = cfg.seed + 1000003 * index
                print(
                    f"dataset {index} / {total}: {annotation_scenario} x "
                    f"{genetic_scenario} (replicate {replicate_in_cell})",
                    flush=True,
                )
                data = simulate_dataset(
                    cfg, annotation_scenario, genetic_scenario, index, seed
                )
                (
                    coefficient,
                    surface,
                    selection,
                    centers,
                    groups,
                    functionals,
                ) = fit_one_dataset(data, cfg, seed + 700001)
                alpha_truth, global_truth = truth_tables(data)
                coefficient_all.append(coefficient)
                surface_all.append(surface)
                if not selection.empty:
                    selection_all.append(selection)
                center_all.append(centers)
                group_all.append(groups)
                functional_all.append(functionals)
                alpha_truth_all.append(alpha_truth)
                global_truth_all.append(global_truth)
                manifest.append(
                    {
                        "replicate": index,
                        "replicate_in_cell": replicate_in_cell,
                        "annotation_scenario": annotation_scenario,
                        "genetic_scenario": genetic_scenario,
                        "seed": seed,
                    }
                )

    coefficient = pd.concat(coefficient_all, ignore_index=True)
    surface = pd.concat(surface_all, ignore_index=True)
    selection = pd.concat(selection_all, ignore_index=True) if selection_all else pd.DataFrame()
    centers = pd.concat(center_all, ignore_index=True)
    groups = pd.concat(group_all, ignore_index=True)
    functionals = pd.concat(functional_all, ignore_index=True)
    alpha_truth = pd.concat(alpha_truth_all, ignore_index=True)
    global_truth = pd.concat(global_truth_all, ignore_index=True)
    recovery = aggregate_recovery(coefficient, surface)
    select_summary = selection_summary(selection)
    select_group_recovery, select_group_summary = selection_group_tables(selection)
    checks = exactness_checks(alpha_truth, global_truth, surface, coefficient)
    manifest_frame = pd.DataFrame(manifest)
    elapsed = time.time() - start_time

    tables = {
        "coefficient_recovery_v5.csv": coefficient,
        "surface_recovery_v5.csv": surface,
        "selection_recovery_v5.csv": selection,
        "fixed_continuation_centers_v5.csv": centers,
        "alpha_truth_v5.csv": alpha_truth,
        "global_truth_v5.csv": global_truth,
        "method_comparison_v5.csv": recovery,
        "selection_summary_v5.csv": select_summary,
        "selection_group_recovery_v5.csv": select_group_recovery,
        "selection_group_summary_v5.csv": select_group_summary,
        "annotation_correlation_groups_v5.csv": groups,
        "functional_diagnostics_v5.csv": functionals,
        "exactness_checks_v5.csv": checks,
        "simulation_manifest_v5.csv": manifest_frame,
    }
    for filename, frame in tables.items():
        frame.to_csv(output_dir / filename, index=False)

    config_dict = asdict(cfg)
    config_dict["python_version"] = sys.version
    config_dict["platform"] = platform.platform()
    config_dict["elapsed_seconds"] = elapsed
    (output_dir / "configuration_v5.json").write_text(
        json.dumps(config_dict, indent=2), encoding="utf-8"
    )
    make_validation_pdf(
        recovery,
        surface,
        coefficient,
        selection,
        output_dir / "validation_v5.pdf",
    )
    write_readme(output_dir, cfg, elapsed)
    print("Done.", flush=True)
    print(recovery.to_string(index=False), flush=True)
    print(f"Output directory: {output_dir.resolve()}", flush=True)
    return tables


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--preset",
        choices=("smoke", "focused", "robustness", "full"),
        default="smoke",
    )
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--seed", type=int, default=20260821)
    parser.add_argument("--no-selection", action="store_true")
    parser.add_argument("--no-hmc", action="store_true")
    parser.add_argument("--run-lowrank", action="store_true")
    parser.add_argument("--m", type=int, default=None)
    parser.add_argument("--p-annotations", type=int, default=None)
    parser.add_argument("--n-eff", type=int, default=None)
    parser.add_argument("--n-per-cell", type=int, default=None)
    parser.add_argument("--n-chains", type=int, default=None)
    parser.add_argument("--n-iter", type=int, default=None)
    parser.add_argument("--burn", type=int, default=None)
    parser.add_argument("--thin", type=int, default=None)
    parser.add_argument(
        "--annotation-scenarios",
        default=None,
        help="Comma-separated subset of null,sparse_independent,correlated_proxy,"
        "rare_binary,later_stick_only,dense",
    )
    parser.add_argument(
        "--genetic-scenarios",
        default=None,
        help="Comma-separated subset of bayesr_exact,major_polygenic,maf_dependent",
    )
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    cfg = apply_preset(Config(seed=args.seed), args.preset)
    overrides = {
        "m": args.m,
        "p_annotations": args.p_annotations,
        "n_eff": args.n_eff,
        "n_per_cell": args.n_per_cell,
        "n_chains": args.n_chains,
        "n_iter": args.n_iter,
        "burn": args.burn,
        "thin": args.thin,
    }
    for key, value in overrides.items():
        if value is not None:
            setattr(cfg, key, value)
    if args.no_selection:
        cfg.run_selection = False
    if args.no_hmc:
        cfg.run_hmc = False
    if args.run_lowrank:
        cfg.run_lowrank = True
    if args.annotation_scenarios:
        cfg.annotation_scenarios = tuple(
            x.strip() for x in args.annotation_scenarios.split(",") if x.strip()
        )
    if args.genetic_scenarios:
        cfg.genetic_scenarios = tuple(
            x.strip() for x in args.genetic_scenarios.split(",") if x.strip()
        )
    if cfg.burn >= cfg.n_iter:
        raise ValueError("burn must be smaller than n_iter")
    output_dir = Path(args.output_dir or f"sbayesrc_alpha_solution_v5_{cfg.preset}")
    run(cfg, output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
