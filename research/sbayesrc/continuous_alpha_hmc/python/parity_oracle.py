#!/usr/bin/env python3
"""Expose deterministic v5 quantities for the independent R parity checks."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from scipy.special import log_ndtr, logsumexp

from sbayesrc_alpha_solution_v5 import (
    CollapsedModel,
    build_designs,
    log_pi_from_theta,
    responsibilities,
)


def read_matrix(path: Path) -> np.ndarray:
    value = np.loadtxt(path, delimiter=",", ndmin=2)
    return np.asarray(value, dtype=float)


def write_matrix(path: Path, value: np.ndarray) -> None:
    np.savetxt(path, np.asarray(value), delimiter=",", fmt="%.17g")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir")
    parser.add_argument("output_dir")
    args = parser.parse_args()
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    annotation = read_matrix(input_dir / "annotation.csv")
    centers = read_matrix(input_dir / "centers.csv")
    theta = read_matrix(input_dir / "theta.csv")
    prior_mean = read_matrix(input_dir / "prior_mean.csv")
    prior_sd = read_matrix(input_dir / "prior_sd.csv")
    summary = read_matrix(input_dir / "summary.csv")
    mixture = read_matrix(input_dir / "mixture.csv").ravel()
    supplied_mode = read_matrix(input_dir / "mode_r.csv")
    supplied_covariance = read_matrix(input_dir / "covariance_r.csv")
    supplied_position = read_matrix(input_dir / "position.csv").ravel()
    supplied_momentum = read_matrix(input_dir / "momentum.csv").ravel()
    hmc_control = read_matrix(input_dir / "hmc_control.csv").ravel()

    v_beta = mixture[0]
    gamma = mixture[1:]
    variance = summary[:, 1, None] ** 2 + v_beta * gamma[None, :]
    log_density = -0.5 * (
        np.log(2.0 * np.pi * variance) + summary[:, 0, None] ** 2 / variance
    )
    A = np.column_stack((np.ones(annotation.shape[0]), annotation))
    designs = build_designs(A, centers)
    model = CollapsedModel(designs, log_density, prior_mean, prior_sd)

    log_pi, q, eta = log_pi_from_theta(theta, designs)
    responsibility, _, _, total_log_likelihood = responsibilities(
        theta, designs, log_density
    )
    log_weight = log_pi + log_density
    marker_log_likelihood = logsumexp(log_weight, axis=1)
    posterior_kernel, gradient = model.logpost_and_grad(theta.ravel())
    prior_kernel = -0.5 * np.sum(((theta - prior_mean) / prior_sd) ** 2)
    prior_constant = -np.log(prior_sd).sum() - 0.5 * theta.size * np.log(
        2.0 * np.pi
    )
    posterior_mode, python_covariance = model.laplace()

    root = np.linalg.cholesky(supplied_covariance)
    supplied_mode_flat = supplied_mode.ravel()
    theta_flat = theta.ravel()
    whitened = np.linalg.solve(root, theta_flat - supplied_mode_flat)
    inverse_whitened = supplied_mode_flat + root @ whitened

    def target(position: np.ndarray):
        flat = supplied_mode_flat + root @ position
        value, gradient_theta = model.logpost_and_grad(flat)
        full_value = value + prior_constant
        return full_value, root.T @ gradient_theta

    step_size = float(hmc_control[0])
    n_step = int(hmc_control[1])
    initial_value, initial_gradient = target(supplied_position)
    position = supplied_position.copy()
    momentum = supplied_momentum.copy()
    momentum += 0.5 * step_size * initial_gradient
    proposed_value = initial_value
    proposed_gradient = initial_gradient
    for leapfrog in range(n_step):
        position += step_size * momentum
        proposed_value, proposed_gradient = target(position)
        if leapfrog != n_step - 1:
            momentum += step_size * proposed_gradient
    momentum += 0.5 * step_size * proposed_gradient
    current_hamiltonian = -initial_value + 0.5 * np.dot(
        supplied_momentum, supplied_momentum
    )
    proposed_hamiltonian = -proposed_value + 0.5 * np.dot(momentum, momentum)
    log_acceptance_ratio = current_hamiltonian - proposed_hamiltonian
    acceptance_probability = np.exp(min(0.0, log_acceptance_ratio))

    outputs = {
        "linear_predictor.csv": eta,
        "continuation_probability.csv": q,
        "log_continuation_probability.csv": log_ndtr(eta),
        "log_continuation_survival_probability.csv": log_ndtr(-eta),
        "log_component_probability.csv": log_pi,
        "component_probability.csv": np.exp(log_pi),
        "component_log_density.csv": log_density,
        "component_log_weight.csv": log_weight,
        "marker_log_likelihood.csv": marker_log_likelihood[:, None],
        "responsibility.csv": responsibility,
        "gradient.csv": gradient.reshape(theta.shape),
        "mode_python.csv": posterior_mode,
        "covariance_python.csv": python_covariance,
        "whitened.csv": whitened[:, None],
        "inverse_whitened.csv": inverse_whitened[:, None],
        "leapfrog_position.csv": position[:, None],
        "leapfrog_momentum.csv": momentum[:, None],
        "scalars.csv": np.array(
            [
                total_log_likelihood,
                prior_kernel,
                prior_kernel + prior_constant,
                posterior_kernel,
                posterior_kernel + prior_constant,
                current_hamiltonian,
                proposed_hamiltonian,
                log_acceptance_ratio,
                acceptance_probability,
            ]
        )[:, None],
    }
    for filename, value in outputs.items():
        write_matrix(output_dir / filename, value)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
