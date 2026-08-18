#!/usr/bin/env python3
"""Small deterministic checks for sbayesrc_alpha_solution_v5.py."""

import unittest

import numpy as np

from sbayesrc_alpha_solution_v5 import (
    CollapsedModel,
    Config,
    alpha_to_theta,
    build_designs,
    marker_log_densities,
    log_pi_from_theta,
    responsibilities,
    simulate_dataset,
    stick_to_pi,
    theta_to_alpha,
)


class TestAlphaSolutionV5(unittest.TestCase):
    def setUp(self):
        self.cfg = Config(m=80, p_annotations=8, n_eff=4000)
        self.data = simulate_dataset(
            self.cfg,
            "sparse_independent",
            "bayesr_exact",
            replicate=1,
            seed=91,
        )

    def test_probabilities_sum_to_one(self):
        pi = stick_to_pi(self.data.q_true)
        self.assertLess(np.max(np.abs(pi.sum(axis=1) - 1.0)), 1e-12)
        self.assertTrue(np.all(pi > 0.0))

    def test_centered_transform_preserves_predictor(self):
        rng = np.random.default_rng(12)
        centers = rng.normal(scale=0.2, size=(3, self.cfg.p_annotations))
        theta = alpha_to_theta(self.data.alpha_true, centers)
        recovered = theta_to_alpha(theta, centers)
        self.assertLess(np.max(np.abs(recovered - self.data.alpha_true)), 1e-12)
        designs = build_designs(self.data.A, centers)
        for stick in range(3):
            lhs = designs[stick] @ theta[stick]
            rhs = self.data.A @ self.data.alpha_true[:, stick]
            self.assertLess(np.max(np.abs(lhs - rhs)), 1e-12)

    def test_collapsed_gradient(self):
        centers = np.zeros((3, self.cfg.p_annotations))
        designs = build_designs(self.data.A, centers)
        d = self.data.A.shape[1]
        mean = np.zeros((3, d))
        mean[:, 0] = (-1.0, 0.0, -0.1)
        sd = np.full((3, d), 0.75)
        sd[:, 0] = 1.5
        model = CollapsedModel(designs, marker_log_densities(self.data), mean, sd)
        rng = np.random.default_rng(3)
        theta = mean + rng.normal(scale=0.1, size=mean.shape)
        value, gradient = model.logpost_and_grad(theta.ravel())
        numerical = np.zeros_like(gradient)
        step = 1e-6
        for j in range(gradient.size):
            plus = theta.ravel().copy()
            minus = theta.ravel().copy()
            plus[j] += step
            minus[j] -= step
            numerical[j] = (
                model.logpost_and_grad(plus)[0]
                - model.logpost_and_grad(minus)[0]
            ) / (2.0 * step)
        self.assertTrue(np.isfinite(value))
        self.assertLess(np.max(np.abs(gradient - numerical)), 2e-5)

    def test_exact_probit_tail_gradient(self):
        designs = [np.ones((2, 1), dtype=float) for _ in range(3)]
        step = 5e-5
        for tail in (-8.0, 8.0):
            summary_beta = np.full(2, 0.3 if tail < 0 else 0.0)
            summary_se = np.full(2, 0.01 if tail < 0 else 1e-8)
            variance = summary_se[:, None] ** 2 + 0.08 * np.array(
                [0.0, 0.01, 0.10, 1.0]
            )[None, :]
            log_density = -0.5 * (
                np.log(2.0 * np.pi * variance)
                + summary_beta[:, None] ** 2 / variance
            )
            model = CollapsedModel(
                designs, log_density, np.zeros((3, 1)), np.full((3, 1), 1.5)
            )
            flat = np.full(3, tail)
            value, analytic = model.logpost_and_grad(flat)
            numerical = np.zeros_like(analytic)
            for j in range(flat.size):
                plus_two = flat.copy()
                plus_one = flat.copy()
                minus_one = flat.copy()
                minus_two = flat.copy()
                plus_two[j] += 2.0 * step
                plus_one[j] += step
                minus_one[j] -= step
                minus_two[j] -= 2.0 * step
                numerical[j] = (
                    -model.logpost_and_grad(plus_two)[0]
                    + 8.0 * model.logpost_and_grad(plus_one)[0]
                    - 8.0 * model.logpost_and_grad(minus_one)[0]
                    + model.logpost_and_grad(minus_two)[0]
                ) / (12.0 * step)
            self.assertTrue(np.isfinite(value))
            self.assertLess(np.max(np.abs(analytic - numerical)), 2e-7)
            change = model.loglik((flat + step).reshape(model.shape)) - model.loglik(
                (flat - step).reshape(model.shape)
            )
            self.assertNotEqual(change, 0.0)

    def test_size_one_axes_are_preserved(self):
        A = np.array([[1.0, 0.25]])
        for n_stick in (3, 1):
            centers = np.zeros((n_stick, 1))
            designs = build_designs(A, centers)
            theta = np.zeros((n_stick, 2))
            log_density = np.linspace(-2.0, -0.5, n_stick + 1)[None, :]
            log_pi, q, eta = log_pi_from_theta(theta, designs)
            r, _, _, _ = responsibilities(theta, designs, log_density)
            self.assertEqual(eta.shape, (1, n_stick))
            self.assertEqual(q.shape, (1, n_stick))
            self.assertEqual(log_pi.shape, (1, n_stick + 1))
            self.assertEqual(r.shape, (1, n_stick + 1))
            self.assertAlmostEqual(float(np.exp(log_pi).sum()), 1.0, places=14)
            self.assertAlmostEqual(float(r.sum()), 1.0, places=14)


if __name__ == "__main__":
    unittest.main()
