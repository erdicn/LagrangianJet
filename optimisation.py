import os
import sys

# Suppress Wayland warnings and force XWayland before Qt loads
os.environ["QT_QPA_PLATFORM"] = "xcb"

import shutil
import subprocess
import warnings
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Process, Manager, Value
from scipy.interpolate import PchipInterpolator, interp1d
from sklearn.neural_network import MLPRegressor
from sklearn.exceptions import ConvergenceWarning

warnings.filterwarnings("ignore", category=ConvergenceWarning)

sys.path.insert(1, "/home/erdi/Uni/LCGAStage/Jets/LagrangianJet/DataVis")
from extracted_data_from_paper import exp_end111mm, exp_init11mm
from myprobs import MyProbs

# -------------------------------------------------------------
# Configuration (Scales up to 64 workers)
# -------------------------------------------------------------
SIMU_BINARY = os.path.abspath("./exe.exe")
N_WORKERS = 8  # Set to 16, 32, or 64 as needed
RADIUS_FILE_TO_METERS = 1e-6

# Diameter span in micrometers for C++ executable
D_GRID_UM = np.array(exp_init11mm.d_prob_dist_bin_centers * 1e6, copy=True)
D_MIN_UM, D_MAX_UM = float(D_GRID_UM.min()), float(D_GRID_UM.max())

# --- 1. Universal Reference Baseline: Gaussian centered at 120 µm ---
GAUSS_MU = 120.0
GAUSS_SIGMA = 50.0
raw_gaussian = np.exp(-0.5 * ((D_GRID_UM - GAUSS_MU) / GAUSS_SIGMA) ** 2)
PDF_INIT_UM = raw_gaussian / np.trapezoid(raw_gaussian, D_GRID_UM)

# Experimental target at 111 mm in SI meters
D_TARGET_M = np.array(exp_end111mm.d_prob_weighted_dist_bin_centers, copy=True)
PDF_VOL_TARGET = np.array(exp_end111mm.d_vol_weighted_prob_dist, copy=True)

VALID_TARGET_MASK = (PDF_VOL_TARGET > 0.0) & np.isfinite(PDF_VOL_TARGET)
D_TARGET_EVAL_M = D_TARGET_M[VALID_TARGET_MASK]
PDF_VOL_TARGET_EVAL = PDF_VOL_TARGET[VALID_TARGET_MASK]

# --- 2. Local Multiscale Parameters (PCHIP Knots + Ricker Wavelets) ---
N_KNOTS = 30  # High density feasible because PCHIP cannot oscillate/overshoot
N_WAVELETS = 8
TOTAL_PARAMS = N_KNOTS + N_WAVELETS
KNOT_X = np.geomspace(D_MIN_UM, D_MAX_UM, N_KNOTS)

# Localized wavelets centered across diameter spectrum
WAVELET_CENTERS = np.geomspace(D_MIN_UM, D_MAX_UM, N_WAVELETS)
WAVELET_SIGMAS = np.gradient(WAVELET_CENTERS) * 0.85

# Precomputed localized Ricker wavelet matrix (N_WAVELETS, len(D_GRID_UM))
WAVELET_BASIS = np.array([
    (1.0 - ((D_GRID_UM - WAVELET_CENTERS[k]) / WAVELET_SIGMAS[k]) ** 2)
    * np.exp(-0.5 * ((D_GRID_UM - WAVELET_CENTERS[k]) / WAVELET_SIGMAS[k]) ** 2)
    for k in range(N_WAVELETS)
])

# -------------------------------------------------------------
# PDF Construction, Projection & File Utilities
# -------------------------------------------------------------
def generate_valid_pdf_um(params):
    """
    PCHIP guarantees that flat parameters yield flat lines without overshoot,
    while isolated spikes yield narrow peaks without non-local wiggles.
    """
    alpha_knots = params[:N_KNOTS]
    alpha_wavelets = params[N_KNOTS:]

    # 1. Monotonicity-preserving Hermite interpolation
    pchip = PchipInterpolator(KNOT_X, alpha_knots)
    macro_perturb = pchip(D_GRID_UM)

    # 2. Localized wavelets only produce oscillations in their immediate band
    micro_perturb = np.dot(alpha_wavelets, WAVELET_BASIS)

    perturbation = np.clip(macro_perturb + micro_perturb, -4.0, 4.0)
    pdf_candidate = PDF_INIT_UM * np.exp(perturbation)

    area = np.trapezoid(pdf_candidate, D_GRID_UM)
    if area <= 0 or np.isnan(area):
        return PDF_INIT_UM
    return pdf_candidate / area


def project_pdf_to_params(target_pdf_um):
    """Maps any physical distribution into the parameter vector."""
    eps = 1e-12
    p_target = np.clip(target_pdf_um, eps, None)
    p_ref = np.clip(PDF_INIT_UM, eps, None)
    log_ratio = np.log(p_target / p_ref)

    interp_log = interp1d(D_GRID_UM, log_ratio, kind="linear", fill_value="extrapolate")
    alpha_knots = np.clip(interp_log(KNOT_X), -3.0, 3.0)
    alpha_wavelets = np.zeros(N_WAVELETS)
    return np.concatenate([alpha_knots, alpha_wavelets])


def load_best_pdf_file():
    for filename in ["distrib_data_best.txt", "distrib_data.txt"]:
        if os.path.exists(filename):
            try:
                with open(filename, "r") as f:
                    lines = [line.strip().replace(",", " ") for line in f if line.strip()]
                d_read = np.array([float(x) for x in lines[0].split()])
                pdf_read = np.array([float(y) for y in lines[1].split()])
                interp_f = interp1d(d_read, pdf_read, kind="linear", bounds_error=False, fill_value=1e-12)
                pdf_interp = interp_f(D_GRID_UM)
                area = np.trapezoid(pdf_interp, D_GRID_UM)
                if area > 0:
                    return pdf_interp / area
            except Exception:
                pass
    return None


def write_distrib_file_atomic(d_vals_um, pdf_vals_um, dest_folder):
    target = os.path.join(dest_folder, "distrib_data.txt")
    tmp = target + ".tmp"
    line1 = " ".join(f"{float(x):.6f}" for x in d_vals_um)
    line2 = " ".join(f"{float(y):.8e}" for y in pdf_vals_um)
    with open(tmp, "w") as f:
        f.write(line1 + "\n" + line2 + "\n")
    os.replace(tmp, target)


# -------------------------------------------------------------
# Cooperative Worker Process with Neural Network Surrogate
# -------------------------------------------------------------
def worker_process(worker_id, shared_memory, shared_best, total_evals, stop_flag, history_queue):
    os.environ["OMP_NUM_THREADS"] = "1"

    worker_dir = os.path.abspath(f"worker_{worker_id}")
    os.makedirs(worker_dir, exist_ok=True)

    local_exe = os.path.join(worker_dir, "exe.exe")
    if not os.path.exists(local_exe):
        try:
            os.symlink(SIMU_BINARY, local_exe)
        except OSError:
            shutil.copy2(SIMU_BINARY, local_exe)

    surrogate_nn = MLPRegressor(
        hidden_layer_sizes=(64, 32),
        activation="relu",
        solver="adam",
        learning_rate_init=0.005,
        max_iter=400,
        tol=1e-3,
        warm_start=True,
        random_state=42 + worker_id,
    )

    local_iter = 0

    def evaluate_simulation(params):
        nonlocal local_iter
        local_iter += 1

        with total_evals.get_lock():
            total_evals.value += 1
            curr_eval = total_evals.value

        candidate_pdf = generate_valid_pdf_um(params)
        write_distrib_file_atomic(D_GRID_UM, candidate_pdf, worker_dir)

        inj_rate = 29.25
        dt = 0.01
        dist_mm = 111
        nb_save = 10000
        base_folder = f"simu_run_{worker_id}_{local_iter}"
        real_out_dir = os.path.join(
            worker_dir,
            f"{base_folder}_inj{inj_rate:.4f}_dt{dt:.4f}_dist{dist_mm:.0f}_n{nb_save}",
        )

        cmd = [
            local_exe,
            base_folder,
            f"injection_rate={inj_rate}",
            f"dt={dt}",
            f"recording_distance_in_mm={dist_mm}",
            f"nb_droplets_saved={nb_save}",
        ]

        try:
            res = subprocess.run(cmd, cwd=worker_dir, capture_output=True, text=True)
            if res.returncode != 0:
                return 1e6

            sim_res = MyProbs.initFromFolder(
                folder_path=real_out_dir,
                label=f"w{worker_id}_{local_iter}",
                n_bins=120,
                diameters_to_si=RADIUS_FILE_TO_METERS,
                vel_to_si=1.0,
            )

            sim_d_m = sim_res.d_prob_weighted_dist_bin_centers
            sim_pdf_vol = sim_res.d_vol_weighted_prob_dist

            if sim_pdf_vol is None or len(sim_d_m) < 5:
                return 1e6

            try:
                history_queue.put_nowait(
                    (candidate_pdf, sim_d_m * 1e6, sim_pdf_vol)
                )
            except Exception:
                pass

            interp = interp1d(
                sim_d_m,
                sim_pdf_vol,
                kind="linear",
                bounds_error=False,
                fill_value=0.0,
            )
            sim_eval = interp(D_TARGET_EVAL_M)

            active_mask = sim_eval > 0.0
            if np.count_nonzero(active_mask) < 5:
                return 1e6

            log_sim = np.log10(sim_eval[active_mask])
            log_exp = np.log10(PDF_VOL_TARGET_EVAL[active_mask])
            loss = float(np.mean((log_sim - log_exp) ** 2))

            # Tail penalty (d > 350 µm must not exceed 10^-1)
            tail_mask = sim_d_m > 350e-6
            tail_vals = sim_pdf_vol[tail_mask]
            excess_mask = tail_vals > 1e-1

            if np.any(excess_mask):
                log_excess = np.log10(tail_vals[excess_mask]) - np.log10(1e-1)
                tail_penalty = 50.0 * float(np.sum(log_excess**2))
                loss += tail_penalty

            shared_memory.append((params.tolist(), loss))

            if loss < shared_best["loss"]:
                shared_best["loss"] = loss
                shared_best["params"] = params
                shared_best["pdf"] = candidate_pdf
                shared_best["sim_d"] = sim_d_m * 1e6
                shared_best["sim_vol"] = sim_pdf_vol
                shared_best["iter"] = curr_eval
                print(f"[Worker {worker_id:02d} | Eval {curr_eval:04d}] New Best Log MSE: {loss:.5f}")

            return loss

        except Exception:
            return 1e6

    # 1. Warm-up Phase: Configured per worker role
    if worker_id == 0:
        pdf_exp11 = (exp_init11mm.d_prob_density / 1e6)
        pdf_exp11 = pdf_exp11 / np.trapezoid(pdf_exp11, D_GRID_UM)
        seed_params = project_pdf_to_params(pdf_exp11)

    elif worker_id == 1:
        pdf_best = load_best_pdf_file()
        if pdf_best is not None:
            seed_params = project_pdf_to_params(pdf_best)
        else:
            pdf_exp11 = (exp_init11mm.d_prob_density / 1e6)
            pdf_exp11 = pdf_exp11 / np.trapezoid(pdf_exp11, D_GRID_UM)
            seed_params = project_pdf_to_params(pdf_exp11)

    else:
        rng = np.random.default_rng(2024 + worker_id)
        rand_mu = rng.normal(120.0, 3.0)
        rand_sigma = rng.uniform(20.0, 55.0)
        raw_g = np.exp(-0.5 * ((D_GRID_UM - rand_mu) / rand_sigma) ** 2)
        pdf_rand = raw_g / np.trapezoid(raw_g, D_GRID_UM)
        seed_params = project_pdf_to_params(pdf_rand)

    evaluate_simulation(seed_params)

    for _ in range(3):
        if stop_flag.value:
            return
        jitter = np.concatenate([
            np.random.normal(0, 0.12, N_KNOTS),
            np.random.normal(0, 0.05, N_WAVELETS),
        ])
        evaluate_simulation(seed_params + jitter)

    # 2. Neural Network Accelerated Active-Optimization Phase
    while not stop_flag.value:
        dataset = list(shared_memory)
        if len(dataset) < 10:
            continue

        X_train = np.array([item[0] for item in dataset])
        y_train = np.array([min(item[1], 10.0) for item in dataset])

        try:
            surrogate_nn.fit(X_train, y_train)
        except Exception:
            pass

        # --- Sparse, Heavy-Tailed Candidate Generation ---
        current_best = shared_best["params"]
        n_candidates = 350

        # Cauchy noise allows large targeted spikes; sparsity mask keeps flat zones flat
        raw_noise_knots = np.random.standard_cauchy(size=(n_candidates, N_KNOTS)) * 0.07
        sparse_mask_knots = np.random.binomial(1, p=0.35, size=(n_candidates, N_KNOTS))
        noise_macro = np.clip(raw_noise_knots * sparse_mask_knots, -0.6, 0.6)

        raw_noise_wavelets = np.random.normal(0, 0.08, size=(n_candidates, N_WAVELETS))
        sparse_mask_wavelets = np.random.binomial(1, p=0.30, size=(n_candidates, N_WAVELETS))
        noise_micro = raw_noise_wavelets * sparse_mask_wavelets

        candidates = current_best + np.hstack([noise_macro, noise_micro])

        predicted_losses = surrogate_nn.predict(candidates)
        best_candidate = candidates[np.argmin(predicted_losses)]

        evaluate_simulation(best_candidate)


# -------------------------------------------------------------
# Main Visualizer and Process Orchestrator
# -------------------------------------------------------------
if __name__ == "__main__":
    plt.ion()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    ax1.plot(
        D_GRID_UM,
        PDF_INIT_UM,
        "k--",
        lw=1.5,
        zorder=4,
        label="Gaussian Base (120±50 µm)",
    )
    (line_cand,) = ax1.plot(
        D_GRID_UM,
        PDF_INIT_UM,
        "blue",
        lw=2.2,
        zorder=5,
        label="Current Best PDF",
    )
    ax1.set_title("Input Number PDF to Executable ($L=11$ mm)")
    ax1.set_xlabel("Diameter [µm]")
    ax1.set_ylabel("Probability Density [1/µm]")
    ax1.legend(loc="upper right")
    ax1.grid(True, ls=":", alpha=0.5)

    ax2.plot(
        D_TARGET_M * 1e6,
        PDF_VOL_TARGET,
        "r-",
        lw=2.2,
        zorder=4,
        label="Exp Target ($L=111$ mm)",
    )
    (line_sim,) = ax2.plot(
        [],
        [],
        color="darkblue",
        marker=".",
        lw=1.5,
        zorder=5,
        label="Current Best Sim",
    )
    ax2.set_yscale("log")
    ax2.set_title("Volume-Weighted PDF at $L=111$ mm")
    ax2.set_xlabel("Diameter [µm]")
    ax2.set_ylabel("Vol-Weighted PDF")
    ax2.set_ylim(bottom=1e-9, top=np.nanmax(PDF_VOL_TARGET) * 3)
    ax2.set_xlim(left=0, right=np.max(D_TARGET_M * 1e6) * 1.05)
    ax2.legend(loc="upper right")
    ax2.grid(True, ls=":", alpha=0.5)
    plt.tight_layout()

    manager = Manager()
    shared_memory = manager.list()
    history_queue = manager.Queue(maxsize=5000)
    shared_best = manager.dict(
        {
            "loss": 1e6,
            "params": np.zeros(TOTAL_PARAMS),
            "pdf": PDF_INIT_UM,
            "sim_d": np.array([]),
            "sim_vol": np.array([]),
            "iter": 0,
        }
    )
    total_evals = Value("i", 0)
    stop_flag = Value("b", False)

    workers = [
        Process(
            target=worker_process,
            args=(
                i,
                shared_memory,
                shared_best,
                total_evals,
                stop_flag,
                history_queue,
            ),
        )
        for i in range(N_WORKERS)
    ]

    print(f"Launching {N_WORKERS} workers with PCHIP and localized wavelets...")
    for w in workers:
        w.start()

    try:
        while any(w.is_alive() for w in workers):
            drain_count = 0
            while not history_queue.empty() and drain_count < 100:
                try:
                    cand_pdf_hist, sim_d_hist, sim_vol_hist = history_queue.get_nowait()

                    ax1.plot(
                        D_GRID_UM,
                        cand_pdf_hist,
                        color="cornflowerblue",
                        alpha=0.04,
                        lw=0.7,
                        zorder=1,
                    )

                    if len(sim_d_hist) > 0:
                        ax2.plot(
                            sim_d_hist,
                            sim_vol_hist,
                            color="gray",
                            alpha=0.05,
                            lw=0.7,
                            zorder=1,
                        )

                    drain_count += 1
                except Exception:
                    break

            if shared_best["loss"] < 1e5:
                line_cand.set_ydata(shared_best["pdf"])
                ax1.set_ylim(
                    0,
                    max(
                        np.max(shared_best["pdf"]) * 1.18,
                        np.max(PDF_INIT_UM) * 1.18,
                    ),
                )

                if len(shared_best["sim_d"]) > 0:
                    line_sim.set_data(
                        shared_best["sim_d"], shared_best["sim_vol"]
                    )

                fig.suptitle(
                    f"Collective Evals: {total_evals.value:04d} | Replay Pool Size: {len(shared_memory)} | "
                    f"Best Support-Masked Log MSE: {shared_best['loss']:.5f}",
                    fontsize=12,
                )

            fig.canvas.draw_idle()
            fig.canvas.flush_events()
            fig.canvas.start_event_loop(0.4)

    except KeyboardInterrupt:
        print("\nStopping worker processes...")
        stop_flag.value = True
        for w in workers:
            w.terminate()

    for w in workers:
        w.join()

    best_pdf = shared_best["pdf"]
    write_distrib_file_atomic(D_GRID_UM, best_pdf, ".")
    shutil.copy2("distrib_data.txt", "distrib_data_best.txt")
    print(f"Optimization complete. Best PDF saved to distrib_data_best.txt (Loss: {shared_best['loss']:.5f})")

    plt.ioff()
    plt.show()