# Continuous Bidirectional Pore Evolution Diffusion Workflow

This folder prepares RTSPHEM interface-image trajectories for diffusion/BBDM
models. It uses `exp_*/interface_images/timestep_####.png` images, removes the
white border, keeps states with global `Porosity <= 0.80` by default, and forms
directed pairs inside each RTM experiment.

## Build Dataset

```powershell
python ReactiveTransport/diffusion_inverse/build_interface_evolution_dataset.py `
  --batch-root outputs/rtm_batches/batch_20260508_074936 `
  --output-root outputs/diffusion_inverse/interface_bidirectional_p080 `
  --porosity-max 0.80 `
  --image-size 256 `
  --max-pairs-per-exp 2000
```

Use `--max-pairs-per-exp 0` only for small batches; directed pairs grow as
`N * (N - 1)` per experiment.

## Train Time Estimator

```powershell
python ReactiveTransport/diffusion_inverse/fit_time_estimator.py `
  --pairs-csv outputs/diffusion_inverse/interface_bidirectional_p080/manifests/pairs.csv `
  --model-path outputs/diffusion_inverse/interface_bidirectional_p080/time_estimator.json `
  --group-metrics-csv outputs/diffusion_inverse/interface_bidirectional_p080/qc/time_estimator_by_da_pe.csv
```

The estimator predicts target `Time_s` from `Da`, `Pe`, layout, source progress,
target progress, signed progress delta, and source time.

## Prepare BBDM Training

```powershell
python ReactiveTransport/diffusion_inverse/train_bidirectional_bbdm.py `
  --dataset-root outputs/diffusion_inverse/interface_bidirectional_p080 `
  --output-config outputs/diffusion_inverse/interface_bidirectional_p080/bbdm/bidirectional_bbdm_config.yaml
```

This writes an RTSPHEM-side config. The official BBDM code stays external and is
passed through `--bbdm-root` when launching a training command.

## Train Local Neural Baselines

Use the CUDA PyTorch environment that was found during the full run:

```powershell
& C:\Users\imgw\Documents\Codex\NMR-agent\.venv\Scripts\python.exe `
  ReactiveTransport/diffusion_inverse/train_torch_unet_baseline.py `
  --train-csv outputs/diffusion_inverse/interface_bidirectional_p080_full/splits/train.csv `
  --test-csv outputs/diffusion_inverse/interface_bidirectional_p080_full/splits/test.csv `
  --checkpoint outputs/diffusion_inverse/interface_bidirectional_p080_full/models/unet_64_50k_2000.pt `
  --generated-dir outputs/diffusion_inverse/interface_bidirectional_p080_full/generated/unet_64_50k_2000/test `
  --predictions-csv outputs/diffusion_inverse/interface_bidirectional_p080_full/generated/unet_64_50k_2000/test_predictions.csv `
  --image-size 64 --batch-size 32 --max-train-pairs 50000 --max-steps 2000 --epochs 2 --device cuda
```

```powershell
& C:\Users\imgw\Documents\Codex\NMR-agent\.venv\Scripts\python.exe `
  ReactiveTransport/diffusion_inverse/train_torch_ddpm_baseline.py `
  --train-csv outputs/diffusion_inverse/interface_bidirectional_p080_full/splits/train.csv `
  --test-csv outputs/diffusion_inverse/interface_bidirectional_p080_full/splits/test.csv `
  --checkpoint outputs/diffusion_inverse/interface_bidirectional_p080_full/models/ddpm_64_50k_4000.pt `
  --generated-dir outputs/diffusion_inverse/interface_bidirectional_p080_full/generated/ddpm_64_50k_4000/test `
  --predictions-csv outputs/diffusion_inverse/interface_bidirectional_p080_full/generated/ddpm_64_50k_4000/test_predictions.csv `
  --image-size 64 --batch-size 32 --max-train-pairs 50000 --max-steps 4000 `
  --epochs 3 --timesteps 100 --sample-steps 20 --device cuda
```

## Train V2 With OpenAI Improved-Diffusion

Clone the official repository outside RTSPHEM:

```powershell
git clone https://github.com/openai/improved-diffusion.git C:\Users\imgw\Documents\Codex\openai-improved-diffusion
```

Run the RTSPHEM adapter. The model and diffusion process come from the external
`improved_diffusion` package; the local script only adapts source/target RTM
pairs into a conditional image-to-image task.

```powershell
& C:\Users\imgw\Documents\Codex\NMR-agent\.venv\Scripts\python.exe `
  ReactiveTransport/diffusion_inverse/train_openai_improved_diffusion.py `
  --improved-diffusion-root C:\Users\imgw\Documents\Codex\openai-improved-diffusion `
  --train-csv outputs/diffusion_inverse/interface_bidirectional_p080_full/splits/train.csv `
  --test-csv outputs/diffusion_inverse/interface_bidirectional_p080_full/splits/test.csv `
  --checkpoint outputs/diffusion_inverse/interface_bidirectional_p080_full/models/openai_improved_64_50k_4000.pt `
  --generated-dir outputs/diffusion_inverse/interface_bidirectional_p080_full/generated/openai_improved_64_50k_4000/test `
  --predictions-csv outputs/diffusion_inverse/interface_bidirectional_p080_full/generated/openai_improved_64_50k_4000/test_predictions.csv `
  --image-size 64 --num-channels 32 --num-res-blocks 2 --attention-resolutions 16,8 `
  --batch-size 32 --max-train-pairs 50000 --max-steps 4000 --epochs 3 `
  --diffusion-steps 100 --timestep-respacing ddim20 --use-ddim --device cuda
```

## Selected Forward/Backward Test Cases

Create the requested three-geometry by three-parameter test manifest:

```powershell
python ReactiveTransport/diffusion_inverse/make_test_case_comparisons.py `
  --dataset-root outputs/diffusion_inverse/interface_bidirectional_p080_full `
  --output-dir outputs/diffusion_inverse/interface_bidirectional_p080_full/selected_test_cases
```

Run the tuned OpenAI improved-diffusion checkpoint on those selected cases:

```powershell
& C:\Users\imgw\Documents\Codex\NMR-agent\.venv\Scripts\python.exe `
  ReactiveTransport/diffusion_inverse/train_openai_improved_diffusion.py `
  --improved-diffusion-root C:\Users\imgw\Documents\Codex\openai-improved-diffusion `
  --train-csv outputs/diffusion_inverse/interface_bidirectional_p080_full/splits/train.csv `
  --test-csv outputs/diffusion_inverse/interface_bidirectional_p080_full/selected_test_cases/selected_pairs.csv `
  --checkpoint outputs/diffusion_inverse/interface_bidirectional_p080_full/tuned/openai_improved_ch64_1000/checkpoint.pt `
  --generated-dir outputs/diffusion_inverse/interface_bidirectional_p080_full/selected_test_cases/generated/openai_improved_ch64_1000 `
  --predictions-csv outputs/diffusion_inverse/interface_bidirectional_p080_full/selected_test_cases/openai_improved_ch64_1000_predictions.csv `
  --batch-size 8 --predict-only --use-ddim --timestep-respacing ddim20 --device cuda
```

Assemble real-vs-generated figures and metrics:

```powershell
python ReactiveTransport/diffusion_inverse/make_test_case_comparisons.py `
  --dataset-root outputs/diffusion_inverse/interface_bidirectional_p080_full `
  --output-dir outputs/diffusion_inverse/interface_bidirectional_p080_full/selected_test_cases `
  --generated-dir outputs/diffusion_inverse/interface_bidirectional_p080_full/selected_test_cases/generated/openai_improved_ch64_1000
```

## Query A Trajectory

```powershell
python ReactiveTransport/diffusion_inverse/sample_interface_trajectory.py `
  --dataset-root outputs/diffusion_inverse/interface_bidirectional_p080 `
  --exp-id exp_230 `
  --source-timestep 8 `
  --target-progresses 0.05,0.20,0.60 `
  --time-estimator outputs/diffusion_inverse/interface_bidirectional_p080/time_estimator.json `
  --output-dir outputs/diffusion_inverse/interface_bidirectional_p080/queries/exp_230_t0008
```

Without generated BBDM images, the query grid falls back to the nearest RTM
reference image so the requested progress/time mapping can still be inspected.

## Evaluate Generated Images

```powershell
python ReactiveTransport/diffusion_inverse/evaluate_interface_trajectory.py `
  --pairs-csv outputs/diffusion_inverse/interface_bidirectional_p080/splits/test.csv `
  --generated-dir outputs/diffusion_inverse/interface_bidirectional_p080/generated/test `
  --output-csv outputs/diffusion_inverse/interface_bidirectional_p080/qc/generated_test_metrics.csv
```
