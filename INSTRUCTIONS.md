# Model Instructions

This file summarizes the key steps for compiling and running the malaria in pregnancy model. For further detail on parameters and advanced usage see `Instructions.docx`.

## Compiling the executable
1. Open Microsoft Visual Studio (the Community Edition works).
2. Create a new **Empty Project** and choose a project name and location.
3. Copy the contents of the `code` directory from this repository into your new project directory.
4. In **Solution Explorer**, right-click **Source Files** and choose **Add &rarr; Existing Item…**. Select `main.cpp` from the copied `code` folder.
5. Select **Release** from the build configuration dropdown (avoids slower debug builds).
6. Build the project. The resulting executable will appear in the project `Release` folder (e.g. `C:\<ProjectName>\Release\<ProjectName>.exe`).
7. Move the compiled `.exe` to the folder where you will run the model.

## Running the model
The simplest approach is to create a batch file. `compiled/run_example.bat` shows the format:

```
<NameOfExecutable>.exe <root_path> <parameter_dir> <output_name>
```

Arguments:
1. **Executable path** – path to the compiled `.exe` and a root path for model inputs and outputs.
2. **Parameter directory** – name of the directory containing the parameter text files.
3. **Output name** – name (and optional path) for the generated output files.

Additional parameter overrides can be supplied as pairs after these three arguments, e.g.

```
<NameOfExecutable>.exe <root> <dir> <out> EIR 5 HB_eval_time 200
```

## Parameter files
A parameter directory must contain the following text files:

- `sim_params` – main simulation parameters you may wish to modify.
- `preg_params` – pregnancy-specific parameters (posterior medians; usually kept as provided).
- `fertility_rates` – 7&times;20 matrix of age- and parity-specific fertility rates.
- `HB_non_inf_file` – mean haemoglobin levels by gestational age and gravidity.
- `HB_inf_file` – impact of malaria on haemoglobin by gestational age and prior infection history.

## Interpreting output
Running the model produces files prefixed by your chosen output name. Of particular interest:

- `<output_name>.cout` – record of all options supplied, useful for troubleshooting.
- `<output_name>_hb_summary.txt` – summary of anaemia-related outputs with columns such as `Grav_Cat`, `EIR`, `primi_prev`, `HB_diff`, `moderate_anaemia`, and `severe_anaemia`.

To generate infection histories instead of summary outputs, set `inf_history` to `1` in `sim_params` (this overrides any ANC strategies).

For further detail on parameters and advanced usage see `Instructions.docx`.
