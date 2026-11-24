#!/usr/bin/env bash
# Install the frozen Helixer release in an isolated environment on Linux.
set -euo pipefail
root=$1
software=$root/software
env=$software/env
mkdir -p "$software" "$root/provenance"
[[ -x $env/bin/python ]] || conda create -y -p "$env" -c conda-forge --override-channels python=3.10 hdf5=1.14 rust pkg-config
export PATH="$env/bin:$PATH" HDF5_DIR="$env" CARGO_BUILD_JOBS=4
export LD_LIBRARY_PATH="$env/lib:${LD_LIBRARY_PATH:-}"
[[ ! -d $software/wheels ]] || export PIP_FIND_LINKS="$software/wheels"
for item in Helixer:v0.3.8 HelixerPost:v0.4.0 GeenuFF:v0.4.0; do
  repo=${item%:*}; tag=${item#*:}
  [[ -d $software/$repo/.git ]] || git clone --depth 1 --branch "$tag" "https://github.com/usadellab/$repo.git" "$software/$repo"
  git -C "$software/$repo" describe --tags --exact-match | grep -Fx "$tag"
  git -C "$software/$repo" rev-parse HEAD > "$root/provenance/$repo.commit"
done
python -m pip install 'pip @ https://files.pythonhosted.org/packages/3f/06/3e4420302363fbb71be3c7f35ce042f54bfca330e3e66b3a2e3c4cca579f/pip-26.2.2-py3-none-any.whl#sha256=fd7f4dfa8190e014013b9f08f3e973f4e4c26ffd179a204895533a7254813a70'
if [[ -d $software/wheels ]]; then
  python -m pip install --no-index --no-deps "$software/wheels/"*.whl
fi
python -m pip install -r "$software/Helixer/requirements.3.10.txt" 'tensorflow[and-cuda]==2.15.1'
python -m pip install "$software/GeenuFF"
python -m pip install --no-deps "$software/Helixer"
python -m pip check
lockfile=$(dirname "$0")/../metadata/helixer_Cargo.lock
[[ ! -f $lockfile ]] || cp "$lockfile" "$software/HelixerPost/Cargo.lock"
(cd "$software/HelixerPost" && { [[ -f Cargo.lock ]] || cargo generate-lockfile; } && cargo build --release --locked)
cp "$software/HelixerPost/Cargo.lock" "$root/provenance/HelixerPost.Cargo.lock"
cp "$software/HelixerPost/target/release/helixer_post_bin" "$env/bin/"
helixer_post_bin --version
model=$software/models/fungi/fungi_v0.3_a_0100.h5
mkdir -p "$(dirname "$model")"
[[ -s $model ]] || curl -fLsS --retry 3 'https://zenodo.org/records/10836346/files/fungi_v0.3_a_0100.h5?download=1' -o "$model"
printf 'c62d9615e722c2f927a0de3f4fcf5b70  %s\n' "$model" | md5sum -c -
python -m pip freeze > "$root/provenance/python_packages.txt"
sha256sum "$software/models/fungi/fungi_v0.3_a_0100.h5" > "$root/provenance/model.sha256"
printf 'Helixer 0.3.8; HelixerPost 0.4.0; GeenuFF 0.4.0\n' > "$root/SOFTWARE_READY"
