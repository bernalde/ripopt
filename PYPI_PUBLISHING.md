# PyPI Publishing Setup

ripopt ships two Python packages to PyPI, both gated on a `v*` git tag:

| PyPI name      | Source                            | Workflow                              | Contents                                         |
|----------------|-----------------------------------|---------------------------------------|--------------------------------------------------|
| `pyomo-ripopt` | `pyomo-ripopt/`                   | `.github/workflows/publish-pyomo.yml` | Pyomo plugin + bundled `ripopt` solver binary    |
| `ripopt`       | `ripopt-py/`                      | `.github/workflows/publish-ripopt-py.yml` | Direct PyO3 bindings (abi3 wheel, Python 3.9+) |

Both workflows use **PyPI Trusted Publishing** (OIDC, no API tokens). That
requires one-time configuration on PyPI, described below.

---

## One-time PyPI setup (per project)

Do this **before** cutting the first release that should land on PyPI.

### 1. Create the project on PyPI

If the project name is still free, reserve it first so nobody else takes it:

```bash
# Check availability
curl -sI https://pypi.org/pypi/ripopt/json      | head -1   # expect 404
curl -sI https://pypi.org/pypi/pyomo-ripopt/json | head -1   # expect 404
```

PyPI doesn't let you "reserve" a name without uploading something, so either:
- (a) Upload a `0.0.0` placeholder manually once (`python -m build && twine
  upload dist/*`) with a temporary token, then configure trusted publishing
  against the now-existing project, or
- (b) Configure a **pending** trusted publisher (PyPI UI → Your Account →
  Publishing → "Add a new pending publisher"). The first successful OIDC
  publish creates the project and binds the publisher.

Option (b) is cleaner and avoids having to mint / rotate a manual token.

### 2. Configure the trusted publisher on PyPI

For each project (`pyomo-ripopt` and `ripopt`), go to:

- **Existing project:** https://pypi.org/manage/project/\<name\>/settings/publishing/
- **Pending (new) project:** https://pypi.org/manage/account/publishing/

Click **Add** and fill in exactly:

| Field                  | `pyomo-ripopt`              | `ripopt`                       |
|------------------------|-----------------------------|--------------------------------|
| PyPI project name      | `pyomo-ripopt`              | `ripopt`                       |
| Owner                  | `jkitchin`                  | `jkitchin`                     |
| Repository name        | `ripopt`                    | `ripopt`                       |
| Workflow filename      | `publish-pyomo.yml`         | `publish-ripopt-py.yml`        |
| Environment            | `pypi`                      | `pypi`                         |

The workflow filename must be **exactly** the basename under
`.github/workflows/` — not a path, not `publish-pyomo`, not with a leading
slash. Mismatches here produce `invalid-publisher` errors that look
cryptic in the Action log.

### 3. Create the `pypi` GitHub environment

In the GitHub repo:

1. Settings → Environments → **New environment** → name it `pypi`.
2. (Optional but recommended) Under **Deployment branches and tags**, pick
   **Selected branches and tags** → Add deployment tag rule → `v*`. That
   restricts OIDC publish to tag pushes matching the workflow trigger.
3. No secrets need to be stored in this environment — Trusted Publishing
   uses short-lived OIDC tokens, not long-lived API keys.

Both workflows already reference this environment:

```yaml
publish:
  environment: pypi
  permissions:
    id-token: write   # required for OIDC
```

### 4. Dry-run before tagging

Before the next real release, sanity-check the setup:

```bash
# pyomo-ripopt: builds the ripopt binary, bundles it, builds a wheel
gh workflow run publish-pyomo.yml --ref main        # will fail the publish step
                                                     # (no tag) but exercises the build matrix

# ripopt-py: maturin wheel build across 5 platforms
gh workflow run publish-ripopt-py.yml --ref main
```

Watch for two things:

- **All five wheel jobs succeed.** The `macos-13` runner dispatch failure
  that killed v0.6.2 should be fixed now — both workflows cross-compile
  x86_64 macOS from `macos-14` (Apple silicon).
- **The publish step is skipped**, not failed. A skipped publish on a
  non-tag branch is the expected outcome of a dry-run.

---

## Per-release flow

Once Trusted Publishing is configured, publishing is automatic:

1. Bump versions per `RELEASE_CHECKLIST.md` §5 (includes
   `pyomo-ripopt/pyproject.toml`, `ripopt-py/Cargo.toml`,
   `ripopt-py/pyproject.toml`).
2. Commit, tag `vX.Y.Z`, push tag.
3. GitHub Actions builds wheels + sdist (pyomo) / wheels only (ripopt-py)
   across all platforms in parallel.
4. Both workflows' final `publish` job uploads to PyPI via OIDC.
5. Within a minute of workflow completion:

   ```bash
   pip install pyomo-ripopt==X.Y.Z
   pip install ripopt==X.Y.Z
   ```

   should both resolve on linux (x86_64/aarch64), macOS (x86_64/arm64), and
   windows (x64).

---

## Troubleshooting

**`invalid-publisher: valid token, but no corresponding publisher`**
Trusted publisher config on PyPI doesn't match the workflow. Triple-check:
repo owner, repo name, workflow filename (case-sensitive, basename only),
and environment name (`pypi`).

**macOS x86_64 wheel fails with `configuration '…-us-default' is not supported`**
GitHub runner dispatch issue on the old `macos-13` image. Both workflows
now cross-compile x86_64 from `macos-14`, which avoids it. If it recurs,
simply re-run the failed job from the Actions UI.

**`publish` job runs on a non-tag push**
The workflow trigger is `on.push.tags: ['v*']`, so the whole workflow is
skipped on branch pushes. If it runs anyway, check you haven't accidentally
added a `workflow_dispatch` path that defaults to publishing.

**Wheel builds but `pip install` errors with `No matching distribution`**
PyPI needed a few seconds to index the upload; wait 30s and retry. If it
persists, check the `publish` job log for an HTTP error during upload.

**`abi3` wheel is wrongly tagged `cp312-cp312` instead of `cp39-abi3`**
`ripopt-py/Cargo.toml` must have `pyo3 = { features = ["abi3-py39", …] }`
and `ripopt-py/pyproject.toml` must have
`features = ["pyo3/extension-module", "pyo3/abi3-py39"]` under
`[tool.maturin]`. Both are required; maturin picks the tag from the Cargo
features, not the pyproject alone.

---

## Why two separate workflows?

`pyomo-ripopt` is a pure-Python package that bundles a **prebuilt `ripopt`
binary** (produced by `cargo build --release --bin ripopt`). It needs a
matrix across platforms only to produce the right binary per wheel.

`ripopt-py` is a **native PyO3 extension** compiled by maturin against
abi3-py39, so one wheel per (platform, arch) covers Python 3.9 onward.

Keeping them in separate workflows means a build failure in one doesn't
block the other, and the YAML for each stays simple.
