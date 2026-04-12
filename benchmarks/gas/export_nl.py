"""
Export gas pipeline network NLP benchmark problems to NL format.

Requires: pip install gas_net  (from https://github.com/Sakshi21299/gas_networks)

Produces 4 NL files:
  gaslib11_steady.nl   — GasLib-11 steady-state (n=204, m=200)
  gaslib11_dynamic.nl  — GasLib-11 24h dynamic   (n=2542, m=2492)
  gaslib40_steady.nl   — GasLib-40 steady-state (n=1694, m=1682)
  gaslib40_dynamic.nl  — GasLib-40 24h dynamic   (n=21058, m=20908)

All problems are equality-constrained NLPs from PDE-discretized gas pipeline
networks. See https://github.com/jkitchin/ripopt/issues/12 for context.
"""

import json
from pathlib import Path

from gas_net.util.import_data import import_data_from_excel
from gas_net.model_nlp import buildNonLinearModel
from gas_net.modelling_library.fix_and_init_vars import init_network_default

# Point this at the gas_networks data directory
DATA = Path(__file__).resolve().parent.parent.parent / "saudi-aramco" / "gas_networks" / "gas_net" / "data"
OUT = Path(__file__).resolve().parent

NETWORKS = {
    "GasLib_11": {
        "dir": DATA / "data_files" / "GasLib_11",
        "fix": lambda m: (
            m.wSource["entry01", :].fix(),
            m.wSource["entry02", :].fix(),
            m.pSource["entry03", :].fix(),
        ),
    },
    "Gaslib_40": {
        "dir": DATA / "data_files" / "Gaslib_40",
        "fix": lambda m: (
            m.pSource["source_1", :].fix(),
            m.pSource["source_2", :].fix(),
            m.wSource["source_3", :].fix(),
        ),
    },
}


def export_network(name, cfg, opt):
    tag = name.lower().replace("_", "")  # gaslib11, gaslib40

    nd, id_ = import_data_from_excel(
        str(cfg["dir"] / "networkData.xlsx"),
        str(cfg["dir"] / "inputData.xlsx"),
    )

    # Steady-state
    opt_ss = {**opt, "dynamic": False, "T": opt["T0"] + opt["dt"] / 3600}
    m_ss = buildNonLinearModel(nd, id_, opt_ss)
    m_ss = init_network_default(m_ss, p_default=55e5)
    cfg["fix"](m_ss)
    fname = OUT / f"{tag}_steady.nl"
    m_ss.write(str(fname), format="nl")
    print(f"  {fname.name}")

    # Dynamic (24h horizon)
    nd, id_ = import_data_from_excel(
        str(cfg["dir"] / "networkData.xlsx"),
        str(cfg["dir"] / "inputData.xlsx"),
    )
    opt_dyn = {**opt, "dynamic": True, "T": 24}
    m_dyn = buildNonLinearModel(nd, id_, opt_dyn)
    m_dyn = init_network_default(m_dyn, p_default=55e5)
    cfg["fix"](m_dyn)
    # Fix initial time step pressures from steady-state solution
    t0 = m_dyn.Times.first()
    for p, vol in m_dyn.Pipes_VolExtrR_interm.data():
        m_dyn.interm_p[p, vol, t0] = m_ss.interm_p[p, vol, t0]
    m_dyn.interm_p[:, :, t0].fix()
    fname = OUT / f"{tag}_dynamic.nl"
    m_dyn.write(str(fname), format="nl")
    print(f"  {fname.name}")


def main():
    opt = json.load(open(DATA / "Options.json"))
    for name, cfg in NETWORKS.items():
        print(f"{name}:")
        export_network(name, cfg, opt)
    print("Done.")


if __name__ == "__main__":
    main()
