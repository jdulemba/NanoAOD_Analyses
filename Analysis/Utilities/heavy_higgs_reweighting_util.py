"""
This script is copied from https://gitlab.cern.ch/ljeppe/heavy-higgs-reweighting/-/blob/master/util.py
"""

import subprocess

def get_rel_width(key):
    keyparts = key.split("_")

    if len(keyparts) < 4:
        raise ValueError()

    mass = float(keyparts[1][1:])
    abswidth = float(keyparts[2][1:].replace("p", "."))
    relwidth = abswidth/mass*100

    return relwidth

def name_passes(name, ignore, require):
    relwidth = get_rel_width(name)
    for s in ignore:
        if "relwidth" in s:
            w = float(s.split("relwidth")[1])
            if w == relwidth:
                return False
        elif s in name:
            return False
    for s in require:
        if "relwidth" in s:
            w = float(s.split("relwidth")[1])
            if w != relwidth:
                return False
        elif not s in name:
            return False
    return True

def get_xrootd_path(file):
    if "store/user" in file:
        raise ValueError(f"Not an official sample, so no XrootD: {file}")
    
    return "root://xrootd-cms.infn.it//" + file[file.index("/store/mc"):]

def get_files_from_das(das_name):
    res = subprocess.run(["dasgoclient", "-query", f"file dataset={das_name}"],
        capture_output=True, text=True, check=True)
    files = res.stdout.split("\n")
    return [f for f in files if f != ""]

from dataclasses import dataclass
from typing import Union, Optional

@dataclass(frozen=True)
class Scenario:
    parity: str  # 'A': pseudoscalar, 'H': scalar
    mass: Union[int, float, None]  # Higgs boson mass in GeV
    width: Union[int, float, None]  # Higgs boson realtive decay width
    part: str  # 'res': resonance, 'int'/'intneg'/'intpos': interference
    channel: Optional[str] = None # 'll' or 'lj'
    suffix: Optional[str] = None # Anything, but normally the systematic

    @classmethod
    def fromstr(cls, s: str):
        fields = s.split("_")
        parity, mass, width, part = fields[:4]
        if parity not in ("A", "H"):
            raise ValueError(f"Invalid parity: {parity}")
        if mass[0] != "m":
            raise ValueError(f"Mass is missing prefix 'm': {mass}")
        if mass == "mnat":
            mass = None
        elif "p" in mass:
            mass = float(mass[1:].replace("p", "."))
        else:
            mass = int(mass[1:])
        if not width.startswith("relw") and not width == "wnat":
            raise ValueError(f"Width is missing prefix 'relw': {width}")
        if width == "wnat":
            width = None
        elif "p" in width:
            width = float(width[4:].replace("p", "."))
        else:
            width = int(width[4:])
        if part not in ["res", "int", "intpos", "intneg"]:
            raise ValueError(f"Invalid part {part}")
        if len(fields) > 4:
            channel = fields[4]
            if channel not in ["ll", "lj"]:
                raise ValueError(f"Invalid channel {channel}")
        else:
            channel = "ll"
        if len(fields) > 5:
            suffix = "_".join(fields[5:])
        else:
            suffix = None
        return cls(
            parity=parity,
            mass=mass,
            width=width,
            part=part,
            suffix=suffix,
            channel=channel
        )

    @classmethod
    def from_sample(cls, s: str):
        fields = s.split("_")
        sname, mass, width, part, syst = fields[:5]
        
        if 'Hscalar' in sname:
            parity = 'H'
        elif "Hpseudo" in sname:
            parity = 'A'
        else:
            raise ValueError(f"Invalid parity: {sname}")

        if '2L2Nu' in sname:
            channel = 'll'
        elif "1L1Nu2J" in sname:
            channel = 'lj'
        else:
            raise ValueError(f"Invalid channel: {sname}") 

        if mass[0] != "m":
            raise ValueError(f"Mass is missing prefix 'm': {mass}")
        if "p" in mass:
            mass = float(mass[1:].replace("p", "."))
        else:
            mass = int(mass[1:])
        
        if not width.startswith("w"):
            raise ValueError(f"Width is missing prefix 'w': {width}")
        if "p" in width:
            width = float(width[1:].replace("p", "."))
        else:
            width = int(width[1:])

        width = width / mass * 100

        if part not in ["res", "int", "intpos", "intneg"]:
            raise ValueError(f"Invalid part {part}")
        
        if syst.startswith("mtop"):
            suffix = syst
        else:
            suffix = None

        return cls(
            parity=parity,
            mass=mass,
            width=width,
            part=part,
            channel=channel,
            suffix=suffix
        )

    def __str__(self):
        if self.mass is None:
            mass = "mnat"
        else:
            mass = "m" + str(self.mass).replace(".", "p")
        if self.width is None:
            width = "wnat"
        else:
            width = "relw" + str(self.width).replace(".", "p")
        fields = [
            self.parity,
            mass,
            width,
            self.part,
            self.channel
        ]
        if self.suffix is not None:
            fields.append(self.suffix)
        return "_".join(fields)

    def get_model(self):
        c = 'dilep' if self.channel == 'll' else 'semilep'
        p = 'scalar' if self.parity == 'H' else 'pseudo'
        return f"heavyhiggs-{c}-{p}-{self.part}"

    def get_param_dict(self):
        if self.parity == "H":
            param_dict = {"MH0": self.mass, "H0Width": self.mass*self.width/100}
        elif self.parity == "A":
            param_dict = {"MA0": self.mass, "A0Width": self.mass*self.width/100}
        else:
            raise ValueError(f"Invalid parity: {self.parity}")

        if self.suffix is not None:
            if self.suffix.startswith("mtop"):
                mtop = float(self.suffix[4:].replace("p","."))
                param_dict["MT"] = mtop
                param_dict["ymt"] = mtop

        return param_dict

    def to_sample(self, relwidth=False):
        parity = "pseudo" if self.parity == 'A' else 'scalar'
        channel = "1L1Nu2J" if self.channel == "lj" else "2L2Nu"
        if relwidth:
            width = str(self.width).replace(".", "p")
        else:
            width = str(self.width * self.mass / 100).replace(".", "p")
        return f"H{parity}ToTTTo{channel}_m{int(self.mass)}_w{width}_{self.part}_TuneCP5_13TeV-madgraph_pythia8"


