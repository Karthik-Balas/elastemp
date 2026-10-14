"""vasprun.xml-based convergence validation, tested with small but real vasprun.xml files."""

import pytest

from elastemp.base.strain_analysis import check_convergence
from elastemp.validation.convergence import validate_vasp_run


def _vasprun(n_scsteps, *, nelm=60, nsw=0, ibrion=-1, n_ionic=1, incar_extra="", energy=-71.78294695):
    """A minimal vasprun.xml with the elements VASP writes and the validator reads."""
    scstep = '  <scstep><energy><i name="e_fr_energy"> -71.7 </i></energy></scstep>\n'
    calc = (
        " <calculation>\n"
        + scstep * n_scsteps
        + f'  <energy><i name="e_fr_energy"> {energy} </i><i name="e_0_energy"> {energy} </i></energy>\n'
        " </calculation>\n"
    )
    return (
        '<?xml version="1.0" encoding="ISO-8859-1"?>\n<modeling>\n'
        f' <incar>\n  <i type="int" name="IBRION"> {ibrion} </i>\n{incar_extra} </incar>\n'
        " <parameters>\n"
        f'  <separator name="electronic">\n   <i type="int" name="NELM"> {nelm} </i>\n  </separator>\n'
        f'  <separator name="ionic">\n   <i type="int" name="NSW"> {nsw} </i>\n'
        f'   <i type="int" name="IBRION"> {ibrion} </i>\n  </separator>\n'
        '  <separator name="response functions">\n   <i type="int" name="NELM"> 1 </i>\n  </separator>\n'
        " </parameters>\n" + calc * n_ionic + "</modeling>\n"
    )


def _write(tmp_path, text):
    (tmp_path / "vasprun.xml").write_text(text)
    return tmp_path


def test_converged_static_run(tmp_path):
    result = validate_vasp_run(_write(tmp_path, _vasprun(11)))
    assert result.converged and result.electronic and result.ionic
    assert result.energy == pytest.approx(-71.78294695)
    assert check_convergence(tmp_path / "OSZICAR")[0] is True


def test_scf_that_hit_nelm_is_rejected(tmp_path):
    # VASP still writes a final energy (and an OSZICAR F= line) when NELM is reached.
    result = validate_vasp_run(_write(tmp_path, _vasprun(60, nelm=60)))
    assert not result.converged and not result.electronic
    assert result.energy is None
    assert "NELM = 60" in result.reason


def test_nelm_is_read_from_the_electronic_block(tmp_path):
    # The response-function block also has NELM (= 1); using it would reject every run.
    assert validate_vasp_run(_write(tmp_path, _vasprun(27, nelm=60))).converged


def test_relaxation_needs_ionic_convergence_only_when_asked(tmp_path):
    _write(tmp_path, _vasprun(10, nsw=5, ibrion=2, n_ionic=5))  # used all NSW steps
    assert validate_vasp_run(tmp_path).converged
    strict = validate_vasp_run(tmp_path, require_ionic=True)
    assert not strict.converged and "Ionic" in strict.reason


def test_missing_xml_is_rejected(tmp_path):
    result = validate_vasp_run(tmp_path)
    assert not result.converged and "Missing" in result.reason
    ok, detail = check_convergence(tmp_path / "OSZICAR")
    assert ok is False and "Missing" in detail


def test_truncated_xml_is_rejected(tmp_path):
    # What a job killed at walltime leaves behind: vasprun.xml cut off mid-file.
    full = _vasprun(11)
    result = validate_vasp_run(_write(tmp_path, full[: len(full) // 2]))
    assert not result.converged and "Truncated" in result.reason


def test_fortran_numbers_without_exponent_letter_are_accepted(tmp_path):
    # DFPT output can contain numbers such as 3.37026300-199 that pymatgen cannot parse.
    text = _vasprun(27).replace("</modeling>", " <dynmat><v> 3.37026300-199 -1.06E-02 </v></dynmat>\n</modeling>")
    assert validate_vasp_run(_write(tmp_path, text)).converged


def test_run_without_ionic_step_is_rejected(tmp_path):
    text = _vasprun(11, n_ionic=0)
    result = validate_vasp_run(_write(tmp_path, text))
    assert not result.converged and "No completed ionic step" in result.reason
